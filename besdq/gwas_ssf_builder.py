"""Two-pass GWAS-SSF builder: parallel Pass 1 (filtering) + serial Pass 2 (index)."""

import datetime
import gzip
import json
import sqlite3
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from .annotation_reader import TraitConfig
from .gwas_ssf_fast_reader import _count_data_lines, read_gwas_ssf_candidates
from .gwas_ssf_reader import GwasSsfRow
from .significance_filter import apply_significance_filter


@dataclass
class _TraitResult:
    """Intermediate result from Pass 1 for one trait."""
    trait: TraitConfig
    rows: List[GwasSsfRow] = field(default_factory=list)
    cis_rows: List[GwasSsfRow] = field(default_factory=list)
    n_total_read: int = 0
    n_retained: int = 0
    estimated_trait_var: Optional[float] = None
    n_cis: int = 0
    n_sig_trans_peaks: int = 0
    n_sig_trans_peaks_approx: bool = False  # True when distance-based fallback was used
    n_sug_trans: int = 0


def _pass1_worker(args: tuple) -> _TraitResult:
    """Pass 1: stream and filter one GWAS-SSF file."""
    (
        file_path, trait_id, trait_name, trait_chr, trait_bp,
        sample_size, trait_var, gene, context, study_metadata,
        cis_radius, sig_threshold, sug_threshold, plink2_pfile,
        sig_radius, clump_r2, clump_kb, force_fallback,
    ) = args

    trait = TraitConfig(
        file_path=file_path,
        trait_id=trait_id,
        trait_name=trait_name,
        trait_chr=trait_chr,
        trait_bp=trait_bp,
        sample_size=sample_size,
        trait_var=trait_var,
        gene=gene,
        context=context,
        study_metadata=study_metadata,
    )

    result = _TraitResult(trait=trait)

    # Count total rows cheaply (wc -l rather than reading the whole file)
    result.n_total_read = _count_data_lines(file_path)

    # Stream only the candidates we'll keep (cis unconditional + p < sug_threshold)
    cis_start = trait_bp - cis_radius if (trait_chr and trait_bp is not None) else None
    cis_end = trait_bp + cis_radius if (trait_chr and trait_bp is not None) else None
    candidates = list(read_gwas_ssf_candidates(
        file_path,
        cis_chr=trait_chr,
        cis_start=cis_start,
        cis_end=cis_end,
        p_threshold=sug_threshold,
        force_fallback=force_fallback,
    ))

    filter_result = apply_significance_filter(
        candidates,
        trait_chr=trait_chr,
        trait_bp=trait_bp,
        cis_radius=cis_radius,
        sig_threshold=sig_threshold,
        sug_threshold=sug_threshold,
    )

    cis_rows = list(filter_result.cis)
    sug_trans_rows = list(filter_result.sug_trans)
    retained = cis_rows + sug_trans_rows

    # LD clumping for significant trans candidates
    sig_trans_peaks = 0
    sig_trans_peaks_approx = False
    if filter_result.sig_trans_candidates and plink2_pfile:
        try:
            from .ld_clumping import clump_trans_peaks
            clump_result = clump_trans_peaks(
                filter_result.sig_trans_candidates,
                plink2_pfile=plink2_pfile,
                sig_radius=sig_radius,
                clump_r2=clump_r2,
                clump_kb=clump_kb,
                all_rows=candidates,
            )
            retained.extend(clump_result.rows)
            sig_trans_peaks = clump_result.peak_count
        except ImportError:
            retained.extend(filter_result.sig_trans_candidates)
            from .peak_count import count_peaks_by_distance
            sig_trans_peaks = count_peaks_by_distance(
                filter_result.sig_trans_candidates, sig_radius
            )
            sig_trans_peaks_approx = True
    elif filter_result.sig_trans_candidates:
        retained.extend(filter_result.sig_trans_candidates)
        from .peak_count import count_peaks_by_distance
        sig_trans_peaks = count_peaks_by_distance(
            filter_result.sig_trans_candidates, sig_radius
        )
        sig_trans_peaks_approx = True

    # Deduplicate by snp_key
    seen: set = set()
    deduped = []
    for row in retained:
        key = row.snp_key
        if key not in seen:
            seen.add(key)
            deduped.append(row)

    result.rows = deduped
    result.cis_rows = cis_rows
    result.n_retained = len(deduped)
    result.n_cis = len(cis_rows)
    result.n_sig_trans_peaks = sig_trans_peaks
    result.n_sig_trans_peaks_approx = sig_trans_peaks_approx
    result.n_sug_trans = len(sug_trans_rows)

    # Estimate trait_var from cis SNPs: median(se^2 * n * 2 * eaf * (1-eaf))
    if cis_rows and sample_size:
        ses_arr = np.array([r.se for r in cis_rows], dtype=np.float64)
        eafs_arr = np.array([r.eaf for r in cis_rows], dtype=np.float64)
        valid = (eafs_arr > 0.01) & (eafs_arr < 0.99) & (ses_arr > 0)
        if valid.sum() >= 10:
            result.estimated_trait_var = float(
                np.median(ses_arr[valid] ** 2 * sample_size * 2
                          * eafs_arr[valid] * (1.0 - eafs_arr[valid]))
            )

    return result


def _snp_sort_key(snp_key: str) -> Tuple:
    """Sort SNP keys by chromosome (numeric then alpha) then position."""
    parts = snp_key.split(':')
    chr_str = parts[0]
    bp = int(parts[1]) if len(parts) > 1 else 0
    try:
        return (0, int(chr_str), bp)
    except ValueError:
        return (1, chr_str, bp)


# Columns written to intermediate TSV.gz files by Stage 1
INTERMEDIATE_TSV_COLUMNS = ['chr', 'bp', 'a1', 'a2', 'rsid', 'eaf', 'beta', 'se', 'p']


def read_intermediate_tsv(tsv_gz_path: str) -> List[GwasSsfRow]:
    """Read filtered rows from an intermediate TSV.gz produced by Stage 1."""
    rows = []
    with gzip.open(tsv_gz_path, 'rt') as fh:
        header = fh.readline().rstrip('\n').split('\t')
        col_idx = {n: i for i, n in enumerate(header)}
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if not line.strip():
                continue
            try:
                rows.append(GwasSsfRow(
                    chr=parts[col_idx['chr']],
                    bp=int(parts[col_idx['bp']]),
                    a1=parts[col_idx['a1']],
                    a2=parts[col_idx['a2']],
                    rsid=parts[col_idx['rsid']] or None,
                    eaf=float(parts[col_idx['eaf']]),
                    beta=float(parts[col_idx['beta']]),
                    se=float(parts[col_idx['se']]),
                    p=float(parts[col_idx['p']]) if 'p' in col_idx else 0.0,
                ))
            except (ValueError, IndexError, KeyError):
                continue
    return rows


class GwasSsfIndexBuilder:
    """Build a BESD-compatible SQLite index from GWAS-SSF files."""

    def __init__(self, db_path: str):
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)

    def build(
        self,
        traits: List[TraitConfig],
        workers: int = 1,
        cis_radius: int = 1_000_000,
        sig_threshold: float = 5e-8,
        sug_threshold: float = 1e-4,
        plink2_pfile: Optional[str] = None,
        sig_radius: int = 500_000,
        clump_r2: float = 0.01,
        clump_kb: int = 10_000,
        force_fallback: bool = False,
    ) -> None:
        """Build the index database from annotation-driven trait configs."""
        if self.db_path.exists():
            self.db_path.unlink()

        # ------- Pass 1: parallel per-file filtering -------
        _log("Pass 1: filtering traits…")
        worker_args = [
            (
                t.file_path, t.trait_id, t.trait_name, t.trait_chr, t.trait_bp,
                t.sample_size, t.trait_var, t.gene, t.context, t.study_metadata,
                cis_radius, sig_threshold, sug_threshold, plink2_pfile,
                sig_radius, clump_r2, clump_kb, force_fallback,
            )
            for t in traits
        ]

        trait_results: Dict[str, _TraitResult] = {}

        if workers == 1:
            for args in worker_args:
                r = _pass1_worker(args)
                trait_results[r.trait.trait_id] = r
                _log(f"  {r.trait.trait_id}: {r.n_retained}/{r.n_total_read} rows retained")
        else:
            with ProcessPoolExecutor(max_workers=workers) as pool:
                futures = {pool.submit(_pass1_worker, args): args[1] for args in worker_args}
                for fut in as_completed(futures):
                    r = fut.result()
                    trait_results[r.trait.trait_id] = r
                    _log(f"  {r.trait.trait_id}: {r.n_retained}/{r.n_total_read} rows retained")

        build_params = {
            'cis_radius': cis_radius,
            'sig_threshold': sig_threshold,
            'sug_threshold': sug_threshold,
            'sig_radius': sig_radius,
        }

        # ------- Pass 2: serial index construction -------
        conn = sqlite3.connect(str(self.db_path))
        cursor = conn.cursor()
        self._create_schema(cursor)
        self._run_pass2(conn, cursor, traits, trait_results, build_params)
        conn.close()

    def build_from_intermediates(
        self,
        workdir: str,
    ) -> None:
        """Build SQLite index from Stage 1 intermediate file pairs in workdir.

        Discovers all GCST*.yaml (skips *.failed.yaml), loads paired TSV.gz,
        and writes a clean SQLite index at self.db_path.
        """
        import yaml

        workdir_path = Path(workdir)
        if not workdir_path.exists():
            raise FileNotFoundError(f"Workdir not found: {workdir}")

        yaml_files = sorted(
            p for p in workdir_path.glob("*.yaml")
            if not p.name.endswith(".failed.yaml")
        )
        if not yaml_files:
            raise ValueError(f"No completed intermediate YAML files found in {workdir}")

        traits: List[TraitConfig] = []
        trait_results: Dict[str, _TraitResult] = {}
        build_params: Optional[Dict] = None

        for yaml_path in yaml_files:
            with open(yaml_path, 'r') as fh:
                meta = yaml.safe_load(fh) or {}

            trait_id = meta.get('trait_id', yaml_path.stem)
            tsv_gz_path = yaml_path.parent / f"{trait_id}.tsv.gz"
            if not tsv_gz_path.exists():
                _log(f"  WARNING: missing TSV for {trait_id}, skipping")
                continue

            trait_cfg = TraitConfig(
                file_path=str(tsv_gz_path),
                trait_id=trait_id,
                trait_name=meta.get('trait_name', trait_id),
                trait_chr=meta.get('trait_chr') or None,
                trait_bp=meta.get('trait_bp') or None,
                sample_size=None,
                trait_var=meta.get('trait_var') or None,
                gene=meta.get('gene') or None,
                context=meta.get('context') or None,
                study_metadata={},
            )
            traits.append(trait_cfg)

            rows = read_intermediate_tsv(str(tsv_gz_path))
            result = _TraitResult(
                trait=trait_cfg,
                rows=rows,
                n_total_read=meta.get('n_total_read', len(rows)),
                n_retained=len(rows),
                estimated_trait_var=None,
                n_cis=meta.get('n_cis', 0),
                n_sig_trans_peaks=meta.get('n_sig_trans_peaks', 0),
                n_sig_trans_peaks_approx=bool(meta.get('n_sig_trans_peaks_approx', False)),
                n_sug_trans=meta.get('n_sug_trans', 0),
            )
            trait_results[trait_id] = result

            if build_params is None:
                build_params = {
                    'cis_radius': meta.get('cis_radius', 1_000_000),
                    'sig_threshold': meta.get('sig_threshold', 5e-8),
                    'sug_threshold': meta.get('sug_threshold', 1e-4),
                    'sig_radius': meta.get('sig_radius', 500_000),
                }

        if not traits:
            raise ValueError("No valid traits loaded from intermediates")

        _log(f"Stage 2: {len(traits)} traits loaded from intermediates")

        if self.db_path.exists():
            self.db_path.unlink()

        conn = sqlite3.connect(str(self.db_path))
        cursor = conn.cursor()
        self._create_schema(cursor)
        self._run_pass2(conn, cursor, traits, trait_results, build_params or {})
        conn.close()
        _log(f"Stage 2 complete. Database written to {self.db_path}")

    def _run_pass2(
        self,
        conn: sqlite3.Connection,
        cursor: sqlite3.Cursor,
        traits: List[TraitConfig],
        trait_results: Dict[str, _TraitResult],
        build_params: Dict,
    ) -> None:
        """Serial index construction: ESI, EPI, probe_data, and besd_meta."""
        _log("Pass 2: building index…")

        cis_radius = build_params.get('cis_radius', 1_000_000)
        sig_threshold = build_params.get('sig_threshold', 5e-8)
        sug_threshold = build_params.get('sug_threshold', 1e-4)
        sig_radius = build_params.get('sig_radius', 500_000)

        # Consolidate SNP universe across all traits, sorted stably
        all_snp_keys: set = set()
        for r in trait_results.values():
            for row in r.rows:
                all_snp_keys.add(row.snp_key)

        sorted_keys = sorted(all_snp_keys, key=_snp_sort_key)
        snp_key_to_idx: Dict[str, int] = {k: i for i, k in enumerate(sorted_keys)}
        _log(f"  ESI size: {len(sorted_keys)} unique SNPs")

        # Write ESI
        snp_data: Dict[str, dict] = {}
        for r in trait_results.values():
            for row in r.rows:
                key = row.snp_key
                if key not in snp_data:
                    snp_data[key] = {
                        'chr': row.chr,
                        'bp': row.bp,
                        'a1': row.a1,
                        'a2': row.a2,
                        'rsid': row.rsid,
                        'eaf': row.eaf,
                    }

        for idx, key in enumerate(sorted_keys):
            d = snp_data[key]
            cursor.execute(
                "INSERT INTO esi (row_idx, chr, snp_id, genetic_dist, bp, a1, a2, freq) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                (idx, d['chr'], d['rsid'], None, d['bp'], d['a1'], d['a2'], d['eaf']),
            )

        _log("  ESI written")

        # Write EPI and probe_data
        study_meta_written = False
        for epi_idx, trait_cfg in enumerate(traits):
            tid = trait_cfg.trait_id
            r = trait_results.get(tid)
            rows = r.rows if r else []

            # Resolve trait_var: user-supplied > estimated from cis SNPs > None
            resolved_trait_var = trait_cfg.trait_var
            if resolved_trait_var is None and r is not None:
                resolved_trait_var = r.estimated_trait_var
            if resolved_trait_var is None:
                _log(f"  WARNING: no trait_var for {tid} — queries will return original units")

            cursor.execute(
                "INSERT INTO epi (row_idx, trait_id, trait_name, trait_chr, trait_bp, "
                "trait_var, gene, context, n_source_snps, n_cis, n_sig_trans_peaks, "
                "n_sig_trans_peaks_approx, n_sug_trans) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (
                    epi_idx, trait_cfg.trait_id, trait_cfg.trait_name,
                    trait_cfg.trait_chr, trait_cfg.trait_bp,
                    resolved_trait_var,
                    trait_cfg.gene, trait_cfg.context,
                    r.n_total_read if r else None,
                    r.n_cis if r else None,
                    r.n_sig_trans_peaks if r else None,
                    int(r.n_sig_trans_peaks_approx) if r else None,
                    r.n_sug_trans if r else None,
                ),
            )

            if rows:
                snp_indices = np.array(
                    [snp_key_to_idx[row.snp_key] for row in rows], dtype=np.int32
                )
                zscores = np.array(
                    [row.beta / row.se if row.se > 0 else 0.0 for row in rows],
                    dtype=np.float64,
                ).astype(np.float16)
                # VectorN: store SE directly (original units, float16)
                se_vector = np.array(
                    [row.se for row in rows], dtype=np.float64
                ).astype(np.float16)
            else:
                snp_indices = np.array([], dtype=np.int32)
                zscores = np.array([], dtype=np.float16)
                se_vector = np.array([], dtype=np.float16)

            cursor.execute(
                "INSERT INTO probe_data (probe_idx, snp_count, snp_indices, zscores, "
                "n_scalar, se_vector) VALUES (?, ?, ?, ?, ?, ?)",
                (
                    epi_idx, len(rows),
                    snp_indices.tobytes(), zscores.tobytes(),
                    None, se_vector.tobytes(),
                ),
            )

            # Store study metadata once (first trait with metadata)
            if not study_meta_written and trait_cfg.study_metadata:
                cursor.execute(
                    "INSERT INTO besd_meta (key, value) VALUES (?, ?)",
                    ('study_metadata', json.dumps(trait_cfg.study_metadata)),
                )
                study_meta_written = True

        # Store basic metadata
        for key, value in [
            ('n_snps', str(len(sorted_keys))),
            ('n_traits', str(len(traits))),
            ('source', 'gwas-ssf'),
            ('cis_radius', str(cis_radius)),
            ('sig_threshold', str(sig_threshold)),
            ('sug_threshold', str(sug_threshold)),
            ('sig_radius', str(sig_radius)),
        ]:
            cursor.execute(
                "INSERT INTO besd_meta (key, value) VALUES (?, ?)", (key, value)
            )

        # Indexes
        cursor.execute("CREATE INDEX idx_esi_chr_bp ON esi(chr, bp)")
        cursor.execute("CREATE INDEX idx_esi_snp_id ON esi(snp_id)")
        cursor.execute("CREATE INDEX idx_epi_trait_chr_bp ON epi(trait_chr, trait_bp)")
        cursor.execute("CREATE INDEX idx_epi_trait_id ON epi(trait_id)")

        conn.commit()
        _log(f"Pass 2 complete. Database written to {self.db_path}")

    def _create_schema(self, cursor: sqlite3.Cursor) -> None:
        cursor.execute("""
            CREATE TABLE besd_meta (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL
            )
        """)
        cursor.execute("""
            CREATE TABLE esi (
                row_idx INTEGER PRIMARY KEY,
                chr TEXT NOT NULL,
                snp_id TEXT,
                genetic_dist REAL,
                bp INTEGER NOT NULL,
                a1 TEXT,
                a2 TEXT,
                freq REAL
            )
        """)
        cursor.execute("""
            CREATE TABLE epi (
                row_idx INTEGER PRIMARY KEY,
                trait_id TEXT NOT NULL,
                trait_name TEXT NOT NULL,
                trait_chr TEXT,
                trait_bp INTEGER,
                trait_var REAL,
                gene TEXT,
                context TEXT,
                n_source_snps INTEGER,
                n_cis INTEGER,
                n_sig_trans_peaks INTEGER,
                n_sig_trans_peaks_approx INTEGER,
                n_sug_trans INTEGER
            )
        """)
        cursor.execute("""
            CREATE TABLE probe_data (
                probe_idx INTEGER PRIMARY KEY,
                snp_count INTEGER NOT NULL,
                snp_indices BLOB NOT NULL,
                zscores BLOB NOT NULL,
                n_scalar INTEGER,
                se_vector BLOB
            )
        """)


def _log(msg: str) -> None:
    ts = datetime.datetime.now().strftime('%H:%M:%S')
    print(f"[{ts}] {msg}", file=sys.stderr)
