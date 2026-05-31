"""Stage 1: filter GWAS-SSF source files and write intermediate TSV+YAML pairs."""

import datetime
import gzip as _gzip
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import yaml

from .gwas_ssf_builder import (
    INTERMEDIATE_TSV_COLUMNS,
    _TraitResult,
    _pass1_worker,
)
from .gwas_ssf_reader import GwasSsfRow
from .annotation_reader import TraitConfig


# ---------------------------------------------------------------------------
# Discovery TSV reader
# ---------------------------------------------------------------------------

@dataclass
class DiscoveryRow:
    pmid: str
    gcst_id: str
    file_path: str
    trait_name: str
    gene: str
    trait_chr: Optional[str]
    trait_bp: Optional[int]
    context: str
    sample_size: Optional[int] = None


def read_discovery_tsv(tsv_path: str) -> List[DiscoveryRow]:
    """Read discovery TSV produced by discover-study."""
    path = Path(tsv_path)
    if not path.exists():
        raise FileNotFoundError(f"Discovery TSV not found: {tsv_path}")

    rows: List[DiscoveryRow] = []
    with open(path, 'r') as fh:
        header = fh.readline().rstrip('\n').split('\t')
        col_idx = {n.strip(): i for i, n in enumerate(header)}

        required = ('pmid', 'gcst_id', 'file_path', 'trait_name')
        missing = [c for c in required if c not in col_idx]
        if missing:
            raise ValueError(f"Discovery TSV missing required columns: {missing}")

        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')

            def _get(col: str) -> str:
                i = col_idx.get(col)
                return parts[i].strip() if i is not None and i < len(parts) else ''

            chr_raw = _get('trait_chr')
            bp_raw = _get('trait_bp')
            trait_chr = chr_raw if chr_raw else None
            trait_bp = int(bp_raw) if bp_raw and bp_raw.lstrip('-').isdigit() else None

            ss_raw = _get('sample_size')
            sample_size = int(ss_raw) if ss_raw and ss_raw.isdigit() else None

            rows.append(DiscoveryRow(
                pmid=_get('pmid'),
                gcst_id=_get('gcst_id'),
                file_path=_get('file_path'),
                trait_name=_get('trait_name'),
                gene=_get('gene'),
                trait_chr=trait_chr,
                trait_bp=trait_bp,
                context=_get('context'),
                sample_size=sample_size,
            ))
    return rows


# ---------------------------------------------------------------------------
# Source locator helpers
# ---------------------------------------------------------------------------

def _is_url(path: str) -> bool:
    return path.startswith(('http://', 'https://', 'ftp://'))


def _download_file(url: str, dest: Path) -> None:
    try:
        import requests
    except ImportError:
        raise ImportError(
            "requests library required for URL downloads. Install: pip install requests"
        )
    resp = requests.get(url, stream=True, timeout=120)
    if resp.status_code != 200:
        raise RuntimeError(f"Download failed for {url}: HTTP {resp.status_code}")
    with open(dest, 'wb') as fh:
        for chunk in resp.iter_content(chunk_size=1 << 16):
            fh.write(chunk)


# ---------------------------------------------------------------------------
# Intermediate file writing
# ---------------------------------------------------------------------------

def _write_intermediate_pair(
    rows: List[GwasSsfRow],
    meta: dict,
    tsv_gz_path: Path,
    yaml_path: Path,
) -> None:
    """Atomically write filtered rows to TSV.gz and metadata to YAML."""
    tsv_tmp = Path(str(tsv_gz_path) + '.tmp')
    yaml_tmp = Path(str(yaml_path) + '.tmp')
    try:
        with _gzip.open(tsv_tmp, 'wt') as fh:
            fh.write('\t'.join(INTERMEDIATE_TSV_COLUMNS) + '\n')
            for row in rows:
                fh.write('\t'.join([
                    str(row.chr),
                    str(row.bp),
                    row.a1 or '',
                    row.a2 or '',
                    row.rsid or '',
                    str(row.eaf),
                    str(row.beta),
                    str(row.se),
                    str(row.p),
                ]) + '\n')
        with open(yaml_tmp, 'w') as fh:
            yaml.dump(meta, fh, default_flow_style=False, sort_keys=True)
        os.replace(tsv_tmp, tsv_gz_path)
        os.replace(yaml_tmp, yaml_path)
    except Exception:
        tsv_tmp.unlink(missing_ok=True)
        yaml_tmp.unlink(missing_ok=True)
        raise


def _write_failed_yaml(path: Path, gcst_id: str, error: str) -> None:
    data = {
        'gcst_id': gcst_id,
        'error': error,
        'timestamp': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    }
    with open(path, 'w') as fh:
        yaml.dump(data, fh, default_flow_style=False)


# ---------------------------------------------------------------------------
# Per-trait processing (top-level for ProcessPoolExecutor pickling)
# ---------------------------------------------------------------------------

def _process_trait_args(args: tuple) -> Tuple[str, str]:
    """Worker: process one trait and return (gcst_id, status)."""
    (
        gcst_id, pmid, trait_name, trait_chr, trait_bp, gene, context,
        file_path, sample_size, workdir_str,
        cis_radius, sig_threshold, sug_threshold, plink2_pfile,
        sig_radius, clump_r2, clump_kb, force_fallback,
    ) = args

    work_pmid_dir = Path(workdir_str) / pmid
    tsv_gz_path = work_pmid_dir / f"{gcst_id}.tsv.gz"
    yaml_path = work_pmid_dir / f"{gcst_id}.yaml"
    failed_path = work_pmid_dir / f"{gcst_id}.failed.yaml"

    # Already completed
    if tsv_gz_path.exists() and yaml_path.exists():
        return (gcst_id, 'skipped')

    # Delete failed marker so we retry
    if failed_path.exists():
        failed_path.unlink()

    source_is_url = _is_url(file_path)
    local_path = file_path
    downloaded = False

    try:
        if source_is_url:
            tmp_dl = work_pmid_dir / f"{gcst_id}.download.tmp"
            # Set downloaded=True before attempt so finally always cleans up partial writes
            local_path = str(tmp_dl)
            downloaded = True
            _download_file(file_path, tmp_dl)

        pass1_args = (
            local_path,
            gcst_id,
            trait_name,
            trait_chr,
            trait_bp,
            sample_size,
            None,       # trait_var — estimated in pass1
            gene,
            context,
            {},         # study_metadata
            cis_radius,
            sig_threshold,
            sug_threshold,
            plink2_pfile,
            sig_radius,
            clump_r2,
            clump_kb,
            force_fallback,
        )
        result = _pass1_worker(pass1_args)

        meta = {
            'trait_id': gcst_id,
            'trait_name': trait_name,
            'pmid': pmid,
            'trait_chr': trait_chr,
            'trait_bp': trait_bp,
            'gene': gene or None,
            'context': context or None,
            'trait_var': result.estimated_trait_var,
            'n_cis': result.n_cis,
            'n_sig_trans_peaks': result.n_sig_trans_peaks,
            'n_sig_trans_peaks_approx': bool(result.n_sig_trans_peaks_approx),
            'n_sug_trans': result.n_sug_trans,
            'n_total_read': result.n_total_read,
            'cis_radius': cis_radius,
            'sig_threshold': sig_threshold,
            'sug_threshold': sug_threshold,
            'sig_radius': sig_radius,
        }
        _write_intermediate_pair(result.rows, meta, tsv_gz_path, yaml_path)
        return (gcst_id, 'completed')

    except Exception as e:
        try:
            _write_failed_yaml(failed_path, gcst_id, str(e))
        except Exception:
            pass
        return (gcst_id, 'failed')

    finally:
        if downloaded:
            Path(local_path).unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# Stage 1 builder
# ---------------------------------------------------------------------------

class Stage1Builder:
    """Run Stage 1: filter source files and write intermediate pairs."""

    def build(
        self,
        discovery_tsv: str,
        workdir: str,
        parallel_downloads: int = 1,
        cis_radius: int = 1_000_000,
        sig_threshold: float = 5e-8,
        sug_threshold: float = 1e-4,
        plink2_pfile: Optional[str] = None,
        sig_radius: int = 500_000,
        clump_r2: float = 0.01,
        clump_kb: int = 10_000,
        force_fallback: bool = False,
    ) -> Dict[str, int]:
        """Process all traits in discovery_tsv. Returns counts dict."""
        rows = read_discovery_tsv(discovery_tsv)
        if not rows:
            raise ValueError(f"No traits found in {discovery_tsv}")

        # Create per-PMID subdirectories up front
        pmids = {r.pmid for r in rows}
        for pmid in pmids:
            Path(workdir, pmid).mkdir(parents=True, exist_ok=True)

        worker_args = [
            (
                r.gcst_id, r.pmid, r.trait_name, r.trait_chr, r.trait_bp,
                r.gene, r.context, r.file_path, r.sample_size, workdir,
                cis_radius, sig_threshold, sug_threshold, plink2_pfile,
                sig_radius, clump_r2, clump_kb, force_fallback,
            )
            for r in rows
        ]

        counts = {'completed': 0, 'skipped': 0, 'failed': 0}

        if parallel_downloads == 1:
            for args in worker_args:
                gcst_id, status = _process_trait_args(args)
                counts[status] += 1
                _log(f"  {gcst_id}: {status}")
        else:
            with ProcessPoolExecutor(max_workers=parallel_downloads) as pool:
                futures = {
                    pool.submit(_process_trait_args, args): args[0]
                    for args in worker_args
                }
                for fut in as_completed(futures):
                    gcst_id, status = fut.result()
                    counts[status] += 1
                    _log(f"  {gcst_id}: {status}")

        _log(
            f"Stage 1 complete — "
            f"{counts['completed']} completed, "
            f"{counts['skipped']} skipped, "
            f"{counts['failed']} failed"
        )
        return counts


def _log(msg: str) -> None:
    ts = datetime.datetime.now().strftime('%H:%M:%S')
    print(f"[{ts}] {msg}", file=sys.stderr)
