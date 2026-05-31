"""EBI GWAS Catalog discovery: query publications by PMID, emit annotation TSV."""

import sys
from typing import Dict, List, Optional

try:
    import requests as _requests
    _REQUESTS_AVAILABLE = True
except ImportError:
    _REQUESTS_AVAILABLE = False


_EBI_API_BASE = "https://www.ebi.ac.uk/gwas/api/v2"
_EBI_PAGE_SIZE = 100

DISCOVERY_TSV_COLUMNS = [
    "pmid", "gcst_id", "file_path", "trait_name",
    "gene", "trait_chr", "trait_bp", "context",
]


def _require_requests():
    if not _REQUESTS_AVAILABLE:
        raise ImportError("requests library is required. Install: pip install requests")


def _fetch_page(pmid: str, page: int, size: int = _EBI_PAGE_SIZE) -> dict:
    _require_requests()
    url = f"{_EBI_API_BASE}/publications/{pmid}/studies/"
    params = {"size": size, "page": page}
    try:
        resp = _requests.get(url, params=params, timeout=30)
    except _requests.RequestException as e:
        raise RuntimeError(f"EBI GWAS API request failed: {e}") from e
    if resp.status_code != 200:
        raise RuntimeError(
            f"EBI GWAS API returned HTTP {resp.status_code} for PMID {pmid}: "
            f"{resp.text[:300]}"
        )
    return resp.json()


def query_ebi_studies(pmid: str) -> List[Dict]:
    """Fetch all studies for a PMID from EBI GWAS Catalog v2 API, paginating as needed."""
    all_studies: List[Dict] = []
    page = 0
    while True:
        data = _fetch_page(pmid, page)
        embedded = data.get("_embedded", {})
        studies = embedded.get("studies", embedded.get("gwasStudies", []))
        all_studies.extend(studies)
        page_info = data.get("page", {})
        total_pages = page_info.get("totalPages", 1)
        current = page_info.get("number", 0)
        if current + 1 >= total_pages or not studies:
            break
        page += 1
    return all_studies


def _harmonised_url(study: Dict) -> str:
    summary_stats = study.get("summaryStatistics", "")
    accession_id = study.get("accessionId", "")
    if not summary_stats or not accession_id:
        return ""
    return f"{summary_stats}/harmonised/{accession_id}.h.tsv.gz"


def _parse_gene_symbol(trait_name: str) -> str:
    parts = trait_name.split() if trait_name else []
    return parts[0] if parts else ""


def load_gene_annotation(path: str) -> Dict[str, Dict]:
    """Load gene annotation TSV (columns: gene_name, chr, bp) keyed by gene_name."""
    from pathlib import Path
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Gene annotation file not found: {path}")
    gene_map: Dict[str, Dict] = {}
    with open(p, 'r') as fh:
        header = fh.readline().rstrip('\n').split('\t')
        col_idx = {n.strip(): i for i, n in enumerate(header)}
        for col in ('gene_name', 'chr', 'bp'):
            if col not in col_idx:
                raise ValueError(f"Gene annotation TSV missing required column: '{col}'")
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if not line.strip():
                continue
            name = parts[col_idx['gene_name']].strip() if col_idx['gene_name'] < len(parts) else ''
            if name:
                gene_map[name] = {
                    'chr': parts[col_idx['chr']].strip() if col_idx['chr'] < len(parts) else '',
                    'bp': parts[col_idx['bp']].strip() if col_idx['bp'] < len(parts) else '',
                }
    return gene_map


def discover_study(
    pmid: str,
    parse_gene: bool = False,
    gene_annotation_path: Optional[str] = None,
) -> List[Dict]:
    """Query EBI GWAS Catalog for all studies for pmid. Return list of row dicts."""
    if parse_gene and not gene_annotation_path:
        raise ValueError("--parse-gene requires --gene-annotation <path>")

    gene_map: Dict[str, Dict] = {}
    if parse_gene and gene_annotation_path:
        gene_map = load_gene_annotation(gene_annotation_path)

    studies = query_ebi_studies(pmid)
    rows = []
    for study in studies:
        gcst_id = study.get("accessionId", "")
        file_path = _harmonised_url(study)
        trait_name = study.get("reportedTrait", "")
        gene = ""
        trait_chr = ""
        trait_bp = ""

        if parse_gene:
            symbol = _parse_gene_symbol(trait_name)
            if symbol:
                if symbol in gene_map:
                    gene = symbol
                    trait_chr = gene_map[symbol].get("chr", "")
                    trait_bp = gene_map[symbol].get("bp", "")
                else:
                    print(
                        f"WARNING: gene '{symbol}' from '{trait_name}' not found in annotation",
                        file=sys.stderr,
                    )

        rows.append({
            "pmid": pmid,
            "gcst_id": gcst_id,
            "file_path": file_path,
            "trait_name": trait_name,
            "gene": gene,
            "trait_chr": trait_chr,
            "trait_bp": trait_bp,
            "context": "",
        })
    return rows


def write_discovery_tsv(rows: List[Dict], output: Optional[str] = None) -> None:
    """Write discovery rows as TSV to stdout or a file path."""
    import sys
    fh = open(output, 'w') if output else sys.stdout
    try:
        fh.write('\t'.join(DISCOVERY_TSV_COLUMNS) + '\n')
        for row in rows:
            fh.write(
                '\t'.join(str(row.get(col, "")) for col in DISCOVERY_TSV_COLUMNS) + '\n'
            )
    finally:
        if output:
            fh.close()
