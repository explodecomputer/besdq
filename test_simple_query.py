#!/usr/bin/env python3
"""
Test simple_query.py tool.
"""

import tempfile
import struct
import subprocess
import sys
from pathlib import Path
import numpy as np


def create_test_files(tmpdir):
    """Create test BESD files."""
    tmpdir = Path(tmpdir)

    # Create .esi file
    esi_path = tmpdir / "test.esi"
    with open(esi_path, "w") as f:
        f.write("1\trs1\t0.0\t1000\tA\tG\t0.45\n")
        f.write("1\trs2\t0.001\t2000\tC\tT\t0.30\n")
        f.write("1\trs3\t0.002\t3000\tG\tA\t0.55\n")

    # Create .epi file
    epi_path = tmpdir / "test.epi"
    with open(epi_path, "w") as f:
        f.write("1\tENSG00000001\t0.0\t1500\tGENE1\t+\n")
        f.write("1\tENSG00000002\t0.001\t2500\tGENE2\t-\n")

    # Create sparse BESD file
    besd_path = tmpdir / "test.besd"
    n_probes = 2
    n_snps = 3

    with open(besd_path, "wb") as f:
        # Header
        magic = 0x3f800000  # Sparse
        f.write(struct.pack("<IIII", magic, 0, n_probes, n_snps))

        # Placeholder for offset table
        offset_table_pos = f.tell()
        for _ in range(n_probes):
            f.write(struct.pack("<QQ", 0, 0))

        # Probe blocks
        probe_offsets = []
        for probe_idx in range(n_probes):
            start_pos = f.tell()

            snp_count = probe_idx + 1
            snp_indices = np.array([i for i in range(snp_count)], dtype=np.uint32)
            betas = np.array([0.05 * (i + 1) for i in range(snp_count)], dtype=np.float32)
            ses = np.array([0.01 * (i + 1) for i in range(snp_count)], dtype=np.float32)

            f.write(struct.pack("<I", snp_count))
            f.write(snp_indices.tobytes())
            f.write(betas.tobytes())
            f.write(ses.tobytes())

            end_pos = f.tell()
            probe_offsets.append((start_pos, end_pos))

        # Update offset table
        f.seek(offset_table_pos)
        for start, end in probe_offsets:
            f.write(struct.pack("<QQ", start, end))

    return tmpdir / "test"


def main():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        prefix = create_test_files(tmpdir)

        print("=" * 60)
        print("Testing simple_query.py")
        print("=" * 60)

        # Test query
        cmd = [
            sys.executable, "simple_query.py",
            "--beqtl-summary", str(prefix),
            "--snp-chr", "1",
            "--from-snp-kb", "0.5",
            "--to-snp-kb", "2.5",
            "--probe-chr", "1",
            "--from-probe-kb", "1.0",
            "--to-probe-kb", "2.7",
            "--out", str(tmpdir / "result"),
        ]

        print(f"\nRunning: {' '.join(cmd)}\n")
        result = subprocess.run(cmd, cwd=Path(__file__).parent)

        if result.returncode == 0:
            # Check output
            output_file = tmpdir / "result.txt"
            if output_file.exists():
                print(f"\nOutput file contents:")
                with open(output_file) as f:
                    print(f.read())
                print("✅ Test passed!")
            else:
                print("❌ Output file not created")
        else:
            print(f"❌ Command failed with exit code {result.returncode}")


if __name__ == "__main__":
    main()
