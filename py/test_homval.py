#!/usr/bin/env python3
"""
test_homval.py -- compare hits_to_edf.py / hits_to_topcat.py against brute-force
expected golden files.

Run make_expected_edf.py and make_expected_tcat.py first to regenerate fixtures.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

TEST_DIR = Path(__file__).resolve().parent.parent / "test"
PY_DIR = Path(__file__).resolve().parent

EDF_TRUTHS = ("superfamily", "superfamilyx")
TOPCAT_TRUTHS = ("topsf", "topfold")


def run_hits_to_edf(truth: str, out_path: Path) -> None:
    cmd = [
        sys.executable,
        str(PY_DIR / "hits_to_edf.py"),
        "--hits",
        str(TEST_DIR / "test_hits.tsv"),
        "--lookup",
        str(TEST_DIR / "test.lookup"),
        "--derived-info",
        str(TEST_DIR / "derived_info.json"),
        "--truth",
        truth,
        "--score",
        "--output",
        str(out_path),
    ]
    subprocess.run(cmd, check=True)


def run_hits_to_topcat(truth: str, out_path: Path) -> None:
    cmd = [
        sys.executable,
        str(PY_DIR / "hits_to_topcat.py"),
        "--hits",
        str(TEST_DIR / "test_hits.tsv"),
        "--lookup",
        str(TEST_DIR / "test.lookup"),
        "--derived-info",
        str(TEST_DIR / "derived_info.json"),
        "--truth",
        truth,
        "--score",
        "--output",
        str(out_path),
    ]
    subprocess.run(cmd, check=True)


def read_edf(path: Path) -> tuple[dict, list[str]]:
    hdr: dict[str, str] = {}
    body: list[str] = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                content = line[1:].strip()
                if content.startswith("SEPQ0.1="):
                    for part in content.split():
                        k, v = part.split("=", 1)
                        hdr[k] = v
                elif "=" in content:
                    k, v = content.split("=", 1)
                    hdr[k] = v
            elif line.startswith("score\t"):
                continue
            elif line.strip():
                body.append(line.strip())
    return hdr, body


def read_tcat(path: Path) -> tuple[dict, list[str]]:
    hdr: dict[str, str] = {}
    body: list[str] = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                content = line[1:].strip()
                if content.startswith("TEPQ0.001="):
                    for part in content.split():
                        k, v = part.split("=", 1)
                        hdr[k] = v
                elif "=" in content:
                    k, v = content.split("=", 1)
                    hdr[k] = v
            elif line.startswith("domain\t"):
                continue
            elif line.strip():
                body.append(line.strip())
    return hdr, body


def compare_edf_truth(truth: str) -> bool:
    expected = TEST_DIR / f"expected.{truth}.edf"
    if not expected.exists():
        print(f"{truth}: missing {expected} (run make_expected_edf.py)", file=sys.stderr)
        return False

    got_path = TEST_DIR / f"_test_output.{truth}.edf"
    run_hits_to_edf(truth, got_path)

    exp_hdr, exp_body = read_edf(expected)
    got_hdr, got_body = read_edf(got_path)

    ok = True
    if got_body != exp_body:
        print(f"{truth}: EDF body mismatch", file=sys.stderr)
        print("  got:", got_body, file=sys.stderr)
        print("  expected:", exp_body, file=sys.stderr)
        ok = False
    else:
        print(f"{truth}: EDF body OK ({len(got_body)} rows)")

    for key in ("truth", "fields", "score_direction", "ndom", "N_possible_tp", "N_possible_fp"):
        if got_hdr.get(key) != exp_hdr.get(key):
            print(
                f"{truth}: header {key} mismatch got={got_hdr.get(key)!r} "
                f"expected={exp_hdr.get(key)!r}",
                file=sys.stderr,
            )
            ok = False

    for key in ("SEPQ0.1", "SEPQ1", "SEPQ10", "Sum3", "SFFP"):
        if got_hdr.get(key) != exp_hdr.get(key):
            print(
                f"{truth}: {key} mismatch got={got_hdr.get(key)!r} "
                f"expected={exp_hdr.get(key)!r}",
                file=sys.stderr,
            )
            ok = False
        else:
            print(f"{truth}: {key} OK ({got_hdr.get(key)})")

    return ok


def compare_topcat_truth(truth: str) -> bool:
    expected = TEST_DIR / f"expected.{truth}.tcat"
    if not expected.exists():
        print(f"{truth}: missing {expected} (run make_expected_tcat.py)", file=sys.stderr)
        return False

    got_path = TEST_DIR / f"_test_output.{truth}.tcat"
    run_hits_to_topcat(truth, got_path)

    exp_hdr, exp_body = read_tcat(expected)
    got_hdr, got_body = read_tcat(got_path)

    ok = True
    if got_body != exp_body:
        print(f"{truth}: TCAT body mismatch", file=sys.stderr)
        print("  got:", got_body, file=sys.stderr)
        print("  expected:", exp_body, file=sys.stderr)
        ok = False
    else:
        print(f"{truth}: TCAT body OK ({len(got_body)} rows)")

    for key in ("truth", "fields", "score_direction", "N_possible_tp"):
        if got_hdr.get(key) != exp_hdr.get(key):
            print(
                f"{truth}: header {key} mismatch got={got_hdr.get(key)!r} "
                f"expected={exp_hdr.get(key)!r}",
                file=sys.stderr,
            )
            ok = False

    for key in ("TEPQ0.001", "TEPQ0.01", "TEPQ0.1", "Top3"):
        if got_hdr.get(key) != exp_hdr.get(key):
            print(
                f"{truth}: {key} mismatch got={got_hdr.get(key)!r} "
                f"expected={exp_hdr.get(key)!r}",
                file=sys.stderr,
            )
            ok = False
        else:
            print(f"{truth}: {key} OK ({got_hdr.get(key)})")

    return ok


def main() -> int:
    derived = TEST_DIR / "derived_info.json"
    if not derived.exists():
        print("missing derived_info.json; run build_derived_info.py first", file=sys.stderr)
        return 1

    all_ok = True
    for truth in EDF_TRUTHS:
        if not compare_edf_truth(truth):
            all_ok = False
    for truth in TOPCAT_TRUTHS:
        if not compare_topcat_truth(truth):
            all_ok = False

    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
