#!/usr/bin/env python3
"""Recompute the finite certificate for the planar Lennard--Jones theorem.

MATHEMATICAL_SCOPE.txt describes the analytic premises and finite assertions.
"""
if not __debug__:
    raise RuntimeError("Assertions must be enabled; do not use -O or -OO.")

import argparse
import hashlib
import json
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
UPSTREAM_COMMIT = "fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_sources():
    """Identify every source file and confirm the unmodified upstream copy."""
    records = []
    for line in (ROOT / "SOURCE_SHA256SUMS.txt").read_text().splitlines():
        if not line.strip():
            continue
        expected, relative = line.split("  ", 1)
        source = (ROOT / relative).resolve()
        if not source.is_relative_to(ROOT):
            raise ValueError("A source-manifest path escapes the certificate directory")
        actual = digest(source)
        if actual != expected:
            raise ValueError(f"Source digest mismatch: {relative}")
        records.append({"path": relative, "sha256": actual})
    manifest = json.loads((ROOT / "upstream/manifest.json").read_text())
    if manifest["commit"] != UPSTREAM_COMMIT:
        raise ValueError("Unexpected upstream commit")
    for record in manifest["files"]:
        source = ROOT / "upstream/atomic" / record["relative_path"]
        if digest(source) != record["sha256"]:
            raise ValueError(f"Pinned upstream file changed: {record['relative_path']}")
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=ROOT / "verification_output",
        help="A new directory for regenerated data and logs; it must not already exist.")
    args = parser.parse_args()
    output = args.output.expanduser().resolve()
    if output.exists():
        parser.error(f"Output directory already exists: {output}. Choose a new --output path.")

    try:
        import flint
    except ImportError as exc:
        raise SystemExit("Install the pinned dependency with: python3 -m pip install -r requirements.txt") from exc
    if flint.__version__ != "0.9.0":
        raise SystemExit(f"Expected python-flint 0.9.0; found {flint.__version__}")

    source_records = verify_sources()
    output.mkdir(parents=True)
    logs = output / "logs"
    logs.mkdir()
    constants = output / "lattice_constants.json"
    reference = output / "reference_coefficients.json"
    moments = output / "signed_moments.json"
    signs = output / "signs.json"

    steps = [
        ("upstream_finite_certificate",
         [sys.executable, "-I", "-B", "verification/check_certificate.py"],
         ROOT / "upstream/atomic", None),
        ("lattice_constants",
         [sys.executable, "-B", str(ROOT / "src/lattice_constants.py"),
          "--output", str(constants)], ROOT, constants),
        ("reference_coefficients",
         [sys.executable, "-B", str(ROOT / "src/reference_coefficients.py"),
          "--lattice-constants", str(constants), "--output", str(reference)], ROOT, reference),
        ("signed_moments",
         [sys.executable, "-B", str(ROOT / "src/verify_signed_moments.py"),
          "--reference", str(reference), "--output", str(moments)], ROOT, moments),
        ("global_sign_inequalities",
         [sys.executable, "-B", str(ROOT / "src/verify_signs.py"),
          "--lattice-constants", str(constants), "--reference", str(reference),
          "--output", str(signs)], ROOT, signs),
    ]

    summary = {
        "scope": "Finite inequalities and arithmetic only; analytic premises are stated in MATHEMATICAL_SCOPE.txt",
        "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(),
        "python_flint": flint.__version__,
        "flint": flint.__FLINT_VERSION__,
        "working_precision_bits": 160,
        "upstream_commit": UPSTREAM_COMMIT,
        "source_manifest_sha256": digest(ROOT / "SOURCE_SHA256SUMS.txt"),
        "source_files": source_records,
        "steps": [],
        "all_finite_checks_passed": False,
    }
    started = time.monotonic()
    for name, command, cwd, generated in steps:
        before = time.monotonic()
        completed = subprocess.run(command, cwd=cwd, text=True,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        elapsed = time.monotonic() - before
        log = logs / f"{name}.log"
        log.write_text(completed.stdout)
        record = {
            "step": name, "exit_code": completed.returncode,
            "elapsed_seconds": elapsed, "log": str(log.relative_to(output)),
            "log_sha256": digest(log),
        }
        if completed.returncode == 0 and generated is not None:
            record["output"] = generated.name
            record["output_sha256"] = digest(generated)
        summary["steps"].append(record)
        print(f"{name}: {'PASS' if completed.returncode == 0 else 'FAIL'} ({elapsed:.2f}s)", flush=True)
        if completed.returncode != 0:
            summary["elapsed_seconds"] = time.monotonic() - started
            (output / "verification_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
            print(completed.stdout, file=sys.stderr)
            return completed.returncode

    reference_data = json.loads(reference.read_text())
    signs_data = json.loads(signs.read_text())
    moments_data = json.loads(moments.read_text())
    assert reference_data["finite_checks_passed"]
    assert moments_data["finite_checks_passed"]
    assert signs_data["finite_checks_passed"]
    summary["all_finite_checks_passed"] = True
    summary["elapsed_seconds"] = time.monotonic() - started
    summary["results"] = {
        "absolute_mixture_mass": reference_data["absolute_mixture_mass"],
        "inherited_finite_error": reference_data["inherited_finite_error"],
        "inherited_tail_error": reference_data["inherited_tail_error"],
        "curvature_intervals": signs_data["curvature_intervals"],
        "other_taylor_intervals": signs_data["other_taylor_intervals"],
        "fourier_origin_derivative": signs_data["origin"]["Fourier_derivative_lower_bound"],
        "tail": signs_data["tail"],
        "normalized_polytope_fourth_coordinate": moments_data["integrated_pairs"][1]["normalized_exponential"],
    }
    (output / "verification_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"PASS: {signs_data['curvature_intervals']} curvature and "
          f"{signs_data['other_taylor_intervals']} other Taylor intervals; "
          "tail comparisons cover s >= 14.5.")
    print("Analytic premises: MATHEMATICAL_SCOPE.txt")
    print("Execution record:", output / "verification_summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
