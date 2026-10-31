# Planar Lennard–Jones energy certificate

Reproducible verification code accompanying **The bulk energy of the planar
Lennard–Jones potential**, by Matthew J. Colbrook.

The certificate verifies the finite inequalities for the signed Gaussian
interpolation construction in the paper. It regenerates the lattice constants,
interpolation coefficients, signed moments, and interval bounds from the supplied
source files. The reference calculation passes all five stages, with **231
curvature intervals and 583 other Taylor intervals**, followed by the estimates
for the unbounded range.

## Quick start

Use Python 3.12. From the repository root:

```sh
python -m pip install -r requirements.txt
python run_verification.py
```

For an isolated installation, first create a virtual environment with
`python -m venv .venv` and activate it using `.venv\Scripts\Activate.ps1`
in PowerShell or `source .venv/bin/activate` in a Unix shell.

The tested environment is Python 3.12.14, python-flint 0.9.0, and FLINT 3.6.0.
Ball arithmetic uses 160-bit working precision. The lattice-sum calculation and
the upstream finite checker use only Python's standard library. After installing
the dependency, verification runs entirely offline.

The driver creates `verification_output/`. For another run, choose a new output
directory:

```sh
python run_verification.py --output verification_output_second
```

Keep Python assertions enabled: do not use `-O` or `-OO`. The driver verifies the
source hashes before starting, rejects an existing output directory, and stops at
the first failed stage. Run time is approximately 20–30 seconds in the tested
environment and depends on the machine.

## Verification stages

| Stage | Calculation |
| --- | --- |
| 1 | Execute the unchanged upstream exact finite checker. |
| 2 | Enclose the triangular lattice sums, optimal spacing, and energy using rational arithmetic and Euler–Maclaurin summation. |
| 3 | Solve the two 20-dimensional systems for the integrated finite approximation; bound the signed mass and coefficient errors. |
| 4 | Recompute the signed moments and verify that the integrated data lie outside the source's Gaussian parameter polytope. |
| 5 | Prove the origin estimates and compact sign inequalities, check the complete interval covering, and verify the numerical tail comparisons. |

A successful run prints `PASS` for every stage and sets
`all_finite_checks_passed` to `true` in `verification_summary.json`. Separate JSON
files contain the numerical enclosures and every certified rational interval;
`logs/` contains each stage's output. The summary includes the software versions,
precision, source hashes, output hashes, and elapsed times.

The driver regenerates all numerical inputs. The bundled [`reference/`](reference/)
directory is a completed execution for comparison and is not used as input.

## Repository contents

```text
.
├── README.md                    Instructions and overview
├── CITATION.cff                 Citation metadata
├── MATHEMATICAL_SCOPE.txt       Analytic premises and verified assertions
├── requirements.txt            Pinned Python dependency
├── run_verification.py         Complete verification driver
├── SOURCE_SHA256SUMS.txt        Source and documentation hashes
├── src/                        Four verification programs
├── reference/                  Verified outputs and execution logs
└── upstream/                   Pinned interpolation source and finite checker
    ├── manifest.json           Provenance, source URLs, and hashes
    ├── LICENSE                 Upstream Apache 2.0 license
    └── atomic/                 Unmodified upstream files
```

Generated runs, virtual environments, and Python caches are excluded by
`.gitignore`. The `.gitattributes` file preserves the exact bytes of the pinned
upstream files and reference outputs; the remaining text files use LF line
endings so their source hashes remain stable across platforms.

## Mathematical scope and provenance

The finite approximation, interpolation basis, damping parameters, correction
estimates, and tail estimates come from OpenAI's *An atomic certificate for
triangular-lattice universal optimality* (26 September 2026), pinned to revision
`fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`. The included
[`upstream/manifest.json`](upstream/manifest.json) identifies each redistributed
file and its public source. The driver checks these files against their pinned
hashes. The [upstream license](upstream/LICENSE) applies to those materials.

The new calculation verifies the inequalities for the signed Lennard–Jones
integral. The analytic interpolation, convergence, Fourier comparison, and
thermodynamic-limit arguments are supplied by the paper and its cited sources.
[`MATHEMATICAL_SCOPE.txt`](MATHEMATICAL_SCOPE.txt) gives the detailed correspondence
between those arguments and the finite arithmetic assertions checked by the
programs.

## Citation and contact

Please cite the accompanying manuscript:

> Matthew J. Colbrook. *The bulk energy of the planar Lennard–Jones potential*.
> Manuscript, 2026.

Machine-readable metadata are provided in [`CITATION.cff`](CITATION.cff).

Matthew J. Colbrook  
Department of Applied Mathematics and Theoretical Physics, University of Cambridge  
[m.colbrook@damtp.cam.ac.uk](mailto:m.colbrook@damtp.cam.ac.uk)
