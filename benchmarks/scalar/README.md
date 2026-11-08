# Scalar polynomial benchmark

This folder contains the scalar-center programs and their exact scalar checks.
The main benchmark has 80 stable instances: transverse dimension `N=1,3`,
reduced order `m=5,7`, and 20 seeds per group. The additional 40 order-three
instances remain available for the full checks. See the
[main guide](../../README.md#1-scalar-polynomial-systems) for the system and domains.

## Sources to read

| Sources | Purpose |
|---|---|
| `run_batch_exact.m` | Compare the numerical formulations and retain candidate records |
| `run_jet_exact.m`, `jet_sos_center.m` | Center the prescribed jet-reduced program and export its Gram point |
| `gen_instance.m`, `gen_instance_exact.m` | Sample the scalar family, then freeze rational coefficients and reconstruct exact criticality |
| `exact_jet.m`, `exact_data.m` | Reconstruct the invariant jet from the field and select exact system coefficients |
| `gen_cm_cert.m`, `cm_verify_exact.m` | Construct an anchored certificate and check its univariate matrix positivity exactly |
| `gen_sos_cert.m`, `gen_sos_cert_z.m`, `jp_cert.m` | Original-coordinate, jet-coordinate and partially quadratic formulations |
| `classify_candidate.m` | Propose numerical invalidity points and check derivative signs exactly |
| `cm_verify_anchors.m`, `cm_verify_defects.m`, `cm_quadratic_obstruction.m` | Check retained anchors, floating-candidate transverse defects and necessary quadratic conditions |
| `run_scalar_regions.m` | Search enlarged strips and attempt exact reconstruction |

The shared reconstruction and SOS verification routines are in
[`code/shared/`](../../code/shared/); cross-benchmark symbolic checks are in
[`code/theory/`](../../code/theory/).

## Run from the repository root

```matlab
run_reproduction("solve", Benchmark="scalar", N=1, M=5, Seeds=1, ...
    OutputDir=tempname);
run_reproduction("anchors");
run_reproduction("defects");
```

The entry point loads dependencies and prepares an external workspace for new
searches. `anchors` and `defects` check the complete retained scalar cohorts.
The `symbolic`, `verify` and `witnesses` actions also include this benchmark;
those three actions cover all applicable retained families.

## Inputs and evidence

- [`data/benchmark/`](../../data/benchmark/) holds the exact field coefficients,
  anchored records and floating candidates.
- [`data/jet_inputs/`](../../data/jet_inputs/) holds the retained central points.
- [`data/certificates/`](../../data/certificates/) and
  [`data/certificates_regions/`](../../data/certificates_regions/) hold the
  rational certificates on the original and enlarged strips.

## Complementary controls

`controls/` contains four relevant sources, rather than another main benchmark:

- `run_example1.m` computes the retained scalar example used by the
  kernel-ratio workflow; `make_fig_data.m` regenerates numerical Gram-kernel
  diagnostics from that example and the full scalar groups.
- `review_dd_scaling.m` runs the small DD-scaling experiment.
  `cm_dd_evidence.m` checks the retained raw records, applied scalings and exact
  witness signs; the `tables` action includes that evidence check.

```matlab
run_reproduction("solve", Benchmark="scalar-strips", OutputDir=tempname);
run_reproduction("solve", Benchmark="dd-scaling", OutputDir=tempname);
run_reproduction("solve", Benchmark="kernel-ratios", OutputDir=tempname);
```

Kernel-ratio regeneration requires all four main scalar groups and 20 seeds.
The historical source under `data/review_dd_scaling/` is authenticated evidence;
the current DD producer is `controls/review_dd_scaling.m`.
