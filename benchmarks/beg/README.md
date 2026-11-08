# BEG mean-field model

The Blume-Emery-Griffiths benchmark compares the degree-five surrogate and the
true analytic field at 15 rational parameter values: 11 stable and four unstable.
See the [main guide](../../README.md#4-beg-analytic-model) for the equations,
tricritical parameterization and certified radii.

| Source | Purpose |
|---|---|
| `run_beg_batch.m` | Build the surrogate and true-field numerical programs, center their prescribed jet faces and export central points |
| `cm_symbolic_beg.m` | Derive the invariant graph and reduced quintic term symbolically, then check the retained parameter values |
| `cm_beg_bounds.m` | Reconstruct Taylor coefficients, rational complex-polydisc bounds and analytic remainders, and check the physical domain |

The shared model, rational reconstruction and certificate/witness verification
are in [`code/shared/`](../../code/shared/). The same independently reconstructed
analytic bounds support both symbolic checks and true-field verification.

## Run from the repository root

```matlab
run_reproduction("solve", Benchmark="beg", ...
    Instances="a_2_1.json", OutputDir=tempname);
```

`a_2_1` denotes `a=2/1`. The entry point loads shared helpers and executes new
searches in an external workspace. To check the retained reference results:

```matlab
run_reproduction("symbolic");
run_reproduction("verify");
run_reproduction("witnesses");
```

These actions include BEG and cover all applicable retained families.
`Benchmark` selects new numerical searches or rounding; it does not restrict
the complete checks above to this model.

## Inputs and evidence

[`data/beg/`](../../data/beg/) contains the exact `instances/`, floating
`candidates/` and rational `certificates/`.
[`results/beg/`](../../results/beg/) holds classifications and centering summaries.
The 40 retained certificate files include surrogate and true-field proofs at
several radii; the main cohort counts 11 surrogate and 11 true-field problems.
The symbolic checks reconstruct all 45 stored analytic bounds independently.

`witnesses` checks the retained surrogate or true-field derivative signs at
stored points. A new numerical solve does not automatically search for new
invalidity witnesses.
