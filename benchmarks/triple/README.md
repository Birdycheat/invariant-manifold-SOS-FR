# Three polynomial centers

This benchmark has 40 frozen rational instances with three central variables
and transverse dimension one or two. Independent weighted AM-GM checks give
32 stable and eight unstable systems. See the
[main guide](../../README.md#3-three-centers) for the equations and domain.

`run_triple_batch.m` contains the numerical driver and its local construction,
centering and export helpers. It compares the formulations, including the
original-coordinate layerwise reduction for transverse dimension one, and
exports the bases actually used with their matching Gram matrices.

The exact quartic graph correction, reduced leading field and system classes
are reconstructed by [`code/theory/`](../../code/theory/). The common model,
rational reconstruction, certificate verification and retained-witness checks
are in [`code/shared/`](../../code/shared/).

## Run from the repository root

```matlab
run_reproduction("solve", Benchmark="triple", ...
    Instances="tr_ny1_s01.json", OutputDir=tempname);
```

The entry point loads shared helpers and prepares a fresh external workspace.
To check the retained reference results:

```matlab
run_reproduction("symbolic");
run_reproduction("verify");
run_reproduction("witnesses");
```

These actions include this benchmark and cover all applicable retained families.
`Benchmark` selects new numerical searches or rounding; it does not restrict
the complete checks above to a single family.

## Inputs and evidence

[`data/triple/`](../../data/triple/) contains the frozen `instances/`, floating
`candidates/` and 32 exact `certificates/`.
[`results/triple/`](../../results/triple/) retains classifications and rounding
summaries. The original-coordinate transverse-dimension-one programs are a
separate representation, rather than another certificate cohort.

New searches consume these exact instance files rather than replaying the
former random-generator stream. `witnesses` checks retained points; a new solve
does not automatically search for new witnesses.
