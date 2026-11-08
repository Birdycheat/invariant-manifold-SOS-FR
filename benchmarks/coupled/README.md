# Two coupled centers

This analytic benchmark has 40 frozen rational instances, with transverse
dimension one or two. Independent exact leading-field checks give 34 stable
and six unstable systems. See the
[main guide](../../README.md#2-two-coupled-centers) for the equations, domain and
analytic enclosure.

| Source | Purpose |
|---|---|
| `run_coupled_batch.m` | Compare the numerical formulations, center the reduced programs and retain their outputs |
| `coupled_sos.m` | Build original/adapted-coordinate SOS programs and the four polynomial enclosure vertices |
| `write_center.m` | Export the central point with its ordered bases and matching Gram blocks |

The shared symbolic routines reconstruct the quartic graph correction and
classify the leading field. Common rational reconstruction, SOS identities and
witness verification are under [`code/shared/`](../../code/shared/); symbolic
jets and classes are under [`code/theory/`](../../code/theory/).

## Run from the repository root

```matlab
run_reproduction("solve", Benchmark="coupled", ...
    Instances="cp_ny1_s01.json", OutputDir=tempname);
```

`Enc=1` selects the reference enclosure; `Enc=2` selects the supplied
higher-order enclosure. The entry point loads shared helpers and executes new
searches in an external workspace.

```matlab
run_reproduction("symbolic");
run_reproduction("verify");
run_reproduction("witnesses");
```

These checks include this benchmark and cover all applicable retained families.
`Benchmark` selects new numerical searches or rounding; it does not restrict
the complete checks above to a single family.

## Inputs and evidence

[`data/coupled/`](../../data/coupled/) contains the frozen `instances/`,
floating `candidates/`, exact `certificates/`, and additional
`certificates_regions/`. The main cohort has 34 true-field certificates; the
jet-four and composite constructions bring the retained file count to 100.
[`results/coupled/`](../../results/coupled/) retains classifications and summaries.

New searches use the supplied exact instances. The `witnesses` action checks
retained points; a new solve does not automatically search for new witnesses.
