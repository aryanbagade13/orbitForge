# Earlier orbitForge work

This directory preserves the original Jupiter experiments and fixed-step
engine for reference. The active experiment is
`experiments/interplanetary_baseline.py` in the project root.

## Experiments

| File | Purpose |
| --- | --- |
| `experiments/deterministic_flyby.py` | Original Jupiter close-up and day-5 burn comparison |
| `experiments/jupiter_flyby_scenario.py` | Reusable original Jupiter setup and measurements |
| `experiments/manoeuvre_sensitivity.py` | Day-5 y-burn sweep from -50 to +50 m/s |
| `experiments/optimise_jupiter_flyby.py` | Early burn-time/y-burn candidate grid, not a completed optimiser |

Run these from the project root with the project interpreter:

```text
python -m legacy.experiments.deterministic_flyby
python -m legacy.experiments.manoeuvre_sensitivity
python -m legacy.experiments.optimise_jupiter_flyby
```

In PyCharm, an old run configuration may still point at its previous file
location. Update it to the new module path above, with the project root as
the working directory. Each plotting experiment opens its figures as before.

## Engine

`engine/integrators.py` contains your original RK4 and fixed-step propagator.
`engine/simulation.py` contains the 12-component Jupiter/spacecraft dynamics.
`engine/mission.py` applies scheduled burns with fixed-step propagation.
`engine/interplanetary.py` retains the early fixed-step ephemeris wrapper.

Shared gravity, body and burn models stay in the current `gravity_assist`
package. Legacy wrappers explicitly declare their result layout; the main
`MissionResult` now defaults to a six-component spacecraft state.

Regression tests live in `tests/legacy/` and remain part of the complete test
suite. Nothing here is automatically imported by the current mission engine.
The earlier code remains useful for understanding and checking the physics.
