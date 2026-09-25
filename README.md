# orbitForge

A Python project for studying interplanetary transfers, gravity assists and
spacecraft manoeuvres, with a desktop viewer for saved simulation results.

The current simulator propagates and evaluates proposed trajectories using JPL
planetary data. The viewer compares two baseline routes and animates the
spacecraft and planets. Neither baseline meets the Jupiter/Saturn encounter
requirements; they demonstrate the simulation pipeline, not a solved mission.

General transfer and arrival-cost helpers are in development. A Lambert solver,
launch-window search and a simulated Saturn capture mission are not yet implemented.

## Getting started

Run these commands from the repository root. The desktop launcher below is for
macOS; the Python simulation can also be run independently of the UI.

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m pip install -r ui/requirements.txt
```

If you already have a project virtual environment, use it instead of creating
another one. The UI dependency is optional when running simulations alone.

### Generate the baseline results

```sh
.venv/bin/python -m experiments.interplanetary_baseline --no-show
```

The first run downloads the JPL ephemeris kernel; later runs use its local
cache. Results are written to `outputs/interplanetary_baseline/`:

- `trajectory_samples.npz`: spacecraft states and aligned planet positions.
- `summary.json`: mission evaluations, frame information and convergence results.
- `trajectory.png` and `encounter_distances.png`: static plots.

Generated outputs are ignored by Git, so a fresh checkout needs to run the
experiment before opening the viewer. Omit `--no-show` to display the plots.
In PyCharm, use the project interpreter and run the experiment as a module.

### Open the desktop viewer

```sh
.venv/bin/python ui/desktop.py
```

On macOS, you can also double-click `ui/orbitForge.app` in Finder. Keep the
launcher in that folder: it uses this repository's `.venv` and saved outputs.
It is not a self-contained distributable application.

The viewer supports trajectory comparison, 3D rotation, an XY projection and
zoom. Use the time slider or playback controls to inspect the spacecraft state,
or jump to departure, the day-180 correction and the arrival epoch. Full-flight
results remain separate from measurements at the selected time.

The desktop app reads local files and needs no browser tab or web server.
It does not run simulations or modify mission parameters. Restart it after
regenerating outputs. See the [viewer guide](ui/README.md) for controls and
custom data folders.

## Simulation workflow

```text
MissionDefinition + MissionCandidate + planetary ephemerides
                         ↓
               spacecraft propagation
                         ↓
                   MissionResult
                         ↓
               mission evaluation
                         ↓
            saved results → desktop viewer
```

`gravity_assist.evaluation.assess_candidate()` propagates one proposal and
reports `feasible`, `infeasible` or `unsafe`. It does not search for a better
proposal. The baseline experiment also compares results at tighter numerical
settings to check convergence.

## Project layout

| Path | Purpose |
| --- | --- |
| `experiments/interplanetary_baseline.py` | Main runnable simulation and output generation |
| `gravity_assist/ephemerides.py` | JPL planetary states and checked interpolation |
| `gravity_assist/mission_definition.py` | Mission and encounter requirements |
| `gravity_assist/mission_candidate.py` | Proposed burns and arrival time |
| `gravity_assist/adaptive.py` | Adaptive spacecraft integration |
| `gravity_assist/interplanetary.py` | Gravity and clearance checks |
| `gravity_assist/evaluation.py` | Arrival and encounter measurements |
| `gravity_assist/transfer.py` | Transfer preparation and Hohmann reference helpers, in development |
| `gravity_assist/arrival.py` | Idealised flyby/insertion cost helpers, in development |
| `ui/` | Desktop viewer, local app launcher and viewer data checks |
| `tests/` | Simulation and shared-model tests |
| `legacy/` | Earlier fixed-step RK4 engine and Jupiter experiments |
| `docs/interplanetary_model.md` | Model assumptions and baseline results |

The current simulation package does not depend on the archived engine.
See the [legacy guide](legacy/README.md) to revisit those experiments.

## Model and limits

The baseline uses Sun, Earth, Jupiter-system and Saturn-system gravity,
six-component spacecraft states, adaptive DOP853 integration and impulsive
burns. Public units are km, seconds and km/s. Simulation states use a Solar
System barycentric origin and ICRS axes; the viewer subtracts the Sun's saved
position and displays positions in AU. Its speed readout remains barycentric.

The baseline begins after Earth departure, excluding launch and departure
costs. Giant-planet system barycentres and point masses are approximations.
Rings, moon encounters and realistic orbit-insertion constraints are not
included. The new arrival helpers estimate an idealised capture burn; they do
not yet add capture to the propagated mission or establish a safe orbit.

See the [model guide](docs/interplanetary_model.md) for the detailed assumptions,
measured baseline results and interpretation of the convergence checks.

## Checks

Run the simulation tests:

```sh
.venv/bin/python -m pytest -q tests
```

To omit the archived regression cases, add `--ignore=tests/legacy`.
Run the viewer's data-conversion checks separately:

```sh
.venv/bin/python -m unittest discover -s ui -p 'test_*.py' -v
```

## Next development

The next step is a general targeted transfer calculation, with Earth to Saturn
as an initial case. This will connect future planetary positions, estimate
Earth departure demand and compare flyby versus orbit-insertion arrival costs.
Launch-date and flight-duration searches will follow, then physically compatible
Jupiter gravity assists and validation with the multi-body simulator.

Monte Carlo uncertainty and stochastic disturbances remain later work.
