# Interplanetary propagation and evaluation

The next search algorithm can call `assess_candidate()` to simulate one
proposal and inspect its cost and constraint residuals. No interplanetary
optimiser is implemented by this work.

## The execution path

1. `JplEphemerisProvider` supplies time-dependent barycentric planetary states.
2. `TabulatedEphemerisProvider` precomputes and checks a mission-length table.
3. `make_departure_state` adds position/velocity offsets to Earth's state at
   departure, in the same fixed coordinate axes.
4. `MissionDefinition` sets arrival and encounter requirements.
5. `MissionCandidate` proposes burns and an arrival time.
6. `simulate_adaptive_mission` propagates the spacecraft between exact burn times.
7. `evaluate_mission` measures arrival, closest approaches, order and clearance.

`assess_candidate` combines the final two steps and returns a status of
`feasible`, `infeasible`, or `unsafe`. An unsafe propagation has a failure
reason instead of a successful result. Configuration errors still raise.
No weighted objective or choice of search method has been imposed.

```python
from gravity_assist.evaluation import assess_candidate

assessment = assess_candidate(
    definition, candidate, ephemerides, bodies,
    max_step_s=86400.0,
)
if assessment.evaluation is not None:
    print(assessment.evaluation.total_delta_v_km_s)
    print(assessment.evaluation.constraint_residuals)
else:
    print(assessment.failure_reason)
```

## Requirements and measurements

Each intermediate body in `flyby_body_names` must have a corresponding
`FlybyRequirement`, in the same order. It specifies a time window and maximum
centre distance. The selected closest approach must lie inside the window,
before arrival, and after the previous selected encounter. Repeated visits
should use separate, non-overlapping windows to disambiguate encounters.
This first evaluator selects each window's closest approach; it does not
enumerate alternative sequences among several encounters within one window.

Arrival distance is measured against the destination at the candidate's exact
arrival time, not against its launch-time location or its nearest point at any
time. Arrival relative velocity is reported as a vector. The baseline evaluator checks
encounter distances; capture/insertion and a terminal velocity constraint are
not modelled in the propagated mission. Separate arrival-cost helpers are in
development, as described below.

Every numeric constraint residual is satisfied at <= 0. Positive values are
violations. Units are recorded in keys (`_km` or `_s`); do not sum these raw
numbers into an objective without choosing meaningful scales. Total delta-v
is the sum of burn magnitudes. Opposing burns both cost delta-v.

## Dynamics and time

States use Solar System barycentric position and ICRS axes. Mission time is
elapsed SI seconds from the timezone-aware departure epoch. Astropy handles
UTC/TAI/TDB conversion. Initial offsets and burn vectors use those same fixed
axes; they are not rotating local-orbital coordinates.

The default four-body force approximation includes the Sun, Earth, Jupiter
system and Saturn system. DE432s gives the Earth centre, and the giant-planet
system barycentres. The model treats each system as a point mass. Moon gravity
and the other planets are omitted. This approximation is for learning and
initial transfer experiments, not precision mission navigation.

GMs in `body_catalogue.py` come from the JPL DE431 parameter kernel:
https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de431.tpc
Jupiter and Saturn use BODY5/BODY6 system GMs, not BODY599/BODY699 planet GMs.
The Earth uses BODY399 and Sun BODY10. These are DE431-era parameters, not a
claim to reproduce every perturbation used to generate DE432s. Equivalent
mass is GM/G to preserve the existing `CelestialBody` interface.

Mean planetary radii are sourced from:
https://ssd.jpl.nasa.gov/planets/phys_par.html
Using those radii at giant-planet barycentres only approximates surface
clearance. Atmospheres, rings, oblateness, satellites, and radiation hazards
need additional constraints before close-flyby results can be called safe.

## Numerical choices

The original fixed-step RK4 is archived in `legacy/engine/integrators.py`.
Its mission wrappers and Jupiter experiments are under `legacy/` too; see
`legacy/README.md` for their module run commands. The current package does
not import the archived engine. The current mission path uses
SciPy DOP853 with adaptive error control and dense output. This is an ODE
integrator, not a route optimiser:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html

Each burn starts a new integration segment. Dense queries at burn times return
the post-burn state; velocity is never interpolated across the impulse.
The endpoint is the candidate's arrival time, including any final burn.

`MissionResult` defaults to the six-component `spacecraft` layout. The
explicit `planet_spacecraft` layout remains supported for archived wrappers,
which set it when producing the original twelve-component results.

Ephemeris Hermite interpolation uses both JPL positions and velocities.
It rejects extrapolation and checks every table interval's midpoint against
JPL. The error threshold applies to checked points; it is not a rigorous
all-times error bound. Rebuild at finer spacing as part of convergence checks.

Clearance event crossings can terminate propagation. The evaluator also
refines closest approaches by solving relative_position dot relative_velocity
= 0 on subdivisions of the adaptive mesh, and checks endpoints and burn
boundaries. This avoids relying on the smallest stored distance. Neither
event sign changes nor the finite refinement mesh constitutes a mathematical
guarantee against unresolved encounters. Verify convergence and restrict the
maximum step near actual flybys.

## Reproducing the demonstration

In PyCharm, run `experiments/interplanetary_baseline.py`. It reports progress,
prints both candidate results, saves reports under
`outputs/interplanetary_baseline/`, and opens two plots. The output folder is
ignored by Git. From the project root, the equivalent command is:

```text
.venv/bin/python -m experiments.interplanetary_baseline
```

Add `--no-show` for a non-interactive run or `--output-dir PATH` to choose
the artifact folder. Dependencies are listed in `requirements.txt`. The first
JPL provider construction downloads the kernel to Astropy's cache.

### Saved trajectory data

`trajectory_samples.npz` contains `times_s` with shape `(1200,)`, the
`no_burn` and `correction` spacecraft states with shape `(1200, 6)`, and
`sun_positions_km`, `earth_positions_km`, `jupiter_positions_km`, and
`saturn_positions_km`, each with shape `(1200, 3)`. Every row refers to the
same timestamp, using barycentric ICRS coordinates. The adjacent `summary.json`
records the departure epoch and frame. The desktop viewer loads these arrays
without querying the planetary ephemeris. For Sun-relative display, subtract the
saved Sun position from each body's or spacecraft's position at the same time.

### Desktop playback

After generating the baseline, install `ui/requirements.txt` in the project
virtual environment and run `.venv/bin/python ui/desktop.py`. The macOS
`ui/orbitForge.app` launcher runs the same entry point. See the
[viewer guide](../ui/README.md) for setup and controls.

The viewer displays Sun-relative positions in AU using ICRS axes. Its XY
projection is not an ecliptic-plane view. Planet markers are enlarged and
planet paths cover only the saved interval. Playback linearly interpolates
saved position samples; it does not rerun the integrator. Speed is interpolated
from saved barycentric velocity magnitudes. In particular, the viewer does not
resolve an impulsive velocity change exactly between neighbouring samples.
The saved full-flight evaluation is independent of the selected playback time.

### Explicitly illustrative conditions

- Departure: 2000-01-01 12:00:00 UTC; flight duration: eight Julian years.
- Start one million km from Earth in the direction of Earth's initial
  heliocentric velocity, with 11 km/s relative velocity in the same direction.
  This is a specified post-departure state, not a solved Earth escape.
- Jupiter: centre distance <= 2 million km, within years 0.25 to 6.
- Saturn: centre distance <= 1 million km at arrival, allowed in years 7 to 9.
- Minimum altitude: Sun 1 million km, Earth 200 km, Jupiter/Saturn 100000 km.
  These are example screening requirements, not validated mission limits.
- Trial correction: 50 m/s on day 180 along the fixed departure direction.
  Launch and Earth-departure delta-v/propellant are excluded from both costs.

### Measured results

| Candidate | Post-departure delta-v | Saturn distance at year 8 | Meets mission |
| --- | ---: | ---: | --- |
| No burn | 0 m/s | 3,117,391,866 km | No |
| Day-180 correction | 50 m/s | 3,125,742,400 km | No |

Both candidates miss Jupiter's encounter requirement as well. The minimum
distance within the configured Jupiter window is about 954.442 million km,
at its lower boundary: this is explicitly not a flyby. The closest points
over the full trajectory are reported separately from windowed measurements.

Tightening integration tolerances by ten, halving maximum step to 12 hours,
and halving ephemeris spacing to three hours changed the no-burn final
position by 0.135679 km, Saturn arrival distance by 0.071233 km, and the
windowed Jupiter distance by 0.000855 km. This is a combined numerical
refinement comparison, not an absolute physical-accuracy estimate or proof
of close-encounter accuracy. The full experiment took about 25 seconds on the
development machine; performance elsewhere will vary.

Tests also exercise an analytically known successful encounter sequence,
wrong order, a missed flyby, and a clearance violation between saved points.
Thus feasibility is tested on both passing and failing examples, even though
the real-ephemeris demonstration deliberately has no solved route.

## Boundary before the optimiser

The evaluation interface is ready for your search logic. The next work is to
choose search variables/bounds, a feasible initial guess or a feasibility
search, and scaled constraints. Departure time/state remain fixed within a
`MissionDefinition`; searching launch windows would create different mission
definitions and ephemeris tables. There is no guarantee that a small set of
correction burns can rescue the illustrative initial conditions above.


## Transfer and arrival helpers in development

`gravity_assist/transfer.py` contains small calculations used to prepare a
future targeted transfer:

- `heliocentric_state` subtracts the Sun's position and velocity from a body's
  state at the same time.
- `transfer_endpoints` obtains the departure body's state at time zero and the
  destination's state at the requested flight duration.
- `relative_speed` returns the magnitude of the difference between two velocity
  vectors expressed in the same frame.
- `hohmann_transfer_time` gives a half-ellipse flight-time reference for circular,
  coplanar orbits. It does not solve transfers between actual planetary states.

`gravity_assist/arrival.py` estimates periapsis speeds and a capture burn using
an isolated two-body approximation. `arrival_delta_v` selects `flyby` (zero
arrival burn) or `orbit_insertion` (the estimated braking cost). Distances are
from the body's centre, not altitudes above its surface. Capture assumes an
instantaneous burn at periapsis with aligned incoming and target velocities.

These helpers are not yet wired into the baseline evaluator or desktop UI.
A zero flyby burn is not a feasibility check, and an estimated capture cost is
not a propagated bound orbit. A Lambert solver and launch-window search remain
to be implemented. The provisional future-mission search range is departures
in 2030–2035 with flight durations of 5–10 years; these are proposed bounds,
not validated launch opportunities or a selected mission.
