# orbitForge

> **Status: Work in progress.** A nominal deterministic Jupiter flyby can now
> be propagated, measured, and visualised. The engine also supports physical
> collision checks and optimiser-ready impulsive manoeuvres scheduled at
> arbitrary times. Ephemeris-driven interplanetary propagation and encounter
> evaluation are now available; the interplanetary optimiser is not implemented.
> Seventy-three automated tests currently protect the foundation.

## Project goal

This project aims to design fuel-efficient spacecraft routes to different
planets using planetary gravity assists and timed thruster manoeuvres.

The deterministic part of the project will propagate candidate trajectories,
apply manoeuvres, and optimise quantities such as encounter timing, flyby
geometry, and burn vectors. The principal objective will be to reach a chosen
destination while minimising total velocity change and propellant use.

Monte Carlo simulation will then test how reliably an optimised route succeeds
when its initial state, navigation solution, or thruster execution is
imperfect. A later stochastic-calculus extension will investigate continuous
random disturbances during flight.

## Research questions

> How can gravity assists and spacecraft manoeuvres be combined to reach a
> target planet using minimal propellant?

> How robust are the resulting trajectories to uncertain initial conditions,
> manoeuvre errors, and continuous stochastic disturbances?

## Why this matters

Gravity assists can produce large changes in heliocentric velocity without
requiring the spacecraft to supply the equivalent change using fuel. Their
outcomes are nevertheless sensitive to the incoming position and velocity.
Small errors can alter the flyby altitude, turning angle, outgoing direction,
and eventual planetary encounter.

The project will therefore connect three different numerical tasks:

1. Propagate and optimise deterministic interplanetary trajectories.
2. Evaluate optimised trajectories with Monte Carlo simulation.
3. Extend the dynamics with a physically justified stochastic process.

Monte Carlo sampling is used to measure robustness, not as a replacement for
the deterministic optimiser that searches for fuel-efficient routes.

## Mathematical model

### Deterministic dynamics

The initial flyby model will use the Sun, Jupiter, and a spacecraft of
negligible mass. The Sun and Jupiter affect the spacecraft, while the
spacecraft does not affect either massive body.

The deterministic equations are

```text
dr/dt = v
dv/dt = a(r, t),
```

where `a(r, t)` is the combined gravitational acceleration. These ordinary
differential equations will be integrated using the fourth-order Runge-Kutta
method (RK4).

### Thruster manoeuvres and fuel

The first thruster model will treat a short burn as an instantaneous velocity
change:

```text
v_after = v_before + delta_v.
```

A candidate mission may contain multiple burns separated by coast and flyby
segments. The initial optimisation objective will minimise total `delta_v`;
the rocket equation will later convert this into propellant mass. Arrival,
flight-time, collision-avoidance, and flyby-altitude requirements will be
treated as constraints.

### Monte Carlo uncertainty

Initial position, initial velocity, navigation, and manoeuvre-execution errors
will be sampled from stated probability distributions. Each sample will then
follow the deterministic mission model. The resulting distribution of
destinations, fuel costs, and failed encounters will quantify the route's
robustness, but this is not by itself stochastic calculus.

### Continuous stochastic disturbances

Continuous unmodelled acceleration will later be represented by the stochastic
differential equation

```text
dr_t = v_t dt
dv_t = a(r_t, t) dt + B dW_t,
```

where `W_t` is a three-dimensional Wiener process and `B` controls the scale
and direction of the disturbance. This model will be simulated with an SDE
method such as Euler-Maruyama rather than RK4.

The noise must be physically interpreted and calibrated; it will not be added
solely to make the model stochastic.

## Units

- Distance: kilometres
- Time: seconds
- Velocity: kilometres per second
- Mass: kilograms
- Acceleration: kilometres per second squared

## Current implementation

See [the interplanetary model guide](docs/interplanetary_model.md) for the
departure builder, explicit flyby requirements, adaptive propagator, evaluator,
and a reproducible eight-year Earth-departure experiment with results.
Run `experiments/interplanetary_baseline.py` in PyCharm to generate its plots.
The two trial candidates are explicitly infeasible; no optimised route is claimed.

### Ephemeris-driven mission propagation

`gravity_assist.interplanetary.simulate_mission` now accepts a
`MissionDefinition`, `MissionCandidate`, an `EphemerisProvider`, a sequence of
`CelestialBody` objects, and `dt_s`. The caller supplies masses and radii for
the bodies included in the force model. Those bodies must include the Sun,
destination, required flyby bodies, and bodies with specified clearance limits.

```python
from gravity_assist.interplanetary import simulate_mission

result = simulate_mission(
    definition=definition,
    candidate=candidate,
    ephemerides=ephemerides,
    bodies=bodies,
    dt_s=dt_s,
)
spacecraft_positions_km = result.spacecraft_states[:, :3]
spacecraft_velocities_km_s = result.spacecraft_states[:, 3:]
```

This example assumes the input objects have already been constructed.
The spacecraft state and burn vectors must use the provider's coordinate axes.
The default mission frame matches `JplEphemerisProvider`: ICRS axes with a
Solar System barycentric origin. Departure epochs and declared frames are
checked for agreement; these checks do not transform incorrectly supplied
coordinates. Do not reuse the old heliocentric initial state directly.

The new derivative queries planetary positions at each RK4 stage and sums
their Newtonian gravitational accelerations. Planetary trajectories are
prescribed by the ephemeris; only the spacecraft is integrated.
Scheduled burns work at arbitrary times, including departure and arrival.
A sample at a burn time contains the post-burn velocity.

`MissionResult.state_layout` explicitly distinguishes `spacecraft` (six
columns) from `planet_spacecraft` (the existing twelve-column default).
`spacecraft_states` provides the six spacecraft columns for either layout.
Existing Jupiter experiments retain their original layout and dynamics.

Clearance checks raise `UnsafeTrajectoryError` at evaluated RK4 stages and
stored samples. They are not continuous collision detection: sufficiently
coarse steps can skip an encounter. Fixed-step convergence must be established
for each mission, especially near flybys. DE432s Jupiter/Saturn positions are
planetary-system barycentres; their surface clearances are approximate.

This entry point **propagates a proposal**. For refined encounter evaluation,
use `simulate_adaptive_mission` followed by `evaluate_mission`, or the combined
`assess_candidate` entry point. Interplanetary optimisation remains to be implemented.
Supplying a required flyby name includes that body in the model; it does not
automatically steer the spacecraft toward it. Earth departure and launch costs
remain outside the initial post-departure model.

Verification includes synthetic moving-body gravity, circular-orbit timestep
convergence, scheduled burns and frame/epoch compatibility. A one-day JPL
smoke test with Sun/Earth/Jupiter/Saturn gravity and an off-grid burn completed
at 3600 s and 1800 s steps, with an approximately 1.6e-6 km difference in final
position. This short test is not validation of a multiyear transfer.

### Existing Jupiter demonstration

The repository currently provides a working deterministic flyby model:

- `CelestialBody` and `OrbitalState` data models;
- gravitational acceleration from multiple bodies;
- a reusable RK4 step and fixed-step trajectory propagator;
- packing of the Jupiter and spacecraft states into a 12-component vector;
- a restricted three-body derivative for a fixed Sun, moving Jupiter, and
  massless spacecraft;
- a 50-day nominal Jupiter encounter propagated at 60-second intervals;
- calculation of closest-approach time, centre distance, and surface altitude;
- comparison of incoming and outgoing velocities at equal Jupiter distances;
- calculation of turning angle and heliocentric speed change;
- comparison of the simulated turning angle with two-body hyperbolic theory;
- physical detection of spacecraft intersections with massive bodies;
- validated impulsive-manoeuvre objects with a calculated `delta_v` magnitude;
- application of manoeuvres without mutating the original spacecraft state;
- mission propagation through ordered burns at arbitrary, off-grid times;
- a controlled Jupiter experiment comparing a nominal trajectory with a
  10 m/s correction manoeuvre on mission day 5;
- a three-dimensional close-up containing a correctly scaled Jupiter and the
  nearby spacecraft trajectory;
- automated tests covering gravity, collisions, manoeuvres, mission propagation,
  ephemerides, encounter measurements and feasibility checks.

Deterministic trajectory optimisation, refined asymptotic-state estimation,
Monte Carlo experiments, and SDE integration have not yet been implemented.

## Timed-manoeuvre experiment

The first controlled manoeuvre experiment applies a 10 m/s velocity change in
the positive y-direction five days after the start of the nominal Jupiter
encounter. Both cases use the same initial state, force model, 60-second step,
and 50-day duration.

| Metric | No burn | Day-5 burn | Change |
| --- | ---: | ---: | ---: |
| Closest-approach altitude | 371,210.4 km | 376,455.0 km | +5,244.7 km |
| Turning angle | 130.840° | 130.490° | -0.350° |
| Heliocentric speed change | +2.973557 km/s | +2.931846 km/s | -0.041711 km/s |

The positive y-burn increases the flyby altitude, so Jupiter bends the
spacecraft's path less and produces a smaller heliocentric speed gain. This is
an illustrative sensitivity experiment, not an optimised manoeuvre. Its
purpose is to verify the complete path from a scheduled burn to a measurable
change in the encounter.

## Current project layout

```text
.
├── experiments/
│   └── deterministic_flyby.py
├── gravity_assist/
│   ├── __init__.py
│   ├── collisions.py
│   ├── constants.py
│   ├── forces.py
│   ├── integrators.py
│   ├── manoeuvres.py
│   ├── mission.py
│   ├── models.py
│   └── simulation.py
├── tests/
│   ├── test_collisions.py
│   ├── test_forces.py
│   ├── test_manoeuvres.py
│   ├── test_mission.py
│   └── test_simulation.py
├── .gitignore
└── README.md
```

## Development roadmap

### Phase 1: validate the single-flyby foundation

- Test the celestial-body and orbital-state validation.
- **Completed:** check gravitational acceleration against a hand-calculated
  case.
- Test the general RK4 step using a differential equation with a known result.
- Perform timestep-convergence checks.
- **Completed:** detect body intersections and stop invalid propagation.
- **Completed:** compare the simulated turning angle with two-body hyperbolic
  theory.

### Phase 2: generalise the mission model

- Represent the time-dependent states of multiple planets.
- Introduce documented planetary ephemerides and reference frames.
- Define departure conditions and target-planet arrival requirements.
- Support propagation through multiple coast and flyby segments.

### Phase 3: introduce thruster manoeuvres and fuel

- **Completed:** represent three-dimensional impulsive manoeuvres and schedule
  them at arbitrary times.
- **Completed:** propagate a mission through an ordered list of burns without
  restricting burn times to the RK4 step grid.
- **Completed:** apply a trajectory-correction manoeuvre to the nominal Jupiter
  encounter and measure the resulting change.
- Accumulate the mission's total `delta_v`.
- Convert `delta_v` into propellant mass using the rocket equation.

### Phase 4: optimise deterministic routes

- Begin with a fixed planet sequence and destination.
- Optimise burn timing, burn vectors, and flyby geometry.
- Minimise fuel subject to arrival, duration, altitude, and safety constraints.
- Compare alternative gravity-assist sequences and target planets.

### Phase 5: propagate initial uncertainty with Monte Carlo

- Choose and justify distributions for navigation and manoeuvre errors.
- Propagate many perturbed versions of an optimised route.
- Measure arrival dispersion, fuel use, unsafe flybys, and mission failures.
- Compare expected performance with conservative fuel-reserve requirements.

### Phase 6: add stochastic calculus

- Define a physically meaningful continuous acceleration-noise model.
- Implement Euler-Maruyama and verify its `sqrt(dt)` noise scaling.
- Compare SDE paths with deterministic and Monte Carlo uncertainty models.
- Perform time-step and sample-size convergence checks.

### Phase 7: produce the research tool and analysis

- Present optimised routes to selected destination planets.
- Compare direct transfers with routes using gravity assists.
- Report fuel, flight-time, and robustness trade-offs.
- Document model limitations and sensitivity to uncertainty assumptions.

## Validation principles

Results will not be treated as meaningful until the model passes appropriate
checks. These include deterministic time-step convergence, conservation tests,
Monte Carlo sample-size convergence, and SDE time-step convergence. Any noise
parameters and probability distributions will be stated explicitly.

## Intended final output

The finished project will be a research-scale interplanetary trajectory-design
desktop application. A user will choose a destination or candidate planet
sequence, and the system will search for a low-fuel combination of gravity
assists and thruster manoeuvres. It will animate the resulting mission and
report its total `delta_v`, estimated propellant use, flight time, closest
approaches, and arrival accuracy.

The planned interface will be a native desktop window built separately from
the numerical engine. Mission controls will supply inputs to the tested Python
modules, while interactive three-dimensional views and result panels will
display the output. This separation keeps the physics usable without the user
interface and allows a future optimiser to call the same mission functions.

Monte Carlo and SDE experiments will then show how the nominal solution changes
under uncertainty. The result is intended as an educational and research
prototype, not an operational mission-planning system. Limitations in the
dynamics, ephemerides, manoeuvre model, optimiser, and uncertainty assumptions
will be documented explicitly.

## Running the current experiment

The current code requires Python 3.10 or later, NumPy, Matplotlib, and pytest.
From the project root, install the dependencies if necessary and run:

```bash
python -m pip install numpy matplotlib pytest
python -m pytest -q
python -m experiments.deterministic_flyby
```

The experiment propagates the nominal and day-5 manoeuvred trajectories over
50 days. It reports their closest approaches, turning angles, Jupiter-relative
outgoing speeds, and heliocentric speed changes. It also opens a correctly
scaled three-dimensional close-up of the nominal encounter.
