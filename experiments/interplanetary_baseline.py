"""Reproducible, deliberately unoptimised Earth-departure comparison.

Run as a module: python -m experiments.interplanetary_baseline --no-show
These are illustrative mission conditions, not a launch recommendation or
a solved Earth-Jupiter-Saturn route. Launch/departure cost is excluded.
"""

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import json
from pathlib import Path
from time import perf_counter

import matplotlib.pyplot as plt
import numpy as np

from gravity_assist.adaptive import simulate_adaptive_mission
from gravity_assist.body_catalogue import transfer_bodies
from gravity_assist.departure import make_departure_state
from gravity_assist.ephemerides import JplEphemerisProvider, TabulatedEphemerisProvider
from gravity_assist.evaluation import evaluate_mission
from gravity_assist.manoeuvres import ImpulsiveManoeuvre
from gravity_assist.mission_candidate import MissionCandidate
from gravity_assist.mission_definition import MissionDefinition, FlybyRequirement

DAY = 86400.0
YEAR = 365.25 * DAY
AU_KM = 149597870.7


def make_scenario(provider):
    earth = provider.get_state("Earth", 0.)
    sun = provider.get_state("Sun", 0.)
    # Outward-going post-departure state; offsets use fixed ICRS axes.
    direction = earth.velocity_km_s - sun.velocity_km_s
    direction /= np.linalg.norm(direction)
    initial = make_departure_state(provider, "Earth", 1e6*direction, 11.*direction)
    return MissionDefinition(
        departure_epoch=provider.departure_epoch,
        initial_spacecraft_state=initial,
        flyby_body_names=("Jupiter",), destination_body_name="Saturn",
        arrival_time_bounds_s=(7.*YEAR, 9.*YEAR),
        destination_max_distance_km=1e6,
        minimum_altitudes_km={"Sun": 1e6, "Earth": 200., "Jupiter": 1e5, "Saturn": 1e5},
        reference_frame=provider.frame,
        flyby_requirements=(FlybyRequirement("Jupiter", (.25*YEAR, 6.*YEAR), 2e6),),
    ), direction


def run_experiment(output_dir, show=True):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    started = perf_counter()
    raw = JplEphemerisProvider(datetime(2000, 1, 1, 12, tzinfo=timezone.utc))
    definition, direction = make_scenario(raw)
    bodies = transfer_bodies()
    print("Building and checking six-hour ephemeris table...", flush=True)
    table = TabulatedEphemerisProvider(raw, [b.name for b in bodies], 8.*YEAR)
    candidates = {
        "No burn": MissionCandidate((), 8.*YEAR),
        "Day-180 correction": MissionCandidate(
            (ImpulsiveManoeuvre(180.*DAY, .05*direction),), 8.*YEAR),
    }
    evaluations, trajectories = {}, {}
    for label, candidate in candidates.items():
        print(f"Simulating {label}...", flush=True)
        result = simulate_adaptive_mission(definition, candidate, table, bodies)
        print(f"Evaluating {label} ({len(result.times_s)} adaptive samples)...", flush=True)
        evaluation = evaluate_mission(definition, candidate, result, table, bodies)
        trajectories[label] = result
        evaluations[label] = asdict(evaluation)
        print(f"  Saturn arrival miss: {evaluation.arrival_distance_km:,.0f} km; "
              f"delta-v: {evaluation.total_delta_v_km_s*1000:.1f} m/s; "
              f"feasible: {evaluation.feasible}", flush=True)

    print("Checking stricter integration and three-hour ephemeris spacing...", flush=True)
    fine_table = TabulatedEphemerisProvider(raw, [b.name for b in bodies], 8.*YEAR, spacing_s=10800.)
    fine_result = simulate_adaptive_mission(
        definition, candidates["No burn"], fine_table, bodies,
        rtol=1e-11, position_atol_km=1e-5, velocity_atol_km_s=1e-11,
        max_step_s=43200.,
    )
    fine_evaluation = evaluate_mission(definition, candidates["No burn"], fine_result, fine_table, bodies)
    baseline = trajectories["No burn"]
    convergence = {
        "final_position_difference_km": float(np.linalg.norm(baseline.states[-1,:3]-fine_result.states[-1,:3])),
        "final_velocity_difference_km_s": float(np.linalg.norm(baseline.states[-1,3:]-fine_result.states[-1,3:])),
        "saturn_arrival_distance_difference_km": abs(evaluations["No burn"]["arrival_distance_km"]-fine_evaluation.arrival_distance_km),
        "jupiter_window_distance_difference_km": abs(evaluations["No burn"]["flybys"][0]["distance_km"]-fine_evaluation.flybys[0].distance_km),
    }
    for key, value in convergence.items():
        print(f"  {key}: {value:.9g}", flush=True)

    # Equally spaced samples for plots, independent of adaptive integrator steps.
    times = np.linspace(0., 8.*YEAR, 1200)
    planet_positions = {name: np.array([table.get_state(name,t).position_km for t in times])
                        for name in ["Sun", "Earth", "Jupiter", "Saturn"]}
    samples = {name: np.array([result.spacecraft_state_at(t) for t in times])
               for name, result in trajectories.items()}
    fig = plt.figure(figsize=(11, 8))
    axes = fig.add_subplot(projection="3d")
    sun = planet_positions["Sun"]
    for name in ["Earth", "Jupiter", "Saturn"]:
        xyz = (planet_positions[name]-sun)/AU_KM
        axes.plot(*xyz.T, alpha=.65, linewidth=1, label=name)
        axes.scatter(*xyz[-1], s=25)
    for name, states in samples.items():
        xyz = (states[:,:3]-sun)/AU_KM
        axes.plot(*xyz.T, label=name, linewidth=1.8)
        axes.scatter(*xyz[-1], marker="x", s=50)
    axes.scatter(0,0,0,color="gold",edgecolor="black",s=90,label="Sun")
    all_xyz = np.concatenate([(v-sun)/AU_KM for v in planet_positions.values()] +
                             [(v[:,:3]-sun)/AU_KM for v in samples.values()])
    centre = (all_xyz.max(axis=0)+all_xyz.min(axis=0))/2
    half = np.ptp(all_xyz,axis=0).max()/2*1.05
    axes.set_xlim(centre[0]-half, centre[0]+half)
    axes.set_ylim(centre[1]-half, centre[1]+half)
    axes.set_zlim(centre[2]-half, centre[2]+half)
    axes.set_box_aspect((1,1,1))
    axes.set_xlabel("Sun-relative ICRS x (AU)")
    axes.set_ylabel("Sun-relative ICRS y (AU)")
    axes.set_zlabel("Sun-relative ICRS z (AU)")
    axes.set_title("Eight-year Earth-departure trial: no optimised route\nMarkers show positions at year 8")
    axes.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(output_dir/"trajectory.png", dpi=160)

    fig2, axes2 = plt.subplots(2,1,figsize=(10,7),sharex=True)
    for ax, planet, limit in zip(axes2,["Jupiter","Saturn"],[2e6,1e6]):
        for name, states in samples.items():
            distance = np.linalg.norm(states[:,:3]-planet_positions[planet],axis=1)
            ax.plot(times/YEAR,distance/1e6,label=name)
        ax.axhline(limit/1e6,color="black",linestyle="--",label="Encounter distance limit")
        ax.set_ylabel(f"Distance to {planet}\n(million km)")
        ax.grid(alpha=.25)
        ax.legend(fontsize=8)
    axes2[-1].set_xlabel("Years after departure")
    fig2.suptitle("Distance to the moving planets — these candidates miss the mission targets")
    fig2.tight_layout()
    fig2.savefig(output_dir/"encounter_distances.png",dpi=160)
    np.savez_compressed(output_dir/"trajectory_samples.npz", times_s=times,
                        no_burn=samples["No burn"], correction=samples["Day-180 correction"])
    report = dict(
        departure_epoch=raw.departure_epoch.isoformat(), frame=raw.frame,
        flight_years=8., departure_offset_km=(1e6*direction).tolist(),
        departure_relative_velocity_km_s=(11.*direction).tolist(),
        departure_cost_included=False, ephemeris="DE432s, six-hour Hermite table",
        interpolation_errors=table.validation_errors, evaluations=evaluations,
        convergence=convergence, elapsed_seconds=perf_counter()-started,
        limitations=["Unoptimised illustrative initial conditions", "Four-body point-mass model",
                     "Jupiter/Saturn system barycentres, approximate surface clearances",
                     "No Saturn capture, ring-plane or radiation-belt model",
                     "Convergence differences are not absolute physical accuracy"],
    )
    (output_dir/"summary.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(f"Saved report and figures to {output_dir.resolve()}",flush=True)
    if show:
        plt.show()
    else:
        plt.close("all")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",default="outputs/interplanetary_baseline")
    parser.add_argument("--no-show",action="store_true")
    args = parser.parse_args()
    run_experiment(args.output_dir, show=not args.no_show)
