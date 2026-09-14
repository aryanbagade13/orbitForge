from gravity_assist.manoeuvres import ImpulsiveManoeuvre

import numpy as np

from experiments.jupiter_flyby_scenario import SECONDS_PER_DAY

burn_times_days = [3.0, 5.0, 7.0]
y_burns_m_s = [-20.0, 0.0, 20.0]

candidate_results = []

from experiments.jupiter_flyby_scenario import (
    propagate_jupiter_flyby,
    calculate_flyby_comparison_metrics,
)



def evaluate_candidate(manoeuvre: ImpulsiveManoeuvre):
    result = propagate_jupiter_flyby([manoeuvre])
    metrics = calculate_flyby_comparison_metrics(
        result.times_s,
        result.states
    )

    metrics["total_delta_v_km_s"] = manoeuvre.magnitude_km_s

    return metrics

def run_candidate_grid(burn_times_days, y_burns_m_s):
    candidate_results = []

    for burn_time in burn_times_days:
        for y_burn_m_s in y_burns_m_s:
            manoeuvre = ImpulsiveManoeuvre(
                time_s=burn_time * SECONDS_PER_DAY,
                delta_velocity_km_s=np.array([
                    0.0,
                    y_burn_m_s / 1000.0,
                    0.0,
                ]),
            )

            metrics = evaluate_candidate(manoeuvre)

            candidate_results.append({
                "manoeuvre": manoeuvre,
                "metrics": metrics,
            })

    return candidate_results

if __name__ == "__main__":
    results = run_candidate_grid(
        burn_times_days=[3.0, 5.0, 7.0],
        y_burns_m_s=[-20.0, 0.0, 20.0],
    )

    for candidate in results:
        manoeuvre = candidate["manoeuvre"]
        metrics = candidate["metrics"]

        burn_time_days = manoeuvre.time_s / SECONDS_PER_DAY
        y_burn_m_s = manoeuvre.delta_velocity_km_s[1] * 1000.0
        total_delta_v_m_s = metrics["total_delta_v_km_s"] * 1000.0

        print(f"\nBurn time: {burn_time_days:.1f} days")
        print(f"Y-burn: {y_burn_m_s:+.1f} m/s")
        print(
            "Closest-approach altitude: "
            f"{metrics['closest_approach_altitude_km']:,.1f} km"
        )
        print(
            "Heliocentric speed gain: "
            f"{metrics['heliocentric_speed_change_km_s']:+.6f} km/s"
        )
        print(f"Total delta-v cost: {total_delta_v_m_s:.1f} m/s")
