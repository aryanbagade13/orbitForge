import matplotlib.pyplot as plt
import numpy as np

from legacy.experiments.jupiter_flyby_scenario import (
    SECONDS_PER_DAY,
    calculate_flyby_comparison_metrics,
    propagate_jupiter_flyby,
)
from gravity_assist.manoeuvres import ImpulsiveManoeuvre

BURN_TIME_DAYS = 5.0


def run_sensitivity_sweep(burns_m_s):
    burns_m_s = np.asarray(burns_m_s, dtype=float)
    closest_approach_altitudes_km = []
    turning_angles_degrees = []
    heliocentric_speed_changes_km_s = []

    for burn_m_s in burns_m_s:
        manoeuvres = []
        if burn_m_s != 0.0:
            manoeuvres = [
                ImpulsiveManoeuvre(
                    time_s=BURN_TIME_DAYS * SECONDS_PER_DAY,
                    delta_velocity_km_s=np.array([
                        0.0,
                        burn_m_s / 1000.0,
                        0.0,
                    ]),
                )
            ]

        result = propagate_jupiter_flyby(manoeuvres)
        metrics = calculate_flyby_comparison_metrics(
            result.times_s, result.states
        )

        closest_approach_altitudes_km.append(
            metrics["closest_approach_altitude_km"]
        )
        turning_angles_degrees.append(metrics["turning_angle_degrees"])
        heliocentric_speed_changes_km_s.append(
            metrics["heliocentric_speed_change_km_s"]
        )

    return {
        "burns_m_s": burns_m_s,
        "closest_approach_altitudes_km": np.asarray(
            closest_approach_altitudes_km
        ),
        "turning_angles_degrees": np.asarray(turning_angles_degrees),
        "heliocentric_speed_changes_km_s": np.asarray(
            heliocentric_speed_changes_km_s
        ),
    }


def plot_sensitivity_results(results):
    burns_m_s = results["burns_m_s"]
    figure, axes = plt.subplots(3, 1, figsize=(9, 11), sharex=True)

    axes[0].plot(
        burns_m_s,
        results["closest_approach_altitudes_km"],
        marker="o",
    )
    axes[0].set_ylabel("Closest-approach altitude (km)")

    axes[1].plot(
        burns_m_s,
        results["turning_angles_degrees"],
        marker="o",
        color="tab:orange",
    )
    axes[1].set_ylabel("Turning angle (degrees)")

    axes[2].plot(
        burns_m_s,
        results["heliocentric_speed_changes_km_s"],
        marker="o",
        color="tab:green",
    )
    axes[2].set_ylabel("Heliocentric speed change (km/s)")
    axes[2].set_xlabel("Day-5 y-directed burn (m/s)")

    for axis in axes:
        axis.axvline(0.0, color="black", linewidth=1, alpha=0.5)
        axis.grid(alpha=0.3)

    figure.suptitle("Jupiter flyby sensitivity to a day-5 manoeuvre")
    figure.tight_layout()
    return figure


def print_sensitivity_results(results):
    print(
        "Burn (m/s) | Closest altitude (km) | "
        "Turning angle (deg) | Heliocentric speed change (km/s)"
    )
    for burn_m_s, altitude_km, angle_degrees, speed_change_km_s in zip(
        results["burns_m_s"],
        results["closest_approach_altitudes_km"],
        results["turning_angles_degrees"],
        results["heliocentric_speed_changes_km_s"],
    ):
        print(
            f"{burn_m_s:10.1f} | {altitude_km:21.1f} | "
            f"{angle_degrees:19.3f} | {speed_change_km_s:+32.6f}"
        )


def main():
    burns_m_s = np.arange(-50.0, 51.0, 10.0)
    results = run_sensitivity_sweep(burns_m_s)
    print_sensitivity_results(results)
    plot_sensitivity_results(results)
    plt.show()


if __name__ == "__main__":
    main()
