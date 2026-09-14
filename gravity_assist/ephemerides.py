"""Interface for planetary positions and velocities over mission time."""

from datetime import datetime
from numbers import Real
from typing import Protocol

import numpy as np
from astropy import units as u
from astropy.coordinates import get_body_barycentric_posvel
from astropy.time import Time, TimeDelta
from astropy.utils.data import download_file

from .models import OrbitalState


class EphemerisProvider(Protocol):
    """Contract implemented by a planetary-state provider.

    The departure epoch must be timezone-aware. Mission time is elapsed
    seconds from that instant; concrete providers handle conversion to their
    data source's time standard.

    Providers must document their origin and axis orientation. All returned
    states, and the spacecraft states compared with them, must share those
    conventions. Positions are in km and velocities are in km/s.
    """

    departure_epoch: datetime
    frame: str

    def get_state(self, body_name: str, time_s: float) -> OrbitalState:
        """Return a body's state at the requested elapsed mission time."""
        ...


class JplEphemerisProvider:
    """DE432s geometric states: Solar System barycentric origin, ICRS axes.

    Earth is its physical centre; Jupiter and Saturn refer to their system
    barycentres, not planet centres. Do not use the latter for precision
    surface-clearance calculations. This frame differs from the existing
    heliocentric Jupiter experiment; its initial state cannot be reused here.

    The kernel is downloaded once into Astropy's cache. All state queries
    subsequently use the local file, with no Horizons/network request.
    Mission seconds are elapsed SI seconds (TAI), converted to TDB for JPL.
    DE432s covers approximately 1950-2050; out-of-range queries raise errors.
    """

    frame = "ICRS axes; Solar System barycentric origin"
    ephemeris_name = "de432s"
    kernel_url = (
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/"
        "spk/planets/de432s.bsp"
    )
    supported_bodies = ("sun", "earth", "jupiter", "saturn")

    def __init__(self, departure_epoch: datetime):
        if not isinstance(departure_epoch, datetime):
            raise TypeError("departure_epoch must be a datetime")
        if departure_epoch.utcoffset() is None:
            raise ValueError("departure_epoch must be timezone-aware")
        self.departure_epoch = departure_epoch
        self._epoch = Time(departure_epoch, scale="utc").tai
        self._kernel_path = download_file(self.kernel_url, cache=True)

    def get_state(self, body_name: str, time_s: float) -> OrbitalState:
        """Return geometric position (km) and velocity (km/s), without light time."""
        if not isinstance(body_name, str):
            raise TypeError("body_name must be a string")
        body_name = body_name.strip().lower()
        if body_name not in self.supported_bodies:
            raise ValueError(f"unsupported body: {body_name!r}")
        if isinstance(time_s, (bool, np.bool_)) or not isinstance(time_s, Real):
            raise TypeError("time_s must be a real number")
        if not np.isfinite(time_s):
            raise ValueError("time_s must be finite")
        query_time = (self._epoch + TimeDelta(float(time_s), format="sec")).tdb
        position, velocity = get_body_barycentric_posvel(
            body_name, query_time, ephemeris=self._kernel_path
        )
        return OrbitalState(
            position_km=position.xyz.to_value(u.km),
            velocity_km_s=velocity.xyz.to_value(u.km / u.s),
        )

    def get_states(self, body_name: str, times_s) -> np.ndarray:
        """Vectorised queries, shape (n, 6), for building interpolation tables."""
        times = np.asarray(times_s, dtype=float)
        if times.ndim != 1 or not np.isfinite(times).all():
            raise ValueError("times_s must be a finite one-dimensional array")
        name = body_name.strip().lower()
        if name not in self.supported_bodies:
            raise ValueError(f"unsupported body: {name!r}")
        time = (self._epoch + TimeDelta(times, format="sec")).tdb
        position, velocity = get_body_barycentric_posvel(name, time, ephemeris=self._kernel_path)
        return np.concatenate((position.xyz.to_value(u.km).T,
                               velocity.xyz.to_value(u.km/u.s).T), axis=1)


class TabulatedEphemerisProvider:
    """Hermite interpolation of position using supplied velocities.

    No extrapolation. Midpoint checks measure interpolation error, not the
    absolute physical error of the source ephemeris. Values at other points
    are not guaranteed to obey the midpoint error bound.
    """

    def __init__(self, source, body_names, end_time_s, spacing_s=21600.0,
                 position_tolerance_km=0.1, velocity_tolerance_km_s=1e-5):
        from scipy.interpolate import CubicHermiteSpline

        if not np.isfinite([end_time_s, spacing_s, position_tolerance_km,
                            velocity_tolerance_km_s]).all() or min(
                            end_time_s, spacing_s, position_tolerance_km,
                            velocity_tolerance_km_s) <= 0:
            raise ValueError("table duration, spacing and tolerances must be positive and finite")
        self.departure_epoch = source.departure_epoch
        self.frame = source.frame
        self.end_time_s = float(end_time_s)
        self._splines = {}
        self.validation_errors = {}
        times = np.linspace(0., end_time_s, int(np.ceil(end_time_s/spacing_s))+1)
        midpoints = (times[:-1] + times[1:])/2

        def query(name, values):
            if hasattr(source, "get_states"):
                return source.get_states(name, values)
            return np.array([np.concatenate((s.position_km, s.velocity_km_s))
                             for t in values for s in [source.get_state(name, t)]])

        for name in body_names:
            name = name.strip().lower()
            if name in self._splines:
                raise ValueError("duplicate body in ephemeris table")
            states = query(name, times)
            spline = CubicHermiteSpline(times, states[:, :3], states[:, 3:], extrapolate=False)
            truth = query(name, midpoints)
            position_error = float(np.max(np.linalg.norm(spline(midpoints)-truth[:, :3], axis=1)))
            velocity_error = float(np.max(np.linalg.norm(spline(midpoints, 1)-truth[:, 3:], axis=1)))
            if position_error > position_tolerance_km or velocity_error > velocity_tolerance_km_s:
                raise ValueError(f"ephemeris interpolation too coarse for {name}: "
                                 f"{position_error} km, {velocity_error} km/s; reduce spacing_s")
            self._splines[name] = spline
            self.validation_errors[name] = dict(position_km=position_error,
                                                 velocity_km_s=velocity_error)

    def get_state(self, body_name, time_s):
        if not np.isfinite(time_s) or not 0 <= time_s <= self.end_time_s:
            raise ValueError("query outside ephemeris table coverage")
        name = body_name.strip().lower()
        if name not in self._splines:
            raise ValueError(f"body {name!r} is not tabulated")
        spline = self._splines[name]
        return OrbitalState(spline(time_s), spline(time_s, 1))
