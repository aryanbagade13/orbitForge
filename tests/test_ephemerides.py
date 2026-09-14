from datetime import datetime, timezone

import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import CartesianRepresentation
from astropy.time import Time

import gravity_assist.ephemerides as ephemerides


def test_provider_converts_units_names_and_elapsed_time(monkeypatch):
    monkeypatch.setattr(ephemerides, "download_file", lambda *a, **k: "local.bsp")
    calls = []

    def fake_state(body, time, ephemeris):
        calls.append((body, time, ephemeris))
        return (
            CartesianRepresentation([1000, 2000, 3000] * u.m),
            CartesianRepresentation([86400, 172800, 259200] * u.km / u.day),
        )

    monkeypatch.setattr(ephemerides, "get_body_barycentric_posvel", fake_state)
    # There is a leap second between this UTC epoch and 2017-01-01.
    provider = ephemerides.JplEphemerisProvider(
        datetime(2016, 12, 31, 23, 59, 59, tzinfo=timezone.utc)
    )
    state = provider.get_state(" Earth ", 2.0)
    np.testing.assert_allclose(state.position_km, [1, 2, 3])
    np.testing.assert_allclose(state.velocity_km_s, [1, 2, 3])
    assert calls[0][0] == "earth"
    assert calls[0][2] == "local.bsp"
    assert calls[0][1].scale == "tdb"
    expected = Time("2017-01-01T00:00:00", scale="utc")
    assert abs((calls[0][1] - expected).to_value(u.s)) < 1e-8


def test_provider_rejects_naive_epoch():
    with pytest.raises(ValueError, match="timezone-aware"):
        ephemerides.JplEphemerisProvider(datetime(2000, 1, 1))


@pytest.mark.parametrize("body,time,error", [
    ("unknown", 0.0, ValueError),
    ("Earth", np.nan, ValueError),
    ("Earth", np.inf, ValueError),
    ("Earth", "1", TypeError),
])
def test_provider_rejects_invalid_query(monkeypatch, body, time, error):
    monkeypatch.setattr(ephemerides, "download_file", lambda *a, **k: "local.bsp")
    provider = ephemerides.JplEphemerisProvider(
        datetime(2000, 1, 1, tzinfo=timezone.utc)
    )
    with pytest.raises(error):
        provider.get_state(body, time)
