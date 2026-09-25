"""Checks for the coordinate conversion at the physics/viewer boundary."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from server import AU_KM, load_data


class ViewerDataTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        (self.directory / "summary.json").write_text(json.dumps({"evaluations": {}}))
        self.arrays = {"times_s": np.array([0., 10.])}
        sun = np.array([[100., 200., 300.], [400., 500., 600.]])
        for body in ("sun", "earth", "jupiter", "saturn"):
            self.arrays[body + "_positions_km"] = sun.copy()
        for route in ("no_burn", "correction"):
            self.arrays[route] = np.column_stack((sun + [AU_KM, 0., 0.], [[3., 4., 0.], [0., 0., 12.]]))

    def load(self):
        np.savez(self.directory / "trajectory_samples.npz", **self.arrays)
        return load_data(self.directory)

    def test_subtracts_each_saved_sun_position_and_converts_to_au(self):
        result = self.load()
        np.testing.assert_allclose(result["positions"]["no_burn"], [[1, 0, 0], [1, 0, 0]])
        np.testing.assert_allclose(result["positions"]["sun_positions_km"], np.zeros((2, 3)))
        self.assertEqual(result["speeds_km_s"]["no_burn"], [5., 12.])

    def test_rejects_misaligned_rows(self):
        self.arrays["correction"] = self.arrays["correction"][:1]
        with self.assertRaisesRegex(ValueError, "correction"):
            self.load()

    def test_rejects_nan_positions(self):
        self.arrays["earth_positions_km"][0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, "earth_positions_km"):
            self.load()

    def test_rejects_nonincreasing_time(self):
        self.arrays["times_s"] = np.array([0., 0.])
        with self.assertRaisesRegex(ValueError, "increasing"):
            self.load()


if __name__ == "__main__":
    unittest.main()
