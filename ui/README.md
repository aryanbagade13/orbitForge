# orbitForge trajectory viewer

A local desktop viewer for `experiments/interplanetary_baseline.py` outputs.
The macOS window uses pywebview and the system WebKit renderer. It reads the
saved simulation directly, without a browser tab, web server or internet access.

On this Mac, double-click `ui/orbitForge.app` in Finder to launch. This is a
local app launcher: keep it in the `ui` folder. It uses the repository's `.venv`
and saved outputs; it is not a self-contained app for distribution to other Macs.

From the orbitForge repository, using its virtual environment:

```sh
.venv/bin/python -m pip install -r ui/requirements.txt
.venv/bin/python ui/desktop.py
```

Close the window to quit. If running this folder outside the repository,
supply the saved output folder:

```sh
/path/to/Orbital_Dynamics/.venv/bin/python desktop.py --data-dir /path/to/Orbital_Dynamics/outputs/interplanetary_baseline
```

Missing outputs? Run the baseline experiment first. Restart the viewer after
regenerating the output files. `server.py` remains an optional browser preview
for development; it is not used by the desktop app.

## Controls

- Select a candidate; use Compare trajectories to overlay the other route.
- Drag to rotate and scroll to zoom. Camera buttons offer keyboard alternatives.
- XY shows the ICRS XY projection. The 3D view uses an orthographic camera.
- Scrub or play the mission; event buttons jump to departure, the correction or arrival.
- Read the selected-time position and speed alongside the saved full-flight evaluation.

## Small, separate layers

`server.py` contains the NPZ validation and optional development server.
It subtracts the saved Sun position from each position and converts km to AU.
`app.js` projects the resulting 3D coordinates onto a canvas, interpolates playback
and updates the measurements. `index.html` and `style.css` define the interface.
`desktop.py` bundles the local interface and data into a native window.
The physics engine is unchanged.

Playback linearly interpolates the saved samples; it is a visualisation, not
another integration. Displayed speed interpolates the magnitudes of the saved
barycentric velocity samples. Marker sizes are illustrative. Planet paths show
only the saved eight-year interval, so outer-planet paths are incomplete arcs.
The correction is small and the two trajectories can overlap at the full-flight
scale. The grid is in the ICRS XY plane, not the ecliptic.

The baseline trials miss their mission targets. The result panel reads their
saved evaluations; the viewer neither optimises nor changes a mission.
