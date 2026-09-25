"use strict";

// The viewer only reads saved results. Camera and playback never run the solver.
const $ = id => document.getElementById(id);
const canvas = $("orbit"), ctx = canvas.getContext("2d");
const AU = 149597870.7, YEAR = 365.25 * 86400;
const colours = {no_burn: "#b9ecc6", correction: "#b1a0e0", earth_positions_km: "#79b2d2", jupiter_positions_km: "#cea079", saturn_positions_km: "#c9c29a"};
const names = {no_burn: "No burn", correction: "Day-180 correction", earth_positions_km: "Earth", jupiter_positions_km: "Jupiter", saturn_positions_km: "Saturn"};
const planets = ["earth_positions_km", "jupiter_positions_km", "saturn_positions_km"];
let data, width = 1, height = 1, selected = "no_burn", playing = false;
let time = 0, yaw = -0.35, tilt = -0.45, zoom = 1, view = "3d", bounds = 12;
let centre = [0, 0, 0], previousFrame = null, labelBoxes = [];
const norm = p => Math.hypot(...p);

function sampleIndex(t) {
  const times = data.times_s;
  let low = 0, high = times.length - 1;
  while (high - low > 1) {
    const mid = (low + high) >> 1;
    if (times[mid] <= t) low = mid; else high = mid;
  }
  return [low, Math.max(0, Math.min(1, (t - times[low]) / (times[high] - times[low])))];
}

function interpolate(values, index, mix) {
  return values[index].map((v, axis) => v + (values[index + 1][axis] - v) * mix);
}

function project(position) {
  const [x, y, z] = position.map((v, axis) => v - centre[axis]);
  const u = x * Math.cos(yaw) - y * Math.sin(yaw);
  const v = x * Math.sin(yaw) + y * Math.cos(yaw);
  const scale = Math.min(width - 65, height - 140) / (bounds * 2) * zoom;
  return [width / 2 + u * scale, height * 0.55 - (v * Math.cos(tilt) - z * Math.sin(tilt)) * scale];
}

function line(points, colour, opacity = 1, weight = 1, dash = []) {
  ctx.beginPath(); ctx.strokeStyle = colour; ctx.globalAlpha = opacity; ctx.lineWidth = weight; ctx.setLineDash(dash);
  points.forEach((p, i) => { const [x, y] = project(p); i ? ctx.lineTo(x, y) : ctx.moveTo(x, y); });
  ctx.stroke(); ctx.setLineDash([]); ctx.globalAlpha = 1;
}

function marker(position, colour, radius, label, spacecraft = false) {
  const [x, y] = project(position);
  ctx.beginPath(); ctx.fillStyle = colour;
  if (spacecraft) { ctx.moveTo(x, y - 5); ctx.lineTo(x + 4, y + 4); ctx.lineTo(x - 4, y + 4); ctx.closePath(); }
  else ctx.arc(x, y, radius, 0, Math.PI * 2);
  ctx.fill();
  if ($("labels").checked && label) {
    ctx.font = "10px ui-monospace, monospace"; ctx.fillStyle = colour;
    const textWidth = ctx.measureText(label).width;
    const lx = Math.max(12, Math.min(width - textWidth - 12, x + radius + 9));
    let ly = Math.max(55, Math.min(height - 40, y - radius - 3));
    for (let attempt = 0; attempt < 12; attempt++) {
      if (!labelBoxes.some(b => lx < b[0] + b[2] + 5 && lx + textWidth + 5 > b[0] && ly - 10 < b[1] + 4 && ly + 4 > b[1] - 10)) break;
      ly += 15;
    }
    labelBoxes.push([lx, ly, textWidth]);
    ctx.fillText(label, lx, ly);
  }
}

function render() {
  if (!data) return;
  ctx.clearRect(0, 0, width, height); labelBoxes = [];
  // Reference grid lies in the ICRS XY plane, not the ecliptic.
  const gridExtent = Math.ceil(bounds / 5) * 5;
  for (let n = -gridExtent; n <= gridExtent; n += 5) {
    line([[n, -gridExtent, 0], [n, gridExtent, 0]], "#718e7a", n === 0 ? 0.13 : 0.055);
    line([[-gridExtent, n, 0], [gridExtent, n, 0]], "#718e7a", n === 0 ? 0.13 : 0.055);
  }
  const [i, mix] = sampleIndex(time);
  if ($("paths").checked) planets.forEach(key => line(data.positions[key], colours[key], 0.38, 1));
  const other = selected === "no_burn" ? "correction" : "no_burn";
  if ($("compare").checked) {
    line(data.positions[other], colours[other], 0.7, 1, [4, 5]);
    marker(interpolate(data.positions[other], i, mix), colours[other], 2, "", true);
  }
  line(data.positions[selected], colours[selected], 0.28, 1.3);
  const current = interpolate(data.positions[selected], i, mix);
  line([...data.positions[selected].slice(0, i + 1), current], colours[selected], 1, 1.8);
  marker([0, 0, 0], "#efdca2", 5, "Sun");
  planets.forEach(key => marker(interpolate(data.positions[key], i, mix), colours[key], key.startsWith("earth") ? 3 : 4, names[key]));
  marker(current, colours[selected], 4, "Spacecraft", true);
  // An unobtrusive orientation triad shows the projection of the actual axes.
  const origin = project([0, 0, 0]);
  [[1, 0, 0], [0, 1, 0], [0, 0, 1]].forEach((axis, a) => {
    const p = project(axis), dx = p[0] - origin[0], dy = p[1] - origin[1];
    const length = Math.hypot(dx, dy) || 1;
    ctx.strokeStyle = "#667a6d"; ctx.fillStyle = "#829388"; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(42, height - 70); ctx.lineTo(42 + dx / length * 20, height - 70 + dy / length * 20); ctx.stroke();
    ctx.font = "8px monospace"; ctx.fillText("XYZ"[a], 39 + dx / length * 29, height - 67 + dy / length * 29);
  });
  $("scale").textContent = "Grid: 5 AU";
}

function updateTelemetry() {
  if (!data) return;
  const [i, mix] = sampleIndex(time), position = interpolate(data.positions[selected], i, mix);
  $("distance").innerHTML = norm(position).toFixed(3) + " <small>AU</small>";
  const speeds = data.speeds_km_s[selected];
  $("speed").innerHTML = (speeds[i] + (speeds[i + 1] - speeds[i]) * mix).toFixed(2) + " <small>km/s</small>";
  ["pos-x", "pos-y", "pos-z"].forEach((id, axis) => $(id).textContent = (position[axis] >= 0 ? "+" : "") + position[axis].toFixed(4));
  $("day").textContent = "Day " + Math.floor(time / 86400);
  $("view-year").textContent = "Year " + (time / YEAR).toFixed(2);
  $("view-date").textContent = new Date(Date.parse(data.report.departure_epoch) + time * 1000).toLocaleDateString("en-GB", {day:"2-digit", month:"short", year:"numeric", timeZone:"UTC"}) + " · UTC";
  $("time").value = time;
  const fraction = (time - data.times_s[0]) / (data.times_s.at(-1) - data.times_s[0]);
  $("time").style.setProperty("--progress", fraction * 100 + "%");
  $("chart-cursor").setAttribute("x1", fraction * 220); $("chart-cursor").setAttribute("x2", fraction * 220);
  render();
}

function updateRoute() {
  document.querySelectorAll("[data-route]").forEach(button => {
    const active = button.dataset.route === selected;
    button.classList.toggle("active", active); button.setAttribute("aria-pressed", active);
  });
  const evaluation = data.report.evaluations[names[selected]];
  $("arrival").innerHTML = (evaluation.arrival_distance_km / AU).toFixed(3) + " <small>AU</small>";
  $("feasibility").textContent = evaluation.feasible ? "Feasible" : "Infeasible";
  $("result-note").textContent = evaluation.feasible ? "Saved candidate meets the evaluated constraints." : "Mission constraints not met. Unoptimised baseline.";
  const values = data.positions[selected].map(norm), maximum = Math.max(...values), minimum = Math.min(...values);
  $("distance-path").setAttribute("d", values.map((v, i) => `${i ? "L" : "M"}${(i / (values.length - 1) * 220).toFixed(2)},${(34 - (v - minimum) / (maximum - minimum || 1) * 29).toFixed(2)}`).join(" "));
  $("distance-path").style.stroke = colours[selected];
  updateTelemetry();
}

function stopPlayback() { playing = false; $("play").textContent = "▶"; $("play").setAttribute("aria-label", "Play trajectory"); }
function seek(value) { if (!data) return; stopPlayback(); time = Math.max(data.times_s[0], Math.min(data.times_s.at(-1), value)); updateTelemetry(); }
function setView(mode) {
  view = mode; yaw = mode === "xy" ? 0 : -0.35; tilt = mode === "xy" ? 0 : -0.45;
  $("top-view").classList.toggle("selected", mode === "xy"); $("top-view").setAttribute("aria-pressed", mode === "xy");
  $("perspective").classList.toggle("selected", mode === "3d"); $("perspective").setAttribute("aria-pressed", mode === "3d"); render();
}

$("time").addEventListener("input", event => seek(Number(event.target.value)));
$("play").addEventListener("click", () => {
  if (!data) return;
  if (playing) return stopPlayback();
  if (time >= data.times_s.at(-1)) time = data.times_s[0];
  playing = true; previousFrame = null; $("play").textContent = "Ⅱ"; $("play").setAttribute("aria-label", "Pause trajectory");
});
document.querySelectorAll("[data-route]").forEach(button => button.addEventListener("click", () => { if (!data) return; selected = button.dataset.route; updateRoute(); }));
["paths", "labels", "compare"].forEach(id => $(id).addEventListener("change", render));
$("departure").onclick = () => seek(data?.times_s[0]);
$("burn").onclick = () => { if (!data) return; selected = "correction"; updateRoute(); seek(180 * 86400); };
$("end").onclick = () => seek(data?.times_s.at(-1));
$("perspective").onclick = () => setView("3d"); $("top-view").onclick = () => setView("xy");
$("reset").onclick = () => { zoom = 1; setView("3d"); };
function changeZoom(factor) { zoom = Math.min(12, Math.max(0.35, zoom * factor)); render(); }
$("zoom-in").onclick = () => changeZoom(1.25); $("zoom-out").onclick = () => changeZoom(0.8);
$("rotate-left").onclick = () => { yaw -= 0.2; render(); }; $("rotate-right").onclick = () => { yaw += 0.2; render(); };
canvas.addEventListener("wheel", event => { event.preventDefault(); changeZoom(Math.exp(-event.deltaY * 0.001)); }, {passive:false});
let drag = null;
canvas.addEventListener("pointerdown", event => { drag = [event.clientX, event.clientY]; canvas.setPointerCapture(event.pointerId); });
canvas.addEventListener("pointermove", event => {
  if (!drag) return;
  yaw += (event.clientX - drag[0]) * 0.006;
  if (view === "3d") tilt = Math.max(-1.5, Math.min(1.5, tilt + (event.clientY - drag[1]) * 0.006));
  drag = [event.clientX, event.clientY]; render();
});
canvas.addEventListener("pointerup", () => drag = null); canvas.addEventListener("pointercancel", () => drag = null);
new ResizeObserver(() => {
  const box = canvas.getBoundingClientRect(), dpr = window.devicePixelRatio || 1;
  width = box.width; height = box.height; canvas.width = Math.round(width * dpr); canvas.height = Math.round(height * dpr);
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0); render();
}).observe(canvas);

function frame(now) {
  if (playing && data) {
    // One whole mission per 45 seconds at 1×, independent of frame rate.
    if (previousFrame !== null) time += Math.min(now - previousFrame, 100) / 45000 * (data.times_s.at(-1) - data.times_s[0]) * Number($("rate").value);
    if (time >= data.times_s.at(-1)) { time = data.times_s.at(-1); stopPlayback(); }
    updateTelemetry();
  }
  previousFrame = now; requestAnimationFrame(frame);
}
requestAnimationFrame(frame);

const trajectorySource = window.orbitForgeData ? Promise.resolve(window.orbitForgeData) : fetch("/api/trajectory").then(response => { if (!response.ok) throw new Error("Could not read saved trajectory."); return response.json(); });
trajectorySource.then(result => {
  data = result;
  // Fit the complete flight, including both candidates and all planet paths.
  const all = Object.values(data.positions).flat();
  const minimum = [Infinity, Infinity, Infinity], maximum = [-Infinity, -Infinity, -Infinity];
  all.forEach(point => point.forEach((v, axis) => { minimum[axis] = Math.min(minimum[axis], v); maximum[axis] = Math.max(maximum[axis], v); }));
  centre = minimum.map((v, axis) => (v + maximum[axis]) / 2);
  bounds = Math.max(...all.map(point => norm(point.map((v, axis) => v - centre[axis])))) * 1.1;
  $("time").min = data.times_s[0]; $("time").max = data.times_s.at(-1); time = data.times_s[0];
  $("time").disabled = false; $("play").disabled = false;
  $("time-ticks").replaceChildren(...Array.from({length:9}, (_, i) => { const span = document.createElement("span"); span.textContent = ((data.times_s[0] + (data.times_s.at(-1) - data.times_s[0]) * i / 8) / YEAR).toFixed(0) + "Y"; return span; }));
  $("sample-count").textContent = data.times_s.length.toLocaleString() + " saved samples · interpolated playback";
  $("model-details").textContent = data.report.ephemeris + ". " + data.report.limitations.join(". ") + ". Departure cost excluded.";
  $("data-note").textContent = "DE432s ephemerides / saved baseline / " + data.times_s.length.toLocaleString() + " samples";
  updateRoute();
}).catch(error => { $("error").hidden = false; $("error").textContent = "Unable to load the trajectory: " + error.message + " Restart the viewer with a valid baseline output folder."; $("view-date").textContent = "Data unavailable"; });
