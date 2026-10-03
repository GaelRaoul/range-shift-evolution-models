/* Browser-only Python engine for the three-model simulator.
   Pyodide 0.25.1 is intentionally pinned because it ships NumPy 1.26.4,
   matching the NumPy version used for the archived scientific code. */

const PYODIDE_VERSION = "0.25.1";
const PYODIDE_BASE = `https://cdn.jsdelivr.net/pyodide/v${PYODIDE_VERSION}/full/`;

importScripts(`${PYODIDE_BASE}pyodide.js`);

let pyodide = null;
let ready = false;

function status(message) {
  self.postMessage({ type: "status", message });
}

async function fetchText(name) {
  const response = await fetch(new URL(name, self.location.href));
  if (!response.ok) throw new Error(`Could not load ${name}: HTTP ${response.status}`);
  return await response.text();
}

async function initialise() {
  status("Loading the browser Python engine…");
  pyodide = await loadPyodide({ indexURL: PYODIDE_BASE });

  status("Loading NumPy and SciPy…");
  await pyodide.loadPackage(["numpy", "scipy"]);

  status("Loading the three simulation models…");
  const files = ["sexual_model.py", "unified_models.py", "browser_bridge.py"];
  for (const name of files) {
    const source = await fetchText(name);
    pyodide.FS.writeFile(`/home/pyodide/${name}`, source);
  }

  await pyodide.runPythonAsync(`
import sys
if "/home/pyodide" not in sys.path:
    sys.path.insert(0, "/home/pyodide")
import numpy, scipy
import browser_bridge
`);

  const environmentJson = pyodide.runPython(`
import json, sys, numpy, scipy
json.dumps({
    "python": sys.version.split()[0],
    "numpy": numpy.__version__,
    "scipy": scipy.__version__,
    "pyodide": "${PYODIDE_VERSION}",
})
`);

  ready = true;
  self.postMessage({ type: "ready", environment: JSON.parse(environmentJson) });
}

const readyPromise = initialise().catch((error) => {
  self.postMessage({ type: "fatal", error: error?.stack || error?.message || String(error) });
  throw error;
});

self.onmessage = async (event) => {
  const { id, action, payload } = event.data || {};
  if (action !== "compare") return;

  try {
    await readyPromise;
    if (!ready) throw new Error("Python engine is not ready");

    status("Running simulations in your browser…");
    pyodide.globals.set("_browser_payload_json", JSON.stringify(payload));
    const resultJson = await pyodide.runPythonAsync(`
import browser_bridge
browser_bridge.run_compare_json(_browser_payload_json)
`);
    pyodide.globals.delete("_browser_payload_json");
    self.postMessage({ type: "result", id, result: JSON.parse(resultJson) });
  } catch (error) {
    self.postMessage({ type: "error", id, error: error?.stack || error?.message || String(error) });
  }
};
