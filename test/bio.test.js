/**
 * Headless tests for Protein Flow Explorer's analysis routines.
 *
 *   node test/bio.test.js
 *
 * The app is an IIFE with no module exports, so we run it inside a stubbed
 * DOM and capture the internals via a probe injection.
 */
const fs = require("fs");
const path = require("path");
const vm = require("vm");

const ROOT = path.join(__dirname, "..");
const noop = () => {};
const el = () => ({
  value: "", textContent: "", hidden: false, disabled: false,
  style: {}, innerHTML: "", dataset: {}, remove: noop,
  classList: { add: noop, remove: noop, toggle: noop, contains: () => false },
  addEventListener: noop, appendChild: noop,
  querySelector: () => el(), querySelectorAll: () => [],
  setAttribute: noop, removeAttribute: noop, focus: noop, click: noop,
});

global.document = {
  getElementById: () => el(), querySelectorAll: () => [],
  querySelector: () => ({ dataset: { filter: "all" } }),
  addEventListener: noop, readyState: "loading", createElement: () => el(),
};
global.window = {
  location: { search: "", pathname: "/" },
  addEventListener: noop, history: { replaceState: noop }, prompt: noop,
};
global.navigator = { clipboard: null };
global.performance = { now: () => 0 };
global.requestAnimationFrame = noop;
global.setTimeout = () => 0;
global.clearTimeout = noop;
global.fetch = () => Promise.reject(new Error("offline in tests"));
global.$3Dmol = {
  createViewer: () => ({
    setStyle: noop, removeAllModels: noop, removeAllShapes: noop,
    removeAllLabels: noop, removeAllSurfaces: noop, zoomTo: noop, render: noop,
    setBackgroundColor: noop, resize: noop, spin: noop, querySelector: () => el(),
    addStyle: noop, addModel: () => ({}), superpose: () => true,
    setClickable: noop, pngURI: () => "", addSphere: noop, addCylinder: noop,
    addLabel: noop, addSurface: noop, removeModel: noop, getIntersection: () => null,
  }),
  SurfaceType: { VDW: 1 },
};

vm.runInThisContext(fs.readFileSync(path.join(ROOT, "data/proteins.js"), "utf8"),
  { filename: "proteins.js" });

const src = fs.readFileSync(path.join(ROOT, "app.js"), "utf8");
const probe = src.replace(/\n\}\)\(\);\s*$/,
  "\n  global.__PROBE__ = { plddtToColor, midPoint, computeDistance, computeAngle,\n" +
  "    computeDihedral, cross, dot, progressiveMsa, conservationScores, nwAlign,\n" +
  "    resolveColor, chainColorFunc, STRUCT_CHAIN_PALETTE };\n})();");
vm.runInThisContext(probe, { filename: "app.js" });
const P = global.__PROBE__;

let pass = 0, fail = 0;
const near = (a, b, tol) => Math.abs(a - b) <= tol;
function check(name, cond, detail = "") {
  if (cond) { pass++; console.log(`  PASS  ${name}${detail ? "  " + detail : ""}`); }
  else { fail++; console.log(`  FAIL  ${name}${detail ? "  " + detail : ""}`); }
}

console.log("\npLDDT thresholds (AlphaFold convention)");
check("pLDDT 95 -> deep blue",  P.plddtToColor(95) === 0x0053d6);
check("pLDDT 85 -> cyan",       P.plddtToColor(85) === 0x65cbf3);
check("pLDDT 60 -> yellow",     P.plddtToColor(60) === 0xffdb13);
check("pLDDT 20 -> orange",     P.plddtToColor(20) === 0xff7d45);

console.log("\nGeometry");
check("distance 3-4-0 = 5", near(P.computeDistance({x:0,y:0,z:0},{x:3,y:4,z:0}), 5, 1e-9));
check("angle 90deg",         near(P.computeAngle({x:1,y:0,z:0},{x:0,y:0,z:0},{x:0,y:1,z:0}), 90, 1e-6));
check("dihedral trans = 180", near(Math.abs(P.computeDihedral({x:1,y:0,z:0},{x:0,y:0,z:0},{x:0,y:1,z:0},{x:-1,y:1,z:0})), 180, 1));
check("dihedral cis = 0",     near(P.computeDihedral({x:1,y:0,z:0},{x:0,y:0,z:0},{x:0,y:1,z:0},{x:1,y:1,z:0}), 0, 1));

console.log("\nChain palette (regression: black-protein bug)");
const pal = P.STRUCT_CHAIN_PALETTE;
const lum = pal.map((v) => {
  const r = (v >> 16) & 255, g = (v >> 8) & 255, b = v & 255;
  return 0.2126 * r + 0.7152 * g + 0.0722 * b;
});
check("palette has >= 8 colours", pal.length >= 8, `(${pal.length})`);
check("all distinct",            new Set(pal).size === pal.length);
check("all visible on white",    lum.every((l) => l > 40),
      `(min lum ${Math.min(...lum).toFixed(0)})`);

console.log("\nColour routing");
check("chain -> colorfunc", typeof P.resolveColor("chain", 1).colorfunc === "function");
check("spectrum passes through", P.resolveColor("spectrum", 1).color === "spectrum");
check("colorfunc returns palette colour",
      pal.includes(P.chainColorFunc({ chain: "A" })));

console.log("\nMSA / conservation");
const fam = [
  "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGE",
  "MKWVTFISLLLLFSSAYSRGVFRRDAHKSEIAHRFKDLGE",
  "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGE",
  "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE",
  "MKWVTFLSLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGE",
  "MQWVTFLLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGE",
];
const m = P.progressiveMsa(fam);
check("all rows aligned",  m.aligned.every((r) => r.length === m.consensus.length));
const sc = P.conservationScores(m);
check("scores in [0,1]",   sc.every((v) => v >= 0 && v <= 1));
check("no NaN",            sc.every((v) => !Number.isNaN(v)));
const identical = P.nwAlign("ACDEFGHIK", "ACDEFGHIK");
check("identical seqs ungapped", !identical.a.includes("-") && !identical.b.includes("-"));

console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
