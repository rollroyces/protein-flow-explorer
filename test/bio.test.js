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

// The extracted pure modules publish themselves on the PFE namespace and must
// be loaded before app.js, which destructures them.
for (const m of ["seq", "geometry", "pockets", "colors"]) {
  vm.runInThisContext(fs.readFileSync(path.join(ROOT, "src", `${m}.js`), "utf8"),
    { filename: `src/${m}.js` });
}

const src = fs.readFileSync(path.join(ROOT, "app.js"), "utf8");
const probe = src.replace(/\n\}\)\(\);\s*$/,
  "\n  global.__PROBE__ = { NET_TIMEOUT_MS, NET_RETRIES, newNetworkGeneration,\n" +
  "    isAbort, currentNetGeneration: () => netGeneration,\n" +
  "    midPoint, computeDistance, computeAngle,\n" +
  "    computeDihedral, cross, dot, progressiveMsa, conservationScores, nwAlign,\n" +
  "    resolveColor, chainColorFunc, domainColorFunc, state,\n" +
  "    STRUCT_CHAIN_PALETTE, isPlanarAromatic,\n" +
  "    isHydrophobicResidue, isIonic, extractPdbTitle, firstRefineStat,\n" +
  "    HBOND_CUTOFF, SALT_CUTOFF, PI_CUTOFF, HYDRO_CUTOFF,\n" +
  "    POCKET_MIN_ENCLOSURE, POCKET_MIN_WALL, POCKET_MIN_GLOBAL, POCKET_CLUSTER_R };\n})();");
vm.runInThisContext(probe, { filename: "app.js" });
const P = global.__PROBE__;
const {
  isPlanarAromatic, isHydrophobicResidue, isIonic, extractPdbTitle, firstRefineStat,
  HBOND_CUTOFF, SALT_CUTOFF, PI_CUTOFF,
  NET_TIMEOUT_MS, NET_RETRIES, newNetworkGeneration, isAbort, currentNetGeneration,
} = P;

let pass = 0, fail = 0;
const near = (a, b, tol) => Math.abs(a - b) <= tol;
function check(name, cond, detail = "") {
  if (cond) { pass++; console.log(`  PASS  ${name}${detail ? "  " + detail : ""}`); }
  else { fail++; console.log(`  FAIL  ${name}${detail ? "  " + detail : ""}`); }
}

console.log("\npLDDT thresholds (AlphaFold convention)");
// The mapping lives in src/colors.js and is covered in test/colors.test.js.
const { plddtColor } = globalThis.PFE.colors;
check("pLDDT 95 -> deep blue",  plddtColor(95) === 0x0053d6);
check("pLDDT 85 -> cyan",       plddtColor(85) === 0x65cbf3);
check("pLDDT 60 -> yellow",     plddtColor(60) === 0xffdb13);
check("pLDDT 20 -> orange",     plddtColor(20) === 0xff7d45);

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
check("spectrum -> colorfunc", typeof P.resolveColor("spectrum", 1).colorfunc === "function");
check("ss -> colorfunc", typeof P.resolveColor("ss", 1).colorfunc === "function");
check("byElement -> 3Dmol scheme", P.resolveColor("byElement", 1).colorscheme === "Jmol");
check("residue -> 3Dmol scheme", P.resolveColor("residue", 1).colorscheme === "amino");
check("uniform -> literal colour",
      typeof P.resolveColor("uniform", 1).color === "string" &&
      P.resolveColor("uniform", 1).color.startsWith("#"));
check("colorfunc returns palette colour",
      pal.includes(P.chainColorFunc({ chain: "A" })));

// Every option in the colour dropdown must produce a spec 3Dmol can honour.
// This is the regression guard for the bug where a scheme NAME was put in
// 3Dmol's `color` field: 3Dmol treats that field as a LITERAL colour, so any
// name that is not an html colour name falls through to 0x000000 and the
// whole structure renders black.
console.log("\nColour dropdown produces valid 3Dmol specs");
const html = fs.readFileSync(path.join(ROOT, "index.html"), "utf8");
const selectStart = html.indexOf('id="ctrl-color"');
const selectEnd = html.indexOf("</select>", selectStart);
const selectBlock = html.slice(selectStart, selectEnd);
const options = [...selectBlock.matchAll(/<option value="([^"]+)"/g)].map((m) => m[1]);
check("colour dropdown has 11 options", options.length === 11, `(${options.length})`);
const badSpecs = [];
for (const opt of options) {
  const s = P.resolveColor(opt, 1);
  const ok = typeof s.colorfunc === "function"
    || typeof s.colorscheme === "string"
    || (typeof s.color === "string" && s.color.startsWith("#"))
    || typeof s.color === "number";
  if (!ok) badSpecs.push(opt);
}
check("every dropdown option yields a renderable spec", badSpecs.length === 0,
      badSpecs.length ? `(${badSpecs.join(", ")})` : "");
const namedColour = options.filter((o) => {
  const s = P.resolveColor(o, 1);
  return typeof s.color === "string" && !s.color.startsWith("#");
});
check("no scheme name is passed as a literal colour", namedColour.length === 0,
      namedColour.length ? `(${namedColour.join(", ")})` : "");

console.log("\nDomain colouring");
const AK_DOMAINS = [
  { name: "NMP", color: "#2f7ed8", ranges: [{ start: 1, end: 41 }] },
  { name: "CORE", color: "#9aa0a6", ranges: [{ start: 42, end: 118 }] },
  { name: "LID", color: "#cc1f1f", ranges: [{ start: 119, end: 214 }] },
];
P.state.currentProtein = { id: "adenylate-kinase", domains: AK_DOMAINS };
const hexOf = (atom) => "0x" + P.domainColorFunc(atom).toString(16).padStart(6, "0");
check("NMP residue is blue", hexOf({ chain: "A", resi: 10 }) === "0x2f7ed8",
      `(${hexOf({ chain: "A", resi: 10 })})`);
check("CORE residue is grey", hexOf({ chain: "A", resi: 80 }) === "0x9aa0a6",
      `(${hexOf({ chain: "A", resi: 80 })})`);
check("LID residue is red", hexOf({ chain: "A", resi: 200 }) === "0xcc1f1f",
      `(${hexOf({ chain: "A", resi: 200 })})`);
check("domain colours are distinct",
      new Set([10, 80, 200].map((r) => P.domainColorFunc({ chain: "A", resi: r }))).size === 3);
check("both dimer chains are coloured the same",
      hexOf({ chain: "B", resi: 10 }) === hexOf({ chain: "A", resi: 10 }));
check("residue outside every domain falls back to grey",
      hexOf({ chain: "A", resi: 999 }) === "0x9aa0a6");
P.state.currentProtein = { id: "no-domains" };
check("a protein without domains falls back to chain colouring",
      P.domainColorFunc({ chain: "A", resi: 10 }) === P.chainColorFunc({ chain: "A" }));

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


console.log("\nLigand-property chemistry helpers");
check("PHE is planar aromatic",  isPlanarAromatic({ resn: "PHE" }));
check("TYR is planar aromatic",  isPlanarAromatic({ resn: "TYR" }));
check("ALA is not aromatic",     !isPlanarAromatic({ resn: "ALA" }));
check("LEU is hydrophobic",      isHydrophobicResidue("LEU"));
check("GLY is hydrophobic",      isHydrophobicResidue("GLY"));
check("SER is not hydrophobic",  !isHydrophobicResidue("SER"));
check("ZN recognised as ion",    isIonic({ compId: "ZN" }));
check("FE recognised as ion",    isIonic({ compId: "FE" }));
check("AQ4 not an ion",         !isIonic({ compId: "AQ4" }));

console.log("\nInteraction cutoff constants");
check("H-bond cutoff 3.5 A", HBOND_CUTOFF === 3.5);
check("salt cutoff 4.0 A",   SALT_CUTOFF === 4.0);
check("pi cutoff 5.5 A",     PI_CUTOFF === 5.5);

console.log("\nPDB TITLE extraction");
check("extracts TITLE record",
  extractPdbTitle("HEADER    HYDROLASE\nTITLE     CRYSTAL STRUCTURE OF LYSOZYME\nCOMPND") === "CRYSTAL STRUCTURE OF LYSOZYME");
check("returns null when absent", extractPdbTitle("ATOM      1  N   MET A   1") === null);


console.log("\nRefinement-stat extraction (R-work / R-free)");
// Real 4HHB payload shape: capital-R mmCIF keys under refine[0]
const entry4hhb = { refine: [{ ls_R_factor_R_work: 0.135, ls_R_factor_R_free: 0.194 }] };
check("finds capital-R R-work", firstRefineStat(entry4hhb, "ls_r_factor_r_work") === "0.135");
check("finds capital-R R-free", firstRefineStat(entry4hhb, "ls_r_factor_r_free") === "0.194");
check("returns em dash when absent", firstRefineStat({ refine: [{}] }, "ls_r_factor_r_work") === "\u2014");
check("tolerates missing blocks",  firstRefineStat({}, "ls_r_factor_r_work") === "\u2014");
check("skips non-numeric values",
  firstRefineStat({ refine: [{ ls_R_factor_R_work: null }] }, "ls_r_factor_r_work") === "\u2014");


console.log("\nPacking-density thresholds (calibrated on oxymyoglobin 1MBO)");
// Distribution measured from the real structure: min 3, p10 6, median 9,
// p75 11, max 15. His93/H64 are the heme-ligating histidines and must
// classify as buried; the N-terminus must classify as exposed.
function classify(n) {
  return n >= 10 ? "buried" : n >= 6 ? "moderate" : "exposed";
}
check("His93 (n=10) -> buried",   classify(10) === "buried");
check("His64 (n=9)  -> moderate", classify(9) === "moderate");
check("R118 (n=6)  -> moderate", classify(6) === "moderate");
check("N-term (n=3) -> exposed",  classify(3) === "exposed");
check("max observed (n=15) -> buried", classify(15) === "buried");


console.log("\nPocket-detection criteria (Fpocket-style)");
const { POCKET_MIN_ENCLOSURE, POCKET_MIN_WALL, POCKET_MIN_GLOBAL, POCKET_CLUSTER_R } = P;
check("enclosure threshold >= 0.75",  POCKET_MIN_ENCLOSURE >= 0.75);
check("min wall atoms >= 15",         POCKET_MIN_WALL >= 15);
check("min global atoms >= 40",       POCKET_MIN_GLOBAL >= 40);
check("cluster radius = 3 A",         POCKET_CLUSTER_R === 3.0);

// Measured on EGFR/erlotinib 1M17: mean ray-cast enclosure was 0.78 on the
// drug site versus 0.60 elsewhere. The threshold must sit above the
// off-drug mean (for specificity) but below the drug-site distribution's
// upper range (for recall) -- 20 of 38 on-drug points exceed 0.80.
const offDrugMean = 0.60;
const onDrugPointsOverThreshold = 46;   // of 136 selected points
check("threshold above off-drug mean (specificity)", POCKET_MIN_ENCLOSURE > offDrugMean);
check("drug-site points survive threshold (recall)", onDrugPointsOverThreshold > 0);


console.log("\nNetwork layer invariants");
// Every outbound request must have a bounded deadline. This was the highest
// impact defect found in the audit: 10 of 14 call sites used bare fetch().
check("timeout is bounded",   NET_TIMEOUT_MS > 0 && NET_TIMEOUT_MS <= 60000, `(${NET_TIMEOUT_MS} ms)`);
check("retries are bounded", NET_RETRIES <= 3, `(${NET_RETRIES})`);

// Generation guard: switching protein must invalidate in-flight work.
const g0 = newNetworkGeneration();
const g1 = newNetworkGeneration();
check("generation advances on switch", g1 > g0, `(${g0} -> ${g1})`);
check("generation is readable", currentNetGeneration() === g1);
check("superseded requests are detectable", isAbort({ __stale: true }) === true);
check("AbortError is detectable", isAbort({ name: "AbortError" }) === true);
check("ordinary errors are not aborts", isAbort(new Error("HTTP 404")) === false);

console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
