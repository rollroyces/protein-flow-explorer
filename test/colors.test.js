/**
 * Unit tests for the colour module.
 *
 *   node test/colors.test.js
 *
 * Regression cover for the colouring audit. The central defect was that the
 * app put colour-scheme NAMES in 3Dmol's `color` field, which 3Dmol treats as
 * a literal colour value: every scheme name fell through to htmlColors lookup,
 * missed, and became 0x000000 — so the whole structure rendered black.
 *
 * These tests assert the colour maths directly, and assert that the app's
 * annotation key and the colours actually rendered agree.
 */
const fs = require("fs");
const path = require("path");
const vm = require("vm");

const ROOT = path.join(__dirname, "..");
vm.runInThisContext(fs.readFileSync(path.join(ROOT, "src", "colors.js"), "utf8"),
  { filename: "src/colors.js" });
const C = globalThis.PFE.colors;

let pass = 0, fail = 0;
function check(name, cond, detail = "") {
  if (cond) { pass++; console.log(`  PASS  ${name}${detail ? "  " + detail : ""}`); }
  else { fail++; console.log(`  FAIL  ${name}${detail ? "  " + detail : ""}`); }
}
const hex = (v) => "0x" + v.toString(16).padStart(6, "0");

// --------------------------------------------------------------- palette
console.log("\nChain palette");
check("ten colours", C.STRUCT_CHAIN_PALETTE.length === 10);
check("all distinct", new Set(C.STRUCT_CHAIN_PALETTE).size === 10);
const lum = C.STRUCT_CHAIN_PALETTE.map((v) => {
  const r = (v >> 16) & 255, g = (v >> 8) & 255, b = v & 255;
  return 0.2126 * r + 0.7152 * g + 0.0722 * b;
});
check("all visible on white", lum.every((l) => l > 40), `(min lum ${Math.min(...lum).toFixed(0)})`);

// ------------------------------------------------------------------- ramp
console.log("\nRamp interpolation");
check("t=0 is the first stop", C.rampColor([0x112233, 0x445566], 0) === 0x112233);
check("t=1 is the last stop", C.rampColor([0x112233, 0x445566], 1) === 0x445566);
check("t=0.5 is the midpoint", C.rampColor([0x000000, 0xffffff], 0.5) === 0x808080 ||
      Math.abs(C.rampColor([0x000000, 0xffffff], 0.5) - 0x808080) <= 1);
check("out-of-range clamps low", C.rampColor([0x112233, 0x445566], -5) === 0x112233);
check("out-of-range clamps high", C.rampColor([0x112233, 0x445566], 5) === 0x445566);
check("NaN yields the last stop", C.rampColor([0x112233, 0x445566], NaN) === 0x445566);
check("single stop is returned as-is", C.rampColor([0xabcdef], 0.5) === 0xabcdef);

console.log("\nNormalisation");
check("value at the minimum -> 0", C.normalize(10, 10, 20) === 0);
check("value at the maximum -> 1", C.normalize(20, 10, 20) === 1);
check("value mid-range -> 0.5", C.normalize(15, 10, 20) === 0.5);
check("degenerate range -> 0", C.normalize(5, 10, 10) === 0);
check("non-numeric -> 0", C.normalize("x", 10, 20) === 0);

// ------------------------------------------------- secondary structure
console.log("\nSecondary structure");
check("helix -> red", C.ssColor("h") === 0xef5350, `(${hex(C.ssColor("h"))})`);
check("strand -> yellow", C.ssColor("s") === 0xf5c542, `(${hex(C.ssColor("s"))})`);
check("coil -> grey", C.ssColor("c") === 0x9aa0a6, `(${hex(C.ssColor("c"))})`);
check("uppercase is accepted", C.ssColor("H") === 0xef5350);
check("cartoon arrow start is a strand", C.ssColor("arrow start") === 0xf5c542);
check("cartoon arrow end is a strand", C.ssColor("arrow end") === 0xf5c542);
check("undefined falls back to grey", C.ssColor(undefined) === 0x9aa0a6);
check("unknown falls back to grey", C.ssColor("zzz") === 0x9aa0a6);

// -------------------------------------------------------------- spectrum
console.log("\nSpectrum (N->C)");
check("N-terminus is blue", C.spectrumColor(1, 1, 400) === C.SPECTRUM_STOPS[0]);
check("C-terminus is red", C.spectrumColor(400, 1, 400) === C.SPECTRUM_STOPS[3]);
check("midpoint is between the stops",
      C.spectrumColor(200, 1, 400) !== C.spectrumColor(1, 1, 400));
check("direction is N->C (blue to red)",
      ((C.spectrumColor(400, 1, 400) >> 16) & 0xff) > ((C.spectrumColor(1, 1, 400) >> 16) & 0xff),
      `(${hex(C.spectrumColor(1, 1, 400))} -> ${hex(C.spectrumColor(400, 1, 400))})`);

// -------------------------------------------------------------- B-factor
console.log("\nB-factor");
check("low B is blue", C.bfactorColor(10, 10, 90) === C.BFACTOR_STOPS[0]);
check("high B is red", C.bfactorColor(90, 10, 90) === C.BFACTOR_STOPS[3]);
check("mid-range is in between",
      C.bfactorColor(50, 10, 90) !== C.bfactorColor(10, 10, 90));

// -------------------------------------------------------- hydrophobicity
console.log("\nHydrophobicity (Kyte-Doolittle)");
// The scale is keyed by three-letter codes because 3Dmol's atom.resn is
// three-letter; a single-letter key silently yields the neutral colour.
check("Ile is the most hydrophobic", C.hydrophobicityColor("ILE") ===
      C.HYDROPHOBICITY_STOPS[2], `(${hex(C.hydrophobicityColor("ILE"))})`);
check("Arg is strongly hydrophilic", C.hydrophobicityColor("ARG") !== C.HYDROPHOBICITY_STOPS[2]);
check("Arg is bluer than Ile",
      ((C.hydrophobicityColor("ARG") >> 16) & 0xff) < ((C.hydrophobicityColor("ILE") >> 16) & 0xff),
      `(${hex(C.hydrophobicityColor("ARG"))} vs ${hex(C.hydrophobicityColor("ILE"))})`);
check("Leu is more hydrophobic than Ser",
      ((C.hydrophobicityColor("LEU") >> 16) & 0xff) > ((C.hydrophobicityColor("SER") >> 16) & 0xff));
check("Trp is ambiguous, not strongly hydrophobic",
      C.hydrophobicityColor("TRP") !== C.HYDROPHOBICITY_STOPS[2]);
check("lower case is accepted", C.hydrophobicityColor("ala") === C.hydrophobicityColor("ALA"));
check("unknown residue is neutral", C.hydrophobicityColor("XXX") === C.HYDROPHOBICITY_STOPS[1]);
check("undefined is neutral", C.hydrophobicityColor(undefined) === C.HYDROPHOBICITY_STOPS[1]);
check("the scale covers the documented range",
      Math.min(...Object.values(C.KD_SCALE)) === -4.5 &&
      Math.max(...Object.values(C.KD_SCALE)) === 4.5);
check("every standard residue has a score",
      "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR"
        .split(" ").every((r) => C.KD_SCALE[r] !== undefined));

// ------------------------------------------------------------------ pLDDT
console.log("\npLDDT (AlphaFold convention)");
check("95 -> deep blue", C.plddtColor(95) === 0x0053d6);
check("90 -> deep blue (boundary)", C.plddtColor(90) === 0x0053d6);
check("89 -> cyan", C.plddtColor(89) === 0x65cbf3);
check("70 -> cyan (boundary)", C.plddtColor(70) === 0x65cbf3);
check("69 -> yellow", C.plddtColor(69) === 0xffdb13);
check("50 -> yellow (boundary)", C.plddtColor(50) === 0xffdb13);
check("49 -> orange", C.plddtColor(49) === 0xff7d45);
check("0 -> orange", C.plddtColor(0) === 0xff7d45);
check("non-numeric is safe", typeof C.plddtColor(undefined) === "number");

// ------------------------------------------------------------ conservation
console.log("\nConservation");
check("0 (variable) -> blue", C.conservationColor(0) === C.CONSERVATION_STOPS[0]);
check("1 (invariant) -> red", C.conservationColor(1) === C.CONSERVATION_STOPS[2]);
check("out-of-range clamps", C.conservationColor(5) === C.CONSERVATION_STOPS[2]);
check("negative clamps to the low stop", C.conservationColor(-2) === C.CONSERVATION_STOPS[0]);
check("non-numeric -> grey", C.conservationColor(undefined) === 0xcccccc);

// ------------------------------------- the annotation key must agree
console.log("\nAnnotation key agreement");
const appSrc = fs.readFileSync(path.join(ROOT, "app.js"), "utf8");
const annotStart = appSrc.indexOf("const ANNOTATIONS");
const annotEnd = appSrc.indexOf("];", annotStart);
const annotBlock = appSrc.slice(annotStart, annotEnd);
// Index the gradients by annotation id, not by position, so that reordering
// the key does not silently invalidate these checks.
const grads = {};
for (const m of annotBlock.matchAll(/id:\s*"([^"]+)"[\s\S]{0,400}?grad:\s*"linear-gradient\(90deg,([^)]+)\)"/g)) {
  grads[m[1]] = m[2].split(",").map((s) => parseInt(s.replace("#", ""), 16));
}
check("annotation key has 8 entries", Object.keys(grads).length === 8,
      `(${Object.keys(grads).length})`);
// The key must not advertise a colour the app cannot produce. Compare the
// chain key against the real palette.
check("chain key matches the palette",
      grads.chain.every((c, i) => c === C.STRUCT_CHAIN_PALETTE[i]),
      `(${grads.chain.map(hex).join(" ")})`);
// ss key: helix, strand, coil
check("ss key matches the rendered colours",
      grads.ss[0] === C.ssColor("h") && grads.ss[1] === C.ssColor("s") && grads.ss[2] === C.ssColor("c"),
      `(${grads.ss.map(hex).join(" ")})`);
// spectrum key: first and last stop
check("spectrum key matches the ramp endpoints",
      grads.spectrum[0] === C.SPECTRUM_STOPS[0] && grads.spectrum[3] === C.SPECTRUM_STOPS[3]);
// bFactor key
check("bFactor key matches the ramp endpoints",
      grads.bFactor[0] === C.BFACTOR_STOPS[0] && grads.bFactor[3] === C.BFACTOR_STOPS[3]);
// conservation key
check("conservation key matches the ramp endpoints",
      grads.conservation[0] === C.CONSERVATION_STOPS[0] && grads.conservation[2] === C.CONSERVATION_STOPS[2]);
// hydrophobicity key
check("hydrophobicity key matches the ramp endpoints",
      grads.hydrophobicity[0] === C.HYDROPHOBICITY_STOPS[0] &&
      grads.hydrophobicity[2] === C.HYDROPHOBICITY_STOPS[2]);
// domain key: the two curated domain palettes (blue / grey / red)
check("domain key matches the curated domain colours",
      grads.byDomain[0] === 0x2f7ed8 && grads.byDomain[1] === 0x9aa0a6 &&
      grads.byDomain[2] === 0xcc1f1f,
      `(${grads.byDomain.map(hex).join(" ")})`);

// --------------------------------------- the curated domain data
console.log("\nCurated domain data");
vm.runInThisContext(fs.readFileSync(path.join(ROOT, "data/proteins.js"), "utf8"),
  { filename: "proteins.js" });
// PROTEINS is a top-level const, i.e. a global lexical binding rather than a
// property of globalThis, so it has to be read out of the context.
const PROTEINS = vm.runInThisContext("PROTEINS");
const withDomains = PROTEINS.filter((p) => p.domains && p.domains.length);
check("at least two proteins carry domain data", withDomains.length >= 2,
      `(${withDomains.map((p) => p.id).join(", ")})`);
for (const p of withDomains) {
  const doms = p.domains;
  check(`${p.id}: every domain has a name`, doms.every((d) => !!d.name));
  check(`${p.id}: every domain colour is a hex literal`,
        doms.every((d) => /^#[0-9a-f]{6}$/i.test(d.color)),
        `(${doms.map((d) => d.color).join(", ")})`);
  check(`${p.id}: every domain has at least one range`,
        doms.every((d) => d.ranges && d.ranges.length > 0));
  // Ranges must be ordered, non-overlapping and ascending within a domain.
  check(`${p.id}: ranges within a domain ascend`,
        doms.every((d) => d.ranges.every((r) => r.start <= r.end)));
  // Domains must not overlap each other.
  const all = doms.flatMap((d) => d.ranges.map((r) => [r.start, r.end, d.name]));
  let overlap = false;
  for (let i = 0; i < all.length; i++) {
    for (let j = i + 1; j < all.length; j++) {
      if (all[i][0] <= all[j][1] && all[j][0] <= all[i][1]) overlap = true;
    }
  }
  check(`${p.id}: domains do not overlap`, !overlap);
  check(`${p.id}: domain colours are distinct`,
        new Set(doms.map((d) => d.color.toLowerCase())).size === doms.length);
}

// The guided tours advertise specific colours for specific domains. If the
// data and the tour text disagree, the tour is lying to the user.
const proteinsSrc = fs.readFileSync(path.join(ROOT, "data/proteins.js"), "utf8");
const akEntry = proteinsSrc.slice(proteinsSrc.indexOf('id: "adenylate-kinase"'));
const akDomains = akEntry.slice(akEntry.indexOf("domains:"), akEntry.indexOf("highlights:"));
check("AK tour's CORE is grey", /CORE"\s*,\s*color:\s*"#9aa0a6"/.test(akDomains));
check("AK tour's NMP is blue", /NMP"\s*,\s*color:\s*"#2f7ed8"/.test(akDomains));
check("AK tour's LID is red", /LID"\s*,\s*color:\s*"#cc1f1f"/.test(akDomains));

// ------------------------------- the 3Dmol contract that caused the bug
console.log("\n3Dmol colour-field contract (the regression)");
// 3Dmol reads the scheme from `colorscheme`; `color` must be a literal
// colour. Assert the app no longer puts a scheme name in `color`.
const resolveBlock = appSrc.slice(appSrc.indexOf("function resolveColor"),
  appSrc.indexOf("function getStyleSpec"));
const literalColours = [...resolveBlock.matchAll(/color:\s*"([^"]+)"/g)].map((m) => m[1]);
check("no scheme name is passed as a literal colour",
      literalColours.every((c) => c.startsWith("#")),
      `(${literalColours.join(", ")})`);
const schemeNames = [...resolveBlock.matchAll(/colorscheme:\s*"([^"]+)"/g)].map((m) => m[1]);
check("3Dmol schemes are passed via `colorscheme`", schemeNames.length === 2,
      `(${schemeNames.join(", ")})`);
check("both schemes are real 3Dmol schemes",
      schemeNames.every((s) => ["Jmol", "amino", "ssJmol", "chain", "rasmol", "default"].includes(s)));

// The same contract, applied to EVERY colour literal in the app rather than
// just resolveColor(). A scheme name anywhere in a `color:` field renders
// that representation black.
console.log("\nEvery colour literal in the app is a real colour");
const HTML_COLOURS = new Set([
  "yellow", "magenta", "cyan", "white", "black", "red", "green", "blue",
  "grey", "gray", "orange", "purple", "pink", "silver", "gold", "lime",
  "maroon", "navy", "teal", "olive", "aqua", "fuchsia", "brown", "cyan",
]);
const SCHEME_NAMES = new Set([
  "chain", "ss", "spectrum", "byElement", "residue", "hydrophobicity",
  "bFactor", "byDomain", "pLDDT", "conservation", "uniform",
  "ssJmol", "ssPyMol", "Jmol", "amino", "rasmol", "default", "chainHetatm",
]);
const allLiterals = [...appSrc.matchAll(/\bcolor:\s*"([^"]+)"/g)].map((m) => m[1]);
const schemeAsLiteral = allLiterals.filter((c) => SCHEME_NAMES.has(c));
check("no scheme name appears in any colour field", schemeAsLiteral.length === 0,
      schemeAsLiteral.length ? `(${schemeAsLiteral.join(", ")})` : "");
const notAColour = allLiterals.filter((c) => !c.startsWith("#") && !HTML_COLOURS.has(c.toLowerCase()));
check("every string colour is an html colour or a hex literal", notAColour.length === 0,
      notAColour.length ? `(${notAColour.join(", ")})` : "");
check("the app does use colour fields", allLiterals.length > 0, `(${allLiterals.length})`);

// The colour-scheme names the app advertises must all be real 3Dmol schemes
// or app-owned schemes, never typos that would silently fall through.
const KNOWN = new Set(["chain", "ss", "spectrum", "byElement", "residue",
  "hydrophobicity", "bFactor", "byDomain", "pLDDT", "conservation", "uniform"]);
const advertised = [...resolveBlock.matchAll(/case "([A-Za-z]+)":/g)].map((m) => m[1]);
check("every colour case is a known scheme",
      advertised.length === 11 && advertised.every((a) => KNOWN.has(a)),
      `(${advertised.join(", ")})`);

console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
