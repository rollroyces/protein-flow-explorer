/**
 * Unit tests for the annotation key module.
 *
 *   node test/annotations.test.js
 *
 * The point of this file is the cross-checks at the bottom. The app's
 * annotation key documents the colours each scheme uses; the colours live
 * in src/colors.js and src/variants.js. Those two things used to be
 * maintained separately, and nothing noticed when they disagreed — the
 * failure that let every scheme except "By chain" render the structure
 * black while the key went on advertising colours the app never produced.
 *
 * So every gradient in the key is compared here against the palette the
 * app really renders. Where the key deliberately shows fewer colours than
 * the palette holds, that gap is pinned by an explicit assertion rather
 * than left to drift.
 */
const fs = require("fs");
const path = require("path");
const vm = require("vm");

const ROOT = path.join(__dirname, "..");
const srcPath = path.join(ROOT, "src", "annotations.js");
const src = fs.readFileSync(srcPath, "utf8");
vm.runInThisContext(src, { filename: "src/annotations.js" });
const A = globalThis.PFE.annotations;

// The cross-checks need the palettes the app really renders. The module
// itself stays standalone on purpose; only the test loads all three.
vm.runInThisContext(fs.readFileSync(path.join(ROOT, "src", "colors.js"), "utf8"),
  { filename: "src/colors.js" });
vm.runInThisContext(fs.readFileSync(path.join(ROOT, "src", "variants.js"), "utf8"),
  { filename: "src/variants.js" });
const C = globalThis.PFE.colors;
const V = globalThis.PFE.variants;

let pass = 0, fail = 0;
function check(name, cond, detail = "") {
  if (cond) { pass++; console.log(`  PASS  ${name}${detail ? "  " + detail : ""}`); }
  else { fail++; console.log(`  FAIL  ${name}${detail ? "  " + detail : ""}`); }
}
const hex = (v) => "0x" + v.toString(16).padStart(6, "0");
const toRgb = (s) => parseInt(s.replace("#", ""), 16);
const stops = (id) => A.stopsFor(id).map(toRgb);
const show = (id) => A.stopsFor(id).join(" ");

// ------------------------------------------------------------ module shape
console.log("\nModule shape");
check("publishes on the PFE namespace", !!(globalThis.PFE && globalThis.PFE.annotations));
check("the exported object is frozen", Object.isFrozen(A));
check("the annotation list is frozen", Object.isFrozen(A.ANNOTATIONS));
check("every entry is frozen", A.ANNOTATIONS.every((a) => Object.isFrozen(a)));
check("every entry's stops are frozen", A.ANNOTATIONS.every((a) => Object.isFrozen(a.stops)));
check("it did not clobber the other modules",
      !!(globalThis.PFE.colors && globalThis.PFE.variants));
check("no fetch call", !/\bfetch\s*\(/.test(src));
check("no DOM access", !/\bdocument\s*\./.test(src) && !/\bwindow\s*\./.test(src));
check("no timers", !/\bsetTimeout\s*\(/.test(src) && !/\bsetInterval\s*\(/.test(src));
check("no XMLHttpRequest", !/XMLHttpRequest/.test(src));
check("not an ES module", !/^\s*(import|export)\s/m.test(src));
check("no require()", !/\brequire\s*\(/.test(src));
check("it is a classic IIFE", /\(function \(global\) \{/.test(src) &&
      /\}\)\(typeof globalThis !== "undefined" \? globalThis : this\);/.test(src));

// -------------------------------------------------------------- the key
console.log("\nThe key itself");
check("nine schemes", A.ANNOTATIONS.length === 9, `(${A.ANNOTATIONS.length})`);
check("the ids are unique",
      new Set(A.ANNOTATIONS.map((a) => a.id)).size === A.ANNOTATIONS.length);
check("every entry names itself",
      A.ANNOTATIONS.every((a) => typeof a.name === "string" && a.name.length > 0));
check("every entry says what it quantifies",
      A.ANNOTATIONS.every((a) => typeof a.quantity === "string" && a.quantity.length > 0));
check("every entry says what it is for",
      A.ANNOTATIONS.every((a) => typeof a.use === "string" && a.use.length > 0));
check("every entry has a gradient",
      A.ANNOTATIONS.every((a) => typeof a.grad === "string" && a.grad.length > 0));
check("the known scheme ids are all present",
      ["chain", "byDomain", "ss", "spectrum", "bFactor", "pLDDT", "conservation",
        "variant", "hydrophobicity"].every((id) => !!A.annotationFor(id)));

// ------------------------------------------------------------ availability
console.log("\nAvailability");
check("five flags are documented", A.AVAILABILITY_FLAGS.length === 5,
      `(${A.AVAILABILITY_FLAGS.join(", ")})`);
check("the flags are the documented five",
      ["crystal", "alphafold", "conservation", "domains", "variants"]
        .every((f) => A.AVAILABILITY_FLAGS.includes(f)));
check("every `needs` names a known flag",
      A.ANNOTATIONS.every((a) => !a.needs || A.AVAILABILITY_FLAGS.includes(a.needs)));
check("every flag is used by exactly one entry",
      A.AVAILABILITY_FLAGS.every((f) =>
        A.ANNOTATIONS.filter((a) => a.needs === f).length === 1));
const none = A.availability({});
check("nothing is conditional when no flag is set",
      A.ANNOTATIONS.filter((a) => !a.needs).every((a) => none[a.id] === true));
check("conditional entries are off when no flag is set",
      A.ANNOTATIONS.filter((a) => a.needs).every((a) => none[a.id] === false));
const all = A.availability({
  crystal: true, alphafold: true, conservation: true, domains: true, variants: true,
});
check("everything is on when every flag is set",
      A.ANNOTATIONS.every((a) => all[a.id] === true));
check("availability returns one key per annotation",
      Object.keys(all).length === A.ANNOTATIONS.length &&
      A.ANNOTATIONS.every((a) => Object.prototype.hasOwnProperty.call(all, a.id)));
check("availability ignores a missing flag object",
      Object.keys(A.availability()).length === A.ANNOTATIONS.length);
check("availability ignores a nonsense flag object",
      Object.keys(A.availability("nope")).length === A.ANNOTATIONS.length);
check("one flag switches exactly one entry",
      A.ANNOTATIONS.filter((a) => a.needs === "crystal").length === 1 &&
      A.ANNOTATIONS.filter((a) =>
        A.isAvailable(a.id, { crystal: true }) !== A.isAvailable(a.id, {})).length === 1);
check("isAvailable agrees with availability",
      A.ANNOTATIONS.every((a) =>
        A.isAvailable(a.id, { variants: true }) === A.availability({ variants: true })[a.id]));
check("an unknown scheme is not available", A.isAvailable("nope", all) === false);
check("an unknown flag name is reported",
      A.unknownFlags({ crytal: true, domains: true }).length === 1 &&
      A.unknownFlags({ crytal: true, domains: true })[0] === "crytal");
check("no unknown flag names when the caller is correct",
      A.unknownFlags({ crystal: true, alphafold: true, conservation: true,
        domains: true, variants: true }).length === 0);
check("a missing flag object has no unknown names", A.unknownFlags(null).length === 0);

// -------------------------------------------------------------- gradients
console.log("\nGradients");
check("the css string is built from the stops",
      A.gradientCss(["#2f7ed8", "#ef553b"]) === "linear-gradient(90deg,#2f7ed8,#ef553b)");
check("rgb integers are accepted",
      A.gradientCss([0x2f7ed8, 0xef553b]) === "linear-gradient(90deg,#2f7ed8,#ef553b)");
check("a single stop is not a gradient", A.gradientCss(["#2f7ed8"]) === "");
check("no stops is not a gradient", A.gradientCss([]) === "");
check("a non-array is not a gradient", A.gradientCss("nope") === "");
check("an invalid stop is not a gradient", A.gradientCss(["#2f7ed8", "#xyz"]) === "");
const parsed = A.parseGradient("linear-gradient(90deg,#ef5350,#f5c542,#9aa0a6)");
check("the gradient parses into its stops",
      parsed.valid && parsed.stops.length === 3 &&
      parsed.stops.join(",") === "#ef5350,#f5c542,#9aa0a6", `(${parsed.stops.join(" ")})`);
check("the gradient parses into its angle", parsed.angle === "90deg");
check("a one-stop gradient is invalid", A.parseGradient("linear-gradient(90deg,#ef5350)").valid === false);
check("a non-hex stop is invalid", A.parseGradient("linear-gradient(90deg,#ef5350,blue)").valid === false);
check("a three-digit hex is invalid", A.parseGradient("linear-gradient(90deg,#fff,#000)").valid === false);
check("a non-gradient string is invalid", A.parseGradient("#2f7ed8").valid === false);
check("a non-string is invalid", A.parseGradient(null).valid === false);
check("every entry's gradient parses back to its stops",
      A.ANNOTATIONS.every((a) => {
        const p = A.parseGradient(a.grad);
        return p.valid && p.stops.join(",") === a.stops.join(",");
      }));
check("every entry's stops are hex colours",
      A.ANNOTATIONS.every((a) => a.stops.every((s) => /^#[0-9a-f]{6}$/.test(s))));
check("stopsFor copies, so a caller cannot corrupt the key",
      (() => {
        const s = A.stopsFor("ss");
        s[0] = "#000000";
        return A.stopsFor("ss")[0] === "#ef5350";
      })());
check("stopsFor of an unknown scheme is empty", A.stopsFor("nope").length === 0);
check("annotationFor of an unknown scheme is null", A.annotationFor("nope") === null);
check("annotationFor ignores a non-string", A.annotationFor(42) === null);

// =========================================================================
// The cross-checks: the key against the palettes the app actually renders.
// =========================================================================

// ------------------------------------------------------------------- chain
console.log("\nCross-check: chain");
const chainStops = stops("chain");
check("the chain key is the start of the chain palette",
      chainStops.every((c, i) => c === C.STRUCT_CHAIN_PALETTE[i]),
      `(${show("chain")})`);
check("the chain key shows every palette colour",
      chainStops.length === C.STRUCT_CHAIN_PALETTE.length,
      `(key ${chainStops.length} of ${C.STRUCT_CHAIN_PALETTE.length})`);
check("the chain key invents no colour of its own",
      chainStops.every((c) => C.STRUCT_CHAIN_PALETTE.includes(c)));

// --------------------------------------------------- secondary structure
console.log("\nCross-check: secondary structure");
const ssStops = stops("ss");
check("the ss key is helix, strand, coil",
      ssStops.length === 3 &&
      ssStops[0] === C.ssColor("h") &&
      ssStops[1] === C.ssColor("s") &&
      ssStops[2] === C.ssColor("c"), `(${show("ss")})`);
check("the ss key invents no colour of its own",
      ssStops.every((c) => [C.ssColor("h"), C.ssColor("s"), C.ssColor("c")].includes(c)));

// --------------------------------------------------------------- spectrum
console.log("\nCross-check: spectrum");
check("the spectrum key is the spectrum ramp",
      stops("spectrum").join(",") === C.SPECTRUM_STOPS.join(","), `(${show("spectrum")})`);
check("the spectrum key matches the ramp endpoints",
      stops("spectrum")[0] === C.SPECTRUM_STOPS[0] &&
      stops("spectrum")[stops("spectrum").length - 1] === C.SPECTRUM_STOPS[C.SPECTRUM_STOPS.length - 1]);

// --------------------------------------------------------------- B-factor
console.log("\nCross-check: B-factor");
check("the bFactor key is the b-factor ramp",
      stops("bFactor").join(",") === C.BFACTOR_STOPS.join(","), `(${show("bFactor")})`);
check("the bFactor key matches the ramp endpoints",
      stops("bFactor")[0] === C.BFACTOR_STOPS[0] &&
      stops("bFactor")[stops("bFactor").length - 1] === C.BFACTOR_STOPS[C.BFACTOR_STOPS.length - 1]);

// ----------------------------------------------------------- conservation
console.log("\nCross-check: conservation");
check("the conservation key is the conservation ramp",
      stops("conservation").join(",") === C.CONSERVATION_STOPS.join(","),
      `(${show("conservation")})`);
check("the conservation key matches the ramp endpoints",
      stops("conservation")[0] === C.CONSERVATION_STOPS[0] &&
      stops("conservation")[stops("conservation").length - 1] ===
      C.CONSERVATION_STOPS[C.CONSERVATION_STOPS.length - 1]);

// --------------------------------------------------------- hydrophobicity
console.log("\nCross-check: hydrophobicity");
check("the hydrophobicity key is the hydrophobicity ramp",
      stops("hydrophobicity").join(",") === C.HYDROPHOBICITY_STOPS.join(","),
      `(${show("hydrophobicity")})`);
check("the hydrophobicity key matches the ramp endpoints",
      stops("hydrophobicity")[0] === C.HYDROPHOBICITY_STOPS[0] &&
      stops("hydrophobicity")[stops("hydrophobicity").length - 1] ===
      C.HYDROPHOBICITY_STOPS[C.HYDROPHOBICITY_STOPS.length - 1]);

// ----------------------------------------------------------------- pLDDT
console.log("\nCross-check: pLDDT");
const plddtStops = stops("pLDDT");
check("the pLDDT key is the four confidence bands",
      plddtStops[0] === C.plddtColor(20) &&
      plddtStops[1] === C.plddtColor(60) &&
      plddtStops[2] === C.plddtColor(85) &&
      plddtStops[3] === C.plddtColor(95), `(${show("pLDDT")})`);
check("the pLDDT key is exactly four stops", plddtStops.length === 4);
check("the pLDDT key covers every band, not just the probes",
      [10, 49, 50, 69, 70, 89, 90, 100].every((p) => plddtStops.includes(C.plddtColor(p))));
check("the pLDDT key runs low confidence to high",
      plddtStops[0] === C.plddtColor(20) && plddtStops[3] === C.plddtColor(95));

// --------------------------------------------------------------- variant
console.log("\nCross-check: variant significance");
const variantStops = stops("variant");
// Every ClinVar class the key documents, in the order it documents them.
// The key used to document only five of the ten the scheme renders, which
// left drug_response -- the class T790M and L858R actually are -- unfindable
// in the key. The gap is closed, and this assertion now prevents it opening
// again.
const SHOWN_CLASSES = V.SIGNIFICANCE_ORDER.slice();
check("the variant key is the documented significance classes",
      variantStops.length === SHOWN_CLASSES.length &&
      SHOWN_CLASSES.every((c, i) => variantStops[i] === V.SIGNIFICANCE_COLORS[c]),
      `(${show("variant")})`);
check("the variant key shows the classes in severity order",
      SHOWN_CLASSES.join(",") ===
      V.SIGNIFICANCE_ORDER.filter((c) => SHOWN_CLASSES.includes(c)).join(","));
check("the variant key invents no colour of its own",
      variantStops.every((c) => Object.values(V.SIGNIFICANCE_COLORS).includes(c)));
check("the key documents every significance class the scheme renders",
      V.SIGNIFICANCE_ORDER.every((c) => variantStops.includes(V.SIGNIFICANCE_COLORS[c])),
      `(${V.SIGNIFICANCE_ORDER.filter((c) => !variantStops.includes(V.SIGNIFICANCE_COLORS[c])).join(", ") || "none missing"})`);
check("no significance class is undocumented",
      V.SIGNIFICANCE_ORDER.filter((c) => !SHOWN_CLASSES.includes(c)).length === 0);
check("every significance class still has a colour",
      V.SIGNIFICANCE_ORDER.every((c) => typeof V.SIGNIFICANCE_COLORS[c] === "number"));

// -------------------------------------------------------------- by domain
console.log("\nCross-check: by domain");
// The domain scheme renders the per-domain colours curated in
// data/proteins.js, so the key must show those colours and no others.
vm.runInThisContext(fs.readFileSync(path.join(ROOT, "data", "proteins.js"), "utf8"),
  { filename: "data/proteins.js" });
const PROTEINS = vm.runInThisContext("PROTEINS");
const domainColors = new Set();
for (const p of PROTEINS) {
  for (const d of p.domains || []) domainColors.add(d.color.toLowerCase());
}
check("the curated dataset has domain colours", domainColors.size > 0,
      `(${[...domainColors].join(" ")})`);
check("the byDomain key shows the curated domain colours",
      stops("byDomain").every((c) => domainColors.has("#" + hex(c).slice(2))) &&
      [...domainColors].every((c) => stops("byDomain").includes(toRgb(c))),
      `(${show("byDomain")})`);
check("the byDomain key invents no colour of its own",
      stops("byDomain").every((c) => domainColors.has("#" + hex(c).slice(2))));

// ------------------------------------------ how app.js obtains the key
console.log("\nHow app.js obtains the key");
// app.js must not keep its own copy of the key: a second copy is exactly how
// the key drifted away from the palettes in the first place.
const appSrc = fs.readFileSync(path.join(ROOT, "app.js"), "utf8");
const hasLocalCopy = /const ANNOTATIONS\s*=\s*\[/.test(appSrc);
check("app.js holds no local copy of the key", !hasLocalCopy,
      hasLocalCopy ? "(app.js still declares its own array)" : "");
check("app.js reads the key from the module",
      /PFE\.annotations\.ANNOTATIONS/.test(appSrc));

// ------------------------------------------------------------- validation
console.log("\nValidation of the real key");
const realValidation = A.validateAll();
check("the shipped key is valid", realValidation.ok,
      realValidation.ok ? "" : `(${realValidation.errors.join("; ")})`);
check("every entry validates on its own",
      A.ANNOTATIONS.every((a) => A.validate(a).ok));
check("assertValid accepts the shipped key", (() => {
  try { A.assertValid(); return true; } catch (e) { return false; }
})());

console.log("\nValidation rejects malformed entries");
const base = A.annotationFor("ss");
const clone = (patch) => ({ ...base, ...patch });
const bad = [
  ["a missing id", clone({ id: "" })],
  ["a missing name", clone({ name: undefined })],
  ["a missing quantity", clone({ quantity: "  " })],
  ["a missing use", clone({ use: undefined })],
  ["a one-stop gradient", clone({ grad: "linear-gradient(90deg,#ef5350)" })],
  ["a non-hex stop", clone({ grad: "linear-gradient(90deg,#ef5350,blue)" })],
  ["a three-digit hex", clone({ grad: "linear-gradient(90deg,#fff,#000)" })],
  ["a gradient that is not one", clone({ grad: "#ef5350" })],
  ["a missing gradient", clone({ grad: undefined })],
  ["an unknown needs flag", clone({ needs: "crytal" })],
  ["stops that disagree with the gradient", clone({ stops: ["#111111", "#222222"] })],
  ["stops that are not colours", clone({ stops: ["#111111", "blue"] })],
  ["stops that are too few", clone({ stops: ["#111111"] })],
  ["an entry that is not an object", null],
];
for (const [name, entry] of bad) {
  const res = A.validate(entry);
  check(`${name} is rejected`, res.ok === false && res.errors.length > 0,
        res.ok ? "(accepted)" : `(${res.errors[0]})`);
}
check("a malformed entry fails assertValid", (() => {
  try { A.assertValid([clone({ grad: "linear-gradient(90deg,#ef5350)" })]); return false; }
  catch (e) { return /malformed/.test(e.message); }
})());
check("a duplicate id is reported", (() => {
  const res = A.validateAll([base, { ...base }]);
  return !res.ok && res.errors.some((e) => /duplicate id/.test(e));
})());
check("validateAll defaults to the shipped key", A.validateAll().ok === true);
check("validateAll accepts an empty list", A.validateAll([]).ok === true);

console.log("\nDocumentation");
// Every exported function must carry a jsdoc block with @param and
// @returns, immediately above it, so the convention cannot rot.
const EXPORTED_FUNCTIONS = [
  "gradientCss", "parseGradient", "annotationFor", "stopsFor", "availability",
  "isAvailable", "unknownFlags", "validate", "validateAll", "assertValid",
];
for (const fn of EXPORTED_FUNCTIONS) {
  const at = src.indexOf(`function ${fn}(`);
  const before = at === -1 ? "" : src.slice(0, at);
  const doc = before.match(/\/\*\*[\s\S]*?\*\/\s*$/);
  check(`${fn} is exported`, typeof A[fn] === "function");
  check(`${fn} has jsdoc with @param and @returns`,
        !!doc && /@param/.test(doc[0]) && /@returns/.test(doc[0]));
}
check("the api exports nothing undocumented",
      Object.keys(A).every((k) => typeof A[k] !== "function" || EXPORTED_FUNCTIONS.includes(k)));
check("the api exports the key and the flags",
      Array.isArray(A.ANNOTATIONS) && Array.isArray(A.AVAILABILITY_FLAGS));

console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
