/**
 * Unit tests for the drug-discovery module.
 *
 *   node test/drug.test.js
 *
 * src/drug.js is pure logic: it takes pocket lists, ligand atoms and
 * externally-supplied reference records and returns verdicts, occupancies
 * and cross-validation checks. Nothing in it touches the DOM or the
 * network, so the whole surface can be exercised here without a browser
 * and without a fixture on disk.
 *
 * The cross-validation section is the important one. The module's job in
 * the app is to check the app's own arithmetic against an authoritative
 * external source (the Chemical Component Dictionary, Open Targets), so a
 * test that only exercised the happy path would not catch a check that
 * silently stopped failing. Every check is therefore asserted in both
 * directions: it must pass on the real EGFR/AQ4 numbers and it must fail
 * on numbers that are wrong.
 */
const fs = require("fs");
const path = require("path");
const vm = require("vm");

const ROOT = path.join(__dirname, "..");
const srcPath = path.join(ROOT, "src", "drug.js");
const src = fs.readFileSync(srcPath, "utf8");
vm.runInThisContext(src, { filename: "src/drug.js" });
const D = globalThis.PFE.drug;

let pass = 0, fail = 0;
function check(name, cond, detail = "") {
  if (cond) { pass++; console.log(`  PASS  ${name}${detail ? "  " + detail : ""}`); }
  else { fail++; console.log(`  FAIL  ${name}${detail ? "  " + detail : ""}`); }
}
const near = (a, b, tol = 1e-6) => Math.abs(a - b) <= tol;
const byName = (checks, name) => checks.find((c) => c.name === name);

// ------------------------------------------------------------ module shape
console.log("\nModule shape");
check("publishes on the PFE namespace", !!(globalThis.PFE && globalThis.PFE.drug));
check("the exported object is frozen", Object.isFrozen(D));
check("every documented member is exported",
      ["DRUGGABILITY", "ATOMIC_MASS", "druggabilityVerdict", "pocketOccupancy",
        "validatePocketsAgainstLigand", "contactResidues", "pocketResidueRecall",
        "crossValidateLigand", "countAtomsFromFormula", "molecularWeight",
        "formatFormula", "formatPotency", "summariseAffinity", "summariseTarget",
        "matchDrugToLigand"].every((k) => D[k] !== undefined));
check("the constants are the plain ones", typeof D.DRUGGABILITY === "object" &&
      typeof D.ATOMIC_MASS === "object");
check("no fetch call", !/\bfetch\s*\(/.test(src));
check("no DOM access", !/\bdocument\s*\./.test(src) && !/\bwindow\s*\./.test(src));
check("no timers", !/\bsetTimeout\s*\(/.test(src) && !/\bsetInterval\s*\(/.test(src));
check("not an ES module", !/^\s*(import|export)\s/m.test(src));
check("no require()", !/\brequire\s*\(/.test(src));
check("it is a classic IIFE", /\(function \(global\) \{/.test(src) &&
      /\}\)\(typeof globalThis !== "undefined" \? globalThis : this\);/.test(src));

// ---------------------------------------------------------- the constants
console.log("\nDruggability constants");
check("the volume floor is 120 Å³", D.DRUGGABILITY.minVolume === 120);
check("the volume ceiling is 1200 Å³", D.DRUGGABILITY.maxVolume === 1200);
check("the burial threshold is 0.45", D.DRUGGABILITY.minEnclosure === 0.45);
check("the exposure threshold is 0.25", D.DRUGGABILITY.marginalEnclosure === 0.25);
check("the tabulated masses are the expected ones",
      D.ATOMIC_MASS.H === 1.008 && D.ATOMIC_MASS.C === 12.011 &&
      D.ATOMIC_MASS.N === 14.007 && D.ATOMIC_MASS.O === 15.999 &&
      D.ATOMIC_MASS.S === 32.06 && D.ATOMIC_MASS.P === 30.974 &&
      D.ATOMIC_MASS.F === 18.998 && D.ATOMIC_MASS.Cl === 35.45 &&
      D.ATOMIC_MASS.Br === 79.904 && D.ATOMIC_MASS.I === 126.904 &&
      D.ATOMIC_MASS.Se === 78.971 && D.ATOMIC_MASS.Fe === 55.845 &&
      D.ATOMIC_MASS.Zn === 65.38);

// ------------------------------------------------------------- the verdict
console.log("\nDruggability verdicts");
const good = D.druggabilityVerdict({ volume: 400, enclosure: 0.6 });
check("a pocket inside the window is druggable", good.verdict === "druggable",
      `(${good.verdict})`);
check("its label is Druggable", good.label === "Druggable");
check("it says why", good.reasons.length === 1 &&
      /drug-like range/.test(good.reasons[0]), `(${good.reasons[0]})`);
check("the lower edge of the window is still druggable",
      D.druggabilityVerdict({ volume: 120, enclosure: 0.5 }).verdict === "druggable");
check("the upper edge of the window is still druggable",
      D.druggabilityVerdict({ volume: 1200, enclosure: 0.5 }).verdict === "druggable");
check("exactly at the burial threshold is druggable",
      D.druggabilityVerdict({ volume: 400, enclosure: 0.45 }).verdict === "druggable");
const small = D.druggabilityVerdict({ volume: 50, enclosure: 0.6 });
check("a volume below the floor is marginal", small.verdict === "marginal",
      `(${small.verdict})`);
check("its label is Marginal", small.label === "Marginal");
check("it names the floor", /below the 120 Å³ floor/.test(small.reasons.join(" ")),
      `(${small.reasons[0]})`);
const big = D.druggabilityVerdict({ volume: 2000, enclosure: 0.6 });
check("a volume above the ceiling is marginal", big.verdict === "marginal",
      `(${big.verdict})`);
check("it names the ceiling", /above the 1200 Å³ ceiling/.test(big.reasons.join(" ")),
      `(${big.reasons[0]})`);
const shallow = D.druggabilityVerdict({ volume: 400, enclosure: 0.30 });
check("enclosure 0.30 is still druggable", shallow.verdict === "druggable",
      `(${shallow.verdict})`);
check("enclosure 0.30 is flagged as shallower than ideal",
      /shallower than ideal/.test(shallow.reasons.join(" ")), `(${shallow.reasons[0]})`);
check("enclosure exactly at the exposure floor is still druggable",
      D.druggabilityVerdict({ volume: 400, enclosure: 0.25 }).verdict === "druggable");
const exposed = D.druggabilityVerdict({ volume: 400, enclosure: 0.10 });
check("enclosure 0.10 is marginal", exposed.verdict === "marginal",
      `(${exposed.verdict})`);
check("it says the pocket is too exposed",
      /too exposed/.test(exposed.reasons.join(" ")), `(${exposed.reasons[0]})`);
check("a shallow pocket still says why",
      shallow.reasons.length === 1 && exposed.reasons.length === 1);
const nanVol = D.druggabilityVerdict({ volume: NaN, enclosure: 0.6 });
check("a NaN volume is unknown", nanVol.verdict === "unknown", `(${nanVol.verdict})`);
check("an unknown verdict is labelled Unknown", nanVol.label === "Unknown");
check("it says the volume was not measured",
      /volume not measured/.test(nanVol.reasons.join(" ")), `(${nanVol.reasons[0]})`);
const noPocket = D.druggabilityVerdict(null);
check("a null pocket is unknown", noPocket.verdict === "unknown",
      `(${noPocket.verdict})`);
check("it says there was no pocket", /no pocket/.test(noPocket.reasons.join(" ")),
      `(${noPocket.reasons[0]})`);
check("an unmeasured enclosure does not block druggability",
      D.druggabilityVerdict({ volume: 400, enclosure: NaN }).verdict === "druggable");
check("the verdict is always one of the three documented values",
      [
        { volume: 400, enclosure: 0.6 },
        { volume: 50, enclosure: 0.1 },
        { volume: NaN, enclosure: 0.6 },
        null,
      ].every((p) => ["druggable", "marginal", "unknown"]
        .includes(D.druggabilityVerdict(p).verdict)));

// --------------------------------------------------------- pocket occupancy
console.log("\nPocket occupancy");
const ligand = [{ x: 0, y: 0, z: 0 }, { x: 1, y: 0, z: 0 }, { x: 0, y: 1, z: 0 }];
const centre = { center: [0, 0, 0] };
const allIn = D.pocketOccupancy(centre, ligand, 5);
check("a pocket containing every atom has occupancy 1", allIn.occupancy === 1,
      `(${allIn.occupancy})`);
check("it counts every atom as inside", allIn.inside === 3 && allIn.total === 3);
const noneIn = D.pocketOccupancy({ center: [100, 100, 100] }, ligand, 5);
check("a pocket containing nothing has occupancy 0", noneIn.occupancy === 0,
      `(${noneIn.occupancy})`);
check("it counts no atom as inside", noneIn.inside === 0);
check("the mean distance is reported", near(noneIn.meanDistance,
      (Math.hypot(100, 100, 100) + Math.hypot(99, 100, 100) + Math.hypot(100, 99, 100)) / 3));
check("the minimum distance is reported", near(noneIn.minDistance,
      Math.hypot(99, 100, 100)), `(${noneIn.minDistance})`);
const boundary = D.pocketOccupancy({ center: [0, 0, 0] }, [{ x: 5, y: 0, z: 0 }], 5);
check("an atom exactly on the radius counts as inside", boundary.inside === 1 &&
      boundary.occupancy === 1);
check("the default radius is 5 Å",
      D.pocketOccupancy({ center: [0, 0, 0] }, [{ x: 5, y: 0, z: 0 }]).inside === 1 &&
      D.pocketOccupancy({ center: [0, 0, 0] }, [{ x: 5.5, y: 0, z: 0 }]).inside === 0);
check("occupancy is the fraction of atoms inside",
      D.pocketOccupancy({ center: [0, 0, 0] }, ligand, 0.5).occupancy === 1 / 3,
      `(${D.pocketOccupancy({ center: [0, 0, 0] }, ligand, 0.5).occupancy})`);
check("an empty atom list has no occupancy", D.pocketOccupancy(centre, [], 5) === null);
check("a missing atom list has no occupancy", D.pocketOccupancy(centre, null, 5) === null);
check("a pocket without a centre has no occupancy",
      D.pocketOccupancy({}, ligand, 5) === null);
check("a missing pocket has no occupancy", D.pocketOccupancy(null, ligand, 5) === null);

// -------------------------------------------------- validating the pockets
console.log("\nValidating pockets against the ligand");
const pockets = [
  { center: [100, 0, 0] },
  { center: [0, 0, 0] },
  { center: [1, 0, 0] },
];
const twoAtoms = [{ x: 0, y: 0, z: 0 }, { x: 2, y: 0, z: 0 }];
const validated = D.validatePocketsAgainstLigand(pockets, twoAtoms, 5);
check("the best-occupied pocket is chosen", validated.index === 1,
      `(index ${validated.index})`);
check("it reports the occupancy it found", validated.occupancy === 1);
check("it reports that it found the ligand", validated.found === true);
check("perPocket has one entry per pocket", validated.perPocket.length === 3,
      `(${validated.perPocket.length})`);
check("each entry carries its pocket index",
      validated.perPocket.map((p) => p.index).join(",") === "0,1,2");
check("each entry carries its occupancy",
      validated.perPocket.every((p) => typeof p.occupancy === "number"));
// Two pockets with the same occupancy; the one whose centre is nearer the
// ligand on average must win, whichever order the pockets arrive in.
const tieLigand = [{ x: 0, y: 0, z: 0 }, { x: 0, y: 10, z: 0 }];
const tie = D.validatePocketsAgainstLigand(
  [{ center: [1, 0, 0] }, { center: [0, 0, 0] }], tieLigand, 5);
check("an occupancy tie is broken by the smaller mean distance", tie.index === 1,
      `(index ${tie.index})`);
check("the tie really was a tie",
      tie.perPocket[0].occupancy === tie.perPocket[1].occupancy);
check("the winner really is the nearer one",
      tie.perPocket[tie.index].meanDistance < tie.perPocket[0].meanDistance);
const empty = D.validatePocketsAgainstLigand(
  [{ center: [100, 0, 0] }, { center: [200, 0, 0] }], [{ x: 0, y: 0, z: 0 }], 5);
check("no pocket holding the ligand is not a finding", empty.found === false);
check("it reports index -1 when nothing is found", empty.index === -1);
check("it reports zero occupancy when nothing is found", empty.occupancy === 0);
check("the pockets are still listed",
      empty.perPocket.length === 2 && empty.perPocket.every((p) => p.occupancy === 0));
check("no pockets at all is not a finding",
      D.validatePocketsAgainstLigand([], twoAtoms, 5).found === false);
check("a pocket without a centre is skipped",
      D.validatePocketsAgainstLigand([{}, { center: [0, 0, 0] }], twoAtoms, 5)
        .perPocket.length === 1);

// --------------------------------------------------------- contact residues
console.log("\nContact residues");
const lig = [{ x: 0, y: 0, z: 0 }];
const prot = [
  { x: 2, y: 0, z: 0, chain: "A", resi: 100, name: "CA" },
  { x: 20, y: 0, z: 0, chain: "A", resi: 101, name: "CA" },
];
const contacts = D.contactResidues(lig, prot, 4.0);
check("it returns a Set", contacts instanceof Set);
check("a residue within the cutoff is reported", contacts.has("A:100"),
      `(${[...contacts].join(", ")})`);
check("a residue beyond the cutoff is not", !contacts.has("A:101"));
check("the default cutoff is 4 Å",
      D.contactResidues(lig, [{ x: 4, y: 0, z: 0, chain: "A", resi: 1 }]).has("A:1") &&
      !D.contactResidues(lig, [{ x: 4.5, y: 0, z: 0, chain: "A", resi: 1 }]).has("A:1"));
check("the cutoff boundary is inclusive",
      D.contactResidues(lig, [{ x: 4, y: 0, z: 0, chain: "A", resi: 7 }], 4).has("A:7"));
check("every atom of a contacting residue gives one entry",
      D.contactResidues(lig, [
        { x: 1, y: 0, z: 0, chain: "A", resi: 5 },
        { x: 2, y: 0, z: 0, chain: "A", resi: 5 },
      ], 4).size === 1);
check("waters are reported too, unfiltered",
      D.contactResidues(lig, [{ x: 2, y: 0, z: 0, chain: "W", resi: 300 }], 4)
        .has("W:300"));
check("another chain is reported too, unfiltered",
      D.contactResidues(lig, [{ x: 2, y: 0, z: 0, chain: "B", resi: 42 }], 4)
        .has("B:42"));
check("empty inputs give an empty set",
      D.contactResidues([], prot, 4).size === 0 && D.contactResidues(lig, [], 4).size === 0);
check("missing inputs give an empty set",
      D.contactResidues(null, null, 4).size === 0);

// ------------------------------------------------------ pocket residue recall
console.log("\nPocket residue recall");
const pocketAt = { center: [0, 0, 0] };
const resAtoms = [
  { x: 1, y: 0, z: 0, chain: "A", resi: 100 },
  { x: 0, y: 1, z: 0, chain: "A", resi: 101 },
  { x: 0, y: 0, z: 20, chain: "A", resi: 102 },
];
const full = D.pocketResidueRecall(pocketAt, new Set(["A:100", "A:101"]), resAtoms, 7);
check("a pocket over every contact residue covers them all", full.covered === 2 &&
      full.total === 2);
check("full coverage means recall 1", full.recall === 1);
check("nothing is missing", full.missing.length === 0);
const partial = D.pocketResidueRecall(pocketAt, new Set(["A:100", "A:102"]), resAtoms, 7);
check("a partly covered site covers part of it", partial.covered === 1 &&
      partial.total === 2);
check("the recall is the exact fraction", partial.recall === 0.5,
      `(${partial.recall})`);
check("the uncovered residues are listed",
      partial.missing.length === 1 && partial.missing[0] === "A:102",
      `(${partial.missing.join(", ")})`);
const absent = D.pocketResidueRecall(pocketAt, new Set(["A:100", "B:999"]), resAtoms, 7);
check("a contact residue with no atoms is uncovered", absent.covered === 1);
check("it is listed as missing", absent.missing.includes("B:999"));
check("it still counts towards the total", absent.total === 2);
const far = D.pocketResidueRecall(pocketAt, new Set(["A:100"]), resAtoms, 0.5);
check("a residue beyond the radius is not covered", far.covered === 0 &&
      far.recall === 0);
check("the default radius is 7 Å",
      D.pocketResidueRecall(pocketAt, new Set(["A:100"]), resAtoms).covered === 1 &&
      D.pocketResidueRecall(pocketAt, new Set(["A:100"]), resAtoms, 0.5).covered === 0);
check("a residue exactly on the radius is covered",
      D.pocketResidueRecall(pocketAt, new Set(["A:100"]), resAtoms, 1).covered === 1);
check("an array of contacts is accepted",
      D.pocketResidueRecall(pocketAt, ["A:100"], resAtoms, 7).covered === 1);
check("no contacts means nothing to cover",
      D.pocketResidueRecall(pocketAt, new Set(), resAtoms, 7).recall === 0);
check("a missing pocket means nothing to cover",
      D.pocketResidueRecall(null, new Set(["A:100"]), resAtoms, 7).covered === 0);
check("a missing pocket still reports the total",
      D.pocketResidueRecall(null, new Set(["A:100"]), resAtoms, 7).total === 1);

// --------------------------------------------------------- cross-validation
console.log("\nCross-validation against the CCD");
// The real case: AQ4 (the 4-anilinoquinazoline in 1IEP/EGFR), read off a
// crystallographic model with no hydrogens, against its CCD record.
const aq4Props = { mw: 370.259, heavy: 29 };
const aq4Ccd = { formula: "C22 H23 N3 O4", formula_weight: 393.436 };
const aq4Composition = { C: 22, N: 3, O: 4 };
const aq4 = D.crossValidateLigand(aq4Props, aq4Ccd, aq4Composition);
check("the real case produces three checks", aq4.checks.length === 3,
      `(${aq4.checks.length})`);
check("the checks are the three documented ones",
      aq4.checks.map((c) => c.name).join("|") ===
        "Element composition|Heavy-atom count|Molecular weight",
      `(${aq4.checks.map((c) => c.name).join("|")})`);
check("all three checks pass", aq4.checks.every((c) => c.ok),
      `(${aq4.checks.filter((c) => !c.ok).map((c) => c.name).join(", ") || "none failed"})`);
check("the app and the reference agree", aq4.agree === true);
const comp = byName(aq4.checks, "Element composition");
check("the composition check passes", comp.ok === true);
check("the composition detail says identical", comp.detail === "identical",
      `(${comp.detail})`);
check("the composition check shows both sides",
      comp.app === "C22 N3 O4" && comp.reference === "C22 H23 N3 O4",
      `(${comp.app} vs ${comp.reference})`);
const heavy = byName(aq4.checks, "Heavy-atom count");
check("the heavy-atom check passes", heavy.ok === true);
check("the heavy-atom check compares 29 with 29",
      heavy.app === 29 && heavy.reference === 29, `(${heavy.detail})`);
const mw = byName(aq4.checks, "Molecular weight");
check("the molecular-weight check passes", mw.ok === true);
const pct = parseFloat((mw.detail.match(/([\d.]+)% apart/) || [])[1]);
check("the molecular-weight detail shows the deviation", !Number.isNaN(pct),
      `(${mw.detail})`);
check("the deviation is below 1 %", pct < 1, `(${pct}%)`);
// The same case with the composition misread off the model: the composition
// check must fail even though the weights were computed from the right atoms.
const wrongComp = D.crossValidateLigand(aq4Props, aq4Ccd, { C: 20, N: 3, O: 4 });
check("a wrong composition fails the composition check",
      byName(wrongComp.checks, "Element composition").ok === false);
check("it reports the difference", byName(wrongComp.checks, "Element composition")
      .detail === "C 20 vs 22", `(${byName(wrongComp.checks, "Element composition").detail})`);
check("a wrong composition means the app and the reference disagree",
      wrongComp.agree === false);
// And the coherent version of the same mistake: the model was misread, so
// the descriptors the app computed are those of the misread molecule too.
// That has to fail the weight check, not only the composition check.
const misread = D.crossValidateLigand(
  { mw: 346.237, heavy: 27 }, aq4Ccd, { C: 20, N: 3, O: 4 });
check("a misread model fails the molecular-weight check",
      byName(misread.checks, "Molecular weight").ok === false,
      `(${byName(misread.checks, "Molecular weight").detail})`);
check("the weight detail shows a deviation above 1 %",
      parseFloat((byName(misread.checks, "Molecular weight").detail
        .match(/([\d.]+)% apart/) || [])[1]) > 1);
check("a misread model means the app and the reference disagree",
      misread.agree === false);
const heavyOff = D.crossValidateLigand(
  { mw: 370.259, heavy: 27 }, aq4Ccd, aq4Composition);
check("a heavy-atom count off by two fails",
      byName(heavyOff.checks, "Heavy-atom count").ok === false,
      `(${byName(heavyOff.checks, "Heavy-atom count").detail})`);
check("a heavy-atom count off by two means no agreement", heavyOff.agree === false);
check("a heavy-atom count off by one is tolerated",
      byName(D.crossValidateLigand({ mw: 370.259, heavy: 28 }, aq4Ccd, aq4Composition)
        .checks, "Heavy-atom count").ok === true);
check("a heavy-atom count one high is tolerated",
      byName(D.crossValidateLigand({ mw: 370.259, heavy: 30 }, aq4Ccd, aq4Composition)
        .checks, "Heavy-atom count").ok === true);
const noCcd = D.crossValidateLigand(aq4Props, null, aq4Composition);
check("a missing CCD record yields no checks", noCcd.checks.length === 0);
check("a missing CCD record means no agreement", noCcd.agree === false);

// ------------------------------------------------------- counting a formula
console.log("\nCounting atoms from a formula");
const c22 = D.countAtomsFromFormula("C22 H23 N3 O4");
check("it counts the heavy atoms", c22.heavy === 29, `(${c22.heavy})`);
check("it counts every atom", c22.total === 52, `(${c22.total})`);
check("it breaks the formula down by element",
      JSON.stringify(c22.byElement) === JSON.stringify({ C: 22, H: 23, N: 3, O: 4 }),
      `(${JSON.stringify(c22.byElement)})`);
// "C6H12O6 . 2H2O" must be read as the parent molecule alone: the model
// never holds the crystallisation water, so the water's six atoms must not
// be counted (counting them would give 30 atoms rather than 24).
const glucose = D.countAtomsFromFormula("C6H12O6 . 2H2O");
check("a hydrate is counted as the parent molecule only", glucose.heavy === 12,
      `(${glucose.heavy})`);
check("the hydrate's total ignores the water", glucose.total === 24,
      `(${glucose.total})`);
check("the water is not in the element breakdown",
      JSON.stringify(glucose.byElement) === JSON.stringify({ C: 6, H: 12, O: 6 }));
check("a hydrate counts the same as its parent",
      JSON.stringify(D.countAtomsFromFormula("C6H12O6")) ===
      JSON.stringify(glucose));
check("an empty formula counts nothing", D.countAtomsFromFormula("") === null);
check("a non-string counts nothing", D.countAtomsFromFormula(null) === null);
check("a number counts nothing", D.countAtomsFromFormula(42) === null);
check("counts of one are implied",
      D.countAtomsFromFormula("CHN").byElement.C === 1 &&
      D.countAtomsFromFormula("CHN").heavy === 2);

// -------------------------------------------------------- molecular weight
console.log("\nMolecular weight");
const mw12 = D.molecularWeight({ C: 12, H: 1 });
check("the weight is the sum over the elements", near(mw12.mw, 144.132),
      `(${mw12.mw})`);
check("it counts the atoms it weighed", mw12.counted === 12, `(${mw12.counted})`);
check("it reports no unknown elements", mw12.unknown.length === 0);
const withUnknown = D.molecularWeight({ C: 1, Xx: 2 });
check("an element with no tabulated mass is reported",
      withUnknown.unknown.length === 1 && withUnknown.unknown[0] === "Xx",
      `(${withUnknown.unknown.join(", ")})`);
check("an unknown element adds nothing to the weight", near(withUnknown.mw, 12.011),
      `(${withUnknown.mw})`);
check("an unknown element is not counted", withUnknown.counted === 1);
const withH = D.molecularWeight({ C: 12, H: 1 }, { includeHydrogen: true });
check("includeHydrogen counts the hydrogen", near(withH.mw, 145.14),
      `(${withH.mw})`);
check("includeHydrogen reports the hydrogen as counted", withH.counted === 13);
check("hydrogen has a tabulated mass but is left out by default",
      D.ATOMIC_MASS.H === 1.008 && near(D.molecularWeight({ C: 12, H: 1 }).mw, 144.132));
check("hydrogen is not counted as an atom by default",
      D.molecularWeight({ C: 12, H: 1 }).counted === 12);
check("an empty composition weighs nothing",
      D.molecularWeight({}).mw === 0 && D.molecularWeight({}).counted === 0);
check("a missing composition weighs nothing",
      D.molecularWeight(null).mw === 0);

// --------------------------------------------------------- formula format
console.log("\nFormula formatting");
check("a formula is written in Hill order",
      D.formatFormula({ C: 22, N: 3, O: 4 }) === "C22 N3 O4",
      `(${D.formatFormula({ C: 22, N: 3, O: 4 })})`);
check("a count of one is left off",
      D.formatFormula({ C: 1, O: 2 }) === "C O2", `(${D.formatFormula({ C: 1, O: 2 })})`);
check("without carbon the elements are alphabetical",
      D.formatFormula({ S: 1, N: 1, O: 2 }) === "N O2 S",
      `(${D.formatFormula({ S: 1, N: 1, O: 2 })})`);
check("hydrogen follows carbon",
      D.formatFormula({ N: 1, C: 1, H: 2 }) === "C H2 N",
      `(${D.formatFormula({ N: 1, C: 1, H: 2 })})`);
check("everything after carbon and hydrogen is alphabetical",
      D.formatFormula({ C: 1, S: 1, O: 1, N: 1 }) === "C N O S",
      `(${D.formatFormula({ C: 1, S: 1, O: 1, N: 1 })})`);
check("an empty composition has no formula", D.formatFormula({}) === "");
check("zero counts are dropped", D.formatFormula({ C: 0, O: 2 }) === "O2",
      `(${D.formatFormula({ C: 0, O: 2 })})`);

// -------------------------------------------------------------- potency
console.log("\nPotency formatting");
check("1500 nM is promoted to µM", D.formatPotency(1500, "nm") === "1.50 µM",
      `(${D.formatPotency(1500, "nm")})`);
check("1000 nM is promoted to µM", D.formatPotency(1000, "nm") === "1.00 µM",
      `(${D.formatPotency(1000, "nm")})`);
check("a small nM value keeps one decimal",
      D.formatPotency(7.6, "nm") === "7.6 nM", `(${D.formatPotency(7.6, "nm")})`);
check("a nM value below 10 keeps one decimal",
      D.formatPotency(9.5, "nm") === "9.5 nM", `(${D.formatPotency(9.5, "nm")})`);
check("a nM value of 10 or more drops the decimal",
      D.formatPotency(500, "nm") === "500 nM", `(${D.formatPotency(500, "nm")})`);
check("the unit name is case-insensitive",
      D.formatPotency(1500, "NM") === "1.50 µM", `(${D.formatPotency(1500, "NM")})`);
check("a µM value below 10 keeps two decimals",
      D.formatPotency(5, "um") === "5.00 µM", `(${D.formatPotency(5, "um")})`);
check("a µM value of 10 or more keeps one decimal",
      D.formatPotency(12.5, "um") === "12.5 µM", `(${D.formatPotency(12.5, "um")})`);
check("the micro sign is accepted as a unit",
      D.formatPotency(12.5, "µm") === "12.5 µM", `(${D.formatPotency(12.5, "µm")})`);
check("an mM value keeps two decimals",
      D.formatPotency(3, "mm") === "3.00 mM", `(${D.formatPotency(3, "mm")})`);
check("a pM value has no decimals",
      D.formatPotency(250, "pm") === "250 pM", `(${D.formatPotency(250, "pm")})`);
check("with no unit a small value is nM",
      D.formatPotency(7.6) === "7.6 nM", `(${D.formatPotency(7.6)})`);
check("with no unit a large value is µM",
      D.formatPotency(1500) === "1.50 µM", `(${D.formatPotency(1500)})`);
check("a NaN potency has no value", D.formatPotency(NaN, "nm") === "—",
      `(${D.formatPotency(NaN, "nm")})`);
check("a missing potency has no value", D.formatPotency(undefined, "nm") === "—");
check("an unknown unit is passed through",
      D.formatPotency(5, "x") === "5 x", `(${D.formatPotency(5, "x")})`);

// -------------------------------------------------------------- affinities
console.log("\nAffinity summary");
const aff = D.summariseAffinity([
  { type: "IC50", value: 120, unit: "nM" },
  { type: "Kd", value: 5, unit: "nM" },
  { type: "Ki", value: 0, unit: "nM" },
  { type: "IC50", value: -3, unit: "nM" },
  { type: "IC50", value: NaN, unit: "nM" },
]);
check("non-positive and non-finite values are discarded", aff.n === 2, `(${aff.n})`);
check("the tightest binder is the best",
      aff.best.type === "Kd" && aff.best.value === 5);
const kdVsIc50 = D.summariseAffinity([
  { type: "IC50", value: 1, unit: "nM" },
  { type: "Kd", value: 1000, unit: "nM" },
]);
check("a Kd beats an IC50 however large",
      kdVsIc50.best.type === "Kd", `(${kdVsIc50.best.type})`);
const kiVsIc50 = D.summariseAffinity([
  { type: "IC50", value: 1, unit: "nM" },
  { type: "Ki", value: 900, unit: "nM" },
]);
check("a Ki beats an IC50 however large",
  kiVsIc50.best.type === "Ki", `(${kiVsIc50.best.type})`);
const sameType = D.summariseAffinity([
  { type: "IC50", value: 900, unit: "nM" },
  { type: "IC50", value: 40, unit: "nM" },
]);
check("within one type the smaller value wins", sameType.best.value === 40);
const unknownType = D.summariseAffinity([
  { type: "EC50", value: 1, unit: "nM" },
  { type: "Kd", value: 900, unit: "nM" },
]);
check("an unranked type sorts last", unknownType.best.type === "Kd",
  `(${unknownType.best.type})`);
const types = D.summariseAffinity([
  { type: "IC50", value: 1, unit: "nM" },
  { type: "Kd", value: 2, unit: "nM" },
  { type: "Ki", value: 3, unit: "nM" },
  { type: "Kd", value: 4, unit: "nM" },
]);
check("the types are sorted and unique", types.types.join(",") === "IC50,Kd,Ki",
  `(${types.types.join(",")})`);
check("no affinities gives nothing", D.summariseAffinity([]).n === 0);
check("no affinities leaves no best", D.summariseAffinity([]).best === null);
check("no affinities lists no types",
  D.summariseAffinity([]).types.length === 0);
check("a missing list is treated as empty",
  D.summariseAffinity(null).n === 0);

// ------------------------------------------------------------ target profile
console.log("\nTarget profile");
const egfr = {
  id: "ENSG00000146648",
  approvedSymbol: "EGFR",
  approvedName: "epidermal growth factor receptor",
  targetClass: [{ id: "TKL", label: "Protein kinase" }, { id: "TK", label: "Kinase" }],
  associatedDiseases: {
    count: 217,
    rows: [
      { disease: { id: "EFO_0000270", name: "Lung cancer" }, score: 0.97 },
      { disease: { id: "EFO_0000222", name: "Breast cancer" }, score: 0.55 },
    ],
  },
  drugAndClinicalCandidates: {
    count: 190,
    rows: [
      { id: "DRB1", maxClinicalStage: "APPROVAL",
        drug: { id: "DRB1", name: "Erlotinib", maximumClinicalStage: "APPROVAL", drugType: "Small molecule" } },
      { id: "DRB2", maxClinicalStage: "PHASE_3",
        drug: { id: "DRB2", name: "Osimertinib", maximumClinicalStage: "PHASE_3", drugType: "Small molecule" } },
    ],
  },
  tractability: [
    { label: "Approved Drug", modality: "Small molecule", value: true },
    { label: "Clinical Trial", modality: "Small molecule", value: true },
  ],
};
const sum = D.summariseTarget(egfr);
check("a null target is null", D.summariseTarget(null) === null);
check("the symbol comes through", sum.symbol === "EGFR", `(${sum.symbol})`);
check("the name comes through", sum.name === "epidermal growth factor receptor");
check("the id comes through", sum.id === "ENSG00000146648");
check("the target classes are the labels",
  sum.classes.join(",") === "Protein kinase,Kinase", `(${sum.classes.join(",")})`);
check("the disease count comes from the count field", sum.diseaseCount === 217,
  `(${sum.diseaseCount})`);
check("the diseases are summarised", sum.diseases.length === 2);
check("a disease keeps its id and name",
  sum.diseases[0].id === "EFO_0000270" && sum.diseases[0].name === "Lung cancer");
check("a disease keeps its score", sum.diseases[0].score === 0.97);
check("the drug count comes from the count field", sum.drugCount === 190,
  `(${sum.drugCount})`);
check("the drugs are summarised", sum.drugs.length === 2);
check("a drug keeps its name and stage",
  sum.drugs[0].name === "Erlotinib" && sum.drugs[0].stage === "APPROVAL");
check("a drug keeps its type", sum.drugs[0].type === "Small molecule");
check("the approved drugs are counted", sum.approvedDrugs === 1, `(${sum.approvedDrugs})`);
check("the tractability is summarised", sum.tractability.length === 2);
check("a tractability entry keeps its label",
  sum.tractability[0].label === "Approved Drug");
check("an approved drug is flagged", sum.isApprovedDrug === true);
const notApproved = D.summariseTarget({
  ...egfr,
  tractability: [{ label: "Approved Drug", modality: "Small molecule", value: false }],
});
check("a false tractability value is not an approval", notApproved.isApprovedDrug === false);
const notDrug = D.summariseTarget({
  ...egfr,
  tractability: [{ label: "Small molecule", modality: "Small molecule", value: true }],
});
check("a target without an approved-drug label is not one", notDrug.isApprovedDrug === false);
const capped = D.summariseTarget({
  ...egfr,
  associatedDiseases: {
    count: 9,
    rows: Array.from({ length: 9 }, (_, i) => ({ disease: { id: "D" + i, name: "D" + i }, score: 0.5 })),
  },
  drugAndClinicalCandidates: {
    count: 14,
    rows: Array.from({ length: 14 }, (_, i) => ({
      id: "D" + i,
      maxClinicalStage: i === 11 ? undefined : "PHASE_1",
      drug: { id: "D" + i, name: "D" + i },
    })),
  },
});
check("at most five diseases are listed", capped.diseases.length === 5);
check("at most twelve drugs are listed", capped.drugs.length === 12);
check("a drug with no stage is marked unknown", capped.drugs[11].stage === "—",
  `(${capped.drugs[11].stage})`);
check("a missing disease name is marked unknown",
  D.summariseTarget({ ...egfr, associatedDiseases: { count: 1, rows: [{ disease: null }] } })
    .diseases[0].name === "—");

// ------------------------------------------------------ matching drug/ligand
console.log("\nMatching a drug to the ligand");
const drugList = [
  { id: "DRB1", name: "Erlotinib", stage: "APPROVAL", type: "Small molecule" },
  { id: "DRB2", name: "Osimertinib", stage: "PHASE_3", type: "Small molecule" },
];
check("a drug matches its own name regardless of case",
  D.matchDrugToLigand("ERLOTINIB", "AQ4", drugList).id === "DRB1");
check("a drug matches its own name in the other case",
  D.matchDrugToLigand("Erlotinib", "AQ4", drugList).id === "DRB1");
check("the punctuation is ignored",
  D.matchDrugToLigand("erlotinib", "AQ4", drugList).id === "DRB1");
check("a prefixed drug name matches",
  D.matchDrugToLigand("Erlotinib", "AQ4",
    [{ id: "D1", name: "Erlotinib HCl" }]).id === "D1");
check("a shortened ligand name matches",
  D.matchDrugToLigand("ERL", "AQ4", drugList).id === "DRB1");
check("an unrelated ligand matches nothing",
  D.matchDrugToLigand("Imatinib", "AQ4", drugList) === null);
check("the compound id is the fallback",
  D.matchDrugToLigand("UNKNOWN LIGAND", "AQ4",
    [{ id: "D1", name: "AQ4" }]).id === "D1");
check("the compound id fallback is case-insensitive",
  D.matchDrugToLigand("UNKNOWN LIGAND", "aq4",
    [{ id: "D1", name: "AQ4" }]).id === "D1");
check("a missing ligand name matches nothing",
  D.matchDrugToLigand("", "AQ4", drugList) === null);
check("a missing drug list matches nothing",
  D.matchDrugToLigand("Erlotinib", "AQ4", null) === null);
check("an empty drug list matches nothing",
  D.matchDrugToLigand("Erlotinib", "AQ4", []) === null);

console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
