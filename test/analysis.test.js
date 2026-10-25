/**
 * Unit tests for the extracted pure analysis modules.
 *
 *   node test/analysis.test.js
 *
 * These modules (src/seq.js, src/geometry.js, src/pockets.js) have no DOM,
 * network or viewer dependency, so they load with a bare require() -- no
 * stubbed document, no IIFE probe injection. That property is the whole
 * point of the extraction, and it is asserted explicitly below.
 */
const path = require("path");

// Load the modules the same way a bundler-free browser would: they publish
// themselves on the global PFE namespace.
require(path.join(__dirname, "..", "src", "seq.js"));
require(path.join(__dirname, "..", "src", "geometry.js"));
require(path.join(__dirname, "..", "src", "pockets.js"));

const { seq, geometry, pockets } = globalThis.PFE;

let pass = 0, fail = 0;
const near = (a, b, tol) => Math.abs(a - b) <= tol;
function check(name, cond, detail = "") {
  if (cond) { pass++; console.log(`  PASS  ${name}${detail ? "  " + detail : ""}`); }
  else { fail++; console.log(`  FAIL  ${name}${detail ? "  " + detail : ""}`); }
}

console.log("\nModule isolation");
check("no DOM required", typeof document === "undefined");
check("no window required", typeof window === "undefined");
check("PFE namespace present", !!(seq && geometry && pockets));
check("namespaces are frozen", Object.isFrozen(seq) && Object.isFrozen(geometry) && Object.isFrozen(pockets));

// ---------------------------------------------------------------- sequences
console.log("\nResidue code mapping");
check("ALA -> A",  seq.AA3TO1.ALA === "A");
check("MSE -> M",  seq.AA3TO1.MSE === "M");
check("SEC -> U",  seq.AA3TO1.SEC === "U");
check("TYR -> Y",  seq.AA3TO1.TYR === "Y");

// A two-chain PDB with one helix and one sheet, plus an altloc duplicate.
// Fields are placed at exact PDB column offsets so the parser's fixed-width
// substring() reads land where they are meant to.
function place(fields) {
  const out = [];
  let pos = 0;
  for (const [start, text] of fields) {
    while (pos < start) { out.push(" "); pos++; }
    out.push(text); pos += text.length;
  }
  return out.join("");
}
// ATOM record: name 13-16, altLoc 17, resName 18-20, chain 22, resSeq 23-26,
// iCode 27, x 31-38 (all 1-indexed).
function pdbLine(chain, resi, resn, atom, x, y, z, altloc) {
  return place([
    [0, "ATOM  "], [6, String(100 + resi)], [11, " "], [12, atom],
    [16, altloc], [17, resn], [20, " "], [21, chain], [22, String(resi)],
    [26, "    "], [30, String(x)], [38, String(y)], [46, String(z)],
    [54, "  1.00  0.00"], [76, atom.charAt(0)],
  ]);
}
// HELIX: initChainID col 20, initSeqNum 22-25, endSeqNum 34-37.
const helixLine = place([
  [0, "HELIX "], [7, "  1"], [12, " 1"], [15, "THR"], [19, "A"],
  [21, "   2"], [27, "GLY"], [31, "A"], [33, "   3"], [38, "  1"],
]);
// SHEET: initChainID col 22, initSeqNum 23-26, endSeqNum 34-37.
const sheetLine = place([
  [0, "SHEET "], [7, "  1"], [11, "A 1"], [17, "PHE"], [21, "B"],
  [22, "   1"], [27, "LEU"], [31, "B"], [33, "   2"], [38, "  1"],
]);
const miniPdb = [
  helixLine,
  sheetLine,
  pdbLine("A", 1, "THR", "CA", 0, 0, 0, " "),
  pdbLine("A", 2, "THR", "CA", 3, 0, 0, " "),
  pdbLine("A", 3, "GLY", "CA", 6, 0, 0, " "),
  pdbLine("A", 3, "GLY", "CA", 6, 0, 0, "B"),   // altloc duplicate of residue 3
  pdbLine("B", 1, "PHE", "CA", 0, 10, 0, " "),
  pdbLine("B", 2, "LEU", "CA", 3, 10, 0, " "),
].join("\n");

const parsed = seq.parseSequenceFromPDB(miniPdb, "pdb");
console.log("\nPDB sequence parsing");
check("two chains found", parsed.chains.length === 2, `(${parsed.chains.length})`);
check("chain ids A and B", parsed.chains[0].id === "A" && parsed.chains[1].id === "B");
check("chain A has 3 residues", parsed.chains[0].residues.length === 3,
      `(${parsed.chains[0].residues.length})`);
check("altloc duplicate dropped", parsed.chains[0].residues.length === 3);
check("one-letter codes", parsed.chains[0].residues[0].one === "T" &&
      parsed.chains[0].residues[2].one === "G");
check("resSeq parsed", parsed.chains[0].residues[2].resi === 3);
check("HELIX marks helix", parsed.chains[0].residues[1].ss === "H",
      `("${parsed.chains[0].residues[1].ss}")`);
check("SHEET marks strand", parsed.chains[1].residues[0].ss === "E",
      `("${parsed.chains[1].residues[0].ss}")`);
check("residues outside assignment have no ss", parsed.chains[0].residues[0].ss === "");
check("unknown residue maps to X",
      seq.parseSequenceFromPDB(pdbLine("A", 1, "UNK", "CA", 0, 0, 0, " "), "pdb")
        .chains[0].residues[0].one === "X");

console.log("\nFASTA parsing");
const fasta = seq.parseFastaSequences(">sp1 desc\nMKWV\nTFIS\n>sp2\nMKWV\n>sp3\nQQQQ\n");
check("three records", fasta.length === 3, `(${fasta.length})`);
check("multi-line record joined", fasta[0] === "MKWVTFIS", `("${fasta[0]}")`);
check("second record", fasta[1] === "MKWV");
check("third record", fasta[2] === "QQQQ");

// ---------------------------------------------------------------- alignment
console.log("\nNeedleman-Wunsch alignment");
const same = seq.nwAlign("ACDEFGHIK", "ACDEFGHIK");
check("identical sequences stay ungapped", !same.a.includes("-") && !same.b.includes("-"));
check("identical sequences preserved", same.a === "ACDEFGHIK" && same.b === "ACDEFGHIK");

const shifted = seq.nwAlign("ACDEFGHIK", "ACDEFXGHIK");
check("aligned outputs have equal length", shifted.a.length === shifted.b.length,
      `(${shifted.a.length} vs ${shifted.b.length})`);
check("insertion is gapped", shifted.a.includes("-") || shifted.b.includes("-"));

const short = seq.nwAlign("ACDE", "ACDEFGHI");
check("gap introduced for the tail", short.a.includes("-") && !short.b.includes("-"),
      `("${short.a}")`);

// Regression: the traceback used to run off the end of the row whenever the
// two sequences differed enough in length. At j=0 with i>0 it read an
// unwritten traceback cell, took the "left" branch, and indexed b[-1] --
// undefined -- so j counted down forever and the alignment string grew until
// V8 threw "Invalid string length". It affected any alignment of sequences
// with very different lengths, e.g. a 1210-residue UniProt sequence against
// a 312-residue crystallographic construct.
console.log("\nAlignment termination (regression)");
{
  const long_ = seq.nwAlign("M".repeat(1200), "M".repeat(312));
  check("very different lengths terminate", long_.a.length === long_.b.length,
        `(len ${long_.a.length}, ${long_.b.length})`);
  check("the short sequence survives intact",
        long_.b.replace(/-/g, "").length === 312);
  const empty = seq.nwAlign("ACDE", "");
  check("an empty second sequence terminates", empty.a.length === empty.b.length);
  const emptyFirst = seq.nwAlign("", "ACDE");
  check("an empty first sequence terminates", emptyFirst.a.length === emptyFirst.b.length);
  // A gap must never be aligned to a gap.
  let bothGapped = false;
  for (let k = 0; k < short.a.length; k++) {
    if (short.a[k] === "-" && short.b[k] === "-") bothGapped = true;
  }
  check("a gap is never aligned to a gap", !bothGapped);
  // Many random pairs of different lengths, to shake out any remaining
  // non-termination.
  let terminated = 0;
  const AAs = "ACDEFGHIKLMNPQRSTVWY";
  let seed = 24680;
  const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff; };
  for (let t = 0; t < 200; t++) {
    const la = 1 + Math.floor(rnd() * 30), lb = 1 + Math.floor(rnd() * 30);
    let a = "", b = "";
    for (let k = 0; k < la; k++) a += AAs[Math.floor(rnd() * 20)];
    for (let k = 0; k < lb; k++) b += AAs[Math.floor(rnd() * 20)];
    const r = seq.nwAlign(a, b);
    if (r.a.length === r.b.length) terminated++;
  }
  check("200 random pairs of differing length all terminate", terminated === 200,
        `(${terminated}/200)`);
}

console.log("\nProgressive MSA");
const fam = [
  "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGE",
  "MKWVTFISLLLLFSSAYSRGVFRRDAHKSEIAHRFKDLGE",
  "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE",
  "MQWVTFLLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGE",
];
const msa = seq.progressiveMsa(fam);
check("all rows same length as consensus",
      msa.aligned.every((r) => r.length === msa.consensus.length));
check("row count preserved", msa.aligned.length === fam.length, `(${msa.aligned.length})`);
check("consensus non-empty", msa.consensus.length > 0);
check("single sequence is a valid MSA", seq.progressiveMsa(["ACDE"]).aligned.length === 1);

console.log("\nConservation scoring");
// progressiveMsa emits rows as arrays of single characters, so the fixture
// must use that shape too: ["A","C","D"], not ["ACD"].
const sc = seq.conservationScores(msa);
check("scores in [0,1]", sc.every((v) => v >= 0 && v <= 1));
check("no NaN", sc.every((v) => !Number.isNaN(v)));
const identicalScores = seq.conservationScores(
  { consensus: "ACD", aligned: [["A", "C", "D"], ["A", "C", "D"], ["A", "C", "D"]] });
check("fully conserved column scores 1", identicalScores.every((v) => v === 1),
      `(${JSON.stringify(identicalScores)})`);
const gapped = seq.conservationScores(
  { consensus: "AC-D", aligned: [["A", "C", "-", "D"], ["A", "C", "-", "D"]] });
check("gap column scores 0", gapped[2] === 0, `(${gapped[2]})`);
const mixed = seq.conservationScores(
  { consensus: "ACD", aligned: [["A", "C", "D"], ["A", "C", "D"], ["A", "C", "E"]] });
check("variable column scores < 1", mixed[2] < 1, `(${mixed[2]})`);
check("conserved columns still score 1", mixed[0] === 1 && mixed[1] === 1);
check("column score is the match fraction",
      near(mixed[2], 2 / 3, 1e-9), `(${mixed[2]})`);

console.log("\nResidue indexing");
const idx = seq.buildResidueIndex(miniPdb);
check("index has 5 residues", idx.map.size === 5, `(${idx.map.size})`);
check("key format chain:resi", idx.map.has("A:1") && idx.map.has("B:2"));
check("seqIdx is monotonic", idx.map.get("A:1").seqIdx === 0 && idx.map.get("B:2").seqIdx === 4);
check("sequence reconstructed", idx.sequence === "TTGFL", `("${idx.sequence}")`);
check("resn preserved", idx.map.get("A:1").resn === "THR");

// ---------------------------------------------------------------- geometry
console.log("\nVector helpers");
const a = { x: 1, y: 0, z: 0 }, b = { x: 0, y: 1, z: 0 };
check("midpoint averages", geometry.midPoint({ x: 0, y: 0, z: 0 }, { x: 4, y: 2, z: 6 }).x === 2);
check("cross of x and y is z", geometry.cross(a, b).z === 1);
check("cross is anti-commutative", geometry.cross(a, b).z === -geometry.cross(b, a).z);
check("dot of x and y is 0", geometry.dot(a, b) === 0);
check("dot of x with itself is 1", geometry.dot(a, a) === 1);

console.log("\nDistances, angles, dihedrals");
check("distance 3-4-0 = 5", near(geometry.computeDistance({x:0,y:0,z:0},{x:3,y:4,z:0}), 5, 1e-9));
check("distance is symmetric",
      near(geometry.computeDistance({x:1,y:2,z:3},{x:4,y:6,z:3}),
           geometry.computeDistance({x:4,y:6,z:3},{x:1,y:2,z:3}), 1e-9));
check("right angle = 90", near(geometry.computeAngle({x:1,y:0,z:0},{x:0,y:0,z:0},{x:0,y:1,z:0}), 90, 1e-6));
check("straight angle = 180", near(geometry.computeAngle({x:1,y:0,z:0},{x:0,y:0,z:0},{x:-1,y:0,z:0}), 180, 1e-6));
check("equilateral angle = 60",
      near(geometry.computeAngle({x:1,y:0,z:0},{x:0,y:0,z:0},{x:0.5,y:Math.sqrt(3)/2,z:0}), 60, 1e-6));
check("trans dihedral = 180",
      near(Math.abs(geometry.computeDihedral({x:1,y:0,z:0},{x:0,y:0,z:0},{x:0,y:1,z:0},{x:-1,y:1,z:0})), 180, 1));
check("cis dihedral = 0",
      near(geometry.computeDihedral({x:1,y:0,z:0},{x:0,y:0,z:0},{x:0,y:1,z:0},{x:1,y:1,z:0}), 0, 1));

console.log("\nPer-residue RMSD");
const fakeModel = (atoms) => ({ selectedAtoms: () => atoms });
const caA = [
  { chain: "A", resi: 1, x: 0, y: 0, z: 0 },
  { chain: "A", resi: 2, x: 1, y: 0, z: 0 },
];
const caB = [
  { chain: "A", resi: 1, x: 0, y: 0, z: 0 },
  { chain: "A", resi: 2, x: 1, y: 0, z: 0 },
];
const caBshift = [
  { chain: "A", resi: 1, x: 3, y: 0, z: 0 },
  { chain: "A", resi: 2, x: 4, y: 0, z: 0 },
];
check("identical models give zero RMSD",
      geometry.computePerResidueRmsd(fakeModel(caA), fakeModel(caB)).every((r) => r.value === 0));
const shiftedRmsd = geometry.computePerResidueRmsd(fakeModel(caA), fakeModel(caBshift));
check("shifted model gives per-residue offset", shiftedRmsd.every((r) => near(r.value, 3, 1e-9)),
      `(${shiftedRmsd.map((r) => r.value).join(", ")})`);
check("unmatched residues are skipped",
      geometry.computePerResidueRmsd(fakeModel(caA), fakeModel([caB[0]])).length === 1);

// ---------------------------------------------------------------- pockets
console.log("\nSpatial hash");
const atoms = [
  { x: 0, y: 0, z: 0 },
  { x: 12, y: 0, z: 0 },
  { x: 25, y: 0, z: 0 },
];
const hash = pockets.buildAtomHash(atoms);
check("three distinct cells", hash.size === 3, `(${hash.size})`);
check("cell keying uses 10 A cells",
      hash.get("0|0|0") && hash.get("1|0|0") && hash.get("2|0|0"));
check("negative coordinates floor correctly",
      pockets.buildAtomHash([{ x: -1, y: 0, z: 0 }]).has("-1|0|0"));

console.log("\nDensity gate counting");
check("counts only atoms inside the radius", pockets.countAtomsNear(hash, 0, 0, 0, 3) === 1);
check("wider radius reaches the next atom", pockets.countAtomsNear(hash, 0, 0, 0, 13) === 2);
check("empty region counts zero", pockets.countAtomsNear(hash, 500, 500, 500, 6) === 0);

console.log("\nRay casting and enclosure");
// The single atom sits far from the query point so that no ray samples
// within the 2.6 A hit tolerance of it.
const lone = pockets.buildAtomHash([{ x: 0, y: 0, z: 0 }]);
check("ray toward an atom hits", pockets.rayHits(lone, 0, 0, 0, 1, 0, 0, 9) === true);
check("ray away from every atom misses", pockets.rayHits(lone, 100, 0, 0, -1, 0, 0, 9) === false);
check("isolated point is not enclosed", near(pockets.enclosureFraction(lone, 100, 0, 0), 0, 1e-9));

// A dense spherical shell: 32 rays must all find material.
const shell = [];
const N = 400;
for (let i = 0; i < N; i++) {
  const y = 1 - (i / (N - 1)) * 2;
  const r = Math.sqrt(Math.max(0, 1 - y * y));
  const th = Math.PI * (3 - Math.sqrt(5)) * i;
  shell.push({ x: Math.cos(th) * r * 6, y: y * 6, z: Math.sin(th) * r * 6 });
}
const shellHash = pockets.buildAtomHash(shell);
check("shell fully encloses the centre", pockets.enclosureFraction(shellHash, 0, 0, 0) === 1,
      `(${pockets.enclosureFraction(shellHash, 0, 0, 0)})`);
check("point outside the shell is not enclosed",
      pockets.enclosureFraction(shellHash, 40, 40, 40) < 0.25);

// ------------------------------------------------------- conservation profile
// conservationProfile() reports three numbers per column:
//   identity - of the sequences that HAVE the column, the fraction that match
//              the consensus (the number the colour map shades by)
//   coverage - the fraction of ALL sequences that have the column at all
//   strength - the fraction of ALL sequences that match the consensus
//              (== identity * coverage)
// The split is the whole point: a column present in a single sequence scores
// identity 1 -- nothing disagrees with it -- but strength 1/n.
console.log("\nConservation profile");
const prof = (consensus, aligned) => seq.conservationProfile({ consensus, aligned });

const profSingle = prof("ACD", ["ACD"]);
check("single row: columns counted", profSingle.columns === 3, `(${profSingle.columns})`);
check("single row: sequences counted", profSingle.sequences === 1, `(${profSingle.sequences})`);
check("single row: identity all 1", profSingle.identity.every((v) => v === 1),
      `(${JSON.stringify(profSingle.identity)})`);
check("single row: coverage all 1", profSingle.coverage.every((v) => v === 1),
      `(${JSON.stringify(profSingle.coverage)})`);
check("single row: strength all 1", profSingle.strength.every((v) => v === 1),
      `(${JSON.stringify(profSingle.strength)})`);

const profTwins = prof("ACD", ["ACD", "ACD"]);
check("identical rows: identity all 1", profTwins.identity.every((v) => v === 1),
      `(${JSON.stringify(profTwins.identity)})`);
check("identical rows: coverage all 1", profTwins.coverage.every((v) => v === 1));
check("identical rows: strength all 1", profTwins.strength.every((v) => v === 1));
check("identical rows: sequences counted", profTwins.sequences === 2, `(${profTwins.sequences})`);

// One column differs. Both rows are the consensus length, so the column is
// present in every sequence and only the disagreement changes the numbers.
const profHalf = prof("ACDT", ["ACDT", "ACDA"]);
check("one of two rows differs: identity is 1/2", near(profHalf.identity[3], 0.5, 1e-9),
      `(${profHalf.identity[3]})`);
check("one of two rows differs: coverage is 1", near(profHalf.coverage[3], 1, 1e-9),
      `(${profHalf.coverage[3]})`);
check("one of two rows differs: strength is 1/2", near(profHalf.strength[3], 0.5, 1e-9),
      `(${profHalf.strength[3]})`);
check("one of two rows differs: other columns stay 1",
      profHalf.identity.slice(0, 3).every((v) => v === 1) &&
      profHalf.coverage.slice(0, 3).every((v) => v === 1) &&
      profHalf.strength.slice(0, 3).every((v) => v === 1));

// The same idea at 3-of-4 rows: identity and strength are 3/4, coverage is 1.
const profThreeOfFour = prof("AAAT", ["AAAT", "AAAT", "AAAT", "AAAA"]);
check("three of four rows match: identity is 3/4", near(profThreeOfFour.identity[3], 0.75, 1e-9),
      `(${profThreeOfFour.identity[3]})`);
check("three of four rows match: coverage is 1", near(profThreeOfFour.coverage[3], 1, 1e-9));
check("three of four rows match: strength is 3/4", near(profThreeOfFour.strength[3], 0.75, 1e-9));

// The key case: the third column exists in only one of the two sequences.
// Nothing disagrees with it, so identity is 1 -- but it is not conserved by
// the alignment, so coverage and strength are 1/2.
const profKey = prof("AAA", ["AAA", "AA-"]);
check("sparse column: identity is 1 (nothing disagrees)", profKey.identity[2] === 1,
      `(${profKey.identity[2]})`);
check("sparse column: coverage is 1/2", near(profKey.coverage[2], 0.5, 1e-9),
      `(${profKey.coverage[2]})`);
check("sparse column: strength is 1/2", near(profKey.strength[2], 0.5, 1e-9),
      `(${profKey.strength[2]})`);
check("sparse column: identity and strength are different numbers",
      profKey.identity[2] !== profKey.strength[2]);
check("sparse column: the other columns are fully conserved",
      profKey.identity[0] === 1 && profKey.identity[1] === 1 &&
      profKey.coverage[0] === 1 && profKey.coverage[1] === 1 &&
      profKey.strength[0] === 1 && profKey.strength[1] === 1);

// A gap in the consensus zeroes the column, even though the rows carry a
// character there.
const profGapCons = prof("AC-D", ["ACXD", "ACXD"]);
check("gap in the consensus: identity is 0", profGapCons.identity[2] === 0,
      `(${profGapCons.identity[2]})`);
check("gap in the consensus: coverage is 0", profGapCons.coverage[2] === 0,
      `(${profGapCons.coverage[2]})`);
check("gap in the consensus: strength is 0", profGapCons.strength[2] === 0,
      `(${profGapCons.strength[2]})`);
check("gap in the consensus: the other columns stay 1",
      profGapCons.identity[0] === 1 && profGapCons.identity[1] === 1 &&
      profGapCons.identity[3] === 1 && profGapCons.coverage[3] === 1 &&
      profGapCons.strength[3] === 1);

// The consensus names a residue that no sequence has.
const profAbsent = prof("ACD", ["A-D", "A-D"]);
check("column absent everywhere: identity is 0", profAbsent.identity[1] === 0,
      `(${profAbsent.identity[1]})`);
check("column absent everywhere: coverage is 0", profAbsent.coverage[1] === 0,
      `(${profAbsent.coverage[1]})`);
check("column absent everywhere: strength is 0", profAbsent.strength[1] === 0,
      `(${profAbsent.strength[1]})`);
check("column absent everywhere: the other columns stay 1",
      profAbsent.identity[0] === 1 && profAbsent.identity[2] === 1 &&
      profAbsent.coverage[2] === 1 && profAbsent.strength[2] === 1);

// One sequence in ten carries the residue: identity 1, strength 0.1.
const profLone = prof("A", ["A", "-", "-", "-", "-", "-", "-", "-", "-", "-"]);
check("one sequence of ten: identity is 1", profLone.identity[0] === 1,
      `(${profLone.identity[0]})`);
check("one sequence of ten: coverage is 0.1", near(profLone.coverage[0], 0.1, 1e-9),
      `(${profLone.coverage[0]})`);
check("one sequence of ten: strength is 0.1", near(profLone.strength[0], 0.1, 1e-9),
      `(${profLone.strength[0]})`);
check("one sequence of ten: sequences counted", profLone.sequences === 10,
      `(${profLone.sequences})`);

console.log("\nConservation profile (mixed alignments)");
// Six rows with four different gap patterns plus one substitution. Column 2
// is present in four rows, three of which match the consensus.
const profMixed = prof("ACDEFGH",
  ["ACDEFGH", "AC-EFGH", "ACD-FGH", "A--EFGH", "ACDEFG-", "ACXEFGH"]);
check("mixed: columns counted", profMixed.columns === 7, `(${profMixed.columns})`);
check("mixed: sequences counted", profMixed.sequences === 6, `(${profMixed.sequences})`);
check("mixed: fully present column is 1/1/1",
      profMixed.identity[0] === 1 && profMixed.coverage[0] === 1 && profMixed.strength[0] === 1);
check("mixed: coverage counts the gaps", near(profMixed.coverage[1], 5 / 6, 1e-9),
      `(${profMixed.coverage[1]})`);
check("mixed: gap-and-mismatch column has identity 3/4", near(profMixed.identity[2], 0.75, 1e-9),
      `(${profMixed.identity[2]})`);
check("mixed: gap-and-mismatch column has coverage 4/6", near(profMixed.coverage[2], 4 / 6, 1e-9),
      `(${profMixed.coverage[2]})`);
check("mixed: gap-and-mismatch column has strength 1/2", near(profMixed.strength[2], 0.5, 1e-9),
      `(${profMixed.strength[2]})`);
check("mixed: strength equals identity times coverage in every column",
      profMixed.strength.every((v, i) => near(v, profMixed.identity[i] * profMixed.coverage[i], 1e-9)));
check("mixed: every identity is within 0..1", profMixed.identity.every((v) => v >= 0 && v <= 1));
check("mixed: every coverage is within 0..1", profMixed.coverage.every((v) => v >= 0 && v <= 1));
check("mixed: every strength is within 0..1", profMixed.strength.every((v) => v >= 0 && v <= 1));
check("mixed: no NaN in any array",
      [...profMixed.identity, ...profMixed.coverage, ...profMixed.strength]
        .every((v) => !Number.isNaN(v)));

// A mismatch and a gap in the same column: identity 1/2, coverage 2/3,
// strength 1/3.
const profMixed2 = prof("ACD", ["AXD", "A-D", "ACD"]);
check("mismatch plus gap: identity is 1/2", near(profMixed2.identity[1], 0.5, 1e-9),
      `(${profMixed2.identity[1]})`);
check("mismatch plus gap: coverage is 2/3", near(profMixed2.coverage[1], 2 / 3, 1e-9),
      `(${profMixed2.coverage[1]})`);
check("mismatch plus gap: strength is 1/3", near(profMixed2.strength[1], 1 / 3, 1e-9),
      `(${profMixed2.strength[1]})`);
check("mismatch plus gap: strength equals identity times coverage in every column",
      profMixed2.strength.every((v, i) => near(v, profMixed2.identity[i] * profMixed2.coverage[i], 1e-9)));

// conservationScores() keeps its old behaviour: it is the identity series.
const profScores = seq.conservationScores({
  consensus: "ACDEFGH",
  aligned: ["ACDEFGH", "AC-EFGH", "ACD-FGH", "A--EFGH", "ACDEFG-", "ACXEFGH"],
});
const profIdentity = seq.conservationProfile({
  consensus: "ACDEFGH",
  aligned: ["ACDEFGH", "AC-EFGH", "ACD-FGH", "A--EFGH", "ACDEFG-", "ACXEFGH"],
}).identity;
check("conservationScores length matches identity length",
      profScores.length === profIdentity.length, `(${profScores.length} vs ${profIdentity.length})`);
check("conservationScores equals identity element by element",
      profScores.every((v, i) => v === profIdentity[i]),
      `(${JSON.stringify(profScores)})`);
check("conservationScores still scores a gap column 0",
      seq.conservationScores({ consensus: "AC-D", aligned: [["A", "C", "-", "D"], ["A", "C", "-", "D"]] })[2] === 0);
check("conservationScores still scores a sparse column 1",
      seq.conservationScores({ consensus: "AAA", aligned: ["AAA", "AA-"] })[2] === 1);

// A single row with gaps must not throw.
let profGappyThrew = false;
let profGappy = null;
try {
  profGappy = seq.conservationProfile({ consensus: "A--C", aligned: ["A--C"] });
} catch (err) {
  profGappyThrew = true;
}
check("single row with gaps does not throw", !profGappyThrew && profGappy !== null);
check("single row with gaps: columns counted", profGappy && profGappy.columns === 4);
check("single row with gaps: sequences counted", profGappy && profGappy.sequences === 1);
check("single row with gaps: the gap columns are 0",
      profGappy && profGappy.identity[1] === 0 && profGappy.identity[2] === 0 &&
      profGappy.coverage[1] === 0 && profGappy.coverage[2] === 0 &&
      profGappy.strength[1] === 0 && profGappy.strength[2] === 0);
check("single row with gaps: the residue columns are 1",
      profGappy && profGappy.identity[0] === 1 && profGappy.identity[3] === 1 &&
      profGappy.coverage[0] === 1 && profGappy.coverage[3] === 1 &&
      profGappy.strength[0] === 1 && profGappy.strength[3] === 1);

// End to end: an MSA built by progressiveMsa() goes through the profile.
const profBuiltMsa = seq.progressiveMsa(["ACDEFG", "ACDEFG", "ACDXYZ"]);
let profBuilt = null;
let profBuiltThrew = false;
try {
  profBuilt = seq.conservationProfile(profBuiltMsa);
} catch (err) {
  profBuiltThrew = true;
}
check("progressive MSA runs through the profile without throwing",
      !profBuiltThrew && profBuilt !== null);
check("progressive MSA: identity has the consensus length",
      profBuilt && profBuilt.identity.length === profBuiltMsa.consensus.length,
      `(${profBuilt && profBuilt.identity.length} vs ${profBuiltMsa.consensus.length})`);
check("progressive MSA: coverage has the consensus length",
      profBuilt && profBuilt.coverage.length === profBuiltMsa.consensus.length);
check("progressive MSA: strength has the consensus length",
      profBuilt && profBuilt.strength.length === profBuiltMsa.consensus.length);
check("progressive MSA: sequences counted", profBuilt && profBuilt.sequences === 3,
      `(${profBuilt && profBuilt.sequences})`);
check("progressive MSA: every value is within 0..1",
      profBuilt && [...profBuilt.identity, ...profBuilt.coverage, ...profBuilt.strength]
        .every((v) => v >= 0 && v <= 1));
check("progressive MSA: conservationScores equals the identity series",
      profBuilt && seq.conservationScores(profBuiltMsa).every((v, i) => v === profBuilt.identity[i]));

// --------------------------------------------------- reference numbering aids
// parseDbrefs() reads the file's own statement of how its own residue
// numbers relate to the reference database's numbering:
//   DBREF  1M17 A  671   998  UNP    P00533   EGFR_HUMAN     695   1022
// says PDB residues 671-998 are UniProt P00533 residues 695-1022. That
// offset is what makes a conservation score meaningful on a domain
// fragment: a crystal does not start at residue 1 of the protein, so the
// two numberings cannot be assumed to agree.
console.log("\nDBREF parsing");
const dbrefLines = [
  "DBREF  1M17 A  671   998  UNP    P00533   EGFR_HUMAN     695   1022",
  "DBREF  2DN2 A    1   141  UNP    P69905   HBA_HUMAN        1    141",
  "DBREF  2DN2 B    1   146  UNP    P68871   HBB_HUMAN        1    146",
];
const dbrefs = seq.parseDbrefs(dbrefLines.join("\n"));
check("three DBREF records parsed", dbrefs.length === 3, `(${dbrefs.length})`);

const dbref1m17 = dbrefs[0];
check("1M17: chain is A", dbref1m17.chain === "A", `("${dbref1m17.chain}")`);
check("1M17: accession is P00533", dbref1m17.acc === "P00533", `("${dbref1m17.acc}")`);
check("1M17: database is UNP", dbref1m17.db === "UNP", `("${dbref1m17.db}")`);
check("1M17: pdbBeg is 671", dbref1m17.pdbBeg === 671, `(${dbref1m17.pdbBeg})`);
check("1M17: pdbEnd is 998", dbref1m17.pdbEnd === 998, `(${dbref1m17.pdbEnd})`);
check("1M17: dbBeg is 695", dbref1m17.dbBeg === 695, `(${dbref1m17.dbBeg})`);
check("1M17: dbEnd is 1022", dbref1m17.dbEnd === 1022, `(${dbref1m17.dbEnd})`);
check("1M17: the two numberings differ by 24",
      dbref1m17.dbBeg - dbref1m17.pdbBeg === 24,
      `(${dbref1m17.dbBeg - dbref1m17.pdbBeg})`);

const dbref2dn2a = dbrefs[1], dbref2dn2b = dbrefs[2];
check("2DN2: the two records name different chains",
      dbref2dn2a.chain === "A" && dbref2dn2b.chain === "B",
      `(${dbref2dn2a.chain}, ${dbref2dn2b.chain})`);
check("2DN2 chain A: accession is P69905", dbref2dn2a.acc === "P69905",
      `("${dbref2dn2a.acc}")`);
check("2DN2 chain B: accession is P68871", dbref2dn2b.acc === "P68871",
      `("${dbref2dn2b.acc}")`);
check("2DN2 chain A: database is UNP", dbref2dn2a.db === "UNP", `("${dbref2dn2a.db}")`);
check("2DN2 chain B: database is UNP", dbref2dn2b.db === "UNP", `("${dbref2dn2b.db}")`);
check("2DN2 chain A: pdb numbering is 1-141",
      dbref2dn2a.pdbBeg === 1 && dbref2dn2a.pdbEnd === 141,
      `(${dbref2dn2a.pdbBeg}-${dbref2dn2a.pdbEnd})`);
check("2DN2 chain A: db numbering is 1-141",
      dbref2dn2a.dbBeg === 1 && dbref2dn2a.dbEnd === 141,
      `(${dbref2dn2a.dbBeg}-${dbref2dn2a.dbEnd})`);
check("2DN2 chain B: pdb numbering is 1-146",
      dbref2dn2b.pdbBeg === 1 && dbref2dn2b.pdbEnd === 146,
      `(${dbref2dn2b.pdbBeg}-${dbref2dn2b.pdbEnd})`);
check("2DN2 chain B: db numbering is 1-146",
      dbref2dn2b.dbBeg === 1 && dbref2dn2b.dbEnd === 146,
      `(${dbref2dn2b.dbBeg}-${dbref2dn2b.dbEnd})`);

// A file with no DBREF records at all.
const noDbrefs = seq.parseDbrefs([
  "HEADER    OXIDOREDUCTASE",
  "ATOM      1  N   ALA A   1      11.104  13.207  10.000  1.00  0.00           N",
].join("\n"));
check("a file with no DBREF lines gives an empty array", noDbrefs.length === 0,
      `(${noDbrefs.length})`);

// Malformed records are skipped; the good ones survive.
const dbrefMixed = seq.parseDbrefs([
  "HEADER    OXIDOREDUCTASE",
  "DBREF  2DN2 A    1   141  UNP    P69905   HBA_HUMAN        1    141",
  "DBREF  2DN2 A    1   141",
  "DBREF  2DN2 B   XX   146  UNP    P68871   HBB_HUMAN        1    146",
  "DBREF  2DN2 C    1   146  UNP    P68872   HBB_HUMAN       ZZ    146",
  "REMARK 999 DBREF  1M17 A  671   998  UNP    P00533   EGFR_HUMAN     695   1022",
  "ATOM      1  N   ALA A   1      11.104  13.207  10.000  1.00  0.00           N",
].join("\n"));
check("a malformed DBREF line is skipped and the good ones survive",
      dbrefMixed.length === 1 && dbrefMixed[0].acc === "P69905" &&
      dbrefMixed[0].chain === "A",
      `(${dbrefMixed.length}: ${dbrefMixed.map((d) => d.chain + "/" + d.acc).join(", ")})`);
check("a DBREF line with a non-numeric pdbBeg is skipped",
      !dbrefMixed.some((d) => d.chain === "B"));
check("a DBREF line with a non-numeric dbBeg is skipped",
      !dbrefMixed.some((d) => d.chain === "C"));
check("a DBREF mention inside another record is ignored",
      !dbrefMixed.some((d) => d.acc === "P00533"));

// Input that is not a PDB file.
check("a null input gives an empty array", seq.parseDbrefs(null).length === 0);
check("an undefined input gives an empty array", seq.parseDbrefs(undefined).length === 0);
check("a number input gives an empty array", seq.parseDbrefs(42).length === 0);
check("an object input gives an empty array", seq.parseDbrefs({}).length === 0);

// canonicalColumnMap() answers "which alignment column is residue n of the
// reference in?". The alignment grows as sequences are added, so the answer
// is not the identity, and gap columns belong to no residue at all.
console.log("\nReference mapping helpers");
const ccmPlain = seq.canonicalColumnMap(["A", "B", "C"], 3);
check("no gaps: residue i maps to column i",
      ccmPlain.length === 3 && ccmPlain[0] === 0 && ccmPlain[1] === 1 && ccmPlain[2] === 2,
      `(${JSON.stringify(ccmPlain)})`);

const ccmGap = seq.canonicalColumnMap(["A", "-", "B", "C"], 3);
check("an early gap shifts the later residues",
      JSON.stringify(ccmGap) === JSON.stringify([0, 2, 3]), `(${JSON.stringify(ccmGap)})`);

const ccmLead = seq.canonicalColumnMap(["-", "-", "A"], 1);
check("two leading gaps push the residue to column 2",
      ccmLead.length === 1 && ccmLead[0] === 2, `(${JSON.stringify(ccmLead)})`);

const ccmShort = seq.canonicalColumnMap(["A", "B", "C"], 5);
check("a residue beyond the row is reported as -1",
      ccmShort.length === 5 && ccmShort[3] === -1 && ccmShort[4] === -1,
      `(${JSON.stringify(ccmShort)})`);
check("the residues that are in the row are still mapped",
      ccmShort[0] === 0 && ccmShort[1] === 1 && ccmShort[2] === 2,
      `(${JSON.stringify(ccmShort)})`);

const ccmStr = seq.canonicalColumnMap("A-BC", 3);
check("a string row gives the same answer as an array row",
      JSON.stringify(ccmStr) === JSON.stringify(ccmGap), `(${JSON.stringify(ccmStr)})`);

check("length 0 gives an empty array", seq.canonicalColumnMap(["A"], 0).length === 0);
check("a negative length gives an empty array", seq.canonicalColumnMap(["A"], -3).length === 0);

// Every value is -1 or a column index into the row.
const ccmRows = [[["A", "B", "C"], 3], [["A", "-", "B", "C"], 3], [["-", "-", "A"], 1],
                 [["A", "B", "C"], 5], [["A", "-", "-"], 2]];
check("every value is -1 or a valid column index into the row",
      ccmRows.every(([row, len]) => {
        const m = seq.canonicalColumnMap(row, len);
        return m.length === len &&
          m.every((v) => v === -1 || (Number.isInteger(v) && v >= 0 && v < row.length));
      }));
check("gap columns are not attributed to any residue",
      !ccmGap.includes(1) && !ccmLead.includes(0) && !ccmLead.includes(1));
check("a trailing gap does not shift the residues before it",
      JSON.stringify(seq.canonicalColumnMap(["A", "B", "-"], 2)) === JSON.stringify([0, 1]));

// The row shape progressiveMsa() produces is an array of single characters.
const ccmMsa = seq.progressiveMsa(["ACDEF", "ACDEF"]);
const ccmMsaMap = seq.canonicalColumnMap(ccmMsa.aligned[0], 5);
check("a progressiveMsa row maps every residue",
      ccmMsaMap.length === 5 && ccmMsaMap.every((v) => v >= 0), `(${JSON.stringify(ccmMsaMap)})`);
check("the map has one entry per reference residue",
      ccmMsaMap.length === "ACDEF".length, `(${ccmMsaMap.length})`);

// sequenceOffset() finds where a chain's residues come from in the
// reference: the chain is usually a fragment, and the reference may still
// carry an initiator methionine the crystal has already lost.
const offExact = seq.sequenceOffset("ABCDEFG", "ABCDEFG");
check("an exact match has offset 0", offExact.offset === 0, `(${offExact.offset})`);
check("an exact match has identity 1", offExact.identity === 1, `(${offExact.identity})`);

const offFrag = seq.sequenceOffset("MMMMMABCDEFGMMMMM", "ABCDEFG");
check("a fragment in the middle is found", offFrag.offset === 5, `(${offFrag.offset})`);
check("a fragment in the middle has identity 1", offFrag.identity === 1, `(${offFrag.identity})`);

const offMet = seq.sequenceOffset("MVLSPADKT", "VLSPADKT");
check("the lost initiator methionine shifts the offset to 1", offMet.offset === 1,
      `(${offMet.offset})`);
check("the initiator-methionine case has identity 1", offMet.identity === 1,
      `(${offMet.identity})`);

const offMismatch = seq.sequenceOffset("MMMMMABCDEFGMMMMM", "ABCDEFX");
check("a single mismatch still finds the right offset", offMismatch.offset === 5,
      `(${offMismatch.offset})`);
check("a single mismatch lowers the identity below 1",
      offMismatch.identity < 1 && near(offMismatch.identity, 6 / 7, 1e-9),
      `(${offMismatch.identity})`);

// A chain from somewhere else entirely: it must still report its best
// (least-bad) placement rather than nothing at all.
const offUnrelated = seq.sequenceOffset(
  "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGE", "GAMPSKGRQKWVE");
check("an unrelated chain scores well below 1", offUnrelated.identity < 0.5,
      `(${offUnrelated.identity})`);
check("an unrelated chain still reports an offset", offUnrelated.offset >= 0,
      `(${offUnrelated.offset})`);
check("an unrelated chain's offset is the least-bad placement", offUnrelated.offset === 9,
      `(${offUnrelated.offset})`);
check("an unrelated chain's offset stays inside the canonical sequence",
      offUnrelated.offset + "GAMPSKGRQKWVE".length <=
        "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGE".length,
      `(${offUnrelated.offset})`);

// Nothing in common at all: there is no placement to report.
const offNoOverlap = seq.sequenceOffset("MVLSPADKTNVTA", "QQQQQQQQQQQQ");
check("a chain with no common residue scores 0", offNoOverlap.identity === 0,
      `(${offNoOverlap.identity})`);
check("a chain with no common residue scores below 1", offNoOverlap.identity < 1,
      `(${offNoOverlap.identity})`);

// Empty or missing input.
check("an empty canonical gives no offset",
      seq.sequenceOffset("", "ABCDEFG").offset === -1);
check("an empty chain gives no offset",
      seq.sequenceOffset("ABCDEFG", "").offset === -1);
check("an empty canonical gives identity 0",
      seq.sequenceOffset("", "ABCDEFG").identity === 0);
check("an empty chain gives identity 0",
      seq.sequenceOffset("ABCDEFG", "").identity === 0);
const offMissing = seq.sequenceOffset(undefined, undefined);
check("missing arguments give no offset", offMissing.offset === -1, `(${offMissing.offset})`);
check("missing arguments give identity 0", offMissing.identity === 0, `(${offMissing.identity})`);
const offNull = seq.sequenceOffset(null, null);
check("null arguments give no offset", offNull.offset === -1, `(${offNull.offset})`);
check("null arguments give identity 0", offNull.identity === 0, `(${offNull.identity})`);

// Equal scores go to the lowest offset: the scan runs upwards and only
// replaces the best on a strict improvement.
const offTie = seq.sequenceOffset("XXAXXXA", "A");
check("a tie goes to the lowest offset", offTie.offset === 2, `(${offTie.offset})`);
check("the tied offset still has identity 1", offTie.identity === 1, `(${offTie.identity})`);
const offTie2 = seq.sequenceOffset("ABA", "A");
check("a tie with the lowest offset first keeps offset 0", offTie2.offset === 0,
      `(${offTie2.offset})`);

// Every offset that was found must sit inside the canonical sequence.
const offsetFixtures = [
  ["ABCDEFG", "ABCDEFG"],
  ["MMMMMABCDEFGMMMMM", "ABCDEFG"],
  ["MVLSPADKT", "VLSPADKT"],
  ["MMMMMABCDEFGMMMMM", "ABCDEFX"],
  ["MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGE", "GAMPSKGRQKWVE"],
  ["ABA", "A"],
  ["XXAXXXA", "A"],
];
check("every offset found is inside the canonical sequence",
      offsetFixtures.every(([canon, chain]) => {
        const r = seq.sequenceOffset(canon, chain);
        return r.offset >= 0 && r.offset < canon.length && r.offset + chain.length <= canon.length;
      }));
check("every identity is a fraction between 0 and 1",
      offsetFixtures.every(([canon, chain]) => {
        const r = seq.sequenceOffset(canon, chain);
        return r.identity >= 0 && r.identity <= 1;
      }));
check("a chain longer than the canonical is aligned at offset 0",
      seq.sequenceOffset("ABC", "ABCDEFG").offset === 0,
      `(${seq.sequenceOffset("ABC", "ABCDEFG").offset})`);

console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
