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

console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
