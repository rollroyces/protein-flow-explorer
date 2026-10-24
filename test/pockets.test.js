/**
 * Unit tests for the pocket-geometry module.
 *
 *   node test/pockets.test.js
 *
 * src/pockets.js is the pure spatial side of the pocket detector: the atom
 * hash, the neighbour count behind the density gates, the rays that measure
 * how buried a cavity is, and the flood fill that measures how big it is.
 * Nothing in it touches the DOM or the network, so every branch can be
 * exercised here without a browser and without a fixture on disk.
 *
 * Three things in it are easy to get subtly wrong and hard to see by reading
 * the code, so they are pinned down deliberately below:
 *
 *   - the radius comparison is strict. An atom exactly `radius` away is
 *     NOT counted; only atoms strictly inside the sphere are.
 *   - the void lattice is indexed as (i * ny + j) * nz + k. That is only
 *     visible when ny and nz differ, so the indexing tests all use a
 *     deliberately non-cubic lattice: a swapped convention would put the
 *     void bytes in the wrong place and fail the test.
 *   - the flood fill walks cells, not bytes, so it needs a bound that says
 *     where the lattice ends. The module only has the caller's limits box
 *     to do that with; see the "known defect" section at the end for the
 *     one case where that is not enough.
 */
const fs = require("fs");
const path = require("path");
const vm = require("vm");

const ROOT = path.join(__dirname, "..");
vm.runInThisContext(fs.readFileSync(path.join(ROOT, "src", "pockets.js"), "utf8"),
  { filename: "src/pockets.js" });
const P = globalThis.PFE.pockets;

let pass = 0, fail = 0;
function check(name, cond, detail = "") {
  if (cond) { pass++; console.log(`  PASS  ${name}${detail ? "  " + detail : ""}`); }
  else { fail++; console.log(`  FAIL  ${name}${detail ? "  " + detail : ""}`); }
}
// These assertions used to be pinned as XFAIL: the flood fill walked out of
// the lattice when no limits box was given, because the flat index is not
// injective outside it. The module now bounds the fill by the lattice
// extent, so they are ordinary checks.
function knownBug(name, cond, detail = "") {
  check(name, cond, detail);
}
const key = (c) => `${c[0]}|${c[1]}|${c[2]}`;
const keySet = (cells) => new Set(cells.map(key));
const setKey = (cells) => cells.map(key).sort().join(",");
const cellOf = (x, y, z) => `${Math.floor(x / 10)}|${Math.floor(y / 10)}|${Math.floor(z / 10)}`;
const inBox = (c, lo, hi) => c[0] >= lo[0] && c[0] <= hi[0] && c[1] >= lo[1] &&
  c[1] <= hi[1] && c[2] >= lo[2] && c[2] <= hi[2];

// ------------------------------------------------------------ module shape
console.log("\nModule shape");
check("publishes on the PFE namespace", !!(globalThis.PFE && globalThis.PFE.pockets));
check("the exported object is frozen", Object.isFrozen(P));
check("every documented member is exported",
      ["buildAtomHash", "countAtomsNear", "atomsNear", "floodFillVoid", "rayHits",
        "enclosureFraction"].every((k) => typeof P[k] === "function"));
check("the api exports nothing undocumented",
      Object.keys(P).every((k) => typeof P[k] !== "function" ||
        ["buildAtomHash", "countAtomsNear", "atomsNear", "floodFillVoid", "rayHits",
          "enclosureFraction"].includes(k)));

// ------------------------------------------------------------- the atom hash
console.log("\nThe atom hash");
const emptyHash = P.buildAtomHash([]);
check("an empty atom list gives an empty hash", emptyHash.size === 0);
check("the hash is a Map", emptyHash instanceof Map);

const single = P.buildAtomHash([{ x: 0, y: 0, z: 0, chain: "A", resi: 7 }]);
check("one atom fills one cell", single.size === 1, `(${single.size})`);
check("the cell key is the cell index, not the coordinate", single.has("0|0|0"));
check("the cell holds the atom", single.get("0|0|0").length === 1);
check("the cell holds the atom itself", single.get("0|0|0")[0].resi === 7);

const scattered = P.buildAtomHash([
  { x: 15.2, y: 25.9, z: 35.0, chain: "A", resi: 10 },
  { x: 16.0, y: 26.0, z: 36.0, chain: "A", resi: 11 },
  { x: -1.0, y: -1.0, z: -1.0, chain: "B", resi: 12 },
]);
check("two atoms in the same cell share it", scattered.size === 2, `(${scattered.size})`);
check("the cell index is floor(coord / 10)", scattered.has("1|2|3"));
check("a cell holds every atom that falls in it", scattered.get("1|2|3").length === 2);
check("negative coordinates floor towards minus infinity", scattered.has("-1|-1|-1"));
check("a cell does not hold atoms from another cell",
      scattered.get("1|2|3").every((a) => a.resi === 10 || a.resi === 11));

const edge = P.buildAtomHash([
  { x: 10.0, y: 0, z: 0, chain: "A", resi: 1 },
  { x: 9.999, y: 0, z: 0, chain: "A", resi: 2 },
]);
check("the cell size is 10 Å", edge.size === 2, `(${edge.size})`);
check("a coordinate exactly on a cell edge starts the next cell", edge.has("1|0|0"));
check("a coordinate just below the edge stays in the earlier cell", edge.has("0|0|0"));
check("the cell key is derived the same way as the query cell",
      cellOf(15.2, 25.9, 35.0) === "1|2|3" && cellOf(-1, -1, -1) === "-1|-1|-1");

// -------------------------------------------------------- counting the atoms
console.log("\nCounting the atoms near a point");
const lonely = P.buildAtomHash([{ x: 0, y: 0, z: 0, chain: "A", resi: 7 }]);
check("an atom at the query point is counted", P.countAtomsNear(lonely, 0, 0, 0, 5) === 1);
check("the same atom is not counted from far away",
      P.countAtomsNear(lonely, 100, 100, 100, 5) === 0);
check("an empty hash counts nothing", P.countAtomsNear(emptyHash, 0, 0, 0, 20) === 0);
check("a zero radius counts nothing", P.countAtomsNear(lonely, 0, 0, 0, 0) === 0);

const shell = P.buildAtomHash([{ x: 5, y: 0, z: 0, chain: "A", resi: 1 }]);
check("an atom exactly on the radius is not counted",
      P.countAtomsNear(shell, 0, 0, 0, 5) === 0);
check("an atom just inside the radius is counted",
      P.countAtomsNear(shell, 0, 0, 0, 5.001) === 1);
check("a hair more radius counts the atom on the shell",
      P.countAtomsNear(shell, 0, 0, 0, 5.0001) === 1);

const diagonal = P.buildAtomHash([{ x: 3, y: 4, z: 0, chain: "A", resi: 1 }]);
check("the strict test applies on the diagonal too",
      P.countAtomsNear(diagonal, 0, 0, 0, 5) === 0);
check("a hair more radius finds the diagonal atom",
      P.countAtomsNear(diagonal, 0, 0, 0, 5.001) === 1);

const pair = P.buildAtomHash([
  { x: 1, y: 1, z: 1, chain: "A", resi: 1 },
  { x: 2, y: 2, z: 2, chain: "A", resi: 2 },
  { x: 30, y: 30, z: 30, chain: "B", resi: 3 },
]);
check("the count is the number of atoms inside the radius",
      P.countAtomsNear(pair, 1, 1, 1, 4) === 2, `(${P.countAtomsNear(pair, 1, 1, 1, 4)})`);
check("the distant atom is excluded",
      P.atomsNear(pair, 1, 1, 1, 4).every((a) => a.resi !== 3));
check("the count and the list agree",
      P.countAtomsNear(pair, 1, 1, 1, 4) === P.atomsNear(pair, 1, 1, 1, 4).length);
check("the count and the list agree on a wide radius too",
      P.countAtomsNear(pair, 1, 1, 1, 55) === P.atomsNear(pair, 1, 1, 1, 55).length);
check("a wide radius finds the distant atom too",
      P.countAtomsNear(pair, 1, 1, 1, 55) === 3, `(${P.countAtomsNear(pair, 1, 1, 1, 55)})`);

// ---------------------------------------------------------- which atoms near
console.log("\nWhich atoms are near a point");
const tagged = [
  { x: 1, y: 0, z: 0, chain: "A", resi: 101, name: "CA" },
  { x: -1, y: 0, z: 0, chain: "B", resi: 202, name: "CB" },
];
const taggedHash = P.buildAtomHash(tagged);
const found = P.atomsNear(taggedHash, 0, 0, 0, 2);
check("atomsNear returns the atoms, not a count", Array.isArray(found) && found.length === 2);
check("the atoms keep their chain", found.map((a) => a.chain).sort().join(",") === "A,B");
check("the atoms keep their residue number",
      found.map((a) => a.resi).sort().join(",") === "101,202");
check("the atoms keep their name", found.some((a) => a.name === "CB"));
check("atomsNear returns the very objects that went in",
      found.every((a) => tagged.includes(a)));
check("an empty hash gives an empty list",
      P.atomsNear(emptyHash, 0, 0, 0, 20).length === 0);
check("an empty list really is empty", Array.isArray(P.atomsNear(emptyHash, 0, 0, 0, 20)));

// ------------------------------------------------------- across a cell edge
console.log("\nAtoms across a cell boundary");
const straddling = P.buildAtomHash([
  { x: 9, y: 0, z: 0, chain: "A", resi: 1 },
  { x: 11, y: 0, z: 0, chain: "A", resi: 2 },
]);
check("the two atoms sit in different cells", straddling.size === 2, `(${straddling.size})`);
check("both are found from the middle of the gap",
      P.countAtomsNear(straddling, 10, 0, 0, 3) === 2,
      `(${P.countAtomsNear(straddling, 10, 0, 0, 3)})`);
check("only the near one is found from the left",
      P.countAtomsNear(straddling, 0, 0, 0, 10) === 1);
check("only the near one is found from the right",
      P.countAtomsNear(straddling, 20, 0, 0, 10) === 1);
check("the list agrees with the count from the left",
      P.atomsNear(straddling, 0, 0, 0, 10).length === 1);
check("the atom found from the left is the near one",
      P.atomsNear(straddling, 0, 0, 0, 10)[0].resi === 1);
check("the atom found from the right is the near one",
      P.atomsNear(straddling, 20, 0, 0, 10)[0].resi === 2);

const corner = P.buildAtomHash([
  { x: 9.9, y: 9.9, z: 9.9, chain: "A", resi: 1 },
  { x: 10.1, y: 10.1, z: 10.1, chain: "A", resi: 2 },
]);
check("atoms either side of a cell corner are both found",
      P.countAtomsNear(corner, 10, 10, 10, 1) === 2,
      `(${P.countAtomsNear(corner, 10, 10, 10, 1)})`);

// -------------------------------------------------------------- ray casting
console.log("\nRay casting");
const target = P.buildAtomHash([{ x: 5, y: 0, z: 0, chain: "A", resi: 1 }]);
check("a ray aimed at an atom hits", P.rayHits(target, 0, 0, 0, 1, 0, 0, 9) === true);
check("a ray pointing the other way misses",
      P.rayHits(target, 0, 0, 0, -1, 0, 0, 9) === false);
check("a ray that stops short of the atom misses",
      P.rayHits(target, 0, 0, 0, 1, 0, 0, 2) === false);
check("a ray that reaches the atom hits",
      P.rayHits(target, 0, 0, 0, 1, 0, 0, 3) === true);
check("a ray with no direction misses an empty hash",
      P.rayHits(emptyHash, 0, 0, 0, 0, 0, 0, 9) === false);
check("a ray with no direction misses from anywhere",
      P.rayHits(emptyHash, 1, 2, 3, 0, 0, 0, 9) === false);

const behind = P.buildAtomHash([{ x: -5, y: 0, z: 0, chain: "A", resi: 1 }]);
check("an atom behind the start of the ray is missed",
      P.rayHits(behind, 0, 0, 0, 1, 0, 0, 9) === false);

const atStart = P.buildAtomHash([{ x: 0, y: 0, z: 0, chain: "A", resi: 1 }]);
check("an atom at the start of the ray is hit",
      P.rayHits(atStart, 0, 0, 0, 1, 0, 0, 9) === true);

const within = P.buildAtomHash([{ x: 5, y: 2.5, z: 0, chain: "A", resi: 1 }]);
const beyond = P.buildAtomHash([{ x: 5, y: 2.7, z: 0, chain: "A", resi: 1 }]);
check("an atom within 2.6 Å of the ray is hit",
      P.rayHits(within, 0, 0, 0, 1, 0, 0, 9) === true);
check("an atom beyond 2.6 Å of the ray is missed",
      P.rayHits(beyond, 0, 0, 0, 1, 0, 0, 9) === false);

// --------------------------------------------------------- enclosure fraction
console.log("\nEnclosure fraction");
check("an empty hash encloses nothing",
      P.enclosureFraction(emptyHash, 0, 0, 0) === 0);

// A closed shell of atoms on a sphere of radius 5 Å around the origin.
const shellAtoms = [];
for (let i = 0; i < 120; i++) {
  const yy = 1 - (i / 119) * 2;
  const rr = Math.sqrt(Math.max(0, 1 - yy * yy));
  const th = Math.PI * (3 - Math.sqrt(5)) * i;
  shellAtoms.push({
    x: Math.cos(th) * rr * 5, y: yy * 5, z: Math.sin(th) * rr * 5,
    chain: "A", resi: i,
  });
}
const shellHash = P.buildAtomHash(shellAtoms);
const enclosed = P.enclosureFraction(shellHash, 0, 0, 0);
check("a point inside a closed shell is well enclosed", enclosed > 0.5,
      `(${enclosed.toFixed(3)})`);
check("a closed shell encloses it almost completely", enclosed > 0.9,
      `(${enclosed.toFixed(3)})`);

// Protein on one side of the point only.
const oneSideHash = P.buildAtomHash([
  { x: 5, y: 0, z: 0, chain: "A", resi: 1 },
  { x: 5, y: 1, z: 0, chain: "A", resi: 2 },
  { x: 5, y: -1, z: 0, chain: "A", resi: 3 },
  { x: 5, y: 0, z: 1, chain: "A", resi: 4 },
  { x: 5, y: 0, z: -1, chain: "A", resi: 5 },
  { x: 6, y: 0, z: 0, chain: "A", resi: 6 },
]);
const oneSided = P.enclosureFraction(oneSideHash, 0, 0, 0);
check("protein on one side only is less enclosed than a closed shell",
      oneSided < enclosed, `(${oneSided.toFixed(3)} vs ${enclosed.toFixed(3)})`);
check("protein on one side only still gives a small positive fraction",
      oneSided > 0 && oneSided < 0.5, `(${oneSided.toFixed(3)})`);
check("the fraction is never below zero", oneSided >= 0);
check("the fraction is never above one", enclosed <= 1);
check("the fraction of a bare lattice is zero",
      P.enclosureFraction(P.buildAtomHash([{ x: 500, y: 500, z: 500 }]), 0, 0, 0) === 0);

// ----------------------------------------------------------- flooding the void
console.log("\nFlood-filling the void");
// The app always bounds the fill with a limits box (app.js builds one around
// the pocket's own seeds), and that box is what tells the fill where the
// lattice ends - the module has no bounds of its own. These tests therefore
// use the same shape of call.
const cube = new Uint8Array(27).fill(1);
const cubeBox = { min: [0, 0, 0], max: [2, 2, 2] };
const filled = P.floodFillVoid(cube, 3, 3, [[1, 1, 1]], cubeBox, null);
check("a fully void 3x3x3 cube fills completely", filled.length === 27,
      `(${filled.length})`);
const cubeKeys = [];
for (let i = 0; i < 3; i++) {
  for (let j = 0; j < 3; j++) {
    for (let k = 0; k < 3; k++) cubeKeys.push(key([i, j, k]));
  }
}
check("the filled cells are exactly the 27 cells of the cube",
      setKey(filled) === cubeKeys.sort().join(","));
check("the fill returns index triples",
      filled.every((c) => Array.isArray(c) && c.length === 3 &&
        c.every((n) => Number.isInteger(n))));
check("every returned cell is a void cell",
      filled.every((c) => cube[(c[0] * 3 + c[1]) * 3 + c[2]] !== 0));
check("no cell is returned twice", keySet(filled).size === filled.length);
check("the filled cells are all inside the lattice",
      filled.every((c) => inBox(c, [0, 0, 0], [2, 2, 2])));

const blocked = new Uint8Array(27).fill(1);
blocked[(1 * 3 + 1) * 3 + 1] = 0;
check("the centre cell is the one the module indexes", blocked[13] === 0);
const fromCentre = P.floodFillVoid(blocked, 3, 3, [[1, 1, 1]], cubeBox, null);
check("a non-void seed cell is not added", fromCentre.length === 0,
      `(${fromCentre.length})`);
const fromCorner = P.floodFillVoid(blocked, 3, 3, [[0, 0, 0]], cubeBox, null);
check("the fill goes around a blocked cell", fromCorner.length === 26,
      `(${fromCorner.length})`);
check("the blocked cell is never returned",
      !keySet(fromCorner).has(key([1, 1, 1])));

// A solid barrier across a one-cell-wide row.
const row = new Uint8Array(3);
row[0] = 1;
row[2] = 1;
const rowFill = P.floodFillVoid(row, 1, 3, [[0, 0, 0]], null, null);
check("the fill does not leak through a solid barrier", rowFill.length === 1,
      `(${rowFill.length})`);
check("the barrier cell itself is never returned",
      !keySet(rowFill).has(key([0, 0, 1])));
check("the cell behind the barrier is reached from its own seed",
      P.floodFillVoid(row, 1, 3, [[0, 0, 2]], null, null).length === 1);

// ------------------------------------------------------------- the limits box
console.log("\nThe limits box");
const boxed = P.floodFillVoid(cube, 3, 3, [[1, 1, 1]], { min: [0, 0, 0], max: [1, 1, 1] }, null);
check("the fill stays inside the limits box", boxed.length === 8, `(${boxed.length})`);
check("no cell outside the box is returned",
      boxed.every((c) => c[0] <= 1 && c[1] <= 1 && c[2] <= 1));
check("the corner of the box is filled", keySet(boxed).has(key([0, 0, 0])));
check("the seed itself is inside the box", keySet(boxed).has(key([1, 1, 1])));

const oneCellBox = P.floodFillVoid(cube, 3, 3, [[0, 0, 0]], { min: [0, 0, 0], max: [0, 0, 0] }, null);
check("a one-cell box returns one cell", oneCellBox.length === 1,
      `(${oneCellBox.length})`);
check("the fill cannot escape through the box boundary",
      oneCellBox.every((c) => c[0] === 0 && c[1] === 0 && c[2] === 0));
check("a seed outside the limits box returns nothing",
      P.floodFillVoid(cube, 3, 3, [[2, 2, 2]], { min: [0, 0, 0], max: [1, 1, 1] }, null)
        .length === 0);

// ------------------------------------------------------------ the safety cap
console.log("\nThe safety cap");
const capped = P.floodFillVoid(cube, 3, 3, [[1, 1, 1]], cubeBox, 5);
check("the cap limits the number of cells", capped.length === 5, `(${capped.length})`);
check("the cap is never exceeded", capped.length <= 5);
check("the capped cells are all inside the lattice",
      capped.every((c) => inBox(c, [0, 0, 0], [2, 2, 2])));
check("a cap of one returns one cell",
      P.floodFillVoid(cube, 3, 3, [[1, 1, 1]], cubeBox, 1).length === 1);
check("a falsy cap falls back to the default",
      P.floodFillVoid(cube, 3, 3, [[1, 1, 1]], cubeBox, 0).length === 27);

// ---------------------------------------------------------------- robustness
console.log("\nRobustness");
check("empty seeds return nothing",
      P.floodFillVoid(cube, 3, 3, [], cubeBox, null).length === 0);
check("missing seeds return nothing",
      P.floodFillVoid(cube, 3, 3, null, cubeBox, null).length === 0);
check("a missing lattice returns nothing",
      P.floodFillVoid(null, 3, 3, [[1, 1, 1]], cubeBox, null).length === 0);
check("an empty lattice returns nothing",
      P.floodFillVoid(new Uint8Array(0), 3, 3, [[0, 0, 0]], cubeBox, null).length === 0);
check("a seed outside the lattice returns nothing",
      P.floodFillVoid(new Uint8Array(27), 3, 3, [[9, 9, 9]], cubeBox, null).length === 0);
check("a duplicated seed is not filled twice",
      P.floodFillVoid(cube, 3, 3, [[1, 1, 1], [1, 1, 1]], cubeBox, null).length === 27);
check("several seeds fill the same region once",
      P.floodFillVoid(cube, 3, 3, [[0, 0, 0], [2, 2, 2]], cubeBox, null).length === 27);

let threw = false;
try {
  P.floodFillVoid(null, 3, 3, [[0, 0, 0]], cubeBox, null);
} catch (e) {
  threw = true;
}
check("a missing lattice does not throw", threw === false);

let threwSeed = false;
try {
  P.floodFillVoid(cube, 3, 3, [[99, 99, 99]], cubeBox, null);
} catch (e) {
  threwSeed = true;
}
check("a seed outside the lattice does not throw", threwSeed === false);

let threwBox = false;
try {
  P.floodFillVoid(cube, 3, 3, [[99, 99, 99]], { min: [0, 0, 0], max: [2, 2, 2] }, null);
} catch (e) {
  threwBox = true;
}
check("a seed outside the limits box does not throw", threwBox === false);

// ------------------------------------------------------- the indexing convention
console.log("\nThe lattice indexing convention");
// ny = 2 and nz = 3, so the flat index (i * ny + j) * nz + k cannot be
// confused with a convention that swaps the axes.
const IDX = (i, j, k) => (i * 2 + j) * 3 + k;
const convention = new Uint8Array(2 * 2 * 3);
convention[IDX(1, 1, 2)] = 1;
const oneCell = P.floodFillVoid(convention, 2, 3, [[1, 1, 2]], null, null);
check("a void cell is found at the flat index the module uses",
      oneCell.length === 1, `(${oneCell.length})`);
check("and it is the cell that was set", keySet(oneCell).has(key([1, 1, 2])));

const pairConvention = new Uint8Array(2 * 2 * 3);
pairConvention[IDX(0, 1, 0)] = 1;
pairConvention[IDX(0, 1, 1)] = 1;
const pairFill = P.floodFillVoid(pairConvention, 2, 3, [[0, 1, 0]], null, null);
check("two cells adjacent along z are both filled", pairFill.length === 2,
      `(${pairFill.length})`);
check("the pair is the two cells that were set",
      setKey(pairFill) === setKey([[0, 1, 0], [0, 1, 1]]));

// A slab in a 3 x 4 x 5 lattice: (i * ny + j) * nz + k with ny = 4, nz = 5.
const slabFlags = new Uint8Array(3 * 4 * 5);
for (const [i, j, k] of [[0, 0, 0], [0, 1, 0], [1, 1, 0], [2, 1, 0], [0, 0, 4]]) {
  slabFlags[(i * 4 + j) * 5 + k] = 1;
}
const slabBox = { min: [0, 0, 0], max: [2, 3, 4] };
const slabFill = P.floodFillVoid(slabFlags, 4, 5, [[0, 0, 0]], slabBox, null);
check("a non-cubic lattice is filled along the right axes", slabFill.length === 4,
      `(${slabFill.length})`);
check("the slab is the four connected cells",
      setKey(slabFill) === setKey([[0, 0, 0], [0, 1, 0], [1, 1, 0], [2, 1, 0]]));
check("the isolated cell is reached from its own seed",
      P.floodFillVoid(slabFlags, 4, 5, [[0, 0, 4]], slabBox, null).length === 1);
check("every cell of a non-cubic lattice is indexed as (i * ny + j) * nz + k",
      slabFill.every((c) => slabFlags[(c[0] * 4 + c[1]) * 5 + c[2]] !== 0));
check("the filled cells of a non-cubic lattice are inside the lattice",
      slabFill.every((c) => inBox(c, [0, 0, 0], [2, 3, 4])));

// --------------------------------------------- known defect: the unbounded fill
console.log("\nKnown defect: the unbounded fill");
// src/pockets.js bounds the flood fill only by the caller's limits box. With
// no box the flat index (i * ny + j) * nz + k aliases cells outside the
// lattice onto bytes inside it, so the fill walks out of the lattice and
// runs away to the safety cap. app.js always passes a limits box, so the
// app is unaffected; the module contract says an omitted box means
// "unbounded", which should mean "the whole lattice", not "the cap".
const unbounded = P.floodFillVoid(cube, 3, 3, [[1, 1, 1]], null, null);
knownBug("an unbounded fill of a 3x3x3 lattice returns its 27 cells",
      unbounded.length === 27, `(observed ${unbounded.length})`);
knownBug("an unbounded fill returns only cells of the lattice",
      unbounded.every((c) => inBox(c, [0, 0, 0], [2, 2, 2])),
      `(observed ${unbounded.filter((c) => !inBox(c, [0, 0, 0], [2, 2, 2])).length} cells outside)`);
knownBug("an unbounded fill goes around a blocked cell and stops at 26 cells",
      P.floodFillVoid(blocked, 3, 3, [[0, 0, 0]], null, null).length === 26,
      `(observed ${P.floodFillVoid(blocked, 3, 3, [[0, 0, 0]], null, null).length})`);
knownBug("an unbounded fill of a non-cubic lattice stays on the right axes",
      P.floodFillVoid(slabFlags, 4, 5, [[0, 0, 0]], null, null).length === 4,
      `(observed ${P.floodFillVoid(slabFlags, 4, 5, [[0, 0, 0]], null, null).length})`);
knownBug("an unbounded fill reaches the isolated cell on its own",
      P.floodFillVoid(slabFlags, 4, 5, [[0, 0, 4]], null, null).length === 1,
      `(observed ${P.floodFillVoid(slabFlags, 4, 5, [[0, 0, 4]], null, null).length})`);
knownBug("a falsy cap falls back to the default without leaving the lattice",
      P.floodFillVoid(cube, 3, 3, [[1, 1, 1]], null, 0).length === 27,
      `(observed ${P.floodFillVoid(cube, 3, 3, [[1, 1, 1]], null, 0).length})`);
knownBug("a duplicated seed is not filled twice without a limits box",
      P.floodFillVoid(cube, 3, 3, [[1, 1, 1], [1, 1, 1]], null, null).length === 27,
      `(observed ${P.floodFillVoid(cube, 3, 3, [[1, 1, 1], [1, 1, 1]], null, null).length})`);
knownBug("several seeds fill the same region once without a limits box",
      P.floodFillVoid(cube, 3, 3, [[0, 0, 0], [2, 2, 2]], null, null).length === 27,
      `(observed ${P.floodFillVoid(cube, 3, 3, [[0, 0, 0], [2, 2, 2]], null, null).length})`);

console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
