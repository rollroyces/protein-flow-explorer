/* ===================================================================
 * Protein Flow Explorer — molecular geometry module
 *
 * Pure vector maths for the measurement tool: distances, angles,
 * dihedrals and per-residue RMSD. No DOM, no network, no app state.
 *
 * Extracted from app.js during the 2026 architecture audit. The module is
 * loaded as a classic <script> before app.js and publishes itself on the
 * single PFE namespace; app.js destructures the members it needs, so no
 * call site changed when the code moved.
 * =================================================================== */
(function (global) {
  "use strict";

  function midPoint(a, b) {
    return { x: (a.x + b.x) / 2, y: (a.y + b.y) / 2, z: (a.z + b.z) / 2 };
  }

  function cross(a, b) {
    return {
      x: a.y * b.z - a.z * b.y,
      y: a.z * b.x - a.x * b.z,
      z: a.x * b.y - a.y * b.x,
    };
  }

  function dot(a, b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

  function computeDistance(a, b) {
    return Math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2);
  }

  function computeAngle(a, b, c) {
    // angle at b between ba and bc vectors
    const v1 = { x: a.x - b.x, y: a.y - b.y, z: a.z - b.z };
    const v2 = { x: c.x - b.x, y: c.y - b.y, z: c.z - b.z };
    const dot = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    const m1 = Math.sqrt(v1.x ** 2 + v1.y ** 2 + v1.z ** 2);
    const m2 = Math.sqrt(v2.x ** 2 + v2.y ** 2 + v2.z ** 2);
    const cos = Math.max(-1, Math.min(1, dot / (m1 * m2)));
    return Math.acos(cos) * 180 / Math.PI;
  }

  function computeDihedral(a, b, c, d) {
    // 4-atom dihedral using vectors b1=b-a, b2=c-b, b3=d-c
    const b1 = { x: b.x - a.x, y: b.y - a.y, z: b.z - a.z };
    const b2 = { x: c.x - b.x, y: c.y - b.y, z: c.z - b.z };
    const b3 = { x: d.x - c.x, y: d.y - c.y, z: d.z - c.z };
    // n1 = b1 × b2 ; n2 = b2 × b3
    const n1 = cross(b1, b2);
    const n2 = cross(b2, b3);
    const m1 = Math.sqrt(dot(n1, n1));
    // n2 is deliberately left unnormalised: atan2 is invariant under a
    // positive scaling of both arguments, so normalising it changes nothing.
    const b2n = Math.sqrt(dot(b2, b2));
    // m = b2 / |b2|
    const m = { x: b2.x / b2n, y: b2.y / b2n, z: b2.z / b2n };
    // x = n1 / |n1|
    const x = { x: n1.x / m1, y: n1.y / m1, z: n1.z / m1 };
    // y = m × x
    const y = cross(m, x);
    // angle = atan2( n2·y , n2·x )
    const ang = Math.atan2(dot(n2, y), dot(n2, x));
    return ang * 180 / Math.PI;
  }

  function computePerResidueRmsd(modelA, modelB) {
    const aAtoms = modelA.selectedAtoms({ atom: "CA" });
    const bAtoms = modelB.selectedAtoms({ atom: "CA" });
    // Build lookup for B atoms by (chain, resi)
    const bMap = new Map();
    for (const a of bAtoms) {
      bMap.set(`${a.chain}:${a.resi}`, a);
    }
    const result = [];
    for (const a of aAtoms) {
      const b = bMap.get(`${a.chain}:${a.resi}`);
      if (!b) continue;
      const dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
      const d = Math.sqrt(dx * dx + dy * dy + dz * dz);
      result.push({ chain: a.chain, resi: a.resi, value: d });
    }
    return result;
  }

  const api = {
    midPoint,
    cross,
    dot,
    computeDistance,
    computeAngle,
    computeDihedral,
    computePerResidueRmsd,
  };

  global.PFE = global.PFE || {};
  global.PFE.geometry = Object.freeze(api);
})(typeof globalThis !== "undefined" ? globalThis : this);
