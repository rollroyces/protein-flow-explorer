/* ===================================================================
 * Protein Flow Explorer — pocket-detection geometry module
 *
 * Pure spatial-hash geometry behind the pocket detector: the atom
 * hash, the density-gate neighbour count, ray casting and the enclosure
 * fraction. The detector itself (thresholds, clustering, UI) stays in app.js.
 *
 * Extracted from app.js during the 2026 architecture audit. The module is
 * loaded as a classic <script> before app.js and publishes itself on the
 * single PFE namespace; app.js destructures the members it needs, so no
 * call site changed when the code moved.
 * =================================================================== */
(function (global) {
  "use strict";

  const POCKET_CELL = 10.0;      // spatial-hash cell for the density gates

  const POCKET_RAY_R = 9.0;

  const POCKET_RAY_N = 32;

  const POCKET_DIRS = (function () {
    const out = [];
    const ga = Math.PI * (3 - Math.sqrt(5));
    for (let i = 0; i < POCKET_RAY_N; i++) {
      const y = 1 - (i / (POCKET_RAY_N - 1)) * 2;
      const r = Math.sqrt(Math.max(0, 1 - y * y));
      const th = ga * i;
      out.push([Math.cos(th) * r, y, Math.sin(th) * r]);
    }
    return out;
  })();

  function buildAtomHash(atoms) {
    const map = new Map();
    for (const a of atoms) {
      const k = `${Math.floor(a.x / POCKET_CELL)}|${Math.floor(a.y / POCKET_CELL)}|${Math.floor(a.z / POCKET_CELL)}`;
      let arr = map.get(k);
      if (!arr) { arr = []; map.set(k, arr); }
      arr.push(a);
    }
    return map;
  }

  function countAtomsNear(hash, x, y, z, radius) {
    const cx = Math.floor(x / POCKET_CELL), cy = Math.floor(y / POCKET_CELL), cz = Math.floor(z / POCKET_CELL);
    const span = Math.ceil(radius / POCKET_CELL);
    const r2 = radius * radius;
    let c = 0;
    for (let dx = -span; dx <= span; dx++) {
      for (let dy = -span; dy <= span; dy++) {
        for (let dz = -span; dz <= span; dz++) {
          const arr = hash.get(`${cx + dx}|${cy + dy}|${cz + dz}`);
          if (!arr) continue;
          for (let i = 0; i < arr.length; i++) {
            const a = arr[i];
            const ax = a.x - x, ay = a.y - y, az = a.z - z;
            if (ax * ax + ay * ay + az * az < r2) c++;
          }
        }
      }
    }
    return c;
  }

  function rayHits(hash, x, y, z, dx, dy, dz, radius) {
    const steps = Math.ceil(radius / 1.0);
    const hit2 = 2.6 * 2.6;
    for (let s = 1; s <= steps; s++) {
      const t = s * 1.0;
      const px = x + dx * t, py = y + dy * t, pz = z + dz * t;
      const cx = Math.floor(px / POCKET_CELL), cy = Math.floor(py / POCKET_CELL), cz = Math.floor(pz / POCKET_CELL);
      for (let ax2 = cx - 1; ax2 <= cx + 1; ax2++) {
        for (let ay2 = cy - 1; ay2 <= cy + 1; ay2++) {
          for (let az2 = cz - 1; az2 <= cz + 1; az2++) {
            const arr = hash.get(`${ax2}|${ay2}|${az2}`);
            if (!arr) continue;
            for (let i = 0; i < arr.length; i++) {
              const a = arr[i];
              const qx = a.x - px, qy = a.y - py, qz = a.z - pz;
              if (qx * qx + qy * qy + qz * qz < hit2) return true;
            }
          }
        }
      }
    }
    return false;
  }

  function enclosureFraction(hash, x, y, z) {
    let hits = 0;
    for (const d of POCKET_DIRS) {
      if (rayHits(hash, x, y, z, d[0], d[1], d[2], POCKET_RAY_R)) hits++;
    }
    return hits / POCKET_RAY_N;
  }

  const api = {
    buildAtomHash,
    countAtomsNear,
    rayHits,
    enclosureFraction,
  };

  global.PFE = global.PFE || {};
  global.PFE.pockets = Object.freeze(api);
})(typeof globalThis !== "undefined" ? globalThis : this);
