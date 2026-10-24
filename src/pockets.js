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

  /**
   * The atoms within `radius` of a point.
   *
   * countAtomsNear answers "how many"; this answers "which", which the
   * pocket detector needs to name the residues that line a cavity.
   */
  function atomsNear(hash, x, y, z, radius) {
    const cx = Math.floor(x / POCKET_CELL), cy = Math.floor(y / POCKET_CELL), cz = Math.floor(z / POCKET_CELL);
    const span = Math.ceil(radius / POCKET_CELL);
    const r2 = radius * radius;
    const out = [];
    for (let dx = -span; dx <= span; dx++) {
      for (let dy = -span; dy <= span; dy++) {
        for (let dz = -span; dz <= span; dz++) {
          const arr = hash.get(`${cx + dx}|${cy + dy}|${cz + dz}`);
          if (!arr) continue;
          for (let i = 0; i < arr.length; i++) {
            const a = arr[i];
            const ax = a.x - x, ay = a.y - y, az = a.z - z;
            if (ax * ax + ay * ay + az * az < r2) out.push(a);
          }
        }
      }
    }
    return out;
  }

  /**
   * Flood-fill the empty space around a set of seeds.
   *
   * The density gates that select pockets also hollow out large cavities:
   * a voxel in the middle of a big void has too few atoms near it to satisfy
   * the wall count, so the surviving voxels are only the shell near the
   * walls, and the resulting "volume" is far too small. Filling the void
   * region outwards from the pocket's own voxels measures the cavity the way
   * a pocket finder is supposed to.
   *
   * `voidFlags` is the lattice of void voxels produced by the scan: a 1 byte
   * per lattice cell array, indexed as (i * ny + j) * nz + k.
   *
   * @param {Uint8Array} voidFlags
   * @param {number} ny number of lattice cells along y
   * @param {number} nz number of lattice cells along z
   * @param {Array<[number,number,number]>} seeds  lattice indices
   * @param {{min:[number,number,number], max:[number,number,number]}} limits
   * @param {number} maxVoxels safety cap
   * @returns {Array<[number,number,number]>} the filled lattice indices
   */
  function floodFillVoid(voidFlags, ny, nz, seeds, limits, maxVoxels) {
    const out = [];
    if (!voidFlags || !seeds || !seeds.length) return out;
    if (!ny || !nz) return out;
    const cap = maxVoxels || 20000;
    const idx = (i, j, k) => (i * ny + j) * nz + k;
    // The lattice's extent in x is implied by its length. Without this bound
    // the flat index is not injective outside the lattice - a cell at
    // k = -1 still lands on an in-range byte - so a fill with no limits box
    // walks out of the lattice and keeps going until the safety cap.
    const nx = Math.floor(voidFlags.length / (ny * nz));
    const seen = new Set();
    const stack = [];
    for (const s of seeds) {
      const kk = `${s[0]}|${s[1]}|${s[2]}`;
      if (seen.has(kk)) continue;
      seen.add(kk);
      stack.push([s[0], s[1], s[2]]);
    }
    const lo = limits ? limits.min : [-Infinity, -Infinity, -Infinity];
    const hi = limits ? limits.max : [Infinity, Infinity, Infinity];
    while (stack.length && out.length < cap) {
      const c = stack.pop();
      const i = c[0], j = c[1], k = c[2];
      if (i < 0 || i >= nx || j < 0 || j >= ny || k < 0 || k >= nz) continue;
      if (i < lo[0] || i > hi[0] || j < lo[1] || j > hi[1] || k < lo[2] || k > hi[2]) continue;
      if (!voidFlags[idx(i, j, k)]) continue;
      out.push([i, j, k]);
      const nbrs = [[i + 1, j, k], [i - 1, j, k], [i, j + 1, k], [i, j - 1, k],
                    [i, j, k + 1], [i, j, k - 1]];
      for (const n of nbrs) {
        const kk = `${n[0]}|${n[1]}|${n[2]}`;
        if (seen.has(kk)) continue;
        seen.add(kk);
        stack.push(n);
      }
    }
    return out;
  }

  const api = {
    buildAtomHash,
    countAtomsNear,
    atomsNear,
    floodFillVoid,
    rayHits,
    enclosureFraction,
  };

  global.PFE = global.PFE || {};
  global.PFE.pockets = Object.freeze(api);
})(typeof globalThis !== "undefined" ? globalThis : this);
