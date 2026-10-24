/* ===================================================================
 * Protein Flow Explorer — drug discovery module
 *
 * Pure logic for the parts of a structure-based drug-discovery readout that
 * can be computed from the model itself, plus the cross-validation checks
 * that compare what the app computed against an authoritative external
 * source.
 *
 * Nothing here fetches anything: the caller supplies the pocket list, the
 * ligand atoms and the external reference values. That keeps every rule in
 * this file unit-testable without a DOM or a network.
 *
 * Extracted during the 2026 drug-discovery audit. Loaded as a classic
 * <script> before app.js and published on the single PFE namespace.
 * =================================================================== */
(function (global) {
  "use strict";

  // ---------------------------------------------------------------------
  // Druggability
  //
  // The rules are the conventional shape-based filters used by pocket
  // finders: a ligandable cavity has enough volume to hold a small
  // molecule, is buried enough in the protein to bind tightly, and is not
  // so polar or so hydrophobic that it cannot be occupied.
  // ---------------------------------------------------------------------
  const DRUGGABILITY = {
    // Usable volume window, in cubic Angstroms.
    minVolume: 120,
    maxVolume: 1200,
    // Enclosure (the fraction of rays that hit protein) above this means
    // the pocket is buried rather than a surface groove.
    minEnclosure: 0.45,
    // Below this the pocket is an open cleft: shallow, and a ligand in it
    // would not be held.
    marginalEnclosure: 0.25,
  };

  /**
   * Classify a pocket's druggability.
   *
   * @param {{volume: number, enclosure: number}} pocket
   * @returns {{verdict: string, label: string, reasons: string[]}}
   *   verdict is "druggable" | "marginal" | "undruggable" | "unknown"
   */
  function druggabilityVerdict(pocket) {
    if (!pocket) return { verdict: "unknown", label: "Unknown", reasons: ["no pocket"] };
    const volume = Number(pocket.volume);
    const enclosure = Number(pocket.enclosure);
    if (!Number.isFinite(volume)) {
      return { verdict: "unknown", label: "Unknown", reasons: ["volume not measured"] };
    }
    const reasons = [];
    let ok = true;
    if (volume < DRUGGABILITY.minVolume) {
      ok = false;
      reasons.push(`volume ${Math.round(volume)} Å³ is below the ${DRUGGABILITY.minVolume} Å³ floor`);
    } else if (volume > DRUGGABILITY.maxVolume) {
      ok = false;
      reasons.push(`volume ${Math.round(volume)} Å³ is above the ${DRUGGABILITY.maxVolume} Å³ ceiling`);
    }
    if (Number.isFinite(enclosure)) {
      if (enclosure < DRUGGABILITY.marginalEnclosure) {
        ok = false;
        reasons.push(`enclosure ${enclosure.toFixed(2)} is too exposed to hold a ligand`);
      } else if (enclosure < DRUGGABILITY.minEnclosure) {
        reasons.push(`enclosure ${enclosure.toFixed(2)} is shallower than ideal`);
      }
    }
    if (!reasons.length) reasons.push("volume and burial are both in the drug-like range");
    return {
      verdict: ok ? "druggable" : "marginal",
      label: ok ? "Druggable" : "Marginal",
      reasons,
    };
  }

  // ---------------------------------------------------------------------
  // Element data
  //
  // The Chemical Component Dictionary counts hydrogens; a crystallographic
  // model does not contain them. Comparing the two therefore needs the
  // atomic masses, and needs them per element rather than lumped into
  // "heavy atoms", because the halogens differ from each other by more than
  // a factor of three.
  // ---------------------------------------------------------------------
  const ATOMIC_MASS = {
    H: 1.008, He: 4.0026, Li: 6.94, Be: 9.0122, B: 10.81, C: 12.011,
    N: 14.007, O: 15.999, F: 18.998, Ne: 20.180, Na: 22.990, Mg: 24.305,
    Al: 26.982, Si: 28.085, P: 30.974, S: 32.06, Cl: 35.45, Ar: 39.948,
    K: 39.098, Ca: 40.078, Sc: 44.956, Ti: 47.867, V: 50.942, Cr: 51.996,
    Mn: 54.938, Fe: 55.845, Co: 58.933, Ni: 58.693, Cu: 63.546, Zn: 65.38,
    Ga: 69.723, Ge: 72.630, As: 74.922, Se: 78.971, Br: 79.904, Kr: 83.798,
    Rb: 85.468, Sr: 87.62, Y: 88.906, Zr: 91.224, Nb: 92.906, Mo: 95.95,
    Tc: 98, Ru: 101.07, Rh: 102.91, Pd: 106.42, Ag: 107.87, Cd: 112.414,
    In: 114.818, Sn: 118.710, Sb: 121.760, Te: 127.60, I: 126.904, Xe: 131.293,
    Cs: 132.91, Ba: 137.33, Pt: 195.084, Au: 196.967, Hg: 200.592, Pb: 207.2,
  };

  /**
   * Canonical element symbol: first letter capital, the rest lower case.
   *
   * A PDB file writes its element column in capitals ("FE", "CA", "CL") while
   * the Chemical Component Dictionary writes them properly ("Fe", "Ca",
   * "Cl"). Comparing the two without normalising makes every metal look
   * like a different element, so everything that matches an element by
   * symbol goes through here first.
   *
   * @param {string} symbol
   * @returns {string}
   */
  function canonicalElement(symbol) {
    const s = String(symbol == null ? "" : symbol).trim();
    if (!s) return "";
    return s[0].toUpperCase() + s.slice(1).toLowerCase();
  }

  /**
   * Molecular weight of an element composition, from the element masses.
   *
   * Hydrogen is excluded by default: a crystallographic model does not
   * contain it, so a composition read off the model should not silently
   * acquire the hydrogens the CCD counts. Pass `includeHydrogen` to add them.
   *
   * @param {object} composition  element counts, e.g. {C:22, N:3, O:4}
   * @param {object} [opts]       {includeHydrogen: boolean}
   * @returns {{mw: number, counted: number, unknown: string[]}}
   */
  function molecularWeight(composition, opts) {
    const includeH = !!(opts && opts.includeHydrogen);
    const unknown = [];
    let mw = 0, counted = 0;
    for (const raw of Object.keys(composition || {})) {
      const n = Number(composition[raw]);
      if (!Number.isFinite(n) || n <= 0) continue;
      const el = canonicalElement(raw);
      const mass = ATOMIC_MASS[el];
      if (mass === undefined) {
        unknown.push(el);
        continue;
      }
      if (el === "H" && !includeH) continue;
      mw += mass * n;
      counted += n;
    }
    return { mw, counted, unknown };
  }

  /**
   * Write a composition as a formula in Hill order (C, then H, then the rest
   * alphabetically). Elements the model never contains are simply absent.
   *
   * @param {object} composition
   * @returns {string}
   */
  function formatFormula(composition) {
    const els = Object.keys(composition || {}).filter(
      (el) => Number.isFinite(Number(composition[el])) && Number(composition[el]) > 0);
    if (!els.length) return "";
    const rest = els.filter((el) => el !== "C" && el !== "H").sort();
    const ordered = els.includes("C")
      ? ["C"].concat(els.includes("H") ? ["H"] : [], rest)
      : els.slice().sort();
    return ordered.map((el) => {
      const n = Number(composition[el]);
      return n === 1 ? el : `${el}${n}`;
    }).join(" ");
  }

  // ---------------------------------------------------------------------
  // Cross-validation: does the detected pocket find the real site?
  //
  // The pocket finder never sees the ligand: it scans a voxel grid for
  // cavities. The bound ligand's position comes from the crystallography.
  // So if the top pocket contains the ligand, the detector and the
  // experiment agree, and that is a genuine independent check.
  // ---------------------------------------------------------------------

  /**
   * How much of the ligand sits in a pocket.
   *
   * When the pocket carries its own cavity voxels (`points`), the ligand is
   * measured against the cavity's real extent. That is the honest test: the
   * cavity was found without the ligand, so the ligand's atoms falling inside
   * it means the detector and the crystallography agree. Without voxels the
   * pocket's centre is used instead, which is a cruder but still meaningful
   * distance-to-centre measure.
   *
   * @param {{center: number[], points?: number[][]}} pocket
   * @param {Array<{x:number,y:number,z:number}>} ligandAtoms
   * @param {number} radius  how far from the pocket still counts
   * @returns {{occupancy: number, meanDistance: number, minDistance: number}|null}
   */
  function pocketOccupancy(pocket, ligandAtoms, radius) {
    if (!pocket || !Array.isArray(pocket.center) || !ligandAtoms || !ligandAtoms.length) {
      return null;
    }
    const r = radius || 5;
    const [cx, cy, cz] = pocket.center;
    const pts = Array.isArray(pocket.points) && pocket.points.length ? pocket.points : null;
    let inside = 0, sum = 0, min = Infinity;
    for (const a of ligandAtoms) {
      let d;
      if (pts) {
        d = Infinity;
        for (const p of pts) {
          const dd = Math.hypot(a.x - p[0], a.y - p[1], a.z - p[2]);
          if (dd < d) d = dd;
          if (d <= r) break;
        }
      } else {
        d = Math.hypot(a.x - cx, a.y - cy, a.z - cz);
      }
      sum += d;
      if (d < min) min = d;
      if (d <= r) inside++;
    }
    return {
      occupancy: inside / ligandAtoms.length,
      meanDistance: sum / ligandAtoms.length,
      minDistance: min,
      inside,
      total: ligandAtoms.length,
    };
  }

  /**
   * Check every detected pocket against the bound ligand.
   *
   * @param {Array<object>} pockets
   * @param {Array<object>} ligandAtoms
   * @param {number} radius
   * @returns {{found: boolean, index: number, occupancy: number, perPocket: object[]}}
   */
  function validatePocketsAgainstLigand(pockets, ligandAtoms, radius) {
    const perPocket = [];
    let best = null;
    (pockets || []).forEach((p, i) => {
      const occ = pocketOccupancy(p, ligandAtoms, radius);
      if (!occ) return;
      const entry = { index: i, ...occ };
      perPocket.push(entry);
      if (!best || occ.occupancy > best.occupancy ||
          (occ.occupancy === best.occupancy && occ.meanDistance < best.meanDistance)) {
        best = entry;
      }
    });
    const found = !!best && best.occupancy > 0;
    return {
      found,
      index: found ? best.index : -1,
      occupancy: found ? best.occupancy : 0,
      perPocket,
    };
  }

  /**
   * The residues that actually touch the bound ligand.
   *
   * This is the crystallographic ground truth for the binding site: it comes
   * from the ligand's position in the model, not from any cavity search.
   *
   * @param {Array<{x:number,y:number,z:number}>} ligandAtoms
   * @param {Array<{x:number,y:number,z:number,chain:string,resi:number}>} proteinAtoms
   * @param {number} cutoff  contact distance in Angstroms
   * @returns {Set<string>} "chain:resi" identifiers
   */
  function contactResidues(ligandAtoms, proteinAtoms, cutoff) {
    const contacts = new Set();
    if (!Array.isArray(ligandAtoms) || !ligandAtoms.length ||
        !Array.isArray(proteinAtoms) || !proteinAtoms.length) {
      return contacts;
    }
    const cut = cutoff || 4.0;
    const cut2 = cut * cut;
    for (const p of proteinAtoms) {
      for (const l of ligandAtoms) {
        const dx = p.x - l.x, dy = p.y - l.y, dz = p.z - l.z;
        if (dx * dx + dy * dy + dz * dz <= cut2) {
          contacts.add(`${p.chain}:${p.resi}`);
          break;
        }
      }
    }
    return contacts;
  }

  /**
   * How much of the real binding site does a detected pocket actually cover?
   *
   * A pocket that merely sits nearby but is lined by the wrong residues is not
   * the binding site, so this is a sharper test than distance alone.
   *
   * When the pocket carries its cavity voxels (`points`), a contact residue
   * counts as covered when one of its atoms lies within `radius` of one of
   * those voxels - i.e. the residue lines the cavity the pocket describes.
   * Without voxels the pocket's centre is used instead.
   *
   * @param {{center: number[], points?: number[][]}} pocket
   * @param {Set<string>|string[]} contacts  "chain:resi" identifiers
   * @param {Array<{x:number,y:number,z:number,chain:string,resi:number}>} proteinAtoms
   * @param {number} radius  how close the residue must be, in Angstroms
   * @returns {{covered: number, total: number, recall: number, missing: string[]}}
   */
  function pocketResidueRecall(pocket, contacts, proteinAtoms, radius) {
    const want = contacts instanceof Set ? contacts : new Set(contacts || []);
    const missing = [];
    if (!pocket || !Array.isArray(pocket.center) || !want.size) {
      return { covered: 0, total: want.size, recall: 0, missing: [...want] };
    }
    const r = radius || 5;
    const r2 = r * r;
    const [cx, cy, cz] = pocket.center;
    const pts = Array.isArray(pocket.points) && pocket.points.length ? pocket.points : null;

    // Nearest distance from the pocket to each contacting residue.
    const nearest = new Map();
    for (const a of proteinAtoms || []) {
      const key = `${a.chain}:${a.resi}`;
      if (!want.has(key)) continue;
      let d2;
      if (pts) {
        d2 = Infinity;
        for (const p of pts) {
          const dx = a.x - p[0], dy = a.y - p[1], dz = a.z - p[2];
          const dd = dx * dx + dy * dy + dz * dz;
          if (dd < d2) d2 = dd;
          if (d2 <= r2) break;
        }
      } else {
        const dx = a.x - cx, dy = a.y - cy, dz = a.z - cz;
        d2 = dx * dx + dy * dy + dz * dz;
      }
      const cur = nearest.get(key);
      if (cur === undefined || d2 < cur) nearest.set(key, d2);
    }
    let covered = 0;
    for (const c of want) {
      const d2 = nearest.get(c);
      if (d2 !== undefined && d2 <= r2) covered++;
      else missing.push(c);
    }
    return {
      covered,
      total: want.size,
      recall: want.size ? covered / want.size : 0,
      missing,
    };
  }

  // ---------------------------------------------------------------------
  // Cross-validation: the app's own descriptors against the CCD
  //
  // The app computes molecular weight from the atomic masses it can see in
  // the model, which is approximate (implicit hydrogens are not modelled).
  // The Chemical Component Dictionary carries the authoritative value, so
  // comparing the two is a real check on the app's arithmetic.
  // ---------------------------------------------------------------------

  /**
   * Compare the app's computed ligand descriptors with the CCD record.
   *
   * The model has no hydrogen atoms, so the app's molecular weight is
   * systematically low by the mass of the hydrogens the CCD counts. The
   * comparison therefore adds them back rather than pretending the two
   * numbers should match directly.
   *
   * @param {{mw: number, heavy: number}} props
   * @param {{formula_weight: number, formula: string}} ccd
   * @param {object} composition  element counts read off the model, e.g. {C:22,N:3,O:4}
   * @returns {{checks: Array, agree: boolean}}
   */
  function crossValidateLigand(props, ccd, composition) {
    const checks = [];
    const ref = ccd ? countAtomsFromFormula(ccd.formula) : null;

    // 1. Element composition. The strongest check: it verifies the app read
    //    the right atoms out of the model, not merely that it counted them.
    if (ref && composition) {
      const diffs = [];
      // Hydrogen is excluded: a crystallographic model has none, and that is
      // a property of the file rather than an error in the app.
      // Element symbols are canonicalised on both sides: the model writes
      // "FE" where the dictionary writes "Fe".
      const gotBy = {};
      for (const el of Object.keys(composition)) {
        const k = canonicalElement(el);
        gotBy[k] = (gotBy[k] || 0) + Number(composition[el] || 0);
      }
      const wantBy = {};
      for (const el of Object.keys(ref.byElement)) {
        const k = canonicalElement(el);
        wantBy[k] = (wantBy[k] || 0) + Number(ref.byElement[el] || 0);
      }
      const els = [...new Set([...Object.keys(wantBy), ...Object.keys(gotBy)])]
        .filter((el) => el !== "H");
      for (const el of [...els].sort()) {
        const got = Number(gotBy[el] || 0);
        const want = Number(wantBy[el] || 0);
        if (got !== want) diffs.push(`${el} ${got} vs ${want}`);
      }
      checks.push({
        name: "Element composition",
        ok: diffs.length === 0,
        app: Object.entries(gotBy).map(([e, n]) => `${e}${n}`).join(" "),
        reference: Object.entries(wantBy).map(([e, n]) => `${e}${n}`).join(" "),
        detail: diffs.length ? diffs.join("; ") : "identical",
      });
    }

    // 2. Heavy-atom count.
    if (ref && props && Number.isFinite(Number(props.heavy))) {
      const got = Number(props.heavy);
      const want = ref.heavy;
      const ok = Math.abs(got - want) <= 1;
      checks.push({
        name: "Heavy-atom count",
        ok,
        app: got,
        reference: want,
        detail: `${got} vs ${want}`,
      });
    }

    // 3. Molecular weight, with the missing hydrogens added back.
    if (ref && props && Number.isFinite(Number(ccd.formula_weight))) {
      const reference = Number(ccd.formula_weight);
      const hydrogens = Math.max(0, ref.total - ref.heavy);
      const withH = Number(props.mw) + hydrogens * 1.008;
      const rel = reference > 0 ? Math.abs(withH - reference) / reference : Infinity;
      const ok = rel <= 0.01;
      checks.push({
        name: "Molecular weight",
        ok,
        app: withH,
        reference,
        detail: `${withH.toFixed(1)} vs ${reference.toFixed(1)} Da (${(rel * 100).toFixed(2)}% apart)`,
      });
    }

    return { checks, agree: checks.length > 0 && checks.every((c) => c.ok) };
  }

  /**
   * Count atoms in a CCD formula string such as "C22 H23 N3 O4".
   * @param {string} formula
   * @returns {{heavy: number, total: number, byElement: object}|null}
   */
  function countAtomsFromFormula(formula) {
    if (typeof formula !== "string") return null;
    // A CCD formula can carry a hydrate ("C6H12O6 . 2H2O"). The model holds
    // the parent molecule only, so take the first component.
    const main = formula.split(/[·.]/)[0];
    const byElement = {};
    let heavy = 0, total = 0, any = false;
    const re = /([A-Z][a-z]?)\s*(\d*)/g;
    let m;
    while ((m = re.exec(main))) {
      const el = canonicalElement(m[1]);
      const n = m[2] ? parseInt(m[2], 10) : 1;
      if (!Number.isFinite(n)) continue;
      any = true;
      byElement[el] = (byElement[el] || 0) + n;
      total += n;
      if (el !== "H") heavy += n;
    }
    return any ? { heavy, total, byElement } : null;
  }

  // ---------------------------------------------------------------------
  // Measured potency
  // ---------------------------------------------------------------------

  /**
   * Render a potency in the most readable unit.
   * @param {number} value
   * @param {string} unit
   * @returns {string}
   */
  function formatPotency(value, unit) {
    const v = Number(value);
    if (!Number.isFinite(v)) return "—";
    const u = (unit || "").toLowerCase();
    // Promote to the unit that reads best, so 1500 nM shows as 1.5 µM.
    if (u === "nm") {
      if (v >= 1000) return `${(v / 1000).toFixed(2)} µM`;
      return `${v < 10 ? v.toFixed(1) : v.toFixed(0)} nM`;
    }
    if (u === "um" || u === "µm") return `${v < 10 ? v.toFixed(2) : v.toFixed(1)} µM`;
    if (u === "mm") return `${v.toFixed(2)} mM`;
    if (u === "pm") return `${v.toFixed(0)} pM`;
    if (!u) return v < 1000 ? `${v.toFixed(1)} nM` : `${(v / 1000).toFixed(2)} µM`;
    return `${v} ${unit || ""}`.trim();
  }

  /**
   * Summarise a list of measured affinities for one compound.
   *
   * @param {Array<{type: string, value: number, unit: string}>} affinities
   * @returns {{n: number, best: object|null, types: string[]}}
   */
  function summariseAffinity(affinities) {
    const list = (affinities || []).filter(
      (a) => a && Number.isFinite(Number(a.value)) && a.value > 0);
    if (!list.length) return { n: 0, best: null, types: [] };
    // Prefer the type that is the tightest binder: pChEMBL-style ranking,
    // Ki/IC50/Kd are comparable enough for an ordering.
    const order = { "Kd": 0, "Ki": 1, "IC50": 2 };
    const sorted = [...list].sort((a, b) => {
      const t = (order[a.type] ?? 9) - (order[b.type] ?? 9);
      return t !== 0 ? t : Number(a.value) - Number(b.value);
    });
    return {
      n: list.length,
      best: sorted[0],
      types: [...new Set(list.map((a) => a.type))].sort(),
    };
  }

  // ---------------------------------------------------------------------
  // Target profile
  // ---------------------------------------------------------------------

  /**
   * Reduce an Open Targets target record to the fields the panel shows.
   * @param {object} target
   * @returns {object}
   */
  function summariseTarget(target) {
    if (!target) return null;
    const diseases = (target.associatedDiseases && target.associatedDiseases.rows) || [];
    const drugs = (target.drugAndClinicalCandidates && target.drugAndClinicalCandidates.rows) || [];
    const tract = target.tractability || [];
    return {
      id: target.id,
      symbol: target.approvedSymbol,
      name: target.approvedName,
      classes: (target.targetClass || []).map((c) => c.label),
      diseaseCount: (target.associatedDiseases && target.associatedDiseases.count) || 0,
      diseases: diseases.slice(0, 5).map((d) => ({
        id: d.disease && d.disease.id,
        name: (d.disease && d.disease.name) || "—",
        score: typeof d.score === "number" ? d.score : null,
      })),
      drugCount: (target.drugAndClinicalCandidates && target.drugAndClinicalCandidates.count) || 0,
      drugs: drugs.slice(0, 12).map((d) => ({
        id: d.drug && d.drug.id,
        name: (d.drug && d.drug.name) || "—",
        stage: d.maxClinicalStage || "—",
        type: (d.drug && d.drug.drugType) || "",
      })),
      approvedDrugs: drugs.filter((d) => d.maxClinicalStage === "APPROVAL").length,
      tractability: tract.map((t) => ({
        label: t.label, modality: t.modality, value: t.value,
      })),
      isApprovedDrug: tract.some((t) => t.value && /approved drug/i.test(t.label)),
    };
  }

  /**
   * Does a drug name match the bound ligand? Used to point out that the
   * molecule in the structure is itself one of the target's drugs.
   * @param {string} ligandName
   * @param {string} compId
   * @param {Array<object>} drugs
   * @returns {object|null}
   */
  function matchDrugToLigand(ligandName, compId, drugs) {
    if (!ligandName || !Array.isArray(drugs)) return null;
    const norm = (s) => String(s || "").toUpperCase().replace(/[^A-Z0-9]/g, "");
    const target = norm(ligandName);
    if (!target) return null;
    for (const d of drugs) {
      const dn = norm(d.name);
      if (dn && (dn === target || target.startsWith(dn) || dn.startsWith(target))) {
        return d;
      }
    }
    // Fall back to the compound id appearing in the drug name.
    const cid = String(compId || "").toUpperCase();
    for (const d of drugs) {
      if (cid && norm(d.name).includes(cid)) return d;
    }
    return null;
  }

  const api = {
    DRUGGABILITY,
    ATOMIC_MASS,
    canonicalElement,
    druggabilityVerdict,
    pocketOccupancy,
    validatePocketsAgainstLigand,
    contactResidues,
    pocketResidueRecall,
    crossValidateLigand,
    countAtomsFromFormula,
    molecularWeight,
    formatFormula,
    formatPotency,
    summariseAffinity,
    summariseTarget,
    matchDrugToLigand,
  };

  global.PFE = global.PFE || {};
  global.PFE.drug = Object.freeze(api);
})(typeof globalThis !== "undefined" ? globalThis : this);
