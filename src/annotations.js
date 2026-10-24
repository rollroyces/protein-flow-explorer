/* ===================================================================
 * Protein Flow Explorer — annotation key module
 *
 * The single source of truth for the app's annotation key: the list of
 * colour schemes the Presentation tab documents, the colours each scheme
 * actually uses, and the rule that decides whether a scheme is available
 * for the structure on screen.
 *
 * This module exists because of a specific failure. The key used to live
 * as a local constant in app.js, far away from the colour tables in
 * src/colors.js and src/variants.js. Two months ago every scheme except
 * "By chain" rendered the whole structure black — the app had passed a
 * scheme *name* where 3Dmol wanted a literal colour — and the key went
 * on advertising colours the app never produced. A clean lint run and a
 * passing test suite both missed it. The failure mode is drift: the
 * documentation of the colours and the colours themselves are maintained
 * separately, and nothing checks that they agree.
 *
 * So the key lives here, and test/annotations.test.js cross-checks every
 * gradient in it against the palettes the app really renders (PFE.colors,
 * PFE.variants). A future edit that changes one side without the other
 * fails the test suite instead of shipping.
 *
 * No DOM, no network, no app state — unit-testable in isolation. Loaded
 * as a classic <script> before app.js and published on the single PFE
 * namespace. The module is deliberately standalone: it does not import
 * the other src/*.js modules, so it loads on its own and so the
 * cross-module comparison happens in exactly one place (the test).
 * =================================================================== */
(function (global) {
  "use strict";

  // The five conditions that decide whether a colour scheme can be used
  // for the structure on screen. The app computes these from the loaded
  // record; a key entry names the flag it depends on in its `needs` field.
  //
  // The names are a constant so that a typo cannot silently disable a key
  // entry — "crytal" instead of "crystal" would otherwise look like a
  // working, simply-unavailable scheme.
  const AVAILABILITY_FLAGS = Object.freeze([
    "crystal",      // experimental method is X-RAY / NMR / EM / ELECTRON / CRYO
    "alphafold",    // the loaded model is an AlphaFold model
    "conservation", // a conservation / MSA result is loaded
    "domains",      // the loaded protein has annotated domains
    "variants",     // variants are loaded and at least one sits on the structure
  ]);

  // The key itself, one entry per colour scheme the app can render.
  //
  // `stops` is the truth: the hex colours the scheme uses, in order. The
  // CSS gradient string is derived from it below, so the swatch the panel
  // paints and the colours the structure wears cannot drift apart inside
  // this file.
  //
  // `needs` names the availability flag a structure must satisfy for the
  // entry to be usable. An entry without `needs` is always available.
  const KEY = [
    {
      id: "chain", name: "By chain", quantity: "polymer entity / chain",
      // All ten palette colours: the scheme renders all ten, so the key must
      // show all ten. A key that shows five of the colours the structure
      // can take is the same kind of lie as a key that shows none.
      stops: ["#2f7ed8", "#ef553b", "#00a878", "#f2b705", "#9b5de5",
              "#00b4d8", "#e76f51", "#6a994e", "#c77dff", "#ff8fab"],
      use: "Locating a domain or subunit across a complex.",
    },
    {
      id: "byDomain", name: "By domain", quantity: "annotated structural domain",
      stops: ["#2f7ed8", "#9aa0a6", "#cc1f1f"],
      use: "Separating the parts of a multi-domain protein that move independently.",
      needs: "domains",
    },
    {
      id: "ss", name: "Secondary structure", quantity: "helix / strand / coil (HELIX + SHEET records)",
      stops: ["#ef5350", "#f5c542", "#9aa0a6"],
      use: "Reading fold architecture; finding helices that pack into the core.",
    },
    {
      id: "spectrum", name: "Spectrum (N→C)", quantity: "residue index along the chain",
      stops: ["#0053d6", "#00b4d8", "#f5c542", "#cc1f1f"],
      use: "Following a sequence walk; spotting a domain boundary.",
    },
    {
      id: "bFactor", name: "B-factor", quantity: "atomic displacement (Å²)",
      stops: ["#2f7ed8", "#00a878", "#f5c542", "#cc1f1f"],
      use: "Experimental flexibility in the crystal. High = disordered or mobile.",
      needs: "crystal",
    },
    {
      id: "pLDDT", name: "pLDDT (AlphaFold)", quantity: "per-residue prediction confidence (0–100)",
      stops: ["#ff7d45", "#ffdb13", "#65cbf3", "#0053d6"],
      use: "Trust in an AlphaFold model. <50 is usually a disordered tail.",
      needs: "alphafold",
    },
    {
      id: "conservation", name: "Conservation", quantity: "MSA column identity across homologs",
      stops: ["#5b8def", "#f5c542", "#cc1f1f"],
      use: "Finding functional sites — the core is invariant, the surface is not.",
      needs: "conservation",
    },
    {
      id: "variant", name: "Variant significance", quantity: "ClinVar clinical significance",
      // Every ClinVar class the scheme renders, in severity order. Drug
      // response matters most here -- EGFR T790M and L858R are classified
      // that way, as resistance mutations rather than inherited disease
      // alleles -- so it must be findable in the key.
      stops: ["#b91c1c", "#f97316", "#a78bfa", "#ec4899", "#fbbf24",
              "#8b5cf6", "#84cc16", "#22c55e", "#2dd4bf", "#cccccc"],
      use: "Where the known disease variants sit, and how severe they are.",
      needs: "variants",
    },
    {
      id: "hydrophobicity", name: "Hydrophobicity", quantity: "Kyte–Doolittle side-chain score",
      stops: ["#0053d6", "#f5c542", "#cc1f1f"],
      use: "Reading a surface for a ligand-facing hydrophobic face.",
    },
  ];

  /**
   * Canonicalise one colour stop to "#rrggbb".
   *
   * Accepts a hex string in any case, or an RGB integer in the 0xrrggbb
   * range (the form src/colors.js stores its palettes in), so a caller can
   * build a gradient straight from a colour table.
   *
   * @param {string|number} stop  e.g. "#2f7ed8" or 0x2f7ed8
   * @returns {string|null} the lower-case "#rrggbb" form, or null if invalid
   */
  function normaliseStop(stop) {
    if (typeof stop === "number") {
      if (!Number.isFinite(stop) || stop < 0 || stop > 0xffffff) return null;
      return "#" + Math.round(stop).toString(16).padStart(6, "0");
    }
    if (typeof stop === "string") {
      const s = stop.trim().toLowerCase();
      return /^#[0-9a-f]{6}$/.test(s) ? s : null;
    }
    return null;
  }

  /**
   * Build the CSS gradient string for a list of colour stops.
   *
   * @param {Array<string|number>} stops  e.g. ["#2f7ed8", "#ef553b"]
   * @returns {string} e.g. "linear-gradient(90deg,#2f7ed8,#ef553b)"; "" if
   *   fewer than two valid stops were given
   */
  function gradientCss(stops) {
    if (!Array.isArray(stops) || stops.length < 2) return "";
    const parts = [];
    for (const s of stops) {
      const hex = normaliseStop(s);
      if (!hex) return "";
      parts.push(hex);
    }
    return `linear-gradient(90deg,${parts.join(",")})`;
  }

  /**
   * Parse a CSS linear-gradient into its angle and its hex stops.
   *
   * A gradient that contains anything other than "#rrggbb" stops is
   * reported as invalid rather than partially parsed, so a half-understood
   * gradient can never pass validation.
   *
   * @param {string} grad  e.g. "linear-gradient(90deg,#ef5350,#f5c542)"
   * @returns {{angle: string, stops: string[], valid: boolean}}
   */
  function parseGradient(grad) {
    const out = { angle: "", stops: [], valid: false };
    if (typeof grad !== "string") return out;
    const m = grad.match(/^\s*linear-gradient\(([^()]*)\)\s*$/i);
    if (!m) return out;
    const parts = m[1].split(",").map((p) => p.trim()).filter((p) => p.length);
    if (!parts.length) return out;
    // The first token is normally the angle ("90deg"); a gradient written
    // without one is read as all-stops, defaulting to 90deg.
    const stopParts = normaliseStop(parts[0]) ? parts : parts.slice(1);
    out.angle = normaliseStop(parts[0]) ? "90deg" : parts[0];
    const stops = [];
    for (const p of stopParts) {
      const hex = normaliseStop(p);
      if (!hex) return out;              // a non-hex stop invalidates the lot
      stops.push(hex);
    }
    out.stops = stops;
    out.valid = stops.length >= 2;
    return out;
  }

  const ANNOTATIONS = Object.freeze(KEY.map((a) => {
    const entry = {
      id: a.id,
      name: a.name,
      quantity: a.quantity,
      stops: Object.freeze(a.stops.slice()),
      grad: gradientCss(a.stops),
      use: a.use,
    };
    if (a.needs) entry.needs = a.needs;
    return Object.freeze(entry);
  }));

  /**
   * Look up one key entry by id.
   * @param {string} id  e.g. "ss"
   * @returns {object|null} the entry, or null if there is no such scheme
   */
  function annotationFor(id) {
    if (typeof id !== "string") return null;
    for (const a of ANNOTATIONS) if (a.id === id) return a;
    return null;
  }

  /**
   * The hex stops of one key entry.
   * @param {string} id  e.g. "ss"
   * @returns {string[]} a copy of the stops, e.g. ["#ef5350","#f5c542","#9aa0a6"]
   */
  function stopsFor(id) {
    const a = annotationFor(id);
    return a ? a.stops.slice() : [];
  }

  /**
   * Decide which key entries are usable for the structure on screen.
   *
   * An entry with no `needs` is always available; an entry with `needs` is
   * available only when the flag it names is set.
   *
   * @param {object} flags  a plain object keyed by AVAILABILITY_FLAGS
   * @returns {object} a map of annotation id -> boolean
   */
  function availability(flags) {
    const f = flags && typeof flags === "object" ? flags : {};
    const out = {};
    for (const a of ANNOTATIONS) out[a.id] = a.needs ? !!f[a.needs] : true;
    return out;
  }

  /**
   * Whether one key entry is usable for the structure on screen.
   * @param {string} id  the annotation id
   * @param {object} flags  a plain object keyed by AVAILABILITY_FLAGS
   * @returns {boolean}
   */
  function isAvailable(id, flags) {
    const a = annotationFor(id);
    if (!a) return false;
    if (!a.needs) return true;
    return !!(flags && typeof flags === "object" && flags[a.needs]);
  }

  /**
   * Flag names the caller passed that are not part of AVAILABILITY_FLAGS.
   *
   * A typo on the calling side is the mirror image of a typo in `needs`:
   * it silently disables an entry rather than enabling one. Reporting the
   * unknown names lets the app say so instead of showing a greyed-out key.
   *
   * @param {object} flags  a plain object keyed by AVAILABILITY_FLAGS
   * @returns {string[]} the unrecognised names, empty when there are none
   */
  function unknownFlags(flags) {
    if (!flags || typeof flags !== "object") return [];
    return Object.keys(flags).filter((k) => !AVAILABILITY_FLAGS.includes(k));
  }

  /**
   * Check one key entry.
   *
   * @param {object} entry  a key entry
   * @returns {{ok: boolean, errors: string[]}} ok is false when the entry
   *   is malformed; errors says why
   */
  function validate(entry) {
    const errors = [];
    if (!entry || typeof entry !== "object") {
      return { ok: false, errors: ["entry is not an object"] };
    }
    for (const field of ["id", "name", "quantity", "use"]) {
      if (typeof entry[field] !== "string" || !entry[field].trim()) {
        errors.push(`missing or empty "${field}"`);
      }
    }
    const parsed = parseGradient(entry.grad);
    if (!parsed.valid) {
      errors.push('grad is not a linear-gradient with at least two "#rrggbb" stops');
    }
    if (entry.stops !== undefined) {
      if (!Array.isArray(entry.stops) || entry.stops.length < 2) {
        errors.push('stops must be an array of at least two hex colours');
      } else if (entry.stops.some((s) => !normaliseStop(s))) {
        errors.push('stops contains a value that is not a "#rrggbb" colour');
      } else if (parsed.valid) {
        const fromStops = entry.stops.map((s) => normaliseStop(s));
        if (fromStops.join(",") !== parsed.stops.join(",")) {
          errors.push("grad does not match stops");
        }
      }
    }
    if (entry.needs !== undefined && !AVAILABILITY_FLAGS.includes(entry.needs)) {
      errors.push(`needs "${entry.needs}" is not one of ${AVAILABILITY_FLAGS.join(", ")}`);
    }
    return { ok: errors.length === 0, errors };
  }

  /**
   * Check a whole key, including duplicate ids.
   *
   * @param {object[]} [entries]  the entries to check; defaults to ANNOTATIONS
   * @returns {{ok: boolean, errors: string[]}}
   */
  function validateAll(entries) {
    const list = Array.isArray(entries) ? entries : ANNOTATIONS;
    const errors = [];
    const seen = new Set();
    list.forEach((entry, i) => {
      const id = entry && typeof entry.id === "string" ? entry.id : `#${i}`;
      for (const e of validate(entry).errors) errors.push(`${id}: ${e}`);
      if (entry && typeof entry.id === "string") {
        if (seen.has(entry.id)) errors.push(`${entry.id}: duplicate id`);
        seen.add(entry.id);
      }
    });
    return { ok: errors.length === 0, errors };
  }

  /**
   * Throw when the key is malformed.
   *
   * The loud counterpart of validateAll: a broken edit to the key should
   * fail loudly at start-up rather than quietly produce a key with a
   * missing swatch or a silently disabled entry.
   *
   * @param {object[]} [entries]  the entries to check; defaults to ANNOTATIONS
   * @returns {object[]} the validated entries
   * @throws {Error} when any entry is malformed
   */
  function assertValid(entries) {
    const res = validateAll(entries);
    if (!res.ok) {
      throw new Error(`Annotation key is malformed:\n  - ${res.errors.join("\n  - ")}`);
    }
    return Array.isArray(entries) ? entries : ANNOTATIONS;
  }

  const api = {
    ANNOTATIONS,
    AVAILABILITY_FLAGS,
    gradientCss,
    parseGradient,
    annotationFor,
    stopsFor,
    availability,
    isAvailable,
    unknownFlags,
    validate,
    validateAll,
    assertValid,
  };

  global.PFE = global.PFE || {};
  global.PFE.annotations = Object.freeze(api);
})(typeof globalThis !== "undefined" ? globalThis : this);
