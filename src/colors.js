/* ===================================================================
 * Protein Flow Explorer — colour module
 *
 * Pure colour maths and the colour tables that the annotation key
 * documents. No DOM, no network, no app state — unit-testable in isolation.
 *
 * Every ramp below matches the gradient printed in the app's annotation
 * key (ANNOTATIONS in app.js). The app therefore renders these through a
 * `colorfunc` rather than delegating to a 3Dmol scheme, so that:
 *
 *   - the colours on screen are exactly the colours the key promises;
 *   - the same colour is used for cartoon, stick, sphere, line and surface,
 *     instead of only for the one representation 3Dmol happens to support;
 *   - the ramp direction is under our control (3Dmol's sinebow and rwb
 *     gradients run opposite to the direction the key documents).
 *
 * Extracted during the 2026 colouring audit. Loaded as a classic <script>
 * before app.js and published on the single PFE namespace.
 * =================================================================== */
(function (global) {
  "use strict";

  // Colourblind-friendly chain palette. Replaces 3Dmol's built-in 'chain'
  // scheme, whose B/J/O/P/R/T chains render grey or near-black.
  const STRUCT_CHAIN_PALETTE = [
    0x2f7ed8, // blue
    0xef553b, // red-orange
    0x00a878, // teal-green
    0xf2b705, // amber
    0x9b5de5, // violet
    0x00b4d8, // cyan
    0xe76f51, // coral
    0x6a994e, // olive
    0xc77dff, // lilac
    0xff8fab, // pink
  ];

  // Secondary structure. Matches the key: helix red, strand yellow, coil grey.
  // 3Dmol stores atom.ss as lowercase 'h' / 's' / 'c', and transiently as
  // 'arrow start' / 'arrow end' while it builds cartoon arrows.
  const SS_COLORS = {
    h: 0xef5350,     // helix
    s: 0xf5c542,     // strand
    c: 0x9aa0a6,     // coil
    helix: 0xef5350,
    strand: 0xf5c542,
    sheet: 0xf5c542,
    coil: 0x9aa0a6,
    "arrow start": 0xf5c542,
    "arrow end": 0xf5c542,
  };
  const SS_FALLBACK = 0x9aa0a6;

  // Spectrum (N→C): blue at the N-terminus, red at the C-terminus.
  const SPECTRUM_STOPS = [0x0053d6, 0x00b4d8, 0xf5c542, 0xcc1f1f];

  // B-factor: blue = ordered, red = displaced.
  const BFACTOR_STOPS = [0x2f7ed8, 0x00a878, 0xf5c542, 0xcc1f1f];

  // Hydrophobicity: blue = hydrophilic, red = hydrophobic.
  const HYDROPHOBICITY_STOPS = [0x0053d6, 0xf5c542, 0xcc1f1f];

  // Conservation: blue = variable, red = highly conserved.
  const CONSERVATION_STOPS = [0x5b8def, 0xf5c542, 0xcc1f1f];

  // Kyte–Doolittle hydropathy index, keyed by three-letter residue code
  // (3Dmol's atom.resn is three-letter).
  const KD_SCALE = {
    ILE: 4.5, VAL: 4.2, LEU: 3.8, PHE: 2.8, CYS: 2.5, MET: 1.9, ALA: 1.8,
    GLY: -0.4, THR: -0.7, SER: -0.8, TRP: -0.9, TYR: -1.3, PRO: -1.6,
    HIS: -3.2, GLU: -3.5, GLN: -3.5, ASP: -3.5, ASN: -3.5, LYS: -3.9, ARG: -4.5,
  };
  const KD_MIN = -4.5;
  const KD_MAX = 4.5;

  /**
   * Interpolate across a list of RGB stops.
   * @param {number[]} stops  RGB integers, e.g. 0x0053d6
   * @param {number} t        position in [0, 1]
   * @returns {number} RGB integer
   */
  function rampColor(stops, t) {
    if (!Array.isArray(stops) || stops.length === 0) return 0xcccccc;
    if (!Number.isFinite(t)) return stops[stops.length - 1];
    if (t <= 0) return stops[0];
    if (t >= 1) return stops[stops.length - 1];
    const seg = (stops.length - 1) * t;
    const i = Math.min(stops.length - 2, Math.floor(seg));
    const f = seg - i;
    const a = stops[i], b = stops[i + 1];
    const ar = (a >> 16) & 255, ag = (a >> 8) & 255, ab = a & 255;
    const br = (b >> 16) & 255, bg = (b >> 8) & 255, bb = b & 255;
    const r = Math.round(ar + (br - ar) * f);
    const g = Math.round(ag + (bg - ag) * f);
    const bl = Math.round(ab + (bb - ab) * f);
    return (r << 16) | (g << 8) | bl;
  }

  /** Map a value in [min, max] onto [0, 1]. */
  function normalize(value, min, max) {
    if (typeof value !== "number" || !Number.isFinite(value)) return 0;
    if (!Number.isFinite(min) || !Number.isFinite(max) || max <= min) return 0;
    return (value - min) / (max - min);
  }

  /**
   * Secondary-structure colour.
   * @param {string} ss  3Dmol's atom.ss ('h', 's', 'c', or undefined)
   * @returns {number} RGB integer
   */
  function ssColor(ss) {
    if (typeof ss !== "string") return SS_FALLBACK;
    const key = ss.toLowerCase();
    return SS_COLORS[key] !== undefined ? SS_COLORS[key] : SS_FALLBACK;
  }

  /**
   * Spectrum colour along the chain, N-terminus to C-terminus.
   * @param {number} resi
   * @param {number} min  lowest residue index in the structure
   * @param {number} max  highest residue index in the structure
   * @returns {number} RGB integer
   */
  function spectrumColor(resi, min, max) {
    return rampColor(SPECTRUM_STOPS, normalize(resi, min, max));
  }

  /**
   * B-factor colour: blue = well ordered, red = displaced.
   * @param {number} b
   * @param {number} min
   * @param {number} max
   * @returns {number} RGB integer
   */
  function bfactorColor(b, min, max) {
    return rampColor(BFACTOR_STOPS, normalize(b, min, max));
  }

  /**
   * Hydrophobicity colour from the Kyte–Doolittle index.
   * @param {string} resn  three-letter residue name
   * @returns {number} RGB integer
   */
  function hydrophobicityColor(resn) {
    if (typeof resn !== "string") return HYDROPHOBICITY_STOPS[1];
    const v = KD_SCALE[resn.trim().toUpperCase()];
    if (v === undefined) return HYDROPHOBICITY_STOPS[1];
    return rampColor(HYDROPHOBICITY_STOPS, (v - KD_MIN) / (KD_MAX - KD_MIN));
  }

  /**
   * pLDDT colour, following the AlphaFold convention.
   *   >= 90 very high (deep blue) · >= 70 confident (cyan)
   *   >= 50 low (yellow)        · < 50 very low (orange)
   * @param {number} p
   * @returns {number} RGB integer
   */
  function plddtColor(p) {
    if (typeof p !== "number" || !Number.isFinite(p)) return 0xffdb13;
    if (p >= 90) return 0x0053d6;
    if (p >= 70) return 0x65cbf3;
    if (p >= 50) return 0xffdb13;
    return 0xff7d45;
  }

  /**
   * Conservation colour: 0 = variable (blue), 1 = invariant (red).
   * @param {number} c
   * @returns {number} RGB integer
   */
  function conservationColor(c) {
    if (typeof c !== "number" || !Number.isFinite(c)) return 0xcccccc;
    return rampColor(CONSERVATION_STOPS, Math.max(0, Math.min(1, c)));
  }

  const api = {
    STRUCT_CHAIN_PALETTE,
    SS_COLORS,
    SPECTRUM_STOPS,
    BFACTOR_STOPS,
    HYDROPHOBICITY_STOPS,
    CONSERVATION_STOPS,
    KD_SCALE,
    rampColor,
    normalize,
    ssColor,
    spectrumColor,
    bfactorColor,
    hydrophobicityColor,
    plddtColor,
    conservationColor,
  };

  global.PFE = global.PFE || {};
  global.PFE.colors = Object.freeze(api);
})(typeof globalThis !== "undefined" ? globalThis : this);
