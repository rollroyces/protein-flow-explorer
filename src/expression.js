/* ===================================================================
 * Protein Flow Explorer — tissue expression module
 *
 * Pure logic for turning GTEx tissue-expression data into values a panel
 * can display:
 *
 *   - readable labels for GTEx's machine tissue ids
 *   - log-scale normalisation of median TPM
 *   - a colour ramp for an expression value
 *   - a ranking of tissues by expression
 *   - a specificity verdict (enriched / enhanced / low / not detected)
 *   - a grouping of the 54 tissues into readable compartments
 *
 * No DOM, no network, no app state — unit-testable in isolation. The
 * fetching lives in app.js; this module only ever sees rows that have
 * already arrived:
 *
 *   { tissueSiteDetailId: "Skin_Sun_Exposed_Lower_leg", median: 78.3374 }
 *
 * The colour ramp is interpolated locally instead of being borrowed from
 * src/colors.js: every src/ module has to stand on its own, because
 * index.html loads them as plain <script> tags in a fixed order and none
 * of them may depend on another's load order.
 *
 * Extracted during the 2026 biomarker audit. Loaded as a classic <script>
 * before app.js and published on the single PFE namespace.
 * =================================================================== */
(function (global) {
  "use strict";

  // ---- thresholds -------------------------------------------------------
  //
  // Detection. GTEx treats a gene as detected in a tissue when the median
  // TPM is at or above 1; below that the reads are background noise.
  const DETECTION_TPM = 1;

  // Elevation. The standard 4x convention used by GTEx and the Human
  // Protein Atlas: a tissue counts as elevated when its median TPM is at
  // least four times the expression seen elsewhere.
  //
  // The comparison is made against BOTH the median and the mean of the
  // OTHER tissues, and the stricter of the two wins. The median on its own
  // is not enough: a housekeeping gene such as ACTB has a median of
  // ~2200 TPM and a spread of 5000-11600, so four times its median is
  // still reached by tissues that are merely high, and the gene would be
  // called "tissue enhanced" when it is in fact expressed everywhere.
  // Measured against the mean of the other tissues the same gene is flat,
  // which is the correct verdict.
  //
  // The tissues being judged are left out of their own baseline: a tissue
  // compared against a median that includes itself is measured against a
  // number it helped create.
  const ENHANCED_FOLD = 4;

  // Dominance. For "tissue enriched" the leading tissue must also stand
  // clear of the runner-up, so a gene expressed in one tissue (or in one
  // group of tissues) is not confused with a gene that is merely higher in
  // a couple of places.
  const ENRICHED_RATIO = 10;

  // ---- colour ramp ------------------------------------------------------
  //
  // Expression, low to high. Muted slate for "nothing here", then the
  // app's own teal-green, amber and red — the same three colours the
  // b-factor and conservation ramps already use, so the panel reads as one
  // system. Every stop is legible against the dark UI (--bg-1 #161b22).
  const EXPRESSION_STOPS = [
    0x64748b, // muted slate blue — undetectable / background
    0x00a878, // teal-green — detected, low
    0xf5c542, // amber — moderate
    0xcc1f1f, // red — very high
  ];

  // Top of the ramp, in TPM. log10(1001) = 3, so the log scale runs
  // 0 -> 1000 TPM over log values 0 -> 3. Values above the ceiling
  // saturate at the red stop; the log keeps 0.045 TPM and 78 TPM apart,
  // which a linear ramp could not do.
  const EXPRESSION_TPM_MAX = 1000;

  // ---- tissue labels ----------------------------------------------------
  //
  // GTEx's tissue ids are machine names ("Skin_Sun_Exposed_Lower_leg").
  // These are the readable names for every tissue in the v8 median
  // expression dataset (54 of them); anything else falls back to a
  // mechanical underscore-to-space title case.
  const TISSUE_LABELS = {
    Adipose_Subcutaneous: "Adipose (subcutaneous)",
    Adipose_Visceral_Omentum: "Adipose (visceral, omentum)",
    Adrenal_Gland: "Adrenal gland",
    Artery_Aorta: "Artery (aorta)",
    Artery_Coronary: "Artery (coronary)",
    Artery_Tibial: "Artery (tibial)",
    Bladder: "Bladder",
    Brain_Amygdala: "Brain — amygdala",
    Brain_Anterior_cingulate_cortex_BA24: "Brain — anterior cingulate cortex (BA24)",
    Brain_Caudate_basal_ganglia: "Brain — caudate (basal ganglia)",
    Brain_Cerebellar_Hemisphere: "Brain — cerebellar hemisphere",
    Brain_Cerebellum: "Brain — cerebellum",
    Brain_Cortex: "Brain — cortex",
    Brain_Frontal_Cortex_BA9: "Brain — frontal cortex (BA9)",
    Brain_Hippocampus: "Brain — hippocampus",
    Brain_Hypothalamus: "Brain — hypothalamus",
    Brain_Nucleus_accumbens_basal_ganglia: "Brain — nucleus accumbens (basal ganglia)",
    Brain_Putamen_basal_ganglia: "Brain — putamen (basal ganglia)",
    "Brain_Spinal_cord_cervical_c-1": "Brain — spinal cord (cervical C1)",
    Brain_Substantia_nigra: "Brain — substantia nigra",
    Breast_Mammary_Tissue: "Breast (mammary tissue)",
    Cells_Cultured_fibroblasts: "Cultured fibroblasts",
    "Cells_EBV-transformed_lymphocytes": "EBV-transformed lymphocytes",
    Cervix_Ectocervix: "Cervix (ectocervix)",
    Cervix_Endocervix: "Cervix (endocervix)",
    Colon_Sigmoid: "Colon (sigmoid)",
    Colon_Transverse: "Colon (transverse)",
    Esophagus_Gastroesophageal_Junction: "Esophagus (gastroesophageal junction)",
    Esophagus_Mucosa: "Esophagus (mucosa)",
    Esophagus_Muscularis: "Esophagus (muscularis)",
    Fallopian_Tube: "Fallopian tube",
    Heart_Atrial_Appendage: "Heart (atrial appendage)",
    Heart_Left_Ventricle: "Heart (left ventricle)",
    Kidney_Cortex: "Kidney (cortex)",
    Kidney_Medulla: "Kidney (medulla)",
    Liver: "Liver",
    Lung: "Lung",
    Minor_Salivary_Gland: "Minor salivary gland",
    Muscle_Skeletal: "Muscle (skeletal)",
    Nerve_Tibial: "Nerve (tibial)",
    Ovary: "Ovary",
    Pancreas: "Pancreas",
    Pituitary: "Pituitary",
    Prostate: "Prostate",
    Skin_Not_Sun_Exposed_Suprapubic: "Skin (not sun-exposed suprapubic)",
    Skin_Sun_Exposed_Lower_leg: "Skin (sun-exposed lower leg)",
    Small_Intestine_Terminal_Ileum: "Small intestine (terminal ileum)",
    Spleen: "Spleen",
    Stomach: "Stomach",
    Testis: "Testis",
    Thyroid: "Thyroid",
    Uterus: "Uterus",
    Vagina: "Vagina",
    Whole_Blood: "Whole blood",
  };

  // ---- compartments -----------------------------------------------------
  //
  // A long list of 54 tissues is unreadable; these buckets let the UI
  // group them into something a person can scan.
  const COMPARTMENT_LABELS = {
    skin: "Skin",
    brain: "Brain",
    cardiovascular: "Cardiovascular",
    muscle: "Muscle",
    digestive: "Digestive",
    respiratory: "Respiratory",
    renal: "Renal & urinary",
    endocrine: "Endocrine",
    reproductive: "Reproductive",
    adipose: "Adipose",
    immune: "Immune & blood",
    nervous: "Nervous",
    cellular: "Cell lines",
    other: "Other",
  };

  // The order compartments are listed in.
  const COMPARTMENT_ORDER = [
    "skin", "brain", "cardiovascular", "muscle", "digestive", "respiratory",
    "renal", "endocrine", "reproductive", "adipose", "immune", "nervous",
    "cellular", "other",
  ];

  // Where every GTEx v8 tissue sits. An id that is not in this table falls
  // back to "other" (see OTHER_COMPARTMENT).
  const TISSUE_COMPARTMENTS = {
    Adipose_Subcutaneous: "adipose",
    Adipose_Visceral_Omentum: "adipose",
    Adrenal_Gland: "endocrine",
    Artery_Aorta: "cardiovascular",
    Artery_Coronary: "cardiovascular",
    Artery_Tibial: "cardiovascular",
    Bladder: "renal",
    Brain_Amygdala: "brain",
    Brain_Anterior_cingulate_cortex_BA24: "brain",
    Brain_Caudate_basal_ganglia: "brain",
    Brain_Cerebellar_Hemisphere: "brain",
    Brain_Cerebellum: "brain",
    Brain_Cortex: "brain",
    Brain_Frontal_Cortex_BA9: "brain",
    Brain_Hippocampus: "brain",
    Brain_Hypothalamus: "brain",
    Brain_Nucleus_accumbens_basal_ganglia: "brain",
    Brain_Putamen_basal_ganglia: "brain",
    "Brain_Spinal_cord_cervical_c-1": "brain",
    Brain_Substantia_nigra: "brain",
    Breast_Mammary_Tissue: "reproductive",
    Cells_Cultured_fibroblasts: "cellular",
    "Cells_EBV-transformed_lymphocytes": "immune",
    Cervix_Ectocervix: "reproductive",
    Cervix_Endocervix: "reproductive",
    Colon_Sigmoid: "digestive",
    Colon_Transverse: "digestive",
    Esophagus_Gastroesophageal_Junction: "digestive",
    Esophagus_Mucosa: "digestive",
    Esophagus_Muscularis: "digestive",
    Fallopian_Tube: "reproductive",
    Heart_Atrial_Appendage: "cardiovascular",
    Heart_Left_Ventricle: "cardiovascular",
    Kidney_Cortex: "renal",
    Kidney_Medulla: "renal",
    Liver: "digestive",
    Lung: "respiratory",
    Minor_Salivary_Gland: "digestive",
    Muscle_Skeletal: "muscle",
    Nerve_Tibial: "nervous",
    Ovary: "reproductive",
    Pancreas: "endocrine",
    Pituitary: "endocrine",
    Prostate: "reproductive",
    Skin_Not_Sun_Exposed_Suprapubic: "skin",
    Skin_Sun_Exposed_Lower_leg: "skin",
    Small_Intestine_Terminal_Ileum: "digestive",
    Spleen: "immune",
    Stomach: "digestive",
    Testis: "reproductive",
    Thyroid: "endocrine",
    Uterus: "reproductive",
    Vagina: "reproductive",
    Whole_Blood: "immune",
  };

  // Fallback compartment for a tissue id that is not in the table above.
  const OTHER_COMPARTMENT = "other";

  // ---- specificity classes ---------------------------------------------
  //
  // Ordered from the most specific verdict to the least. "insufficient
  // data" covers a payload with fewer than two tissues, where no
  // comparison between tissues is possible at all.
  const SPECIFICITY_CLASSES = [
    "tissue_enriched",
    "tissue_enhanced",
    "low_tissue_specificity",
    "not_detected",
    "insufficient_data",
  ];

  const SPECIFICITY_LABELS = {
    tissue_enriched: "Tissue enriched",
    tissue_enhanced: "Tissue enhanced",
    low_tissue_specificity: "Low tissue specificity",
    not_detected: "Not detected",
    insufficient_data: "Insufficient data",
  };

  /**
   * Interpolate across a list of RGB stops.
   *
   * A local copy of the helper in src/colors.js — the modules deliberately
   * do not depend on one another.
   *
   * @param {number[]} stops  RGB integers, e.g. 0x00a878
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

  /**
   * Coerce a median expression value to a usable TPM number.
   *
   * Missing, non-numeric, negative and non-finite values all count as
   * zero: GTEx reports an undetectable gene as a zero median, so a
   * missing number means the same thing as "not expressed here".
   *
   * @param {*} value  the `median` field of a GTEx expression row
   * @returns {number} a finite, non-negative TPM
   */
  function toTpm(value) {
    if (typeof value !== "number" || !Number.isFinite(value)) return 0;
    if (value < 0) return 0;
    return value;
  }

  /**
   * Median of a list of numbers.
   * @param {number[]} values
   * @returns {number} 0 for an empty list
   */
  function median(values) {
    if (!values.length) return 0;
    const sorted = values.slice().sort((a, b) => a - b);
    const mid = Math.floor(sorted.length / 2);
    if (sorted.length % 2 === 1) return sorted[mid];
    return (sorted[mid - 1] + sorted[mid]) / 2;
  }

  /**
   * Mean of a list of numbers.
   * @param {number[]} values
   * @returns {number} 0 for an empty list
   */
  function mean(values) {
    if (!values.length) return 0;
    let sum = 0;
    for (const v of values) sum += v;
    return sum / values.length;
  }

  /**
   * Human-readable label for a GTEx tissue id.
   *
   *   "Skin_Sun_Exposed_Lower_leg" -> "Skin (sun-exposed lower leg)"
   *
   * Unmapped ids fall back to a mechanical underscore-to-space title case,
   * so a tissue this build has never heard of still prints something
   * readable.
   *
   * @param {string} tissueId  a GTEx tissueSiteDetailId
   * @returns {string}
   */
  function tissueLabel(tissueId) {
    if (typeof tissueId !== "string" || !tissueId) return "";
    const mapped = TISSUE_LABELS[tissueId];
    if (mapped) return mapped;
    return tissueId
      .split("_")
      .filter((part) => part.length > 0)
      .map((part) => part.charAt(0).toUpperCase() + part.slice(1))
      .join(" ");
  }

  /**
   * The compartment a GTEx tissue id belongs to.
   * @param {string} tissueId  a GTEx tissueSiteDetailId
   * @returns {string} a compartment id; "other" when the tissue is unmapped
   */
  function tissueCompartment(tissueId) {
    if (typeof tissueId !== "string" || !tissueId) return OTHER_COMPARTMENT;
    const mapped = TISSUE_COMPARTMENTS[tissueId];
    return mapped || OTHER_COMPARTMENT;
  }

  /**
   * Readable label for a compartment id.
   * @param {string} compartmentId
   * @returns {string}
   */
  function compartmentLabel(compartmentId) {
    if (typeof compartmentId !== "string" || !compartmentId) return COMPARTMENT_LABELS.other;
    return COMPARTMENT_LABELS[compartmentId] || compartmentId;
  }

  /**
   * Log-scale normalisation of a median TPM value.
   *
   * log10(TPM + 1) keeps zero at zero and separates the bottom of the
   * range (0.045 TPM in whole blood) from the top (78 TPM in sun-exposed
   * skin) in a way a linear scale cannot.
   *
   * @param {number} tpm  median transcripts per million
   * @returns {number} the log value; 0 for zero, negative or non-finite input
   */
  function logExpression(tpm) {
    if (typeof tpm !== "number" || !Number.isFinite(tpm) || tpm <= 0) return 0;
    return Math.log10(tpm + 1);
  }

  /**
   * Colour for an expression value, on the log scale.
   *
   * By default the ramp is absolute: it always spans 0 to
   * EXPRESSION_TPM_MAX TPM, so the colour of a given TPM means the same
   * thing for every gene. Passing maxTpm rescales the ramp to that
   * ceiling instead, which is useful when one gene's values sit far above
   * the scale (albumin in liver is ~25,000 TPM).
   *
   * @param {number} tpm      median transcripts per million
   * @param {number} [maxTpm] optional ceiling of the ramp, in TPM
   * @returns {number} RGB integer, e.g. 0x00a878
   */
  function expressionColor(tpm, maxTpm) {
    const log = logExpression(tpm);
    let ceiling = EXPRESSION_TPM_MAX;
    if (typeof maxTpm === "number" && Number.isFinite(maxTpm) && maxTpm > 0) {
      ceiling = maxTpm;
    }
    const logCeiling = Math.log10(ceiling + 1);
    if (!(logCeiling > 0)) return EXPRESSION_STOPS[0];
    const t = Math.max(0, Math.min(1, log / logCeiling));
    return rampColor(EXPRESSION_STOPS, t);
  }

  /**
   * Pull the expression rows out of a GTEx median-expression payload.
   *
   * Accepts the whole response ({ data: [...] }) or a bare array of rows.
   * Rows without a tissue id are dropped; a missing or broken `median`
   * becomes 0.
   *
   * @param {object|Array<object>} payload  a GTEx medianGeneExpression response
   * @returns {Array<object>} the usable rows
   */
  function parseMedianExpression(payload) {
    if (!payload) return [];
    const rows = Array.isArray(payload) ? payload
      : (Array.isArray(payload.data) ? payload.data : []);
    const out = [];
    for (const row of rows) {
      if (!row || typeof row !== "object") continue;
      const tissueId = typeof row.tissueSiteDetailId === "string" ? row.tissueSiteDetailId : "";
      if (!tissueId) continue;
      out.push({ ...row, tissueSiteDetailId: tissueId, median: toTpm(row.median) });
    }
    return out;
  }

  /**
   * Rank tissues by median TPM, highest first.
   *
   * Every entry is display-ready: a readable label, the compartment, the
   * TPM value, its log value and the colour for that value.
   *
   * @param {Array<object>} rows  GTEx median-expression rows
   * @returns {Array<{tissueId: string, label: string, compartment: string,
   *   tpm: number, log: number, color: number, detected: boolean}>}
   */
  function rankTissues(rows) {
    const list = Array.isArray(rows) ? rows : [];
    const entries = [];
    for (const row of list) {
      if (!row || typeof row !== "object") continue;
      const tissueId = typeof row.tissueSiteDetailId === "string" ? row.tissueSiteDetailId : "";
      if (!tissueId) continue;
      const tpm = toTpm(row.median);
      entries.push({
        tissueId,
        label: tissueLabel(tissueId),
        compartment: tissueCompartment(tissueId),
        tpm,
        log: logExpression(tpm),
        color: expressionColor(tpm),
        detected: tpm >= DETECTION_TPM,
      });
    }
    entries.sort((a, b) => (b.tpm - a.tpm) ||
      (a.tissueId < b.tissueId ? -1 : a.tissueId > b.tissueId ? 1 : 0));
    return entries;
  }

  /**
   * Group tissues into compartments, for display.
   *
   * The compartments come out in COMPARTMENT_ORDER, empty compartments are
   * dropped, and the tissues inside each compartment keep the ranking
   * order (highest expression first).
   *
   * @param {Array<object>} rows  GTEx median-expression rows
   * @returns {Array<{compartment: string, label: string, tissues: Array<object>}>}
   */
  function groupByCompartment(rows) {
    const entries = rankTissues(rows);
    const buckets = new Map();
    for (const id of COMPARTMENT_ORDER) buckets.set(id, []);
    for (const entry of entries) {
      const key = entry.compartment || OTHER_COMPARTMENT;
      if (!buckets.has(key)) buckets.set(key, []);
      buckets.get(key).push(entry);
    }
    const ordered = COMPARTMENT_ORDER.concat(
      [...buckets.keys()].filter((k) => !COMPARTMENT_ORDER.includes(k)));
    const out = [];
    for (const key of ordered) {
      const tissues = buckets.get(key);
      if (!tissues || tissues.length === 0) continue;
      out.push({ compartment: key, label: compartmentLabel(key), tissues });
    }
    return out;
  }

  /**
   * Summarise the expression of one gene across tissues.
   *
   * The verdict follows the thresholds documented at the top of this file:
   *
   *   not_detected          no tissue reaches DETECTION_TPM (1 TPM)
   *   tissue_enriched       one elevated tissue, ENRICHED_RATIO (10x) clear
   *                         of the runner-up
   *   tissue_enhanced       one or more tissues at least ENHANCED_FOLD (4x)
   *                         above the rest
   *   low_tissue_specificity  detected, but no tissue is several-fold above
   *                         the others
   *   insufficient_data     fewer than two tissues, so no comparison is
   *                         possible
   *
   * @param {Array<object>} rows  GTEx median-expression rows
   * @returns {object} the summary and the specificity verdict
   */
  function summarizeExpression(rows) {
    const entries = rankTissues(rows);
    const tpms = entries.map((e) => e.tpm);
    const detected = entries.filter((e) => e.detected);
    const medianTpm = median(tpms);
    const meanTpm = mean(tpms);

    // A tissue is elevated when its median TPM is at least ENHANCED_FOLD
    // times both the median and the mean of the OTHER tissues; the
    // stricter of the two comparisons decides.
    const elevated = [];
    for (const entry of entries) {
      if (entry.tpm < DETECTION_TPM) continue;
      const others = entries.length > 1
        ? entries.filter((e) => e !== entry).map((e) => e.tpm)
        : [];
      const baseline = Math.max(median(others), mean(others));
      if (baseline > 0 && entry.tpm >= ENHANCED_FOLD * baseline) elevated.push(entry);
      else if (baseline <= 0 && entry.tpm >= DETECTION_TPM) elevated.push(entry);
    }

    const top = entries.length ? entries[0] : null;
    const runnerUp = entries.length > 1 ? entries[1] : null;
    // The dominance of the leading tissue over the runner-up. A runner-up
    // of zero means the leading tissue is the only one expressed at all,
    // which is as dominant as it gets.
    const dominanceRatio = !top ? null
      : (!runnerUp || runnerUp.tpm <= 0 ? Infinity : top.tpm / runnerUp.tpm);

    let specificity;
    if (entries.length < 2) {
      specificity = "insufficient_data";
    } else if (detected.length === 0) {
      specificity = "not_detected";
    } else if (elevated.length === 1 &&
               dominanceRatio !== null && dominanceRatio >= ENRICHED_RATIO) {
      specificity = "tissue_enriched";
    } else if (elevated.length >= 1) {
      specificity = "tissue_enhanced";
    } else {
      specificity = "low_tissue_specificity";
    }

    return {
      specificity,
      specificityLabel: SPECIFICITY_LABELS[specificity],
      tissueCount: entries.length,
      detectedCount: detected.length,
      undetectedCount: entries.length - detected.length,
      medianTpm,
      meanTpm,
      maxTpm: top ? top.tpm : null,
      maxTissueId: top ? top.tissueId : null,
      maxTissueLabel: top ? top.label : null,
      dominantTissueId: top ? top.tissueId : null,
      dominantTissueLabel: top ? top.label : null,
      dominanceRatio,
      foldVsMedian: top && medianTpm > 0 ? top.tpm / medianTpm : null,
      elevatedCount: elevated.length,
      elevatedTissueIds: elevated.map((e) => e.tissueId),
    };
  }

  const api = {
    // Thresholds
    DETECTION_TPM,
    ENHANCED_FOLD,
    ENRICHED_RATIO,
    EXPRESSION_TPM_MAX,
    // Colour
    EXPRESSION_STOPS,
    // Tables
    TISSUE_LABELS,
    TISSUE_COMPARTMENTS,
    COMPARTMENT_LABELS,
    COMPARTMENT_ORDER,
    SPECIFICITY_CLASSES,
    SPECIFICITY_LABELS,
    // Functions
    tissueLabel,
    tissueCompartment,
    compartmentLabel,
    logExpression,
    expressionColor,
    parseMedianExpression,
    rankTissues,
    groupByCompartment,
    summarizeExpression,
  };

  global.PFE = global.PFE || {};
  global.PFE.expression = Object.freeze(api);
})(typeof globalThis !== "undefined" ? globalThis : this);
