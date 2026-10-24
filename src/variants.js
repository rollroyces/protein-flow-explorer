/* ===================================================================
 * Protein Flow Explorer — variant module
 *
 * Pure logic for turning clinical variant records into things you can put
 * on a structure:
 *
 *   - parsing ClinVar's protein_change ("N1095S, N881S") into positions
 *   - parsing ClinVar's canonical SPDI into a gnomAD variant id
 *   - classifying a ClinVar significance string
 *   - mapping a UniProt position onto a PDB residue number
 *
 * No DOM, no network, no app state — unit-testable in isolation.
 *
 * Extracted during the 2026 biomarker audit. Loaded as a classic <script>
 * before app.js and published on the single PFE namespace.
 * =================================================================== */
(function (global) {
  "use strict";

  // ClinVar significance classes, most severe / most actionable first.
  //
  // ClinVar's germline classification has more classes than most people
  // expect. "drug response" is the one that matters most here: EGFR T790M
  // and L858R are classified as drug response, not pathogenic, because they
  // are resistance mutations in a treated cancer rather than inherited
  // disease alleles.
  const SIGNIFICANCE_ORDER = [
    "pathogenic",
    "likely_pathogenic",
    "drug_response",
    "risk_factor",
    "vus",
    "conflicting",
    "likely_benign",
    "benign",
    "protective",
    "unknown",
  ];

  // Colours used by the "Variant" scheme and the variant spheres. These are
  // the same colours the panel's own legend prints (styles.css .dot-path and
  // friends), so the key and the structure always agree.
  const SIGNIFICANCE_COLORS = {
    pathogenic: 0xb91c1c,
    likely_pathogenic: 0xf97316,
    drug_response: 0xa78bfa,
    risk_factor: 0xec4899,
    vus: 0xfbbf24,
    conflicting: 0x8b5cf6,
    likely_benign: 0x84cc16,
    benign: 0x22c55e,
    protective: 0x2dd4bf,
    unknown: 0xcccccc,
  };

  // ClinVar's somatic clinical-impact tiers (AMP/ASCO/CAP).
  const IMPACT_TIERS = {
    "tier i - strong": "Tier I",
    "tier i - moderate": "Tier I",
    "tier ii": "Tier II",
    "tier iii": "Tier III",
    "tier iv": "Tier IV",
    "tier i": "Tier I",
    "tier ii - strong": "Tier II",
    "tier ii - moderate": "Tier II",
    "tier iii - moderate": "Tier III",
  };

  // ClinVar's oncogenicity classification.
  const ONCOGENICITY = {
    oncogenic: "Oncogenic",
    "likely oncogenic": "Likely oncogenic",
    "oncogenic other": "Oncogenic (other)",
    "likely oncogenic other": "Likely oncogenic (other)",
  };

  // Population-frequency bands. Thresholds follow the usual convention that
  // a variant in more than ~1% of the population is common, and one seen in
  // fewer than ~0.1% is rare.
  const RARITY_BANDS = [
    { max: 0.0001, label: "rare" },
    { max: 0.01, label: "uncommon" },
    { max: 0.05, label: "common" },
    { max: Infinity, label: "very common" },
  ];

  // Three-letter to one-letter residue codes, for the protein changes ClinVar
  // writes in "p.Asn1148Ser" form.
  const AA3_TO_1 = {
    ALA: "A", ARG: "R", ASN: "N", ASP: "D", CYS: "C", GLN: "Q", GLU: "E",
    GLY: "G", HIS: "H", ILE: "I", LEU: "L", LYS: "K", MET: "M", PHE: "F",
    PRO: "P", SER: "S", THR: "T", TRP: "W", TYR: "Y", VAL: "V", SEC: "U",
    PYL: "O", MSE: "M", ASX: "B", GLX: "Z", XLE: "J", TER: "*",
  };

  /**
   * Parse a protein-change string.
   *
   *   "N1148S"          -> { wt: "N", pos: 1148, mut: "S" }
   *   "Asn1148Ser"      -> { wt: "N", pos: 1148, mut: "S" }
   *   "Trp26Ter"        -> { wt: "W", pos: 26,  mut: "*" }
   *   "Leu1034fs"       -> { wt: "L", pos: 1034, mut: "fs" }
   *   "N1095S, N881S"   -> two changes
   *
   * Both one-letter and three-letter forms are accepted, because ClinVar
   * uses one-letter in the `protein_change` field and three-letter in the
   * `variation_name` field.
   *
   * @param {string} text
   * @returns {Array<{wt: string, pos: number, mut: string}>}
   */
  function parseProteinChange(text) {
    if (typeof text !== "string") return [];
    const out = [];
    for (const part of text.split(",")) {
      const s = part.trim().replace(/^p\./i, "").replace(/[()]/g, "");
      const m = s.match(/^([A-Z][a-z]{2}|\*|[A-Z])(\d+)(.*)$/);
      if (!m) continue;
      const pos = parseInt(m[2], 10);
      if (!Number.isFinite(pos) || pos < 1) continue;
      let wt = m[1] === "*" ? "*" : (AA3_TO_1[m[1].toUpperCase()] || m[1]);
      if (wt.length !== 1) continue;            // unknown three-letter code
      const rest = (m[3] || "").trim();
      let mut = "";
      if (rest) {
        const vm = rest.match(/^([A-Z][a-z]{2}|\*|[A-Z])(.*)$/);
        if (vm) {
          const v = vm[1] === "*" ? "*" : (AA3_TO_1[vm[1].toUpperCase()] || vm[1]);
          mut = v.length === 1 ? v : vm[1].toLowerCase();   // fs, del, dup, ext
        } else {
          mut = rest.toLowerCase();
        }
      }
      out.push({ wt, pos, mut });
    }
    return out;
  }

  /**
   * Pull the protein change out of a ClinVar variation name.
   *
   *   "NM_005228.5(EGFR):c.3443A>G (p.Asn1148Ser)" -> one change
   *
   * This is the better source than the `protein_change` field: the field
   * lists the same variant under every transcript ClinVar knows about, each
   * with its own residue numbering, whereas `variation_name` always refers
   * to the canonical (MANE) transcript -- which is the numbering UniProt
   * uses.
   *
   * @param {string} variationName
   * @returns {Array<{wt: string, pos: number, mut: string}>}
   */
  function extractProteinChange(variationName) {
    if (typeof variationName !== "string") return [];
    const m = variationName.match(/\((p\.[^)]+)\)/i);
    if (!m) return [];
    return parseProteinChange(m[1]);
  }

  /**
   * Parse a ClinVar `canonical_spdi` string.
   *
   *   "NC_000007.14:55205426:A:G" -> { chrom: "7", pos: 55205427, ref: "A", alt: "G" }
   *
   * SPDI positions are 0-based; gnomAD variant ids are 1-based.
   *
   * @param {string} spdi
   * @returns {{chrom: string, pos: number, ref: string, alt: string}|null}
   */
  function parseSpdi(spdi) {
    if (typeof spdi !== "string") return null;
    const m = spdi.trim().match(/^([^:]+):(\d+):([^:]*):([^:]*)$/);
    if (!m) return null;
    const acc = m[1].split(".")[0];            // NC_000007.14 -> NC_000007
    const num = acc.match(/(\d+)$/);
    if (!num) return null;
    const chrom = String(parseInt(num[1], 10));
    const pos = parseInt(m[2], 10) + 1;         // SPDI is 0-based
    const ref = (m[3] || "").toUpperCase();
    const alt = (m[4] || "").toUpperCase();
    if (!Number.isFinite(pos) || pos < 1 || !ref || !alt) return null;
    return { chrom, pos, ref, alt };
  }

  /**
   * Build a gnomAD variant id ("7-55205427-A-G") from a ClinVar record.
   * Falls back to the GRCh38 `variation_loc` when the SPDI is missing.
   *
   * @param {object} vs  a ClinVar variation_set entry
   * @returns {string|null}
   */
  function gnomadVariantId(vs) {
    if (!vs) return null;
    const spdi = parseSpdi(vs.canonical_spdi);
    if (spdi) return `${spdi.chrom}-${spdi.pos}-${spdi.ref}-${spdi.alt}`;
    // `variation_loc` carries the chromosome and position but its ref and alt
    // fields are empty in ClinVar's current response, so the id cannot be
    // completed from it. A missing SPDI means no gnomAD lookup.
    return null;
  }

  /**
   * Classify a ClinVar clinical-significance description.
   *
   * @param {string} description
   * @returns {string} one of SIGNIFICANCE_ORDER
   */
  function classifySignificance(description) {
    const t = String(description || "").toLowerCase().trim();
    if (!t) return "unknown";
    // Conflicting must be checked first: "Conflicting interpretations of
    // pathogenicity" contains "pathogenic".
    if (t.includes("conflict")) return "conflicting";
    if (t.includes("drug response")) return "drug_response";
    if (t.includes("risk factor")) return "risk_factor";
    if (t.includes("protective")) return "protective";
    const benign = t.includes("benign");
    const pathogenic = t.includes("pathogenic");
    if (pathogenic && !benign) return t.includes("likely") ? "likely_pathogenic" : "pathogenic";
    if (benign && !pathogenic) return t.includes("likely") ? "likely_benign" : "benign";
    if (t.includes("uncertain")) return "vus";
    return "vus";
  }

  /**
   * ClinVar's somatic clinical-impact tier, e.g. "Tier I - Strong".
   * @param {string} description
   * @returns {string|null}
   */
  function classifyImpactTier(description) {
    const t = String(description || "").toLowerCase().trim();
    if (!t) return null;
    return IMPACT_TIERS[t] || null;
  }

  /**
   * ClinVar's oncogenicity classification, e.g. "Oncogenic".
   * @param {string} description
   * @returns {string|null}
   */
  function classifyOncogenicity(description) {
    const t = String(description || "").toLowerCase().trim();
    if (!t) return null;
    return ONCOGENICITY[t] || null;
  }

  /**
   * Rank a significance class, most severe first.
   * @param {string} sig
   * @returns {number}
   */
  function significanceRank(sig) {
    const i = SIGNIFICANCE_ORDER.indexOf(sig);
    return i === -1 ? SIGNIFICANCE_ORDER.length : i;
  }

  /**
   * Label a population allele frequency.
   * @param {number} af
   * @returns {string}
   */
  function rarityLabel(af) {
    if (typeof af !== "number" || !Number.isFinite(af)) return "unknown";
    for (const band of RARITY_BANDS) if (af <= band.max) return band.label;
    return "unknown";
  }

  /**
   * Build a UniProt position -> PDB residue map from an alignment.
   *
   * The alignment is the two gapped strings produced by nwAlign(uniProtSeq,
   * chainSeq). Only positions where both sides are aligned to each other AND
   * the amino acids are identical are kept, and only positions inside the
   * SIFTS-constructed range, because Needleman-Wunsch can otherwise find an
   * equally-scoring placement for a repeat elsewhere in a long protein.
   *
   * @param {string} alignedUniProt  gapped UniProt sequence
   * @param {string} alignedChain   gapped chain sequence
   * @param {Array<{resi: number, resn: string}>} chainResidues
   * @param {number} minPos  lowest UniProt position in the construct (SIFTS)
   * @param {number} maxPos  highest UniProt position in the construct (SIFTS)
   * @returns {Map<number, {resi: number, resn: string}>}
   */
  function buildUniprotMap(alignedUniProt, alignedChain, chainResidues, minPos, maxPos) {
    const map = new Map();
    if (typeof alignedUniProt !== "string" || typeof alignedChain !== "string") return map;
    const len = Math.min(alignedUniProt.length, alignedChain.length);
    let u = 0, c = 0;
    for (let i = 0; i < len; i++) {
      const ua = alignedUniProt[i], cb = alignedChain[i];
      if (ua !== "-") u++;
      if (cb !== "-") c++;
      if (ua === "-" || cb === "-") continue;
      if (ua !== cb) continue;                     // only exact matches
      if (Number.isFinite(minPos) && u < minPos) continue;
      if (Number.isFinite(maxPos) && u > maxPos) continue;
      const res = chainResidues[c - 1];
      if (!res) continue;
      map.set(u, { resi: res.resi, resn: res.resn });
    }
    return map;
  }

  /**
   * Merge variant records that describe the same protein change, keeping the
   * most severe significance and the lowest allele frequency seen.
   *
   * @param {Array<object>} records  each { pos, wt, mut, significance, af, ... }
   * @returns {Array<object>}
   */
  function mergeByChange(records) {
    const byKey = new Map();
    for (const r of records || []) {
      if (!r || typeof r.pos !== "number") continue;
      const key = `${r.wt}${r.pos}${r.mut}`;
      const prev = byKey.get(key);
      if (!prev) {
        byKey.set(key, { ...r });
        continue;
      }
      if (significanceRank(r.significance) < significanceRank(prev.significance)) {
        prev.significance = r.significance;
        prev.review = r.review || prev.review;
        prev.tier = r.tier || prev.tier;
        prev.oncogenicity = r.oncogenicity || prev.oncogenicity;
      }
      if (typeof r.af === "number" && (typeof prev.af !== "number" || r.af < prev.af)) {
        prev.af = r.af;
        prev.afSource = r.afSource;
      }
      if (!prev.id && r.id) prev.id = r.id;
      if (r.gnomadId && !prev.gnomadId) prev.gnomadId = r.gnomadId;
    }
    return [...byKey.values()];
  }

  const api = {
    SIGNIFICANCE_ORDER,
    SIGNIFICANCE_COLORS,
    IMPACT_TIERS,
    ONCOGENICITY,
    RARITY_BANDS,
    AA3_TO_1,
    parseProteinChange,
    extractProteinChange,
    parseSpdi,
    gnomadVariantId,
    classifySignificance,
    classifyImpactTier,
    classifyOncogenicity,
    significanceRank,
    rarityLabel,
    buildUniprotMap,
    mergeByChange,
  };

  global.PFE = global.PFE || {};
  global.PFE.variants = Object.freeze(api);
})(typeof globalThis !== "undefined" ? globalThis : this);
