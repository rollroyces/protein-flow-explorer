/* ===================================================================
 * Protein Flow Explorer — validation module
 *
 * Pure boundary validation: what the user typed into the search boxes,
 * and what the external APIs answered.
 *
 * Two kinds of check live here, and they fail in different ways:
 *
 *   1. Input-boundary validators (PDB ID, gene symbol, UniProt
 *      accession). Whatever they accept is interpolated into a URL and
 *      into an HTML attribute, so they are a security boundary as well
 *      as a UX one. They match a whitelist and never a denylist, so a
 *      string containing `"><script>` cannot get through: it is simply
 *      not a gene symbol.
 *
 *   2. Payload-shape validators (ClinVar, gnomAD, RCSB, UniProt,
 *      FASTA). Each one checks the fields the app actually reads --
 *      not the fields the API documentation promises -- and names the
 *      first thing it cannot use, so the caller can print a message
 *      instead of crashing on `undefined.description`.
 *
 * The convention is the same for every function here:
 *
 *      { ok: true }                               the value is usable
 *      { ok: true, value: "4HHB" }                ... and normalised
 *      { ok: false, code: "invalid", error: "..." }
 *
 * `error` is written for a human, in the app's own voice. `code` is a
 * stable machine token, so a caller can branch on *why* something
 * failed: gnomAD's "rate limited" and "no such variant" mean very
 * different things to a user, and only the code tells them apart.
 *
 * No DOM, no network, no app state — unit-testable in isolation.
 *
 * Loaded as a classic <script> before app.js and published on the
 * single PFE namespace.
 * =================================================================== */
(function (global) {
  "use strict";

  // Stable machine tokens for the `code` field of a failed validation.
  const CODES = Object.freeze({
    INVALID: "invalid",              // the input is not a valid identifier
    INVALID_SHAPE: "invalid_shape",  // the payload is not shaped the way the app needs
    RATE_LIMITED: "rate_limited",    // gnomAD throttled the request
    NOT_FOUND: "not_found",          // the API has no such record
  });

  // A 4-character PDB ID: a digit followed by three letters or digits.
  const PDB_SHORT = /^[0-9][A-Z0-9]{3}$/;
  // The extended form, "pdb_00001abc": five digits then three characters
  // of RCSB's base-32 alphabet (0-9, a-v). The canonical form is lower
  // case, so the match is case-insensitive and the value is normalised.
  const PDB_EXTENDED = /^pdb_[0-9]{5}[0-9a-v]{3}$/i;
  // Human gene symbols: upper-case letters, digits and hyphens. The
  // case-insensitive twin is checked first, so that "EGFR&foo=bar" is
  // reported as the bad character it is rather than as a case problem.
  const GENE_SYMBOL = /^[A-Z][A-Z0-9]*(?:-[A-Z0-9]+)*$/;
  const GENE_CHARS = /^[A-Za-z0-9]+(?:-[A-Za-z0-9]+)*$/;
  // UniProt's own accession pattern, six or ten characters.
  const UNIPROT_ACCESSION =
    /^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$/;
  // The isoform suffix form, "P00533-2".
  const UNIPROT_ISOFORM = /^(.+)-([0-9]+)$/;

  // Ceilings, so that a mistyped paste of a whole file is rejected
  // rather than treated as an identifier.
  const MAX_IDENTIFIER_LENGTH = 64;
  const MAX_GENE_LENGTH = 32;
  const MAX_ECHO_LENGTH = 32;

  /**
   * Name the type of a value, for use in an error message.
   * @param {*} value
   * @returns {string}
   */
  function describeType(value) {
    if (value === null) return "null";
    if (Array.isArray(value)) return "an array";
    return typeof value;
  }

  /**
   * Render an untrusted value for inclusion in an error message.
   *
   * The message may end up in a DOM node, so the echo is stripped of
   * everything that could close an attribute or open a tag, and
   * clipped. The result can only contain letters, digits and a short
   * list of punctuation, so it is safe to quote.
   *
   * @param {*} value
   * @returns {string}
   */
  function displayValue(value) {
    const s = String(value);
    const clipped = s.length > MAX_ECHO_LENGTH ? s.slice(0, MAX_ECHO_LENGTH) + "…" : s;
    return clipped.replace(/[^A-Za-z0-9 ._:+-]/g, "?");
  }

  /**
   * A plain object: not null, not an array.
   * @param {*} value
   * @returns {boolean}
   */
  /**
   * A genomic coordinate: a number, or a string of digits.
   *
   * ClinVar's ESummary returns `variation_loc[].start` and `.stop` as
   * strings ("55181378"), so a strict number check rejects real records.
   */
  function isCoordinate(value) {
    if (typeof value === "number") return Number.isFinite(value);
    if (typeof value === "string") return /^\d+$/.test(value.trim());
    return false;
  }

  function isPlainObject(value) {
    return !!value && typeof value === "object" && !Array.isArray(value);
  }

  /**
   * A string with at least one non-whitespace character.
   * @param {*} value
   * @returns {boolean}
   */
  function isFilledString(value) {
    return typeof value === "string" && value.trim() !== "";
  }

  /**
   * A failed validation.
   * @param {string} code  one of CODES
   * @param {string} error  a message for a human
   * @returns {{ok: false, code: string, error: string}}
   */
  function fail(code, error) {
    return { ok: false, code, error };
  }

  /**
   * A successful validation, optionally carrying a normalised value.
   * @param {*} [value]
   * @returns {{ok: true, value?: *}}
   */
  function pass(value) {
    return value === undefined ? { ok: true } : { ok: true, value };
  }

  /**
   * Validate a PDB entry identifier.
   *
   * Two forms are accepted:
   *   - the classic 4-character ID ("1M17", "4AKE", "6VXX"), which the
   *     PDB numbers as a digit followed by three letters or digits;
   *   - the extended ID ("pdb_00001abc"), which the PDB adopted for
   *     entries beyond the 4-character range.
   *
   * Case is not significant here -- the app upper-cases the ID before
   * building a URL -- so "1m17" is accepted and normalised to "1M17".
   * Everything else is rejected: an empty string, a string with a space
   * or a URL metacharacter in it, and a non-string.
   *
   * @param {string} value  the text the user typed
   * @returns {{ok: true, value: string}|{ok: false, code: string, error: string}}
   */
  function validatePdbId(value) {
    if (typeof value !== "string") {
      return fail(CODES.INVALID,
        `Invalid PDB ID: expected a string, got ${describeType(value)}.`);
    }
    const id = value.trim();
    if (!id) {
      return fail(CODES.INVALID, "Invalid PDB ID: the value is empty.");
    }
    if (id.length > MAX_IDENTIFIER_LENGTH) {
      return fail(CODES.INVALID,
        `Invalid PDB ID: too long (${id.length} characters, maximum ${MAX_IDENTIFIER_LENGTH}).`);
    }
    if (PDB_EXTENDED.test(id)) {
      return pass(id.toLowerCase());
    }
    const upper = id.toUpperCase();
    if (upper.length !== 4) {
      return fail(CODES.INVALID,
        `Invalid PDB ID "${displayValue(id)}": expected 4 characters, e.g. 4HHB.`);
    }
    if (!/^[0-9A-Z]{4}$/.test(upper)) {
      return fail(CODES.INVALID,
        `Invalid PDB ID "${displayValue(id)}": only letters and digits are allowed, e.g. 4HHB.`);
    }
    if (!PDB_SHORT.test(upper)) {
      return fail(CODES.INVALID,
        `Invalid PDB ID "${displayValue(id)}": a 4-character PDB ID starts with a digit, e.g. 4HHB.`);
    }
    return pass(upper);
  }

  /**
   * Validate a human gene symbol.
   *
   * HGNC symbols are upper-case letters, digits and hyphens ("EGFR",
   * "TP53", "HBB", "CALM1", "SCN5A", "HLA-DRB1"). The symbol is
   * interpolated into a UniProt query string, so it is matched against
   * that whitelist rather than a list of forbidden characters: an
   * `&`, a `#`, a `/`, a `?`, a quote, an angle bracket or a space all
   * fail simply by not being gene symbols.
   *
   * Lower case is rejected rather than folded, because the symbol is
   * defined in upper case and the error can say what was meant.
   *
   * @param {string} value  the text the user typed
   * @returns {{ok: true, value: string}|{ok: false, code: string, error: string}}
   */
  function validateGeneSymbol(value) {
    if (typeof value !== "string") {
      return fail(CODES.INVALID,
        `Invalid gene symbol: expected a string, got ${describeType(value)}.`);
    }
    const symbol = value.trim();
    if (!symbol) {
      return fail(CODES.INVALID, "Invalid gene symbol: the value is empty.");
    }
    if (symbol.length > MAX_GENE_LENGTH) {
      return fail(CODES.INVALID,
        `Invalid gene symbol "${displayValue(symbol)}": too long ` +
        `(${symbol.length} characters, maximum ${MAX_GENE_LENGTH}).`);
    }
    if (!GENE_CHARS.test(symbol)) {
      return fail(CODES.INVALID,
        `Invalid gene symbol "${displayValue(symbol)}": only letters, digits and ` +
        `hyphens are allowed, e.g. EGFR.`);
    }
    if (/[a-z]/.test(symbol)) {
      return fail(CODES.INVALID,
        `Invalid gene symbol "${displayValue(symbol)}": gene symbols are upper case, e.g. EGFR.`);
    }
    if (!GENE_SYMBOL.test(symbol)) {
      return fail(CODES.INVALID,
        `Invalid gene symbol "${displayValue(symbol)}": a gene symbol starts with ` +
        `a letter, e.g. EGFR.`);
    }
    return pass(symbol);
  }

  /**
   * Validate a UniProt accession.
   *
   * The canonical forms are six or ten characters ("P00533", "P69905",
   * "A0A087X1C4") and the isoform form adds a numeric suffix
   * ("P00533-2"). The match uses UniProt's own accession pattern rather
   * than "any six upper-case characters", so "AAAAAA" is rejected: it
   * is not an accession anybody has.
   *
   * @param {string} value  the text the user typed
   * @returns {{ok: true, value: string}|{ok: false, code: string, error: string}}
   */
  function validateUniprotAccession(value) {
    if (typeof value !== "string") {
      return fail(CODES.INVALID,
        `Invalid UniProt accession: expected a string, got ${describeType(value)}.`);
    }
    const accession = value.trim().toUpperCase();
    if (!accession) {
      return fail(CODES.INVALID, "Invalid UniProt accession: the value is empty.");
    }
    if (accession.length > MAX_IDENTIFIER_LENGTH) {
      return fail(CODES.INVALID,
        `Invalid UniProt accession: too long (${accession.length} characters, ` +
        `maximum ${MAX_IDENTIFIER_LENGTH}).`);
    }
    if (UNIPROT_ACCESSION.test(accession)) {
      return pass(accession);
    }
    if (!/^[A-Z0-9-]+$/.test(accession)) {
      return fail(CODES.INVALID,
        `Invalid UniProt accession "${displayValue(accession)}": only letters, digits ` +
        `and hyphens are allowed, e.g. P00533.`);
    }
    const iso = accession.match(UNIPROT_ISOFORM);
    if (iso && UNIPROT_ACCESSION.test(iso[1])) {
      return pass(accession);
    }
    if (/^[A-Z0-9]+-[0-9]+$/.test(accession)) {
      return fail(CODES.INVALID,
        `Invalid UniProt accession "${displayValue(accession)}": the isoform suffix ` +
        `is not a known accession, e.g. P00533-2.`);
    }
    if (/^[A-Z]{6,10}$/.test(accession)) {
      return fail(CODES.INVALID,
        `Invalid UniProt accession "${displayValue(accession)}": not a UniProt ` +
        `accession pattern, e.g. P00533.`);
    }
    return fail(CODES.INVALID,
      `Invalid UniProt accession "${displayValue(accession)}": expected 6 or 10 ` +
      `characters, e.g. P00533 or P00533-2.`);
  }

  /**
   * Validate one ClinVar ESummary record.
   *
   * Required, because the app cannot place the record on a structure
   * without them:
   *   - an identity: `accession_version` or `uid`;
   *   - `germline_classification`, with a non-empty `description`.
   *     This is the *current* field. ClinVar no longer publishes
   *     `clinical_significance`, so a record that still carries only
   *     the legacy field is called out by name rather than being
   *     silently classified as a VUS;
   *   - a non-empty `variation_set`, whose first entry names the
   *     variant on the canonical transcript;
   *   - a protein change, from either `variation_set[0].variation_name`
   *     or the `protein_change` field.
   *
   * Optional, but checked for type when present:
   *   - `clinical_impact_classification` and `oncogenicity_classification`
   *     (somatic tiers, not every record has them);
   *   - `variation_set[0].canonical_spdi` and `variation_set[0].variation_loc`
   *     (the gnomAD lookup falls back to the GRCh38 location);
   *   - `genes[0].symbol`, `title`, `protein_change`.
   *
   * @param {object} record  one entry of a ClinVar ESummary result
   * @returns {{ok: true}|{ok: false, code: string, error: string}}
   */
  function validateClinVarRecord(record) {
    if (!isPlainObject(record)) {
      return fail(CODES.INVALID_SHAPE,
        `ClinVar record must be an object, got ${describeType(record)}.`);
    }
    const id = isFilledString(record.accession_version)
      ? record.accession_version.trim()
      : (isFilledString(record.uid) ? String(record.uid).trim() : "");
    const label = id ? `ClinVar record ${id}` : "ClinVar record";

    if (!id) {
      return fail(CODES.INVALID_SHAPE, `${label} has no accession_version or uid.`);
    }
    if (!isPlainObject(record.germline_classification)) {
      if (isPlainObject(record.clinical_significance)) {
        return fail(CODES.INVALID_SHAPE,
          `${label} has only the legacy clinical_significance field; ClinVar ` +
          `publishes germline_classification.`);
      }
      return fail(CODES.INVALID_SHAPE, `${label} has no germline_classification.`);
    }
    const gc = record.germline_classification;
    if (!isFilledString(gc.description)) {
      return fail(CODES.INVALID_SHAPE, `${label} has no germline_classification.description.`);
    }
    if (gc.review_status !== undefined && typeof gc.review_status !== "string") {
      return fail(CODES.INVALID_SHAPE,
        `${label} has a germline_classification.review_status that is not a string.`);
    }
    if (record.clinical_impact_classification !== undefined) {
      const ic = record.clinical_impact_classification;
      if (!isPlainObject(ic)) {
        return fail(CODES.INVALID_SHAPE,
          `${label} has a clinical_impact_classification that is not an object.`);
      }
      if (ic.description !== undefined && typeof ic.description !== "string") {
        return fail(CODES.INVALID_SHAPE,
          `${label} has a clinical_impact_classification.description that is not a string.`);
      }
    }
    if (record.oncogenicity_classification !== undefined) {
      const oc = record.oncogenicity_classification;
      if (!isPlainObject(oc)) {
        return fail(CODES.INVALID_SHAPE,
          `${label} has an oncogenicity_classification that is not an object.`);
      }
      if (oc.description !== undefined && typeof oc.description !== "string") {
        return fail(CODES.INVALID_SHAPE,
          `${label} has an oncogenicity_classification.description that is not a string.`);
      }
    }
    if (!Array.isArray(record.variation_set)) {
      return fail(CODES.INVALID_SHAPE, `${label} has no variation_set array.`);
    }
    if (record.variation_set.length === 0) {
      return fail(CODES.INVALID_SHAPE, `${label} has an empty variation_set.`);
    }
    const vs = record.variation_set[0];
    if (!isPlainObject(vs)) {
      return fail(CODES.INVALID_SHAPE, `${label} has a variation_set[0] that is not an object.`);
    }
    if (vs.variation_name !== undefined && typeof vs.variation_name !== "string") {
      return fail(CODES.INVALID_SHAPE,
        `${label} has a variation_set[0].variation_name that is not a string.`);
    }
    if (vs.canonical_spdi !== undefined && typeof vs.canonical_spdi !== "string") {
      return fail(CODES.INVALID_SHAPE,
        `${label} has a variation_set[0].canonical_spdi that is not a string.`);
    }
    if (vs.variation_loc !== undefined) {
      if (!Array.isArray(vs.variation_loc)) {
        return fail(CODES.INVALID_SHAPE,
          `${label} has a variation_set[0].variation_loc that is not an array.`);
      }
      for (const loc of vs.variation_loc) {
        if (!isPlainObject(loc)) {
          return fail(CODES.INVALID_SHAPE,
            `${label} has a variation_set[0].variation_loc entry that is not an object.`);
        }
        if (loc.assembly_name !== undefined && typeof loc.assembly_name !== "string") {
          return fail(CODES.INVALID_SHAPE,
            `${label} has a variation_set[0].variation_loc assembly_name that is not a string.`);
        }
        if (loc.chr !== undefined && typeof loc.chr !== "string" && typeof loc.chr !== "number") {
          return fail(CODES.INVALID_SHAPE,
            `${label} has a variation_set[0].variation_loc chr that is not a string or number.`);
        }
        // ClinVar returns the coordinates as STRINGS ("55181378"), not
        // numbers. Requiring a number rejected every real record, which
        // silently emptied the variant map. Accept either, and reject
        // anything that is not an integer.
        if (loc.start !== undefined && !isCoordinate(loc.start)) {
          return fail(CODES.INVALID_SHAPE,
            `${label} has a variation_set[0].variation_loc start that is not a coordinate.`);
        }
        if (loc.stop !== undefined && !isCoordinate(loc.stop)) {
          return fail(CODES.INVALID_SHAPE,
            `${label} has a variation_set[0].variation_loc stop that is not a coordinate.`);
        }
      }
    }
    if (record.protein_change !== undefined && typeof record.protein_change !== "string") {
      return fail(CODES.INVALID_SHAPE, `${label} has a protein_change that is not a string.`);
    }
    if (!isFilledString(vs.variation_name) && !isFilledString(record.protein_change)) {
      return fail(CODES.INVALID_SHAPE,
        `${label} has neither a variation_name nor a protein_change.`);
    }
    if (record.genes !== undefined) {
      if (!Array.isArray(record.genes)) {
        return fail(CODES.INVALID_SHAPE, `${label} has a genes field that is not an array.`);
      }
      const gene = record.genes[0];
      if (gene !== undefined) {
        if (!isPlainObject(gene)) {
          return fail(CODES.INVALID_SHAPE, `${label} has a genes[0] that is not an object.`);
        }
        if (gene.symbol !== undefined && typeof gene.symbol !== "string") {
          return fail(CODES.INVALID_SHAPE, `${label} has a genes[0].symbol that is not a string.`);
        }
      }
    }
    if (record.title !== undefined && typeof record.title !== "string") {
      return fail(CODES.INVALID_SHAPE, `${label} has a title that is not a string.`);
    }
    return pass();
  }

  /**
   * Validate a gnomAD GraphQL response.
   *
   * gnomAD answers a missing variant with an *error*, not an empty
   * result, and answers a throttled request with HTTP 400 and the
   * message "Query rate limit exceeded". The two mean different things
   * to a user -- "we do not have this variant" versus "we were not
   * allowed to ask" -- so they are given different codes.
   *
   * The allele frequency itself is required from at least one of the
   * two datasets, because that is the number the app puts on the
   * structure.
   *
   * @param {object} payload  the parsed GraphQL response body
   * @returns {{ok: true}|{ok: false, code: string, error: string}}
   */
  function validateGnomadResponse(payload) {
    if (!isPlainObject(payload)) {
      return fail(CODES.INVALID_SHAPE,
        `gnomAD response must be an object, got ${describeType(payload)}.`);
    }
    let message = "";
    if (Array.isArray(payload.errors) && payload.errors.length) {
      const first = payload.errors[0];
      if (isPlainObject(first) && typeof first.message === "string") {
        message = first.message.trim();
      }
    }
    if (message) {
      if (/rate limit|too many/i.test(message)) {
        return fail(CODES.RATE_LIMITED, `gnomAD rate limit exceeded: ${message}.`);
      }
      if (/not found|no such variant/i.test(message)) {
        return fail(CODES.NOT_FOUND, `gnomAD has no record for this variant: ${message}.`);
      }
      return fail(CODES.INVALID_SHAPE, `gnomAD returned an error: ${message}.`);
    }
    if (!isPlainObject(payload.data)) {
      return fail(CODES.INVALID_SHAPE, "gnomAD response has no data object.");
    }
    const variant = payload.data.variant;
    if (variant === null || variant === undefined) {
      return fail(CODES.NOT_FOUND, "gnomAD returned no variant data.");
    }
    if (!isPlainObject(variant)) {
      return fail(CODES.INVALID_SHAPE, "gnomAD response has a data.variant that is not an object.");
    }
    if (!isFilledString(variant.variant_id)) {
      return fail(CODES.INVALID_SHAPE, "gnomAD response has no data.variant.variant_id.");
    }
    for (const dataset of ["exome", "genome"]) {
      const d = variant[dataset];
      if (d === undefined || d === null) continue;
      if (!isPlainObject(d)) {
        return fail(CODES.INVALID_SHAPE,
          `gnomAD response has a data.variant.${dataset} that is not an object.`);
      }
      for (const field of ["ac", "an", "af"]) {
        if (d[field] === undefined || d[field] === null) continue;
        if (typeof d[field] !== "number" || !Number.isFinite(d[field])) {
          return fail(CODES.INVALID_SHAPE,
            `gnomAD response has a data.variant.${dataset}.${field} that is not a finite number.`);
        }
      }
    }
    const hasFrequency = ["genome", "exome"].some((dataset) => {
      const d = variant[dataset];
      return isPlainObject(d) && typeof d.af === "number" && Number.isFinite(d.af);
    });
    if (!hasFrequency) {
      return fail(CODES.NOT_FOUND, "gnomAD has no allele frequency for this variant.");
    }
    return pass();
  }

  /**
   * Validate an RCSB polymer_entity response.
   *
   * The app reads three things from it: the UniProt reference in
   * `rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers`,
   * the SIFTS alignment in `rcsb_polymer_entity_align[].aligned_regions`,
   * and the chain list in `entity_poly.pdbx_strand_id`.
   *
   * None of the three is guaranteed -- a synthetic construct has no
   * UniProt reference, and not every entity is aligned -- so the
   * response must carry at least one of them, and each one that is
   * present must have the right shape.
   *
   * @param {object} payload  the parsed RCSB polymer_entity body
   * @returns {{ok: true}|{ok: false, code: string, error: string}}
   */
  function validateRcsbEntity(payload) {
    if (!isPlainObject(payload)) {
      return fail(CODES.INVALID_SHAPE,
        `RCSB polymer entity response must be an object, got ${describeType(payload)}.`);
    }
    const hasIds = isPlainObject(payload.rcsb_polymer_entity_container_identifiers);
    const hasAlign = Array.isArray(payload.rcsb_polymer_entity_align);
    const hasPoly = isPlainObject(payload.entity_poly);
    if (!hasIds && !hasAlign && !hasPoly) {
      return fail(CODES.INVALID_SHAPE,
        "RCSB polymer entity response has no rcsb_polymer_entity_container_identifiers, " +
        "rcsb_polymer_entity_align or entity_poly.");
    }
    if (payload.rcsb_polymer_entity_container_identifiers !== undefined && !hasIds) {
      return fail(CODES.INVALID_SHAPE,
        "RCSB polymer entity rcsb_polymer_entity_container_identifiers is not an object.");
    }
    if (hasIds) {
      const refs = payload.rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers;
      if (refs !== undefined) {
        if (!Array.isArray(refs)) {
          return fail(CODES.INVALID_SHAPE,
            "RCSB polymer entity reference_sequence_identifiers is not an array.");
        }
        for (let i = 0; i < refs.length; i++) {
          const ref = refs[i];
          if (!isPlainObject(ref)) {
            return fail(CODES.INVALID_SHAPE,
              `RCSB polymer entity reference_sequence_identifiers[${i}] is not an object.`);
          }
          if (ref.database_name !== undefined && typeof ref.database_name !== "string") {
            return fail(CODES.INVALID_SHAPE,
              `RCSB polymer entity reference_sequence_identifiers[${i}].database_name ` +
              `is not a string.`);
          }
          if (ref.database_accession !== undefined && typeof ref.database_accession !== "string") {
            return fail(CODES.INVALID_SHAPE,
              `RCSB polymer entity reference_sequence_identifiers[${i}].database_accession ` +
              `is not a string.`);
          }
        }
      }
    }
    if (payload.rcsb_polymer_entity_align !== undefined && !hasAlign) {
      return fail(CODES.INVALID_SHAPE,
        "RCSB polymer entity rcsb_polymer_entity_align is not an array.");
    }
    if (hasAlign) {
      for (let i = 0; i < payload.rcsb_polymer_entity_align.length; i++) {
        const al = payload.rcsb_polymer_entity_align[i];
        if (!isPlainObject(al)) {
          return fail(CODES.INVALID_SHAPE,
            `RCSB polymer entity rcsb_polymer_entity_align[${i}] is not an object.`);
        }
        if (al.reference_database_name !== undefined && typeof al.reference_database_name !== "string") {
          return fail(CODES.INVALID_SHAPE,
            `RCSB polymer entity rcsb_polymer_entity_align[${i}].reference_database_name ` +
            `is not a string.`);
        }
        if (al.reference_database_accession !== undefined &&
            typeof al.reference_database_accession !== "string") {
          return fail(CODES.INVALID_SHAPE,
            `RCSB polymer entity rcsb_polymer_entity_align[${i}].reference_database_accession ` +
            `is not a string.`);
        }
        const regions = al.aligned_regions;
        if (regions === undefined) continue;
        if (!Array.isArray(regions)) {
          return fail(CODES.INVALID_SHAPE,
            `RCSB polymer entity rcsb_polymer_entity_align[${i}].aligned_regions ` +
            `is not an array.`);
        }
        for (let j = 0; j < regions.length; j++) {
          const reg = regions[j];
          if (!isPlainObject(reg)) {
            return fail(CODES.INVALID_SHAPE,
              `RCSB polymer entity rcsb_polymer_entity_align[${i}].aligned_regions[${j}] ` +
              `is not an object.`);
          }
          if (reg.ref_beg_seq_id !== undefined &&
              (typeof reg.ref_beg_seq_id !== "number" || !Number.isFinite(reg.ref_beg_seq_id))) {
            return fail(CODES.INVALID_SHAPE,
              `RCSB polymer entity rcsb_polymer_entity_align[${i}].aligned_regions[${j}]` +
              `.ref_beg_seq_id is not a number.`);
          }
          if (reg.length !== undefined &&
              (typeof reg.length !== "number" || !Number.isFinite(reg.length))) {
            return fail(CODES.INVALID_SHAPE,
              `RCSB polymer entity rcsb_polymer_entity_align[${i}].aligned_regions[${j}]` +
              `.length is not a number.`);
          }
        }
      }
    }
    if (payload.entity_poly !== undefined && !hasPoly) {
      return fail(CODES.INVALID_SHAPE, "RCSB polymer entity entity_poly is not an object.");
    }
    if (hasPoly && payload.entity_poly.pdbx_strand_id !== undefined &&
        typeof payload.entity_poly.pdbx_strand_id !== "string") {
      return fail(CODES.INVALID_SHAPE,
        "RCSB polymer entity entity_poly.pdbx_strand_id is not a string.");
    }
    return pass();
  }

  /**
   * Validate a UniProt JSON response.
   *
   * The app reads exactly two things from it: the gene symbol in
   * `genes[0].geneName.value`, which names the ClinVar query, and the
   * canonical sequence in `sequence.value`, which the variants are
   * mapped onto. Without either one the variant pipeline cannot run,
   * so both are required.
   *
   * @param {object} payload  the parsed UniProt entry body
   * @returns {{ok: true}|{ok: false, code: string, error: string}}
   */
  function validateUniProtResponse(payload) {
    if (!isPlainObject(payload)) {
      return fail(CODES.INVALID_SHAPE,
        `UniProt response must be an object, got ${describeType(payload)}.`);
    }
    if (!isPlainObject(payload.sequence)) {
      return fail(CODES.INVALID_SHAPE, "UniProt response has no sequence object.");
    }
    if (!isFilledString(payload.sequence.value)) {
      return fail(CODES.INVALID_SHAPE, "UniProt response has no sequence.value.");
    }
    if (!Array.isArray(payload.genes)) {
      return fail(CODES.INVALID_SHAPE, "UniProt response has no genes array.");
    }
    const gene = payload.genes[0];
    if (!isPlainObject(gene)) {
      return fail(CODES.INVALID_SHAPE, "UniProt response has no genes[0] object.");
    }
    if (!isPlainObject(gene.geneName)) {
      return fail(CODES.INVALID_SHAPE, "UniProt response has no genes[0].geneName object.");
    }
    if (!isFilledString(gene.geneName.value)) {
      return fail(CODES.INVALID_SHAPE, "UniProt response has no genes[0].geneName.value.");
    }
    return pass();
  }

  /**
   * Validate a FASTA response (UniProt, AlphaFold or any other source).
   *
   * The body has to be FASTA text that `parseFastaSequences` (in
   * src/seq.js) can read: a header line starting with ">", and at
   * least one non-empty amino-acid sequence. An HTML error page is
   * rejected explicitly, because a CDN that answers a missing file
   * with HTML rather than a 404 would otherwise be parsed as a
   * one-letter sequence.
   *
   * @param {string} text  the response body
   * @returns {{ok: true}|{ok: false, code: string, error: string}}
   */
  function validateFasta(text, source) {
    const what = source || "FASTA response";
    if (typeof text !== "string") {
      return fail(CODES.INVALID_SHAPE,
        `${what} must be text, got ${describeType(text)}.`);
    }
    const body = text.trim();
    if (!body) {
      return fail(CODES.INVALID_SHAPE, `${what} is empty.`);
    }
    if (/^</.test(body) || /<!DOCTYPE/i.test(body)) {
      return fail(CODES.INVALID_SHAPE,
        `${what} is an HTML page, not FASTA.`);
    }
    const lines = body.split(/\r?\n/);
    if (!lines[0].startsWith(">")) {
      return fail(CODES.INVALID_SHAPE,
        `${what} has no '>' header line.`);
    }
    let records = 0;
    let residues = 0;
    let inRecord = false;
    for (const line of lines) {
      if (line.startsWith(">")) {
        records++;
        inRecord = true;
        continue;
      }
      if (!inRecord) continue;
      const seq = line.trim();
      if (!seq) continue;
      if (!/^[A-Za-z*-]+$/.test(seq)) {
        return fail(CODES.INVALID_SHAPE,
          `${what} contains a character that is not an amino acid.`);
      }
      residues += seq.length;
    }
    if (records === 0) {
      return fail(CODES.INVALID_SHAPE, `${what} has no FASTA record.`);
    }
    if (residues === 0) {
      return fail(CODES.INVALID_SHAPE, `${what} has an empty sequence.`);
    }
    return pass();
  }

  /**
   * Validate a list of payloads with one validator and collect the
   * failures, so a batch of records can be summarised in one message
   * rather than one message per record.
   *
   * @param {Array} items  the payloads to validate
   * @param {function} validator  a validator returning { ok, code, error }
   * @param {function} [identify]  optional (item, index) => label, used
   *        to say which record failed
   * @returns {{ok: boolean, count: number, failures: Array<{index: number,
   *          id: string|null, code: string, error: string}>}}
   */
  function validateList(items, validator, identify) {
    if (!Array.isArray(items)) {
      return {
        ok: false,
        code: CODES.INVALID,
        error: `Expected a list of records, got ${describeType(items)}.`,
        count: 0,
        failures: [],
      };
    }
    if (typeof validator !== "function") {
      return {
        ok: false,
        code: CODES.INVALID,
        error: "validateList needs a validator function.",
        count: 0,
        failures: [],
      };
    }
    const failures = [];
    for (let i = 0; i < items.length; i++) {
      const result = validator(items[i], i);
      if (result && result.ok) continue;
      const error = result && result.error ? result.error : "invalid record";
      const code = result && result.code ? result.code : CODES.INVALID_SHAPE;
      let id = null;
      if (typeof identify === "function") {
        const label = identify(items[i], i);
        if (label) id = String(label);
      }
      failures.push({ index: i, id, code, error });
    }
    if (failures.length === 0) {
      return { ok: true, count: items.length, failures };
    }
    return {
      ok: false,
      code: failures[0].code,
      error: `${failures.length} of ${items.length} records failed validation: ` +
        failures[0].error,
      count: items.length,
      failures,
    };
  }

  /**
   * Validate a batch of ClinVar records and collect the failures.
   *
   * @param {Array<object>} records  the records of one ESummary result
   * @returns {{ok: boolean, count: number, failures: Array<{index: number,
   *          id: string|null, code: string, error: string}>}}
   */
  function validateClinVarRecords(records) {
    return validateList(records, validateClinVarRecord, (record) => {
      if (!isPlainObject(record)) return null;
      return record.accession_version || record.uid || null;
    });
  }

  const api = {
    CODES,
    validatePdbId,
    validateGeneSymbol,
    validateUniprotAccession,
    validateClinVarRecord,
    validateClinVarRecords,
    validateGnomadResponse,
    validateRcsbEntity,
    validateUniProtResponse,
    validateFasta,
    validateList,
  };

  global.PFE = global.PFE || {};
  global.PFE.validate = Object.freeze(api);
})(typeof globalThis !== "undefined" ? globalThis : this);
