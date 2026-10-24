/**
 * Unit tests for the validation module.
 *
 *   node test/validate.test.js
 *
 * Two boundaries are covered here.
 *
 *   1. The input boundary: what the user typed into the "PDB ID" and
 *      "gene symbol" boxes. Whatever passes is interpolated into a
 *      URL and into an HTML attribute, so the validators are a
 *      security boundary as well as a UX one. The hostile inputs below
 *      are the ones that would corrupt a query string, walk a path, or
 *      inject markup.
 *
 *   2. The payload boundary: what RCSB, UniProt, ClinVar, gnomAD and
 *      AlphaFold answered. Every fixture here is inline: the tests
 *      make no network calls at all.
 *
 * The module is pure, so it loads with a bare require()-style probe
 * and no DOM stub -- that property is asserted explicitly below.
 */
const fs = require("fs");
const path = require("path");
const vm = require("vm");

const ROOT = path.join(__dirname, "..");
vm.runInThisContext(fs.readFileSync(path.join(ROOT, "src", "validate.js"), "utf8"),
  { filename: "src/validate.js" });
const V = globalThis.PFE.validate;

let pass = 0, fail = 0;
function check(name, cond, detail = "") {
  if (cond) { pass++; console.log(`  PASS  ${name}${detail ? "  " + detail : ""}`); }
  else { fail++; console.log(`  FAIL  ${name}${detail ? "  " + detail : ""}`); }
}

// ------------------------------------------------------------------ module
console.log("\nModule isolation");
const SRC = fs.readFileSync(path.join(ROOT, "src", "validate.js"), "utf8");
check("no DOM required", typeof document === "undefined");
check("PFE.validate is published", !!V);
check("PFE.validate is frozen", Object.isFrozen(V));
check("no fetch call in the source", !/\bfetch\s*\(/.test(SRC));
check("no document access in the source", !/\bdocument\b/.test(SRC));
check("no window access in the source", !/\bwindow\b/.test(SRC));
check("no timers in the source", !/setTimeout|setInterval/.test(SRC));
check("no ES module syntax", !/^\s*(import|export)\s/m.test(SRC));
check("every expected member is exported",
  ["CODES", "validatePdbId", "validateGeneSymbol", "validateUniprotAccession",
   "validateClinVarRecord", "validateClinVarRecords", "validateGnomadResponse",
   "validateRcsbEntity", "validateUniProtResponse", "validateAlphaFoldFasta",
   "validateList"].every((k) => typeof V[k] === "function" || V[k] === V.CODES));
check("the error codes are the documented ones",
  V.CODES.INVALID === "invalid" && V.CODES.INVALID_SHAPE === "invalid_shape" &&
  V.CODES.RATE_LIMITED === "rate_limited" && V.CODES.NOT_FOUND === "not_found");

// ---------------------------------------------------------------- PDB IDs
console.log("\nPDB IDs");
check("1M17 is valid", V.validatePdbId("1M17").ok);
check("4AKE is valid", V.validatePdbId("4AKE").ok);
check("6VXX is valid", V.validatePdbId("6VXX").ok);
check("a valid ID is normalised to upper case",
  V.validatePdbId("1m17").ok && V.validatePdbId("1m17").value === "1M17");
check("surrounding whitespace is trimmed",
  V.validatePdbId("  1M17  ").ok && V.validatePdbId("  1M17  ").value === "1M17");
check("the extended ID form is accepted",
  V.validatePdbId("pdb_00001abc").ok && V.validatePdbId("pdb_00001abc").value === "pdb_00001abc");
check("an extended ID is normalised to lower case",
  V.validatePdbId("PDB_00001ABC").ok && V.validatePdbId("PDB_00001ABC").value === "pdb_00001abc");
check("a three-character ID is rejected", !V.validatePdbId("1M1").ok);
check("a five-character ID is rejected", !V.validatePdbId("1M17X").ok);
check("a letter-leading ID is rejected", !V.validatePdbId("ABHH").ok);
check("the error says what was expected", /4 characters/.test(V.validatePdbId("1M17X").error));

console.log("\nPDB IDs: hostile input");
const HOSTILE = [
  ['"><script>alert(1)</script>', "markup injection"],
  ["EGFR&foo=bar", "query-string injection"],
  ["../../etc/passwd", "path traversal"],
  ["javascript:alert(1)", "javascript scheme"],
  ["", "empty string"],
  ["   ", "whitespace only"],
  ["A".repeat(10000), "10,000-character string"],
  [null, "null"],
  [undefined, "undefined"],
  [1234, "a number instead of a string"],
  [{}, "an object instead of a string"],
  [["1M17"], "an array instead of a string"],
];
for (const [value, label] of HOSTILE) {
  check(`PDB ID rejects ${label}`, !V.validatePdbId(value).ok);
  check(`gene symbol rejects ${label}`, !V.validateGeneSymbol(value).ok);
  check(`UniProt accession rejects ${label}`, !V.validateUniprotAccession(value).ok);
}
check("a rejected PDB ID never returns a value", V.validatePdbId("1M 7").value === undefined);
check("a rejected input carries a code", V.validatePdbId("").code === V.CODES.INVALID);
check("a rejected input carries a message", V.validatePdbId("").error.length > 10);

// ------------------------------------------------------------ gene symbols
console.log("\nGene symbols");
check("EGFR is valid", V.validateGeneSymbol("EGFR").ok);
check("TP53 is valid", V.validateGeneSymbol("TP53").ok);
check("HBB is valid", V.validateGeneSymbol("HBB").ok);
check("CALM1 is valid", V.validateGeneSymbol("CALM1").ok);
check("SCN5A is valid", V.validateGeneSymbol("SCN5A").ok);
check("a hyphenated symbol is valid", V.validateGeneSymbol("HLA-DRB1").ok);
check("surrounding whitespace is trimmed",
  V.validateGeneSymbol("  EGFR  ").ok && V.validateGeneSymbol("  EGFR  ").value === "EGFR");
check("lower case is rejected", !V.validateGeneSymbol("egfr").ok);
check("the case error says what was expected",
  /upper case/.test(V.validateGeneSymbol("egfr").error));
check("a space inside the symbol is rejected", !V.validateGeneSymbol("EG FR").ok);
check("a leading digit is rejected", !V.validateGeneSymbol("1ABC").ok);
check("a doubled hyphen is rejected", !V.validateGeneSymbol("ABC--DEF").ok);
check("a trailing hyphen is rejected", !V.validateGeneSymbol("ABC-").ok);
check("an over-long symbol is rejected", !V.validateGeneSymbol("A".repeat(10000)).ok);
check("the length error says the limit", /maximum/.test(V.validateGeneSymbol("A".repeat(10000)).error));

// --------------------------------------------------------- UniProt accessions
console.log("\nUniProt accessions");
check("P00533 is valid", V.validateUniprotAccession("P00533").ok);
check("P69905 is valid", V.validateUniprotAccession("P69905").ok);
check("P04637 is valid", V.validateUniprotAccession("P04637").ok);
check("the ten-character form is valid", V.validateUniprotAccession("A0A087X1C4").ok);
check("the isoform form is valid",
  V.validateUniprotAccession("P00533-2").ok && V.validateUniprotAccession("P00533-2").value === "P00533-2");
check("lower case is folded to upper case",
  V.validateUniprotAccession("p00533").ok && V.validateUniprotAccession("p00533").value === "P00533");
check("a five-character accession is rejected", !V.validateUniprotAccession("P0053").ok);
check("a non-accession of the right length is rejected",
  !V.validateUniprotAccession("AAAAAA").ok);
check("an empty isoform suffix is rejected", !V.validateUniprotAccession("P00533-").ok);
check("a non-numeric isoform suffix is rejected", !V.validateUniprotAccession("P00533-x").ok);
check("an isoform of an unknown accession is rejected",
  !V.validateUniprotAccession("AAAAAA-2").ok);

// --------------------------------------------------------- ClinVar records
console.log("\nClinVar ESummary records");
const CLINVAR = {
  uid: "16609",
  accession_version: "RCV000016609.13",
  title: "EGFR T790M",
  variation_set: [{
    variation_name: "NM_005228.5(EGFR):c.2369C>T (p.Thr790Met)",
    canonical_spdi: "NC_000007.14:55249071:C:T",
    variation_loc: [
      { assembly_name: "GRCh38", chr: "7", start: 55249072, stop: 55249072 },
      { assembly_name: "GRCh37", chr: "7", start: 55259515, stop: 55259515 },
    ],
  }],
  germline_classification: {
    description: "drug response",
    review_status: "criteria provided, single submitter",
  },
  clinical_impact_classification: { description: "Tier I - Moderate" },
  oncogenicity_classification: { description: "Oncogenic" },
  protein_change: "T790M",
  genes: [{ symbol: "EGFR" }],
};
const clone = (o) => JSON.parse(JSON.stringify(o));

check("a complete record is valid", V.validateClinVarRecord(CLINVAR).ok);
check("a record with only the somatic blocks is still valid",
  V.validateClinVarRecord((() => {
    const r = clone(CLINVAR);
    delete r.clinical_impact_classification;
    delete r.oncogenicity_classification;
    return r;
  })()).ok);
check("a record without genes is still valid",
  V.validateClinVarRecord((() => {
    const r = clone(CLINVAR);
    delete r.genes;
    return r;
  })()).ok);
check("a record without variation_loc is still valid",
  V.validateClinVarRecord((() => {
    const r = clone(CLINVAR);
    delete r.variation_set[0].variation_loc;
    return r;
  })()).ok);

check("a null payload is rejected", !V.validateClinVarRecord(null).ok);
check("an array payload is rejected", !V.validateClinVarRecord([CLINVAR]).ok);
check("a string payload is rejected", !V.validateClinVarRecord("{}").ok);
check("an empty object is rejected", !V.validateClinVarRecord({}).ok);
check("the empty-object error names the missing identity",
  /accession_version/.test(V.validateClinVarRecord({}).error));

check("germline_classification present but empty is rejected",
  !V.validateClinVarRecord((() => {
    const r = clone(CLINVAR);
    r.germline_classification = {};
    return r;
  })()).ok);
check("the empty-classification error names the field",
  /germline_classification\.description/.test(
    V.validateClinVarRecord((() => {
      const r = clone(CLINVAR);
      r.germline_classification = {};
      return r;
    })()).error));

check("a record with only the legacy clinical_significance is rejected",
  !V.validateClinVarRecord({
    accession_version: "RCV000000001.1",
    clinical_significance: { description: "Pathogenic", review_status: "no assertion" },
    variation_set: [{ variation_name: "NM_005228.5(EGFR):c.2369C>T (p.Thr790Met)" }],
  }).ok);
check("the legacy-field error names both fields",
  /clinical_significance/.test(V.validateClinVarRecord({
    accession_version: "RCV000000001.1",
    clinical_significance: { description: "Pathogenic" },
  }).error) && /germline_classification/.test(V.validateClinVarRecord({
    accession_version: "RCV000000001.1",
    clinical_significance: { description: "Pathogenic" },
  }).error));

check("a record with no germline_classification is rejected",
  !V.validateClinVarRecord((() => {
    const r = clone(CLINVAR);
    delete r.germline_classification;
    return r;
  })()).ok);
check("a record with an empty variation_set is rejected",
  !V.validateClinVarRecord((() => {
    const r = clone(CLINVAR);
    r.variation_set = [];
    return r;
  })()).ok);
check("the empty-variation_set error names the field",
  /variation_set/.test(V.validateClinVarRecord((() => {
    const r = clone(CLINVAR);
    r.variation_set = [];
    return r;
  })()).error));
check("a record with no variation_set at all is rejected",
  !V.validateClinVarRecord((() => {
    const r = clone(CLINVAR);
    delete r.variation_set;
    return r;
  })()).ok);
check("a record with neither variation_name nor protein_change is rejected",
  !V.validateClinVarRecord((() => {
    const r = clone(CLINVAR);
    r.variation_set[0].variation_name = "";
    r.protein_change = "";
    return r;
  })()).ok);
check("a record with a non-string protein_change is rejected",
  !V.validateClinVarRecord((() => {
    const r = clone(CLINVAR);
    r.protein_change = 42;
    return r;
  })()).ok);
check("a record with a non-object variation_set[0] is rejected",
  !V.validateClinVarRecord((() => {
    const r = clone(CLINVAR);
    r.variation_set = ["nope"];
    return r;
  })()).ok);
check("a record with a non-string variation_loc is rejected",
  !V.validateClinVarRecord((() => {
    const r = clone(CLINVAR);
    r.variation_set[0].variation_loc = "chr7";
    return r;
  })()).ok);
check("a rejected record carries the invalid_shape code",
  V.validateClinVarRecord({}).code === V.CODES.INVALID_SHAPE);

// ----------------------------------------------------------------- gnomAD
console.log("\ngnomAD responses");
const GNOMAD = {
  data: {
    variant: {
      variant_id: "7-55249072-C-T",
      exome: { ac: 4, an: 1461884, af: 2.7e-06 },
      genome: { ac: 12, an: 12345678, af: 9.7e-07 },
    },
  },
};
check("a complete response is valid", V.validateGnomadResponse(GNOMAD).ok);
check("a genome-only response is valid",
  V.validateGnomadResponse({ data: { variant: { variant_id: "1-1-A-T", genome: { af: 0.01 } } } }).ok);
check("an exome-only response is valid",
  V.validateGnomadResponse({ data: { variant: { variant_id: "1-1-A-T", exome: { af: 0.01 } } } }).ok);

check("a null payload is rejected", !V.validateGnomadResponse(null).ok);
check("an array payload is rejected", !V.validateGnomadResponse([]).ok);
check("an empty object is rejected", !V.validateGnomadResponse({}).ok);
check("the empty-object error names the data field",
  /data/.test(V.validateGnomadResponse({}).error));

const gnomadNotFound = { errors: [{ message: "Variant not found" }], data: { variant: null } };
const notFound = V.validateGnomadResponse(gnomadNotFound);
check("an error with null data is rejected", !notFound.ok);
check("a missing variant is the not_found code", notFound.code === V.CODES.NOT_FOUND);
check("the not-found error repeats the API's own message",
  /Variant not found/.test(notFound.error));

const gnomadRateLimited = { errors: [{ message: "Query rate limit exceeded" }] };
const rateLimited = V.validateGnomadResponse(gnomadRateLimited);
check("a rate-limit response is rejected", !rateLimited.ok);
check("a rate limit is its own code", rateLimited.code === V.CODES.RATE_LIMITED);
check("rate limiting and absence are told apart",
  rateLimited.code !== notFound.code);
check("the rate-limit error says so", /rate limit/i.test(rateLimited.error));

check("an unexpected gnomAD error is a shape failure, not an absence",
  V.validateGnomadResponse({ errors: [{ message: "Internal server error" }] }).code
    === V.CODES.INVALID_SHAPE);
check("a response with no variant is rejected", !V.validateGnomadResponse({ data: {} }).ok);
check("a variant with no id is rejected",
  !V.validateGnomadResponse({ data: { variant: { genome: { af: 0.1 } } } }).ok);
check("a variant with no allele frequency is rejected",
  !V.validateGnomadResponse({ data: { variant: { variant_id: "1-1-A-T" } } }).ok);
check("a non-numeric allele frequency is rejected",
  !V.validateGnomadResponse({
    data: { variant: { variant_id: "1-1-A-T", genome: { af: "0.1" } } },
  }).ok);
check("a non-object dataset is rejected",
  !V.validateGnomadResponse({
    data: { variant: { variant_id: "1-1-A-T", genome: "genome" } },
  }).ok);

// ------------------------------------------------------------------- RCSB
console.log("\nRCSB polymer_entity responses");
const RCSB = {
  rcsb_polymer_entity_container_identifiers: {
    reference_sequence_identifiers: [
      { database_name: "UniProt", database_accession: "P00533" },
      { database_name: "GenBank", database_accession: "NG_009458" },
    ],
    entity_poly: { pdbx_description: "Tyrosine kinase" },
  },
  rcsb_polymer_entity_align: [{
    reference_database_name: "UniProt",
    reference_database_accession: "P00533",
    aligned_regions: [{ ref_beg_seq_id: 696, length: 282 }],
  }],
  entity_poly: { pdbx_strand_id: "A,B", type: "polypeptide(L)" },
};
check("a complete entity is valid", V.validateRcsbEntity(RCSB).ok);
check("an entity with only reference identifiers is valid",
  V.validateRcsbEntity({
    rcsb_polymer_entity_container_identifiers: {
      reference_sequence_identifiers: [{ database_name: "UniProt", database_accession: "P00533" }],
    },
  }).ok);
check("an entity with only an alignment is valid",
  V.validateRcsbEntity({
    rsb_polymer_entity_align: undefined,
    rcsb_polymer_entity_align: [{
      reference_database_name: "UniProt",
      reference_database_accession: "P00533",
      aligned_regions: [{ ref_beg_seq_id: 696, length: 282 }],
    }],
  }).ok);
check("an entity with only entity_poly is valid",
  V.validateRcsbEntity({ entity_poly: { pdbx_strand_id: "A" } }).ok);

check("a null payload is rejected", !V.validateRcsbEntity(null).ok);
check("an array payload is rejected", !V.validateRcsbEntity([RCSB]).ok);
check("an empty object is rejected", !V.validateRcsbEntity({}).ok);
check("the empty-object error names the blocks it looked for",
  /rcsb_polymer_entity_container_identifiers/.test(V.validateRcsbEntity({}).error));
check("a non-object container identifiers is rejected",
  !V.validateRcsbEntity({ rcsb_polymer_entity_container_identifiers: "P00533" }).ok);
check("non-array reference_sequence_identifiers is rejected",
  !V.validateRcsbEntity({
    rcsb_polymer_entity_container_identifiers: { reference_sequence_identifiers: {} },
  }).ok);
check("a non-object reference entry is rejected",
  !V.validateRcsbEntity({
    rcsb_polymer_entity_container_identifiers: { reference_sequence_identifiers: ["P00533"] },
  }).ok);
check("a non-string database_accession is rejected",
  !V.validateRcsbEntity({
    rcsb_polymer_entity_container_identifiers: {
      reference_sequence_identifiers: [{ database_name: "UniProt", database_accession: 533 }],
    },
  }).ok);
check("a non-array align block is rejected",
  !V.validateRcsbEntity({ rcsb_polymer_entity_align: {} }).ok);
check("a non-array aligned_regions is rejected",
  !V.validateRcsbEntity({
    rcsb_polymer_entity_align: [{ aligned_regions: {} }],
  }).ok);
check("a non-numeric ref_beg_seq_id is rejected",
  !V.validateRcsbEntity({
    rcsb_polymer_entity_align: [{ aligned_regions: [{ ref_beg_seq_id: "696", length: 282 }] }],
  }).ok);
check("a non-numeric aligned length is rejected",
  !V.validateRcsbEntity({
    rcsb_polymer_entity_align: [{ aligned_regions: [{ ref_beg_seq_id: 696, length: "282" }] }],
  }).ok);
check("a non-string pdbx_strand_id is rejected",
  !V.validateRcsbEntity({ entity_poly: { pdbx_strand_id: ["A", "B"] } }).ok);

// ----------------------------------------------------------------- UniProt
console.log("\nUniProt responses");
const UNIPROT = {
  primaryAccession: "P00533",
  genes: [{ geneName: { value: "EGFR" } }],
  sequence: { length: 1210, value: "MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKLTQLGTFEDHFLSL" },
};
check("a complete entry is valid", V.validateUniProtResponse(UNIPROT).ok);
check("a null payload is rejected", !V.validateUniProtResponse(null).ok);
check("an array payload is rejected", !V.validateUniProtResponse([UNIPROT]).ok);
check("an empty object is rejected", !V.validateUniProtResponse({}).ok);
check("the empty-object error names the sequence",
  /sequence/.test(V.validateUniProtResponse({}).error));
check("a missing sequence.value is rejected",
  !V.validateUniProtResponse({ genes: UNIPROT.genes, sequence: {} }).ok);
check("an empty sequence.value is rejected",
  !V.validateUniProtResponse({ genes: UNIPROT.genes, sequence: { value: "" } }).ok);
check("a missing genes array is rejected",
  !V.validateUniProtResponse({ sequence: UNIPROT.sequence }).ok);
check("an empty genes array is rejected",
  !V.validateUniProtResponse({ genes: [], sequence: UNIPROT.sequence }).ok);
check("a missing geneName object is rejected",
  !V.validateUniProtResponse({ genes: [{}], sequence: UNIPROT.sequence }).ok);
check("a missing geneName.value is rejected",
  !V.validateUniProtResponse({ genes: [{ geneName: {} }], sequence: UNIPROT.sequence }).ok);

// ------------------------------------------------------------------- FASTA
console.log("\nAlphaFold FASTA responses");
check("a normal FASTA is valid",
  V.validateAlphaFoldFasta(">sp|P00533|EGFR_HUMAN\nMRPSGTAGAALLALLAALCPASR\nALEEKKVCQGT\n").ok);
check("a single-line FASTA is valid",
  V.validateAlphaFoldFasta(">sp|P00533|EGFR_HUMAN\nMRPSGTAGAALL\n").ok);
check("a null payload is rejected", !V.validateAlphaFoldFasta(null).ok);
check("an array payload is rejected", !V.validateAlphaFoldFasta([]).ok);
check("a number payload is rejected", !V.validateAlphaFoldFasta(42).ok);
check("an empty body is rejected", !V.validateAlphaFoldFasta("").ok);
check("a whitespace-only body is rejected", !V.validateAlphaFoldFasta("   \n  ").ok);
check("an HTML error page is rejected",
  !V.validateAlphaFoldFasta("<!DOCTYPE html><html><body>404</body></html>").ok);
check("the HTML error says so", /HTML/.test(V.validateAlphaFoldFasta("<html>404</html>").error));
check("text with no header line is rejected",
  !V.validateAlphaFoldFasta("MRPSGTAGAALLALLAALCPASR\n").ok);
check("a header with no sequence is rejected",
  !V.validateAlphaFoldFasta(">sp|P00533|EGFR_HUMAN\n").ok);
check("a sequence containing a digit is rejected",
  !V.validateAlphaFoldFasta(">sp|SEQ\nMRPS1\n").ok);

// ------------------------------------------------------------------ batches
console.log("\nBatch validation");
const batch = V.validateClinVarRecords([CLINVAR, clone(CLINVAR), clone(CLINVAR)]);
check("a clean batch is valid", batch.ok);
check("a clean batch counts its records", batch.count === 3);
check("a clean batch has no failures", batch.failures.length === 0);

const mixed = V.validateClinVarRecords([
  CLINVAR,
  { accession_version: "RCV000000002.1" },
  clone(CLINVAR),
]);
check("a batch with one bad record is not ok", !mixed.ok);
check("the failure index points at the bad record", mixed.failures[0].index === 1);
check("the failure names the bad record", mixed.failures[0].id === "RCV000000002.1");
check("the failure carries a code", mixed.failures[0].code === V.CODES.INVALID_SHAPE);
check("the failure carries a message", mixed.failures[0].error.length > 10);
check("every failure is collected, not just the first",
  V.validateClinVarRecords([{}, {}, {}]).failures.length === 3);
check("the batch summary counts the failures",
  /1 of 3 records failed/.test(mixed.error));

check("an empty batch is valid", V.validateClinVarRecords([]).ok);
check("a non-array batch is rejected", !V.validateClinVarRecords(null).ok);
check("a non-array batch says what it got",
  /Expected a list/.test(V.validateClinVarRecords(null).error));

const list = V.validateList([1, 2, 3], (n) => (n % 2 === 0 ? { ok: true } : { ok: false, code: "odd", error: `${n} is odd` }));
check("validateList accepts a custom validator", !list.ok);
check("validateList counts the items", list.count === 3);
check("validateList collects every failure", list.failures.length === 2);
check("validateList uses the validator's code", list.failures[0].code === "odd");
check("validateList uses the validator's message", list.failures[0].error === "1 is odd");
check("validateList labels a failure when asked",
  V.validateList(["x"], () => ({ ok: false, error: "no" }),
    (item) => item.toUpperCase()).failures[0].id === "X");
check("validateList rejects a non-list", !V.validateList("nope", () => ({ ok: true })).ok);
check("validateList rejects a missing validator", !V.validateList([], null).ok);

// ------------------------------------------------------- message discipline
console.log("\nMessage discipline");
const MESSAGES = [
  V.validatePdbId("1M17X").error,
  V.validateGeneSymbol("EGFR&foo=bar").error,
  V.validateUniprotAccession("P0053").error,
  V.validateClinVarRecord({}).error,
  V.validateGnomadResponse(gnomadRateLimited).error,
  V.validateRcsbEntity({}).error,
  V.validateUniProtResponse({}).error,
  V.validateAlphaFoldFasta("").error,
];
check("every message is a non-empty string",
  MESSAGES.every((m) => typeof m === "string" && m.length > 10));
check("every message ends in a full stop", MESSAGES.every((m) => m.endsWith(".")));
check("no message leaks an angle bracket",
  !MESSAGES.some((m) => /[<>]/.test(m)));
// The message quotes the offending value, so the quoted part itself has
// to be inert: the value is echoed with anything that could close an
// attribute or open a tag replaced.
const echoOf = (m) => {
  const quoted = m.match(/"([^"]*)"/);
  return quoted ? quoted[1] : "";
};
const ECHOED = [
  V.validatePdbId("\"><script>alert(1)</script>").error,
  V.validateGeneSymbol("EGFR&foo=bar").error,
  V.validateGeneSymbol("../../etc/passwd").error,
  V.validateUniprotAccession("javascript:alert(1)").error,
];
check("the quoted echo carries no markup or metacharacter",
  ECHOED.every((m) => !/[<>&"']/.test(echoOf(m))));
check("the raw hostile input is never repeated verbatim",
  ECHOED.every((m) => !m.includes("\"><script>")) &&
  !ECHOED.some((m) => m.includes("EGFR&foo=bar")) &&
  !ECHOED.some((m) => m.includes("../../etc/passwd")) &&
  !ECHOED.some((m) => m.includes("javascript:alert(1)")));

console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
