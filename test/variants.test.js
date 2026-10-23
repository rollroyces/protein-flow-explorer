/**
 * Unit tests for the variant module.
 *
 *   node test/variants.test.js
 *
 * Covers the logic that turns clinical variant records into something you
 * can put on a structure: parsing ClinVar's protein_change, deriving a gnomAD
 * variant id, classifying significance, and mapping canonical UniProt
 * positions onto PDB residue numbers.
 */
const fs = require("fs");
const path = require("path");
const vm = require("vm");

const ROOT = path.join(__dirname, "..");
vm.runInThisContext(fs.readFileSync(path.join(ROOT, "src", "variants.js"), "utf8"),
  { filename: "src/variants.js" });
const V = globalThis.PFE.variants;

let pass = 0, fail = 0;
function check(name, cond, detail = "") {
  if (cond) { pass++; console.log(`  PASS  ${name}${detail ? "  " + detail : ""}`); }
  else { fail++; console.log(`  FAIL  ${name}${detail ? "  " + detail : ""}`); }
}

// --------------------------------------------------------- protein_change
console.log("\nProtein-change parsing");
// Real ClinVar value: "N1095S, N881S, N1103S, N1148S"
const multi = V.parseProteinChange("N1095S, N881S, N1103S, N1148S");
check("four changes parsed", multi.length === 4, `(${multi.length})`);
check("first change is N1095S", multi[0].wt === "N" && multi[0].pos === 1095 && multi[0].mut === "S");
check("second change is N881S", multi[1].wt === "N" && multi[1].pos === 881 && multi[1].mut === "S");
const single = V.parseProteinChange("N1148S");
check("a single change parses", single.length === 1 && single[0].pos === 1148);
check("whitespace is tolerated", V.parseProteinChange(" N1148S , T790M ").length === 2);
check("a leading p. is tolerated", V.parseProteinChange("p.N1148S").length === 1);
check("non-missense changes are skipped", V.parseProteinChange("").length === 0);
check("nonsense input yields nothing", V.parseProteinChange(null).length === 0);
check("stop codons are marked", V.parseProteinChange("R175*")[0].mut === "*");
check("a non-numeric position is rejected", V.parseProteinChange("NxxS").length === 0);
check("position zero is rejected", V.parseProteinChange("N0S").length === 0);
check("three-letter wild type is accepted",
      V.parseProteinChange("Asn1148Ser")[0].wt === "N" &&
      V.parseProteinChange("Asn1148Ser")[0].pos === 1148 &&
      V.parseProteinChange("Asn1148Ser")[0].mut === "S");
check("a stop codon is a star", V.parseProteinChange("Trp26Ter")[0].mut === "*");
check("a frameshift keeps its label", V.parseProteinChange("Leu1034fs")[0].mut === "fs");
check("the protein change can come from a variation name",
      V.extractProteinChange("NM_005228.5(EGFR):c.3443A>G (p.Asn1148Ser)").length === 1);
check("the change from a variation name is the canonical numbering",
      V.extractProteinChange("NM_005228.5(EGFR):c.2369C>T (p.Thr790Met)")[0].pos === 790);
check("a variation name with no protein change yields nothing",
      V.extractProteinChange("NM_005228.5(EGFR):c.2369C>T").length === 0);

// ------------------------------------------------------------------- SPDI
console.log("\nSPDI parsing and gnomAD ids");
// Real ClinVar value for EGFR c.3443A>G (p.Asn1148Ser)
const spdi = V.parseSpdi("NC_000007.14:55205426:A:G");
check("chromosome extracted", spdi.chrom === "7", `("${spdi.chrom}")`);
check("position is 1-based", spdi.pos === 55205427, `(${spdi.pos})`);
check("reference allele", spdi.ref === "A");
check("alternate allele", spdi.alt === "G");
check("a malformed SPDI is rejected", V.parseSpdi("not a spdi") === null);
check("an empty SPDI is rejected", V.parseSpdi("") === null);
check("a missing allele is rejected", V.parseSpdi("NC_000007.14:55205426::G") === null);

check("gnomAD id from SPDI",
      V.gnomadVariantId({ canonical_spdi: "NC_000007.14:55205426:A:G" }) === "7-55205427-A-G",
      `(${V.gnomadVariantId({ canonical_spdi: "NC_000007.14:55205426:A:G" })})`);
check("gnomAD id needs the alleles", V.gnomadVariantId({ canonical_spdi: "" }) === null);
check("gnomAD id is null without a SPDI", V.gnomadVariantId({}) === null);
check("a real ClinVar variation_set yields an id",
      V.gnomadVariantId({
        canonical_spdi: "NC_000007.14:55205426:A:G",
        variation_loc: [{ assembly_name: "GRCh38", chr: "7", start: "55205427" }],
      }) === "7-55205427-A-G");

// ------------------------------------------------------------ significance
console.log("\nSignificance classification");
check("Pathogenic", V.classifySignificance("Pathogenic") === "pathogenic");
check("Likely pathogenic", V.classifySignificance("Likely pathogenic") === "likely_pathogenic");
check("Uncertain significance", V.classifySignificance("Uncertain significance") === "vus");
check("Likely benign", V.classifySignificance("Likely benign") === "likely_benign");
check("Benign", V.classifySignificance("Benign") === "benign");
check("Benign/Likely benign", V.classifySignificance("Benign/Likely benign") === "likely_benign");
check("Pathogenic/Likely pathogenic",
      V.classifySignificance("Pathogenic/Likely pathogenic") === "likely_pathogenic");
check("conflicting beats pathogenic",
      V.classifySignificance("Conflicting interpretations of pathogenicity") === "conflicting");
// ClinVar's germline classes include drug response and risk factor, which
// matter here: EGFR T790M and L858R are classified "drug response", not
// pathogenic, because they are resistance mutations rather than inherited
// disease alleles.
check("drug response is its own class",
      V.classifySignificance("drug response") === "drug_response");
check("risk factor is its own class",
      V.classifySignificance("risk factor") === "risk_factor");
check("protective is its own class",
      V.classifySignificance("protective") === "protective");
check("drug response outranks uncertain",
      V.significanceRank("drug_response") < V.significanceRank("vus"));
check("every class has a colour",
      V.SIGNIFICANCE_ORDER.every((c) => typeof V.SIGNIFICANCE_COLORS[c] === "number"));
check("the classes are all distinct",
      new Set(V.SIGNIFICANCE_ORDER).size === V.SIGNIFICANCE_ORDER.length);

console.log("\nClinical-impact tier and oncogenicity");
check("Tier I - Strong", V.classifyImpactTier("Tier I - Strong") === "Tier I");
check("Tier II", V.classifyImpactTier("Tier II") === "Tier II");
check("Tier III - Moderate", V.classifyImpactTier("Tier III - Moderate") === "Tier III");
check("an unknown tier is null", V.classifyImpactTier("nonsense") === null);
check("a missing tier is null", V.classifyImpactTier("") === null);
check("Oncogenic", V.classifyOncogenicity("Oncogenic") === "Oncogenic");
check("Likely oncogenic", V.classifyOncogenicity("Likely oncogenic") === "Likely oncogenic");
check("an unknown oncogenicity is null", V.classifyOncogenicity("nonsense") === null);
check("an empty description is unknown", V.classifySignificance("") === "unknown");
check("a missing description is unknown", V.classifySignificance(undefined) === "unknown");
check("case is ignored", V.classifySignificance("pathogenic") === "pathogenic");

console.log("\nSignificance ordering");
check("pathogenic is the most severe", V.significanceRank("pathogenic") === 0);
check("benign is the least severe", V.significanceRank("benign") > V.significanceRank("vus"));
check("unknown ranks last", V.significanceRank("unknown") >= V.significanceRank("benign"));
check("every class has a colour",
      V.SIGNIFICANCE_ORDER.every((s) => typeof V.SIGNIFICANCE_COLORS[s] === "number"));

console.log("\nFrequency bands");
check("0 is rare", V.rarityLabel(0) === "rare");
check("1e-6 is rare", V.rarityLabel(0.000001) === "rare");
check("1e-3 is uncommon", V.rarityLabel(0.001) === "uncommon");
check("0.02 is common", V.rarityLabel(0.02) === "common");
check("0.5 is very common", V.rarityLabel(0.5) === "very common");
check("a missing frequency is unknown", V.rarityLabel(undefined) === "unknown");

// --------------------------------------------------------------- the map
console.log("\nUniProt -> PDB residue map");
// A UniProt segment and the chain it covers, aligned with one gap.
const upSeg = "MKWVTFISLL";
const chainSeq = "MKWVTFISLL";
const residues = chainSeq.split("").map((_, i) => ({ resi: 100 + i, resn: "ALA" }));
const aligned = { a: upSeg, b: chainSeq };
const map = V.buildUniprotMap(aligned.a, aligned.b, residues, 1, 1000);
check("every position maps", map.size === upSeg.length, `(${map.size})`);
check("position 1 maps to the first residue", map.get(1).resi === 100);
check("position 10 maps to the tenth residue", map.get(10).resi === 109);
check("the map is empty when the ranges do not overlap",
      V.buildUniprotMap(aligned.a, aligned.b, residues, 5000, 6000).size === 0);
check("mismatched residues are not mapped",
      V.buildUniprotMap("AAAA", "CCCC", residues, 1, 1000).size === 0);
// A gap in the UniProt sequence shifts the chain index but not the UniProt
// index, which is the whole point of the map.
const gapped = V.buildUniprotMap("A-AAAAA", "AAAAAAA", residues, 1, 1000);
check("gapped columns are skipped, the rest still map", gapped.size === 6, `(${gapped.size})`);
check("position 1 maps to the first residue", gapped.get(1).resi === 100);
check("after a gap the chain index has moved on", gapped.get(2).resi === 102,
      `(${gapped.get(2) ? gapped.get(2).resi : "none"})`);

// ------------------------------------------------------------- merging
console.log("\nMerging duplicate changes");
const merged = V.mergeByChange([
  { pos: 790, wt: "T", mut: "M", significance: "vus", af: 0.01 },
  { pos: 790, wt: "T", mut: "M", significance: "pathogenic", af: 0.02 },
  { pos: 858, wt: "L", mut: "R", significance: "likely_pathogenic", af: 0.001 },
]);
check("duplicates collapse", merged.length === 2, `(${merged.length})`);
const t790 = merged.find((m) => m.pos === 790);
check("the most severe significance wins", t790.significance === "pathogenic");
check("the lowest frequency wins", t790.af === 0.01);
check("a different change is kept separate", merged.some((m) => m.pos === 858));
check("records without a position are dropped",
      V.mergeByChange([{ wt: "T", mut: "M" }]).length === 0);
check("an empty list is handled", V.mergeByChange([]).length === 0);

console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
