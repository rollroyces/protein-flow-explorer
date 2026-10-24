/**
 * Unit tests for the tissue expression module.
 *
 *   node test/expression.test.js
 *
 * Covers the logic behind the expression panel: GTEx tissue ids turned
 * into readable labels, median TPM turned into a log value and a colour,
 * a ranking of tissues, a specificity verdict and a grouping of the 54
 * tissues into compartments.
 *
 * The fixtures are real GTEx v8 median-expression values (EGFR, ALB,
 * ACTB, MYH1, INS), inlined so the test never touches the network.
 */
const fs = require("fs");
const path = require("path");
const vm = require("vm");

const ROOT = path.join(__dirname, "..");
vm.runInThisContext(fs.readFileSync(path.join(ROOT, "src", "expression.js"), "utf8"),
  { filename: "src/expression.js" });
const E = globalThis.PFE.expression;

let pass = 0, fail = 0;
function check(name, cond, detail = "") {
  if (cond) { pass++; console.log(`  PASS  ${name}${detail ? "  " + detail : ""}`); }
  else { fail++; console.log(`  FAIL  ${name}${detail ? "  " + detail : ""}`); }
}
const near = (a, b, tol) => Math.abs(a - b) <= tol;
const hex = (v) => "0x" + v.toString(16).padStart(6, "0");

// ---------------------------------------------------------------- fixtures
// Real GTEx v8 median TPM for EGFR (ENSG00000146648.17), all 54 tissues.
const EGFR_ROWS = [
  ["Adipose_Subcutaneous", 39.7322],
  ["Adipose_Visceral_Omentum", 25.8382],
  ["Adrenal_Gland", 14.6419],
  ["Artery_Aorta", 23.4493],
  ["Artery_Coronary", 20.1819],
  ["Artery_Tibial", 33.1978],
  ["Bladder", 23.5296],
  ["Brain_Amygdala", 4.94781],
  ["Brain_Anterior_cingulate_cortex_BA24", 5.58234],
  ["Brain_Caudate_basal_ganglia", 7.9245600000000005],
  ["Brain_Cerebellar_Hemisphere", 2.39642],
  ["Brain_Cerebellum", 4.72511],
  ["Brain_Cortex", 6.2543],
  ["Brain_Frontal_Cortex_BA9", 5.61144],
  ["Brain_Hippocampus", 2.9834],
  ["Brain_Hypothalamus", 3.28776],
  ["Brain_Nucleus_accumbens_basal_ganglia", 6.9520100000000005],
  ["Brain_Putamen_basal_ganglia", 5.38999],
  ["Brain_Spinal_cord_cervical_c-1", 2.0341],
  ["Brain_Substantia_nigra", 4.05847],
  ["Breast_Mammary_Tissue", 32.9422],
  ["Cells_Cultured_fibroblasts", 60.6276],
  ["Cells_EBV-transformed_lymphocytes", 0.148094],
  ["Cervix_Ectocervix", 30.9548],
  ["Cervix_Endocervix", 24.7174],
  ["Colon_Sigmoid", 18.9776],
  ["Colon_Transverse", 15.0498],
  ["Esophagus_Gastroesophageal_Junction", 34.6313],
  ["Esophagus_Mucosa", 29.9136],
  ["Esophagus_Muscularis", 31.8016],
  ["Fallopian_Tube", 28.1227],
  ["Heart_Atrial_Appendage", 4.9930699999999995],
  ["Heart_Left_Ventricle", 3.11803],
  ["Kidney_Cortex", 5.2046],
  ["Kidney_Medulla", 10.6746],
  ["Liver", 16.9142],
  ["Lung", 22.1273],
  ["Minor_Salivary_Gland", 31.625999999999998],
  ["Muscle_Skeletal", 6.7897300000000005],
  ["Nerve_Tibial", 43.1228],
  ["Ovary", 20.7359],
  ["Pancreas", 6.84974],
  ["Pituitary", 4.23149],
  ["Prostate", 27.0509],
  ["Skin_Not_Sun_Exposed_Suprapubic", 75.9298],
  ["Skin_Sun_Exposed_Lower_leg", 78.3374],
  ["Small_Intestine_Terminal_Ileum", 12.8558],
  ["Spleen", 6.2702],
  ["Stomach", 11.0075],
  ["Testis", 5.731380000000001],
  ["Thyroid", 25.0389],
  ["Uterus", 21.3998],
  ["Vagina", 40.6031],
  ["Whole_Blood", 0.0450544],
].map(([tissueSiteDetailId, median]) => ({
  tissueSiteDetailId, median, unit: "TPM", datasetId: "gtex_v8",
  gencodeId: "ENSG00000146648.17", geneSymbol: "EGFR",
}));

const rows = (pairs) => pairs.map(([tissueSiteDetailId, median]) => ({
  tissueSiteDetailId, median, unit: "TPM",
}));

// --------------------------------------------------------- module isolation
console.log("\nModule isolation");
check("no DOM required", typeof document === "undefined");
check("no window required", typeof window === "undefined");
check("PFE namespace present", !!E);
check("the namespace is frozen", Object.isFrozen(E));

const src = fs.readFileSync(path.join(ROOT, "src", "expression.js"), "utf8");
check("no network call", !/\bfetch\s*\(|XMLHttpRequest|\bWebSocket\b/.test(src));
check("no DOM use", !/\bdocument\s*\.|\bwindow\s*\./.test(src));
check("no timers", !/setTimeout\s*\(|setInterval\s*\(|requestAnimationFrame/.test(src));
check("no module syntax", !/^\s*(import|export)\s/m.test(src));

// ------------------------------------------------------------- thresholds
console.log("\nDocumented thresholds");
check("detection threshold is 1 TPM", E.DETECTION_TPM === 1);
check("elevation fold is 4x", E.ENHANCED_FOLD === 4);
check("dominance ratio is 10x", E.ENRICHED_RATIO === 10);
check("the ramp tops out at 1000 TPM", E.EXPRESSION_TPM_MAX === 1000);
check("four ramp stops", E.EXPRESSION_STOPS.length === 4, `(${E.EXPRESSION_STOPS.length})`);
check("the stops are the documented colours",
      E.EXPRESSION_STOPS.join(",") === "6621481,0,168,120,245,197,66,204,31,31".split(",").map(Number).join(",") ||
      E.EXPRESSION_STOPS[0] === 0x64748b && E.EXPRESSION_STOPS[1] === 0x00a878 &&
      E.EXPRESSION_STOPS[2] === 0xf5c542 && E.EXPRESSION_STOPS[3] === 0xcc1f1f,
      `(${E.EXPRESSION_STOPS.map(hex).join(" ")})`);
check("the stops are distinct", new Set(E.EXPRESSION_STOPS).size === E.EXPRESSION_STOPS.length);
check("every stop is a 24-bit colour",
      E.EXPRESSION_STOPS.every((c) => Number.isInteger(c) && c >= 0 && c <= 0xffffff));
check("the ramp runs cool to hot",
      (E.EXPRESSION_STOPS[0] >> 16) < (E.EXPRESSION_STOPS[3] >> 16) &&
      (E.EXPRESSION_STOPS[2] >> 8) > (E.EXPRESSION_STOPS[0] >> 8));

// ------------------------------------------------------------ tissue labels
console.log("\nTissue labels");
check("sun-exposed lower leg",
      E.tissueLabel("Skin_Sun_Exposed_Lower_leg") === "Skin (sun-exposed lower leg)");
check("cerebellar hemisphere",
      E.tissueLabel("Brain_Cerebellar_Hemisphere") === "Brain — cerebellar hemisphere");
check("cultured fibroblasts",
      E.tissueLabel("Cells_Cultured_fibroblasts") === "Cultured fibroblasts");
check("whole blood", E.tissueLabel("Whole_Blood") === "Whole blood");
check("a brain region",
      E.tissueLabel("Brain_Substantia_nigra") === "Brain — substantia nigra");
check("the c-1 tissue id is mapped",
      E.tissueLabel("Brain_Spinal_cord_cervical_c-1") === "Brain — spinal cord (cervical C1)");
check("an unmapped id falls back to title case",
      E.tissueLabel("Some_New_Tissue") === "Some New Tissue");
check("the fallback keeps hyphens",
      E.tissueLabel("Foo_bar-baz") === "Foo Bar-baz");
check("an empty id has no label", E.tissueLabel("") === "");
check("a missing id has no label", E.tissueLabel(undefined) === "");
check("a non-string id has no label", E.tissueLabel(42) === "");

// ------------------------------------------------------------ compartments
console.log("\nCompartments");
check("skin", E.tissueCompartment("Skin_Sun_Exposed_Lower_leg") === "skin");
check("brain", E.tissueCompartment("Brain_Cerebellum") === "brain");
check("muscle", E.tissueCompartment("Muscle_Skeletal") === "muscle");
check("cardiovascular", E.tissueCompartment("Heart_Left_Ventricle") === "cardiovascular");
check("digestive", E.tissueCompartment("Liver") === "digestive");
check("reproductive", E.tissueCompartment("Uterus") === "reproductive");
check("immune", E.tissueCompartment("Whole_Blood") === "immune");
check("endocrine", E.tissueCompartment("Thyroid") === "endocrine");
check("an unmapped tissue falls back to other",
      E.tissueCompartment("Some_New_Tissue") === "other");
check("a missing tissue falls back to other", E.tissueCompartment(undefined) === "other");
check("the compartment label", E.compartmentLabel("immune") === "Immune & blood");
check("an unknown compartment keeps its id",
      E.compartmentLabel("nonsense") === "nonsense");
check("the order lists every compartment in the table",
      E.COMPARTMENT_ORDER.length === 14 &&
      Object.values(E.TISSUE_COMPARTMENTS).every((c) => E.COMPARTMENT_ORDER.includes(c)),
      `(${E.COMPARTMENT_ORDER.length})`);
check("every compartment in the order has a label",
      E.COMPARTMENT_ORDER.every((c) => !!E.COMPARTMENT_LABELS[c]));
// Every tissue GTEx returns must be classified, or the UI would drop it.
const ALL_TISSUES = EGFR_ROWS.map((r) => r.tissueSiteDetailId);
check("all 54 GTEx tissues are labelled",
      ALL_TISSUES.every((t) => E.TISSUE_LABELS[t] && E.TISSUE_LABELS[t].length > 2),
      `(${ALL_TISSUES.length})`);
check("all 54 GTEx tissues have a compartment",
      ALL_TISSUES.every((t) => E.TISSUE_COMPARTMENTS[t] &&
        E.COMPARTMENT_ORDER.includes(E.TISSUE_COMPARTMENTS[t])),
      `(${ALL_TISSUES.length})`);

// ------------------------------------------------------------------- the log
console.log("\nLog normalisation");
check("zero stays zero", E.logExpression(0) === 0);
check("9 TPM -> log10(10) = 1", near(E.logExpression(9), 1, 1e-12));
check("99 TPM -> log10(100) = 2", near(E.logExpression(99), 2, 1e-12));
check("999 TPM -> log10(1000) = 3", near(E.logExpression(999), 3, 1e-12));
check("a negative value is zero", E.logExpression(-5) === 0);
check("NaN is zero", E.logExpression(NaN) === 0);
check("Infinity is zero", E.logExpression(Infinity) === 0);
check("a missing value is zero", E.logExpression(undefined) === 0);
check("null is zero", E.logExpression(null) === 0);
check("a non-number is zero", E.logExpression("78") === 0);
// The whole point of the log: 0.045 and 78 have to be far apart.
check("0.045 TPM and 78 TPM are far apart on the log scale",
      E.logExpression(78.3374) - E.logExpression(0.0450544) > 1.8,
      `(${E.logExpression(78.3374).toFixed(3)} vs ${E.logExpression(0.0450544).toFixed(3)})`);
check("the log is monotone",
      E.logExpression(0.0450544) < E.logExpression(1) &&
      E.logExpression(1) < E.logExpression(78.3374));

// ------------------------------------------------------------------ colours
console.log("\nColour ramp");
check("no expression is the first stop",
      E.expressionColor(0) === E.EXPRESSION_STOPS[0], `(${hex(E.expressionColor(0))})`);
check("a missing value is the first stop",
      E.expressionColor(undefined) === E.EXPRESSION_STOPS[0]);
check("a negative value is the first stop",
      E.expressionColor(-1) === E.EXPRESSION_STOPS[0]);
check("NaN is the first stop", E.expressionColor(NaN) === E.EXPRESSION_STOPS[0]);
check("the top of the scale is the last stop",
      E.expressionColor(E.EXPRESSION_TPM_MAX) === E.EXPRESSION_STOPS[3],
      `(${hex(E.expressionColor(E.EXPRESSION_TPM_MAX))})`);
check("a value above the ramp maximum clamps to the last stop",
      E.expressionColor(50000) === E.EXPRESSION_STOPS[3]);
check("the colour of a mid value is between the stops",
      E.expressionColor(10) !== E.EXPRESSION_STOPS[0] &&
      E.expressionColor(10) !== E.EXPRESSION_STOPS[3]);
check("a higher value never goes back down the ramp",
      E.expressionColor(1) !== E.expressionColor(10) &&
      E.expressionColor(10) !== E.expressionColor(100));
check("every colour is a 24-bit integer",
      [0, 1, 10, 78.3374, 1000, 5000].every((v) => {
        const c = E.expressionColor(v);
        return Number.isInteger(c) && c >= 0 && c <= 0xffffff;
      }));
check("the colours are legible on the dark UI",
      E.EXPRESSION_STOPS.every((c) => {
        const lum = 0.2126 * ((c >> 16) & 255) + 0.7152 * ((c >> 8) & 255) +
                    0.0722 * (c & 255);
        return lum > 60;
      }));
// A relative scale: the ceiling is the caller's choice.
check("a relative scale puts the ceiling at the last stop",
      E.expressionColor(100, 100) === E.EXPRESSION_STOPS[3]);
check("a relative scale keeps a lower value below it",
      E.expressionColor(50, 100) !== E.EXPRESSION_STOPS[3]);
check("a relative scale still starts at the first stop",
      E.expressionColor(0, 100) === E.EXPRESSION_STOPS[0]);
check("a broken ceiling falls back to the absolute scale",
      E.expressionColor(E.EXPRESSION_TPM_MAX, NaN) === E.EXPRESSION_STOPS[3]);

// --------------------------------------------------------- parsing payloads
console.log("\nParsing a GTEx payload");
const parsed = E.parseMedianExpression({ data: EGFR_ROWS });
check("every row survives", parsed.length === EGFR_ROWS.length, `(${parsed.length})`);
check("the median is preserved", parsed[0].median === 39.7322);
check("a bare array is accepted",
      E.parseMedianExpression(EGFR_ROWS).length === EGFR_ROWS.length);
check("a missing median becomes zero",
      E.parseMedianExpression([{ tissueSiteDetailId: "Liver" }])[0].median === 0);
check("a NaN median becomes zero",
      E.parseMedianExpression([{ tissueSiteDetailId: "Liver", median: NaN }])[0].median === 0);
check("a row without a tissue id is dropped",
      E.parseMedianExpression([{ median: 5 }]).length === 0);
check("a null payload yields nothing", E.parseMedianExpression(null).length === 0);
check("a payload without data yields nothing", E.parseMedianExpression({}).length === 0);

// ------------------------------------------------------------------ ranking
console.log("\nRanking");
const ranked = E.rankTissues(EGFR_ROWS);
check("every tissue is ranked", ranked.length === 54, `(${ranked.length})`);
check("the order is descending",
      ranked.every((r, i) => i === 0 || ranked[i - 1].tpm >= r.tpm));
check("the top tissue is sun-exposed skin",
      ranked[0].tissueId === "Skin_Sun_Exposed_Lower_leg" && ranked[0].tpm === 78.3374);
check("the bottom tissue is whole blood",
      ranked[53].tissueId === "Whole_Blood" && ranked[53].tpm === 0.0450544);
check("each entry carries a label", ranked[0].label === "Skin (sun-exposed lower leg)");
check("each entry carries a compartment", ranked[0].compartment === "skin");
check("each entry carries a log value",
      near(ranked[0].log, E.logExpression(78.3374), 1e-12));
check("each entry carries a colour", ranked[0].color === E.expressionColor(78.3374));
check("each entry carries the detection flag",
      ranked[0].detected === true && ranked[53].detected === false);
check("an empty array ranks nothing", E.rankTissues([]).length === 0);
check("a missing array ranks nothing", E.rankTissues(undefined).length === 0);
check("a single tissue ranks one entry",
      E.rankTissues(rows([["Liver", 12]])).length === 1);
check("an unmapped tissue still gets a label",
      E.rankTissues(rows([["Some_New_Tissue", 5]]))[0].label === "Some New Tissue");
check("an unmapped tissue falls into other",
      E.rankTissues(rows([["Some_New_Tissue", 5]]))[0].compartment === "other");
check("a missing median is zero",
      E.rankTissues([{ tissueSiteDetailId: "Liver" }])[0].tpm === 0);
check("a missing median is not detected",
      E.rankTissues([{ tissueSiteDetailId: "Liver" }])[0].detected === false);
check("a NaN median is zero",
      E.rankTissues([{ tissueSiteDetailId: "Liver", median: NaN }])[0].tpm === 0);
check("a negative median is zero",
      E.rankTissues([{ tissueSiteDetailId: "Liver", median: -3 }])[0].tpm === 0);
check("rows without a tissue id are dropped",
      E.rankTissues([{ median: 5 }, { tissueSiteDetailId: "Liver", median: 5 }]).length === 1);

// ------------------------------------------------------------------ grouping
console.log("\nGrouping");
const groups = E.groupByCompartment(EGFR_ROWS);
check("the groups are in the declared order",
      groups.every((g, i) => E.COMPARTMENT_ORDER.indexOf(g.compartment) >=
        (i === 0 ? -1 : E.COMPARTMENT_ORDER.indexOf(groups[i - 1].compartment))));
check("empty compartments are dropped",
      groups.every((g) => g.tissues.length > 0));
check("every tissue is in exactly one group",
      groups.reduce((n, g) => n + g.tissues.length, 0) === 54,
      `(${groups.reduce((n, g) => n + g.tissues.length, 0)})`);
check("the tissues inside a group keep the ranking",
      groups.every((g) => g.tissues.every((t, i) => i === 0 || g.tissues[i - 1].tpm >= t.tpm)));
check("each group carries a label",
      groups.every((g) => typeof g.label === "string" && g.label.length > 0));
check("skin is a group of two",
      groups.find((g) => g.compartment === "skin").tissues.length === 2);
check("the brain group holds every brain tissue",
      groups.find((g) => g.compartment === "brain").tissues.length === 13,
      `(${groups.find((g) => g.compartment === "brain").tissues.length})`);
check("an unmapped tissue lands in other",
      E.groupByCompartment(rows([["Some_New_Tissue", 5]]))[0].compartment === "other");
check("an empty payload groups nothing", E.groupByCompartment([]).length === 0);

// -------------------------------------------------------------- specificity
console.log("\nSpecificity verdict");
const egfr = E.summarizeExpression(EGFR_ROWS);
check("EGFR is tissue enhanced", egfr.specificity === "tissue_enhanced",
      `(${egfr.specificity})`);
check("the verdict has a label", egfr.specificityLabel === "Tissue enhanced");
check("EGFR is detected in 52 of 54 tissues", egfr.detectedCount === 52,
      `(${egfr.detectedCount}/${egfr.tissueCount})`);
check("the summary counts every tissue", egfr.tissueCount === 54);
check("the leading tissue is reported",
      egfr.maxTissueId === "Skin_Sun_Exposed_Lower_leg" && egfr.maxTpm === 78.3374);
check("the leading tissue has a label",
      egfr.maxTissueLabel === "Skin (sun-exposed lower leg)");
check("the cross-tissue median is reported", near(egfr.medianTpm, 14.8458, 0.001),
      `(${egfr.medianTpm})`);
check("the fold over the median is reported", near(egfr.foldVsMedian, 78.3374 / 14.8458, 0.001));
check("the elevated tissues are listed", egfr.elevatedCount === 2 &&
      egfr.elevatedTissueIds.includes("Skin_Sun_Exposed_Lower_leg"));
check("the two skin samples are not dominant enough to be enriched",
      near(egfr.dominanceRatio, 78.3374 / 75.9298, 1e-9));

console.log("\nSpecificity thresholds");
// The detection threshold is inclusive: 1 TPM counts as detected.
check("exactly 1 TPM is detected",
      E.summarizeExpression(rows([["Liver", 1]])).detectedCount === 1);
check("just under 1 TPM is not detected",
      E.summarizeExpression(rows([["Liver", 0.999]])).detectedCount === 0);

// Elevation is 4x the baseline, inclusive. Here the baseline is 10
// (both the median and the mean of the other tissues).
check("4x the baseline is elevated",
      E.summarizeExpression(rows([["A", 40], ["B", 10], ["C", 10]])).elevatedCount === 1);
check("just under 4x is not elevated",
      E.summarizeExpression(rows([["A", 39.9], ["B", 10], ["C", 10]])).elevatedCount === 0);

// Dominance is 10x the runner-up, inclusive.
check("10x the runner-up is enriched",
      E.summarizeExpression(rows([["A", 100], ["B", 10], ["C", 10]])).specificity ===
        "tissue_enriched");
check("just under 10x the runner-up is only enhanced",
      E.summarizeExpression(rows([["A", 99], ["B", 10], ["C", 10]])).specificity ===
        "tissue_enhanced");
check("one tissue far above the rest is enriched",
      E.summarizeExpression(rows([["Pancreas", 2325.26], ["Liver", 0.4]])).specificity ===
        "tissue_enriched");
check("a flat gene has low tissue specificity",
      E.summarizeExpression(rows([["A", 2600], ["B", 2500], ["C", 2400], ["D", 2300]]))
        .specificity === "low_tissue_specificity");
check("a gene below the detection threshold is not detected",
      E.summarizeExpression(rows([["A", 0.9], ["B", 0.5]])).specificity === "not_detected");
check("a gene with no expression at all is not detected",
      E.summarizeExpression(rows([["A", 0], ["B", 0]])).specificity === "not_detected");
check("an empty payload has insufficient data",
      E.summarizeExpression([]).specificity === "insufficient_data");
check("a single tissue has insufficient data",
      E.summarizeExpression(rows([["Liver", 50]])).specificity === "insufficient_data");
check("an empty payload counts no tissues",
      E.summarizeExpression([]).tissueCount === 0 &&
      E.summarizeExpression([]).detectedCount === 0);
check("a single tissue is still counted",
      E.summarizeExpression(rows([["Liver", 50]])).tissueCount === 1);
check("every verdict is a declared class",
      Object.keys(E.SPECIFICITY_LABELS).length === E.SPECIFICITY_CLASSES.length &&
      E.SPECIFICITY_CLASSES.every((c) => E.SPECIFICITY_LABELS[c]),
      `(${E.SPECIFICITY_CLASSES.length} classes)`);
check("the summary's verdict is a declared class",
      E.SPECIFICITY_CLASSES.includes(egfr.specificity));
check("the class list is ordered most specific first",
      E.SPECIFICITY_CLASSES[0] === "tissue_enriched" &&
      E.SPECIFICITY_CLASSES[1] === "tissue_enhanced");

console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
