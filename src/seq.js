/* ===================================================================
 * Protein Flow Explorer — sequence & alignment module
 *
 * Pure sequence-domain logic: PDB sequence parsing, FASTA parsing,
 * Needleman-Wunsch alignment, progressive MSA, conservation scoring and
 * residue indexing. No DOM, no network, no app state — unit-testable in isolation.
 *
 * Extracted from app.js during the 2026 architecture audit. The module is
 * loaded as a classic <script> before app.js and publishes itself on the
 * single PFE namespace; app.js destructures the members it needs, so no
 * call site changed when the code moved.
 * =================================================================== */
(function (global) {
  "use strict";

  const AA3TO1 = {
    ALA: "A", ARG: "R", ASN: "N", ASP: "D", CYS: "C",
    GLN: "Q", GLU: "E", GLY: "G", HIS: "H", ILE: "I",
    LEU: "L", LYS: "K", MET: "M", PHE: "F", PRO: "P",
    SER: "S", THR: "T", TRP: "W", TYR: "Y", VAL: "V",
    MSE: "M", SEC: "U", PYR: "U", // rare / non-standard
  };

  const AA1TO3 = {};
  for (const k in AA3TO1) AA1TO3[AA3TO1[k]] = k;

  const NW_MATCH = 2;

  const NW_MISMATCH = -1;

  const NW_GAP = -2;

  function parseSequenceFromPDB(text, format) {
    // Build per-chain residue list from ATOM records + SS from HELIX/SHEET.
    const chainsMap = {}; // chain -> [residues]
    const ssByResi = {}; // chain -> resi -> 'H' | 'E' | ''

    const lines = text.split(/\r?\n/);
    const seen = {}; // key chain+resi+atom -> true (dedupe altloc)
    for (const line of lines) {
      if (!line.startsWith("ATOM") && !line.startsWith("HETATM")) continue;
      const chain = line.substring(21, 22).trim() || "A";
      // For CA-only sequence display, keep only CA atoms
      const atomName = line.substring(12, 16).trim();
      if (atomName !== "CA" && atomName !== "P" && format === "pdb") continue;
      const resi = parseInt(line.substring(22, 26).trim(), 10);
      const resn = line.substring(17, 20).trim();
      const ins = line.substring(26, 27).trim();
      if (!resi) continue;
      const key = chain + ":" + resi + ":" + ins;
      if (seen[key]) continue;
      seen[key] = true;
      if (!chainsMap[chain]) chainsMap[chain] = [];
      chainsMap[chain].push({ resi, resn, insCode: ins });
    }
    // SS
    for (const line of lines) {
      if (line.startsWith("HELIX")) {
        const chain = line.substring(19, 20).trim() || "A";
        const start = parseInt(line.substring(21, 25).trim(), 10);
        const end = parseInt(line.substring(33, 37).trim(), 10);
        if (!ssByResi[chain]) ssByResi[chain] = {};
        for (let i = start; i <= end; i++) ssByResi[chain][i] = "H";
      } else if (line.startsWith("SHEET")) {
        const chain = line.substring(21, 22).trim() || "A";
        const start = parseInt(line.substring(22, 26).trim(), 10);
        const end = parseInt(line.substring(33, 37).trim(), 10);
        if (!ssByResi[chain]) ssByResi[chain] = {};
        for (let i = start; i <= end; i++) ssByResi[chain][i] = "E";
      }
    }
    // Build ordered list
    const chainIds = Object.keys(chainsMap).sort();
    return {
      chains: chainIds.map((c) => ({
        id: c,
        residues: chainsMap[c].map((r) => ({
          resi: r.resi,
          resn: r.resn,
          one: AA3TO1[r.resn] || "X",
          ss: ssByResi[c]?.[r.resi] || "",
        })),
      })),
    };
  }

  function parseFastaSequences(text) {
    const out = [];
    const lines = text.split(/\r?\n/);
    let cur = "";
    for (const line of lines) {
      if (line.startsWith(">")) {
        if (cur) out.push(cur);
        cur = "";
      } else {
        cur += line.trim();
      }
    }
    if (cur) out.push(cur);
    return out;
  }

  function nwAlign(a, b) {
    const m = a.length, n = b.length;
    // Score matrix
    const F = new Int16Array((m + 1) * (n + 1));
    // Traceback matrix: 0=diag, 1=up, 2=left
    const Tb = new Uint8Array((m + 1) * (n + 1));
    for (let i = 0; i <= m; i++) F[i * (n + 1)] = i * NW_GAP;
    for (let j = 0; j <= n; j++) F[j] = j * NW_GAP;
    for (let i = 1; i <= m; i++) {
      for (let j = 1; j <= n; j++) {
        const s = a[i - 1] === b[j - 1] ? NW_MATCH : NW_MISMATCH;
        const diag = F[(i - 1) * (n + 1) + (j - 1)] + s;
        const up = F[(i - 1) * (n + 1) + j] + NW_GAP;
        const left = F[i * (n + 1) + (j - 1)] + NW_GAP;
        let best = diag, dir = 0;
        if (up > best) { best = up; dir = 1; }
        if (left > best) { best = left; dir = 2; }
        F[i * (n + 1) + j] = best;
        Tb[i * (n + 1) + j] = dir;
      }
    }
    // Traceback — return aligned strings
    let i = m, j = n, ai = "", bj = "";
    while (i > 0 || j > 0) {
      const dir = Tb[i * (n + 1) + j];
      if (i > 0 && j > 0 && dir === 0) {
        ai = a[i - 1] + ai; bj = b[j - 1] + bj; i--; j--;
      } else if (i > 0 && dir === 1) {
        ai = a[i - 1] + ai; bj = "-" + bj; i--;
      } else {
        ai = "-" + ai; bj = b[j - 1] + bj; j--;
      }
    }
    return { a: ai, b: bj };
  }

  function progressiveMsa(seqs) {
    if (seqs.length === 0) return [];
    let consensus = seqs[0];
    let aligned = [consensus.split("")];
    for (let i = 1; i < seqs.length; i++) {
      const result = nwAlign(consensus, seqs[i]);
      // Re-align previous sequences to the new consensus to keep all the same length
      const newAligned = aligned.map((row) => {
        const s = row.join("");
        const r = nwAlign(result.a, s);
        return r.b.split("");
      });
      newAligned.push(result.b.split(""));
      consensus = result.a;
      aligned = newAligned;
    }
    // Pad all rows to consensus length
    const len = consensus.length;
    for (const row of aligned) {
      while (row.length < len) row.push("-");
    }
    return { consensus, aligned };
  }

  function conservationScores(msa) {
    const cons = msa.consensus;
    const n = cons.length;
    const out = new Array(n).fill(1);
    for (let j = 0; j < n; j++) {
      const ref = cons[j];
      if (ref === "-") { out[j] = 0; continue; }
      let matches = 0, total = 0;
      for (const row of msa.aligned) {
        const ch = row[j];
        if (ch === "-") continue;
        total++;
        if (ch === ref) matches++;
      }
      out[j] = total > 0 ? matches / total : 0;
    }
    return out;
  }

  function buildResidueIndex(pdbText) {
    const map = new Map(); // "chain:resi" -> {indexInSequence, resn}
    let seqIdx = -1;
    const seen = new Set();
    const lines = pdbText.split(/\r?\n/);
    for (const line of lines) {
      if (!line.startsWith("ATOM")) continue;
      const atomName = line.substring(12, 16).trim();
      if (atomName !== "CA") continue;
      const chain = line.substring(21, 22).trim() || "A";
      const resi = parseInt(line.substring(22, 26).trim(), 10);
      const resn = line.substring(17, 20).trim();
      const key = chain + ":" + resi;
      if (seen.has(key)) continue;
      seen.add(key);
      seqIdx++;
      map.set(key, { seqIdx, resn });
    }
    return { map, sequence: (() => {
      const seq = [];
      for (const [, v] of map) seq.push(v.seqIdx);  // not used
      // Reconstruct 1-letter sequence in seqIdx order
      const arr = new Array(seqIdx + 1).fill("X");
      for (const [k, v] of map) arr[v.seqIdx] = (AA3TO1[v.resn] || "X");
      return arr.join("");
    })() };
  }

  const api = {
    parseSequenceFromPDB,
    parseFastaSequences,
    nwAlign,
    progressiveMsa,
    conservationScores,
    buildResidueIndex,
    AA3TO1,
    AA1TO3,
  };

  global.PFE = global.PFE || {};
  global.PFE.seq = Object.freeze(api);
})(typeof globalThis !== "undefined" ? globalThis : this);
