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
    // Traceback — return aligned strings.
    //
    // Two things matter here:
    //   - The traceback must terminate. When j reaches 0 with i still > 0 the
    //     only legal move is "up" (consume a[i-1] against a gap), but the
    //     traceback matrix is never written for the j=0 column, so dir reads
    //     back as 0 and the old code fell into the "left" branch with
    //     b[j-1] === undefined. That decremented j past zero forever and the
    //     alignment string grew until V8 threw "Invalid string length". It
    //     only triggers when the two sequences differ a lot in length, e.g.
    //     aligning a 1200-residue UniProt sequence to a 300-residue construct.
    //   - The aligned strings are built with push + join rather than repeated
    //     prepending, which is O(n) instead of O(n^2).
    let i = m, j = n;
    const ai = [], bj = [];
    while (i > 0 || j > 0) {
      const dir = Tb[i * (n + 1) + j];
      if (i > 0 && j > 0 && dir === 0) {
        ai.push(a[i - 1]); bj.push(b[j - 1]); i--; j--;
      } else if (i > 0 && (dir === 1 || j === 0)) {
        ai.push(a[i - 1]); bj.push("-"); i--;
      } else {
        ai.push("-"); bj.push(b[j - 1]); j--;
      }
    }
    ai.reverse(); bj.reverse();
    return { a: ai.join(""), b: bj.join("") };
  }

  /**
   * Align a set of sequences to the first one.
   *
   * The columns are the reference sequence's residues. Every other sequence
   * is aligned to the reference on its own, so an insertion in one homolog
   * does not shift the columns of the sequences after it.
   *
   * Aligning each new sequence to a running consensus instead - the
   * obvious way to write a progressive alignment - re-aligns every earlier
   * row to an already gapped consensus on each pass, and the gap columns
   * accumulate. For a family of distant homologues the consensus grew to
   * several times the length of the reference: aligning E. coli adenylate
   * kinase to nine human adenylate kinases produced a 1,911-column
   * alignment of a 214-residue protein, with a mean identity of 10 %.
   *
   * @param {string[]} seqs  the first one is the reference
   * @returns {{consensus: string, aligned: string[][]}}
   */
  function progressiveMsa(seqs) {
    if (!seqs || !seqs.length) return { consensus: "", aligned: [] };
    const reference = seqs[0] || "";
    const len = reference.length;
    const rows = [];
    for (const s of seqs) {
      const r = nwAlign(reference, s || "");
      const chars = new Array(len).fill("-");
      let ref = 0;
      for (let k = 0; k < r.a.length && ref < len; k++) {
        // A gap in the reference side means this sequence has an extra
        // residue here. It has no reference position of its own, so it
        // does not get a column.
        if (r.a[k] === "-") continue;
        chars[ref] = r.b[k] === undefined ? "-" : r.b[k];
        ref++;
      }
      rows.push(chars);
    }
    return { consensus: reference, aligned: rows };
  }

  function conservationScores(msa) {
    return conservationProfile(msa).identity;
  }

  /**
   * Per-column conservation, in three flavours.
   *
   * `identity` is what the colour map uses: of the sequences that have the
   * column, the fraction that match the consensus. It is the right number
   * for shading, but on its own it is misleading - a position that exists
   * in only one sequence scores a perfect 1.0, because nothing disagrees
   * with it.
   *
   * `coverage` is the fraction of the whole alignment that has the column
   * at all, and `strength` is the fraction of the whole alignment that
   * matches the consensus. Strength is identity multiplied by coverage,
   * and it is the honest measure of a position: 1.0 means every aligned
   * sequence has the same residue there, 0.1 means one sequence has it
   * and nothing else agrees.
   *
   * @param {{consensus: string, aligned: string[]}} msa
   * @returns {{identity: number[], coverage: number[], strength: number[], columns: number, sequences: number}}
   */
  function conservationProfile(msa) {
    const cons = msa.consensus;
    const n = cons.length;
    const rows = msa.aligned || [];
    const totalSeqs = rows.length;
    const identity = new Array(n).fill(1);
    const coverage = new Array(n).fill(0);
    const strength = new Array(n).fill(0);
    for (let j = 0; j < n; j++) {
      const ref = cons[j];
      if (ref === "-") { identity[j] = 0; continue; }
      let matches = 0, present = 0;
      for (const row of rows) {
        const ch = row[j];
        // A row shorter than the consensus has no character here. That is
        // an absent residue, not a mismatching one.
        if (ch === undefined || ch === "-") continue;
        present++;
        if (ch === ref) matches++;
      }
      identity[j] = present > 0 ? matches / present : 0;
      coverage[j] = totalSeqs > 0 ? present / totalSeqs : 0;
      strength[j] = totalSeqs > 0 ? matches / totalSeqs : 0;
    }
    return { identity, coverage, strength, columns: n, sequences: totalSeqs };
  }

  /**
   * Parse the DBREF records of a PDB file.
   *
   * DBREF is the file's own statement of how its residue numbers relate to
   * the reference database's numbering:
   *
   *   DBREF  1M17 A  671   998  UNP    P00533   EGFR_HUMAN     695   1022
   *
   * which says PDB residues 671-998 are UniProt P00533 residues 695-1022.
   * That offset is what makes it possible to put a conservation score on a
   * residue of the structure: a crystal of a domain does not start at
   * residue 1 of the protein, so the two numberings cannot be assumed to
   * agree.
   *
   * @param {string} pdbText
   * @returns {Array<{chain: string, pdbBeg: number, pdbEnd: number, db: string,
   *                  acc: string, dbBeg: number, dbEnd: number}>}
   */
  function parseDbrefs(pdbText) {
    const out = [];
    if (typeof pdbText !== "string") return out;
    for (const line of pdbText.split(/\r?\n/)) {
      if (!line.startsWith("DBREF")) continue;
      // The fields are separated by runs of spaces, so splitting on
      // whitespace is more robust than fixed columns.
      const f = line.trim().split(/\s+/);
      if (f.length < 9) continue;
      const pdbBeg = parseInt(f[3], 10);
      const pdbEnd = parseInt(f[4], 10);
      const dbBeg = parseInt(f[8], 10);
      const dbEnd = parseInt(f[9], 10);
      if (!Number.isFinite(pdbBeg) || !Number.isFinite(dbBeg)) continue;
      out.push({
        chain: f[2],
        pdbBeg, pdbEnd: Number.isFinite(pdbEnd) ? pdbEnd : pdbBeg,
        db: f[5],
        acc: f[6],
        dbBeg, dbEnd: Number.isFinite(dbEnd) ? dbEnd : dbBeg,
      });
    }
    return out;
  }

  /**
   * Map each residue of the reference sequence to its column in the
   * alignment.
   *
   * The alignment grows as sequences are added, so the reference sequence's
   * residue n is not at column n. This counts the non-gap characters of the
   * reference row to find where each residue ended up.
   *
   * @param {string[]|string} alignedRow  the reference sequence's aligned row
   * @param {number} length  the length of the reference sequence
   * @returns {number[]}  index = residue index, value = alignment column
   */
  function canonicalColumnMap(alignedRow, length) {
    const row = typeof alignedRow === "string"
      ? alignedRow.split("")
      : (alignedRow || []);
    const out = new Array(Math.max(0, length | 0)).fill(-1);
    let ci = 0;
    for (let c = 0; c < row.length; c++) {
      if (row[c] === "-") continue;
      if (ci < out.length) out[ci] = c;
      ci++;
    }
    return out;
  }

  /**
   * Map "chain:resi" keys onto reference-database residue numbers.
   *
   * @param {Array<object>} dbrefs  from parseDbrefs
   * @param {string} accession  the reference accession the scores belong to
   * @param {Iterable<string>} residueKeys  "chain:resi" identifiers
   * @returns {Map<string, number>}  "chain:resi" -> reference position (1-based)
   */
  function referencePositionMap(dbrefs, accession, residueKeys) {
    const map = new Map();
    const refs = (dbrefs || []).filter((d) => d.acc === accession);
    if (!refs.length) return map;
    for (const key of residueKeys || []) {
      const sep = String(key).indexOf(":");
      if (sep < 0) continue;
      const chain = key.slice(0, sep);
      const resi = parseInt(key.slice(sep + 1), 10);
      if (!Number.isFinite(resi)) continue;
      for (const d of refs) {
        if (d.chain !== chain) continue;
        if (resi < d.pdbBeg || resi > d.pdbEnd) continue;
        map.set(key, d.dbBeg + (resi - d.pdbBeg));
        break;
      }
    }
    return map;
  }

  /**
   * Where a chain's sequence sits in a reference sequence.
   *
   * A crystal's residues are numbered by the depositor and the reference
   * sequence is numbered by the database, and the two disagree in ways
   * that are easy to get wrong: a construct usually starts part way into
   * the protein, and a UniProt sequence may still carry an initiator
   * methionine that the crystal has already lost. Matching the observed
   * sequence against the reference directly is self-correcting - wherever
   * the chain's residues actually come from is where they match.
   *
   * @param {string} canonical  the reference sequence
   * @param {string} chainSeq  the chain's residues, in order
   * @returns {{offset: number, identity: number}}
   */
  function sequenceOffset(canonical, chainSeq) {
    const out = { offset: -1, identity: 0 };
    if (!canonical || !chainSeq) return out;
    // The chain usually covers only part of the reference, so the offset
    // ranges over the whole reference and the score is the fraction of the
    // chain that matches. Requiring the chain to fit inside the reference
    // would rule out every fragment that starts part way into it.
    const need = Math.min(chainSeq.length, 20, canonical.length);
    let best = -1, bestScore = 0;
    for (let o = 0; o < canonical.length; o++) {
      const len = Math.min(chainSeq.length, canonical.length - o);
      if (len < need) break;   // the overlap only gets shorter from here
      let m = 0;
      for (let i = 0; i < len; i++) {
        if (canonical.charCodeAt(o + i) === chainSeq.charCodeAt(i)) m++;
      }
      const score = m / chainSeq.length;
      if (score > bestScore) { bestScore = score; best = o; }
    }
    if (best < 0) return out;
    out.offset = best;
    out.identity = bestScore;
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
      // Reconstruct the 1-letter sequence in seqIdx order.
      const arr = new Array(seqIdx + 1).fill("X");
      for (const v of map.values()) arr[v.seqIdx] = (AA3TO1[v.resn] || "X");
      return arr.join("");
    })() };
  }

  const api = {
    parseSequenceFromPDB,
    parseFastaSequences,
    nwAlign,
    progressiveMsa,
    conservationScores,
    conservationProfile,
    parseDbrefs,
    canonicalColumnMap,
    referencePositionMap,
    sequenceOffset,
    buildResidueIndex,
    AA3TO1,
    AA1TO3,
  };

  global.PFE = global.PFE || {};
  global.PFE.seq = Object.freeze(api);
})(typeof globalThis !== "undefined" ? globalThis : this);
