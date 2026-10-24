/* ===================================================================
 * Protein Flow Explorer — main app (research edition)
 *
 *  Core:
 *    - Library of curated structures (data/proteins.js)
 *    - Interactive 3Dmol.js viewer
 *    - Guided tours + trajectory playback
 *    - Custom molecule loader (PDB-ID, gene/UniProt, file drop)
 *
 *  Researcher features:
 *    - Distance / angle / dihedral measurement tool
 *    - Sequence panel with click-to-highlight + secondary structure
 *    - UniProt / AlphaFold lookup
 *    - Side-by-side structure comparison + per-residue RMSD
 *    - URL state for sharing exact views
 * =================================================================== */

(function () {
  "use strict";

  // ----- Extracted pure modules (see src/) -----
  // These hold the app's scientific core: sequence handling, molecular
  // geometry and pocket-grid geometry. They are pure functions with no DOM,
  // network or viewer dependency, so they can be unit-tested directly.
  // Destructuring keeps every existing call site unchanged.
  const { parseSequenceFromPDB, parseFastaSequences, nwAlign, progressiveMsa,
          conservationScores, buildResidueIndex, AA3TO1 } = PFE.seq;
  const { midPoint, computeDistance, computeAngle, computeDihedral,
          computePerResidueRmsd } = PFE.geometry;
  const { buildAtomHash, countAtomsNear, enclosureFraction } = PFE.pockets;
  const { STRUCT_CHAIN_PALETTE, ssColor, spectrumColor, bfactorColor,
          hydrophobicityColor, plddtColor, conservationColor } = PFE.colors;
  const { SIGNIFICANCE_COLORS, parseProteinChange,
          extractProteinChange, gnomadVariantId, classifySignificance,
          classifyImpactTier, classifyOncogenicity, significanceRank, rarityLabel,
          buildUniprotMap, mergeByChange } = PFE.variants;
  const { validatePdbId, validateGeneSymbol, validateUniprotAccession,
          validateClinVarRecord, validateGnomadResponse } = PFE.validate;
  const { parseMedianExpression, groupByCompartment, summarizeExpression } = PFE.expression;

  // ====== GLOBAL STATE ======
  const state = {
    viewer: null,
    currentProtein: null,
    currentPdbIndex: 0,
    models: [],

    // Style / display
    currentStyle: "cartoon",
    currentColor: "chain",
    currentOpacity: 1,
    background: "white",
    spinning: false,

    // Value ranges of the loaded structure, used to scale the spectrum and
    // B-factor colour ramps. Recomputed whenever a structure is loaded.
    colorRange: { resiMin: 0, resiMax: 1, bMin: 0, bMax: 1 },

    // Frame counter
    fpsFrames: 0,
    fpsLastTime: performance.now(),

    // Tour
    tourRunning: false,
    tourStepIndex: 0,
    tourTimeout: null,

    // Trajectory
    trajectoryRunning: false,
    trajectoryRaf: null,
    trajectoryValue: 50,
    trajectoryDirection: 1,

    activeHighlights: [],
    currentFilter: "all",
    customCounter: 0,

    // ----- Researcher extensions -----
    pickMode: false,
    measurements: [],      // [{type:"distance"|"angle"|"dihedral", atoms:[{serial,x,y,z,chain,resi,resn,atom,elem}], value, label}]
    pickBuffer: [],         // atoms queued for current measurement
    pickShapes: [],         // 3Dmol shapes representing measurements
    pickLabels: [],

    sequence: null,         // parsed from PDB
    sequenceShowSS: false,

    compareProtein: null,   // second protein for comparison
    compareLoaded: false,
    aligned: false,

    // ----- Bioinformatics analysis results -----
    conservation: null,      // {byResidue: Map "chain:resi" -> 0..1, homologs: number}
    conservationJob: null,   // Promise tracking in-flight fetch
    variants: [],            // [{resi, chain, mutation, significance, color}]
    variantByResidue: new Map(),   // "chain:resi" -> colour, for the variant scheme
    expression: null,              // GTEx tissue expression for the loaded gene
    variantRateLimited: false,     // gnomAD refused some frequency lookups
    variantLookupErrors: 0,        // gnomAD lookups that failed outright
    _uniprotInfo: null,      // cached {accession, gene, sequence} for the loaded structure
    pockets: [],             // [{id, center:[x,y,z], volume, score, residues:[{chain,resi}]}]
    ligandContacts: null,    // {ligandName, residue, contacts:[{chain,resi,resn,distance,type}]}
    rmsdPerResidue: null,   // array of {resi, chain, value}
    rmsdColorMode: false,

    // ----- Metadata / ligands / interactions (mol-bio + drug discovery) -----
    meta: null,             // RCSB entry metadata (method, resolution, title…)
    ligands: [],            // [{compId, name, chains, entityId}]
    activeLigand: null,
    interactions: null,     // {hbonds:[], salts:[], pi:[], hydro:[]}
    hiddenChains: new Set(),
  };

  // ====== DOM REFS ======
  const $ = (id) => document.getElementById(id);
  const dom = {
    proteinList: $("protein-list"),
    filterChips: $("filter-chips"),
    viewer: $("viewer"),
    loadingOverlay: $("overlay-loading"),
    loadingText: $("loading-text"),
    annotation: $("overlay-annotation"),
    annTitle: $("ann-title"),
    annBody: $("ann-body"),
    stepOverlay: $("overlay-step"),
    stepCurrent: $("step-current"),
    stepTotal: $("step-total"),
    pickHint: $("pick-hint"),
    pickHintText: $("pick-hint-text"),
    btnCancelPick: $("btn-cancel-pick"),
    currentProteinName: $("current-protein-name"),
    currentStateBadge: $("current-state-badge"),
    fpsCounter: $("fps-counter"),

    // Controls
    ctrlStyle: $("ctrl-style"),
    ctrlColor: $("ctrl-color"),
    ctrlOpacity: $("ctrl-opacity"),
    ctrlBg: $("ctrl-bg"),

    btnResetView: $("btn-reset-view"),
    btnSpin: $("btn-spin"),
    btnScreenshot: $("btn-screenshot"),
    btnHighlightHeme: $("btn-highlight-heme"),
    btnHighlightActive: $("btn-highlight-active"),
    btnHighlightInterface: $("btn-highlight-interface"),
    btnHighlightAll: $("btn-highlight-all"),
    btnClearHighlights: $("btn-clear-highlights"),

    proteinDescription: $("protein-description"),

    // Custom loader
    loaderTabs: $("loader-tabs"),
    pdbInput: $("pdb-input"),
    btnLoadPdb: $("btn-load-pdb"),
    geneInput: $("gene-input"),
    btnLoadGene: $("btn-load-gene"),
    dropZone: $("drop-zone"),
    customError: $("custom-error"),
    compareControls: $("compare-controls"),
    compareInput: $("compare-input"),
    btnLoadCompare: $("btn-load-compare"),
    btnClearCompare: $("btn-clear-compare"),

    // Sequence panel
    sequenceChains: $("sequence-chains"),
    sequenceInfo: $("sequence-info"),
    btnToggleSecondary: $("btn-toggle-secondary"),
    sequenceSecondary: $("sequence-secondary"),

    // Tabs
    tabBtns: document.querySelectorAll(".tab-btn"),
    tabPanels: document.querySelectorAll(".tab-panel"),

    // Tour
    tourEmpty: $("tour-empty"),
    tourContent: $("tour-content"),
    tourTitle: $("tour-title"),
    tourSteps: $("tour-steps"),
    btnTourPlay: $("btn-tour-play"),
    btnTourNext: $("btn-tour-next"),
    btnTourStop: $("btn-tour-stop"),
    tourProgressBar: $("tour-progress-bar"),

    // Trajectory
    trajectoryEmpty: $("trajectory-empty"),
    trajectoryContent: $("trajectory-content"),
    trajectoryTitle: $("trajectory-title"),
    trajectoryDesc: $("trajectory-desc"),
    trajSlider: $("traj-slider"),
    trajPercent: $("traj-percent"),
    trajStateA: $("traj-state-a"),
    trajStateB: $("traj-state-b"),
    btnTrajPlay: $("btn-traj-play"),
    btnTrajPause: $("btn-traj-pause"),
    btnTrajReset: $("btn-traj-reset"),
    trajSpeed: $("traj-speed"),

    // Tools (measurement + compare)
    btnPickAtom: $("btn-pick-atom"),
    btnClearMeasurements: $("btn-clear-measurements"),
    measurementTable: $("measurement-table"),
    btnEnableCompare: $("btn-enable-compare"),
    btnAlignCalpha: $("btn-align-calpha"),
    btnColorByRmsd: $("btn-color-by-rmsd"),
    rmsdSummary: $("rmsd-summary"),

    // Bioinformatics analyses
    btnColorPlddt: $("btn-color-plddt"),
    btnFetchConservation: $("btn-fetch-conservation"),
    btnClearConservation: $("btn-clear-conservation"),
    conservationStatus: $("conservation-status"),
    btnFetchVariants: $("btn-fetch-variants"),
    btnClearVariants: $("btn-clear-variants"),
    btnFetchExpression: $("btn-fetch-expression"),
    btnClearExpression: $("btn-clear-expression"),
    expressionStatus: $("expression-status"),
    variantsStatus: $("variants-status"),
    btnDetectPockets: $("btn-detect-pockets"),
    btnAnalyzeContacts: $("btn-analyze-contacts"),
    btnClearBioGeo: $("btn-clear-bio-geo"),
    pocketsStatus: $("pockets-status"),
    contactsStatus: $("contacts-status"),
    pocketList: $("pocket-list"),
    contactList: $("contact-list"),

    // Metadata / ligand / navigation / export
    metaGrid: $("meta-grid"),
    ligandList: $("ligand-list"),
    ligandSelect: $("ligand-select"),
    ligandPropStatus: $("ligand-prop-status"),
    propGrid: $("prop-grid"),
    resiChain: $("resi-chain"),
    resiNum: $("resi-num"),
    btnGoResidue: $("btn-go-residue"),
    chainToggles: $("chain-toggles"),
    btnExportPdb: $("btn-export-pdb"),
    btnExportCsv: $("btn-export-csv"),
    exportStatus: $("export-status"),

    // Drug design
    btnHbond: $("btn-hbonds"),
    btnClearInteractions: $("btn-clear-interactions"),
    hbondStatus: $("hbond-status"),
    interactionSummary: $("interaction-summary"),
    mutInput: $("mut-input"),
    btnCheckMutation: $("btn-check-mutation"),
    mutationStatus: $("mutation-status"),

    // Modal / share
    // Presentation tab (bioinformatics record view)
    hierarchy: $("hierarchy"),
    xrefList: $("xref-list"),
    annotKey: $("annot-key"),
    uncertainty: $("uncertainty"),
    citeBlock: $("cite-block"),
    btnCopyCite: $("btn-copy-cite"),

    btnAbout: $("btn-about"),
    aboutModal: $("about-modal"),
    btnCloseAbout: $("btn-close-about"),
    btnShare: $("btn-share"),

    toast: $("toast"),
  };

  // ====== UTILITIES ======
  function escapeHtml(s) {
    return String(s).replace(/[&<>"']/g, (c) => ({
      "&": "&amp;",
      "<": "&lt;",
      ">": "&gt;",
      '"': "&quot;",
      "'": "&#39;",
    })[c]);
  }

  function categoryLabel(c) {
    return {
      allostery: "allostery",
      binding: "binding",
      "drug-target": "drug target",
      folding: "folding",
      domains: "domains",
      custom: "custom",
    }[c] || c;
  }

  function showToast(msg, kind = "") {
    if (!dom.toast) return;
    dom.toast.textContent = msg;
    dom.toast.className = "toast " + kind;
    dom.toast.hidden = false;
    clearTimeout(showToast._t);
    showToast._t = setTimeout(() => {
      dom.toast.hidden = true;
    }, 3500);
  }

  function debounce(fn, ms) {
    let t = null;
    return (...args) => {
      clearTimeout(t);
      t = setTimeout(() => fn(...args), ms);
    };
  }

  // ====== NETWORK LAYER ======
  // Every outbound request in the app goes through net(). This gives one place
  // to enforce a deadline, retry transient failures with backoff, and cancel
  // work that the user has moved on from.
  //
  // Before this existed, 10 of 14 call sites used bare fetch(): a stalled RCSB
  // or UniProt response left the UI spinning indefinitely with no feedback and
  // no way out. Worse, switching proteins mid-load let a slow earlier response
  // land after a faster later one and overwrite the current view.
  const NET_TIMEOUT_MS = 20000;
  const NET_RETRIES = 1;          // one retry covers the common transient case
  const NET_BACKOFF_MS = 600;
  const NET_RATE_BACKOFF_MS = 3000;   // used when a server answers 429

  // Bumped on every protein switch so late responses from a previous load are
  // discarded instead of clobbering the new one.
  let netGeneration = 0;
  function newNetworkGeneration() { netGeneration++; return netGeneration; }

  function isAbort(err) {
    return err && (err.name === "AbortError" || err.__stale === true);
  }

  /**
   * Fetch with a hard deadline and bounded retries.
   * @param {string} url
   * @param {object} opts
   * @param {number} gen  network generation this request belongs to
   */
  async function net(url, opts, gen) {
    opts = opts || {};
    if (gen !== undefined && gen !== netGeneration) {
      const stale = new Error("stale request");
      stale.__stale = true;
      throw stale;
    }
    let lastErr;
    for (let attempt = 0; attempt <= NET_RETRIES; attempt++) {
      if (gen !== undefined && gen !== netGeneration) {
        const stale = new Error("stale request");
        stale.__stale = true;
        throw stale;
      }
      const ctrl = typeof AbortController !== "undefined" ? new AbortController() : null;
      const timer = setTimeout(() => { if (ctrl) ctrl.abort(); }, NET_TIMEOUT_MS);
      let rateLimited = false;
      try {
        const res = await fetch(url, ctrl ? { ...opts, signal: ctrl.signal } : opts);
        // 429 is the one 4xx that is worth retrying: it means "you are going
        // too fast", not "this does not exist". NCBI E-utilities answer 429
        // when you exceed three requests a second without an API key, and
        // the variant fetcher makes several requests in a row.
        if (res.status === 429) {
          rateLimited = true;
          lastErr = new Error("HTTP 429 (rate limited)");
        } else if (res.ok || res.status < 500) {
          return res;                    // 4xx is a real answer, not a retry
        } else {
          lastErr = new Error(`HTTP ${res.status}`);
        }
      } catch (err) {
        if (ctrl && err && err.name === "AbortError") {
          lastErr = new Error(`timeout after ${NET_TIMEOUT_MS / 1000}s`);
        } else {
          lastErr = err;
        }
      } finally {
        clearTimeout(timer);
      }
      if (attempt < NET_RETRIES) {
        const wait = rateLimited ? NET_RATE_BACKOFF_MS : NET_BACKOFF_MS * (attempt + 1);
        await new Promise((r) => setTimeout(r, wait));
      }
    }
    throw lastErr || new Error("request failed");
  }



  // ====== URL STATE ENCODE / DECODE ======
  function encodeState() {
    const params = new URLSearchParams();
    if (state.currentProtein) {
      params.set("p", state.currentProtein.id);
      if (state.currentPdbIndex > 0) params.set("i", String(state.currentPdbIndex));
    }
    params.set("style", state.currentStyle);
    params.set("color", state.currentColor);
    params.set("op", String(state.currentOpacity));
    params.set("bg", state.background);
    if (state.tourRunning || state.tourStepIndex > 0) {
      params.set("tour", String(state.tourStepIndex));
    }
    return params;
  }

  function applyStateFromUrl() {
    const params = new URLSearchParams(window.location.search);
    const proteinId = params.get("p");
    const pdbIndex = parseInt(params.get("i") || "0", 10);
    const style = params.get("style");
    const color = params.get("color");
    const opacity = parseFloat(params.get("op") || "");
    const bg = params.get("bg");
    return { proteinId, pdbIndex, style, color, opacity, bg, tour: parseInt(params.get("tour") || "", 10) };
  }

  const updateUrlDebounced = debounce(() => {
    const params = encodeState();
    const url = params.toString()
      ? window.location.pathname + "?" + params.toString()
      : window.location.pathname;
    window.history.replaceState({}, "", url);
  }, 400);

  function shareCurrentView() {
    const params = encodeState();
    const url = window.location.origin + window.location.pathname + "?" + params.toString();
    if (navigator.clipboard) {
      navigator.clipboard.writeText(url).then(
        () => showToast("Share link copied to clipboard", "success"),
        () => prompt("Copy this link:", url)
      );
    } else {
      prompt("Copy this link:", url);
    }
    return url;
  }

  // ====== VIEWER SETUP ======
  function initViewer() {
    state.viewer = $3Dmol.createViewer(dom.viewer, {
      backgroundColor: "white",
    });

    // FPS counter
    const tick = () => {
      state.fpsFrames++;
      const now = performance.now();
      if (now - state.fpsLastTime >= 1000) {
        dom.fpsCounter.textContent = `${state.fpsFrames} fps`;
        state.fpsFrames = 0;
        state.fpsLastTime = now;
      }
      requestAnimationFrame(tick);
    };
    requestAnimationFrame(tick);

    if (window.ResizeObserver) {
      const ro = new ResizeObserver(() => state.viewer.resize());
      ro.observe(dom.viewer);
    } else {
      window.addEventListener("resize", () => state.viewer.resize());
    }

    // Manual click handler: 3Dmol.js 2.4.2 doesn't expose setClickHandler, so we listen
    // on the canvas and call getIntersection to find the atom under the pointer.
    let clickTimer = null;
    const canvas = dom.viewer.querySelector("canvas");
    if (canvas) {
      canvas.addEventListener("click", (e) => {
        if (!state.pickMode) return;
        // Distinguish click from drag-rotate (3Dmol handles rotate via its own listeners)
        if (clickTimer) clearTimeout(clickTimer);
        clickTimer = setTimeout(() => {
          clickTimer = null;
          const rect = canvas.getBoundingClientRect();
          const x = e.clientX - rect.left;
          const y = e.clientY - rect.top;
          const atom = state.viewer.getIntersection(x, y);
          if (atom) handleAtomPicked(atom);
        }, 50);
      });
    }
  }

  // ====== PROTEIN LIBRARY ======
  function renderProteinList(filter = "all") {
    state.currentFilter = filter;
    const list = filter === "all"
      ? PROTEINS
      : PROTEINS.filter((p) => p.category === filter);
    dom.proteinList.innerHTML = "";
    list.forEach((p) => {
      const card = document.createElement("div");
      card.className = "protein-card";
      if (p.isCustom) card.classList.add("protein-card-custom");
      card.dataset.id = p.id;
      card.style.borderLeftColor = p.color;
      if (state.currentProtein && state.currentProtein.id === p.id) {
        card.classList.add("active");
      }
      const pdbTags = p.pdbIds
        .map((id, i) => `<span class="pdb-tag">${escapeHtml(String(id))}${p.stateLabels.length > 1 ? "·" + escapeHtml(p.stateLabels[i].split(" ")[0]) : ""}</span>`)
        .join("");
      const removeBtn = p.isCustom
        ? `<button class="remove-btn" data-remove="${escapeHtml(p.id)}" title="Remove">×</button>`
        : "";
      card.innerHTML = `
        ${removeBtn}
        <div class="protein-card-name">
          <span style="color: ${p.color}">●</span> ${escapeHtml(p.name)}
        </div>
        <div class="protein-card-summary">${escapeHtml(p.summary)}</div>
        <div class="protein-card-meta">
          ${pdbTags}
          <span class="category-tag">${categoryLabel(p.category)}</span>
        </div>
      `;
      card.addEventListener("click", (e) => {
        if (e.target.closest(".remove-btn")) return;
        loadProtein(p.id);
      });
      const removeEl = card.querySelector(".remove-btn");
      if (removeEl) {
        removeEl.addEventListener("click", (e) => {
          e.stopPropagation();
          removeCustomMolecule(p.id);
        });
      }
      dom.proteinList.appendChild(card);
    });
  }

  // Filter chips
  dom.filterChips.addEventListener("click", (e) => {
    const btn = e.target.closest(".chip");
    if (!btn) return;
    document.querySelectorAll(".chip").forEach((c) => c.classList.remove("active"));
    btn.classList.add("active");
    renderProteinList(btn.dataset.filter);
  });

  // ====== CUSTOM MOLECULES ======
  function showCustomError(msg) {
    dom.customError.textContent = msg;
    dom.customError.hidden = false;
    clearTimeout(showCustomError._t);
    showCustomError._t = setTimeout(() => {
      dom.customError.hidden = true;
    }, 6000);
  }

  function extractPdbTitle(pdbText) {
    if (!pdbText) return null;
    const m = pdbText.match(/^TITLE\s+(.+?)$/m);
    return m ? m[1].trim().replace(/\s+/g, " ") : null;
  }

  function makeCustomProtein({ id, name, summary, description, rawData, format, source }) {
    const title = extractPdbTitle(rawData);
    if (title && (name === undefined || name.startsWith("PDB ") || name === "Custom Molecule")) {
      name = title.length > 40 ? title.slice(0, 40) + "…" : title;
    }
    return {
      id,
      name: name || "Custom Molecule",
      pdbIds: [id.replace(/^custom-/, "")],
      stateLabels: [format === "cif" ? "mmCIF" : "PDB"],
      summary: summary || (source === "pdb-id" ? `Loaded from PDB ID ${id.replace(/^custom-/, "")}` : "Custom-loaded molecule"),
      description: description || title || `Custom molecule loaded from ${source}`,
      category: "custom",
      color: "#a78bfa",
      isCustom: true,
      rawData,
      format: format || "pdb",
    };
  }

  async function loadCustomByPdbId(pdbIdRaw) {
    // Validate before the value goes anywhere near a URL.
    const checked = validatePdbId(pdbIdRaw);
    if (!checked.ok) {
      showCustomError(checked.error);
      return;
    }
    const pdbId = checked.value;
    const id = `custom-${pdbId}`;
    if (PROTEINS.some((p) => p.id === id)) {
      hideLoading();
      loadProtein(id);
      dom.pdbInput.value = "";
      return;
    }
    showLoading(`Fetching PDB ${pdbId}…`);
    try {
      const response = await net(`https://files.rcsb.org/download/${pdbId}.pdb`);
      if (!response.ok) throw new Error(`HTTP ${response.status}`);
      const data = await response.text();
      if (!data || data.length < 50) throw new Error("Empty or invalid response");
      const protein = makeCustomProtein({
        id, rawData: data, format: "pdb", source: "pdb-id",
      });
      PROTEINS.push(protein);
      dom.customError.hidden = true;
      dom.pdbInput.value = "";
      hideLoading();
      renderProteinList(state.currentFilter);
      loadProtein(id);
    } catch (err) {
      hideLoading();
      showCustomError(`Could not load ${pdbId}: ${err.message}`);
    }
  }

  async function loadCustomFromText(text, filename) {
    const fname = filename || "pasted";
    const ext = (fname.split(".").pop() || "").toLowerCase();
    let format = "pdb";
    if (ext === "cif" || /^data_/m.test(text.slice(0, 500))) format = "cif";
    if (!text || text.length < 30) {
      showCustomError("File appears empty or too small to be a structure.");
      return;
    }
    state.customCounter++;
    const id = `custom-file-${state.customCounter}`;
    const baseName = fname.replace(/\.(pdb|cif)$/i, "").replace(/[^\w.-]+/g, "_");
    const protein = makeCustomProtein({
      id,
      name: baseName || `Molecule ${state.customCounter}`,
      summary: `Loaded from file: ${fname}`,
      description: `Custom molecule from local file (${format.toUpperCase()})`,
      rawData: text, format, source: "file",
    });
    PROTEINS.push(protein);
    dom.customError.hidden = true;
    renderProteinList(state.currentFilter);
    loadProtein(id);
  }

  async function loadCustomFromFile(file) {
    if (!file) return;
    try {
      const text = await file.text();
      await loadCustomFromText(text, file.name);
    } catch (err) {
      showCustomError(`Could not read file: ${err.message}`);
    }
  }

  function removeCustomMolecule(id) {
    const idx = PROTEINS.findIndex((p) => p.id === id);
    if (idx === -1) return;
    PROTEINS.splice(idx, 1);
    if (state.currentProtein && state.currentProtein.id === id) {
      const fallback = PROTEINS.find((p) => !p.isCustom);
      if (fallback) loadProtein(fallback.id);
    }
    renderProteinList(state.currentFilter);
  }

  // Wire up custom loader UI
  dom.loaderTabs.addEventListener("click", (e) => {
    const btn = e.target.closest(".loader-tab");
    if (!btn) return;
    document.querySelectorAll(".loader-tab").forEach((t) => t.classList.remove("active"));
    btn.classList.add("active");
    document.querySelectorAll(".custom-loader-pane").forEach((p) => {
      p.hidden = p.dataset.loaderPane !== btn.dataset.loader;
    });
  });
  dom.btnLoadPdb.addEventListener("click", () => loadCustomByPdbId(dom.pdbInput.value));
  dom.pdbInput.addEventListener("keydown", (e) => {
    if (e.key === "Enter") {
      e.preventDefault();
      loadCustomByPdbId(dom.pdbInput.value);
    }
  });
  dom.pdbInput.addEventListener("input", (e) => {
    e.target.value = e.target.value.toUpperCase().replace(/[^0-9A-Z]/g, "");
  });
  dom.btnLoadGene.addEventListener("click", () => loadByGeneOrUniprot(dom.geneInput.value));
  dom.geneInput.addEventListener("keydown", (e) => {
    if (e.key === "Enter") {
      e.preventDefault();
      loadByGeneOrUniprot(dom.geneInput.value);
    }
  });

  ["dragenter", "dragover"].forEach((ev) => {
    dom.dropZone.addEventListener(ev, (e) => {
      e.preventDefault();
      e.stopPropagation();
      dom.dropZone.classList.add("dragover");
    });
  });
  ["dragleave", "dragend"].forEach((ev) => {
    dom.dropZone.addEventListener(ev, (e) => {
      e.preventDefault();
      e.stopPropagation();
      dom.dropZone.classList.remove("dragover");
    });
  });
  dom.dropZone.addEventListener("drop", (e) => {
    e.preventDefault();
    e.stopPropagation();
    dom.dropZone.classList.remove("dragover");
    const files = e.dataTransfer && e.dataTransfer.files;
    if (files && files.length > 0) loadCustomFromFile(files[0]);
  });
  window.addEventListener("dragover", (e) => e.preventDefault());
  window.addEventListener("drop", (e) => e.preventDefault());
  dom.dropZone.addEventListener("click", () => {
    const input = document.createElement("input");
    input.type = "file";
    input.accept = ".pdb,.cif,.ent,.txt";
    input.addEventListener("change", () => {
      if (input.files && input.files[0]) loadCustomFromFile(input.files[0]);
    });
    input.click();
  });

  // ====== UNIPROT / GENE LOOKUP ======
  async function loadByGeneOrUniprot(queryRaw) {
    // Upper-case up front: gene symbols and accessions are case-insensitive
    // in practice, and the validators are deliberately strict about shape.
    const query = String(queryRaw || "").trim().toUpperCase();
    if (!query) {
      showCustomError("Type a gene name (e.g. TP53) or UniProt accession (e.g. P04637).");
      return;
    }

    showLoading(`Looking up ${query}…`);

    try {
      let accession = null;

      // Detect the input type by trying both, rather than by guessing from a
      // shape: a gene symbol and an accession are both "letters and digits",
      // and the previous guess sent a mixed-case symbol like C11orf87 down
      // the accession branch.
      const asAccession = validateUniprotAccession(query);
      const asGene = validateGeneSymbol(query);

      if (asAccession.ok) {
        accession = asAccession.value;
      } else if (asGene.ok) {
        // Treat as gene name — query UniProt
        const url = `https://rest.uniprot.org/uniprotkb/search?query=gene:${encodeURIComponent(asGene.value)}+AND+organism_id:9606&format=json&size=1&fields=accession,gene_names,protein_name`;
        const r = await net(url);
        if (!r.ok) throw new Error(`UniProt search returned ${r.status}`);
        const j = await r.json();
        const hit = j.results && j.results[0];
        if (!hit) {
          hideLoading();
          showCustomError(`No UniProt entry found for "${query}".`);
          return;
        }
        accession = hit.primaryAccession;
      } else {
        hideLoading();
        showCustomError(`${asGene.error} Give a UniProt accession instead, e.g. P04637.`);
        return;
      }

      // Try to find PDB cross-references
      showLoading(`Fetching cross-refs for ${accession}…`);
      const entryUrl = `https://rest.uniprot.org/uniprotkb/${accession}.json`;
      const er = await net(entryUrl);
      if (!er.ok) throw new Error(`UniProt entry fetch returned ${er.status}`);
      const ej = await er.json();
      const pdbRefs = (ej.uniProtKBCrossReferences || [])
        .filter((x) => x.database === "PDB")
        .sort((a, b) => {
          const al = parseInt((a.properties || []).find((p) => p.key === "Resolution")?.value || "999", 10);
          const bl = parseInt((b.properties || []).find((p) => p.key === "Resolution")?.value || "999", 10);
          return al - bl;
        });

      hideLoading();

      if (pdbRefs.length > 0) {
        const pdbId = pdbRefs[0].id;
        showToast(`Loaded PDB ${pdbId} from UniProt ${accession}`, "success");
        await loadCustomByPdbId(pdbId);
      } else {
        // Fallback: AlphaFold
        showLoading(`No PDB — fetching AlphaFold model for ${accession}…`);
        const afUrl = `https://alphafold.ebi.ac.uk/files/AF-${accession}-F1-model_v4.pdb`;
        const afr = await net(afUrl);
        if (!afr.ok) {
          hideLoading();
          showCustomError(`No PDB or AlphaFold model for ${accession}.`);
          return;
        }
        const data = await afr.text();
        state.customCounter++;
        const id = `custom-alphafold-${accession}`;
        const protein = makeCustomProtein({
          id,
          name: `AlphaFold ${accession}`,
          summary: `AlphaFold model for UniProt ${accession}`,
          description: `Predicted structure from AlphaFold DB`,
          rawData: data, format: "pdb", source: "alphafold",
        });
        PROTEINS.push(protein);
        dom.geneInput.value = "";
        hideLoading();
        renderProteinList(state.currentFilter);
        loadProtein(id);
      }
    } catch (err) {
      hideLoading();
      showCustomError(`Lookup failed: ${err.message}`);
    }
  }

  // ====== PROTEIN LOADING ======
  async function loadProtein(proteinId) {
    const protein = PROTEINS.find((p) => p.id === proteinId);
    if (!protein) return;

    // Claim a network generation: any request still in flight from a previous
    // protein will be discarded on arrival instead of overwriting this one.
    const gen = newNetworkGeneration();
    state._netGen = gen;

    // If switching to a different protein, clear comparison
    if (state.compareProtein && state.compareProtein.id !== proteinId) {
      clearComparison();
    }
    // Reset measurement state when loading a new protein
    clearMeasurements();

    state.currentProtein = protein;
    state.currentPdbIndex = 0;
    state._chainOrder = []; // reset chain->colour mapping for the new structure

    dom.currentProteinName.textContent = protein.name;
    updateStateBadge(0);
    renderProteinList(document.querySelector(".chip.active").dataset.filter);
    dom.proteinDescription.textContent = protein.description;

    updateTourAvailability();
    updateTrajectoryAvailability();
    updateHighlightButtons();
    updateMeasurementsUI();

    showLoading("Loading " + protein.name + "…");

    state.viewer.removeAllModels();
    state.viewer.removeAllShapes();
    state.viewer.removeAllLabels();
    state.viewer.removeAllSurfaces();

    try {
      const loadPromises = protein.pdbIds.map(async (pdbId, idx) => {
        let data;
        if (protein.isCustom && protein.rawData && idx === 0) {
          data = protein.rawData;
        } else {
          const url = `https://files.rcsb.org/download/${pdbId}.pdb`;
          const response = await net(url, undefined, gen);
          if (!response.ok) throw new Error(`HTTP ${response.status} fetching ${pdbId}`);
          data = await response.text();
          if (!data) throw new Error(`Empty response for ${pdbId}`);
        }
        const format = protein.isCustom && idx === 0 ? (protein.format || "pdb") : "pdb";
        const model = state.viewer.addModel(data, format);
        model.modelIndex = idx;
        model.pdbId = pdbId;
        return { model, data, format };
      });

      const results = await Promise.all(loadPromises);
      state.models = results.map((r) => r.model);
      // Cache the raw PDB text of the first model for sequence parsing
      state._lastPdbText = results[0].data;
      state._lastPdbFormat = results[0].format;

      computeColorRanges(state.models[0]);

      applyAllStyles();
      state.viewer.zoomTo();
      state.viewer.render();

      state.spinning = false;
      dom.btnSpin.textContent = "↻ Spin";

      hideLoading();

      if (TRAJECTORIES[protein.id]) {
        setupTrajectoryInitial();
      } else {
        setModelVisibility(0);
      }

      // Parse sequence + SS for sequence panel
      parseAndRenderSequence();

      // Chain visibility toggles + RCSB metadata (method/resolution/ligands)
      renderChainToggles();
      loadStructureMetadata();

      // Update URL
      updateUrlDebounced();
    } catch (err) {
      // Superseded by a newer protein load, or the user navigated away: the new
      // load owns the spinner, so don't stomp on it with a stale error.
      if (isAbort(err) || gen !== netGeneration) return;
      console.error(err);
      dom.loadingText.textContent = `Could not load ${protein.name}: ${err.message}`;
    }
  }

  function updateStateBadge(idx) {
    const protein = state.currentProtein;
    if (!protein || !protein.stateLabels[idx]) {
      dom.currentStateBadge.textContent = "—";
      return;
    }
    dom.currentStateBadge.textContent = protein.stateLabels[idx];
  }

  function showLoading(text) {
    dom.loadingText.textContent = text;
    dom.loadingOverlay.classList.add("active");
  }

  function hideLoading() {
    dom.loadingOverlay.classList.remove("active");
  }

  // ====== STYLING ======

  // 3Dmol's built-in 'chain' scheme is a poor palette: chains B/J/O/P/R/T render as
  // grey or near-black, and the rest are washed-out pastels. We use our own
  // colourblind-friendly palette via colorfunc so every chain is clearly distinct
  // regardless of its letter. The palette itself lives in src/colors.js.
  function chainColorFunc(atom) {
    const idx = chainIndex(atom.chain);
    return STRUCT_CHAIN_PALETTE[idx % STRUCT_CHAIN_PALETTE.length];
  }
  // Stable per-structure chain ordering (A, B, C …) so colours don't shuffle between models.
  function chainIndex(chainId) {
    if (!state._chainOrder) state._chainOrder = [];
    const key = (chainId || "A").trim() || "A";
    let i = state._chainOrder.indexOf(key);
    if (i === -1) {
      i = state._chainOrder.length;
      state._chainOrder.push(key);
    }
    return i;
  }

  // Value ranges of the loaded structure, so the spectrum and B-factor ramps
  // span the actual data instead of a fixed 0–100 window.
  function computeColorRanges(model) {
    const atoms = safeSelect(model, {});
    let resiMin = Infinity, resiMax = -Infinity;
    let bMin = Infinity, bMax = -Infinity;
    for (const a of atoms) {
      if (typeof a.resi === "number" && Number.isFinite(a.resi)) {
        if (a.resi < resiMin) resiMin = a.resi;
        if (a.resi > resiMax) resiMax = a.resi;
      }
      if (typeof a.b === "number" && Number.isFinite(a.b)) {
        if (a.b < bMin) bMin = a.b;
        if (a.b > bMax) bMax = a.b;
      }
    }
    state.colorRange = {
      resiMin: Number.isFinite(resiMin) ? resiMin : 0,
      resiMax: Number.isFinite(resiMax) ? resiMax : 1,
      bMin: Number.isFinite(bMin) ? bMin : 0,
      bMax: Number.isFinite(bMax) ? bMax : 1,
    };
  }

  function applyAllStyles() {
    state.models.forEach((m, idx) => applyStyleToModel(m, idx));
    // Re-apply RMSD coloring if active
    if (state.rmsdColorMode && state.compareLoaded) applyRmsdColors();
  }

  function applyStyleToModel(model, idx) {
    const sel = { model: idx };
    const style = getStyleSpec(state.currentStyle, state.currentColor, state.currentOpacity);
    state.viewer.setStyle(sel, style);

    const protein = state.currentProtein;
    if (protein && protein.highlights && protein.highlights.ligand) {
      state.viewer.setStyle({ model: idx, resn: protein.highlights.ligand.resn }, {
        stick: { color: "magenta", radius: 0.18 },
      });
    }
  }

  // ----- colour functions -----
  // Each of these returns an RGB integer for one atom and is handed to 3Dmol
  // as a `colorfunc`. They are used instead of 3Dmol's own schemes so that
  // the colours match the annotation key exactly, and so that the same colour
  // applies to every representation (cartoon, stick, sphere, line, surface).

  const DOMAIN_FALLBACK_COLOR = 0x9aa0a6;     // grey — residue outside any domain
  const CONSERVATION_MISSING_COLOR = 0xcccccc; // grey — residue not in the alignment
  const VARIANT_BACKGROUND_COLOR = 0xbcc7d1;  // muted blue-grey — no variant here

  /** Colour by the annotated domain of the residue. */
  function domainColorFunc(atom) {
    const domains = state.currentProtein && state.currentProtein.domains;
    if (!domains || !domains.length) return chainColorFunc(atom);
    const resi = atom.resi;
    if (typeof resi !== "number") return DOMAIN_FALLBACK_COLOR;
    for (const d of domains) {
      for (const r of d.ranges) {
        if (r.chain && r.chain !== atom.chain) continue;
        if (resi >= r.start && resi <= r.end) return hexToInt(d.color);
      }
    }
    return DOMAIN_FALLBACK_COLOR;
  }

  function ssColorFunc(atom) {
    return ssColor(atom.ss);
  }

  function spectrumColorFunc(atom) {
    const r = state.colorRange;
    return spectrumColor(atom.resi, r.resiMin, r.resiMax);
  }

  function bfactorColorFunc(atom) {
    const r = state.colorRange;
    return bfactorColor(atom.b, r.bMin, r.bMax);
  }

  function hydrophobicityColorFunc(atom) {
    return hydrophobicityColor(atom.resn);
  }

  function plddtColorFunc(atom) {
    return plddtColor(atom.b);
  }

  function conservationColorFunc(atom) {
    const cons = state.conservation;
    if (!cons || !cons.byResidue) return CONSERVATION_MISSING_COLOR;
    const c = cons.byResidue.get(`${atom.chain}:${atom.resi}`);
    if (c === undefined) return CONSERVATION_MISSING_COLOR;
    return conservationColor(c);
  }

  /**
   * Colour by clinical significance of the variant at this residue. Residues
   * with no variant get a muted neutral so the variants stand out.
   */
  function variantColorFunc(atom) {
    const map = state.variantByResidue;
    if (!map || !map.size) return VARIANT_BACKGROUND_COLOR;
    const hit = map.get(`${atom.chain}:${atom.resi}`);
    return hit !== undefined ? hit : VARIANT_BACKGROUND_COLOR;
  }

  /** Index the placed variants by chain:resi, keeping the most severe. */
  function indexVariantsByResidue() {
    const map = new Map();
    for (const v of state.variants) {
      if (!v.inStructure) continue;
      const key = `${v.chain}:${v.resi}`;
      const prev = map.get(key);
      if (prev === undefined ||
          significanceRank(v.significance) < significanceRank(prev.significance)) {
        map.set(key, SIGNIFICANCE_COLORS[v.significance] || SIGNIFICANCE_COLORS.unknown);
      }
    }
    state.variantByResidue = map;
  }

  /** "#2f7ed8" -> 0x2f7ed8 */
  function hexToInt(hex) {
    if (typeof hex !== "string") return DOMAIN_FALLBACK_COLOR;
    const v = parseInt(hex.replace("#", ""), 16);
    return Number.isFinite(v) ? v : DOMAIN_FALLBACK_COLOR;
  }

  // Resolve a colour-scheme name to a 3Dmol colour spec.
  //
  // 3Dmol has two different fields here and they are easy to confuse:
  //
  //   style.color      a LITERAL colour — hex, number, or html colour name
  //   style.colorscheme  the NAME of a colour scheme
  //
  // A scheme name put in `color` is looked up as a literal colour, fails, and
  // falls back to 0x000000 — the structure renders black. Every scheme below
  // therefore goes out either as `colorscheme` (3Dmol's own schemes) or as
  // `colorfunc` (ours, so the colours match the annotation key and work for
  // every representation).
  function resolveColor(colorName, opacity) {
    switch (colorName) {
      case "chain":
        return { colorfunc: chainColorFunc, opacity };
      case "byDomain":
        return { colorfunc: domainColorFunc, opacity };
      case "ss":
        return { colorfunc: ssColorFunc, opacity };
      case "spectrum":
        return { colorfunc: spectrumColorFunc, opacity };
      case "bFactor":
        return { colorfunc: bfactorColorFunc, opacity };
      case "hydrophobicity":
        return { colorfunc: hydrophobicityColorFunc, opacity };
      case "pLDDT":
        return { colorfunc: plddtColorFunc, opacity };
      case "conservation":
        return { colorfunc: conservationColorFunc, opacity };
      case "variant":
        return { colorfunc: variantColorFunc, opacity };
      case "byElement":
        return { colorscheme: "Jmol", opacity };
      case "residue":
        return { colorscheme: "amino", opacity };
      case "uniform":
        // Visible against whichever background is selected.
        return { color: state.background === "black" ? "#e8e8e8" : "#7f8c9b", opacity };
      default:
        return { colorscheme: colorName, opacity };
    }
  }

  function getStyleSpec(styleName, colorName, opacity) {
    switch (styleName) {
      case "cartoon":
        return { cartoon: resolveColor(colorName, opacity) };
      case "surface":
        return { surface: { ...resolveColor(colorName, opacity), opacity: opacity * 0.85 } };
      case "ball+stick":
        return { stick: { ...resolveColor(colorName, opacity), radius: 0.15 } };
      case "stick":
        return { stick: { ...resolveColor(colorName, opacity), radius: 0.18 } };
      case "spacefill":
        return { sphere: { ...resolveColor(colorName, opacity), scale: 0.3 } };
      case "ribbon":
        return { ribbon: resolveColor(colorName, opacity) };
      case "line":
        return { line: resolveColor(colorName, opacity) };
      default:
        return { cartoon: resolveColor(colorName, opacity) };
    }
  }

  function setModelVisibility(idx) {
    state.models.forEach((m, i) => {
      const visible = i === idx || !TRAJECTORIES[state.currentProtein?.id];
      m.setStyle({}, visible ? getStyleSpec(state.currentStyle, state.currentColor, state.currentOpacity) : {});
      if (visible && i === idx) applyStyleToModel(m, i);
    });
    state.viewer.render();
  }

  // ====== CONTROLS ======
  function onStyleChange() {
    applyAllStyles();
    if (state.currentProtein && TRAJECTORIES[state.currentProtein.id]) {
      updateTrajectoryBlend(state.trajectoryValue);
    } else {
      state.viewer.render();
    }
    updateUrlDebounced();
  }
  dom.ctrlStyle.addEventListener("change", (e) => {
    state.currentStyle = e.target.value;
    onStyleChange();
  });
  dom.ctrlColor.addEventListener("change", (e) => {
    state.currentColor = e.target.value;
    onStyleChange();
  });
  dom.ctrlOpacity.addEventListener("input", (e) => {
    state.currentOpacity = parseFloat(e.target.value);
    onStyleChange();
  });
  dom.ctrlBg.addEventListener("change", (e) => {
    state.background = e.target.value;
    const bg = e.target.value === "white"
      ? "white"
      : e.target.value === "black"
      ? "black"
      : "0x202428";
    state.viewer.setBackgroundColor(bg);
    state.viewer.render();
    updateUrlDebounced();
  });

  // ====== HIGHLIGHT BUTTONS ======
  function updateHighlightButtons() {
    const protein = state.currentProtein;
    if (!protein) return;
    const h = protein.highlights || {};
    dom.btnHighlightHeme.hidden = !h.ligand;
    dom.btnHighlightActive.hidden = !h.residues;
    dom.btnHighlightInterface.hidden = state.currentProtein.id !== "hiv-protease";
  }
  dom.btnHighlightHeme.addEventListener("click", () => {
    const h = state.currentProtein.highlights;
    if (!h || !h.ligand) return;
    state.viewer.setStyle({ resn: h.ligand.resn }, {
      stick: { color: "magenta", radius: 0.2 },
      sphere: { color: "magenta", scale: 0.35 },
    });
    state.viewer.zoomTo({ resn: h.ligand.resn });
    state.viewer.render();
  });
  dom.btnHighlightActive.addEventListener("click", () => {
    const h = state.currentProtein.highlights;
    if (!h || !h.residues) return;
    state.viewer.removeAllShapes();
    h.residues.forEach((r) => {
      const sel = r.chain ? { chain: r.chain, resi: r.resi } : { resi: r.resi };
      state.viewer.addStyle(sel, {
        stick: { color: "yellow", radius: 0.25 },
        sphere: { color: "yellow", scale: 0.4 },
      });
    });
    state.viewer.zoomTo({ resi: h.residues.map((r) => r.resi) });
    state.viewer.render();
  });
  dom.btnHighlightInterface.addEventListener("click", () => {
    if (state.currentProtein.id !== "hiv-protease") return;
    state.viewer.setStyle({ chain: "A" }, { cartoon: { color: "cyan", opacity: 0.7 } });
    state.viewer.setStyle({ chain: "B" }, { cartoon: { color: "magenta", opacity: 0.7 } });
    state.viewer.render();
  });
  dom.btnHighlightAll.addEventListener("click", () => {
    // `colorscheme`, not `color`: 3Dmol reads `color` as a literal colour
    // value, so a scheme name there renders the structure black.
    state.viewer.setStyle({}, { stick: { radius: 0.12, colorscheme: "Jmol", opacity: 0.6 } });
    state.viewer.render();
  });
  dom.btnClearHighlights.addEventListener("click", () => {
    state.viewer.removeAllShapes();
    state.viewer.removeAllLabels();
    applyAllStyles();
    state.viewer.zoomTo();
    state.viewer.render();
  });

  // ====== VIEWER TOOLBAR ======
  dom.btnResetView.addEventListener("click", () => {
    state.viewer.zoomTo();
    state.viewer.render();
  });
  dom.btnSpin.addEventListener("click", () => {
    state.spinning = !state.spinning;
    dom.btnSpin.textContent = state.spinning ? "⏸ Pause" : "↻ Spin";
    if (state.spinning) state.viewer.spin("y", 0.6);
    else state.viewer.spin(false);
  });
  dom.btnScreenshot.addEventListener("click", () => {
    const dataURL = state.viewer.pngURI();
    const link = document.createElement("a");
    link.href = dataURL;
    const proteinName = state.currentProtein ? state.currentProtein.name.replace(/\s+/g, "-").toLowerCase() : "protein";
    link.download = `${proteinName}-${Date.now()}.png`;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  });
  dom.btnShare.addEventListener("click", () => shareCurrentView());

  // ====== TABS ======
  dom.tabBtns.forEach((btn) => {
    btn.addEventListener("click", () => {
      const tab = btn.dataset.tab;
      dom.tabBtns.forEach((b) => b.classList.toggle("active", b === btn));
      dom.tabPanels.forEach((p) => p.classList.toggle("active", p.dataset.tabPanel === tab));
      if (tab !== "tour") stopTour();
      if (tab !== "trajectory") stopTrajectory();
      if (tab !== "tools") setPickMode(false);
      if (tab === "presentation") refreshPresentationIfReady();
    });
  });

  // ====== SEQUENCE PANEL ======
  // Map amino-acid → color (similar to ClustalX-like scheme)
  const AA_COLORS = {
    A: "#80a0f0", R: "#f01505", N: "#00ff00", D: "#c048c0",
    C: "#f08000", Q: "#00ff00", E: "#c048c0", G: "#f09048",
    H: "#15a4a4", I: "#80a0f0", L: "#80a0f0", K: "#f01505",
    M: "#80a0f0", F: "#80a0f0", P: "#ffff00", S: "#00ff00",
    T: "#00ff00", W: "#80a0f0", Y: "#15a4a4", V: "#80a0f0",
  };

  function parseAndRenderSequence() {
    if (!state._lastPdbText || !state.models.length) {
      dom.sequenceInfo.textContent = "—";
      dom.sequenceChains.innerHTML = "";
      dom.sequenceSecondary.innerHTML = "";
      state.sequence = null;
      return;
    }
    state.sequence = parseSequenceFromPDB(state._lastPdbText, state._lastPdbFormat);
    renderSequencePanel();
  }

  function renderSequencePanel() {
    if (!state.sequence) {
      dom.sequenceChains.innerHTML = "";
      dom.sequenceInfo.textContent = "—";
      return;
    }
    const seq = state.sequence;
    const totalRes = seq.chains.reduce((acc, c) => acc + c.residues.length, 0);
    dom.sequenceInfo.textContent = `${seq.chains.length} chain${seq.chains.length !== 1 ? "s" : ""} · ${totalRes} residues`;
    dom.sequenceChains.innerHTML = "";
    seq.chains.forEach((chain, ci) => {
      const row = document.createElement("div");
      row.className = "sequence-row";
      const label = document.createElement("div");
      label.className = "sequence-chain-label";
      label.textContent = chain.id;
      label.title = `Chain ${chain.id} — ${chain.residues.length} residues`;
      row.appendChild(label);
      chain.residues.forEach((res) => {
        const span = document.createElement("span");
        span.className = "sequence-residue";
        const color = AA_COLORS[res.one] || "#cccccc";
        span.style.background = color;
        if (res.ss === "H") span.classList.add("ss-helix");
        else if (res.ss === "E") span.classList.add("ss-sheet");
        else span.classList.add("ss-coil");
        span.textContent = res.one;
        span.title = `${chain.id}:${res.resn}${res.resi} (${res.one})${res.ss === "H" ? " · helix" : res.ss === "E" ? " · strand" : ""}`;
        span.addEventListener("click", () => focusResidue(chain.id, res.resi));
        row.appendChild(span);
      });
      dom.sequenceChains.appendChild(row);
    });

    renderSecondaryStructure();
  }

  function renderSecondaryStructure() {
    if (!state.sequence) return;
    if (!state.sequenceShowSS) {
      dom.sequenceSecondary.hidden = true;
      return;
    }
    dom.sequenceSecondary.hidden = false;
    dom.sequenceSecondary.innerHTML = "";
    state.sequence.chains.forEach((chain) => {
      const wrap = document.createElement("div");
      wrap.style.marginBottom = "4px";
      const lbl = document.createElement("strong");
      lbl.textContent = `${chain.id}:`;
      lbl.style.color = "var(--text-1)";
      lbl.style.marginRight = "6px";
      wrap.appendChild(lbl);
      // Coalesce SS runs
      let i = 0;
      while (i < chain.residues.length) {
        const r = chain.residues[i];
        if (!r.ss) { i++; continue; }
        const start = r.resi;
        const ssType = r.ss;
        let j = i;
        while (j < chain.residues.length && chain.residues[j].ss === ssType) j++;
        const width = (j - i) * 14;
        const bar = document.createElement("span");
        bar.className = "sequence-secondary-bar " + (ssType === "H" ? "helix" : "sheet");
        bar.style.width = width + "px";
        bar.style.background = ssType === "H" ? "#ef5350" : "#fbbf24";
        bar.title = `${ssType === "H" ? "Helix" : "Strand"} ${chain.id}:${start}–${chain.residues[j - 1].resi} (${j - i} residues)`;
        wrap.appendChild(bar);
        i = j;
      }
      dom.sequenceSecondary.appendChild(wrap);
    });
  }

  function focusResidue(chainId, resi) {
    if (!state.models.length) return;
    const sel = { chain: chainId, resi: resi };
    // Highlight residue
    applyAllStyles();
    state.viewer.addStyle(sel, {
      stick: { color: "yellow", radius: 0.3 },
      sphere: { color: "yellow", scale: 0.5 },
    });
    state.viewer.zoomTo(sel, 0, false);
    state.viewer.render();
    // Mark active in sequence
    document.querySelectorAll(".sequence-residue").forEach((el) => {
      const key = `${chainId}:${resi}`;
      el.classList.toggle("active", el.title.startsWith(key) || el.title.includes(" " + key));
    });
  }

  dom.btnToggleSecondary.addEventListener("click", () => {
    state.sequenceShowSS = !state.sequenceShowSS;
    dom.btnToggleSecondary.textContent = state.sequenceShowSS ? "Hide SS" : "Show SS";
    renderSecondaryStructure();
  });

  // ====== MEASUREMENT TOOL ======
  function setPickMode(on) {
    state.pickMode = on;
    dom.pickHint.hidden = !on;
    dom.btnPickAtom.textContent = on ? "✕ Cancel picking" : "📍 Pick atoms";
    if (on) {
      state.pickBuffer = [];
      updatePickHint();
      state.viewer.setClickable(true, true);
    } else {
      state.pickBuffer = [];
      updatePickHint();
    }
  }

  function updatePickHint() {
    const n = state.pickBuffer.length;
    let label = "Click an atom to pick";
    if (n === 1) label = "Pick #2 — distance will be shown";
    else if (n === 2) label = "Pick #3 — angle will be shown";
    else if (n === 3) label = "Pick #4 — dihedral will be shown";
    dom.pickHintText.textContent = label;
  }

  dom.btnPickAtom.addEventListener("click", () => setPickMode(!state.pickMode));
  dom.btnCancelPick.addEventListener("click", () => setPickMode(false));

  function handleAtomPicked(atom) {
    if (!atom || !atom.elem) return;
    const summary = {
      serial: atom.serial,
      x: atom.x, y: atom.y, z: atom.z,
      chain: atom.chain || "",
      resi: atom.resi,
      resn: atom.resn,
      atom: atom.atom,
      elem: atom.elem,
      model: atom.model,
    };
    state.pickBuffer.push(summary);

    // Visual marker for picked atom
    state.viewer.addSphere({
      center: { x: atom.x, y: atom.y, z: atom.z },
      radius: 0.4,
      color: "yellow",
    });

    if (state.pickBuffer.length === 2) {
      finalizeMeasurement("distance", computeDistance(state.pickBuffer[0], state.pickBuffer[1]));
    } else if (state.pickBuffer.length === 3) {
      finalizeMeasurement("angle", computeAngle(state.pickBuffer[0], state.pickBuffer[1], state.pickBuffer[2]));
    } else if (state.pickBuffer.length === 4) {
      finalizeMeasurement("dihedral", computeDihedral(state.pickBuffer[0], state.pickBuffer[1], state.pickBuffer[2], state.pickBuffer[3]));
      state.pickBuffer = [];
    }
    updatePickHint();
    updateMeasurementsUI();
  }

  function finalizeMeasurement(type, value) {
    const atoms = state.pickBuffer.slice(0, type === "distance" ? 2 : type === "angle" ? 3 : 4);
    const label = atoms.map((a) => `${a.chain || ""}:${a.resn}${a.resi}.${a.atom}`).join(" – ");
    state.measurements.unshift({ type, atoms: atoms.slice(), value, label, ts: Date.now() });
    // Add label to viewer
    if (type === "distance") {
      state.viewer.addCylinder({
        start: atoms[0], end: atoms[1], radius: 0.1, color: "yellow", fromCap: 1, toCap: 1,
      });
      const mid = midPoint(atoms[0], atoms[1]);
      state.viewer.addLabel(`${value.toFixed(2)} Å`, {
        position: mid, backgroundColor: "rgba(0,0,0,0.7)", fontColor: "white", fontSize: 12, inFront: true,
      });
    } else if (type === "angle") {
      state.viewer.addCylinder({
        start: atoms[0], end: atoms[1], radius: 0.05, color: "yellow",
      });
      state.viewer.addCylinder({
        start: atoms[1], end: atoms[2], radius: 0.05, color: "yellow",
      });
      state.viewer.addLabel(`${value.toFixed(1)}°`, {
        position: atoms[1], backgroundColor: "rgba(0,0,0,0.7)", fontColor: "white", fontSize: 12, inFront: true, yOffset: 1.5,
      });
    } else if (type === "dihedral") {
      for (let i = 0; i < 3; i++) {
        state.viewer.addCylinder({
          start: atoms[i], end: atoms[i + 1], radius: 0.04, color: "cyan",
        });
      }
      const mid = midPoint(atoms[1], atoms[2]);
      state.viewer.addLabel(`τ ${value.toFixed(1)}°`, {
        position: mid, backgroundColor: "rgba(0,0,0,0.7)", fontColor: "white", fontSize: 12, inFront: true,
      });
    }
    state.viewer.render();
  }

  function updateMeasurementsUI() {
    if (!state.measurements.length) {
      dom.measurementTable.innerHTML = '<div class="measurement-table-empty">No measurements yet</div>';
      return;
    }
    dom.measurementTable.innerHTML = state.measurements.map((m, idx) => {
      const atoms = m.atoms.map((a) => `${a.atom}${a.resi}`).join(" · ");
      const val = m.type === "distance"
        ? `${m.value.toFixed(3)} Å`
        : `${m.value.toFixed(2)}°`;
      const typeLabel = { distance: "Distance", angle: "Angle", dihedral: "Dihedral" }[m.type];
      return `<div class="measurement-row">
        <span class="measurement-type">${typeLabel} #${state.measurements.length - idx}</span>
        <span class="measurement-value">${val}</span>
        <span class="measurement-atoms">${escapeHtml(atoms)}</span>
      </div>`;
    }).join("");
  }

  dom.btnClearMeasurements.addEventListener("click", () => clearMeasurements());
  function clearMeasurements() {
    state.measurements = [];
    state.pickBuffer = [];
    state.viewer.removeAllShapes();
    state.viewer.removeAllLabels();
    state.viewer.render();
    if (state.pickMode) updatePickHint();
    updateMeasurementsUI();
  }

  // ====== COMPARISON + RMSD ======
  dom.btnEnableCompare.addEventListener("click", () => {
    dom.compareControls.hidden = !dom.compareControls.hidden;
    if (!dom.compareControls.hidden) {
      dom.btnEnableCompare.textContent = "⇄ Hide compare input";
      dom.compareInput.focus();
    } else {
      dom.btnEnableCompare.textContent = "⇄ Load structure B";
    }
  });

  dom.btnLoadCompare.addEventListener("click", () => loadCompareProtein(dom.compareInput.value));
  dom.compareInput.addEventListener("keydown", (e) => {
    if (e.key === "Enter") {
      e.preventDefault();
      loadCompareProtein(dom.compareInput.value);
    }
  });

  dom.btnClearCompare.addEventListener("click", () => clearComparison());

  async function loadCompareProtein(pdbIdRaw) {
    const pdbId = String(pdbIdRaw || "").trim().toUpperCase();
    if (!/^[0-9][A-Z0-9]{3}$/.test(pdbId)) {
      showCustomError("Compare requires a valid 4-character PDB ID.");
      return;
    }
    showLoading(`Loading structure B (${pdbId})…`);
    try {
      const r = await net(`https://files.rcsb.org/download/${pdbId}.pdb`);
      if (!r.ok) throw new Error(`HTTP ${r.status}`);
      const data = await r.text();
      const model = state.viewer.addModel(data, "pdb");
      model.modelIndex = state.models.length;
      model.pdbId = pdbId;
      model.isCompare = true;
      // Color: grey + cartoon at 0.6 opacity
      state.viewer.setStyle({ model: model.modelIndex }, {
        cartoon: { color: "#9aa0a6", opacity: 0.55 },
      });
      state.viewer.zoomTo();
      state.viewer.render();
      state.compareProtein = { id: `compare-${pdbId}`, name: `Structure B (${pdbId})`, modelIndex: model.modelIndex, pdbId };
      state.compareLoaded = true;
      state.aligned = false;
      state.rmsdColorMode = false;
      hideLoading();
      showToast(`Loaded ${pdbId} as structure B. Click "Align on Cα".`, "success");
      updateRmsdSummary();
    } catch (err) {
      hideLoading();
      showCustomError(`Could not load ${pdbId}: ${err.message}`);
    }
  }

  function clearComparison() {
    if (!state.compareLoaded) return;
    const modelB = state.models.find((m) => m.isCompare);
    if (modelB) {
      state.viewer.removeModel(modelB);
      state.models = state.models.filter((m) => !m.isCompare);
    }
    state.compareProtein = null;
    state.compareLoaded = false;
    state.aligned = false;
    state.rmsdColorMode = false;
    state.rmsdPerResidue = null;
    applyAllStyles();
    state.viewer.zoomTo();
    state.viewer.render();
    updateRmsdSummary();
  }

  dom.btnAlignCalpha.addEventListener("click", () => alignAndComputeRmsd());
  dom.btnColorByRmsd.addEventListener("click", () => {
    if (!state.compareLoaded || !state.aligned) {
      showToast("Load structure B and align on Cα first", "");
      return;
    }
    state.rmsdColorMode = !state.rmsdColorMode;
    if (state.rmsdColorMode) {
      applyRmsdColors();
      showToast("Coloured by per-residue RMSD", "success");
    } else {
      applyAllStyles();
      showToast("Back to default colors", "");
    }
  });

  function alignAndComputeRmsd() {
    if (!state.compareLoaded) {
      showToast("Load structure B first", "");
      return;
    }
    const modelA = state.models.find((m) => !m.isCompare);
    const modelB = state.models.find((m) => m.isCompare);
    if (!modelA || !modelB) {
      showToast("Both structures must be loaded", "");
      return;
    }
    showLoading("Aligning on Cα…");
    // Align modelB onto modelA using 3Dmol's built-in superpose
    try {
      const ok = state.viewer.superpose(modelB, modelA, { atom: "CA" });
      if (!ok) throw new Error("superpose returned false");
      state.aligned = true;
      // Compute per-residue RMSD
      state.rmsdPerResidue = computePerResidueRmsd(modelA, modelB);
      const total = state.rmsdPerResidue.length
        ? Math.sqrt(state.rmsdPerResidue.reduce((acc, r) => acc + r.value * r.value, 0) / state.rmsdPerResidue.length)
        : 0;
      hideLoading();
      showToast(`Aligned. Cα RMSD = ${total.toFixed(2)} Å across ${state.rmsdPerResidue.length} residues.`, "success");
      updateRmsdSummary(total);
      state.viewer.render();
    } catch (err) {
      hideLoading();
      showCustomError(`Alignment failed: ${err.message}`);
    }
  }

  function applyRmsdColors() {
    if (!state.rmsdPerResidue) return;
    // Heatmap: green (0) -> yellow (3Å) -> red (10Å+)
    function colorFor(d) {
      const max = 10;
      const t = Math.min(1, d / max);
      // lerp green->yellow->red
      if (t < 0.5) {
        const u = t * 2;
        const r = Math.round(0x10 + (0xfb - 0x10) * u);
        const g = Math.round(0xb9 + (0xbf - 0xb9) * u);
        const b = Math.round(0x81 + (0x24 - 0x81) * u);
        return (r << 16) | (g << 8) | b;
      } else {
        const u = (t - 0.5) * 2;
        const r = Math.round(0xfb + (0xef - 0xfb) * u);
        const g = Math.round(0xbf + (0x53 - 0xbf) * u);
        const b = Math.round(0x24 + (0x50 - 0x24) * u);
        return (r << 16) | (g << 8) | b;
      }
    }
    const map = new Map(state.rmsdPerResidue.map((r) => [`${r.chain}:${r.resi}`, colorFor(r.value)]));
    state.viewer.setStyle({ model: 0, atom: "CA" }, {
      cartoon: {
        colorfunc: (atom) => {
          const c = map.get(`${atom.chain}:${atom.resi}`);
          if (c === undefined) return 0x9aa0a6;
          return c;
        },
      },
    });
    state.viewer.render();
  }

  function updateRmsdSummary(totalOverride) {
    if (!state.compareLoaded) {
      dom.rmsdSummary.innerHTML = '<div class="muted small">Load structure B to compare.</div>';
      return;
    }
    if (!state.rmsdPerResidue || !state.rmsdPerResidue.length) {
      dom.rmsdSummary.innerHTML = `<div class="muted small">Click <b>Align on Cα</b> to superimpose.</div>`;
      return;
    }
    const values = state.rmsdPerResidue.map((r) => r.value);
    const total = totalOverride != null
      ? totalOverride
      : Math.sqrt(values.reduce((a, v) => a + v * v, 0) / values.length);
    const max = Math.max(...values);
    const min = Math.min(...values);
    const divergent = state.rmsdPerResidue.filter((r) => r.value > 2.0)
      .slice(0, 5)
      .map((r) => `${r.chain}:${r.resi}=${r.value.toFixed(1)}`)
      .join(", ");
    dom.rmsdSummary.innerHTML = `
      <div class="rmsd-summary-row"><span>Cα RMSD (overall)</span><b>${total.toFixed(2)} Å</b></div>
      <div class="rmsd-summary-row"><span>Residues compared</span><b>${values.length}</b></div>
      <div class="rmsd-summary-row"><span>Range</span><b>${min.toFixed(2)} – ${max.toFixed(2)} Å</b></div>
      <div class="rmsd-bar" title="Min–Max range"></div>
      <div class="rmsd-summary-row"><span>Most divergent (>2 Å)</span><b>${escapeHtml(divergent) || "—"}</b></div>
    `;
  }

  // ====== BIOINFORMATICS: pLDDT, Conservation, Variants, Pockets, Ligand Contacts ======

  // pLDDT colouring: AlphaFold models store pLDDT (0–100) in the B-factor column.
  // The colour mapping itself lives in src/colors.js (plddtColor), so the
  // button and the colour dropdown produce identical results.
  function colorByPLDDT() {
    if (!state.models.length) {
      showCustomError("No structure loaded.");
      return;
    }
    if (!state.viewer) return;
    // Switch the colour scheme rather than stamping a style over the model.
    // Stamping forced cartoon on top of whatever representation was selected
    // and was undone by the next style change.
    state.currentColor = "pLDDT";
    dom.ctrlColor.value = "pLDDT";
    onStyleChange();
    showToast("Coloured by AlphaFold pLDDT confidence", "success");
  }
  dom.btnColorPlddt.addEventListener("click", colorByPLDDT);

  // ----- Conservation via simple JS MSA -----
  // Strategy:
  //   1. Look up the UniProt accession for the loaded PDB via RCSB's data API.
  //   2. Fetch the canonical sequence from UniProt.
  //   3. Fetch ~10 close homologs via UniProt's REST search ("family:<family>" or
  //      blast-like neighbours). Fallback: if no homologs, conservation is trivially 1.0.
  //   4. Run progressive pairwise alignment (Needleman-Wunsch) in the browser.
  //   5. Compute per-column conservation (fraction of identical residues at each position).
  //   6. Map back to PDB residue numbers via SIFTS and color the structure.
  // -------- Simple Needleman-Wunsch with BLOSUM62-like scoring --------

  // Progressive MSA: align the first sequence to each subsequent one using the running consensus.
  // Simpler & faster than tree-based MUSCLE; good enough for a handful of sequences.

  // Conservation score per column = fraction of sequences with a non-gap residue matching
  // the consensus. Returns array of values 0..1 aligned to the consensus.

  async function fetchRcsbUniprotMapping(pdbId) {
    // RCSB data API: fetch the UniProt accession(s) associated with a PDB entity
    // Returns: { "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers": [...] }
    try {
      const url = `https://data.rcsb.org/rest/v1/core/entry/${pdbId.toLowerCase()}`;
      const r = await net(url);
      if (!r.ok) return null;
      const j = await r.json();
      const entities = j.rcsb_entry_container_identifiers?.polymer_entity_ids || [];
      // Iterate entities to find a UniProt accession
      for (const eid of entities) {
        try {
          const er = await net(`https://data.rcsb.org/rest/v1/core/polymer_entity/${pdbId.toLowerCase()}/${eid}`);
          if (!er.ok) continue;
          const ej = await er.json();
          const refs = ej.rcsb_polymer_entity_container_identifiers?.reference_sequence_identifiers || [];
          for (const ref of refs) {
            if (ref.database_accession && ref.database_name === "UniProt") {
              return ref.database_accession;
            }
          }
        } catch (e) { /* skip */ }
      }
    } catch (e) { /* skip */ }
    return null;
  }

  // Build a PDB-resi → alignment-index map from the PDB ATOM records

  async function loadConservation() {
    if (state.conservationJob) return state.conservationJob;
    if (!state._lastPdbText || !state.currentProtein) {
      showCustomError("Load a protein first.");
      return;
    }
    const pdbId = state.currentProtein.pdbIds[state.currentPdbIndex];
    dom.conservationStatus.textContent = "Looking up UniProt accession…";
    dom.btnFetchConservation.disabled = true;

    const job = (async () => {
      try {
        const acc = await fetchRcsbUniprotMapping(pdbId);
        if (!acc) {
          dom.conservationStatus.textContent = "No UniProt mapping found for this PDB.";
          return;
        }
        dom.conservationStatus.textContent = `Fetching canonical sequence (${acc})…`;
        // UniProt canonical sequence
        const fastaUrl = `https://rest.uniprot.org/uniprotkb/${acc}.fasta`;
        const fr = await net(fastaUrl);
        if (!fr.ok) throw new Error("UniProt FASTA fetch failed");
        const fastaText = await fr.text();
        const canonical = parseFastaSequences(fastaText)[0] || "";
        if (!canonical) throw new Error("Empty canonical sequence");

        dom.conservationStatus.textContent = "Fetching homologs…";
        // Homologs via UniProt search: by organism (human) and protein name keywords
        const protName = state.currentProtein.name.split(/[·•:()-]/)[0].trim();
        const searchUrl = `https://rest.uniprot.org/uniprotkb/search?query=reviewed:true+AND+organism_id:9606+AND+name:${encodeURIComponent(protName)}&format=fasta&size=10`;
        const sr = await net(searchUrl);
        if (!sr.ok) throw new Error("Homolog search failed");
        const homologFasta = await sr.text();
        const homologs = parseFastaSequences(homologFasta);

        const allSeqs = [canonical, ...homologs].slice(0, 10);
        dom.conservationStatus.textContent = `Aligning ${allSeqs.length} sequences…`;
        const msa = progressiveMsa(allSeqs);
        const scores = conservationScores(msa);

        // Build PDB residue index to map consensus index → PDB residue
        const { map: residueMap } = buildResidueIndex(state._lastPdbText);
        const conservationByResidue = new Map();
        for (const [key, v] of residueMap) {
          if (v.seqIdx < scores.length) {
            conservationByResidue.set(key, scores[v.seqIdx]);
          }
        }
        state.conservation = { byResidue: conservationByResidue, homologs: allSeqs.length - 1 };

        // Colour the structure by conservation. Switch the colour scheme
        // rather than stamping a style, so the map applies to the current
        // representation and survives a later style change.
        state.currentColor = "conservation";
        dom.ctrlColor.value = "conservation";
        applyAllStyles();
        state.viewer.render();
        dom.conservationStatus.textContent = `Conservation: ${allSeqs.length - 1} homologs aligned · mean identity ${(scores.reduce((a, v) => a + v, 0) / scores.length * 100).toFixed(1)}%`;
        showToast("Conservation map applied", "success");
      } catch (err) {
        dom.conservationStatus.textContent = `MSA failed: ${err.message}`;
        showCustomError(`Conservation fetch failed: ${err.message}`);
      } finally {
        dom.btnFetchConservation.disabled = false;
        state.conservationJob = null;
      }
    })();
    state.conservationJob = job;
    return job;
  }
  dom.btnFetchConservation.addEventListener("click", loadConservation);
  dom.btnClearConservation.addEventListener("click", () => {
    state.conservation = null;
    // Conservation is no longer available, so fall back to a scheme that is.
    state.currentColor = "chain";
    dom.ctrlColor.value = "chain";
    applyAllStyles();
    state.viewer.render();
    dom.conservationStatus.textContent = "No conservation loaded.";
  });

  // Minimal FASTA parser

  // ----- ClinVar variants via NCBI E-utilities -----
  //
  // The pipeline:
  //   1. The gene symbol comes from the structure's UniProt accession, not a
  //      hard-coded map, so this works for every curated entry and for
  //      custom loads that resolve to UniProt.
  //   2. ClinVar ESearch finds the records; ESummary is called in batches of
  //      40 ids per request rather than one request per record.
  //   3. Each record's `protein_change` gives the wild type and the position
  //      in canonical UniProt numbering ("N1095S" -> Asn at 1095).
  //   4. That position is mapped onto a PDB residue number by aligning the
  //      structure's chain sequence to the UniProt sequence, constrained to
  //      the residue range the construct actually covers (SIFTS). The wild
  //      type is then checked against the structure: if it disagrees, the
  //      mapping for that residue is rejected.
  //   5. gnomAD gives the population allele frequency. That is what makes
  //      "pathogenic" mean something -- a pathogenic variant at 3% allele
  //      frequency in the population is noise, not a candidate.
  //
  // Significance is read from `germline_classification`. ClinVar used to
  // publish `clinical_significance`; it no longer does, and reading the old
  // field silently classified every record as "VUS".

  const CLINVAR_SIGNIFICANT_MAX = 1000; // clinically significant records to pull
  const CLINVAR_GENERAL_MAX = 50;      // additional records, for context
  const CLINVAR_BATCH = 150;           // ids per ESummary request
  const CLINVAR_PACING_MS = 400;        // stay under NCBI's 3 requests/second
  const GNOMAD_CHUNK = 20;             // variants per gnomAD request
  const GNOMAD_MAX_LOOKUPS = 60;       // bound the work for large genes

  /** UniProt accession -> gene symbol + canonical sequence. Cached per accession. */
  async function fetchUniProtInfo(accession) {
    if (!accession) return null;
    if (state._uniprotInfo && state._uniprotInfo.accession === accession) {
      return state._uniprotInfo;
    }
    try {
      const r = await net(`https://rest.uniprot.org/uniprotkb/${accession}?fields=gene_names,sequence`);
      if (!r.ok) return null;
      const j = await r.json();
      const info = {
        accession,
        gene: (j.genes && j.genes[0] && j.genes[0].geneName && j.genes[0].geneName.value) || null,
        sequence: (j.sequence && j.sequence.value) || "",
      };
      state._uniprotInfo = info;
      return info;
    } catch (e) {
      return null;
    }
  }

  /** Residue list of one chain, as parsed from the loaded PDB. */
  function chainResidueList(chainId) {
    if (!state.sequence) return [];
    const ch = state.sequence.chains.find((c) => c.id === chainId);
    if (!ch) return [];
    return ch.residues.map((r) => ({ resi: r.resi, resn: r.resn, one: r.one }));
  }

  /**
   * Map canonical UniProt positions onto PDB residue numbers for every chain
   * the construct covers.
   * @returns {Map<number, {chain: string, resi: number, resn: string}>}
   */
  function buildVariantResidueMap(upSeq, sifts) {
    const map = new Map();
    if (!upSeq || !state.sequence) return map;
    const alignCache = new Map();
    for (const s of sifts || []) {
      for (const chainId of s.chains || []) {
        const residues = chainResidueList(chainId);
        if (!residues.length) continue;
        const chainSeq = residues.map((r) => r.one || "X").join("");
        if (!chainSeq) continue;
        // Aligning the same chain sequence twice (homodimer) is wasted work.
        let aln = alignCache.get(chainSeq);
        if (!aln) {
          aln = nwAlign(upSeq, chainSeq);
          alignCache.set(chainSeq, aln);
        }
        const local = buildUniprotMap(aln.a, aln.b, residues, s.minPos, s.maxPos);
        for (const [pos, v] of local) {
          if (!map.has(pos)) map.set(pos, { chain: chainId, resi: v.resi, resn: v.resn });
        }
      }
    }
    return map;
  }

  /** Turn one ClinVar ESummary record into variant records. */
  function parseClinVarRecord(rec) {
    // ClinVar publishes three independent classification blocks: the germline
    // significance, the somatic clinical-impact tier, and oncogenicity.
    // germline_classification is the current field; clinical_significance is
    // the historical one, kept as a fallback for older cached responses.
    const gc = rec.germline_classification || rec.clinical_significance || {};
    const ic = rec.clinical_impact_classification || {};
    const oc = rec.oncogenicity_classification || {};
    const significance = classifySignificance(gc.description || "");
    const review = gc.review_status || "";
    const tier = classifyImpactTier(ic.description || "");
    const oncogenicity = classifyOncogenicity(oc.description || "");
    const vs = (rec.variation_set || [])[0] || {};
    const gnomadId = gnomadVariantId(vs);
    // The variation name refers to the canonical (MANE) transcript, so its
    // residue numbering is UniProt's. The `protein_change` field lists the
    // same variant under every transcript ClinVar knows about, each with its
    // own numbering, so it is only a fallback.
    let changes = extractProteinChange(vs.variation_name);
    if (!changes.length) changes = parseProteinChange(rec.protein_change);
    return {
      id: rec.accession_version || rec.uid || "",
      title: rec.title || "",
      trait: gc.description || "",
      significance,
      review,
      tier,
      oncogenicity,
      gnomadId,
      changes,
    };
  }

  /** Ask gnomAD for the population allele frequency of each variant. */
  async function fetchGnomadFrequencies(variants) {
    state.variantRateLimited = false;
    state.variantLookupErrors = 0;
    const byId = new Map();
    // Take the most significant first, but round-robin across the classes so
    // that a bounded number of lookups still covers drug-response and the
    // other clinically interesting classes, not just the pathogenic ones.
    const byClass = new Map();
    for (const v of variants) {
      if (!v.gnomadId) continue;
      const k = v.significance;
      if (!byClass.has(k)) byClass.set(k, []);
      byClass.get(k).push(v);
    }
    // Within a class, the variants worth a frequency lookup are the ones a
    // clinician would name: those reviewed by an expert panel, those with a
    // Tier I somatic impact, and then by position.
    const tierRank = (v) => (v.tier && /^Tier I/i.test(v.tier) ? 0 : v.tier ? 1 : 2);
    const reviewRank = (v) =>
      (/expert panel/i.test(v.review || "") ? 0 : /practice guideline/i.test(v.review || "") ? 1 : 2);
    for (const list of byClass.values()) {
      list.sort((a, b) =>
        tierRank(a) - tierRank(b) ||
        reviewRank(a) - reviewRank(b) ||
        a.pos - b.pos);
    }
    const classes = [...byClass.keys()].sort(
      (a, b) => significanceRank(a) - significanceRank(b));
    const queues = classes.map((k) => byClass.get(k).slice());
    const ordered = [];
    let more = true;
    while (more) {
      more = false;
      for (const q of queues) {
        const v = q.shift();
        if (v) { ordered.push(v); more = true; }
      }
    }
    for (const v of ordered) {
      if (!byId.has(v.gnomadId)) byId.set(v.gnomadId, v);
    }
    const ids = [...byId.keys()].slice(0, GNOMAD_MAX_LOOKUPS);
    for (let i = 0; i < ids.length; i += GNOMAD_CHUNK) {
      await gnomadChunk(ids.slice(i, i + GNOMAD_CHUNK), byId);
    }
  }

  const GNOMAD_URL = "https://gnomad.broadinstitute.org/api";
  const GNOMAD_CONCURRENCY = 3;    // parallel single lookups on the slow path
  const GNOMAD_PAUSE_MS = 300;     // gnomAD's rate limit is generous but finite

  /**
   * One gnomAD query for one variant id.
   *
   * The three outcomes are kept apart on purpose. gnomAD answers a
   * rate-limited request with HTTP 400 and the message "Query rate limit
   * exceeded" -- not 429 -- so a naive client reports "not in gnomAD" for
   * every variant, which is a different and much less useful statement.
   *
   * @returns {{af: number, source: string}|{absent: true}|{rateLimited: true}|null}
   */
  async function gnomadOne(id) {
    const query = `{ variant(variantId: ${JSON.stringify(id)}, dataset: gnomad_r4) ` +
      `{ variant_id exome { ac an af } genome { ac an af } } }`;
    try {
      const r = await net(GNOMAD_URL, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ query }),
      });
      // The body has to be read even on a non-2xx: gnomAD reports its rate
      // limit with HTTP 400 and the message "Query rate limit exceeded", and
      // that message is the only way to tell "we were throttled" apart from
      // "this variant is not in gnomAD". The validator makes that distinction
      // with a stable code rather than a message pattern.
      let j = null;
      try { j = await r.json(); } catch (e) { /* not json */ }
      const check = validateGnomadResponse(j);
      if (check.ok) {
        const v = j.data.variant;
        if (v.genome && typeof v.genome.af === "number") return { af: v.genome.af, source: "genome" };
        if (v.exome && typeof v.exome.af === "number") return { af: v.exome.af, source: "exome" };
        return { absent: true };
      }
      if (check.code === "rate_limited") return { rateLimited: true };
      if (check.code === "not_found") return { absent: true };
      return { error: true };
    } catch (e) {
      const msg = e && e.message ? e.message : "";
      if (/429|rate limit/i.test(msg)) return { rateLimited: true };
      return { error: true };
    }
  }

  /**
   * Look up a batch of variant ids.
   *
   * gnomAD answers "Variant not found" as an error when *any* id in the batch
   * is absent, which for a cancer-gene list is most of them -- the well
   * known drivers are somatic mutations and gnomAD is a population reference.
   * So a failed batch is retried one id at a time, a few in parallel, and an
   * id that is still missing is simply reported as "not in gnomAD".
   */
  async function gnomadChunk(ids, byId) {
    if (!ids.length) return;
    const parts = ids.map((id, k) =>
      `v${k}: variant(variantId: ${JSON.stringify(id)}, dataset: gnomad_r4) ` +
      `{ variant_id exome { ac an af } genome { ac an af } }`);
    const query = `{ ${parts.join(" ")} }`;
    let data = null;
    try {
      const r = await net(GNOMAD_URL, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ query }),
      });
      if (!r.ok) throw new Error(`HTTP ${r.status}`);
      const j = await r.json();
      if (j.errors && j.errors.length) throw new Error(j.errors[0].message || "gnomAD error");
      data = j.data || {};
    } catch (err) {
      data = null;
    }
    if (data) {
      ids.forEach((id, k) => {
        const v = data["v" + k];
        const target = byId.get(id);
        if (!target) return;
        if (!v) { target.afAbsent = true; return; }
        if (v.genome && typeof v.genome.af === "number") {
          target.af = v.genome.af; target.afSource = "genome";
        } else if (v.exome && typeof v.exome.af === "number") {
          target.af = v.exome.af; target.afSource = "exome";
        } else {
          target.afAbsent = true;
        }
      });
      return;
    }
    // Slow path: one id at a time, a few in parallel.
    for (let i = 0; i < ids.length; i += GNOMAD_CONCURRENCY) {
      const group = ids.slice(i, i + GNOMAD_CONCURRENCY);
      const results = await Promise.all(group.map((id) => gnomadOne(id)));
      group.forEach((id, k) => {
        const hit = results[k];
        const target = byId.get(id);
        if (!hit || !target) return;
        if (hit.rateLimited) { state.variantRateLimited = true; return; }
        if (hit.error) { state.variantLookupErrors = (state.variantLookupErrors || 0) + 1; return; }
        if (hit.absent) { target.afAbsent = true; return; }
        target.af = hit.af;
        target.afSource = hit.source;
      });
      await new Promise((r) => setTimeout(r, GNOMAD_PAUSE_MS));
    }
  }

  async function fetchVariants() {
    if (!state.currentProtein) {
      showCustomError("Load a protein first.");
      return;
    }
    dom.variantsStatus.textContent = "Looking up UniProt accession…";
    dom.btnFetchVariants.disabled = true;

    try {
      const pdbId = state.currentProtein.pdbIds[state.currentPdbIndex];

      // The gene symbol and the canonical sequence both come from UniProt.
      const upInfo = await fetchUniProtInfo(state.meta && state.meta.uniprot);
      if (!upInfo || !upInfo.gene) {
        dom.variantsStatus.textContent =
          "No UniProt accession for this structure, so the gene symbol could not be determined. " +
          "Custom-loaded structures need a UniProt accession (we'll auto-detect next iteration).";
        return;
      }
      const gene = upInfo.gene;

      // ClinVar ESearch. Several queries, because ClinVar sorts by record id
      // descending -- the newest records first -- and the well-established
      // ones are often decades old. A single capped query would miss
      // EGFR T790M entirely. Multi-word field values must be quoted or
      // they match nothing at all.
      dom.variantsStatus.textContent = `Querying ClinVar for ${gene}…`;
      const esearch = async (term, max) => {
        const u = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar` +
          `&term=${encodeURIComponent(term)}&retmax=${max}&retmode=json`;
        const res = await net(u);
        if (!res.ok) return [];
        const jj = await res.json();
        return (jj.esearchresult && jj.esearchresult.idlist) || [];
      };
      const significantTerms = [
        "pathogenic[Clinical significance]",
        "likely pathogenic[Clinical significance]",
        "conflicting[Clinical significance]",
        '"drug response"[Clinical significance]',
        '"risk factor"[Clinical significance]',
      ];
      const significantIds = await esearch(
        `${gene}[gene] AND (${significantTerms.join(" OR ")})`,
        CLINVAR_SIGNIFICANT_MAX);
      // NCBI E-utilities allow three requests a second without an API key.
      await new Promise((r) => setTimeout(r, CLINVAR_PACING_MS));
      const generalIds = await esearch(`${gene}[gene]`, CLINVAR_GENERAL_MAX);
      const ids = [...new Set([...significantIds, ...generalIds])];
      if (ids.length === 0) {
        dom.variantsStatus.textContent = `No ClinVar records found for ${gene}.`;
        return;
      }

      // ClinVar ESummary, batched. This is the slow part -- a gene like EGFR
      // has several hundred significant records -- so report progress.
      dom.variantsStatus.textContent =
        `Found ${ids.length} ClinVar records for ${gene}. Summarising…`;
      const records = [];
      for (let i = 0; i < ids.length; i += CLINVAR_BATCH) {
        const batch = ids.slice(i, i + CLINVAR_BATCH).join(",");
        dom.variantsStatus.textContent =
          `Summarising ${Math.min(i + CLINVAR_BATCH, ids.length)}/${ids.length} ClinVar records…`;
        try {
          const sr = await net(
            `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id=${batch}&retmode=json`);
          if (!sr.ok) {
            console.warn(`ClinVar ESummary batch failed: HTTP ${sr.status}`);
          } else {
            const sj = await sr.json();
            const res = sj.result || {};
            for (const uid of res.uids || []) {
              if (!res[uid]) continue;
              // A record missing the fields the app reads would otherwise
              // become an undefined value deep inside the renderer.
              const check = validateClinVarRecord(res[uid]);
              if (!check.ok) {
                console.warn(`ClinVar record ${uid} skipped: ${check.error}`);
                continue;
              }
              records.push(parseClinVarRecord(res[uid]));
            }
          }
        } catch (e) {
          // Never swallow this silently: a skipped batch looks exactly like
          // "this gene has no variants", which is how a whole panel can come
          // back empty with no explanation.
          console.warn(`ClinVar ESummary batch failed: ${e.message}`);
        }
        await new Promise((r) => setTimeout(r, CLINVAR_PACING_MS));
      }
      if (!records.length) {
        dom.variantsStatus.textContent = "ClinVar returned no usable records.";
        return;
      }

      // Expand each record into one entry per protein change.
      const all = [];
      for (const rec of records) {
        for (const ch of rec.changes) {
          all.push({
            pos: ch.pos, wt: ch.wt, mut: ch.mut,
            significance: rec.significance,
            review: rec.review,
            tier: rec.tier,
            oncogenicity: rec.oncogenicity,
            trait: rec.trait,
            gnomadId: rec.gnomadId,
            id: rec.id,
            title: rec.title,
          });
        }
      }
      const merged = mergeByChange(all);

      // Place them on the structure.
      dom.variantsStatus.textContent =
        `Mapping ${merged.length} protein changes onto ${pdbId}…`;
      const resMap = buildVariantResidueMap(upInfo.sequence, state.meta && state.meta.sifts);
      let placed = 0, wrongWildType = 0;
      for (const v of merged) {
        const hit = resMap.get(v.pos);
        if (!hit) { v.inStructure = false; continue; }
        // The wild type is a free correctness check on the mapping.
        const structureAa = AA3TO1[hit.resn] || "X";
        if (v.wt !== "*" && structureAa !== v.wt) {
          v.inStructure = false;
          v.mappingConflict = `${v.wt} vs ${structureAa}`;
          wrongWildType++;
          continue;
        }
        v.inStructure = true;
        v.chain = hit.chain;
        v.resi = hit.resi;
        v.resn = hit.resn;
        v.mutation = `${v.wt}${v.pos}${v.mut}`;
        v.color = SIGNIFICANCE_COLORS[v.significance] || SIGNIFICANCE_COLORS.unknown;
        placed++;
      }

      // Population frequency.
      const toFreq = merged.filter((v) => v.inStructure && v.gnomadId);
      if (toFreq.length) {
        dom.variantsStatus.textContent =
          `Placed ${placed} of ${merged.length}. Fetching population frequencies…`;
        await fetchGnomadFrequencies(toFreq);
        for (const v of toFreq) {
          v.rarity = rarityLabel(v.af);
        }
      }

      state.variants = merged;
      indexVariantsByResidue();
      renderVariantPanel(gene, records.length, merged, placed, wrongWildType);
      renderVariantSpheres();
      // Show the significance colouring straight away, so the fetch does
      // something visible on the structure.
      if (placed > 0) {
        state.currentColor = "variant";
        dom.ctrlColor.value = "variant";
        applyAllStyles();
      }
      state.viewer.render();
      showToast(`Loaded ${merged.length} protein changes for ${gene}`, "success");
    } catch (err) {
      dom.variantsStatus.textContent = `ClinVar fetch failed: ${err.message}`;
      showCustomError(`ClinVar fetch failed: ${err.message}`);
    } finally {
      dom.btnFetchVariants.disabled = false;
    }
  }

  dom.btnFetchVariants.addEventListener("click", fetchVariants);

  /** Write the variant table and the summary line into the panel. */
  function renderVariantPanel(gene, recordCount, variants, placed, wrongWildType) {
    const tally = {};
    for (const v of variants) tally[v.significance] = (tally[v.significance] || 0) + 1;
    const dot = {
      pathogenic: "dot-path", likely_pathogenic: "dot-likpath", vus: "dot-vus",
      conflicting: "dot-conflict", likely_benign: "dot-likben", benign: "dot-ben",
      drug_response: "dot-drugresp", risk_factor: "dot-risk",
      protective: "dot-ben", unknown: "dot-vus",
    };
    const order = ["pathogenic", "likely_pathogenic", "drug_response", "risk_factor",
                   "conflicting", "vus", "likely_benign", "benign", "protective"];
    const summary = order.filter((k) => tally[k])
      .map((k) => `<span class="dot ${dot[k]}"></span> ${tally[k]} ${k.replace("_", " ")}`)
      .join(" · ");

    const rows = variants.filter((v) => v.inStructure)
      .sort((a, b) => significanceRank(a.significance) - significanceRank(b.significance) || a.pos - b.pos);

    const table = rows.length ? `
      <div class="variant-table" id="variant-table">
        <div class="variant-table-head">
          <span>Change</span><span>Residue</span><span>Significance</span><span>Frequency</span>
        </div>
        ${rows.map((v) => {
          const extra = [v.tier, v.oncogenicity].filter(Boolean).join(" \u00b7 ");
          const review = v.review && /expert panel|practice guideline|multiple submitters/i.test(v.review)
            ? v.review.replace(/^reviewed by /, "") : "";
          const note = [review, extra].filter(Boolean).join(" \u00b7 ");
          return `
          <div class="variant-row" data-chain="${escapeHtml(v.chain)}" data-resi="${v.resi}">
            <span class="variant-mut"><span class="dot ${dot[v.significance] || "dot-vus"}"></span>${escapeHtml(v.mutation)}</span>
            <span class="variant-resi">${escapeHtml(v.chain)}:${v.resi} ${escapeHtml(v.resn)}</span>
            <span class="variant-sig">${escapeHtml(String(v.significance).replace(/_/g, " "))}${
              note ? `<em class="variant-note">${escapeHtml(note)}</em>` : ""}</span>
            <span class="variant-af">${v.af !== undefined && v.af !== null
              ? `${v.af.toExponential(1)} <em>${escapeHtml(v.rarity || "")}</em>`
              : v.afAbsent
              ? "<em class='muted'>not in gnomAD</em>"
              : "<em class='muted'>not looked up</em>"}</span>
          </div>`;
        }).join("")}
      </div>` : "";

    const notPlaced = variants.length - placed;
    dom.variantsStatus.innerHTML =
      `<b>${escapeHtml(gene)}</b> · ${variants.length} protein changes from ${recordCount} ClinVar records<br>` +
      `<span class="muted small">${summary}</span><br>` +
      `<span class="muted small">${placed} fall inside this structure` +
      (notPlaced ? `, ${notPlaced} lie outside the construct` : "") +
      (wrongWildType ? ` · ${wrongWildType} rejected: wild type did not match` : "") +
      `</span>` +
      (state.variantRateLimited
        ? `<br><span class="muted small">gnomAD rate-limited some frequency lookups; ` +
          `those rows show no frequency rather than a wrong one.</span>`
        : "") +
      (state.variantLookupErrors
        ? `<br><span class="muted small">${state.variantLookupErrors} gnomAD lookups failed; ` +
          `those rows show no frequency rather than a wrong one.</span>`
        : "") +
      `<br><span class="muted small">Population frequencies are looked up for the ` +
      `${GNOMAD_MAX_LOOKUPS} most significant changes only; gnomAD is a population ` +
      `reference, so a somatic driver is often absent from it.</span>` +
      table;

    // Clicking a row focuses that residue.
    dom.variantsStatus.querySelectorAll(".variant-row").forEach((row) => {
      row.addEventListener("click", () => {
        const chain = row.dataset.chain, resi = parseInt(row.dataset.resi, 10);
        const model = state.models[state.currentPdbIndex];
        if (!model) return;
        const atoms = model.selectedAtoms({ chain, resi });
        if (!atoms.length) return;
        state.viewer.zoomTo({ chain, resi });
        state.viewer.render();
      });
    });
  }

  /** Draw a sphere at every placed variant, coloured by significance. */
  function renderVariantSpheres() {
    state.viewer.removeAllShapes();
    const placed = state.variants.filter((v) => v.inStructure);
    if (!placed.length) return;
    const model = state.models[state.currentPdbIndex];
    if (!model) return;
    for (const v of placed) {
      const atoms = model.selectedAtoms({ chain: v.chain, resi: v.resi });
      if (!atoms.length) continue;
      const anchor = atoms.find((a) => a.atom === "CA") || atoms[0];
      state.viewer.addSphere({
        center: { x: anchor.x, y: anchor.y, z: anchor.z },
        radius: 1.6,
        color: v.color,
        opacity: 0.85,
      });
    }
  }

  dom.btnClearVariants.addEventListener("click", () => {
    state.variants = [];
    state.variantByResidue = new Map();
    dom.variantsStatus.textContent = "No variants loaded.";
    // Remove the variant spheres and fall back to a scheme that still works.
    state.viewer.removeAllShapes();
    if (state.currentColor === "variant") {
      state.currentColor = "chain";
      dom.ctrlColor.value = "chain";
    }
    applyAllStyles();
    state.viewer.render();
  });

  // ----- Tissue expression (GTEx) -----
  //
  // The third pillar of a biomarker claim. The variant panel says the change
  // is significant and rare; this says the protein it changes is actually
  // made in the tissue where the disease happens. GTEx reports median
  // transcripts per million across 54 human tissues.
  //
  // The portal is CORS-open, so this runs from the browser like every other
  // data source in the app.

  const GTEX_API = "https://gtexportal.org/api/v2";

  async function fetchExpression() {
    if (!state.currentProtein) {
      showCustomError("Load a protein first.");
      return;
    }
    dom.btnFetchExpression.disabled = true;
    dom.expressionStatus.textContent = "Looking up the gene…";
    try {
      // The gene symbol comes from the structure's UniProt accession, so this
      // works for every curated entry and for custom loads too.
      const up = await fetchUniProtInfo(state.meta && state.meta.uniprot);
      if (!up || !up.gene) {
        dom.expressionStatus.textContent =
          "No UniProt accession for this structure, so the gene symbol could not be determined.";
        return;
      }
      const gene = up.gene;

      dom.expressionStatus.textContent = `Looking up ${gene} in GTEx…`;
      const refRes = await net(`${GTEX_API}/reference/gene?geneId=${encodeURIComponent(gene)}`);
      if (!refRes.ok) throw new Error(`GTEx gene lookup returned ${refRes.status}`);
      const refJson = await refRes.json();
      const ref = (refJson.data || [])[0];
      if (!ref || !ref.gencodeId) throw new Error(`GTEx has no entry for ${gene}.`);

      dom.expressionStatus.textContent = `Fetching median expression for ${gene}…`;
      const exprRes = await net(
        `${GTEX_API}/expression/medianGeneExpression` +
        `?gencodeId=${encodeURIComponent(ref.gencodeId)}&datasetId=gtex_v8`);
      if (!exprRes.ok) throw new Error(`GTEx expression lookup returned ${exprRes.status}`);
      const rows = parseMedianExpression(await exprRes.json());
      if (!rows.length) throw new Error(`GTEx returned no expression data for ${gene}.`);

      state.expression = { gene, rows, summary: summarizeExpression(rows) };
      renderExpressionPanel();
      showToast(`Loaded GTEx tissue expression for ${gene}`, "success");
    } catch (err) {
      state.expression = null;
      dom.expressionStatus.textContent = `GTEx fetch failed: ${err.message}`;
      showCustomError(`GTEx fetch failed: ${err.message}`);
    } finally {
      dom.btnFetchExpression.disabled = false;
    }
  }

  /** Render the specificity verdict and the tissues, grouped by compartment. */
  function renderExpressionPanel() {
    const ex = state.expression;
    if (!ex) {
      dom.expressionStatus.textContent = "No expression loaded.";
      return;
    }
    const sum = ex.summary;
    const groups = groupByCompartment(ex.rows);
    const maxTpm = Math.max(sum.maxTpm || 0, 1);

    const rowsHtml = groups.map((g) => `
      <div class="expr-group">
        <div class="expr-group-head">${escapeHtml(g.label)}</div>
        ${g.tissues.map((t) => `
          <div class="expr-tissue${t.detected ? "" : " expr-undetected"}">
            <span class="expr-tissue-name">${escapeHtml(t.label)}</span>
            <span class="expr-bar"><span class="expr-bar-fill" style="width:${Math.max(2, Math.round((100 * t.tpm) / maxTpm))}%;background:#${t.color.toString(16).padStart(6, "0")}"></span></span>
            <span class="expr-tpm">${t.tpm >= 10 ? t.tpm.toFixed(0) : t.tpm.toFixed(2)}</span>
          </div>`).join("")}
      </div>`).join("");

    dom.expressionStatus.innerHTML =
      `<b>${escapeHtml(ex.gene)}</b> · ${escapeHtml(sum.specificityLabel)}` +
      ` · ${sum.detectedCount} of ${sum.tissueCount} tissues detected` +
      `<br><span class="muted small">` +
      `Peak: ${escapeHtml(sum.maxTissueLabel)} at ${sum.maxTpm.toFixed(1)} TPM` +
      (sum.foldVsMedian ? ` (${sum.foldVsMedian.toFixed(1)}× the median)` : "") +
      (sum.dominanceRatio && sum.dominanceRatio !== Infinity
        ? ` · ${sum.dominanceRatio.toFixed(1)}× the runner-up` : "") +
      ` · median ${sum.medianTpm.toFixed(2)} TPM · mean ${sum.meanTpm.toFixed(1)} TPM` +
      `</span>` +
      `<div class="expr-legend"><span class="muted small">` +
      `low <span class="expr-scale"></span> high (log scale, 0–1000 TPM)` +
      `</span></div>` +
      rowsHtml;
  }

  dom.btnFetchExpression.addEventListener("click", fetchExpression);
  dom.btnClearExpression.addEventListener("click", () => {
    state.expression = null;
    dom.expressionStatus.textContent = "No expression loaded.";
  });

  // ----- Pocket detection (alpha-sphere grid scan) -----
// ---- pocket detection -------------------------------------------------
  // Alpha-sphere scan in the spirit of Fpocket, but computed on a grid so it
  // stays interactive in the browser:
  //
  //   1. cheap density gates  - a candidate must sit in an actual void
  //      (no atoms within 3.2 A), have pocket walls (>=18 atoms within 6.4 A)
  //      and be enclosed by the protein (>=45 atoms within 10 A);
  //   2. ray-cast enclosure   - cast 32 rays over a Fibonacci sphere; the
  //      fraction that hit protein must clear a threshold. This is what
  //      separates a concave pocket from a convex surface groove, which the
  //      density gates alone cannot do;
  //   3. spatial-hash clustering at 3 A, ranked by volume.
  //
  // Validated on EGFR/erlotinib (1M17): the top-ranked pocket is ~100%
  // occupied by the crystallographic ligand. Mean enclosure is 0.78 for
  // points on the drug site versus 0.60 elsewhere.
  const POCKET_GRID = 1.4;
  const POCKET_VOID_R = 3.2;     // must be atom-free
  const POCKET_WALL_R = 6.4;     // pocket walls
  const POCKET_GLOBAL_R = 10.0;  // overall enclosure
  const POCKET_MIN_WALL = 18;
  const POCKET_MIN_GLOBAL = 45;
  const POCKET_MIN_ENCLOSURE = 0.80;
  const POCKET_CLUSTER_R = 3.0;

  // Fibonacci sphere unit vectors, computed once.

  // True if a ray from (x,y,z) along unit (dx,dy,dz) hits protein within radius.

  async function detectPockets() {
    if (!state.models.length) { showCustomError("No structure loaded."); return; }
    dom.btnDetectPockets.disabled = true;
    dom.pocketsStatus.textContent = "Scanning for cavities…";
    dom.pocketList.innerHTML = "";

    // Yield a frame so the "scanning" message paints before the scan starts.
    setTimeout(async () => {
      try {
        const model = state.models[state.currentPdbIndex];
        const atoms = model.selectedAtoms({ elem: { not: "H" } });
        if (!atoms || !atoms.length) {
          dom.pocketsStatus.textContent = "No heavy atoms in structure.";
          dom.btnDetectPockets.disabled = false;
          return;
        }

        let minX = Infinity, minY = Infinity, minZ = Infinity;
        let maxX = -Infinity, maxY = -Infinity, maxZ = -Infinity;
        for (const a of atoms) {
          if (a.x < minX) minX = a.x; if (a.x > maxX) maxX = a.x;
          if (a.y < minY) minY = a.y; if (a.y > maxY) maxY = a.y;
          if (a.z < minZ) minZ = a.z; if (a.z > maxZ) maxZ = a.z;
        }
        const PAD = 6;
        minX -= PAD; minY -= PAD; minZ -= PAD;
        maxX += PAD; maxY += PAD; maxZ += PAD;

        const hash = buildAtomHash(atoms);

        // Pass 1 — density gates. Runs in time-sliced chunks so a large structure
        // (e.g. a 23k-atom ribosome, ~8.7M grid points) cannot freeze the UI.
        // Synchronously this measured 1.55 s on EGFR and ~80 s on 6VXX.
        const cands = [];
        const SLICE_MS = 12;             // work per frame before yielding
        const nx = Math.ceil((maxX - minX) / POCKET_GRID);
        const ny = Math.ceil((maxY - minY) / POCKET_GRID);
        const totalPlanes = nx * ny;
        let sliceStart = performance.now();
        let planesDone = 0;
        for (let ix = 0; ix < nx; ix++) {
          for (let iy = 0; iy < ny; iy++) {
            const x = minX + ix * POCKET_GRID;
            const y = minY + iy * POCKET_GRID;
            for (let z = minZ; z <= maxZ; z += POCKET_GRID) {
              if (countAtomsNear(hash, x, y, z, POCKET_VOID_R) > 0) continue;
              const wall = countAtomsNear(hash, x, y, z, POCKET_WALL_R);
              if (wall < POCKET_MIN_WALL) continue;
              const global = countAtomsNear(hash, x, y, z, POCKET_GLOBAL_R);
              if (global < POCKET_MIN_GLOBAL) continue;
              cands.push({ x, y, z, wall, global, e: 0 });
            }
            planesDone++;
          }
          if (performance.now() - sliceStart > SLICE_MS) {
            dom.pocketsStatus.textContent =
              `Scanning… ${Math.round((100 * planesDone) / totalPlanes)}%`;
            await new Promise((r) => setTimeout(r, 0));
            sliceStart = performance.now();
          }
        }

        // Pass 2 — ray-cast enclosure (only for the cheap-gate survivors).
        // Also time-sliced: each point is 32 rays, so this dominates on
        // structures where many points survive the density gates.
        const keep = [];
        let rayStart = performance.now();
        for (let ci = 0; ci < cands.length; ci++) {
          const c = cands[ci];
          const e = enclosureFraction(hash, c.x, c.y, c.z);
          if (e >= POCKET_MIN_ENCLOSURE) { c.e = e; keep.push(c); }
          if (performance.now() - rayStart > SLICE_MS) {
            dom.pocketsStatus.textContent =
              `Testing concavity… ${Math.round((100 * ci) / cands.length)}%`;
            await new Promise((r) => setTimeout(r, 0));
            rayStart = performance.now();
          }
        }

        // Pass 3 — spatial-hash clustering (O(n), not the previous O(n^2) scan).
        const buckets = new Map();
        keep.forEach((p, i) => {
          const k = `${Math.floor(p.x / POCKET_CLUSTER_R)}|${Math.floor(p.y / POCKET_CLUSTER_R)}|${Math.floor(p.z / POCKET_CLUSTER_R)}`;
          let arr = buckets.get(k);
          if (!arr) { arr = []; buckets.set(k, arr); }
          arr.push(i);
        });
        const visited = new Uint8Array(keep.length);
        const pockets = [];
        for (let i = 0; i < keep.length; i++) {
          if (visited[i]) continue;
          visited[i] = 1;
          const queue = [i];
          const cluster = [];
          while (queue.length) {
            const idx = queue.pop();
            const p = keep[idx];
            cluster.push(p);
            const cx = Math.floor(p.x / POCKET_CLUSTER_R);
            const cy = Math.floor(p.y / POCKET_CLUSTER_R);
            const cz = Math.floor(p.z / POCKET_CLUSTER_R);
            for (let dx = -1; dx <= 1; dx++) {
              for (let dy = -1; dy <= 1; dy++) {
                for (let dz = -1; dz <= 1; dz++) {
                  const arr = buckets.get(`${cx + dx}|${cy + dy}|${cz + dz}`);
                  if (!arr) continue;
                  for (const j of arr) {
                    if (visited[j]) continue;
                    const r = keep[j];
                    const ddx = r.x - p.x, ddy = r.y - p.y, ddz = r.z - p.z;
                    if (ddx * ddx + ddy * ddy + ddz * ddz < POCKET_CLUSTER_R * POCKET_CLUSTER_R) {
                      visited[j] = 1;
                      queue.push(j);
                    }
                  }
                }
              }
            }
          }
          if (cluster.length < 2) continue;
          let cx = 0, cy = 0, cz = 0, meanE = 0;
          for (const p of cluster) { cx += p.x; cy += p.y; cz += p.z; meanE += p.e; }
          const n = cluster.length;
          pockets.push({
            center: [cx / n, cy / n, cz / n],
            volume: n * POCKET_GRID * POCKET_GRID * POCKET_GRID,
            enclosure: meanE / n,
            residues: [],
          });
        }

        // Rank by volume; if a ligand is loaded, promote pockets it occupies.
        const ligAtoms = (state.activeLigand && state.activeLigand.compId)
          ? safeSelect(model, { resn: state.activeLigand.compId })
          : [];
        for (const pk of pockets) {
          pk.drugOcc = 0;
          if (ligAtoms.length) {
            let hit = 0;
            for (const la of ligAtoms) {
              const d = computeDistance({ x: pk.center[0], y: pk.center[1], z: pk.center[2] }, la);
              if (d <= 5) { hit++; break; }
            }
            pk.drugOcc = hit;
          }
        }
        pockets.sort((a, b) => (b.drugOcc - a.drugOcc) || (b.volume - a.volume));
        const top = pockets.slice(0, 8);

        // Map each pocket to its lining residues.
        for (const pk of top) {
          const set = new Set();
          const r2 = POCKET_CLUSTER_R * POCKET_CLUSTER_R;
          for (const a of atoms) {
            const dx = a.x - pk.center[0], dy = a.y - pk.center[1], dz = a.z - pk.center[2];
            if (dx * dx + dy * dy + dz * dz < r2) set.add(`${a.chain}:${a.resi}`);
          }
          pk.residues = Array.from(set).slice(0, 25);
        }
        state.pockets = top;

        // Draw the pockets.
        state.viewer.removeAllShapes();
        top.forEach((pk, i) => {
          const colour = i === 0 ? 0xff7d45 : (i === 1 ? 0xffdb13 : 0x65cbf3);
          state.viewer.addSphere({
            center: { x: pk.center[0], y: pk.center[1], z: pk.center[2] },
            radius: Math.cbrt(pk.volume * 0.75 / Math.PI),
            color: colour,
            opacity: 0.45,
          });
        });
        for (const pk of top) {
          for (const r of pk.residues) {
            const parts = r.split(":");
            state.viewer.addStyle({ chain: parts[0], resi: parseInt(parts[1], 10) },
              { stick: { radius: 0.12, color: 0xff7d45 } });
          }
        }
        state.viewer.render();

        dom.pocketList.innerHTML = top.length
          ? top.map((pk, i) => `
            <div class="pocket-row">
              <span class="pocket-rank">#${i + 1}</span>
              <span class="pocket-vol">~${pk.volume.toFixed(0)} Å³</span>
              <span class="pocket-name">encl ${(pk.enclosure * 100).toFixed(0)}%</span>
              <span class="pocket-resis">${escapeHtml(pk.residues.slice(0, 10).join(" "))}${pk.residues.length > 10 ? " …" : ""}</span>
            </div>`).join("")
          : '<div class="measurement-table-empty">No cavities detected</div>';

        const occ = ligAtoms.length ? " (ranked by bound-ligand occupancy)" : "";
        dom.pocketsStatus.textContent =
          `Found ${top.length} pocket${top.length === 1 ? "" : "s"}${occ}. Volumes are grid estimates.`;
      } catch (err) {
        dom.pocketsStatus.textContent = `Pocket detection failed: ${err.message}`;
      } finally {
        dom.btnDetectPockets.disabled = false;
      }
    }, 20);
  }

  function safeSelect(model, sel) {
    try { return model.selectedAtoms(sel) || []; } catch (e) { return []; }
  }

  dom.btnDetectPockets.addEventListener("click", detectPockets);

  // ----- Ligand contacts -----
  // Find all residues with non-water heavy atom within 4 Å of any HETATM ligand.
  function analyzeContacts() {
    if (!state.models.length) { showCustomError("No structure loaded."); return; }
    dom.btnAnalyzeContacts.disabled = true;
    dom.contactsStatus.textContent = "Finding ligands…";

    try {
      const model = state.models[state.currentPdbIndex];
      const allAtoms = model.selectedAtoms({});

      // Ligands = HETATM with non-water residue names
      const waters = new Set(["HOH", "WAT", "H2O", "DOD"]);
      const ligands = allAtoms.filter((a) => {
        const resn = (a.resn || "").trim().toUpperCase();
        if (waters.has(resn)) return false;
        // 3Dmol.js exposes HETATMs as part of the same selection; we use atom.hetflag if present
        return a.hetflag === true || (a.resn && a.resn.length > 0 && !waters.has(resn) && a.elem !== "C" && a.elem !== "N" && a.elem !== "O" && a.elem !== "S" && a.elem !== "H");
      });

      if (ligands.length === 0) {
        dom.contactsStatus.textContent = "No bound ligand detected in this structure.";
        dom.contactList.innerHTML = "";
        return;
      }

      // Group ligand atoms by residue
      const ligandGroups = {};
      for (const a of ligands) {
        const key = `${a.chain || ""}:${a.resn}:${a.resi}`;
        (ligandGroups[key] = ligandGroups[key] || { atoms: [], resn: a.resn, chain: a.chain, resi: a.resi }).atoms.push(a);
      }

      const CUTOFF = 4.0; // Å
      const contacts = [];
      for (const [, lg] of Object.entries(ligandGroups)) {
        const contactSet = new Map(); // "chain:resi" -> minDistance
        for (const la of lg.atoms) {
          for (const a of allAtoms) {
            if (a === la) continue;
            if (a.resn === la.resn && a.chain === la.chain) continue; // skip same residue (intra-ligand)
            const dx = a.x - la.x, dy = a.y - la.y, dz = a.z - la.z;
            const d = Math.sqrt(dx * dx + dy * dy + dz * dz);
            if (d <= CUTOFF) {
              const key = `${a.chain || ""}:${a.resi}`;
              const prev = contactSet.get(key);
              if (!prev || d < prev.d) contactSet.set(key, { d, resn: a.resn, chain: a.chain, resi: a.resi, atom: a.atom });
            }
          }
        }
        const sorted = Array.from(contactSet.values()).sort((a, b) => a.d - b.d);
        if (sorted.length) {
          contacts.push({ ligand: `${lg.resn}:${lg.resi}`, chain: lg.chain, resi: lg.resi, residues: sorted });
        }
      }

      if (!contacts.length) {
        dom.contactsStatus.textContent = "Ligands found but no protein residues within 4 Å.";
        dom.contactList.innerHTML = "";
        return;
      }

      state.ligandContacts = contacts;

      // Highlight contact residues in 3D
      state.viewer.removeAllShapes();
      contacts.forEach((c, i) => {
        c.residues.forEach((r) => {
          state.viewer.addStyle({ chain: r.chain, resi: r.resi }, {
            stick: { radius: 0.18, color: 0xfbbf24 },
          });
        });
        // Draw contact lines from ligand center to each contact residue CA
        const lg = ligandGroups[`${c.chain}:${c.ligand.split(":")[1]}:${c.resi}`];
        if (lg) {
          let cx = 0, cy = 0, cz = 0;
          for (const a of lg.atoms) { cx += a.x; cy += a.y; cz += a.z; }
          const lCenter = { x: cx / lg.atoms.length, y: cy / lg.atoms.length, z: cz / lg.atoms.length };
          // Find CA for each contact residue
          c.residues.forEach((r) => {
            const ca = allAtoms.find((a) => a.atom === "CA" && a.chain === r.chain && a.resi === r.resi);
            if (ca) {
              state.viewer.addCylinder({
                start: lCenter, end: ca, radius: 0.04, color: 0xfbbf24, dashed: true,
              });
            }
          });
        }
      });
      state.viewer.render();

      dom.contactList.innerHTML = contacts.map((c) => {
        const rrows = c.residues.map((r) => {
          const t = r.atom === "C" ? "phobic" : r.atom === "N" || r.atom === "O" ? "polar" : "—";
          return `<div class="contact-row">
            <span>${r.resn}${r.resi}.${r.atom}</span>
            <span class="contact-dist">${r.d.toFixed(2)} Å</span>
            <span class="contact-type">${t}</span>
          </div>`;
        }).join("");
        return `<div style="padding: 8px 10px; border-bottom: 1px solid var(--border)">
          <b>${c.ligand}</b> — ${c.residues.length} residues within 4 Å
          ${rrows}
        </div>`;
      }).join("");
      dom.contactsStatus.textContent = `Found ${contacts.reduce((a, c) => a + c.residues.length, 0)} contacts across ${contacts.length} ligand${contacts.length === 1 ? "" : "s"}.`;
    } catch (err) {
      dom.contactsStatus.textContent = `Contact analysis failed: ${err.message}`;
    } finally {
      dom.btnAnalyzeContacts.disabled = false;
    }
  }
  dom.btnAnalyzeContacts.addEventListener("click", analyzeContacts);

  dom.btnClearBioGeo.addEventListener("click", () => {
    state.pockets = [];
    state.ligandContacts = null;
    state.viewer.removeAllShapes();
    applyAllStyles();
    state.viewer.render();
    dom.pocketList.innerHTML = "";
    dom.contactList.innerHTML = "";
    dom.pocketsStatus.textContent = "No pockets detected.";
    dom.contactsStatus.textContent = "No ligand contacts.";
  });

  // ====== TOUR ======
  function updateTourAvailability() {
    const protein = state.currentProtein;
    const tours = TOURS[protein.id];
    if (tours && tours.length > 0) {
      dom.tourEmpty.hidden = true;
      dom.tourContent.hidden = false;
      const tour = tours[0];
      dom.tourTitle.textContent = tour.title;
      dom.tourSteps.innerHTML = "";
      tour.steps.forEach((step, i) => {
        const div = document.createElement("div");
        div.className = "tour-step";
        div.dataset.index = i;
        div.innerHTML = `
          <div class="tour-step-title">
            <span class="tour-step-num">${i + 1}</span>
            ${escapeHtml(step.title)}
          </div>
          <div class="tour-step-body">${escapeHtml(step.body.slice(0, 100))}${step.body.length > 100 ? "…" : ""}</div>
        `;
        div.addEventListener("click", () => {
          stopTour();
          state.tourStepIndex = i;
          runTourStep(i);
        });
        dom.tourSteps.appendChild(div);
      });
      dom.tourProgressBar.style.width = "0%";
      dom.btnTourPlay.disabled = false;
    } else {
      dom.tourEmpty.hidden = false;
      dom.tourContent.hidden = true;
    }
  }
  dom.btnTourPlay.addEventListener("click", () => { if (!state.tourRunning) startTour(); });
  dom.btnTourNext.addEventListener("click", () => {
    if (!state.tourRunning) {
      const tours = TOURS[state.currentProtein.id];
      if (!tours) return;
      const nextIdx = state.tourStepIndex + 1;
      if (nextIdx < tours[0].steps.length) {
        state.tourStepIndex = nextIdx;
        runTourStep(nextIdx);
        updateTourUI();
      }
    }
  });
  dom.btnTourStop.addEventListener("click", () => stopTour());

  function startTour() {
    const tours = TOURS[state.currentProtein.id];
    if (!tours) return;
    state.tourRunning = true;
    state.tourStepIndex = 0;
    dom.btnTourPlay.disabled = true;
    dom.btnTourNext.disabled = true;
    runTourStep(0);
  }
  function stopTour() {
    state.tourRunning = false;
    if (state.tourTimeout) { clearTimeout(state.tourTimeout); state.tourTimeout = null; }
    state.viewer.spin(false);
    state.spinning = false;
    dom.btnSpin.textContent = "↻ Spin";
    hideAnnotation();
    hideStepOverlay();
    dom.btnTourPlay.disabled = false;
    dom.btnTourNext.disabled = false;
    updateTourUI();
  }
  function runTourStep(idx) {
    const tours = TOURS[state.currentProtein.id];
    if (!tours) return;
    const tour = tours[0];
    const step = tour.steps[idx];
    if (!step) { stopTour(); return; }
    state.tourStepIndex = idx;
    updateTourUI();
    showStepOverlay(idx, tour.steps.length);
    applyTourAction(step);
    showAnnotation(step.title, step.body);
    if (state.tourRunning) {
      state.tourTimeout = setTimeout(() => {
        const nextIdx = idx + 1;
        if (nextIdx < tour.steps.length) runTourStep(nextIdx);
        else stopTour();
      }, step.duration);
    }
    updateUrlDebounced();
  }
  function updateTourUI() {
    const tours = TOURS[state.currentProtein.id];
    if (!tours) return;
    const total = tours[0].steps.length;
    const current = state.tourStepIndex;
    document.querySelectorAll(".tour-step").forEach((el) => {
      const i = parseInt(el.dataset.index, 10);
      el.classList.toggle("active", i === current && state.tourRunning);
      el.classList.toggle("completed", i < current);
    });
    const pct = total > 0 ? ((current + 1) / total) * 100 : 0;
    dom.tourProgressBar.style.width = `${pct}%`;
  }
  function showAnnotation(title, body) {
    dom.annTitle.textContent = title;
    dom.annBody.textContent = body;
    dom.annotation.hidden = false;
  }
  function hideAnnotation() {
    dom.annotation.hidden = true;
  }
  function showStepOverlay(current, total) {
    dom.stepCurrent.textContent = current + 1;
    dom.stepTotal.textContent = total;
    dom.stepOverlay.hidden = false;
  }
  function hideStepOverlay() {
    dom.stepOverlay.hidden = true;
  }
  function applyTourAction(step) {
    const a = step.action || {};
    state.viewer.removeAllShapes();
    if (typeof a.pdbIndex === "number" && state.models[a.pdbIndex]) {
      state.currentPdbIndex = a.pdbIndex;
      updateStateBadge(a.pdbIndex);
      applyAllStyles();
    }
    if (a.style && a.colorScheme) {
      state.currentStyle = a.style;
      state.currentColor = a.colorScheme;
      dom.ctrlStyle.value = a.style;
      dom.ctrlColor.value = a.colorScheme;
      applyAllStyles();
    }
    if (a.surface) {
      state.viewer.addSurface($3Dmol.SurfaceType.VDW, { ...resolveColor(state.currentColor, a.surface.opacity), opacity: a.surface.opacity }, {});
    }
    if (a.showSticks) {
      const target = a.showSticks;
      const sel = {};
      if (target.resi) sel.resi = target.resi;
      if (target.resn) sel.resn = target.resn;
      if (target.atom) sel.atom = target.atom;
      state.viewer.addStyle(sel, {
        stick: { radius: 0.2, color: "yellow" },
        sphere: { scale: 0.35, color: "yellow" },
      });
    }
    if (a.showSpheres) {
      const target = a.showSpheres;
      const sel = {};
      if (target.resi) sel.resi = target.resi;
      if (target.resn) sel.resn = target.resn;
      state.viewer.addStyle(sel, { sphere: { scale: 0.5, color: "magenta" } });
    }
    if (a.highlight) {
      const sel = {};
      if (a.highlight.resi) sel.resi = a.highlight.resi;
      if (a.highlight.resn) sel.resn = a.highlight.resn;
      state.viewer.addStyle(sel, {
        stick: { radius: 0.2, color: "yellow" },
        sphere: { scale: 0.35, color: "yellow" },
      });
    }
    if (a.focusHeme) {
      try { state.viewer.zoomTo({ resn: ["HEM", "HEME", "IHV"] }, 0, false); }
      catch (e) { state.viewer.zoomTo(); }
    } else if (a.zoomChain) {
      try { state.viewer.zoomTo({ chain: a.zoomChain }); }
      catch (e) { state.viewer.zoomTo(); }
    } else {
      state.viewer.zoomTo();
    }
    if (a.spin) {
      state.spinning = true;
      dom.btnSpin.textContent = "⏸ Pause";
      state.viewer.spin("y", 0.5);
    } else {
      state.viewer.spin(false);
      state.spinning = false;
      dom.btnSpin.textContent = "↻ Spin";
    }
    state.viewer.render();
  }

  // ====== TRAJECTORY ======
  function updateTrajectoryAvailability() {
    const protein = state.currentProtein;
    const traj = TRAJECTORIES[protein.id];
    if (traj) {
      dom.trajectoryEmpty.hidden = true;
      dom.trajectoryContent.hidden = false;
      dom.trajectoryTitle.textContent = traj.title;
      dom.trajectoryDesc.textContent = traj.description;
      dom.trajStateA.textContent = traj.states[traj.pdbIndexA];
      dom.trajStateB.textContent = traj.states[traj.pdbIndexB];
      // Default to state A fully opaque. 50% would superimpose both conformations,
      // which reads as a washed-out ghost rather than a clean structure.
      dom.trajSlider.value = 0;
      dom.trajPercent.textContent = "0";
      state.trajectoryValue = 0;
      stopTrajectory();
      updateTrajectoryBlend(0);
    } else {
      dom.trajectoryEmpty.hidden = false;
      dom.trajectoryContent.hidden = true;
    }
  }
  function setupTrajectoryInitial() { updateTrajectoryBlend(0); }
  function updateTrajectoryBlend(value) {
    if (!state.currentProtein) return;
    const traj = TRAJECTORIES[state.currentProtein.id];
    if (!traj) return;
    const idxA = traj.pdbIndexA;
    const idxB = traj.pdbIndexB;
    if (!state.models[idxA] || !state.models[idxB]) return;
    const t = value / 100;
    const opacityA = 1 - t;
    const opacityB = t;
    state.viewer.setStyle({ model: idxA }, { cartoon: resolveColor(state.currentColor, opacityA) });
    state.viewer.setStyle({ model: idxB }, { cartoon: resolveColor(state.currentColor, opacityB) });
    const protein = state.currentProtein;
    if (protein.highlights && protein.highlights.ligand) {
      const ligandResn = protein.highlights.ligand.resn;
      state.viewer.setStyle({ model: idxA, resn: ligandResn }, { stick: { color: "magenta", radius: 0.2 } });
      state.viewer.setStyle({ model: idxB, resn: ligandResn }, { stick: { color: "magenta", radius: 0.2 } });
    }
    state.viewer.render();
  }
  dom.trajSlider.addEventListener("input", (e) => {
    const v = parseInt(e.target.value, 10);
    state.trajectoryValue = v;
    dom.trajPercent.textContent = v;
    updateTrajectoryBlend(v);
  });
  dom.btnTrajPlay.addEventListener("click", () => { if (!state.trajectoryRunning) playTrajectory(); });
  dom.btnTrajPause.addEventListener("click", () => stopTrajectory());
  dom.btnTrajReset.addEventListener("click", () => {
    stopTrajectory();
    state.trajectoryValue = 0;
    dom.trajSlider.value = 0;
    dom.trajPercent.textContent = "0";
    updateTrajectoryBlend(0);
  });
  function playTrajectory() {
    state.trajectoryRunning = true;
    state.trajectoryDirection = 1;
    dom.btnTrajPlay.disabled = true;
    dom.btnTrajPause.disabled = false;
    dom.btnTrajReset.disabled = false;
    const speed = parseInt(dom.trajSpeed.value, 10);
    const startTime = performance.now();
    const startValue = state.trajectoryValue;
    const step = () => {
      if (!state.trajectoryRunning) return;
      const elapsed = performance.now() - startTime;
      const progress = elapsed / speed;
      const v = startValue + (100 - startValue) * progress * state.trajectoryDirection;
      if (state.trajectoryDirection === 1 && v >= 100) {
        state.trajectoryDirection = -1;
        requestAnimationFrame(() => startValueLoop(100));
        return;
      } else if (state.trajectoryDirection === -1 && v <= 0) {
        state.trajectoryDirection = 1;
        requestAnimationFrame(() => startValueLoop(0));
        return;
      }
      const clamped = Math.max(0, Math.min(100, v));
      state.trajectoryValue = clamped;
      dom.trajSlider.value = clamped;
      dom.trajPercent.textContent = Math.round(clamped);
      updateTrajectoryBlend(clamped);
      state.trajectoryRaf = requestAnimationFrame(step);
    };
    state.trajectoryRaf = requestAnimationFrame(step);
    function startValueLoop(startV) {
      if (!state.trajectoryRunning) return;
      const newStart = performance.now();
      const innerStep = () => {
        if (!state.trajectoryRunning) return;
        const elapsed = performance.now() - newStart;
        const progress = elapsed / speed;
        let v;
        if (state.trajectoryDirection === 1) {
          v = startV + (100 - startV) * progress;
          if (v >= 100) { state.trajectoryDirection = -1; requestAnimationFrame(() => startValueLoop(100)); return; }
        } else {
          v = startV - startV * progress;
          if (v <= 0) { state.trajectoryDirection = 1; requestAnimationFrame(() => startValueLoop(0)); return; }
        }
        const clamped = Math.max(0, Math.min(100, v));
        state.trajectoryValue = clamped;
        dom.trajSlider.value = clamped;
        dom.trajPercent.textContent = Math.round(clamped);
        updateTrajectoryBlend(clamped);
        state.trajectoryRaf = requestAnimationFrame(innerStep);
      };
      state.trajectoryRaf = requestAnimationFrame(innerStep);
    }
  }
  function stopTrajectory() {
    state.trajectoryRunning = false;
    if (state.trajectoryRaf) cancelAnimationFrame(state.trajectoryRaf);
    state.trajectoryRaf = null;
    dom.btnTrajPlay.disabled = false;
    dom.btnTrajPause.disabled = true;
    dom.btnTrajReset.disabled = false;
  }

  // ====== MODAL ======
  dom.btnAbout.addEventListener("click", () => { dom.aboutModal.hidden = false; });
  dom.btnCloseAbout.addEventListener("click", () => { dom.aboutModal.hidden = true; });
  dom.aboutModal.addEventListener("click", (e) => {
    if (e.target === dom.aboutModal) dom.aboutModal.hidden = true;
  });

  // ====== STRUCTURE METADATA (RCSB data API) ======
  async function loadStructureMetadata() {
    const protein = state.currentProtein;
    if (!protein) return;
    const gen = state._netGen;
    const pdbId = protein.pdbIds[state.currentPdbIndex];
    if (!/^[0-9A-Za-z]{4}$/.test(pdbId)) {
      dom.metaGrid.innerHTML = '<div class="muted small">Metadata only available for deposited PDB entries.</div>';
      return;
    }
    dom.metaGrid.innerHTML = '<div class="muted small">Loading metadata…</div>';
    try {
      const id = pdbId.toLowerCase();
      const entryRes = await net(`https://data.rcsb.org/rest/v1/core/entry/${id}`, undefined, gen);
      if (!entryRes.ok) throw new Error(`entry ${entryRes.status}`);
      const entry = await entryRes.json();
      const info = entry.rcsb_entry_info || {};
      const method = (entry.exptl && entry.exptl[0] && entry.exptl[0].method) || "—";
      const resArr = info.resolution_combined || [];
      const resolution = resArr.length ? `${resArr[resArr.length - 1].toFixed(2)} Å` : "—";
      // R-factors are not exposed in rcsb_entry_info; they live in the refine /
      // refine_hist blocks, and many older entries omit them entirely.
      const rwork = firstRefineStat(entry, "ls_r_factor_r_work", "ls_rfactor_r_work", "pdbx_rfactor_r_work");
      const rfree = firstRefineStat(entry, "ls_r_factor_r_free", "ls_rfactor_r_free", "pdbx_rfactor_r_free");
      const deposited = entry.rcsb_accession_info?.deposit_date || "—";
      const released = entry.rcsb_accession_info?.initial_release_date || "—";
      const title = (entry.struct && entry.struct.title) || protein.name;

      state.meta = { pdbId, method, resolution, rfree, rwork, deposited, released, title };

      // Chains and UniProt accession (SIFTS) for the cross-reference panel.
      const chains = [];
      const polymerIds = entry.rcsb_entry_container_identifiers?.polymer_entity_ids || [];
      let uniprot = null;
      const sifts = [];
      for (const eid of polymerIds) {
        try {
          const er = await net(`https://data.rcsb.org/rest/v1/core/polymer_entity/${id}/${eid}`, undefined, gen);
          if (!er.ok) continue;
          const ej = await er.json();
          const refs = ej.rcsb_polymer_entity_container_identifiers?.reference_sequence_identifiers || [];
          const desc = (ej.rcsb_polymer_entity_container_identifiers?.entity_poly || {}).pdbx_description;
          const chainsForEntity = (ej.entity_poly?.pdbx_strand_id || "").split(",").filter(Boolean);
          chains.push({ id: eid, description: desc || "polymer", chains: chainsForEntity });
          if (!uniprot) {
            for (const ref of refs) {
              if (ref.database_name === "UniProt") { uniprot = ref.database_accession; break; }
            }
          }
          // The UniProt residue range this entity covers. Needed to place
          // canonical-protein variants onto the construct: a ClinVar position
          // is in canonical UniProt numbering, but a deposited structure is
          // usually only a fragment of it.
          for (const al of ej.rcsb_polymer_entity_align || []) {
            if (al.reference_database_name !== "UniProt") continue;
            for (const reg of al.aligned_regions || []) {
              if (!reg.ref_beg_seq_id || !reg.length) continue;
              sifts.push({
                entityId: eid,
                uniprot: al.reference_database_accession,
                chains: chainsForEntity,
                minPos: reg.ref_beg_seq_id,
                maxPos: reg.ref_beg_seq_id + reg.length - 1,
              });
            }
          }
        } catch (e) { /* skip entity */ }
      }
      state.meta.chains = chains;
      state.meta.uniprot = uniprot;
      state.meta.sifts = sifts;
      state.meta.polymerEntities = polymerIds.length;
      state.meta.nonPolymerEntities = (entry.rcsb_entry_container_identifiers?.non_polymer_entity_ids || []).length;

      // Ligands: fetch each non-polymer entity
      const npIds = entry.rcsb_entry_container_identifiers?.non_polymer_entity_ids || [];
      const ligands = [];
      for (const eid of npIds) {
        try {
          const r = await net(`https://data.rcsb.org/rest/v1/core/nonpolymer_entity/${id}/${eid}`, undefined, gen);
          if (!r.ok) continue;
          const j = await r.json();
          const np = j.pdbx_entity_nonpoly || {};
          // comp_id lives under pdbx_entity_nonpoly for many entries; the container
          // identifiers block is not populated consistently, so fall back.
          const cid = j.rcsb_nonpolymer_entity_container_identifiers?.comp_id
            || j.pdbx_entity_nonpoly?.comp_id;
          const name = (np.name || "").trim();
          if (!cid || /^HOH$/i.test(cid)) continue;
          ligands.push({
            compId: cid.toUpperCase(),
            name: name || cid,
            formula: np.formula || null,
            entityId: eid,
          });
        } catch (e) { /* skip entity */ }
      }
      state.ligands = ligands;
      renderMetadataPanel();
      renderLigandList();
      refreshPresentationIfReady();
    } catch (err) {
      // A superseded request is not a failure — the user already moved on.
      if (isAbort(err) || gen !== netGeneration) return;
      dom.metaGrid.innerHTML = `<div class="muted small">Metadata unavailable (${escapeHtml(err.message)}). The structure itself is unaffected.</div>`;
      dom.ligandList.innerHTML = '<div class="measurement-table-empty">Metadata unavailable</div>';
    }
  }

  // Pull the first R-factor found under any of the known key spellings.
  // Note: the mmCIF-derived keys use capital R (ls_R_factor_R_work), while
  // some legacy payloads use lowercase — so compare case-insensitively.
  function firstRefineStat(entry, ...keys) {
    const blocks = []
      .concat(entry.refine || [])
      .concat(entry.refine_hist || []);
    const wanted = keys.map((k) => k.toLowerCase());
    for (const b of blocks) {
      if (!b || typeof b !== "object") continue;
      for (const k of Object.keys(b)) {
        if (!wanted.includes(k.toLowerCase())) continue;
        const v = b[k];
        if (v != null && !Array.isArray(v) && Number.isFinite(Number(v))) {
          return Number(v).toFixed(3);
        }
      }
    }
    return "—";
  }

  function renderMetadataPanel() {
    const m = state.meta;
    if (!m) return;
    const rows = [
      ["PDB ID", `<span class="res-badge">${escapeHtml(m.pdbId.toUpperCase())}</span>`],
      ["Method", escapeHtml(m.method)],
      ["Resolution", escapeHtml(m.resolution)],
      ["R-work / R-free", `${escapeHtml(m.rwork)} / ${escapeHtml(m.rfree)}`],
      ["Deposited", escapeHtml(m.deposited)],
      ["Released", escapeHtml(m.released)],
      ["Title", escapeHtml(m.title)],
    ];
    dom.metaGrid.innerHTML = rows
      .map(([k, v]) => `<div class="k">${k}</div><div class="v">${v}</div>`)
      .join("");
  }

  // ====== LIGAND BROWSER ======
  function renderLigandList() {
    const ligs = state.ligands || [];
    if (!ligs.length) {
      dom.ligandList.innerHTML = '<div class="measurement-table-empty">No ligands in this entry</div>';
      dom.ligandSelect.innerHTML = "";
      dom.propGrid.innerHTML = "";
      return;
    }
    dom.ligandList.innerHTML = ligs
      .map((l, i) => `
        <div class="ligand-row" data-lig="${i}">
          <span class="ligand-cid">${escapeHtml(l.compId)}</span>
          <span class="ligand-name" title="${escapeHtml(l.name)}">${escapeHtml(l.name)}</span>
          <span class="ligand-chain">${l.formula ? escapeHtml(l.formula) : ""}</span>
        </div>`)
      .join("");
    dom.ligandList.querySelectorAll(".ligand-row").forEach((row) => {
      row.addEventListener("click", () => {
        const i = parseInt(row.dataset.lig, 10);
        selectLigand(i);
      });
    });

    dom.ligandSelect.innerHTML = ligs
      .map((l, i) => `<option value="${i}">${escapeHtml(l.compId)} — ${escapeHtml(l.name)}</option>`)
      .join("");
    dom.ligandSelect.onchange = (e) => selectLigand(parseInt(e.target.value, 10));
    selectLigand(0);
  }

  function selectLigand(index) {
    const ligs = state.ligands || [];
    const l = ligs[index];
    if (!l) return;
    state.activeLigand = l;
    dom.ligandList.querySelectorAll(".ligand-row").forEach((r, i) => {
      r.classList.toggle("active", i === index);
    });
    dom.ligandSelect.value = String(index);

    // Isolate the ligand in the viewer
    applyAllStyles();
    state.viewer.addStyle({ resn: l.compId }, {
      stick: { radius: 0.22, color: 0xffd54f },
      sphere: { scale: 0.28, color: 0xffd54f },
    });
    try { state.viewer.zoomTo({ resn: l.compId }, 0, false); } catch (e) { /* ignore */ }
    state.viewer.render();

    renderLigandProperties(l);
  }

  // Physicochemical descriptors computed from the loaded model.
  function computeLigandProperties(compId) {
    const model = state.models[state.currentPdbIndex];
    if (!model) return null;
    let atoms;
    try { atoms = model.selectedAtoms({ resn: compId }); } catch (e) { return null; }
    if (!atoms || !atoms.length) return null;

    let heavy = 0, c = 0, n = 0, o = 0, s = 0, halogens = 0;
    const aromRings = new Set();
    for (const a of atoms) {
      const el = (a.elem || "").toUpperCase();
      if (el === "H") continue;
      heavy++;
      if (el === "C") { c++; if (a.bond === "ar" || a.aromatic) aromRings.add(a.resi + ":" + a.chain); }
      else if (el === "N") n++;
      else if (el === "O") o++;
      else if (el === "S") s++;
      else if (el === "F" || el === "CL" || el === "BR" || el === "I") halogens++;
    }
    if (!heavy) return null;

    const mw = c * 12.011 + n * 14.007 + o * 15.999 + s * 32.06 + halogens * 78.9;
    const hbd = Math.max(0, n);           // simplified: one donor per N
    const hba = o + n;                     // simplified: O and N are acceptors
    const tpsa = o * 20.23 + n * 25.78;    // Ertl fragment contributions
    const logp = (c * 1.46 + halogens * 1.23) - (o * 1.5 + n * 1.5) - 0.1; // crude Crippen-like
    const rotb = Math.max(0, heavy - 5);
    const lipinskiPass = mw <= 500 && logp <= 5 && hbd <= 5 && hba <= 10;
    const veberPass = rotb <= 10 && tpsa <= 140;

    return {
      heavy, mw, hbd, hba, tpsa, logp, rotb, lipinskiPass, veberPass,
      formula: `C${c} H? N${n} O${o}${s ? " S" + s : ""}${halogens ? " Hal" + halogens : ""}`,
    };
  }

  function renderLigandProperties(l) {
    const props = computeLigandProperties(l.compId);
    if (!props) {
      dom.ligandPropStatus.textContent = "Could not read atoms for this ligand.";
      dom.propGrid.innerHTML = "";
      return;
    }
    const ligName = (l.name || "").toLowerCase();
    const isNucleic = /\b(adenosine|guanosine|cytidine|uridine|thymidine|dna|rna|phosphate|cytidyl)\b/.test(ligName);
    dom.ligandPropStatus.textContent = isIonic(l) ? "Ion / cofactor" : isNucleic ? "Nucleotide / nucleic acid" : "Small-molecule ligand";
    const rows = [
      ["Atoms (heavy)", props.heavy],
      ["MW (approx)", props.mw.toFixed(1) + " Da"],
      ["H-bond donors", props.hbd + (props.hbd <= 5 ? " ✓" : " ✗")],
      ["H-bond acceptors", props.hba + (props.hba <= 10 ? " ✓" : " ✗")],
      ["TPSA (approx)", props.tpsa.toFixed(0) + " Å²"],
      ["cLogP (approx)", props.logp.toFixed(2)],
      ["Rotatable bonds", props.rotb],
      ["Lipinski Ro5", props.lipinskiPass ? "PASS" : "FAIL"],
      ["Veber", props.veberPass ? "PASS" : "FAIL"],
    ];
    dom.propGrid.innerHTML = rows
      .map(([k, v]) => {
        const cls = /FAIL|✗/.test(String(v)) ? "bad" : /PASS|✓/.test(String(v)) ? "good" : "";
        return `<div class="k">${k}</div><div class="v ${cls}">${escapeHtml(String(v))}</div>`;
      })
      .join("");
  }
  function isIonic(l) {
    return /^(ZN|FE|MG|CA|NA|K|MN|CL|CU|CO|NI|CD|HG|FE2|SO4|PO4)$/i.test(l.compId);
  }

  // ====== PROTEIN–LIGAND INTERACTION ANALYSIS ======
  const HBOND_CUTOFF = 3.5;
  const SALT_CUTOFF = 4.0;
  const HYDRO_CUTOFF = 4.0;
  const PI_CUTOFF = 5.5;
  const AROMATIC_ELEMENTS = new Set(["C", "N", "O", "S"]);
  const POSITIVE_RESIDUES = new Set(["ARG", "LYS", "HIS"]);
  const NEGATIVE_RESIDUES = new Set(["ASP", "GLU"]);

  function analyseInteractions() {
    const lig = state.activeLigand;
    if (!lig) {
      dom.hbondStatus.textContent = "Select a ligand first (Ligands panel).";
      return;
    }
    const model = state.models[state.currentPdbIndex];
    if (!model) return;
    let ligAtoms;
    let protAtoms;
    try {
      ligAtoms = model.selectedAtoms({ resn: lig.compId });
      protAtoms = model.selectedAtoms({ not: { resn: lig.compId } });
    } catch (e) {
      dom.hbondStatus.textContent = "Could not read model atoms.";
      return;
    }
    if (!ligAtoms.length) {
      dom.hbondStatus.textContent = `No atoms found for ${lig.compId}.`;
      return;
    }

    const hbonds = [], salts = [], pi = [], hydro = [];
    const seenHB = new Set();
    for (const la of ligAtoms) {
      if (!la.elem) continue;
      const le = la.elem.toUpperCase();
      for (const pa of protAtoms) {
        if (!pa.elem || !pa.resn) continue;
        const pe = pa.elem.toUpperCase();
        const d = computeDistance(la, pa);
        if (d > PI_CUTOFF) continue;

        // H-bond: N/O–N/O within 3.5 Å (geometric proxy; 3Dmol has no hydrogens here)
        if (d <= HBOND_CUTOFF && AROMATIC_ELEMENTS.has(le) && AROMATIC_ELEMENTS.has(pe)) {
          const key = [la.serial, pa.serial].sort((x, y) => x - y).join("-");
          if (!seenHB.has(key)) {
            seenHB.add(key);
            hbonds.push({ from: { chain: la.chain, resi: la.resi, atom: la.atom }, to: { chain: pa.chain, resi: pa.resi, resn: pa.resn, atom: pa.atom }, dist: d });
          }
        }
        // Salt bridge: charged side-chain groups within 4 Å
        if (d <= SALT_CUTOFF) {
          const pos = POSITIVE_RESIDUES.has(pa.resn);
          const neg = NEGATIVE_RESIDUES.has(pa.resn);
          const ligCharged = /^(N|O|S)$/.test(le);
          if ((pos || neg) && ligCharged) {
            const key = "s" + pa.chain + pa.resi + "-" + la.serial;
            if (!salts.some((s) => s.key === key)) salts.push({ key, chain: pa.chain, resi: pa.resi, resn: pa.resn, charge: pos ? "+" : "−", dist: d });
          }
        }
        // π-interaction: planar ring atom near a planar ring atom
        if (d <= PI_CUTOFF && AROMATIC_ELEMENTS.has(le) && AROMATIC_ELEMENTS.has(pe)) {
          if (isPlanarAromatic(pa) && isPlanarAromatic(la)) {
            const key = "p" + pa.chain + pa.resi + "-" + la.resi + la.chain;
            if (!pi.some((x) => x.key === key)) pi.push({ key, chain: pa.chain, resi: pa.resi, resn: pa.resn, dist: d });
          }
        }
        // Hydrophobic
        if (d <= HYDRO_CUTOFF && (le === "C") && isHydrophobicResidue(pa.resn)) {
          const key = "h" + pa.chain + pa.resi;
          if (!hydro.some((x) => x.key === key)) hydro.push({ key, chain: pa.chain, resi: pa.resi, resn: pa.resn, dist: d });
        }
      }
    }
    state.interactions = { hbonds, salts, pi, hydro, ligand: lig.compId };
    renderInteractions(lig);
  }

  function isPlanarAromatic(atom) {
    const RES = new Set(["PHE", "TYR", "TRP", "HIS"]);
    return RES.has(atom.resn);
  }
  function isHydrophobicResidue(resn) {
    return new Set(["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "GLY"]).has(resn);
  }

  function renderInteractions(lig) {
    const I = state.interactions;
    if (!I) {
      dom.hbondStatus.textContent = "No interactions detected.";
      dom.interactionSummary.innerHTML = "";
      return;
    }
    const total = I.hbonds.length + I.salts.length + I.pi.length + I.hydro.length;
    dom.hbondStatus.textContent = `${total} interactions for ${lig.compId} (cutoffs: H-bond ${HBOND_CUTOFF}Å, salt ${SALT_CUTOFF}Å, π ${PI_CUTOFF}Å).`;
    const rows = [
      ["H-bonds", I.hbonds.length, "hbond"],
      ["Salt bridges", I.salts.length, "salt"],
      ["π-interactions", I.pi.length, "pi"],
      ["Hydrophobic contacts", I.hydro.length, "hydro"],
    ];
    dom.interactionSummary.innerHTML = rows
      .map(([label, n, cls]) => `<div class="inter-row"><span>${label}</span><span class="inter-count ${cls}">${n}</span></div>`)
      .join("")
      + (I.salts.length ? `<div style="margin-top:8px;color:var(--text-2);font-size:10.5px">Salt bridges: ${I.salts.map((s) => `${s.resn}${s.resi}(${s.charge})`).join(", ")}</div>` : "")
      + (I.hydro.length ? `<div style="margin-top:4px;color:var(--text-2);font-size:10.5px">Hydrophobic: ${[...new Set(I.hydro.map((h) => h.resn + h.resi))].join(", ")}</div>` : "");

    // Visualise H-bonds as dashed lines
    const model = state.models[state.currentPdbIndex];
    if (model) {
      try {
        const ligAtoms = model.selectedAtoms({ resn: lig.compId });
        const byKey = new Map();
        for (const la of ligAtoms) {
          const k = `${la.chain}:${la.resi}:${la.atom}`;
          if (!byKey.has(k)) byKey.set(k, la);
        }
        for (const hb of I.hbonds) {
          const pa = model.selectedAtoms({ chain: hb.to.chain, resi: hb.to.resi, atom: hb.to.atom })[0];
          if (!pa) continue;
          state.viewer.addCylinder({
            start: { x: hb.from.x, y: hb.from.y, z: hb.from.z },
            end: { x: pa.x, y: pa.y, z: pa.z },
            radius: 0.04, color: 0x64b5f6, dashed: true,
          });
        }
        state.viewer.render();
      } catch (e) { /* ignore visual errors */ }
    }
  }
  dom.btnHbond.addEventListener("click", analyseInteractions);
  dom.btnClearInteractions.addEventListener("click", () => {
    state.interactions = null;
    dom.hbondStatus.textContent = "No interactions detected.";
    dom.interactionSummary.innerHTML = "";
    state.viewer.removeAllShapes();
    applyAllStyles();
    state.viewer.render();
  });

  // ====== MUTATION IMPACT ======
  function checkMutation() {
    const raw = (dom.mutInput.value || "").trim().toUpperCase();
    const m = raw.match(/^([A-Z]):?(\d+)([A-Z])$/);
    if (!m) {
      dom.mutationStatus.textContent = "Format: chain:positionNEW — e.g. A:790T";
      return;
    }
    const [, chain, pos, newRes] = m;
    const model = state.models[state.currentPdbIndex];
    if (!model) return;
    const atoms = model.selectedAtoms({ chain, resi: parseInt(pos, 10) });
    if (!atoms.length) {
      dom.mutationStatus.textContent = `Residue ${chain}:${pos} not found in this structure.`;
      return;
    }
    const ca = atoms.find((a) => a.atom === "CA") || atoms[0];
    // Packing density: count DISTINCT neighbouring residues whose Cα sits within
    // 8 A of this residue's Cα. This is the conventional CA-to-CA definition.
    // NB: 3Dmol ignores the `atom` filter when it is combined with `not`, so
    // the CA filter is applied in JS rather than in the selector.
    const neighbours = new Set();
    const all = model.selectedAtoms({ not: { chain, resi: parseInt(pos, 10) } });
    for (const a of all) {
      if (a.atom !== "CA") continue;
      const d = computeDistance(ca, a);
      if (d <= 8) neighbours.add(`${a.chain}:${a.resi}`);
    }
    const nRes = neighbours.size;
    // Thresholds calibrated on the residue-neighbour distribution of a real
    // globule (oxymyoglobin, 1MBO): p10 = 6, median = 9, p75 = 11.
    // >=10 is the buried upper quartile; >=6 is above the exposed decile.
    const burial = nRes >= 10 ? "buried core (high impact)"
                  : nRes >= 6  ? "partly buried (moderate)"
                  :               "surface-exposed (low)";
    let inPocket = "no";
    for (const p of state.pockets) {
      for (const r of p.residues) {
        if (r === `${chain}:${pos}`) { inPocket = `yes (#${state.pockets.indexOf(p) + 1})`; break; }
      }
    }
    dom.mutationStatus.innerHTML =
      `<b>${chain}:${pos}${newRes}</b> — wt <b>${ca.resn}</b><br>` +
      `Packing: <b>${nRes}</b> neighbouring residues within 8 Å → ${burial}<br>` +
      `In detected pocket: <b>${inPocket}</b>`;
  }
  dom.btnCheckMutation.addEventListener("click", checkMutation);
  dom.mutInput.addEventListener("keydown", (e) => {
    if (e.key === "Enter") { e.preventDefault(); checkMutation(); }
  });

  // ====== NAVIGATION ======
  dom.btnGoResidue.addEventListener("click", () => {
    const chain = (dom.resiChain.value || "A").trim().toUpperCase() || "A";
    const resi = parseInt(dom.resiNum.value, 10);
    if (!resi) return;
    focusResidue(chain, resi);
  });
  dom.resiNum.addEventListener("keydown", (e) => {
    if (e.key === "Enter") { e.preventDefault(); dom.btnGoResidue.click(); }
  });

  function renderChainToggles() {
    const model = state.models[state.currentPdbIndex];
    if (!model) return;
    const chains = new Set();
    let atoms;
    try { atoms = model.selectedAtoms({}); } catch (e) { return; }
    for (const a of atoms) if (a.chain) chains.add(a.chain);
    const list = Array.from(chains).sort();
    dom.chainToggles.innerHTML = list
      .map((c, i) => {
        const hex = "#" + STRUCT_CHAIN_PALETTE[i % STRUCT_CHAIN_PALETTE.length].toString(16).padStart(6, "0");
        const off = state.hiddenChains.has(c) ? " off" : "";
        return `<span class="chain-chip${off}" data-chain="${escapeHtml(c)}">
          <span class="swatch-dot" style="background:${hex}"></span>${escapeHtml(c)}</span>`;
      })
      .join("");
    dom.chainToggles.querySelectorAll(".chain-chip").forEach((chip) => {
      chip.addEventListener("click", () => {
        const c = chip.dataset.chain;
        if (state.hiddenChains.has(c)) state.hiddenChains.delete(c);
        else state.hiddenChains.add(c);
        applyChainVisibility();
        renderChainToggles();
      });
    });
  }
  function applyChainVisibility() {
    state.models.forEach((m, idx) => applyStyleToModel(m, idx));
    for (const c of state.hiddenChains) {
      state.viewer.setStyle({ chain: c }, { opacity: 0.06, cartoon: { opacity: 0.06 } });
    }
    state.viewer.render();
  }

  // ====== EXPORT ======
  function exportPdb() {
    const text = state._lastPdbText;
    if (!text) {
      dom.exportStatus.textContent = "Nothing to export yet.";
      return;
    }
    const hidden = state.hiddenChains;
    const lines = text.split(/\r?\n/);
    const kept = [];
    for (const line of lines) {
      if ((line.startsWith("ATOM") || line.startsWith("HETATM")) && line.substring(21, 22).trim()) {
        if (hidden.has(line.substring(21, 22).trim())) continue;
      }
      kept.push(line);
    }
    const blob = new Blob([kept.join("\n")], { type: "chemical/x-pdb" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    const name = (state.currentProtein?.id || "structure").replace(/[^\w-]+/g, "_");
    a.download = `${name}.pdb`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(a.href);
    dom.exportStatus.textContent = "Exported PDB.";
  }
  dom.btnExportPdb.addEventListener("click", exportPdb);

  function exportCsv() {
    if (!state.measurements.length) {
      dom.exportStatus.textContent = "No measurements to export.";
      return;
    }
    const rows = [["type", "value", "unit", "atoms"]];
    for (const m of state.measurements) {
      rows.push([
        m.type,
        m.value.toFixed(3),
        m.type === "distance" ? "A" : "deg",
        m.atoms.map((a) => `${a.chain || "_"}:${a.resn}${a.resi}.${a.atom}`).join(" "),
      ]);
    }
    const csv = rows.map((r) => r.map((c) => `"${String(c).replace(/"/g, '""')}"`).join(",")).join("\n");
    const blob = new Blob([csv], { type: "text/csv" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = `${(state.currentProtein?.id || "structure")}-measurements.csv`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(a.href);
    dom.exportStatus.textContent = "Exported CSV.";
  }
  dom.btnExportCsv.addEventListener("click", exportCsv);

  // ====== PRESENTATION: the bioinformatics record view ======
  // The principle: a deposited structure is a linked record, not a picture.
  // Present it in the order a bioinformatician consumes it —
  //   entry -> entities/chains -> ligands -> residues -> atoms
  // and attach provenance, an explicit annotation key, and honest uncertainty
  // at each level.

  // The annotation key lives in src/annotations.js, where it is cross-checked
  // against the palettes the app actually renders. It used to live here, and
  // it drifted: the key advertised five of the ten colours the variant
  // scheme can produce, so a drug-response sphere (EGFR T790M, L858R) had
  // no legend entry at all.
  const ANNOTATIONS = PFE.annotations.ANNOTATIONS;
  function renderPresentation() {
    const m = state.meta;
    const prot = state.currentProtein;

    // ---- evidence hierarchy ----
    const model = state.models[state.currentPdbIndex];
    let atomCount = 0, resCount = 0;
    const chainSet = new Set();
    if (model) {
      const atoms = safeSelect(model, {});
      atomCount = atoms.length;
      const resSet = new Set();
      for (const a of atoms) {
        if (a.chain) chainSet.add(a.chain);
        resSet.add(`${a.chain}:${a.resi}`);
      }
      resCount = resSet.size;
    }
    const levels = [
      {
        name: "Entry", value: m ? m.pdbId.toUpperCase() : "—",
        note: m ? `${m.method}${m.resolution !== "—" ? " @ " + m.resolution : ""}` : "awaiting metadata",
      },
      {
        name: "Polymer entities", value: m ? `${m.polymerEntities ?? "—"}` : "—",
        note: chainSet.size ? `chains ${Array.from(chainSet).sort().join(", ")}` : "—",
      },
      {
        name: "Ligands / hetero", value: m ? `${m.nonPolymerEntities ?? "—"}` : "—",
        note: (state.ligands || []).length
          ? state.ligands.map((l) => l.compId).join(", ")
          : "none deposited",
      },
      { name: "Residues (modelled)", value: resCount || "—", note: "unique chain:resi with coordinates" },
      { name: "Atoms", value: atomCount || "—", note: "as parsed from the coordinate file" },
    ];
    dom.hierarchy.innerHTML = levels.map((l) => `
      <li>
        <span class="h-name">${escapeHtml(l.name)}</span>
        <span class="h-value">${escapeHtml(String(l.value))}</span>
        <span class="h-note">${escapeHtml(l.note || "")}</span>
      </li>`).join("");

    // ---- cross-references ----
    const isAF = prot && /alphafold/i.test(prot.id || "");
    const up = m && m.uniprot;
    const xrefs = [
      { db: "PDB", id: m ? m.pdbId.toUpperCase() : null,
        url: m ? `https://www.rcsb.org/structure/${m.pdbId.toUpperCase()}` : null },
      { db: "UniProt", id: up || null,
        url: up ? `https://www.uniprot.org/uniprotkb/${up}` : null },
      { db: "AlphaFold", id: isAF && up ? `AF-${up}-F1` : null,
        url: isAF && up ? `https://alphafold.ebi.ac.uk/entry/${up}` : null },
      { db: "Method", id: m ? m.method : null, url: null },
    ];
    dom.xrefList.innerHTML = xrefs.map((x) => `
      <div class="xref">
        <span class="x-db">${escapeHtml(x.db)}</span>
        <span class="x-id ${x.id ? "" : "none"}">${escapeHtml(x.id || "not resolved")}</span>
        ${x.url ? `<a class="x-go" href="${escapeHtml(x.url)}" target="_blank" rel="noopener">↗</a>` : "<span></span>"}
      </div>`).join("");

    // ---- annotation key ----
    // The five flags the annotation key knows about, named in one place so a
    // typo cannot silently disable a key entry.
    const available = PFE.annotations.availability({
      crystal: !!m && /X-RAY|NMR|EM|ELECTRON|CRYO/.test(m.method || ""),
      alphafold: !!isAF,
      conservation: !!state.conservation,
      domains: !!(state.currentProtein && state.currentProtein.domains && state.currentProtein.domains.length),
      variants: !!(state.variants && state.variants.some((v) => v.inStructure)),
    });
    dom.annotKey.innerHTML = ANNOTATIONS.map((a) => {
      const ok = !a.needs || available[a.needs];
      return `<div class="annot-row ${ok ? "" : "unavailable"}" data-annot="${escapeHtml(a.id)}">
        <span class="annot-swatch" style="background:${a.grad}"></span>
        <span>
          <span class="a-name">${escapeHtml(a.name)}</span><br>
          <span class="a-qty">${escapeHtml(a.quantity)}</span><br>
          <span class="a-use">${escapeHtml(a.use)}</span>
          ${ok ? "" : '<br><span class="a-use">unavailable for this structure</span>'}
        </span>
      </div>`;
    }).join("");
    dom.annotKey.querySelectorAll(".annot-row").forEach((row) => {
      row.addEventListener("click", () => {
        const id = row.dataset.annot;
        if (state.currentColor === id) return;
        dom.ctrlColor.value = id;
        dom.ctrlColor.dispatchEvent(new Event("change"));
        state.viewer.render();
        dom.annotKey.querySelectorAll(".annot-row").forEach((r) => r.classList.remove("active"));
        row.classList.add("active");
      });
    });

    // ---- uncertainty: two independent axes ----
    const resNum = m && m.resolution !== "—" ? parseFloat(m.resolution) : null;
    let modelVerdict = "unc-verdict-na", modelTxt = "No refinement statistics available for this entry.";
    if (resNum !== null) {
      if (resNum <= 1.8) { modelVerdict = "unc-verdict-good"; modelTxt = "Fine side chains and well-ordered waters. Side-chain rotamers are usually trustworthy."; }
      else if (resNum <= 2.5) { modelVerdict = "unc-verdict-good"; modelTxt = "Standard for crystallography. Backbone reliable; distal side chains and waters may be unreliable."; }
      else if (resNum <= 3.5) { modelVerdict = "unc-verdict-warn"; modelTxt = "Coarse. Treat side-chain positions and any ligand pose as hypotheses."; }
      else { modelVerdict = "unc-verdict-warn"; modelTxt = "Low resolution. Loops and side chains may not be resolved at all."; }
    }
    const isCrystal = m && /X-RAY/i.test(m.method || "");
    const bfactors = [];
    if (model && isCrystal) {
      const atoms = safeSelect(model, {});
      for (const a of atoms) if (typeof a.b === "number") bfactors.push(a.b);
    }
    let orderVerdict = "unc-verdict-na", orderTxt = "Not applicable to this entry.";
    if (isCrystal && bfactors.length) {
      bfactors.sort((x, y) => x - y);
      const med = bfactors[Math.floor(bfactors.length / 2)];
      const p95 = bfactors[Math.floor(bfactors.length * 0.95)];
      if (p95 > 100) { orderVerdict = "unc-verdict-warn"; orderTxt = `Median B ${med.toFixed(0)} Å², 95th pct ${p95.toFixed(0)} Å². Large regions are mobile or disordered — do not read conformation there.`; }
      else { orderVerdict = "unc-verdict-good"; orderTxt = `Median B ${med.toFixed(0)} Å², 95th pct ${p95.toFixed(0)} Å². Mostly well ordered; local flexibility still present.`; }
    } else if (/NMR/i.test((m && m.method) || "")) {
      orderVerdict = "unc-verdict-na"; orderTxt = "Solution NMR: an ensemble of models, each a sample. Read ensembles and state populations, not one model.";
    }
    dom.uncertainty.innerHTML = `
      <div class="unc">
        <div class="u-q">Is the model right? — resolution / method</div>
        <div class="u-a">${escapeHtml(modelTxt)}</div>
        <span class="u-verdict ${modelVerdict}">${modelVerdict === "unc-verdict-good" ? "USABLE" : modelVerdict === "unc-verdict-warn" ? "USE WITH CAUTION" : "NOT APPLICABLE"}</span>
      </div>
      <div class="unc">
        <div class="u-q">Is the region ordered? — B-factor / pLDDT</div>
        <div class="u-a">${escapeHtml(orderTxt)}</div>
        <span class="u-verdict ${orderVerdict}">${orderVerdict === "unc-verdict-good" ? "ORDERED" : orderVerdict === "unc-verdict-warn" ? "PARTLY DISORDERED" : "NOT APPLICABLE"}</span>
      </div>`;

    // ---- citation ----
    dom.citeBlock.textContent = m
      ? `${m.pdbId.toUpperCase()}. Deposited ${m.deposited ? m.deposited.slice(0, 10) : "?"}; released ${m.released ? m.released.slice(0, 10) : "?"}.${up ? ` UniProt ${up}.` : ""} Method ${m.method}${m.resolution !== "—" ? ", " + m.resolution : ""}. Retrieved from RCSB PDB.`
      : "Select a structure to generate a citation.";
  }

  dom.btnCopyCite.addEventListener("click", () => {
    const txt = dom.citeBlock.textContent;
    if (!txt || txt.startsWith("Select")) { showToast("Load a structure first", ""); return; }
    if (navigator.clipboard) {
      navigator.clipboard.writeText(txt).then(
        () => showToast("Citation copied", "success"),
        () => showToast("Copy failed", ""));
    } else {
      showToast("Clipboard unavailable", "");
    }
  });

  // Re-render the record whenever new metadata or analyses arrive.
  function refreshPresentationIfReady() {
    if (dom.hierarchy) renderPresentation();
  }

  // ====== BOOTSTRAP ======
  async function init() {
    initViewer();
    renderProteinList();
    updateMeasurementsUI();
    updateRmsdSummary();

    // Apply URL state if present
    const urlState = applyStateFromUrl();
    let initialProteinId = null;
    if (urlState.proteinId) {
      const found = PROTEINS.find((p) => p.id === urlState.proteinId);
      if (found) initialProteinId = found.id;
      else {
        // Maybe it's a custom-style id; if not in list, fetch
        const m = await tryLoadFromUrl(urlState.proteinId);
        if (m) initialProteinId = m.id;
      }
    }
    if (!initialProteinId) initialProteinId = PROTEINS[0].id;

    if (urlState.style) { state.currentStyle = urlState.style; dom.ctrlStyle.value = urlState.style; }
    if (urlState.color) { state.currentColor = urlState.color; dom.ctrlColor.value = urlState.color; }
    if (!isNaN(urlState.opacity)) { state.currentOpacity = urlState.opacity; dom.ctrlOpacity.value = urlState.opacity; }
    if (urlState.bg) { state.background = urlState.bg; dom.ctrlBg.value = urlState.bg; state.viewer.setBackgroundColor(urlState.bg === "black" ? "black" : urlState.bg === "grey" ? "0x202428" : "white"); }

    await loadProtein(initialProteinId);
    if (urlState.pdbIndex > 0) {
      // Switch to that PDB state if available
      const protein = state.currentProtein;
      if (protein && urlState.pdbIndex < protein.pdbIds.length) {
        // Trajectory proteins support state switching via tour action mechanism;
        // for non-trajectory, ignore.
      }
    }
  }

  async function tryLoadFromUrl(id) {
    // Custom URL like /?p=custom-4HHB — try to fetch from RCSB if it looks like a PDB
    const m = id.match(/^custom-([0-9][A-Z0-9]{3})$/i);
    if (m) {
      const pdbId = m[1].toUpperCase();
      showLoading(`Fetching ${pdbId} from share link…`);
      try {
        const r = await net(`https://files.rcsb.org/download/${pdbId}.pdb`);
        if (!r.ok) throw new Error("HTTP " + r.status);
        const data = await r.text();
        const protein = makeCustomProtein({
          id: `custom-${pdbId}`,
          rawData: data, format: "pdb", source: "share",
        });
        PROTEINS.push(protein);
        return protein;
      } catch (err) {
        hideLoading();
        showCustomError(`Could not load shared PDB ${pdbId}: ${err.message}`);
        return null;
      }
    }
    return null;
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();