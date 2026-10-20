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
 *    - Bioinformatics analyses: pLDDT confidence colouring,
 *      MSA-derived conservation, ClinVar variant counts,
 *      druggable-pocket detection, ligand-contact analysis.
 * =================================================================== */

(function () {
  "use strict";

  // ====== GLOBAL STATE ======
  const state = {
    viewer: null,
    currentProtein: null,
    currentPdbIndex: 0,
    models: [],

    currentStyle: "cartoon",
    currentColor: "chain",
    currentOpacity: 1,
    background: "white",
    spinning: false,

    fpsFrames: 0,
    fpsLastTime: performance.now(),

    tourRunning: false,
    tourStepIndex: 0,
    tourTimeout: null,

    trajectoryRunning: false,
    trajectoryRaf: null,
    trajectoryValue: 50,
    trajectoryDirection: 1,

    activeHighlights: [],
    currentFilter: "all",
    customCounter: 0,

    // ----- Researcher extensions -----
    pickMode: false,
    measurements: [],
    pickBuffer: [],
    sequence: null,
    sequenceShowSS: false,
    compareProtein: null,
    compareLoaded: false,
    aligned: false,
    rmsdPerResidue: null,
    rmsdColorMode: false,

    // ----- Bioinformatics analysis results -----
    conservation: null,
    conservationJob: null,
    variants: [],
    pockets: [],
    ligandContacts: null,
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

    sequenceChains: $("sequence-chains"),
    sequenceInfo: $("sequence-info"),
    btnToggleSecondary: $("btn-toggle-secondary"),
    sequenceSecondary: $("sequence-secondary"),

    tabBtns: document.querySelectorAll(".tab-btn"),
    tabPanels: document.querySelectorAll(".tab-panel"),

    tourEmpty: $("tour-empty"),
    tourContent: $("tour-content"),
    tourTitle: $("tour-title"),
    tourSteps: $("tour-steps"),
    btnTourPlay: $("btn-tour-play"),
    btnTourNext: $("btn-tour-next"),
    btnTourStop: $("btn-tour-stop"),
    tourProgressBar: $("tour-progress-bar"),

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

    btnPickAtom: $("btn-pick-atom"),
    btnClearMeasurements: $("btn-clear-measurements"),
    measurementTable: $("measurement-table"),
    btnEnableCompare: $("btn-enable-compare"),
    btnAlignCalpha: $("btn-align-calpha"),
    btnColorByRmsd: $("btn-color-by-rmsd"),
    rmsdSummary: $("rmsd-summary"),

    btnColorPlddt: $("btn-color-plddt"),
    btnFetchConservation: $("btn-fetch-conservation"),
    btnClearConservation: $("btn-clear-conservation"),
    conservationStatus: $("conservation-status"),
    btnFetchVariants: $("btn-fetch-variants"),
    btnClearVariants: $("btn-clear-variants"),
    variantsStatus: $("variants-status"),
    btnDetectPockets: $("btn-detect-pockets"),
    btnAnalyzeContacts: $("btn-analyze-contacts"),
    btnClearBioGeo: $("btn-clear-bio-geo"),
    pocketsStatus: $("pockets-status"),
    contactsStatus: $("contacts-status"),
    pocketList: $("pocket-list"),
    contactList: $("contact-list"),

    btnAbout: $("btn-about"),
    aboutModal: $("about-modal"),
    btnCloseAbout: $("btn-close-about"),
    btnShare: $("btn-share"),

    toast: $("toast"),
  };

  // ====== UTILITIES ======
  function escapeHtml(s) {
    return String(s).replace(/[&<>"']/g, (c) => ({
      "&": "&", "<": "<", ">": ">", '"': """, "'": "&#39;",
    })[c]);
  }
  function categoryLabel(c) {
    return { allostery: "allostery", binding: "binding", "drug-target": "drug target",
             folding: "folding", domains: "domains", custom: "custom" }[c] || c;
  }
  function showToast(msg, kind = "") {
    if (!dom.toast) return;
    dom.toast.textContent = msg;
    dom.toast.className = "toast " + kind;
    dom.toast.hidden = false;
    clearTimeout(showToast._t);
    showToast._t = setTimeout(() => { dom.toast.hidden = true; }, 3500);
  }
  function debounce(fn, ms) {
    let t = null;
    return (...args) => { clearTimeout(t); t = setTimeout(() => fn(...args), ms); };
  }

  // ====== AA codes ======
  const AA3TO1 = {
    ALA: "A", ARG: "R", ASN: "N", ASP: "D", CYS: "C",
    GLN: "Q", GLU: "E", GLY: "G", HIS: "H", ILE: "I",
    LEU: "L", LYS: "K", MET: "M", PHE: "F", PRO: "P",
    SER: "S", THR: "T", TRP: "W", TYR: "Y", VAL: "V",
    MSE: "M", SEC: "U", PYR: "U",
  };

  // ====== URL STATE ======
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
    return {
      proteinId: params.get("p"),
      pdbIndex: parseInt(params.get("i") || "0", 10),
      style: params.get("style"),
      color: params.get("color"),
      opacity: parseFloat(params.get("op") || ""),
      bg: params.get("bg"),
      tour: parseInt(params.get("tour") || "", 10),
    };
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
    state.viewer = $3Dmol.createViewer(dom.viewer, { backgroundColor: "white" });

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

    const canvas = dom.viewer.querySelector("canvas");
    if (canvas) {
      canvas.addEventListener("click", (e) => {
        if (!state.pickMode) return;
        const rect = canvas.getBoundingClientRect();
        const x = e.clientX - rect.left;
        const y = e.clientY - rect.top;
        const atom = state.viewer.getIntersection(x, y);
        if (atom) handleAtomPicked(atom);
      });
    }
  }

  // ====== PROTEIN LIBRARY ======
  function renderProteinList(filter = "all") {
    state.currentFilter = filter;
    const list = filter === "all" ? PROTEINS : PROTEINS.filter((p) => p.category === filter);
    dom.proteinList.innerHTML = "";
    list.forEach((p) => {
      const card = document.createElement("div");
      card.className = "protein-card";
      if (p.isCustom) card.classList.add("protein-card-custom");
      card.dataset.id = p.id;
      card.style.borderLeftColor = p.color;
      if (state.currentProtein && state.currentProtein.id === p.id) card.classList.add("active");
      const pdbTags = p.pdbIds
        .map((id, i) => `<span class="pdb-tag">${escapeHtml(String(id))}${p.stateLabels.length > 1 ? "·" + escapeHtml(p.stateLabels[i].split(" ")[0]) : ""}</span>`)
        .join("");
      const removeBtn = p.isCustom
        ? `<button class="remove-btn" data-remove="${escapeHtml(p.id)}" title="Remove">×</button>` : "";
      card.innerHTML = `
        ${removeBtn}
        <div class="protein-card-name">
          <span style="color: ${p.color}">●</span> ${escapeHtml(p.name)}
        </div>
        <div class="protein-card-summary">${escapeHtml(p.summary)}</div>
        <div class="protein-card-meta">
          ${pdbTags}
          <span class="category-tag">${categoryLabel(p.category)}</span>
        </div>`;
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
    showCustomError._t = setTimeout(() => { dom.customError.hidden = true; }, 6000);
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
    const pdbId = String(pdbIdRaw || "").trim().toUpperCase();
    if (!/^[0-9][A-Z0-9]{3}$/.test(pdbId)) {
      showCustomError(`Invalid PDB ID "${pdbIdRaw}". Use 4 characters, e.g. 4HHB.`);
      return;
    }
    const id = `custom-${pdbId}`;
    if (PROTEINS.some((p) => p.id === id)) {
      hideLoading();
      loadProtein(id);
      dom.pdbInput.value = "";
      return;
    }
    showLoading(`Fetching PDB ${pdbId}…`);
    try {
      const response = await fetch(`https://files.rcsb.org/download/${pdbId}.pdb`);
      if (!response.ok) throw new Error(`HTTP ${response.status}`);
      const data = await response.text();
      if (!data || data.length < 50) throw new Error("Empty or invalid response");
      const protein = makeCustomProtein({ id, rawData: data, format: "pdb", source: "pdb-id" });
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
    if (e.key === "Enter") { e.preventDefault(); loadCustomByPdbId(dom.pdbInput.value); }
  });
  dom.pdbInput.addEventListener("input", (e) => {
    e.target.value = e.target.value.toUpperCase().replace(/[^0-9A-Z]/g, "");
  });
  dom.btnLoadGene.addEventListener("click", () => loadByGeneOrUniprot(dom.geneInput.value));
  dom.geneInput.addEventListener("keydown", (e) => {
    if (e.key === "Enter") { e.preventDefault(); loadByGeneOrUniprot(dom.geneInput.value); }
  });

  ["dragenter", "dragover"].forEach((ev) => {
    dom.dropZone.addEventListener(ev, (e) => {
      e.preventDefault(); e.stopPropagation();
      dom.dropZone.classList.add("dragover");
    });
  });
  ["dragleave", "dragend"].forEach((ev) => {
    dom.dropZone.addEventListener(ev, (e) => {
      e.preventDefault(); e.stopPropagation();
      dom.dropZone.classList.remove("dragover");
    });
  });
  dom.dropZone.addEventListener("drop", (e) => {
    e.preventDefault(); e.stopPropagation();
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
    const query = String(queryRaw || "").trim();
    if (!query) {
      showCustomError("Type a gene name (e.g. TP53) or UniProt accession (e.g. P04637).");
      return;
    }
    showLoading(`Looking up ${query}…`);
    try {
      let accession = null;
      const isUniprotAcc = /^[A-Z][A-Z0-9]{5,9}$/.test(query);
      if (isUniprotAcc) {
        accession = query;
      } else {
        const url = `https://rest.uniprot.org/uniprotkb/search?query=gene:${encodeURIComponent(query)}+AND+organism_id:9606&format=json&size=1&fields=accession,gene_names,protein_name`;
        const r = await fetch(url);
        if (!r.ok) throw new Error(`UniProt search returned ${r.status}`);
        const j = await r.json();
        const hit = j.results && j.results[0];
        if (!hit) {
          hideLoading();
          showCustomError(`No UniProt entry found for "${query}".`);
          return;
        }
        accession = hit.primaryAccession;
      }
      showLoading(`Fetching cross-refs for ${accession}…`);
      const entryUrl = `https://rest.uniprot.org/uniprotkb/${accession}.json`;
      const er = await fetch(entryUrl);
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
        showLoading(`No PDB — fetching AlphaFold model for ${accession}…`);
        const afUrl = `https://alphafold.ebi.ac.uk/files/AF-${accession}-F1-model_v4.pdb`;
        const afr = await fetch(afUrl);
        if (!afr.ok) {
          hideLoading();
          showCustomError(`No PDB or AlphaFold model for ${accession}.`);
          return;
        }
        const data = await afr.text();
        state.customCounter++;
        const id = `custom-alphafold-${accession}`;
        const protein = makeCustomProtein({
          id, name: `AlphaFold ${accession}`,
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
    if (state.compareProtein && state.compareProtein.id !== proteinId) clearComparison();
    clearMeasurements();

    state.currentProtein = protein;
    state.currentPdbIndex = 0;

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
          const response = await fetch(url);
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
      state._lastPdbText = results[0].data;
      state._lastPdbFormat = results[0].format;

      applyAllStyles();
      state.viewer.zoomTo();
      state.viewer.render();

      state.spinning = false;
      dom.btnSpin.textContent = "↻ Spin";
      hideLoading();

      if (TRAJECTORIES[protein.id]) setupTrajectoryInitial();
      else setModelVisibility(0);

      parseAndRenderSequence();
      updateUrlDebounced();
    } catch (err) {
      console.error(err);
      dom.loadingText.textContent = `Could not load ${protein.name}. Check your connection.`;
    }
  }

  function updateStateBadge(idx) {
    const protein = state.currentProtein;
    if (!protein || !protein.stateLabels[idx]) { dom.currentStateBadge.textContent = "—"; return; }
    dom.currentStateBadge.textContent = protein.stateLabels[idx];
  }
  function showLoading(text) { dom.loadingText.textContent = text; dom.loadingOverlay.classList.add("active"); }
  function hideLoading() { dom.loadingOverlay.classList.remove("active"); }

  // ====== STYLING ======
  function applyAllStyles() {
    state.models.forEach((m, idx) => applyStyleToModel(m, idx));
    if (state.rmsdColorMode && state.compareLoaded) applyRmsdColors();
  }
  function applyStyleToModel(model, idx) {
    const sel = { model: idx };
    const style = getStyleSpec(state.currentStyle, state.currentColor, state.currentOpacity);
    state.viewer.setStyle(sel, style);
    const protein = state.currentProtein;
    if (protein && protein.highlights && protein.highlights.ligand) {
      state.viewer.setStyle({ model: idx, resn: protein.highlights.ligand.resn }, { stick: { color: "magenta", radius: 0.18 } });
    }
  }
  function getStyleSpec(styleName, colorName, opacity) {
    const base = { opacity };
    switch (styleName) {
      case "cartoon": return { cartoon: { ...base, color: colorName } };
      case "surface": return { surface: { ...base, color: colorName, opacity: 0.85 } };
      case "ball+stick": return { stick: { ...base, color: colorName, radius: 0.15 } };
      case "stick": return { stick: { ...base, color: colorName, radius: 0.18 } };
      case "spacefill": return { sphere: { ...base, color: colorName, scale: 0.3 } };
      case "ribbon": return { ribbon: { ...base, color: colorName } };
      case "line": return { line: { ...base, color: colorName } };
      default: return { cartoon: { ...base, color: colorName } };
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

  function onStyleChange() {
    applyAllStyles();
    if (state.currentProtein && TRAJECTORIES[state.currentProtein.id]) updateTrajectoryBlend(state.trajectoryValue);
    else state.viewer.render();
    updateUrlDebounced();
  }
  dom.ctrlStyle.addEventListener("change", (e) => { state.currentStyle = e.target.value; onStyleChange(); });
  dom.ctrlColor.addEventListener("change", (e) => { state.currentColor = e.target.value; onStyleChange(); });
  dom.ctrlOpacity.addEventListener("input", (e) => { state.currentOpacity = parseFloat(e.target.value); onStyleChange(); });
  dom.ctrlBg.addEventListener("change", (e) => {
    state.background = e.target.value;
    const bg = e.target.value === "white" ? "white" : e.target.value === "black" ? "black" : "0x202428";
    state.viewer.setBackgroundColor(bg);
    state.viewer.render();
    updateUrlDebounced();
  });

  // ====== HIGHLIGHTS ======
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
    state.viewer.setStyle({ resn: h.ligand.resn }, { stick: { color: "magenta", radius: 0.2 }, sphere: { color: "magenta", scale: 0.35 } });
    state.viewer.zoomTo({ resn: h.ligand.resn });
    state.viewer.render();
  });
  dom.btnHighlightActive.addEventListener("click", () => {
    const h = state.currentProtein.highlights;
    if (!h || !h.residues) return;
    state.viewer.removeAllShapes();
    h.residues.forEach((r) => {
      const sel = r.chain ? { chain: r.chain, resi: r.resi } : { resi: r.resi };
      state.viewer.addStyle(sel, { stick: { color: "yellow", radius: 0.25 }, sphere: { color: "yellow", scale: 0.4 } });
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
    state.viewer.setStyle({}, { stick: { radius: 0.12, color: "byElement", opacity: 0.6 } });
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
  dom.btnResetView.addEventListener("click", () => { state.viewer.zoomTo(); state.viewer.render(); });
  dom.btnSpin.addEventListener("click", () => {
    state.spinning = !state.spinning;
    dom.btnSpin.textContent = state.spinning ? "⏸ Pause" : "↻ Spin";
    if (state.spinning) state.viewer.spin("y", 0.6); else state.viewer.spin(false);
  });
  dom.btnScreenshot.addEventListener("click", () => {
    const dataURL = state.viewer.pngURI();
    const link = document.createElement("a");
    link.href = dataURL;
    const proteinName = state.currentProtein ? state.currentProtein.name.replace(/\s+/g, "-").toLowerCase() : "protein";
    link.download = `${proteinName}-${Date.now()}.png`;
    document.body.appendChild(link); link.click(); document.body.removeChild(link);
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
    });
  });

  // ====== SEQUENCE PANEL ======
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
  function parseSequenceFromPDB(text, format) {
    const chainsMap = {};
    const ssByResi = {};
    const lines = text.split(/\r?\n/);
    const seen = {};
    for (const line of lines) {
      if (!line.startsWith("ATOM") && !line.startsWith("HETATM")) continue;
      const chain = line.substring(21, 22).trim() || "A";
      const atomName = line.substring(12, 16).trim();
      if (atomName !== "CA" && atomName !== "P" && format === "pdb") continue;
      const resi = parseInt(line.substring(22, 26).trim(), 10);
      const resn = line.substring(17, 20).trim();
      if (!resi) continue;
      const key = chain + ":" + resi;
      if (seen[key]) continue;
      seen[key] = true;
      if (!chainsMap[chain]) chainsMap[chain] = [];
      chainsMap[chain].push({ resi, resn });
    }
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
    const chainIds = Object.keys(chainsMap).sort();
    return {
      chains: chainIds.map((c) => ({
        id: c,
        residues: chainsMap[c].map((r) => ({
          resi: r.resi, resn: r.resn, one: AA3TO1[r.resn] || "X",
          ss: ssByResi[c]?.[r.resi] || "",
        })),
      })),
    };
  }
  function renderSequencePanel() {
    if (!state.sequence) { dom.sequenceChains.innerHTML = ""; dom.sequenceInfo.textContent = "—"; return; }
    const seq = state.sequence;
    const totalRes = seq.chains.reduce((acc, c) => acc + c.residues.length, 0);
    dom.sequenceInfo.textContent = `${seq.chains.length} chain${seq.chains.length !== 1 ? "s" : ""} · ${totalRes} residues`;
    dom.sequenceChains.innerHTML = "";
    seq.chains.forEach((chain) => {
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
    if (!state.sequenceShowSS) { dom.sequenceSecondary.hidden = true; return; }
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
    applyAllStyles();
    state.viewer.addStyle(sel, {
      stick: { color: "yellow", radius: 0.3 },
      sphere: { color: "yellow", scale: 0.5 },
    });
    state.viewer.zoomTo(sel, 0, false);
    state.viewer.render();
    document.querySelectorAll(".sequence-residue").forEach((el) => {
      el.classList.toggle("active", el.title.includes(`${chainId}:${resi}`));
    });
  }
  dom.btnToggleSecondary.addEventListener("click", () => {
    state.sequenceShowSS = !state.sequenceShowSS;
    dom.btnToggleSecondary.textContent = state.sequenceShowSS ? "Hide SS" : "Show SS";
    renderSecondaryStructure();
  });

  // ====== MEASUREMENT ======
  function setPickMode(on) {
    state.pickMode = on;
    dom.pickHint.hidden = !on;
    dom.btnPickAtom.textContent = on ? "✕ Cancel picking" : "📍 Pick atoms";
    if (on) { state.pickBuffer = []; updatePickHint(); }
    else { state.pickBuffer = []; updatePickHint(); }
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
      serial: atom.serial, x: atom.x, y: atom.y, z: atom.z,
      chain: atom.chain || "", resi: atom.resi, resn: atom.resn,
      atom: atom.atom, elem: atom.elem, model: atom.model,
    };
    state.pickBuffer.push(summary);
    state.viewer.addSphere({ center: { x: atom.x, y: atom.y, z: atom.z }, radius: 0.4, color: "yellow" });

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
    if (type === "distance") {
      state.viewer.addCylinder({ start: atoms[0], end: atoms[1], radius: 0.1, color: "yellow", fromCap: 1, toCap: 1 });
      const mid = midPoint(atoms[0], atoms[1]);
      state.viewer.addLabel(`${value.toFixed(2)} Å`, { position: mid, backgroundColor: "rgba(0,0,0,0.7)", fontColor: "white", fontSize: 12, inFront: true });
    } else if (type === "angle") {
      state.viewer.addCylinder({ start: atoms[0], end: atoms[1], radius: 0.05, color: "yellow" });
      state.viewer.addCylinder({ start: atoms[1], end: atoms[2], radius: 0.05, color: "yellow" });
      state.viewer.addLabel(`${value.toFixed(1)}°`, { position: atoms[1], backgroundColor: "rgba(0,0,0,0.7)", fontColor: "white", fontSize: 12, inFront: true, yOffset: 1.5 });
    } else if (type === "dihedral") {
      for (let i = 0; i < 3; i++) {
        state.viewer.addCylinder({ start: atoms[i], end: atoms[i + 1], radius: 0.04, color: "cyan" });
      }
      const mid = midPoint(atoms[1], atoms[2]);
      state.viewer.addLabel(`τ ${value.toFixed(1)}°`, { position: mid, backgroundColor: "rgba(0,0,0,0.7)", fontColor: "white", fontSize: 12, inFront: true });
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
      const val = m.type === "distance" ? `${m.value.toFixed(3)} Å` : `${m.value.toFixed(2)}°`;
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
  function midPoint(a, b) { return { x: (a.x + b.x) / 2, y: (a.y + b.y) / 2, z: (a.z + b.z) / 2 }; }
  function computeDistance(a, b) { return Math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2); }
  function computeAngle(a, b, c) {
    const v1 = { x: a.x - b.x, y: a.y - b.y, z: a.z - b.z };
    const v2 = { x: c.x - b.x, y: c.y - b.y, z: c.z - b.z };
    const dot = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    const m1 = Math.sqrt(v1.x ** 2 + v1.y ** 2 + v1.z ** 2);
    const m2 = Math.sqrt(v2.x ** 2 + v2.y ** 2 + v2.z ** 2);
    const cos = Math.max(-1, Math.min(1, dot / (m1 * m2)));
    return Math.acos(cos) * 180 / Math.PI;
  }
  function computeDihedral(a, b, c, d) {
    const b1 = { x: b.x - a.x, y: b.y - a.y, z: b.z - a.z };
    const b2 = { x: c.x - b.x, y: c.y - b.y, z: c.z - b.z };
    const b3 = { x: d.x - c.x, y: d.y - c.y, z: d.z - c.z };
    const n1 = cross(b1, b2);
    const n2 = cross(b2, b3);
    const m1 = Math.sqrt(dot(n1, n1));
    const m2 = Math.sqrt(dot(n2, n2));
    const b2n = Math.sqrt(dot(b2, b2));
    const m = { x: b2.x / b2n, y: b2.y / b2n, z: b2.z / b2n };
    const x = { x: n1.x / m1, y: n1.y / m1, z: n1.z / m1 };
    const y = cross(m, x);
    return Math.atan2(dot(n2, y), dot(n2, x)) * 180 / Math.PI;
  }
  function cross(a, b) { return { x: a.y * b.z - a.z * b.y, y: a.z * b.x - a.x * b.z, z: a.x * b.y - a.y * b.x }; }
  function dot(a, b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

  // ====== COMPARISON + RMSD ======
  dom.btnEnableCompare.addEventListener("click", () => {
    dom.compareControls.hidden = !dom.compareControls.hidden;
    dom.btnEnableCompare.textContent = dom.compareControls.hidden ? "⇄ Load structure B" : "⇄ Hide compare input";
    if (!dom.compareControls.hidden) dom.compareInput.focus();
  });
  dom.btnLoadCompare.addEventListener("click", () => loadCompareProtein(dom.compareInput.value));
  dom.compareInput.addEventListener("keydown", (e) => {
    if (e.key === "Enter") { e.preventDefault(); loadCompareProtein(dom.compareInput.value); }
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
      const r = await fetch(`https://files.rcsb.org/download/${pdbId}.pdb`);
      if (!r.ok) throw new Error(`HTTP ${r.status}`);
      const data = await r.text();
      const model = state.viewer.addModel(data, "pdb");
      model.modelIndex = state.models.length;
      model.pdbId = pdbId;
      model.isCompare = true;
      state.viewer.setStyle({ model: model.modelIndex }, { cartoon: { color: "#9aa0a6", opacity: 0.55 } });
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
    if (state.rmsdColorMode) { applyRmsdColors(); showToast("Coloured by per-residue RMSD", "success"); }
    else { applyAllStyles(); showToast("Back to default colors", ""); }
  });

  function alignAndComputeRmsd() {
    if (!state.compareLoaded) { showToast("Load structure B first", ""); return; }
    const modelA = state.models.find((m) => !m.isCompare);
    const modelB = state.models.find((m) => m.isCompare);
    if (!modelA || !modelB) { showToast("Both structures must be loaded", ""); return; }
    showLoading("Aligning on Cα…");
    try {
      const ok = state.viewer.superpose(modelB, modelA, { atom: "CA" });
      if (!ok) throw new Error("superpose returned false");
      state.aligned = true;
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
  function computePerResidueRmsd(modelA, modelB) {
    const aAtoms = modelA.selectedAtoms({ atom: "CA" });
    const bAtoms = modelB.selectedAtoms({ atom: "CA" });
    const bMap = new Map();
    for (const a of bAtoms) bMap.set(`${a.chain}:${a.resi}`, a);
    const result = [];
    for (const a of aAtoms) {
      const b = bMap.get(`${a.chain}:${a.resi}`);
      if (!b) continue;
      const dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
      result.push({ chain: a.chain, resi: a.resi, value: Math.sqrt(dx * dx + dy * dy + dz * dz) });
    }
    return result;
  }
  function applyRmsdColors() {
    if (!state.rmsdPerResidue) return;
    function colorFor(d) {
      const t = Math.min(1, d / 10);
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
          return c === undefined ? 0x9aa0a6 : c;
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
          <div class="tour-step-body">${escapeHtml(step.body.slice(0, 100))}${step.body.length > 100 ? "…" : ""}</div>`;
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
  function hideAnnotation() { dom.annotation.hidden = true; }
  function showStepOverlay(current, total) {
    dom.stepCurrent.textContent = current + 1;
    dom.stepTotal.textContent = total;
    dom.stepOverlay.hidden = false;
  }
  function hideStepOverlay() { dom.stepOverlay.hidden = true; }
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
      state.viewer.addSurface($3Dmol.SurfaceType.VDW, { opacity: a.surface.opacity, color: state.currentColor }, {});
    }
    if (a.showSticks) {
      const target = a.showSticks;
      const sel = {};
      if (target.resi) sel.resi = target.resi;
      if (target.resn) sel.resn = target.resn;
      if (target.atom) sel.atom = target.atom;
      state.viewer.addStyle(sel, { stick: { radius: 0.2, color: "yellow" }, sphere: { scale: 0.35, color: "yellow" } });
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
      state.viewer.addStyle(sel, { stick: { radius: 0.2, color: "yellow" }, sphere: { scale: 0.35, color: "yellow" } });
    }
    if (a.focusHeme) {
      try { state.viewer.zoomTo({ resn: ["HEM", "HEME", "IHV"] }, 0, false); }
      catch (e) { state.viewer.zoomTo(); }
    } else if (a.zoomChain) {
      try { state.viewer.zoomTo({ chain: a.zoomChain }); }
      catch (e) { state.viewer.zoomTo(); }
    } else { state.viewer.zoomTo(); }
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
      dom.trajSlider.value = 50;
      dom.trajPercent.textContent = "50";
      state.trajectoryValue = 50;
      stopTrajectory();
      updateTrajectoryBlend(50);
    } else {
      dom.trajectoryEmpty.hidden = false;
      dom.trajectoryContent.hidden = true;
    }
  }
  function setupTrajectoryInitial() { updateTrajectoryBlend(50); }
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
    state.viewer.setStyle({ model: idxA }, { cartoon: { color: state.currentColor, opacity: opacityA } });
    state.viewer.setStyle({ model: idxB }, { cartoon: { color: state.currentColor, opacity: opacityB } });
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
    state.trajectoryValue = 50;
    dom.trajSlider.value = 50;
    dom.trajPercent.textContent = "50";
    updateTrajectoryBlend(50);
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

  // ====== BIOINFORMATICS: pLDDT, Conservation, Variants, Pockets, Ligand Contacts ======
  function plddtToColor(p) {
    if (p >= 90) return 0x0053d6;
    if (p >= 70) return 0x65cbf3;
    if (p >= 50) return 0xffdb13;
    return 0xff7d45;
  }
  function colorByPLDDT() {
    if (!state.models.length) { showCustomError("No structure loaded."); return; }
    if (!state.viewer) return;
    state.viewer.setStyle({}, {
      cartoon: { colorfunc: (atom) => plddtToColor(atom.b) },
    });
    state.viewer.render();
    showToast("Coloured by AlphaFold pLDDT confidence", "success");
  }
  dom.btnColorPlddt.addEventListener("click", colorByPLDDT);

  const NW_MATCH = 2;
  const NW_MISMATCH = -1;
  const NW_GAP = -2;
  function nwAlign(a, b) {
    const m = a.length, n = b.length;
    const F = new Int16Array((m + 1) * (n + 1));
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
    let i = m, j = n, ai = "", bj = "";
    while (i > 0 || j > 0) {
      const dir = Tb[i * (n + 1) + j];
      if (i > 0 && j > 0 && dir === 0) { ai = a[i - 1] + ai; bj = b[j - 1] + bj; i--; j--; }
      else if (i > 0 && dir === 1) { ai = a[i - 1] + ai; bj = "-" + bj; i--; }
      else { ai = "-" + ai; bj = b[j - 1] + bj; j--; }
    }
    return { a: ai, b: bj };
  }

  function progressiveMsa(seqs) {
    if (seqs.length === 0) return [];
    let consensus = seqs[0];
    let aligned = [consensus.split("")];
    for (let i = 1; i < seqs.length; i++) {
      const result = nwAlign(consensus, seqs[i]);
      const newAligned = aligned.map((row) => {
        const s = row.join("");
        const r = nwAlign(result.a, s);
        return r.b.split("");
      });
      newAligned.push(result.b.split(""));
      consensus = result.a;
      aligned = newAligned;
    }
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

  async function fetchRcsbUniprotMapping(pdbId) {
    try {
      const url = `https://data.rcsb.org/rest/v1/core/entry/${pdbId.toLowerCase()}`;
      const r = await fetch(url);
      if (!r.ok) return null;
      const j = await r.json();
      const entities = j.rcsb_entry_container_identifiers?.polymer_entity_ids || [];
      for (const eid of entities) {
        try {
          const er = await fetch(`https://data.rcsb.org/rest/v1/core/polymer_entity/${pdbId.toLowerCase()}/${eid}`);
          if (!er.ok) continue;
          const ej = await er.json();
          const refs = ej.rcsb_polymer_entity_container_identifiers?.reference_sequence_identifiers || [];
          for (const ref of refs) {
            if (ref.database_accession && ref.database_name === "UniProt") {
              return ref.database_accession;
            }
          }
        } catch (e) {}
      }
    } catch (e) {}
    return null;
  }

  function buildResidueIndex(pdbText) {
    const map = new Map();
    let seqIdx = -1;
    const seen = new Set();
    const arr = [];
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
      arr[seqIdx] = AA3TO1[resn] || "X";
    }
    return { map, sequence: arr.join("") };
  }

  function parseFastaSequences(text) {
    const out = [];
    const lines = text.split(/\r?\n/);
    let cur = "";
    for (const line of lines) {
      if (line.startsWith(">")) { if (cur) out.push(cur); cur = ""; }
      else cur += line.trim();
    }
    if (cur) out.push(cur);
    return out;
  }

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
        const fastaUrl = `https://rest.uniprot.org/uniprotkb/${acc}.fasta`;
        const fr = await fetch(fastaUrl);
        if (!fr.ok) throw new Error("UniProt FASTA fetch failed");
        const fastaText = await fr.text();
        const seqs = parseFastaSequences(fastaText);
        const canonical = seqs[0] || "";
        if (!canonical) throw new Error("Empty canonical sequence");

        dom.conservationStatus.textContent = "Fetching homologs…";
        const protName = state.currentProtein.name.split(/[·•:()\-]/)[0].trim();
        const searchUrl = `https://rest.uniprot.org/uniprotkb/search?query=reviewed:true+AND+organism_id:9606+AND+name:${encodeURIComponent(protName)}&format=fasta&size=10`;
        const sr = await fetch(searchUrl);
        const homologs = sr.ok ? parseFastaSequences(await sr.text()) : [];

        const allSeqs = [canonical, ...homologs].slice(0, 10);
        dom.conservationStatus.textContent = `Aligning ${allSeqs.length} sequences…`;
        const msa = progressiveMsa(allSeqs);
        const scores = conservationScores(msa);

        const { map: residueMap } = buildResidueIndex(state._lastPdbText);
        const conservationByResidue = new Map();
        for (const [key, v] of residueMap) {
          if (v.seqIdx < scores.length) conservationByResidue.set(key, scores[v.seqIdx]);
        }
        state.conservation = { byResidue: conservationByResidue, homologs: allSeqs.length - 1 };

        state.viewer.setStyle({}, {
          cartoon: { colorfunc: (atom) => {
            const c = conservationByResidue.get(`${atom.chain}:${atom.resi}`);
            if (c === undefined) return 0xcccccc;
            const t = c;
            const r = Math.round(0x5b + (0xcc - 0x5b) * t);
            const g = Math.round(0x8d + (0x1f - 0x8d) * t);
            const b = Math.round(0xef + (0x1f - 0xef) * t);
            return (r << 16) | (g << 8) | b;
          } },
        });
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
    applyAllStyles();
    state.viewer.render();
    dom.conservationStatus.textContent = "No conservation loaded.";
  });

  async function fetchVariants() {
    if (!state.currentProtein) { showCustomError("Load a protein first."); return; }
    dom.variantsStatus.textContent = "Querying ClinVar…";
    dom.btnFetchVariants.disabled = true;
    try {
      const prot = state.currentProtein;
      const geneGuess =
        prot.id === "p53" ? "TP53" :
        prot.id === "egfr-kinase" ? "EGFR" :
        prot.id === "kras-g12c" ? "KRAS" :
        prot.id === "hiv-protease" ? "PR" :
        prot.id === "calmodulin" ? "CALM1" :
        prot.id === "insulin" ? "INS" :
        prot.id === "hemoglobin" ? "HBB" :
        prot.id === "lysozyme" ? "LYZ" :
        prot.id === "myoglobin" ? "MB" :
        null;
      if (!geneGuess) {
        dom.variantsStatus.textContent = "No gene-symbol mapping for this entry. Try TP53/EGFR/KRAS.";
        return;
      }
      const esearchUrl = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=${encodeURIComponent(geneGuess + "[gene]")}&retmax=200&retmode=json`;
      const r = await fetch(esearchUrl);
      if (!r.ok) throw new Error(`ClinVar esearch ${r.status}`);
      const j = await r.json();
      const ids = (j.esearchresult && j.esearchresult.idlist) || [];
      if (ids.length === 0) {
        dom.variantsStatus.textContent = `No ClinVar records found for ${geneGuess}.`;
        return;
      }
      dom.variantsStatus.textContent = `Found ${ids.length} ClinVar records. Summarising…`;
      const variants = [];
      const batch = ids.slice(0, 40);
      for (const id of batch) {
        try {
          const surl = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id=${id}&retmode=json`;
          const sr = await fetch(surl);
          if (!sr.ok) continue;
          const sj = await sr.json();
          const r2 = sj.result && sj.result[id];
          if (!r2) continue;
          const trait = (r2.clinical_significance && r2.clinical_significance.description) || "";
          const sig = trait.toLowerCase().includes("pathogenic") && !trait.toLowerCase().includes("benign")
            ? (trait.toLowerCase().includes("likely") ? "likely_pathogenic" : "pathogenic")
            : trait.toLowerCase().includes("benign")
            ? (trait.toLowerCase().includes("likely") ? "likely_benign" : "benign")
            : "vus";
          variants.push({ id, sig, trait });
        } catch (e) {}
      }
      const tally = { pathogenic: 0, likely_pathogenic: 0, vus: 0, likely_benign: 0, benign: 0 };
      for (const v of variants) tally[v.sig] = (tally[v.sig] || 0) + 1;
      state.variants = variants;
      dom.variantsStatus.innerHTML =
        `<b>${geneGuess}</b> · ${variants.length} variants from ${ids.length} ClinVar records:<br>` +
        `<span style="color:#b91c1c">● ${tally.pathogenic} pathogenic</span> · ` +
        `<span style="color:#f97316">● ${tally.likely_pathogenic} likely path.</span> · ` +
        `<span style="color:#fbbf24">● ${tally.vus} VUS</span> · ` +
        `<span style="color:#84cc16">● ${tally.likely_benign} likely ben.</span> · ` +
        `<span style="color:#22c55e">● ${tally.benign} benign</span><br>` +
        `<span class="muted small">Summary count by clinical significance. Per-residue mapping requires UniProt SIFTS lookup (next iteration).</span>`;
      showToast(`Loaded ${variants.length} ClinVar entries for ${geneGuess}`, "success");
    } catch (err) {
      dom.variantsStatus.textContent = `ClinVar fetch failed: ${err.message}`;
      showCustomError(`ClinVar fetch failed: ${err.message}`);
    } finally {
      dom.btnFetchVariants.disabled = false;
    }
  }
  dom.btnFetchVariants.addEventListener("click", fetchVariants);
  dom.btnClearVariants.addEventListener("click", () => {
    state.variants = [];
    dom.variantsStatus.textContent = "No variants loaded.";
    applyAllStyles();
    state.viewer.render();
  });

  function detectPockets() {
    if (!state.models.length) { showCustomError("No structure loaded."); return; }
    dom.btnDetectPockets.disabled = true;
    dom.pocketsStatus.textContent = "Scanning alpha-spheres…";
    try {
      const model = state.models[state.currentPdbIndex];
      const atoms = model.selectedAtoms({});
      if (!atoms.length) { dom.pocketsStatus.textContent = "No atoms in structure."; return; }
      let minX = Infinity, minY = Infinity, minZ = Infinity;
      let maxX = -Infinity, maxY = -Infinity, maxZ = -Infinity;
      for (const a of atoms) {
        if (a.x < minX) minX = a.x; if (a.x > maxX) maxX = a.x;
        if (a.y < minY) minY = a.y; if (a.y > maxY) maxY = a.y;
        if (a.z < minZ) minZ = a.z; if (a.z > maxZ) maxZ = a.z;
      }
      const PAD = 5;
      minX -= PAD; minY -= PAD; minZ -= PAD;
      maxX += PAD; maxY += PAD; maxZ += PAD;
      const GRID = 1.4;
      const PROBE_RADIUS_MAX = 6.0;
      const SCORE_THRESHOLD = 5;
      const cellSize = PROBE_RADIUS_MAX + 1;
      const grid = new Map();
      function cellKey(x, y, z) {
        return `${Math.floor(x / cellSize)}|${Math.floor(y / cellSize)}|${Math.floor(z / cellSize)}`;
      }
      for (const a of atoms) {
        const k = cellKey(a.x, a.y, a.z);
        let arr = grid.get(k);
        if (!arr) { arr = []; grid.set(k, arr); }
        arr.push(a);
      }
      function countNearby(x, y, z, radius) {
        const cx = Math.floor(x / cellSize), cy = Math.floor(y / cellSize), cz = Math.floor(z / cellSize);
        const span = Math.ceil(radius / cellSize);
        let count = 0;
        for (let dx = -span; dx <= span; dx++) {
          for (let dy = -span; dy <= span; dy++) {
            for (let dz = -span; dz <= span; dz++) {
              const arr = grid.get(`${cx + dx}|${cy + dy}|${cz + dz}`);
              if (!arr) continue;
              const r2 = radius * radius;
              for (let i = 0; i < arr.length; i++) {
                const a = arr[i];
                const ddx = a.x - x, ddy = a.y - y, ddz = a.z - z;
                if (ddx * ddx + ddy * ddy + ddz * ddz < r2) count++;
              }
            }
          }
        }
        return count;
      }
      const alphaCentres = [];
      for (let x = minX; x <= maxX; x += GRID) {
        for (let y = minY; y <= maxY; y += GRID) {
          for (let z = minZ; z <= maxZ; z += GRID) {
            const count = countNearby(x, y, z, PROBE_RADIUS_MAX);
            if (count >= SCORE_THRESHOLD * 2) {
              alphaCentres.push({ x, y, z, r: PROBE_RADIUS_MAX, score: count });
            }
          }
        }
      }
      const visited = new Set();
      const pockets = [];
      for (let i = 0; i < alphaCentres.length; i++) {
        if (visited.has(i)) continue;
        const seed = alphaCentres[i];
        const queue = [i];
        const cluster = [];
        visited.add(i);
        while (queue.length) {
          const idx = queue.shift();
          cluster.push(alphaCentres[idx]);
          for (let j = 0; j < alphaCentres.length; j++) {
            if (visited.has(j)) continue;
            const dx = alphaCentres[j].x - alphaCentres[idx].x;
            const dy = alphaCentres[j].y - alphaCentres[idx].y;
            const dz = alphaCentres[j].z - alphaCentres[idx].z;
            if (dx * dx + dy * dy + dz * dz < 9) { visited.add(j); queue.push(j); }
          }
        }
        if (cluster.length >= 3) {
          let cx = 0, cy = 0, cz = 0, score = 0;
          for (const p of cluster) { cx += p.x; cy += p.y; cz += p.z; score += p.score; }
          const n = cluster.length;
          pockets.push({
            center: [cx / n, cy / n, cz / n],
            volume: (4 / 3) * Math.PI * Math.pow(cluster[0].r, 3) * n * 0.4,
            score: score / n,
            residues: [],
          });
        }
      }
      for (const p of pockets) {
        const resiSet = new Set();
        const radius = Math.cbrt(p.volume * 0.75 / Math.PI) + 4;
        const r2 = radius * radius;
        for (const a of atoms) {
          const dx = a.x - p.center[0], dy = a.y - p.center[1], dz = a.z - p.center[2];
          if (dx * dx + dy * dy + dz * dz < r2) resiSet.add(`${a.chain}:${a.resi}`);
        }
        p.residues = Array.from(resiSet).slice(0, 25);
      }
      pockets.sort((a, b) => b.volume - a.volume);
      const top = pockets.slice(0, 8);
      state.pockets = top;

      state.viewer.removeAllShapes();
      top.forEach((p, i) => {
        const color = i === 0 ? 0xff7d45 : (i === 1 ? 0xffdb13 : 0x65cbf3);
        state.viewer.addSphere({
          center: { x: p.center[0], y: p.center[1], z: p.center[2] },
          radius: Math.cbrt(p.volume * 0.75 / Math.PI),
          color, opacity: 0.45,
        });
      });
      const allResis = new Set();
      for (const p of top) for (const r of p.residues) allResis.add(r);
      for (const r of allResis) {
        const [chain, resi] = r.split(":");
        state.viewer.addStyle({ chain, resi: parseInt(resi, 10) }, { stick: { radius: 0.12, color: 0xff7d45 } });
      }
      state.viewer.render();

      dom.pocketList.innerHTML = top.map((p, i) => `
        <div class="pocket-row">
          <span class="pocket-rank">#${i + 1}</span>
          <span class="pocket-vol">${p.volume.toFixed(0)} ų</span>
          <span class="pocket-name">score ${p.score.toFixed(1)}</span>
          <span class="pocket-resis">${p.residues.slice(0, 12).join(" ")}${p.residues.length > 12 ? " …" : ""}</span>
        </div>
      `).join("") || '<div class="measurement-table-empty">No pockets found</div>';
      dom.pocketsStatus.textContent = `Found ${top.length} pocket${top.length === 1 ? "" : "s"} (volumes in ų).`;
    } catch (err) {
      dom.pocketsStatus.textContent = `Pocket detection failed: ${err.message}`;
    } finally {
      dom.btnDetectPockets.disabled = false;
    }
  }
  dom.btnDetectPockets.addEventListener("click", detectPockets);

  function analyzeContacts() {
    if (!state.models.length) { showCustomError("No structure loaded."); return; }
    dom.btnAnalyzeContacts.disabled = true;
    dom.contactsStatus.textContent = "Finding ligands…";
    try {
      const model = state.models[state.currentPdbIndex];
      const allAtoms = model.selectedAtoms({});
      const waters = new Set(["HOH", "WAT", "H2O", "DOD"]);
      const ligands = allAtoms.filter((a) => {
        const resn = (a.resn || "").trim().toUpperCase();
        if (waters.has(resn)) return false;
        return a.hetflag === true;
      });
      if (ligands.length === 0) {
        dom.contactsStatus.textContent = "No bound ligand detected in this structure.";
        dom.contactList.innerHTML = "";
        return;
      }
      const ligandGroups = {};
      for (const a of ligands) {
        const key = `${a.chain || ""}:${a.resn}:${a.resi}`;
        (ligandGroups[key] = ligandGroups[key] || { atoms: [], resn: a.resn, chain: a.chain, resi: a.resi }).atoms.push(a);
      }
      const CUTOFF = 4.0;
      const contacts = [];
      for (const [lname, lg] of Object.entries(ligandGroups)) {
        const contactSet = new Map();
        for (const la of lg.atoms) {
          for (const a of allAtoms) {
            if (a === la) continue;
            if (a.resn === la.resn && a.chain === la.chain) continue;
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
      state.viewer.removeAllShapes();
      contacts.forEach((c) => {
        c.residues.forEach((r) => {
          state.viewer.addStyle({ chain: r.chain, resi: r.resi }, { stick: { radius: 0.18, color: 0xfbbf24 } });
        });
        const lg = ligandGroups[`${c.chain}:${c.ligand.split(":")[1]}:${c.resi}`];
        if (lg) {
          let cx = 0, cy = 0, cz = 0;
          for (const a of lg.atoms) { cx += a.x; cy += a.y; cz += a.z; }
          const lCenter = { x: cx / lg.atoms.length, y: cy / lg.atoms.length, z: cz / lg.atoms.length };
          c.residues.forEach((r) => {
            const ca = allAtoms.find((a) => a.atom === "CA" && a.chain === r.chain && a.resi === r.resi);
            if (ca) {
              try {
                state.viewer.addCylinder({
                  start: lCenter, end: ca, radius: 0.04, color: 0xfbbf24, dashed: true,
                });
              } catch (e) {}
            }
          });
        }
      });
      state.viewer.render();
      dom.contactList.innerHTML = contacts.map((c) => {
        const rrows = c.residues.map((r) => {
          const t = r.atom === "C" ? "phobic" : (r.atom === "N" || r.atom === "O" ? "polar" : "—");
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

  // ====== BOOTSTRAP ======
  async function init() {
    initViewer();
    renderProteinList();
    updateMeasurementsUI();
    updateRmsdSummary();

    const urlState = applyStateFromUrl();
    let initialProteinId = null;
    if (urlState.proteinId) {
      const found = PROTEINS.find((p) => p.id === urlState.proteinId);
      if (found) initialProteinId = found.id;
      else {
        const m = await tryLoadFromUrl(urlState.proteinId);
        if (m) initialProteinId = m.id;
      }
    }
    if (!initialProteinId) initialProteinId = PROTEINS[0].id;

    if (urlState.style) { state.currentStyle = urlState.style; dom.ctrlStyle.value = urlState.style; }
    if (urlState.color) { state.currentColor = urlState.color; dom.ctrlColor.value = urlState.color; }
    if (!isNaN(urlState.opacity)) { state.currentOpacity = urlState.opacity; dom.ctrlOpacity.value = urlState.opacity; }
    if (urlState.bg) {
      state.background = urlState.bg;
      dom.ctrlBg.value = urlState.bg;
      state.viewer.setBackgroundColor(urlState.bg === "black" ? "black" : urlState.bg === "grey" ? "0x202428" : "white");
    }

    await loadProtein(initialProteinId);
  }

  async function tryLoadFromUrl(id) {
    const m = id.match(/^custom-([0-9][A-Z0-9]{3})$/i);
    if (m) {
      const pdbId = m[1].toUpperCase();
      showLoading(`Fetching ${pdbId} from share link…`);
      try {
        const r = await fetch(`https://files.rcsb.org/download/${pdbId}.pdb`);
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