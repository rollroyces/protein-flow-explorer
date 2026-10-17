/* ===================================================================
 * Protein Flow Explorer — main app
 *
 *  - Library of curated structures (data/proteins.js)
 *  - Interactive 3Dmol.js viewer
 *  - Guided tours with annotations + camera moves
 *  - Conformational-change playback (slider between two PDB states)
 *  - Representation & color controls
 *  - Custom molecules (PDB-ID input + .pdb/.cif drag-and-drop)
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
    pdbInput: $("pdb-input"),
    btnLoadPdb: $("btn-load-pdb"),
    dropZone: $("drop-zone"),
    customError: $("custom-error"),

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

    // Modal
    btnAbout: $("btn-about"),
    aboutModal: $("about-modal"),
    btnCloseAbout: $("btn-close-about"),
  };

  // ===================================================================
  //  VIEWER SETUP
  // ===================================================================
  function initViewer() {
    state.viewer = $3Dmol.createViewer(dom.viewer, {
      backgroundColor: "white",
    });

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
  }

  // ===================================================================
  //  PROTEIN LIBRARY
  // ===================================================================
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

  function escapeHtml(s) {
    return String(s).replace(/[&<>"']/g, (c) => ({
      "&": "&amp;",
      "<": "&lt;",
      ">": "&gt;",
      '"': "&quot;",
      "'": "&#39;",
    })[c]);
  }

  dom.filterChips.addEventListener("click", (e) => {
    const btn = e.target.closest(".chip");
    if (!btn) return;
    document.querySelectorAll(".chip").forEach((c) => c.classList.remove("active"));
    btn.classList.add("active");
    renderProteinList(btn.dataset.filter);
  });

  // ===================================================================
  //  CUSTOM MOLECULES (PDB-ID input + file drop)
  // ===================================================================
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
      const protein = makeCustomProtein({
        id,
        rawData: data,
        format: "pdb",
        source: "pdb-id",
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
      rawData: text,
      format,
      source: "file",
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

  dom.btnLoadPdb.addEventListener("click", () => {
    loadCustomByPdbId(dom.pdbInput.value);
  });
  dom.pdbInput.addEventListener("keydown", (e) => {
    if (e.key === "Enter") {
      e.preventDefault();
      loadCustomByPdbId(dom.pdbInput.value);
    }
  });
  dom.pdbInput.addEventListener("input", (e) => {
    e.target.value = e.target.value.toUpperCase().replace(/[^0-9A-Z]/g, "");
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
    if (files && files.length > 0) {
      loadCustomFromFile(files[0]);
    }
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

  // ===================================================================
  //  PROTEIN LOADING
  // ===================================================================
  async function loadProtein(proteinId) {
    const protein = PROTEINS.find((p) => p.id === proteinId);
    if (!protein) return;

    state.currentProtein = protein;
    state.currentPdbIndex = 0;
    state.activeHighlights = [];

    dom.currentProteinName.textContent = protein.name;
    updateStateBadge(0);
    renderProteinList(document.querySelector(".chip.active").dataset.filter);

    dom.proteinDescription.textContent = protein.description;

    updateTourAvailability();
    updateTrajectoryAvailability();
    updateHighlightButtons();

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
        return model;
      });

      const models = await Promise.all(loadPromises);
      state.models = models;

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
    } catch (err) {
      console.error(err);
      dom.loadingText.textContent = `Could not load ${protein.name}. Check your connection.`;
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

  // ===================================================================
  //  STYLING
  // ===================================================================
  function applyAllStyles() {
    state.models.forEach((m, idx) => applyStyleToModel(m, idx));
  }

  function applyStyleToModel(model, idx) {
    const sel = { model: idx };
    const style = getStyleSpec(state.currentStyle, state.currentColor, state.currentOpacity);
    state.viewer.setStyle(sel, style);

    const protein = state.currentProtein;
    if (protein && protein.highlights && protein.highlights.ligand) {
      const ligandSel = { model: idx, resn: protein.highlights.ligand.resn };
      state.viewer.setStyle(ligandSel, { stick: { color: "magenta", radius: 0.18 } });
    }
  }

  function getStyleSpec(styleName, colorName, opacity) {
    const base = { opacity };
    switch (styleName) {
      case "cartoon":
        return { cartoon: { ...base, color: colorName } };
      case "surface":
        return { surface: { ...base, color: colorName, opacity: 0.85 } };
      case "ball+stick":
        return { stick: { ...base, color: colorName, radius: 0.15 } };
      case "stick":
        return { stick: { ...base, color: colorName, radius: 0.18 } };
      case "spacefill":
        return { sphere: { ...base, color: colorName, scale: 0.3 } };
      case "ribbon":
        return { ribbon: { ...base, color: colorName } };
      case "line":
        return { line: { ...base, color: colorName } };
      default:
        return { cartoon: { ...base, color: colorName } };
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

  // ===================================================================
  //  CONTROLS
  // ===================================================================
  dom.ctrlStyle.addEventListener("change", (e) => {
    state.currentStyle = e.target.value;
    applyAllStyles();
    if (state.currentProtein && TRAJECTORIES[state.currentProtein.id]) {
      updateTrajectoryBlend(state.trajectoryValue);
    } else {
      state.viewer.render();
    }
  });

  dom.ctrlColor.addEventListener("change", (e) => {
    state.currentColor = e.target.value;
    applyAllStyles();
    if (state.currentProtein && TRAJECTORIES[state.currentProtein.id]) {
      updateTrajectoryBlend(state.trajectoryValue);
    } else {
      state.viewer.render();
    }
  });

  dom.ctrlOpacity.addEventListener("input", (e) => {
    state.currentOpacity = parseFloat(e.target.value);
    applyAllStyles();
    if (state.currentProtein && TRAJECTORIES[state.currentProtein.id]) {
      updateTrajectoryBlend(state.trajectoryValue);
    } else {
      state.viewer.render();
    }
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
  });

  // ===================================================================
  //  HIGHLIGHT ACTIONS
  // ===================================================================
  function updateHighlightButtons() {
    const protein = state.currentProtein;
    if (!protein) return;
    const h = protein.highlights || {};
    dom.btnHighlightHeme.hidden = !h.ligand;
    dom.btnHighlightActive.hidden = !h.residues;
    if (state.currentProtein.id === "hiv-protease") {
      dom.btnHighlightInterface.hidden = false;
    } else {
      dom.btnHighlightInterface.hidden = true;
    }
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
    state.viewer.setStyle({}, { stick: { radius: 0.12, color: "byElement", opacity: 0.6 } });
    state.viewer.render();
  });

  dom.btnClearHighlights.addEventListener("click", () => {
    clearAllHighlights();
    applyAllStyles();
    state.viewer.zoomTo();
    state.viewer.render();
  });

  function clearAllHighlights() {
    state.viewer.removeAllShapes();
    state.viewer.removeAllLabels();
  }

  // ===================================================================
  //  VIEWER TOOLBAR
  // ===================================================================
  dom.btnResetView.addEventListener("click", () => {
    state.viewer.zoomTo();
    state.viewer.render();
  });

  dom.btnSpin.addEventListener("click", () => {
    state.spinning = !state.spinning;
    dom.btnSpin.textContent = state.spinning ? "⏸ Pause" : "↻ Spin";
    if (state.spinning) {
      state.viewer.spin("y", 0.6);
    } else {
      state.viewer.spin(false);
    }
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

  // ===================================================================
  //  TABS
  // ===================================================================
  dom.tabBtns.forEach((btn) => {
    btn.addEventListener("click", () => {
      const tab = btn.dataset.tab;
      dom.tabBtns.forEach((b) => b.classList.toggle("active", b === btn));
      dom.tabPanels.forEach((p) => p.classList.toggle("active", p.dataset.tabPanel === tab));
      if (tab !== "tour") stopTour();
      if (tab !== "trajectory") stopTrajectory();
    });
  });

  // ===================================================================
  //  TOURS
  // ===================================================================
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

  dom.btnTourPlay.addEventListener("click", () => {
    if (state.tourRunning) return;
    startTour();
  });

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
    if (state.tourTimeout) {
      clearTimeout(state.tourTimeout);
      state.tourTimeout = null;
    }
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
    clearAllHighlights();
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
      state.viewer.addSurface(
        $3Dmol.SurfaceType.VDW,
        { opacity: a.surface.opacity, color: state.currentColor },
        {}
      );
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
      const hemeSel = { resn: ["HEM", "HEME", "IHV"] };
      try { state.viewer.zoomTo(hemeSel, 0, false); }
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

  // ===================================================================
  //  TRAJECTORIES
  // ===================================================================
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

  function setupTrajectoryInitial() {
    updateTrajectoryBlend(50);
  }

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

    state.viewer.setStyle({ model: idxA }, {
      cartoon: { color: state.currentColor, opacity: opacityA },
    });
    state.viewer.setStyle({ model: idxB }, {
      cartoon: { color: state.currentColor, opacity: opacityB },
    });

    const protein = state.currentProtein;
    if (protein.highlights && protein.highlights.ligand) {
      const ligandResn = protein.highlights.ligand.resn;
      state.viewer.setStyle({ model: idxA, resn: ligandResn }, {
        stick: { color: "magenta", radius: 0.2 },
      });
      state.viewer.setStyle({ model: idxB, resn: ligandResn }, {
        stick: { color: "magenta", radius: 0.2 },
      });
    }

    state.viewer.render();
  }

  dom.trajSlider.addEventListener("input", (e) => {
    const v = parseInt(e.target.value, 10);
    state.trajectoryValue = v;
    dom.trajPercent.textContent = v;
    updateTrajectoryBlend(v);
  });

  dom.btnTrajPlay.addEventListener("click", () => {
    if (state.trajectoryRunning) return;
    playTrajectory();
  });

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
          if (v >= 100) {
            state.trajectoryDirection = -1;
            requestAnimationFrame(() => startValueLoop(100));
            return;
          }
        } else {
          v = startV - startV * progress;
          if (v <= 0) {
            state.trajectoryDirection = 1;
            requestAnimationFrame(() => startValueLoop(0));
            return;
          }
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

  // ===================================================================
  //  MODAL
  // ===================================================================
  dom.btnAbout.addEventListener("click", () => { dom.aboutModal.hidden = false; });
  dom.btnCloseAbout.addEventListener("click", () => { dom.aboutModal.hidden = true; });
  dom.aboutModal.addEventListener("click", (e) => {
    if (e.target === dom.aboutModal) dom.aboutModal.hidden = true;
  });

  // ===================================================================
  //  BOOTSTRAP
  // ===================================================================
  function init() {
    initViewer();
    renderProteinList();
    if (PROTEINS.length > 0) loadProtein(PROTEINS[0].id);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();