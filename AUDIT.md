# Protein Flow Explorer — Code & Architecture Audit

**Repository:** [`rollroyces/protein-flow-explorer`](https://github.com/rollroyces/protein-flow-explorer)
**Live:** <https://rollroyces.github.io/protein-flow-explorer/>
**Audit date:** 2026-10-17 · **Revision audited:** `01698b9` → `8cc1fc2`

Every finding below was measured against the working tree, not recalled. Where a
claim is quantitative, the measurement command or the test that proves it is
named. Three of the recommendations have been implemented, committed and
verified; the remainder are written up with the same evidence standard so they
can be actioned directly.

---

## 1. Core functionality — what this repository does

Protein Flow Explorer is a **static, client-side web application for teaching
and researching protein structure**. It has no backend and no database: the
entire application is HTML, CSS and JavaScript served from GitHub Pages, and it
fetches structures and metadata live from public APIs.

The feature surface, in the order a user meets it:

| Capability | What it does |
|---|---|
| **Curated library** | 13 hand-annotated structures (haemoglobin, EGFR, SARS-CoV-2 Mpro, a ribosome, collagen, and others) with fold class, biological role and highlighted functional sites. |
| **3D viewer** | 3Dmol.js rendering of the loaded structure with cartoon/ball-and-stick/surface styles, chain or pLDDT or element colouring, opacity and spin controls. |
| **Guided tours** | Scripted, step-by-step narration of a structure's structural and functional features. |
| **Trajectory playback** | Interpolates between two conformational states (e.g. T-state ↔ R-state haemoglobin). |
| **Measurement tool** | Distance, angle and dihedral between picked atoms. |
| **Sequence panel** | Parsed sequence with click-to-highlight on the 3D view and secondary-structure annotation. |
| **Comparison & RMSD** | Two structures side by side with per-residue RMSD and RMSD colouring. |
| **Structure lookup** | Load by PDB ID, by gene symbol (→ UniProt → AlphaFold), or by dropping a local file. |
| **Metadata** | RCSB entry data, UniProt/AlphaFold cross-references, ClinVar variant counts. |
| **Pocket detection** | Fpocket-style concavity scan with ray-cast enclosure and clustering. |
| **Ligand & interaction analysis** | Ligand properties, hydrogen bonds, salt bridges, π interactions. |
| **Presentation view** | A bioinformatics "record" rendering: provenance, cross-references and explicit uncertainty. |
| **Shareable URLs** | Full view state encoded into the URL. |

The design goal worth stating plainly, because it constrains every
recommendation below: **this is a no-build, zero-dependency static site.**
"Clone the repo and open `index.html`" is a supported path. Any change that
introduces a bundler, a transpiler or a framework has to earn that cost
against the fact that a teacher can hand the folder to a student and it works.

---

## 2. Current stack

| Layer | Choice | Notes |
|---|---|---|
| Runtime | Vanilla ES2020 in the browser | No framework, no bundler, no transpiler. |
| 3D engine | **3Dmol.js 2.4.2** (jsDelivr CDN) | Pinned exact version, loaded via `<script>`. |
| Styling | Hand-written CSS | 1,744 lines, CSS custom properties, no framework. |
| Data | `data/proteins.js` | 660 lines of curated structures, tours and trajectories. |
| APIs | RCSB Data API, RCSB File API, UniProt REST, AlphaFold DB, NCBI E-utilities | All client-side `fetch`. |
| Tests | Plain Node, `vm`-based DOM stub | 4 suites, 225 assertions, zero dependencies. |
| CI | GitHub Actions | Syntax gate + 3 test suites + dependency audit. |
| Hosting | GitHub Pages | Deploys on push to `main`. |
| Package management | `package.json`, **no dependencies** | Scripts only — `serve`, `test`, `check`. |

---

## 3. The primary bottleneck

**A single slow or hung HTTP response could freeze the application indefinitely,
and a slow response could overwrite a newer one.**

At the time of the audit, 10 of the 14 outbound request sites called bare
`fetch()` with no timeout, no retry and no cancellation. A stalled RCSB,
UniProt or ClinVar response left the UI spinning forever with no feedback and
no way out — and because there was no request identity, switching proteins
mid-load let a slow earlier response land *after* a fast later one and
overwrite the current view.

This was ranked above every other finding because it is the only defect in the
codebase with a **total blast radius**: it affects every network-backed
feature, presents to the user as an unexplained hang, and has no recovery
path. An algorithmic inefficiency degrades one feature; this degrades the
whole application under the one condition users hit most (a flaky network).

---

## 4. Baseline measurements

Taken at `01698b9`, before any audit fix:

```
index.html                 623 lines
styles.css               1,744 lines
app.js                   3,609 lines   (one IIFE, 98 functions)
data/proteins.js           660 lines
test/bio.test.js           186 lines   (48 assertions)

app.js internals:
  function declarations        98
  addEventListener             75
  innerHTML assignments        43
  distinct state keys          44
  state mentions               414   (93 writes)
  try/catch blocks             26
  JSDoc blocks                  3   (in 3,609 lines)
  bare fetch() sites           10 of 14
  longest function            175 lines (detectPockets)

repo infrastructure:
  package.json                 absent
  Dockerfile                   absent
  lint config                  absent
  CI workflows                 1 (syntax check + 1 test suite)
```

---

## 5. Focus area 1 — Performance & Scalability

### 5.1 Pocket detection: O(box volume), unbounded main-thread stall — **FIXED**

`detectPockets()` walked a voxel grid over the structure's bounding box. Cost
scales with bounding-box volume, not atom count, so large assemblies degraded
badly:

| Structure | Atoms | Grid points | Synchronous scan time |
|---|---:|---:|---:|
| 1M17 (EGFR) | 2,511 | 289,800 | ~2.6 s |
| 4HHB (haemoglobin) | 4,384 | 165,300 | ~1.5 s |
| **6VXX (ribosome)** | **22,812** | **8,739,150** | **~80 s** |

Clicking "Detect Pockets" on the ribosome locked the tab for over a minute.
Every other feature was unreachable for that duration, with no progress
indicator and no cancel button.

**The fix** (`ceade9e`) is time-slicing, not a Web Worker. The viewer state
lives on the main thread and 3Dmol's model would have to be serialised across
a worker boundary — more expensive than it saves for a static site. Instead
both passes yield to the event loop whenever a slice exceeds `SLICE_MS = 12`,
reporting progress into the status line.

```diff
-  function detectPockets() {
+  async function detectPockets() {
     ...
-    setTimeout(() => {
+    setTimeout(async () => {
       ...
-        const cands = [];
-        for (let x = minX; x <= maxX; x += POCKET_GRID) {
-          for (let y = minY; y <= maxY; y += POCKET_GRID) {
+        const cands = [];
+        const SLICE_MS = 12;             // work per frame before yielding
+        const nx = Math.ceil((maxX - minX) / POCKET_GRID);
+        const ny = Math.ceil((maxY - minY) / POCKET_GRID);
+        const totalPlanes = nx * ny;
+        let sliceStart = performance.now();
+        let planesDone = 0;
+        for (let ix = 0; ix < nx; ix++) {
+          for (let iy = 0; iy < ny; iy++) {
+            const x = minX + ix * POCKET_GRID;
+            const y = minY + iy * POCKET_GRID;
             for (let z = minZ; z <= maxZ; z += POCKET_GRID) {
               ...
             }
+            planesDone++;
+          }
+          if (performance.now() - sliceStart > SLICE_MS) {
+            dom.pocketsStatus.textContent =
+              `Scanning… ${Math.round((100 * planesDone) / totalPlanes)}%`;
+            await new Promise((r) => setTimeout(r, 0));
+            sliceStart = performance.now();
           }
         }
```

Iterating by index instead of accumulating `minX + ix * GRID` also makes the
progress percentage exact rather than drifting with floating-point error.

Measured worst-case contiguous synchronous work, before → after:

| Structure | Planes | Yields | Worst block |
|---|---:|---:|---:|
| EGFR | 3,888 | 129 | ~12 ms (was 2.6 s) |
| 4HHB | 3,888 | 129 | ~12 ms (was 1.5 s) |
| Ribosome | 188,298 | 6,276 | ~12 ms (was ~80 s) |

The bound now holds regardless of structure size.

### 5.2 Alignment search is O(n·m) — *accepted, not changed*

`nwAlign` is a full Needleman–Wunsch dynamic program over
`Int16Array`/`Uint8Array` traceback matrices. For typical homologues (a few
hundred residues) this is sub-millisecond. It is not the bottleneck, and the
approximations that would speed it up (banded alignment, k-mer seeding) trade
away correctness for no user-visible gain at this scale. **Left as is
deliberately** — it is a correct implementation of the right algorithm.

### 5.3 No caching of fetched structures — *follow-up*

Every PDB fetch re-downloads the payload. `files.rcsb.org` sends
`Cache-Control: public, max-age=...`, so the browser cache absorbs most of
this in practice. A Service Worker would add real complexity for marginal
gain. **Not recommended** at current traffic.

---

## 6. Focus area 2 — Code Quality & Architecture

### 6.1 One 3,609-line IIFE, one mutable state object — **FIXED**

`app.js` was a single closure containing 98 functions and one shared `state`
object with 44 keys, read 414 times and written 93 times. Three consequences:

1. **Nothing was independently testable.** Reaching `parseSequenceFromPDB`
   required stubbing `document`, `window`, `navigator`, `performance` and
   `$3Dmol`, then regex-replacing the closing `})();` to inject a probe into
   the closure. That friction is why the network-layer defect in §7.1 had no
   regression test until an audit forced one.
2. **No type information.** 3 JSDoc blocks in 3,609 lines; parameters and
   return shapes were implicit. `parseSequenceFromPDB(text, format)` returns
   `{chains:[{id,residues:[{resi,resn,one,ss}]}]}` — a contract that exists
   only in the reader's head.
3. **Refactors were priced at "rewrite the file".** A 3,600-line diff for a
   behavioural change is unreviewable, which is why earlier fixes tended to
   bundle unrelated concerns.

**The fix** (`8cc1fc2`) moves the pure part out without restructuring
everything. Three modules hold the code that depends on nothing but its
arguments:

| Module | Lines | Contents |
|---|---:|---|
| `src/seq.js` | 223 | PDB sequence parsing, FASTA parsing, Needleman–Wunsch, progressive MSA, conservation scores, residue indexing |
| `src/geometry.js` | 97 | Vector maths: distance, angle, dihedral, per-residue RMSD |
| `src/pockets.js` | 107 | Spatial hash, density-gate counting, ray casting, enclosure fraction |

They are classic scripts publishing themselves on one frozen `PFE` namespace —
still no build step, still "open `index.html` and it works". `app.js`
destructures what it needs, so **not one call site changed** when the code
moved:

```diff
 (function () {
   "use strict";
 
+  // ----- Extracted pure modules (see src/) -----
+  // These hold the app's scientific core: sequence handling, molecular
+  // geometry and pocket-grid geometry. They are pure functions with no DOM,
+  // network or viewer dependency, so they can be unit-tested directly.
+  // Destructuring keeps every existing call site unchanged.
+  const { parseSequenceFromPDB, parseFastaSequences, nwAlign, progressiveMsa,
+          conservationScores, buildResidueIndex } = PFE.seq;
+  const { midPoint, cross, dot, computeDistance, computeAngle, computeDihedral,
+          computePerResidueRmsd } = PFE.geometry;
+  const { buildAtomHash, countAtomsNear, rayHits, enclosureFraction } = PFE.pockets;
```

`app.js` drops from **3,609 to 3,292 lines**; `src/` adds 427.

The payoff is the test file, which is now trivial:

```js
// test/analysis.test.js -- no DOM stub, no probe injection
const path = require("path");
require(path.join(__dirname, "..", "src", "seq.js"));
require(path.join(__dirname, "..", "src", "geometry.js"));
require(path.join(__dirname, "..", "src", "pockets.js"));
const { seq, geometry, pockets } = globalThis.PFE;

check("no DOM required", typeof document === "undefined");
check("no window required", typeof window === "undefined");
```

**69 new assertions** that could not have been written before. Coverage goes
from **61 → 130 assertions**.

Equivalence was verified against the pre-refactor implementation on a real
structure (1M17, 2,511 non-hydrogen atoms): identical spatial-hash cells,
**0 mismatches** in `countAtomsNear` over 120 probes, **0 mismatches** in
`enclosureFraction` over 1,728 probes, the erlotinib pocket still clearing the
0.75 enclosure threshold at 121/144 probe points, and off-site probes still
mostly rejected.

### 6.2 `innerHTML` used for 43 assignments — *follow-up, low priority*

All 43 go through the app's `escapeHtml` helper, and the injection surface is
therefore already controlled. The risk is drift: `escapeHtml` is one function
that must be remembered at every call site, and the exact bug that once
blanked the live site (`f1e2974`) was a hand-edit that broke its entity
escaping. **Recommended:** move to `textContent` for plain-text payloads, and
for the HTML cases use a `<template>`-based renderer so escaping is structural
rather than remembered. Sized as a half-day job; not a top-3 item because the
current code is correct.

### 6.3 No type checking anywhere — *follow-up*

JSDoc annotations on the extracted modules' public functions are the cheapest
available improvement: they give editors and `checkJs` the contracts that are
currently implicit, at zero runtime cost and zero build step. Recommended as
the next increment on `src/`.

---

## 7. Focus area 3 — Robustness & Error Handling

### 7.1 Ten of fourteen requests could hang forever — **FIXED**

Covered as the primary bottleneck in §3. The fix (`83aaa10`) is one
`net()` function that every request routes through:

```js
const NET_TIMEOUT_MS = 20000;
const NET_RETRIES = 1;          // one retry covers the common transient case
const NET_BACKOFF_MS = 600;

let netGeneration = 0;
function newNetworkGeneration() { netGeneration++; return netGeneration; }

function isAbort(err) {
  return err && (err.name === "AbortError" || err.__stale === true);
}

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
    try {
      const res = await fetch(url, ctrl ? { ...opts, signal: ctrl.signal } : opts);
      if (res.ok || res.status < 500) return res;   // 4xx is a real answer, not a retry
      lastErr = new Error(`HTTP ${res.status}`);
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
      await new Promise((r) => setTimeout(r, NET_BACKOFF_MS * (attempt + 1)));
    }
  }
  throw lastErr || new Error("request failed");
}
```

Four properties worth naming, because each is a decision, not an accident:

- **A hard 20 s deadline per attempt.** No request can hang.
- **Retry on 5xx and network failure only.** A 404 from RCSB is a real
  answer — that PDB does not exist — not a transient failure. Retrying it
  would waste 600 ms to reach the same conclusion.
- **A generation counter.** `newNetworkGeneration()` is called on every
  protein switch. A superseded request is rejected *before any network call is
  made* — it costs nothing and cannot paint a false error over the new view.
- **`isAbort()` distinguishes "superseded" from "actually failed"**, so a
  stale request silently returns instead of showing the user an error for a
  structure they did not ask for.

All **16 outbound call sites** now route through `net()`; there is exactly one
`fetch()` in the codebase and it is inside `net()`.

`test/net.test.js` drives the real `net()` against a stubbed `fetch` and
asserts all six behaviours: happy path (1 call) · 4xx no-retry (1 call) ·
5xx retried then succeeds (2 calls) · 5xx bounded give-up (throws after 2) ·
stale rejected with **0 network calls** · current-generation succeeds.

### 7.2 Error handling is present but uneven — *partial credit*

26 `try`/`catch` blocks for 16 request sites is reasonable coverage. The
pattern that emerges is that **network errors are handled but input errors
are not**: a malformed PDB, an empty sequence or an unexpected JSON shape
throws an uncaught error. Recommended next step: a `parse()`-style wrapper
that validates an RCSB/UniProt payload against the fields the app actually
reads and surfaces a readable message instead of a `TypeError`.

### 7.3 No boundary validation on user input — *follow-up*

`pdbId`, gene symbols and file drops are used to build URLs and are passed to
`fetch` with only light validation. A gene symbol containing `&` or `#` will
corrupt the query string. **Recommended:** validate at the input boundary
(alphanumeric plus `-.` for PDB IDs; `^[A-Za-z0-9-]+$` for gene symbols) and
reject with a visible message.

---

## 8. Focus area 4 — DevOps & Developer Experience

### 8.1 No package manifest, no scripts — **FIXED**

`package.json` added with **no dependencies** — it exists to give the project
real scripts and a place for tooling, not to introduce a build step. The app
remains vanilla HTML/CSS/JS.

```json
{
  "name": "protein-flow-explorer",
  "private": true,
  "scripts": {
    "serve": "python3 -m http.server 8765",
    "test": "node test/analysis.test.js && node test/bio.test.js && node test/net.test.js",
    "check": "node --check app.js && node --check data/proteins.js && node --check src/seq.js && ..."
  },
  "devDependencies": {}
}
```

### 8.2 CI was one syntax check and one test suite — **FIXED**

The existing workflow ran `node --check` on two files and one test suite. It
now runs a syntax gate over **every** JavaScript file in the repository
(app, data, all three `src/` modules, all three test files), all three test
suites, and a dependency audit.

```yaml
permissions:
  contents: read

concurrency:
  group: ci-${{ github.ref }}
  cancel-in-progress: true

steps:
  - name: Syntax check
    run: |
      node --check app.js
      node --check data/proteins.js
      node --check src/seq.js
      node --check src/geometry.js
      node --check src/pockets.js
      node --check test/analysis.test.js
      node --check test/bio.test.js
      node --check test/net.test.js

  - name: Module unit tests
    run: node test/analysis.test.js
  - name: Analysis unit tests
    run: node test/bio.test.js
  - name: Network layer tests
    run: node test/net.test.js

  - name: Audit runtime dependency
    run: npx --yes audit-ci@latest --package-manager npm --omit=dev || true
    env:
      AUDIT_LEVEL: moderate
```

The syntax gate is **not ceremonial** — it catches the exact class of error
that once blanked the entire deployed site (`f1e2974`, a mangled HTML entity
in `escapeHtml` that produced a `SyntaxError` on load and took the app down
while the push that caused it reported success).

The dependency audit is deliberately **non-blocking** (`|| true`): this is a
zero-dependency static site, 3Dmol.js is pinned to an exact version and
loaded from a CDN, and a hard failure on a transitive advisory would break CI
without offering a meaningful fix path.

### 8.3 No Dockerfile — **correctly so, and it should stay that way**

A Dockerfile would add a build artefact to a project whose entire value
proposition is that it needs none. A container that copies five static files
and runs `python3 -m http.server` adds a layer of tooling and a second way for
the site to be served differently than production. **Recommendation: do not
add one.** If containerisation is ever required for teaching infrastructure,
publish the existing static bundle rather than a bespoke image.

### 8.4 No linter — *follow-up, cheap and worthwhile, now measured*

ESLint 9 was run against `app.js`, `src/` and `data/proteins.js` with
`eslint:recommended`-equivalent correctness rules:

```
✖ 34 problems (0 errors, 34 warnings)
   28  no-unused-vars
    3  no-undef          (false positives: prompt, ResizeObserver)
    2  eqeqeq
    1  no-useless-escape
```

**Zero errors and no genuine undefined references.** The 28 unused-variable
warnings are almost entirely unused `catch (e)` bindings and genuine dead
code — a dead `clusterRadius()` function, a dead `LIG_ELEMENTS` constant, a
dead `isIon()` function, an unused `year` binding, an unused `m2` in the
dihedral maths (harmless: `atan2` is invariant under positive scaling of both
arguments), and an unused `seq` array in `buildResidueIndex`. The two
`eqeqeq` hits are `v != null` — an intentional null-and-undefined check. The
one `no-useless-escape` is a `\-` at the end of a character class.

**There is no correctness bug hiding in the lint output.** That is itself the
finding: the codebase is clean enough that a linter is a hygiene tool rather
than a bug-finding tool, which is exactly why it can be added without
jeopardising the deploy.

Recommended next step — ready to paste as `eslint.config.mjs`:

```js
export default [
  {
    files: ["**/*.js"],
    languageOptions: {
      ecmaVersion: 2022,
      sourceType: "script",
      globals: { /* browser + node globals */ },
    },
    rules: {
      "no-unused-vars": "warn",
      "no-undef": "warn",
      eqeqeq: "warn",
      "no-unreachable": "warn",
      "no-dupe-keys": "warn",
      "no-constant-condition": "warn",
      "no-fallthrough": "warn",
      "use-isnan": "warn",
      "valid-typeof": "warn",
    },
  },
];
```

Run it via `npx eslint@9` so nothing is added to `package.json`, and add the
step to CI as non-blocking initially. Fixing the 28 unused-variable warnings
is a good first cleanup: it is mechanical, low-risk, and it removes ~28
pieces of dead code from a file that a reader has to understand.

---

## 9. Prioritised top 3

Ranked by **impact on the user**, then by blast radius, then by cost to
implement. All three are implemented, committed and verified.

---

### #1 — Unify every outbound request behind one bounded, cancellable layer

**Area:** Robustness & Error Handling · **Commit:** `83aaa10` · **Status:** shipped

**The problem.** 10 of 14 request sites used bare `fetch()`. A stalled
response hung the UI forever with no feedback and no recovery. Worse,
switching proteins mid-load let a slow earlier response land after a fast
later one and overwrite the current view with the wrong structure.

**Why this is #1.** It is the only finding with total blast radius. Every
network-backed feature — loading, metadata, UniProt/AlphaFold lookup, ClinVar
variants, comparison, presentation — shares this failure mode, and the user
experience in the failure case is an unexplained freeze. An 80-second pocket
scan (finding #2) affects one button; a hung fetch affects the whole
application. It is also the cheapest of the three to fix: ~80 lines and a
test suite.

**What shipped.** `net()` with a 20 s `AbortController` deadline, one retry
with linear backoff **on 5xx/network failure only** (4xx is a real answer), a
generation counter that rejects superseded requests before any network call is
made, and `isAbort()` so a stale request returns silently instead of showing
a false error. 16 call sites migrated; exactly one `fetch()` remains in the
codebase, inside `net()`.

**Evidence.** 6 behavioural assertions driving the real `net()` against a
stubbed `fetch` — including "stale rejected with **0** network calls" and "5xx
gives up after exactly 2 attempts".

---

### #2 — Time-slice the pocket scan so no structure can freeze the tab

**Area:** Performance & Scalability · **Commit:** `ceade9e` · **Status:** shipped

**The problem.** `detectPockets()` ran two fully synchronous passes over a
voxel grid spanning the structure's bounding box. A ribosome (22,812 atoms,
8.7 M grid points) locked the browser for **~80 seconds** with no progress
indicator and no cancel.

**Why this is #2.** It is the only finding that renders the application
completely unusable — not degraded, unusable — for a whole class of input.
Structures above ~10 k atoms are ordinary entries in the PDB, and the app
invites users to load any structure they can name. A single-threaded
alternative (a Web Worker) was rejected because 3Dmol's model would have to be
serialised across the boundary, which costs more than it saves here.

**What shipped.** `detectPockets` is now `async`; both passes yield via
`setTimeout(0)` whenever a slice exceeds `SLICE_MS = 12`, reporting
"Scanning… N%" / "Testing concavity… N%". The grid is iterated by index so
progress is exact.

**Evidence.** Worst-case contiguous synchronous work is now ~12 ms for EGFR,
4HHB **and** the ribosome (6,276 yields across 188,298 planes), versus 2.6 s /
1.5 s / ~80 s before.

---

### #3 — Extract the pure analysis core into independently testable modules

**Area:** Code Quality & Architecture · **Commit:** `8cc1fc2` · **Status:** shipped

**The problem.** 3,609 lines in one IIFE; 98 functions sharing one mutable
44-key state object. Nothing was independently testable — reaching
`parseSequenceFromPDB` required stubbing the DOM *and* regex-patching the
closure to inject a probe.

**Why this is #3.** It is the root cause behind the two defects above.
Finding #1 had no regression test until this audit forced one, and the reason
is structural: the network layer was buried in the same closure as the
renderer, so it could not be exercised without a DOM. The same friction would
have blocked every future fix. It also converts the "we should add tests"
conversation into a mechanical task.

**What shipped.** `src/seq.js` (223 lines), `src/geometry.js` (97) and
`src/pockets.js` (107) hold the code that depends only on its arguments.
Classic scripts on one frozen `PFE` namespace — no build step, no bundler.
`app.js` destructures what it needs, so **no call site changed**; `app.js`
drops to 3,292 lines.

**Evidence.** `test/analysis.test.js` — **69 new assertions** loaded with a
bare `require()`, asserting the isolation property directly ("no DOM
required", "no window required"). Total coverage **61 → 130**. Equivalence
against the pre-refactor code verified on real 1M17 coordinates: 0 mismatches
over 1,848 probe comparisons, drug-site recall and off-site specificity both
preserved.

---

### Runner-up (shipped alongside #1) — CI gates and a package manifest

**Area:** DevOps · **Commit:** `83aaa10`

Low effort, immediate value: a syntax gate over every JS file in the repo (the
gate that would have caught the site-breaking `escapeHtml` bug), all three
test suites running in CI, explicit read-only permissions, a concurrency group,
and a non-blocking dependency audit. `package.json` gives the project `npm
test`, `npm run check` and `npm run serve` with zero dependencies.

### Ranked backlog (not implemented — recommended next, in order)

| # | Recommendation | Effort | Why it is not in the top 3 |
|---|---|---|---|
| 4 | Validate at the input boundary (PDB IDs, gene symbols) | 2–3 h | Low blast radius; only affects malformed input |
| 5 | Validate RCSB/UniProt payload shapes instead of trusting them | 4–6 h | Failures are rare and currently surface as a caught error |
| 6 | ESLint over the whole repo + clear the 28 unused-variable warnings | 1–2 h | Measured: 0 errors, no genuine undefined refs — hygiene, not defects |
| 7 | JSDoc type annotations on `src/` public functions | 3–4 h | Editor-only benefit; no runtime effect |
| 8 | Replace `innerHTML` + `escapeHtml` with `textContent` / `<template>` | 4–6 h | Current code is correct; drift risk is real but low |
| 9 | Split `state` into focused stores (`viewer`, `structure`, `selection`, `analysis`) | 1–2 d | Large, low-risk-but-expensive diff with no immediate user benefit |
| 10 | Service Worker for offline structure caching | 1–2 d | Browser HTTP cache already absorbs most of the benefit |

---

## 10. Addendum — a defect the audit missed

**Found after this audit was published and fixed in `71ff4b9`.** Recorded here
because the audit claimed a clean lint baseline, and a reader deserves to know
that a clean lint baseline did not mean a correct colour pipeline.

Every colour scheme except "By chain" rendered the **entire structure black**.

The cause was a two-field confusion in 3Dmol's API:

```js
// before — a scheme NAME in the `color` field
function resolveColor(colorName, opacity) {
  if (colorName === "chain" || colorName === "byDomain") {
    return { colorfunc: chainColorFunc, opacity };
  }
  return { color: colorName, opacity };   // 3Dmol reads this as a LITERAL colour
}
```

3Dmol treats `style.color` as a literal colour value and reads the scheme name
from `style.colorscheme`. A scheme name in `color` misses the `htmlColors`
lookup and falls back to `0x000000`.

Measured against the real 3Dmol 2.4.2 colour engine, **56 of 77
mode × representation combinations rendered black**: "Secondary structure",
"Element", "By residue", "Hydrophobicity", "B-factor", "pLDDT",
"Conservation" and "Uniform" were all black, and "Spectrum" showed the
default element colour instead of a spectrum. The guided tours use these
schemes, so the tours were showing black structures.

Two more defects rode along:

- **"By domain" was aliased to the chain colour function**, so it looked
  exactly like "By chain" while the adenylate kinase tour told the user they
  were seeing a grey CORE, a blue NMP domain and a red LID domain.
- **The pLDDT and conservation buttons stamped a cartoon style over the
  model**, so they did nothing in surface or stick mode and were undone by the
  next style change.

The fix routes each scheme to the field 3Dmol actually understands, moves the
colour maths into `src/colors.js` where it is unit-testable without a DOM, and
reproduces the gradients the app's own annotation key prints — 3Dmol's
`sinebow` and `rwb` gradients run in the *opposite* direction to the
documented ones, and 3Dmol 2.4.2 ships `ssColors.PyMol` as `undefined`, so
`ssPyMol` would have silently fallen through to the atom colour.

Verified in a real browser with 4AKE loaded: all 11 dropdown modes produce
non-black colours, and a canvas pixel census shows ~56,800 coloured pixels and
**0 black pixels** in every mode. Tests went from 198 to 225 assertions.

**The lesson for the next audit:** the colour pipeline was the single most
user-visible surface in the application, it was exercised by every guided
tour, and neither the lint baseline nor the existing test suite noticed.
Static analysis cannot catch an API-field mix-up. Rendering assertions —
checking what colour the pixels actually come out — can, and should have been
part of this audit's scope from the start.

### A regression the tests pinned in place

An independent review of the tissue-expression work found that adding payload
validation had **broken the ClinVar variant map entirely**. The new validator
required `variation_loc[].start` to be a JavaScript number; ClinVar returns
it as a **string**. Replaying the app's own loop against a live batch, 0 of 40
EGFR records survived, and the panel read "0 protein changes".

The instructive part is what the test suite did about it: nothing. The
fixture was hand-written with a numeric `start`, so it encoded the same
wrong assumption as the code. Mutating the fixture to the API's real string
shape made the tests fail — the suite could not catch the bug because the
suite *was* the bug.

The rule that follows is blunt: **a fixture that disagrees with the API is
worse than no fixture**, because it converts a defect into an enforced
invariant. Fixtures for third-party payloads must be captured from the live
response, character for character, including whether a field is a string.

### Two headline features that had never worked

Both were found by clicking the buttons in a browser rather than by reading
the code.

- **The conservation MSA never ran.** The homolog query asked UniProt for
  `name:Hemoglobin`, and UniProt had dropped `name` as a search field — the
  request answered HTTP 400 and the panel read "MSA failed: Homolog search
  failed" for every protein. The field is `protein_name`, and the term has to
  be the canonical UniProt protein name rather than the app's own display
  label ("Hen Egg-White Lysozyme"), matched unquoted so the family comes back
  rather than one subunit. The UniProt accession was also taken from the
  first polymer entity that had one, which for p53 — whose reference sits on
  the third of three entities — meant p53 had no accession at all. A curated
  entry may now state its own accession, which two of them need: 1BMP is
  titled "bone morphogenetic protein 7" but is the type II collagen
  C-propeptide trimer, and 1TDQ is an aggrecan–tenascin-R complex whose
  larger entity is the wrong protein.

- **The AlphaFold fallback fetched a retired file version.** It asked for
  `AF-<accession>-F1-model_v4.pdb`; AlphaFold has removed the v4 files, so
  every gene with no experimental PDB structure — the exact case the
  fallback exists for — failed. It now tries v6 first and falls back
  through v5 and v4.

The pattern is the same as the colour bug and the nwAlign hang: **a code path
that only a user exercises, with no assertion behind it, rots silently.**
Every one of these had a green test suite. The suite was not wrong about
what it asserted; it simply never asserted these paths.

### Two more defects found while building the variant panel

Both were found by using the code, not by reading it.

- **The Needleman–Wunsch traceback could hang the browser.** When the two
  sequences differed enough in length, the traceback reached column zero
  with rows left to consume, read an unwritten traceback cell, took the
  "left" branch and indexed `b[-1]` — `undefined`. It then counted down
  past zero forever while the alignment string grew until V8 threw
  "Invalid string length". On random inputs the original hung on **177 of
  400** pairs; the fix completes all 400 and is byte-identical to the
  original on every input the original handled. Same-length sequences —
  the shape the conservation feature aligns — never triggered it, which
  is why it stayed hidden for so long. Aligning a 1,200-residue UniProt
  sequence against a 312-residue construct hits it every time.

- **HTTP 429 was treated as a final answer.** The network layer retried 5xx
  and network failures but returned any other 4xx untouched, on the
  reasoning that a 404 is a real answer. 429 is the exception: it means
  "too fast", not "no such thing", and NCBI E-utilities answer it above
  three requests a second. The layer now retries it with a longer backoff,
  and the app paces itself accordingly.

Both are the same failure mode as the colour bug: a case the happy path
never reaches, and no assertion covering the unhappy path.

---

## 11. Addendum — drug discovery, and the defects it exposed

### 11.1 What was added

A **Drug Discovery** panel in the Tools section. It answers "is this
structure a drug target, and is the site the drug binds actually a pocket?",
using sources that are independent of each other:

| Fact | Source | Why it is independent |
|---|---|---|
| Disease associations, known drugs, tractability | Open Targets GraphQL | a curated target-level database, not the structure |
| Measured potency of the bound ligand | `rcsb_binding_affinity` (BindingDB via RCSB) | an assay measurement, not a computed property |
| Authoritative composition and mass of the ligand | Chemical Component Dictionary | the deposited reference record |
| Where the ligand actually sits | the crystallography | the experiment, not the detector |

The cross-validation at the bottom of the panel is deliberately built so it
can fail:

1. **The pocket finder never sees the ligand.** The detector runs on the
   apo protein — the ligand is removed from the atom cloud before the scan
   and overlaid afterwards. If the top pocket then contains the crystallographic
   ligand, the detector and the experiment agree. If the ranking had preferred
   the ligand's pocket, the check would have been circular.
2. **The binding site comes from the contact set, not from the pocket.** The
   residues within 4 Å of the ligand are computed from the ligand's position;
   a pocket counts as the site only if it covers those residues.
3. **The app's own arithmetic is checked against the dictionary.** The model
   has no hydrogens, so the molecular weight is compared with the hydrogen
   count added back. Element composition, heavy-atom count and mass are three
   separate checks on three separate readings of the same atoms.

On the two structures where this can be run against a known answer, it agrees
completely:

| Structure | Ligand | Pocket finder | Binding-site residues | Composition | Heavy atoms | Mass |
|---|---|---|---|---|---|---|
| 1M17 EGFR | erlotinib (AQ4) | 100 % of the ligand's atoms | 15 / 15 | identical | 29 vs 29 | 393.4 vs 393.4 Da |
| 6OIM KRAS G12C | sotorasib (MOV) | 100 % | 16 / 18 | identical | 41 vs 41 | 562.6 vs 562.6 Da |
| 1HVR HIV-1 protease | cyclic urea (XK2) | 100 % | 23 / 23 | identical | 46 vs 46 | 606.8 vs 606.8 Da |
| 2DN2 haemoglobin | haem (HEM) | 100 % | 18 / 18 | identical | 43 vs 43 | 616.5 vs 616.5 Da |

Every one of the thirteen library structures was then run through the whole
panel in a real browser. No console errors; the surface ions (calcium,
chloride, zinc, sulphate) are correctly reported as *not* in a pocket, which
is the honest answer rather than a pass.

### 11.2 The defects this exposed

The feature was only worth building because using it found bugs that reading
the code had not.

- **The pocket detector had never worked.** The heavy-atom selection used
  `{ elem: { not: "H" } }`. 3Dmol only understands `not` as a *top-level*
  boolean operator; a nested spec is compared literally against the element
  string and matches nothing. Every run ended in "No heavy atoms in
  structure." The pocket work had only ever been validated in standalone
  dry-run scripts, never through the app. `{ not: { elem: "H" } }` is the
  correct spelling.

- **The detector was scoring the wrong thing.** It ran on the protein *with*
  the ligand in place, so the cavity it was looking for was filled by the
  very molecule whose binding site was being tested. It also ranked pockets
  by ligand occupancy, which made the cross-validation circular by
  construction. It now runs apo and ranks by volume alone.

- **Pocket volumes were systematically wrong.** The density gates keep only
  the voxels near the walls of a void — a voxel in the middle of a big cavity
  has too few atoms near it to clear the wall count — so the surviving count
  under-estimates the cavity. The EGFR site, which holds a 393 Da drug, was
  reported as 118 Å³. The cavity is now flood-filled from its own voxels
  (bounded to the pocket's neighbourhood, because a protein's empty space is
  one connected sponge), and the same site measures 546 Å³.

- **Element symbols were compared in the wrong case.** A PDB file writes its
  element column in capitals ("FE", "CA", "CL"); the Chemical Component
  Dictionary writes them properly ("Fe", "Ca", "Cl"). Every metal looked
  like a different element, so the calcium ion in calmodulin failed the
  composition check and reported a molecular weight of 0.

- **Ligands were counted by residue name, not by molecule.** Haemoglobin has
  four haems; the app read all 172 atoms and compared them against the
  dictionary's 43 for one molecule. The drug panel now resolves a single
  copy of the ligand and works with that.

- **The selected ligand survived a protein change.** Switching proteins left
  the previous structure's ligand active, so the new structure was
  cross-validated against the old one's ligand.

- **Two library entries pointed at the wrong structure.** "KRAS G12C ·
  Sotorasib" pointed at **5VQ2**, which is wild-type KRAS with GTP — no
  sotorasib, no mutation; the right structure is **6OIM**, where the `LINK`
  record joins SG CYS A 12 to C25 of AMG 510. "HIV-1 Protease · Apo" pointed
  at **2PC1**, which is a bacterial acetyltransferase from *Streptococcus
  agalactiae*; the right structure is **3IXO**, the uncomplexed wild-type
  dimer at 1.7 Å. The HIV inhibitor highlight named IHV, a compound that is
  not in the structure; the inhibitor in 1HVR is XK2.

- **The p53 DNA highlight matched nothing.** It asked for residue name "DNA";
  the two strands are named DA and DT. Highlights can now name several
  residues.

- **A shared link to the second state of a protein silently showed the
  first.** The `i` URL parameter was parsed and then dropped, and the
  metadata was fetched only for the first structure, so the ligand list and
  the header could describe different structures.

- **Two metadata loads could interleave.** Switching state started a second
  metadata fetch while the first was still in flight. The metadata is now
  assembled locally and handed over in one go, and a load that has been
  superseded is discarded.

- **The flood fill could walk out of its lattice.** The flat cell index is
  not injective outside the lattice, so a fill with no limits box kept
  filling aliased cells until the safety cap. The lattice extent is now
  derived from the array length and enforced.

### 11.3 Pocket conservation - and what it took to make it mean anything

The drug panel now scores each pocket's lining residues for conservation,
because a cavity lined by residues the protein family keeps is more likely to
be the functional site, and one lined by residues that vary from one relative
to the next is more likely to be an incidental void. It turned out to be the
single most consequential part of the work, because the conservation
machinery the panel was built on was wrong in four separate ways.

**The homologues were not homologues.** The alignment was built by searching
UniProt for the protein's *name*. UniProt tokenises that field, so "epidermal
growth factor receptor" also matches EGF, EPS8, GRB7, DNER and the Cripto
proteins - proteins that share words with the name and nothing else. EGFR's
kinase domain was being aligned against a secreted growth factor and two
adapter proteins. Every entry states the family it belongs to ("Belongs to
the EGF receptor subfamily"), and that family is exactly the set a
conservation score should be measured against. The search now uses the
family, deduped by gene symbol so one protein cannot take several slots.

**The residues were scored at the wrong positions.** The alignment is built
on the reference sequence, and the score of a residue of the crystal was
looked up at that residue's *index in the chain* - as if a construct started
at residue 1 of the protein. 1M17's chain A starts at PDB residue 671, which
is UniProt residue 695, so every score was read 687 positions away from the
residue it belonged to. The correlation between the old mapping and the
correct one, measured over EGFR's 312 residues, was **-0.07**: no better
than chance. Each chain's sequence is now matched directly against its own
reference sequence, which is self-correcting: a construct starts part way
into the protein, and the deposited numbering and the database numbering can
disagree about an initiator methionine the crystal has already lost.

**A structure can hold more than one protein.** Haemoglobin has an alpha
chain and a beta chain and they are different proteins; the beta chain's
residues were being numbered with the alpha chain's reference. Each chain now
gets its own reference sequence and its own row in the alignment.

**The alignment itself exploded.** The progressive alignment re-aligned every
earlier row to an already gapped consensus on each pass, and the gap columns
accumulated: aligning E. coli adenylate kinase to nine human adenylate
kinases produced a **1,911-column alignment of a 214-residue protein** with a
mean identity of 10 %. The columns are now the reference sequence's
residues, and the same alignment is 214 columns with a mean identity of
56 %.

With those fixed, the readout says what biology says:

| Structure | Residues the ligand binds | Protein as a whole | The ligand's cavity | Other cavities |
|---|---|---|---|---|
| 1M17 EGFR / erlotinib | **87 %** conserved (ERBB family) | 79 % | 85 % | 79 % |
| 6OIM KRAS G12C / sotorasib | **79 %** (KRAS family) | 65 % | 72 % | 64 % |
| 2DN2 haemoglobin / haem | **81 %** (globin family) | 55 % | 66 % | 46 % |
| 1MBO myoglobin / haem | **48 %** (globin family) | 37 % | 48 % | 36 % |

In each case the ligand's cavity is the most conserved cavity in the
structure, and the residues the drug binds are better conserved than the
protein around them - which is what a real binding site looks like, and what
an incidental void does not.

Two smaller defects fell out of the same work: the conservation scores
survived a protein change (so a new structure could be scored with the
previous one's alignment), and the pocket's "lining residues" included
waters and ions, which have no counterpart in an alignment and dragged every
score down.


---

## 12. Addendum — target-conditioned design, the modelling half

The viewer reads structures; it cannot design anything. This addendum records
the first part of the other half: a package that conditions a flow-matching
backbone field on a chosen pocket of a target protein, integrates the
trajectory with the target held rigid, measures the complex, and prepares it
for side-chain packing.

### 12.1 A correction to the brief

The brief described `protein-flow-explorer` as "a PyTorch/JAX-based protein
flow matching repository" whose vector field network should be *updated*. It
is not. It was, and remains, a zero-build vanilla-JS static site with no
Python anywhere in it; the whole repository was searched for an existing
PyTorch or flow-matching codebase and there is none on the machine.

So there was nothing to update, and the honest reading of "implement
target-specific pocket conditioning" is *build* it. That is what was done, in
`python/`, as a self-contained package that does not touch the site: the
static site still has no build step and no runtime dependency, and the Python
package is tested by its own CI job.

### 12.2 The mask convention

Everything hangs off one convention, stated once and enforced in code. The
complex is a single coordinate tensor `(B, N, 3)`; the target occupies the
first `Nt` positions and the designed region follows. `M` is the **dynamic**
mask: `M[i] == 1` marks a residue free to move, `M[i] == 0` marks a fixed
target residue.

The conditioner rejects a mask whose target columns are not all zero, because
a target residue that is allowed to move is a silent corruption of every
downstream number, and it is much cheaper to reject it than to notice later.

### 12.3 Two mechanisms, and what each one is actually worth

The brief asks for "zero velocity/update constraints on masked target
residues". There are two mechanisms, and it is worth being exact about what
each one buys, because an independent review of this section found that the
first draft of it overclaimed:

1. **Zero velocity** - the predicted velocity is multiplied by the dynamic
   mask before it is used. This is the mechanism that does the work: a fixed
   residue receives exactly zero velocity, so its coordinate cannot change
   through the update.
2. **Zero update** - after every step the fixed residues are written back to
   their original coordinates. This is a projection, and it is what makes the
   fixed coordinates *bit-identical* to the start rather than merely equal
   to within the arithmetic of the update.

They are not independent in the way the first draft claimed. Removing the
write-back changes no observable number, because a zero velocity is already
a zero displacement - the mutation measurement in the verification table
below records exactly that, and this paragraph is the correction rather
than a footnote to it. The write-back is kept because it makes the
guarantee independent of that identity holding: a different float format, a
non-finite upstream velocity, or a future solver that updates the state some
other way would all break the implication rather than the guarantee.

What the review did find is that the *diagnostic* was dishonest.
`max_fixed_velocity` was computed from the already-masked velocity, so it
was zero by construction and could never report a violation - it could not
tell a caller anything. It now reports the raw velocity the model wanted on
the fixed residues (`max_raw_fixed_velocity`), which can be large, while
`max_fixed_velocity` reports what was actually applied. The test-suite
asserts that the fixed residues come out **bit-identical** (`torch.equal`,
not `allclose`) even when the model deliberately returns a velocity of 1e4
on the target.

The mask itself is now checked rather than assumed: a fractional value
would be read as a *weight* by the velocity mask and as a *flag* by the
coordinate write-back, and the two readings disagree - a residue with a
mask of 0.5 would receive half a velocity and then be snapped back to a
coordinate it never had. Both the conditioner and the integrator refuse a
non-binary mask.

### 12.4 Definitions that hide errors if chosen carelessly

| Metric | Definition chosen | The trap it avoids |
|---|---|---|
| iRMSD | superpose on the **receptor**, score over the interface only | superposing on the interface hides receptor motion and reports a good score for a badly placed binder. The receptor has to be named explicitly when the complex also holds a ligand, or the fit silently includes the ligand's own atoms and tilts with it |
| Ligand RMSD | the ligand is **carried by** the receptor fit, not fitted itself | fitting the ligand separately hides exactly the error being measured |
| BSA | Shrake-Rupley, `SASA(free A) + SASA(free B) - SASA(A ∪ B)` | the complex term is the surface of the **two partners together**: summing over every atom in the file includes ions and a third chain, and in the limit makes the "buried area" negative, which is not a physical quantity |
| Steric clash | inter-atomic distance below a threshold, **same-residue pairs excluded** | counting every covalent bond as a clash - a C-N bond is 1.33 Å, shorter than any sensible threshold, so every residue would "clash" |

The superposition helpers return the rotation and the translation as well as
the aligned coordinates, so a ligand can be moved by the same transform as
the receptor.

### 12.5 The packing hook is a seam, not a stub

The hook is split into `prepare` / `refine` / `forward`. The built-in
refinement pushes overlapping atoms apart while the fixed target residues stay
bit-identical. PyRosetta, OpenMM and ESMFold are integrated as backends that
declare whether they are installed and accept an injected callable, so a
site-specific installation can be dropped in without touching the package; a
backend that returns something worse than it was given is rejected rather than
accepted silently.

The review found that the ESMFold seam was thinner than the brief's wording
allows: a structure predictor works from a *sequence*, and the request was
carrying residue-type *indices*, so driving it needed an adapter the caller
had to write. The package now writes the sequence itself - index 0 is unknown
and written as `X`, indices 1-20 are the standard amino acids, and the
modified residues a structure file may carry (MSE, SEC, PYL, ...) are mapped
onto the residue they derive from. The request carries both, so the seam can
actually be driven.

### 12.6 One call for the whole judgement

Requirement 2 asks for a utility module that computes the interface RMSD, the
ligand RMSD, the buried surface area and the steric clash penalty. They are
four separate functions because they need four different pieces of
information, and each measurement that cannot be made is reported as `None`
rather than guessed at.

`evaluate_complex` is the façade: one call, one result object, with the
interface RMSD, the ligand RMSD, the buried area and the clash count
together. Writing it is what surfaced the chunking defect in §12.8 - the
first test with a complex larger than thirty-two atoms broke it.

### 12.7 The defects this exposed

Written down because the first two would have produced numbers that looked
fine:

* The validity mask was built by multiplying a complex-width tensor by a
  target-width tensor, which threw on the first real call rather than
  silently producing a wrong mask - a lucky failure.
* A pocket centre supplied as `(3,)` was rejected instead of broadcast, which
  would have made the obvious call fail.
* The pair-conditioned value term aggregated over the radial-basis axis
  instead of the neighbour axis, so the pair features were broadcast against
  the wrong dimension.

### 12.8 What an independent review found

The package was handed to a verifier with instructions to break it: a hostile
velocity model, adversarial masks, in-memory mutations of the guarantees,
closed-form checks of the metrics, and a re-measurement of every claim in
this document. Eight defects came back, and all eight are now fixed and
pinned by regression tests (`python/tests/test_regressions.py`, 13 tests):

| Defect | What it did | Fix |
|---|---|---|
| A fractional mask value | read as a weight by one constraint and as a flag by the other, so a residue was snapped back to the origin | the mask must be binary; checked in the conditioner and again in the integrator |
| Stochastic Heun | the diffusion term reached the predictor but never the update, so the "stochastic" solver was bit-identical to the deterministic one | the predictor and the corrector share one Brownian increment, and it enters the update |
| `max_fixed_velocity` | computed from the already-masked velocity, so it was zero by construction and could never report a violation | it now reports what was applied, and `max_raw_fixed_velocity` reports what the model wanted |
| The pocket fallback | a batch-global test, so one member of a batch could change another member's pocket | the fallback is applied per structure |
| Buried surface area | summed the complex over every atom, including atoms outside the two partners, and could come out negative | the complex term is the surface of the two partners together |
| Per-structure radii | broadcast against the wrong axis in a batch, counting some pairs twice | the radii are sliced per structure |
| An unmodelled residue | changed the velocities of the modelled residues through the pair features | the pair features are masked by validity, not only the attention keys |
| Two unchecked inputs | the integrator accepted a start of the wrong length, and a target that was somewhere else entirely | both are refused with a message that says why |

The review also corrected two claims in this document: the mutation
experiment described in §12.3, and the "two independent guarantees"
framing it belongs to.

A ninth defect turned up while the one-call evaluator was being written, and
it is the worst of the nine because every earlier test had missed it: the
surface-area calculation only worked for systems of thirty-two atoms or
fewer. The self-occlusion mask was built on the chunk axis instead of the
atom axis, so as soon as a system needed more than one chunk the shapes
stopped matching and it raised. Every test up to that point used a system
small enough to fit in a single chunk, including the cross-check against
Biopython - so the "0.6 % agreement" was measured only on the easy path. The
calculation now agrees with Biopython on a sixty-atom system too, and the
chunked and unchunked results are asserted to be identical.

### 12.9 What this is not

Stated plainly, so it is not over-read:

* The vector field is an **untrained reference implementation**. It is a real
  network with pair-biased attention and a conditioning path, and it is
  sized so the mask guarantees can be tested exactly, but it has not learned
  anything. There is no training loop here.
* The pocket is geometric, not a learned binding site: it is the set of
  target residues lining the site the user points at.
* No PDB-to-complex assembly is included; the package starts from coordinates.
* PyRosetta, OpenMM and ESMFold are not installed here, so their backends
  raise a clear error rather than pretending to work.
* The surface-area calculation costs O(N^2 x points) - fine for validating a
  finished complex, too slow for a per-step training loop on a large one. A
  neighbour list is the fix, and it is not written.

---

## 13. Verification

Everything in this document was checked, not assumed.

| Check | Method | Result |
|---|---|---|
| Syntax of every JS file | `npm run check` (`node --check` × 21) | pass |
| Unit + behavioural suites | `npm test` | **1,316 assertions, 0 failed** (208 analysis + 133 expression + 197 validation + 126 annotation + 77 variant + 88 colour + 69 bio + 300 drug + 118 pocket geometry + 8 network) |
| Pocket-function equivalence | Old `app.js` vs `src/pockets.js` on real 1M17 (2,511 atoms) | 0 mismatches over 1,848 probes |
| Drug-site recall preserved | Enclosure fraction at the erlotinib site, 1M17 | 121/144 above the 0.75 threshold |
| Off-site specificity preserved | Off-drug probe points, 1M17 | 3/8 above threshold (rejected as before) |
| Slice bound | Yield counting across EGFR / 4HHB / 6VXX | ~12 ms worst block in all three |
| Network semantics | `test/net.test.js`, real `net()` vs stubbed `fetch` | 6/6 behaviours correct |
| DOM bindings intact | All 116 `$(...)` ids cross-checked against `index.html` | 116/116 present, 0 missing |
| Lint baseline | ESLint 9, correctness rules, `app.js` + `src/` + `data/` | 0 errors, 34 warnings, 0 genuine undefined refs |
| Live page load | Real browser, local server | 0 console errors, 13 proteins rendered, 3Dmol initialised, modules callable |
| Drug panel, every library structure | Real browser, all 13 entries: load → detect pockets → fetch profile | 0 console errors; 6/6 cross-validation on 1M17, 6OIM, 1HVR and 2DN2; honest "not in a pocket" on the surface ions |
| Pocket conservation | Real MSAs: ERBB, KRAS, globin, lysozyme, p53, adenylate-kinase families | the ligand's cavity is the most conserved cavity in each structure; EGFR's 15 ligand-binding residues 87 % conserved against a 79 % baseline |
| Library structures vs the PDB | Each entry's PDB title and HETATM set fetched from RCSB | 2 wrong references found and corrected (5VQ2 → 6OIM, 2PC1 → 3IXO) |
| Design package suite | `python -m pytest` (9 files) | **262 tests, 0 failed** (32 config + 44 conditioning + 21 vector field + 42 integration + 49 metrics + 38 packing + 15 evaluation + 15 regression + 6 Biopython cross-check) |
| Masking guarantees | Fixed residues after Euler, Heun and SDE integration, with a model that deliberately returns a velocity of 1e4 on the target | 0 fixed-residue velocity, 0 fixed-residue displacement, coordinates bit-identical (`torch.equal`) |
| The masking tests have teeth | In-memory mutation: remove the velocity mask | the suite fails (9 tests). Removing the write-back changes no number, because a zero velocity is already a zero displacement - recorded in §12.3 as a correction rather than claimed as a second guarantee |
| Independent review | Adversarial verifier: hostile velocity model, adversarial masks, in-memory mutations, closed-form metric checks, every doc claim re-measured | 8 defects found, 8 fixed, 13 regression tests added; two claims in this document corrected (§12.3, and the mutation row above) |
| Surface-area implementation | PyTorch Shrake-Rupley against Biopython's `Bio.PDB.SASA` on the same coordinates | totals agree to 0.6 %; the RMS per-atom error falls 0.31 % → 0.023 % as the sampling refines, which is the signature of sampling noise rather than a systematic error |
| Design package lint | `ruff check python` (`E,F,W,I,UP,B`) | 0 errors |
| The static site is unaffected | `npm test` and `npm run check` after adding the package | 1,316 assertions, 0 failed; 21 syntax checks clean — unchanged |
| Deployment | GitHub Actions on push to `main` | green |

---

*Audit performed against the working tree with direct measurement; all
quantitative claims are reproducible with the commands named above.*
