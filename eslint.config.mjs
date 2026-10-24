/**
 * ESLint configuration for Protein Flow Explorer.
 *
 * The project is deliberately dependency-free at runtime, so this config adds
 * no npm dependency: run it with `npx eslint@9 .` and it is fetched on
 * demand. The globals are listed explicitly rather than pulled from the
 * `globals` package for the same reason.
 *
 * Only correctness rules are enabled here. Style is left alone -- this
 * codebase's formatting is deliberate and a style linter would be noise.
 */

// Browser globals used by the app.
const browserGlobals = {
  window: "readonly", document: "readonly", navigator: "readonly",
  location: "readonly", history: "readonly", performance: "readonly",
  fetch: "readonly", console: "readonly", alert: "readonly", prompt: "readonly",
  confirm: "readonly", setTimeout: "readonly", clearTimeout: "readonly",
  setInterval: "readonly", clearInterval: "readonly",
  requestAnimationFrame: "readonly", cancelAnimationFrame: "readonly",
  requestIdleCallback: "readonly", cancelIdleCallback: "readonly",
  queueMicrotask: "readonly", structuredClone: "readonly",
  URL: "readonly", URLSearchParams: "readonly", URLPattern: "readonly",
  Blob: "readonly", File: "readonly", FileReader: "readonly",
  FormData: "readonly", Headers: "readonly", Request: "readonly",
  Response: "readonly", AbortController: "readonly", AbortSignal: "readonly",
  TextDecoder: "readonly", TextEncoder: "readonly", crypto: "readonly",
  Image: "readonly", XMLHttpRequest: "readonly", Event: "readonly",
  CustomEvent: "readonly", MouseEvent: "readonly", KeyboardEvent: "readonly",
  ResizeObserver: "readonly", IntersectionObserver: "readonly",
  MutationObserver: "readonly", DOMParser: "readonly", Node: "readonly",
  localStorage: "readonly", sessionStorage: "readonly",
  getComputedStyle: "readonly", matchMedia: "readonly",
  Intl: "readonly", WebGLRenderingContext: "readonly",
};

// Node globals, used by the test harness.
const nodeGlobals = {
  console: "readonly", process: "readonly", Buffer: "readonly",
  __dirname: "readonly", __filename: "readonly", global: "writable",
  module: "writable", exports: "writable", require: "readonly",
  setTimeout: "readonly", clearTimeout: "readonly",
  setInterval: "readonly", clearInterval: "readonly",
  queueMicrotask: "readonly", structuredClone: "readonly",
  TextDecoder: "readonly", TextEncoder: "readonly", URL: "readonly",
  URLSearchParams: "readonly", fetch: "readonly", performance: "readonly",
  AbortController: "readonly", AbortSignal: "readonly",
};

// Globals the app itself defines for 3Dmol.js and the curated dataset.
const appGlobals = {
  $3Dmol: "readonly", PROTEINS: "readonly", TOURS: "readonly",
  TRAJECTORIES: "readonly", PFE: "readonly",
};

export default [
  {
    ignores: ["node_modules/**", ".git/**", "**/*.min.js"],
  },
  {
    // The app, the data and the browser-side modules.
    files: ["app.js", "src/**/*.js", "data/**/*.js"],
    languageOptions: {
      ecmaVersion: 2022,
      sourceType: "script",
      globals: { ...browserGlobals, ...appGlobals },
    },
    rules: {
      // --- correctness ---------------------------------------------------
      "no-undef": "error",
      "no-unused-vars": ["warn", { args: "none", caughtErrors: "none", varsIgnorePattern: "^_" }],
      "no-const-assign": "error",
      "no-dupe-args": "error",
      "no-dupe-class-members": "error",
      "no-dupe-keys": "error",
      "no-duplicate-case": "error",
      "no-dupe-else-if": "error",
      "no-fallthrough": "error",
      "no-irregular-whitespace": "error",
      "no-self-assign": "error",
      "no-sparse-arrays": "error",
      "no-unreachable": "error",
      "no-unsafe-finally": "error",
      "no-unsafe-negation": "error",
      "use-isnan": "error",
      "valid-typeof": "error",
      eqeqeq: ["warn", "smart"],
      "no-cond-assign": "error",
      "no-constant-condition": ["error", { checkLoops: false }],
      "no-empty": ["warn", { allowEmptyCatch: true }],
      "no-func-assign": "error",
      "no-obj-calls": "error",
      "no-sparse-arrays": "error",
      "no-prototype-builtins": "error",
      "no-script-url": "error",
      "no-shadow-restricted-names": "error",
      "no-template-curly-in-string": "error",
      "no-unmodified-loop-condition": "error",
      "no-unreachable-loop": "error",
      "no-useless-backreference": "error",
      "no-useless-catch": "error",
      "no-useless-escape": "warn",
      "no-useless-rename": "error",
      "no-with": "error",
      "require-yield": "error",
      yoda: "warn",
    },
  },
  {
    // The curated dataset is a classic script that also exports itself for
    // the Node test harness.
    files: ["data/**/*.js"],
    languageOptions: {
      ecmaVersion: 2022,
      sourceType: "script",
      globals: { ...browserGlobals, ...appGlobals, module: "writable", exports: "writable" },
    },
    rules: {
      "no-undef": "error",
      "no-unused-vars": ["warn", { args: "none", caughtErrors: "none", varsIgnorePattern: "^_" }],
      "no-dupe-keys": "error",
      "no-unreachable": "error",
    },
  },
  {
    // The test harness runs under Node.
    files: ["test/**/*.js"],
    languageOptions: {
      ecmaVersion: 2022,
      sourceType: "script",
      globals: { ...nodeGlobals, ...appGlobals },
    },
    rules: {
      "no-undef": "error",
      "no-unused-vars": ["warn", { args: "none", caughtErrors: "none", varsIgnorePattern: "^_" }],
      "no-const-assign": "error",
      "no-dupe-keys": "error",
      "no-duplicate-case": "error",
      "no-fallthrough": "error",
      "no-self-assign": "error",
      "no-sparse-arrays": "error",
      "no-unreachable": "error",
      "no-useless-escape": "warn",
      "use-isnan": "error",
      "valid-typeof": "error",
      eqeqeq: ["warn", "smart"],
      "no-empty": ["warn", { allowEmptyCatch: true }],
    },
  },
];
