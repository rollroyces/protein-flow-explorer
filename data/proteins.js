/**
 * Curated protein dataset with educational context.
 * Each protein entry contains:
 *   id          - unique key
 *   name        - display name
 *   pdbIds      - array of PDB IDs (multiple = conformational states)
 *   stateLabels - label per PDB ID (e.g., "T-state (deoxy)")
 *   summary     - one-line role
 *   description - 2-3 sentence educational description
 *   category    - "folding" | "domains" | "binding" | "allostery" | "drug-target"
 *   color       - accent color for cards
 *   highlights  - { chain, residues } for important functional sites
 */

const PROTEINS = [
  {
    id: "hemoglobin",
    name: "Hemoglobin",
    pdbIds: ["2DN2", "1HHO"],
    stateLabels: ["T-state (deoxy)", "R-state (oxy)"],
    summary: "Oxygen transporter in red blood cells",
    description:
      "Hemoglobin is a tetramer of two α and two β chains. Each subunit carries a heme group with an iron atom that reversibly binds O₂. The T→R transition is the textbook example of allosteric cooperativity.",
    category: "allostery",
    color: "#ef4444",
    highlights: {
      chains: ["A", "B", "C", "D"],
      residues: [
        { chain: "A", resi: 87, label: "α His87 (proximal)" },
        { chain: "A", resi: 58, label: "α His58 (distal)" },
      ],
    },
  },
  {
    id: "lysozyme",
    name: "Hen Egg-White Lysozyme",
    pdbIds: ["1LYZ"],
    stateLabels: ["Native"],
    summary: "Classic antibacterial enzyme",
    description:
      "Lysozyme cleaves the peptidoglycan layer of bacterial cell walls. It has a deep active-site cleft containing Glu35 and Asp52, which act as the general acid/base catalysts. It is the model system for studying protein folding and stability.",
    category: "binding",
    color: "#10b981",
    highlights: {
      residues: [
        { resi: 35, label: "Glu35 (acid catalyst)" },
        { resi: 52, label: "Asp52 (nucleophile)" },
        { resi: 62, label: "Trp62 (substrate stacker)" },
      ],
    },
  },
  {
    id: "hiv-protease",
    name: "HIV-1 Protease",
    pdbIds: ["1HVR", "2PC1"],
    stateLabels: ["With inhibitor", "Apo (open)"],
    summary: "Antiretroviral drug target",
    description:
      "HIV-1 protease is a homodimeric aspartic protease. Its C2-symmetric active site uses two Asp25 residues (one from each monomer) to cleave the viral polyprotein. Most approved HIV drugs target this enzyme.",
    category: "drug-target",
    color: "#8b5cf6",
    highlights: {
      residues: [
        { resi: 25, label: "Asp25 (catalytic)" },
        { chain: "A", resi: 25, label: "Asp25 A" },
        { chain: "B", resi: 25, label: "Asp25 B" },
      ],
      ligand: { resn: "IHV", label: "Inhibitor (drug)" },
    },
  },
  {
    id: "adenylate-kinase",
    name: "Adenylate Kinase",
    pdbIds: ["4AKE", "1AKE"],
    stateLabels: ["Open (inactive)", "Closed (active)"],
    summary: "Phosphoryl transfer enzyme with large domain motion",
    description:
      "Adenylate kinase balances cellular ATP/ADP/AMP by transferring a phosphoryl group. It undergoes one of the largest known substrate-induced domain closures: the LID and NMP domains close over the bound nucleotides during catalysis.",
    category: "allostery",
    color: "#f59e0b",
    highlights: {
      residues: [
        { resi: 10, label: "Gly10 (P-loop)" },
        { resi: 30, label: "Lys30 (catalytic)" },
      ],
    },
  },
  {
    id: "myoglobin",
    name: "Sperm Whale Myoglobin",
    pdbIds: ["1MBO"],
    stateLabels: ["Oxy"],
    summary: "Oxygen storage in muscle",
    description:
      "Myoglobin is a monomeric oxygen-storage protein with a single heme group. It was the first protein whose structure was solved by X-ray crystallography (Kendrew, 1958). It serves as a model for studying heme-protein interactions and folding.",
    category: "binding",
    color: "#3b82f6",
    highlights: {
      residues: [
        { resi: 64, label: "His64 (distal)" },
        { resi: 93, label: "His93 (proximal)" },
      ],
      ligand: { resn: "HEM", label: "Heme group" },
    },
  },
  {
    id: "calmodulin",
    name: "Calmodulin",
    pdbIds: ["1CLL", "1CDL"],
    stateLabels: ["Extended (Ca²⁺-free)", "Compact (target-bound)"],
    summary: "Calcium sensor with dramatic conformational change",
    description:
      "Calmodulin has two globular domains connected by a flexible α-helix. Upon binding calcium, the helices reorganize and the protein wraps around target kinases, phosphatases, and channels, regulating many calcium-dependent pathways.",
    category: "allostery",
    color: "#ec4899",
  },
  {
    id: "insulin",
    name: "Insulin Hexamer",
    pdbIds: ["4INS"],
    stateLabels: ["T-state storage"],
    summary: "Hormone that regulates blood glucose",
    description:
      "Insulin is stored in the pancreas as a zinc-coordinated hexamer of three dimers. Active insulin is a dimer of one A-chain and one B-chain connected by disulfide bonds. The hexameric form protects the active dimer from degradation.",
    category: "domains",
    color: "#06b6d4",
    highlights: {
      ligand: { resn: "ZN", label: "Zinc ion" },
    },
  },
  {
    id: "t4-lysozyme",
    name: "T4 Lysozyme Mutants",
    pdbIds: ["1L63"],
    stateLabels: ["Wild-type-like"],
    summary: "Model system for protein stability and folding",
    description:
      "T4 lysozyme has been used extensively as a model for studying how point mutations affect protein stability, dynamics, and folding. The Matthews laboratory has solved hundreds of mutant structures to map the folding energy landscape.",
    category: "folding",
    color: "#84cc16",
  },
];

/**
 * Pre-defined tours: each tour is a sequence of steps that animate the viewer.
 * Step types:
 *   camera      - { pdbIndex, zoom, spin }  switches to a state and frames camera
 *   represent   - { style, colorScheme, opacity }  change rendering style
 *   highlight   - { atoms: "selector", color, label }  spotlight a feature
 *   label       - { html, position }  add a 2D text callout on screen
 *   show        - { atoms }  add a sticks/highlight representation
 *   hide        - removes highlights/labels
 */
const TOURS = {
  hemoglobin: [
    {
      title: "Hemoglobin: Oxygen's Molecular Shuttle",
      steps: [
        {
          title: "Quaternary Structure",
          body:
            "Hemoglobin is a tetramer of two α chains (red, green) and two β chains (yellow, blue). The four subunits work cooperatively — when one binds O₂, the others become more receptive.",
          action: { pdbIndex: 0, style: "cartoon", colorScheme: "chain", spin: true },
          duration: 8000,
        },
        {
          title: "The Heme Cofactor",
          body:
            "Each chain cradles a heme group — a porphyrin ring with an iron atom at its center. The heme is the actual oxygen-binding site. We will zoom into one subunit to see it.",
          action: { pdbIndex: 0, style: "cartoon", colorScheme: "chain", zoomChain: "A", focusHeme: true, spin: false },
          duration: 9000,
        },
        {
          title: "The Iron-Oxygen Bond",
          body:
            "The iron atom is held by His F8 ('proximal') below, and His E7 ('distal') gates access from above. In the oxy state, O₂ binds to the iron between the two histidines.",
          action: { pdbIndex: 1, style: "ball+stick", colorScheme: "byElement", focusHeme: true, highlight: { resn: "HEM" } },
          duration: 10000,
        },
        {
          title: "T → R Conformational Change",
          body:
            "When O₂ binds, the iron moves into the heme plane, pulling the F-helix with it. This shifts the entire quaternary structure from the tense (T) to the relaxed (R) state — Cooperativity at the molecular level.",
          action: { pdbIndex: 1, style: "cartoon", colorScheme: "spectrum", spin: true },
          duration: 9000,
        },
      ],
    },
  ],

  lysozyme: [
    {
      title: "Lysozyme: How an Antibacterial Enzyme Works",
      steps: [
        {
          title: "The α/β Fold",
          body:
            "Lysozyme has a compact α+β fold: a small N-terminal β-sheet domain and a larger C-terminal α-helical domain. The deep cleft between them is the active site.",
          action: { style: "cartoon", colorScheme: "ss" },
          duration: 7000,
        },
        {
          title: "The Active Site Cleft",
          body:
            "The cleft binds six consecutive sugar units of the bacterial peptidoglycan. The colored spectrum shows the binding groove running across the protein.",
          action: { style: "cartoon", colorScheme: "spectrum", highlight: { resi: [35, 52, 62, 63, 101, 34] } },
          duration: 8000,
        },
        {
          title: "The Catalytic Dyad",
          body:
            "Glu35 (acid) and Asp52 (nucleophile) sit on opposite sides of the cleft. They work together to cleave the glycosidic bond — the same chemistry used by every lysozyme from phage to humans.",
          action: { style: "cartoon", colorScheme: "ss", showSticks: { resi: [35, 52] } },
          duration: 9000,
        },
        {
          title: "Substrate Binding Mode",
          body:
            "The enzyme distorts its substrate toward the transition-state geometry. Trp62 and Trp63 stack against the sugar rings, holding them in place for catalysis.",
          action: { style: "cartoon", colorScheme: "ss", showSticks: { resi: [35, 52, 62, 63, 101] } },
          duration: 8000,
        },
      ],
    },
  ],

  "hiv-protease": [
    {
      title: "HIV-1 Protease: A Master Drug Target",
      steps: [
        {
          title: "A C2-Symmetric Dimer",
          body:
            "HIV-1 protease is a homodimer — two identical subunits (cyan, magenta) related by a 2-fold axis. The active site sits exactly on this symmetry axis.",
          action: { style: "cartoon", colorScheme: "chain" },
          duration: 7000,
        },
        {
          title: "The Dimer Interface",
          body:
            "Most of the protein is buried at the dimer interface. Drug resistance mutations often cluster here — they destabilize the interface or reshape the active site.",
          action: { style: "cartoon", colorScheme: "chain", surface: { opacity: 0.35 } },
          duration: 8000,
        },
        {
          title: "Two Aspartates, One Active Site",
          body:
            "Asp25 from each monomer forms the catalytic dyad. They activate a water molecule that attacks the peptide bond of the substrate.",
          action: { style: "cartoon", colorScheme: "ss", showSticks: { resi: [25] } },
          duration: 9000,
        },
        {
          title: "Drug Binding",
          body:
            "The bound inhibitor (shown as spheres) mimics the transition state of the substrate. Modern HIV drugs like darunavir are engineered to fit tightly and overcome resistance.",
          action: { style: "cartoon", colorScheme: "ss", showSpheres: { resn: "IHV" } },
          duration: 9000,
        },
      ],
    },
  ],

  "adenylate-kinase": [
    {
      title: "Adenylate Kinase: The Molecular Hinge",
      steps: [
        {
          title: "Three Domains",
          body:
            "Adenylate kinase has a CORE (gray), a NMP-binding domain (blue), and a LID domain (red). The substrate binds at the interface of all three.",
          action: { pdbIndex: 0, style: "cartoon", colorScheme: "byDomain", spin: true },
          duration: 8000,
        },
        {
          title: "The Open, Inactive State",
          body:
            "Without substrate, the enzyme is open and flexible. The active site cleft is exposed to solvent — it cannot catalyze phosphoryl transfer in this state.",
          action: { pdbIndex: 0, style: "cartoon", colorScheme: "byDomain", spin: false },
          duration: 7000,
        },
        {
          title: "Induced Fit",
          body:
            "When substrate binds, the NMP and LID domains close over it. The closure excludes water and positions the substrates precisely for phosphoryl transfer.",
          action: { pdbIndex: 1, style: "cartoon", colorScheme: "byDomain" },
          duration: 8000,
        },
        {
          title: "The Conformational Movie",
          body:
            "Toggle between open and closed states to see the magnitude of this motion. The NMP domain rotates ~40° and the LID ~70° — among the largest domain movements known.",
          action: { pdbIndex: 1, style: "cartoon", colorScheme: "byDomain" },
          duration: 7000,
        },
      ],
    },
  ],

  myoglobin: [
    {
      title: "Myoglobin: The First Protein Structure",
      steps: [
        {
          title: "Eight α-Helices",
          body:
            "Myoglobin was the first protein whose 3D structure was solved by X-ray crystallography (Kendrew et al., 1958, Nobel Prize 1962). It is built from 8 α-helices (A–H) packed around a central heme.",
          action: { style: "cartoon", colorScheme: "ss" },
          duration: 8000,
        },
        {
          title: "The Heme Pocket",
          body:
            "The heme group sits in a hydrophobic pocket. His F8 (proximal) coordinates the iron from below; His E7 (distal) provides a hydrogen bond to bound O₂.",
          action: { style: "cartoon", colorScheme: "ss", showSticks: { resn: "HEM" } },
          duration: 8000,
        },
        {
          title: "The Globin Fold",
          body:
            "The 3-over-3 α-helical sandwich ('globin fold') is shared by hundreds of proteins, including hemoglobin. It is one of the most ancient and successful protein designs in nature.",
          action: { style: "surface", colorScheme: "ss", surface: { opacity: 0.55 } },
          duration: 8000,
        },
      ],
    },
  ],

  calmodulin: [
    {
      title: "Calmodulin: A Calcium-Triggered Switch",
      steps: [
        {
          title: "Two Lobes, One Hinge",
          body:
            "Calmodulin has an N-terminal lobe (blue) and a C-terminal lobe (red), each with two EF-hand calcium-binding sites, connected by a long flexible α-helix.",
          action: { pdbIndex: 0, style: "cartoon", colorScheme: "byDomain", spin: true },
          duration: 8000,
        },
        {
          title: "Extended State",
          body:
            "In the calcium-bound but target-free state, calmodulin adopts an extended dumbbell shape. This exposes hydrophobic patches on each lobe that will grip the target.",
          action: { pdbIndex: 0, style: "cartoon", colorScheme: "byDomain" },
          duration: 7000,
        },
        {
          title: "Wrapped Around Its Target",
          body:
            "When calmodulin binds a target peptide, the central helix unwinds and the protein collapses into a compact globular shape, wrapping around the substrate.",
          action: { pdbIndex: 1, style: "cartoon", colorScheme: "byDomain" },
          duration: 8000,
        },
      ],
    },
  ],

  insulin: [
    {
      title: "Insulin: From Hexamer to Hormone",
      steps: [
        {
          title: "Three Dimers, Six Chains",
          body:
            "Insulin is stored as a hexamer — three dimers arranged with 3-fold symmetry around two central zinc ions. Each dimer (cyan/magenta) is the active hormone unit.",
          action: { style: "cartoon", colorScheme: "chain", spin: true },
          duration: 8000,
        },
        {
          title: "The A and B Chains",
          body:
            "Each monomer has a 21-residue A chain and a 30-residue B chain, linked by two inter-chain and one intra-chain disulfide bond (yellow).",
          action: { style: "cartoon", colorScheme: "chain", showSticks: { atom: "SG" } },
          duration: 8000,
        },
        {
          title: "Zinc Coordination",
          body:
            "Three His-B10 residues (one per dimer) coordinate two zinc ions at the central axis. This zinc coordination is what drives hexamer assembly in storage granules.",
          action: { style: "cartoon", colorScheme: "chain", showSpheres: { resn: "ZN" }, showSticks: { resi: 10 } },
          duration: 8000,
        },
      ],
    },
  ],

  "t4-lysozyme": [
    {
      title: "T4 Lysozyme: A Folding Sandbox",
      steps: [
        {
          title: "Classic α+β Fold",
          body:
            "T4 lysozyme is the textbook model for protein stability studies. It has a small N-terminal β-sheet and a large C-terminal α-helical domain.",
          action: { style: "cartoon", colorScheme: "ss" },
          duration: 7000,
        },
        {
          title: "Core vs. Surface",
          body:
            "Buried residues form a tightly packed hydrophobic core — every mutation to a core residue has measurable effects on stability, often 1–5 kcal/mol per methylene group.",
          action: { style: "surface", colorScheme: "hydrophobicity", surface: { opacity: 0.7 } },
          duration: 8000,
        },
        {
          title: "Cavity-Forming Mutations",
          body:
            "Trp → Ala or Ile → Val mutations create cavities. The classic Leu99→Ala mutation creates a 150 Å³ cavity that has been probed by 100+ ligand-binding experiments.",
          action: { style: "cartoon", colorScheme: "ss", highlight: { resi: [99] } },
          duration: 8000,
        },
      ],
    },
  ],
};

/**
 * Conformational-change trajectories: pairs of states that can be played back.
 * Each trajectory has a slow 'morph' view that fades between the two structures.
 */
const TRAJECTORIES = {
  hemoglobin: {
    title: "Hemoglobin T → R Allosteric Transition",
    states: ["T-state (deoxy)", "R-state (oxy)"],
    description:
      "Watch the textbook allosteric transition: as O₂ binds, the iron moves into the heme plane and the entire tetramer rearranges from T (tense) to R (relaxed).",
    pdbIndexA: 0,
    pdbIndexB: 1,
  },
  "adenylate-kinase": {
    title: "Adenylate Kinase Domain Closure",
    states: ["Open", "Closed"],
    description:
      "One of the largest known domain motions in a single subunit. The NMP and LID domains close ~40–70° around the substrate to enable phosphoryl transfer.",
    pdbIndexA: 0,
    pdbIndexB: 1,
  },
  calmodulin: {
    title: "Calmodulin: Extended → Compact",
    states: ["Extended", "Compact"],
    description:
      "The long central helix partially unwinds and the protein collapses into a globular shape that grips its target peptide.",
    pdbIndexA: 0,
    pdbIndexB: 1,
  },
  "hiv-protease": {
    title: "HIV-1 Protease: Apo vs Inhibitor-Bound",
    states: ["Apo (open)", "Drug-bound (closed)"],
    description:
      "The flap regions close down over the inhibitor. Drug-resistant mutations reshape these flaps to weaken drug binding.",
    pdbIndexA: 1,
    pdbIndexB: 0,
  },
};

if (typeof module !== "undefined") {
  module.exports = { PROTEINS, TOURS, TRAJECTORIES };
}