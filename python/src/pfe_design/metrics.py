"""Geometry evaluation metrics for generated complexes.

Four families of measurement, all computed on plain ``torch`` tensors so that
they can be used inside a loss, inside a validation loop, or standalone:

* **iRMSD** - interface RMSD after superposing the complex on the receptor.
* **Ligand RMSD** - the displacement of the ligand after the same
  superposition, which is the number a docking or binding study would quote.
* **BSA** - buried surface area on complex formation, by the Shrake-Rupley
  spherical-point method.
* **Steric clash penalty** - the overlap of atomic spheres that is too severe
  to be physical.

Conventions
-----------
The superposition helpers return the aligned coordinates *and* the rotation
and translation, so a ligand can be carried along with the receptor fit
instead of being fitted separately (which would hide the receptor motion).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import torch
from torch import Tensor

from .config import MetricConfig

__all__ = [
    "ClashReport",
    "BSAResult",
    "ComplexEvaluation",
    "VDW_RADII",
    "vdw_radii",
    "superpose",
    "apply_transform",
    "interface_rmsd",
    "ligand_rmsd",
    "sasa",
    "buried_surface_area",
    "clash_report",
    "clash_penalty",
]

# Van der Waals radii in Angstrom. Every element that Biopython's
# Bio.PDB.SASA.ATOMIC_RADII also defines has the same value here, so the two
# implementations can be compared directly (see the test-suite).
VDW_RADII = {
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "P": 1.80,
    "H": 1.20,
    "SE": 1.90,
    "F": 1.47,
    "CL": 1.75,
    "BR": 1.85,
    "I": 1.98,
    "MG": 1.73,
    "ZN": 1.39,
    "CA": 2.31,
    "NA": 2.27,
    "K": 2.75,
    "MN": 1.61,
    "CO": 1.52,
    "NI": 1.63,
    "CU": 1.40,
}
_DEFAULT_RADIUS = 1.70


@dataclass(frozen=True)
class ClashReport:
    """Detailed result of a steric clash analysis of a single structure."""

    n_clashes: int
    """Number of clashing atom pairs."""

    penalty: float
    """Sum of the squared overlaps of the clashing pairs."""

    max_overlap: float
    """Largest single overlap (Angstrom)."""

    pairs: Tensor
    """``(k, 2)`` indices of the clashing pairs."""

    n_atoms: int
    """Number of atoms that were analysed."""


@dataclass(frozen=True)
class BSAResult:
    """Buried surface area on complex formation."""

    bsa: float
    """Buried surface area, ``SASA(free A) + SASA(free B) - SASA(complex)``."""

    sasa_free_a: float
    sasa_free_b: float
    sasa_complex: float
    per_atom_complex: Tensor
    """Per-atom solvent-accessible surface area in the complex."""


def vdw_radii(elements: Sequence[str]) -> Tensor:
    """Look up van der Waals radii for a list of element symbols."""
    return torch.tensor(
        [VDW_RADII.get(str(e).strip().upper(), _DEFAULT_RADIUS) for e in elements],
        dtype=torch.float64,
    )


# --------------------------------------------------------------------- fitting


def superpose(mobile: Tensor, target: Tensor) -> tuple[Tensor, Tensor, Tensor]:
    """Superpose ``mobile`` onto ``target`` (Kabsch).

    Parameters
    ----------
    mobile, target:
        Coordinates of identical shape ``(N, 3)`` or ``(B, N, 3)``.

    Returns
    -------
    (aligned, rotation, translation)
        ``aligned`` has the shape of ``mobile``; ``rotation`` is ``(3, 3)`` or
        ``(B, 3, 3)``; ``translation`` is ``(3,)`` or ``(B, 3)``.
    """
    if mobile.shape != target.shape:
        raise ValueError(
            f"mobile and target must have the same shape; got {tuple(mobile.shape)} "
            f"and {tuple(target.shape)}"
        )
    if mobile.dim() == 2:
        m = mobile.unsqueeze(0)
        t = target.unsqueeze(0)
        squeeze = True
    elif mobile.dim() == 3:
        m, t, squeeze = mobile, target, False
    else:
        raise ValueError(f"expected a 2-D or 3-D coordinate tensor, got {mobile.dim()} dimensions")

    if m.shape[1] < 3:
        raise ValueError("superposition needs at least three points")

    work = torch.float64
    md = m.to(work)
    td = t.to(work)
    mc = md.mean(dim=1, keepdim=True)
    tc = td.mean(dim=1, keepdim=True)
    p = md - mc
    q = td - tc

    h = p.transpose(1, 2) @ q
    u, _s, vh = torch.linalg.svd(h)
    det = torch.det(vh.transpose(1, 2) @ u.transpose(1, 2))
    ones = torch.ones_like(det)
    d = torch.diag_embed(torch.stack([ones, ones, torch.sign(det)], dim=-1))
    rotation = vh.transpose(1, 2) @ d @ u.transpose(1, 2)
    translation = tc[:, 0, :] - torch.bmm(rotation, mc.transpose(1, 2))[:, :, 0]

    aligned = torch.bmm(md, rotation.transpose(1, 2)) + translation.unsqueeze(1)
    aligned = aligned.to(m.dtype)
    rotation = rotation.to(m.dtype)
    translation = translation.to(m.dtype)

    if squeeze:
        return aligned[0], rotation[0], translation[0]
    return aligned, rotation, translation


def apply_transform(coords: Tensor, rotation: Tensor, translation: Tensor) -> Tensor:
    """Apply a rotation and translation to coordinates."""
    if coords.dim() == 2:
        return coords @ rotation.transpose(-1, -2) + translation
    return torch.bmm(coords, rotation.transpose(-1, -2)) + translation.unsqueeze(-2)


@dataclass(frozen=True)
class ComplexEvaluation:
    """The four geometry measurements of a complex, computed in one call.

    Attributes
    ----------
    irmsd:
        Interface RMSD after superposing on the receptor, in Angstrom.
        ``None`` when there is no interface to measure.
    ligand_rmsd:
        Ligand RMSD after the same superposition, in Angstrom.  ``None`` when
        no ligand was given, or when the ligand mask is empty.
    bsa:
        Buried surface area on complex formation, in square Angstrom.
        ``None`` when the two partner groups were not given.
    n_clashes:
        Number of clashing atom pairs.
    clash_penalty:
        Sum of the squared overlaps of the clashing pairs.
    max_overlap:
        Largest single overlap, in Angstrom.
    """

    irmsd: float | None
    ligand_rmsd: float | None
    bsa: float | None
    n_clashes: int
    clash_penalty: float
    max_overlap: float

    @property
    def is_clean(self) -> bool:
        """True when no steric clash was found."""
        return self.n_clashes == 0


def evaluate_complex(
    mobile: Tensor,
    target: Tensor,
    interface_mask: Tensor,
    receptor_mask: Tensor | None = None,
    ligand_mask: Tensor | None = None,
    group_a: Tensor | None = None,
    group_b: Tensor | None = None,
    radii: Tensor | None = None,
    residue_index: Tensor | None = None,
    cfg: MetricConfig | None = None,
) -> ComplexEvaluation:
    """Measure a complex: iRMSD, ligand RMSD, buried area and clashes at once.

    This is the façade over the individual metrics, for the case where a
    generated complex has to be judged as a whole.  Each measurement that
    needs an extra piece of information (a ligand, the two partner groups)
    is reported as ``None`` when that information was not supplied, rather
    than being guessed at.

    Parameters
    ----------
    mobile, target:
        Coordinates of the complex and of the reference, ``(N, 3)``.
    interface_mask:
        1 for the residues that form the interface, ``(N,)``.
    receptor_mask:
        The atoms used for the superposition.  Without it the fit uses
        everything that is not the interface, which for a complex with a
        ligand includes the ligand.
    ligand_mask:
        1 for the ligand's atoms; needed for the ligand RMSD.
    group_a, group_b:
        The two partner groups; needed for the buried surface area.
    radii:
        Atomic radii, needed for the clash analysis in ``vdw`` mode.
    residue_index:
        Residue index per atom; same-residue pairs are then not counted as
        clashes.
    cfg:
        Metric settings.

    Returns
    -------
    ComplexEvaluation
    """
    if mobile.dim() != 2:
        raise ValueError(
            f"evaluate_complex takes one complex, with shape (n_residues, 3); "
            f"got {tuple(mobile.shape)}. Use interface_rmsd for a batch."
        )
    cfg = cfg or MetricConfig()

    irmsd_value = interface_rmsd(mobile, target, interface_mask, receptor_mask=receptor_mask)
    irmsd = float(irmsd_value) if bool(torch.isfinite(irmsd_value).all()) else None

    ligand_value: float | None = None
    if ligand_mask is not None:
        fit = receptor_mask
        if fit is None:
            # Without an explicit receptor, the fit uses everything that is
            # neither the interface nor the ligand - the ligand must not be
            # part of the fit it is measured against.
            fit = ~(
                _as_batch(interface_mask, 1, mobile.shape[0])[0]
                | _as_batch(ligand_mask, 1, mobile.shape[0])[0]
            )
        measured = ligand_rmsd(mobile, target, fit, ligand_mask)
        ligand_value = float(measured) if bool(torch.isfinite(measured).all()) else None

    bsa_value: float | None = None
    if group_a is not None and group_b is not None:
        if radii is None:
            radii = torch.full((mobile.shape[0],), 1.70, dtype=torch.float64)
        bsa_value = buried_surface_area(mobile, radii, group_a, group_b, cfg).bsa

    report = clash_report(mobile, cfg, radii, residue_index)
    return ComplexEvaluation(
        irmsd=irmsd,
        ligand_rmsd=ligand_value,
        bsa=bsa_value,
        n_clashes=report.n_clashes,
        clash_penalty=report.penalty,
        max_overlap=report.max_overlap,
    )


def _as_batch(mask: Tensor, batch: int, n: int) -> Tensor:
    if mask.dim() == 1:
        if mask.shape[0] != n:
            raise ValueError(
                f"mask must have shape ({batch}, {n}) or ({n},); got {tuple(mask.shape)}"
            )
        mask = mask.unsqueeze(0).expand(batch, n)
    if mask.shape != (batch, n):
        raise ValueError(
            f"mask must have shape ({batch}, {n}) or ({n},); got {tuple(mask.shape)}"
        )
    return mask.to(torch.bool)


def _fit_subset(
    mobile: Tensor, target: Tensor, selection: Tensor
) -> tuple[Tensor, Tensor, Tensor]:
    """Superpose the selected atoms and carry the whole structure with them."""
    if int(selection.sum()) < 3:
        identity = torch.eye(3, dtype=mobile.dtype, device=mobile.device)
        zero = torch.zeros(3, dtype=mobile.dtype, device=mobile.device)
        return mobile.clone(), identity, zero
    aligned, rotation, translation = superpose(mobile[selection], target[selection])
    moved = apply_transform(mobile, rotation, translation)
    return moved, rotation, translation


def interface_rmsd(
    mobile: Tensor,
    target: Tensor,
    interface_mask: Tensor,
    valid_mask: Tensor | None = None,
    superposition: str = "receptor",
    receptor_mask: Tensor | None = None,
) -> Tensor:
    """Interface RMSD after superposing on the receptor (or the interface).

    Parameters
    ----------
    mobile, target:
        Coordinates of the complex, ``(N, 3)`` or ``(B, N, 3)``.
    interface_mask:
        1 for the residues that form the interface, ``(N,)`` or ``(B, N)``.
    valid_mask:
        Optional mask of modelled residues; residues outside it are ignored.
    superposition:
        ``"receptor"`` fits the non-interface residues, ``"interface"`` fits
        the interface itself.
    receptor_mask:
        Optional explicit set of atoms for the fit, ``(N,)`` or ``(B, N)``.
        Without it the fit uses everything that is not interface, which for a
        complex that contains a ligand silently includes the ligand's own
        atoms - pass the receptor explicitly in that case.

    Returns
    -------
    Tensor
        Per-structure iRMSD, ``(B,)`` (or a 0-dim tensor for unbatched input).
    """
    if mobile.shape != target.shape:
        raise ValueError("mobile and target must have the same shape")
    batched = mobile.dim() == 3
    m = mobile.unsqueeze(0) if not batched else mobile
    t = target.unsqueeze(0) if not batched else target
    batch, n, _ = m.shape

    iface = _as_batch(interface_mask, batch, n)
    valid = _as_batch(valid_mask, batch, n) if valid_mask is not None else torch.ones(
        (batch, n), dtype=torch.bool, device=m.device
    )
    iface = iface & valid

    if superposition == "receptor":
        if receptor_mask is None:
            selection = valid & ~iface
        else:
            selection = _as_batch(receptor_mask, batch, n) & valid & ~iface
    elif superposition == "interface":
        selection = iface
    else:
        raise ValueError(
            f"superposition must be 'receptor' or 'interface'; got {superposition!r}"
        )

    values = []
    for b in range(batch):
        if int(iface[b].sum()) == 0:
            # There is no interface, so there is nothing to measure. A zero
            # here would be a perfect score for a measurement that was never
            # made, so it is reported as not-a-number.
            values.append(torch.full((), float("nan"), dtype=m.dtype, device=m.device))
            continue
        moved, _r, _t = _fit_subset(m[b], t[b], selection[b])
        delta = moved[iface[b]] - t[b][iface[b]]
        rmsd = torch.sqrt((delta ** 2).sum(dim=-1).mean())
        values.append(rmsd)
    out = torch.stack(values)
    return out if batched else out[0]


def ligand_rmsd(
    mobile: Tensor,
    target: Tensor,
    receptor_mask: Tensor,
    ligand_mask: Tensor,
    valid_mask: Tensor | None = None,
) -> Tensor:
    """Ligand RMSD after superposing the receptor.

    The ligand is carried along with the receptor fit; it is not fitted
    itself, so a ligand that does not follow the receptor is penalised.
    """
    if mobile.shape != target.shape:
        raise ValueError("mobile and target must have the same shape")
    batched = mobile.dim() == 3
    m = mobile.unsqueeze(0) if not batched else mobile
    t = target.unsqueeze(0) if not batched else target
    batch, n, _ = m.shape

    receptor = _as_batch(receptor_mask, batch, n)
    ligand = _as_batch(ligand_mask, batch, n)
    valid = _as_batch(valid_mask, batch, n) if valid_mask is not None else torch.ones(
        (batch, n), dtype=torch.bool, device=m.device
    )
    receptor = receptor & valid
    ligand = ligand & valid

    values = []
    for b in range(batch):
        if int(ligand[b].sum()) == 0:
            # No ligand, so no ligand RMSD - not a perfect zero.
            values.append(torch.full((), float("nan"), dtype=m.dtype, device=m.device))
            continue
        moved, _r, _t = _fit_subset(m[b], t[b], receptor[b])
        delta = moved[ligand[b]] - t[b][ligand[b]]
        values.append(torch.sqrt((delta ** 2).sum(dim=-1).mean()))
    out = torch.stack(values)
    return out if batched else out[0]


# ----------------------------------------------------------------------- BSA


def _sphere_points(n_points: int, device: torch.device) -> Tensor:
    """Fibonacci-spiral points spread evenly over the unit sphere."""
    idx = torch.arange(n_points, dtype=torch.float64, device=device)
    phi = torch.acos(1.0 - 2.0 * (idx + 0.5) / n_points)
    theta = math.pi * (1.0 + 5.0 ** 0.5) * idx
    return torch.stack(
        [torch.cos(theta) * torch.sin(phi), torch.sin(theta) * torch.sin(phi), torch.cos(phi)],
        dim=-1,
    )


def sasa(
    coords: Tensor,
    radii: Tensor,
    probe: float = 1.4,
    n_points: int = 100,
    chunk: int = 32,
) -> Tensor:
    """Solvent-accessible surface area per atom (Shrake-Rupley).

    The cost is O(n_atoms^2 * n_points): every sphere point is tested against
    every atom. That is fine for validating a finished complex, and too slow
    to put inside a per-step training loop for a large one - for that, use a
    neighbour list.

    Parameters
    ----------
    coords:
        Atomic coordinates, ``(N, 3)``.
    radii:
        Atomic radii, ``(N,)``.
    probe:
        Radius of the solvent probe, ``1.4`` Angstrom for water.
    n_points:
        Number of sphere points per atom; more points converge faster.
    chunk:
        Number of atoms processed at once (memory/speed trade-off).

    Returns
    -------
    Tensor
        Per-atom SASA in square Angstrom.
    """
    coords = coords.to(torch.float64)
    radii = radii.to(torch.float64)
    n = coords.shape[0]
    if radii.shape != (n,):
        raise ValueError(f"radii must have shape ({n},); got {tuple(radii.shape)}")
    if n == 0:
        return torch.zeros(0, dtype=torch.float64)

    unit = _sphere_points(n_points, coords.device)
    expanded = radii + probe
    area_per_point = 4.0 * math.pi * expanded ** 2 / n_points

    out = torch.zeros(n, dtype=torch.float64, device=coords.device)
    atom_ids = torch.arange(n, device=coords.device)
    for start in range(0, n, chunk):
        end = min(start + chunk, n)
        points = (
            coords[start:end, None, :] + expanded[start:end, None, None] * unit[None, :, :]
        )  # (c, P, 3)
        distances = torch.cdist(points.reshape(-1, 3), coords).reshape(
            end - start, n_points, n
        )
        # An atom must not occlude itself. The chunk index c stands for the
        # atom (start + c), while the last axis of `distances` runs over all
        # n atoms - so the self-hit mask is built on that axis, not on the
        # chunk axis. (Getting this wrong only shows up once there is more
        # than one chunk, which is why the system must be tested at a size
        # larger than the chunk.)
        chunk_index = torch.arange(end - start, device=coords.device)
        self_hit = (chunk_index[:, None, None] + start) == atom_ids[None, None, :]
        distances = distances.masked_fill(self_hit, float("inf"))
        occluded = (distances < expanded[None, None, :]).any(dim=-1)
        out[start:end] = (~occluded).sum(dim=-1) * area_per_point[start:end]
    return out


def _biopython_structure(coords: Tensor, elements: Sequence[str]):
    """Build a one-atom-per-residue Bio.PDB structure for the reference SASA.

    Biopython resolves an atom's radius from its element, so the atoms are
    given the elements they really are.  The radii table used here matches
    Biopython's own for every element the two share.
    """
    try:
        from Bio.PDB.Atom import Atom
        from Bio.PDB.Chain import Chain
        from Bio.PDB.Model import Model
        from Bio.PDB.Residue import Residue
        from Bio.PDB.Structure import Structure
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise RuntimeError(
            "Biopython is not installed; install it with `pip install biopython` "
            "or use the default PyTorch implementation"
        ) from exc

    structure = Structure("complex")
    model = Model(0)
    structure.add(model)
    chain = Chain("A")
    model.add(chain)
    for index, (xyz, element) in enumerate(zip(coords, elements, strict=False)):
        name = str(element).strip().upper() or "C"
        residue = Residue((" ", index + 1, " "), name, "")
        atom = Atom(
            name=name,
            coord=[float(v) for v in xyz],
            bfactor=0.0,
            occupancy=1.0,
            altloc=" ",
            fullname=name,
            serial_number=index + 1,
            element=name,
        )
        residue.add(atom)
        chain.add(residue)
    return structure


def sasa_biopython(
    coords: Tensor,
    elements: Sequence[str],
    probe: float = 1.4,
    n_points: int = 100,
) -> Tensor:
    """Per-atom solvent-accessible surface area, from Biopython.

    This is the reference implementation, kept beside the vectorised PyTorch
    one so the two can be compared and so that a user who wants the reference
    number rather than the fast one can have it.  Biopython is an optional
    dependency; without it this raises a clear error rather than failing
    obscurely.

    Parameters
    ----------
    coords:
        Atomic coordinates, ``(N, 3)``.
    elements:
        The element symbol of each atom, ``(N,)``.
    probe:
        Radius of the solvent probe, ``1.4`` Angstrom for water.
    n_points:
        Number of sphere points per atom.

    Returns
    -------
    Tensor
        Per-atom SASA in square Angstrom.
    """
    try:
        from Bio.PDB.SASA import ShrakeRupley
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise RuntimeError(
            "Biopython is not installed; install it with `pip install biopython` "
            "or use the default PyTorch implementation"
        ) from exc

    if coords.shape[0] == 0:
        # Biopython refuses an empty structure; the PyTorch path returns an
        # empty tensor, and the two should agree.
        return torch.zeros(0, dtype=torch.float64)

    structure = _biopython_structure(coords, elements)
    calculator = ShrakeRupley(probe_radius=probe, n_points=n_points)
    calculator.compute(structure, level="R")
    values = [
        float(atom.sasa) for residue in structure.get_residues() for atom in residue
    ]
    return torch.tensor(values, dtype=torch.float64)


def buried_surface_area(
    complex_coords: Tensor,
    complex_radii: Tensor,
    group_a: Tensor,
    group_b: Tensor,
    cfg: MetricConfig | None = None,
    implementation: str = "torch",
    elements: Sequence[str] | None = None,
) -> BSAResult:
    """Buried surface area on complex formation.

    Parameters
    ----------
    complex_coords:
        Coordinates of every atom in the complex, ``(N, 3)``.
    complex_radii:
        Radii of every atom, ``(N,)``.
    group_a, group_b:
        Indices of the atoms belonging to each partner.
    cfg:
        Metric settings (probe radius, number of sphere points).
    implementation:
        ``"torch"`` (the default) uses the vectorised PyTorch calculation;
        ``"biopython"`` uses Biopython's reference implementation and then
        needs ``elements``.
    elements:
        The element symbol of each atom; required for the Biopython
        implementation, which resolves radii from the element.

    Returns
    -------
    BSAResult
    """
    cfg = cfg or MetricConfig()
    complex_coords = complex_coords.to(torch.float64)
    complex_radii = complex_radii.to(torch.float64)
    idx_a = torch.as_tensor(group_a, dtype=torch.long).reshape(-1)
    idx_b = torch.as_tensor(group_b, dtype=torch.long).reshape(-1)
    if bool(torch.isin(idx_a, idx_b).any()):
        raise ValueError(
            "the two partner groups share atoms; a group cannot bury itself. "
            "Give the burial of A against B, with the partners disjoint."
        )

    # The complex term is the surface of the two partners together.  Atoms
    # outside both groups (ions, a third chain) are not part of the burial of
    # A against B, and including them would make the difference meaningless -
    # in the limit, a large third group makes the "buried area" negative,
    # which is not a physical quantity.
    union = torch.unique(torch.cat([idx_a, idx_b]))
    if implementation == "biopython":
        if elements is None:
            raise ValueError(
                "the Biopython implementation resolves radii from the element, "
                "so the elements of every atom must be given"
            )
        element_list = [str(e) for e in elements]

        def area_for(indices: Tensor) -> Tensor:
            return sasa_biopython(
                complex_coords[indices],
                [element_list[i] for i in indices.tolist()],
                probe=cfg.bsa_probe_radius,
                n_points=cfg.bsa_n_points,
            )

        per_atom_union = area_for(union)
        free_a = area_for(idx_a).sum()
        free_b = area_for(idx_b).sum()
    else:
        per_atom_union = sasa(
            complex_coords[union],
            complex_radii[union],
            probe=cfg.bsa_probe_radius,
            n_points=cfg.bsa_n_points,
        )
        free_a = sasa(
            complex_coords[idx_a],
            complex_radii[idx_a],
            probe=cfg.bsa_probe_radius,
            n_points=cfg.bsa_n_points,
        ).sum()
        free_b = sasa(
            complex_coords[idx_b],
            complex_radii[idx_b],
            probe=cfg.bsa_probe_radius,
            n_points=cfg.bsa_n_points,
        ).sum()
    sasa_complex = per_atom_union.sum()
    bsa = free_a + free_b - sasa_complex

    # The per-atom values are reported for the whole complex, with the atoms
    # outside the two partners left at zero.
    per_atom = torch.zeros(complex_coords.shape[0], dtype=torch.float64)
    per_atom[union] = per_atom_union

    return BSAResult(
        bsa=float(bsa),
        sasa_free_a=float(free_a),
        sasa_free_b=float(free_b),
        sasa_complex=float(sasa_complex),
        per_atom_complex=per_atom,
    )


# --------------------------------------------------------------------- clashes


def _overlap_matrix(
    coords: Tensor,
    cfg: MetricConfig,
    radii: Tensor | None = None,
) -> Tensor:
    """Overlap (Angstrom) of every atom pair; 0 where there is no overlap."""
    n = coords.shape[0]
    if n == 0:
        return torch.zeros((0, 0), dtype=torch.float64)
    distances = torch.cdist(coords.to(torch.float64), coords.to(torch.float64))
    if cfg.clash_mode == "vdw":
        if radii is None:
            raise ValueError("vdw clash mode needs atomic radii")
        r = radii.to(torch.float64)
        limit = (r[:, None] + r[None, :]) - float(cfg.vdw_tolerance)
    else:
        limit = torch.full((n, n), float(cfg.clash_distance), dtype=torch.float64)
    # The overlap of two coincident atoms is the full threshold, not zero: the
    # worst possible steric overlap must be the most visible one. The
    # self-pairs are excluded by the upper-triangle mask, not here.
    return limit - distances


def clash_report(
    coords: Tensor,
    cfg: MetricConfig | None = None,
    radii: Tensor | None = None,
    residue_index: Tensor | None = None,
) -> ClashReport:
    """Analyse steric clashes in a single structure.

    Parameters
    ----------
    coords:
        Atomic coordinates, ``(N, 3)``.
    cfg:
        Metric settings (threshold, mode, tolerance).
    radii:
        Atomic radii; required in ``vdw`` clash mode.
    residue_index:
        Optional residue index per atom.  When given, pairs inside the same
        residue are ignored, because a covalent bond is shorter than any
        clash threshold and would otherwise be reported every time.
    """
    cfg = cfg or MetricConfig()
    overlap = _overlap_matrix(coords, cfg, radii)
    n = coords.shape[0]

    tri = torch.triu(torch.ones((n, n), dtype=torch.bool), diagonal=1)
    if residue_index is not None and cfg.exclude_intra_residue_pairs:
        res = torch.as_tensor(residue_index, dtype=torch.long).reshape(-1)
        same = res[:, None] == res[None, :]
        tri = tri & ~same

    hits = tri & (overlap > 0)
    pairs = torch.nonzero(hits, as_tuple=False)
    if pairs.numel() == 0:
        return ClashReport(
            n_clashes=0,
            penalty=0.0,
            max_overlap=0.0,
            pairs=pairs.reshape(0, 2),
            n_atoms=n,
        )
    values = overlap[hits]
    return ClashReport(
        n_clashes=int(pairs.shape[0]),
        penalty=float((values ** 2).sum()),
        max_overlap=float(values.max()),
        pairs=pairs,
        n_atoms=n,
    )


def clash_penalty(
    coords: Tensor,
    cfg: MetricConfig | None = None,
    radii: Tensor | None = None,
    residue_index: Tensor | None = None,
) -> Tensor:
    """Total clash penalty of each structure in a batch, ``(B,)``."""
    cfg = cfg or MetricConfig()
    batched = coords.dim() == 3
    c = coords.unsqueeze(0) if not batched else coords
    batch = c.shape[0]
    values = []
    for b in range(batch):
        res = None
        if residue_index is not None:
            res = residue_index[b] if residue_index.dim() == 2 else residue_index
        atom_radii = None
        if radii is not None:
            # Radii may be given per structure or shared by all of them.
            atom_radii = radii[b] if radii.dim() == 2 else radii
        report = clash_report(c[b], cfg, atom_radii, res)
        values.append(torch.tensor(report.penalty, dtype=torch.float64))
    out = torch.stack(values)
    return out if batched else out[0]


def clash_penalty_fn(coords: Tensor, cfg: MetricConfig | None = None) -> Tensor:
    """Clash penalty of a batch of structures, differentiable in ``coords``.

    Unlike :func:`clash_penalty` this keeps the graph, so it can be used as a
    loss term during refinement.

    It is deliberately narrower than :func:`clash_report`: a fixed distance
    threshold, no residue indices, and every pair counted.  A covalent bond
    is shorter than any sensible threshold, so this function scores a bonded
    pair as a clash - use :func:`clash_report` for anything that is read by a
    person.
    """
    cfg = cfg or MetricConfig()
    distances = torch.cdist(coords.to(torch.float64), coords.to(torch.float64))
    limit = float(cfg.clash_distance)
    overlap = (limit - distances).clamp(min=0.0)
    n = coords.shape[0]
    if n > 1:
        eye = torch.eye(n, dtype=torch.bool, device=coords.device)
        overlap = overlap.masked_fill(eye, 0.0)
    # Count each pair once.
    tri = torch.triu(torch.ones((n, n), dtype=torch.bool), diagonal=1)
    overlap = overlap * tri
    return (overlap ** 2).sum()
