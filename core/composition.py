"""Resolve a YAML composition spec into per-component molecule counts.

Three modes are supported:

``counts``
    Explicit per-component integer counts. Used by the smoke configs and by
    anyone who already knows exactly how many of each they want.

``molality``
    Salt molality (mol salt / kg solvent), per-solvent mass fractions, and a
    total-solvent-molecule budget. The salt count is computed so that the
    requested molality is satisfied against the resulting total solvent mass.

``molarity``
    Salt molarity (mol salt / L solution), per-solvent mass fractions, and a
    salt-molecule budget. Uses the target solution density (taken from
    ``model.density_g_cm3``) to convert the requested molarity into a volume,
    then derives solvent counts so that the total mass fills that volume at
    the requested density. After NPT relaxes the box to the target density
    the final concentration matches the requested molarity to within
    integer-rounding error on the molecule counts.

The resolver does *not* know molecular weights ahead of time — they have to
be supplied by the caller after the single-molecule structure step has run,
because some components are SMILES (MW from RDKit) and some are ion presets
(MW from atom-mass sum). The molarity mode additionally needs the target
solution density, passed in as a keyword argument.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from utils.constants import AVOGADRO


@dataclass
class ResolvedComponent:
    """Per-component result of composition resolution."""

    name: str
    kind: str         # "solvent" | "cation" | "anion"
    smiles: str | None
    preset: str | None
    charge_source: str
    count: int


def resolve_composition(
    composition_cfg: dict[str, Any],
    components_cfg: list[dict[str, Any]],
    mw_lookup: dict[str, float],
    *,
    density_g_cm3: float | None = None,
) -> list[ResolvedComponent]:
    """Compute the integer molecule count for each component.

    Parameters
    ----------
    composition_cfg:
        The ``composition:`` section of the YAML config.
    components_cfg:
        The ``components:`` section (list of dicts) of the YAML config.
    mw_lookup:
        Map ``component_name → molecular weight (g/mol)``. Must contain an
        entry for every name in ``components_cfg``.
    density_g_cm3:
        Target solution density (g/cm^3); only the ``molarity`` mode reads
        this, but it should be passed in by the caller whenever a config
        might use that mode. Typically taken from ``model.density_g_cm3``.

    Returns
    -------
    list of ResolvedComponent, in the same order as ``components_cfg``.
    """
    mode = str(composition_cfg.get("mode", "counts")).lower()

    if mode == "counts":
        return _resolve_counts(composition_cfg, components_cfg)
    if mode == "molality":
        return _resolve_molality(composition_cfg, components_cfg, mw_lookup)
    if mode == "molarity":
        if density_g_cm3 is None or density_g_cm3 <= 0:
            raise ValueError(
                "composition.mode='molarity' requires a positive target "
                "solution density. Pass density_g_cm3=... (typically from "
                "model.density_g_cm3)."
            )
        return _resolve_molarity(
            composition_cfg, components_cfg, mw_lookup,
            density_g_cm3=density_g_cm3,
        )

    raise ValueError(f"Unsupported composition.mode: {mode!r}")


def _resolve_counts(
    composition_cfg: dict[str, Any],
    components_cfg: list[dict[str, Any]],
) -> list[ResolvedComponent]:
    counts = composition_cfg.get("counts") or {}
    resolved: list[ResolvedComponent] = []
    for comp in components_cfg:
        name = comp["name"]
        if name not in counts:
            raise ValueError(
                f"composition.counts is missing an entry for component {name!r}"
            )
        resolved.append(_make_resolved(comp, int(counts[name])))
    return resolved


def _resolve_molality(
    composition_cfg: dict[str, Any],
    components_cfg: list[dict[str, Any]],
    mw_lookup: dict[str, float],
) -> list[ResolvedComponent]:
    salt = composition_cfg.get("salt") or {}
    cation_name = salt.get("cation")
    anion_name = salt.get("anion")
    molality = float(salt.get("molality", 0.0))
    if not cation_name or not anion_name:
        raise ValueError("composition.salt must specify cation and anion names")
    if molality <= 0:
        raise ValueError("composition.salt.molality must be > 0")

    solvents_cfg = composition_cfg.get("solvents") or {}
    if not solvents_cfg:
        raise ValueError("composition.solvents is empty")

    total_solvent_target = int(composition_cfg.get("total_solvent_molecules", 0))
    if total_solvent_target <= 0:
        raise ValueError("composition.total_solvent_molecules must be > 0")

    mass_fracs: dict[str, float] = {}
    for sname, sspec in solvents_cfg.items():
        if sname not in mw_lookup:
            raise ValueError(
                f"solvent {sname!r} listed in composition.solvents but not in "
                "components: — add it to components: first"
            )
        mf = sspec.get("mass_fraction")
        if mf is None:
            raise ValueError(f"solvent {sname!r}: missing 'mass_fraction'")
        mass_fracs[sname] = float(mf)

    total_mf = sum(mass_fracs.values())
    if total_mf <= 0:
        raise ValueError("composition.solvents mass_fractions sum to zero")
    if abs(total_mf - 1.0) > 1e-6:
        # Silently renormalise so users can pass non-normalised ratios.
        mass_fracs = {k: v / total_mf for k, v in mass_fracs.items()}

    # N_s ∝ (w_s / MW_s); normalise so they sum to total_solvent_target.
    weighted = {s: mass_fracs[s] / mw_lookup[s] for s in mass_fracs}
    denom = sum(weighted.values())
    solvent_counts: dict[str, int] = {
        s: int(round(total_solvent_target * weighted[s] / denom)) for s in weighted
    }
    # Push rounding drift onto the most-abundant solvent so mass fractions
    # stay closest to the requested values.
    drift = total_solvent_target - sum(solvent_counts.values())
    if drift != 0:
        biggest = max(solvent_counts, key=lambda s: solvent_counts[s])
        solvent_counts[biggest] += drift

    # Solvent mass in the box (g for the actual number of molecules present).
    total_solvent_mass_g = sum(
        solvent_counts[s] * mw_lookup[s] for s in solvent_counts
    )
    # N_salt = m × (m_solvent_g / 1000)
    n_salt = int(round(molality * total_solvent_mass_g / 1000.0))
    if n_salt <= 0:
        raise ValueError(
            "Computed salt molecule count is zero. Increase "
            "total_solvent_molecules or molality "
            f"(m={molality}, solvent_mass_g={total_solvent_mass_g:.3f})."
        )

    resolved: list[ResolvedComponent] = []
    for comp in components_cfg:
        name = comp["name"]
        if name == cation_name or name == anion_name:
            count = n_salt
        elif name in solvent_counts:
            count = solvent_counts[name]
        else:
            raise ValueError(
                f"Component {name!r} appears in components: but is referenced "
                "neither in composition.salt.{cation,anion} nor in "
                "composition.solvents.*"
            )
        resolved.append(_make_resolved(comp, count))
    return resolved


def _resolve_molarity(
    composition_cfg: dict[str, Any],
    components_cfg: list[dict[str, Any]],
    mw_lookup: dict[str, float],
    *,
    density_g_cm3: float,
) -> list[ResolvedComponent]:
    """Convert (molarity, target_salt_molecules, density) → integer counts.

    Algorithm for a 1:1 salt:
        V_cm3       = N_salt / (molarity * N_A) * 1000   # solution volume
        m_solution  = density * V_cm3                    # total mass
        m_salt      = N_salt * (MW_cation + MW_anion) / N_A
        m_solvent   = m_solution - m_salt
        N_solvent_s = round(m_solvent * mass_fraction[s] * N_A / MW[s])

    Once Packmol sizes the box against the same density, the resulting
    concentration recovers the requested molarity to within rounding.
    """
    salt = composition_cfg.get("salt") or {}
    cation_name = salt.get("cation")
    anion_name = salt.get("anion")
    molarity = float(salt.get("molarity", 0.0))
    if not cation_name or not anion_name:
        raise ValueError("composition.salt must specify cation and anion names")
    if molarity <= 0:
        raise ValueError("composition.salt.molarity must be > 0")
    if cation_name not in mw_lookup or anion_name not in mw_lookup:
        raise ValueError(
            f"salt cation/anion {cation_name!r}/{anion_name!r} not found in "
            "components — declare them in components: first"
        )

    n_salt = int(composition_cfg.get("target_salt_molecules", 0))
    if n_salt <= 0:
        raise ValueError("composition.target_salt_molecules must be > 0")

    solvents_cfg = composition_cfg.get("solvents") or {}
    if not solvents_cfg:
        raise ValueError("composition.solvents is empty")

    mass_fracs: dict[str, float] = {}
    for sname, sspec in solvents_cfg.items():
        if sname not in mw_lookup:
            raise ValueError(
                f"solvent {sname!r} listed in composition.solvents but not in "
                "components: — add it to components: first"
            )
        mf = sspec.get("mass_fraction")
        if mf is None:
            raise ValueError(f"solvent {sname!r}: missing 'mass_fraction'")
        mass_fracs[sname] = float(mf)

    total_mf = sum(mass_fracs.values())
    if total_mf <= 0:
        raise ValueError("composition.solvents mass_fractions sum to zero")
    if abs(total_mf - 1.0) > 1e-6:
        mass_fracs = {k: v / total_mf for k, v in mass_fracs.items()}

    # Solution volume in cm^3 from molarity.
    v_cm3 = n_salt / (molarity * AVOGADRO) * 1000.0
    m_solution_g = density_g_cm3 * v_cm3

    mw_salt = mw_lookup[cation_name] + mw_lookup[anion_name]
    m_salt_g = n_salt * mw_salt / AVOGADRO
    m_solvent_g = m_solution_g - m_salt_g

    if m_solvent_g <= 0:
        raise ValueError(
            f"Computed solvent mass is non-positive (got {m_solvent_g:.3e} g). "
            f"Solution mass {m_solution_g:.3e} g is smaller than salt mass "
            f"{m_salt_g:.3e} g — molarity {molarity} is too high relative to "
            f"density {density_g_cm3} g/cm^3."
        )

    solvent_counts: dict[str, int] = {}
    for s, mf in mass_fracs.items():
        n_raw = m_solvent_g * mf * AVOGADRO / mw_lookup[s]
        solvent_counts[s] = max(int(round(n_raw)), 1)

    resolved: list[ResolvedComponent] = []
    for comp in components_cfg:
        name = comp["name"]
        if name == cation_name or name == anion_name:
            count = n_salt
        elif name in solvent_counts:
            count = solvent_counts[name]
        else:
            raise ValueError(
                f"Component {name!r} appears in components: but is referenced "
                "neither in composition.salt.{cation,anion} nor in "
                "composition.solvents.*"
            )
        resolved.append(_make_resolved(comp, count))
    return resolved


def _make_resolved(comp: dict[str, Any], count: int) -> ResolvedComponent:
    return ResolvedComponent(
        name=comp["name"],
        kind=str(comp.get("kind", "solvent")),
        smiles=comp.get("smiles"),
        preset=comp.get("preset"),
        charge_source=str(comp.get("charge_source", "openff_am1bcc")),
        count=count,
    )


__all__ = ["ResolvedComponent", "resolve_composition"]
