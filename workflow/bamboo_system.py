"""Build BAMBOO-ready systems for the viscosity workflow."""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from core.bamboo_packing import PackmolBuilder, PackmolComponent
from core.bamboo_structure import generate_molecule, load_ion_preset, write_xyz
from core.composition import ResolvedComponent, resolve_composition
from core.data_builder import AtomRecord, build_bamboo_data
from forcefields.openff_charges import get_openff_charges


@dataclass(frozen=True)
class BambooSystemResult:
    """Artifacts produced by the BAMBOO system builder."""

    data_file: Path
    build_report: Path
    box_lengths: tuple[float, float, float]
    elements: list[str]
    n_molecules: int
    n_atoms: int


def build_bamboo_system(
    cfg: dict[str, Any],
    workdir: Path,
    *,
    config_path: Path | None = None,
) -> BambooSystemResult:
    """Generate single-component or multi-component BAMBOO ``in.data``."""
    case_cfg = cfg.get("case") or {}
    model_cfg = cfg.get("model") or {}
    components_cfg = _components_from_config(cfg)
    composition_cfg = _composition_from_config(cfg, components_cfg)

    workdir.mkdir(parents=True, exist_ok=True)

    print("\n-- Stage 1: build BAMBOO system -------------------------------")
    print(f"  case    : {case_cfg.get('name', '(unnamed)')}")
    print(f"  workdir : {workdir}")

    per_component = _build_per_component(components_cfg, workdir)
    mw_lookup = {p["name"]: p["mw"] for p in per_component}
    density_g_cm3 = float(model_cfg.get("density_g_cm3", 1.0))
    resolved = resolve_composition(
        composition_cfg,
        components_cfg,
        mw_lookup,
        density_g_cm3=density_g_cm3,
    )

    print("\n  resolved composition:")
    for item in resolved:
        print(f"    {item.name:<10s} count={item.count}")

    pack_components = [
        PackmolComponent(
            name=p["name"],
            xyz_file=p["xyz_path"],
            count=r.count,
            n_atoms_per_molecule=len(p["xyz_atoms"]),
            molecular_weight=p["mw"],
        )
        for p, r in zip(per_component, resolved)
    ]

    packer = PackmolBuilder(workdir)
    system_xyz, box_lengths = packer.build_box(
        pack_components,
        density_g_cm3=density_g_cm3,
        density_scale=float(model_cfg.get("packmol_density_scale", 0.85)),
        tolerance=float(model_cfg.get("packmol_tolerance_A", 2.0)),
        seed=int(model_cfg.get("packmol_seed", 192911)),
        max_attempts=int(model_cfg.get("packmol_max_attempts", 5)),
        seed_step=int(model_cfg.get("packmol_seed_step", 17)),
        nloop=model_cfg.get("packmol_nloop"),
        full_box=bool(model_cfg.get("packmol_full_box", False)),
        margin=model_cfg.get("packmol_margin_A"),
        allow_imperfect=not bool(model_cfg.get("packmol_strict", False)),
        box_aspect_ratio=(
            tuple(model_cfg["box_aspect_ratio"])
            if "box_aspect_ratio" in model_cfg
            else None
        ),
        box_lengths=(
            tuple(model_cfg["box_lengths_A"])
            if "box_lengths_A" in model_cfg
            else None
        ),
    )

    system_coords = _load_xyz_coords(system_xyz)
    atom_records, n_molecules = _assemble_atom_records(
        per_component,
        resolved,
        system_coords,
    )

    bamboo_cfg = cfg.get("bamboo") or {}
    data_name = str(bamboo_cfg.get("data_file", "in.data"))
    data_path, elements = build_bamboo_data(
        atom_records,
        box_lengths=box_lengths,
        output=workdir / data_name,
        title=f"LAMMPS data file - BAMBOO viscosity system ({case_cfg.get('name', 'case')})",
        neutrality_tolerance=float(bamboo_cfg.get("neutrality_tolerance", 1.0e-3)),
    )

    Lx, Ly, Lz = box_lengths
    report = {
        "case_name": case_cfg.get("name"),
        "workdir": str(workdir),
        "config_path": str(config_path) if config_path is not None else None,
        "n_atoms": len(atom_records),
        "n_molecules": n_molecules,
        "n_components": len(resolved),
        "components": [
            {
                "name": p["name"],
                "kind": r.kind,
                "preset": p["cfg"].get("preset"),
                "smiles": p["cfg"].get("smiles"),
                "count": r.count,
                "n_atoms_per_molecule": len(p["xyz_atoms"]),
                "molecular_weight_g_mol": p["mw"],
                "net_charge_per_molecule": float(sum(p["charges"])),
                "charge_method": p["charge_method"],
            }
            for p, r in zip(per_component, resolved)
        ],
        "elements_in_pair_coeff": elements,
        "box_lengths_A": [Lx, Ly, Lz],
        "box_volume_A3": Lx * Ly * Lz,
        "density_target_g_cm3": density_g_cm3,
        "in_data": str(data_path),
        "system_xyz": str(system_xyz),
    }
    report_path = workdir / str(bamboo_cfg.get("build_report_file", "build_report.json"))
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False))

    print(f"  in.data : {data_path}")
    print(f"  report  : {report_path}")
    return BambooSystemResult(
        data_file=data_path,
        build_report=report_path,
        box_lengths=box_lengths,
        elements=elements,
        n_molecules=n_molecules,
        n_atoms=len(atom_records),
    )


def load_bamboo_build_report(workdir: Path, cfg: dict[str, Any]) -> dict[str, Any]:
    """Load the build report written by :func:`build_bamboo_system`."""
    bamboo_cfg = cfg.get("bamboo") or {}
    report_path = workdir / str(bamboo_cfg.get("build_report_file", "build_report.json"))
    if not report_path.exists():
        raise FileNotFoundError(
            f"BAMBOO build report missing: {report_path}. "
            "Run build_system before write_input."
        )
    payload = json.loads(report_path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"BAMBOO build report must be a JSON object: {report_path}")
    return payload


def _components_from_config(cfg: dict[str, Any]) -> list[dict[str, Any]]:
    components = cfg.get("components")
    if isinstance(components, list) and components:
        return [dict(item) for item in components]

    # Backward-compatible single-molecule shorthand.
    case_cfg = cfg.get("case") or {}
    smiles = case_cfg.get("smiles")
    name = case_cfg.get("name", "molecule")
    if not smiles:
        raise ValueError(
            "Viscosity config must define either components: or case.smiles."
        )
    charge_source = str(case_cfg.get("charge_source", "openff_am1bcc"))
    return [
        {
            "name": name,
            "kind": "solvent",
            "smiles": smiles,
            "charge_source": charge_source,
        }
    ]


def _composition_from_config(
    cfg: dict[str, Any],
    components_cfg: list[dict[str, Any]],
) -> dict[str, Any]:
    composition = cfg.get("composition")
    if isinstance(composition, dict) and composition:
        return composition

    model_cfg = cfg.get("model") or {}
    if len(components_cfg) != 1:
        raise ValueError(
            "Multi-component viscosity configs must include a composition: block."
        )
    name = str(components_cfg[0]["name"])
    return {
        "mode": "counts",
        "counts": {
            name: int(model_cfg.get("n_molecules", 500)),
        },
    }


def _build_per_component(
    components_cfg: list[dict[str, Any]],
    workdir: Path,
) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for comp_cfg in components_cfg:
        name = comp_cfg["name"]
        preset = comp_cfg.get("preset")
        smiles = comp_cfg.get("smiles")
        charge_source = str(comp_cfg.get("charge_source", "openff_am1bcc")).lower()
        print(f"  [{name}] generating geometry and charges")

        if preset:
            xyz_atoms, charges, mw = load_ion_preset(preset)
            method = f"preset:{preset}"
        else:
            if not smiles:
                raise ValueError(
                    f"component {name!r}: must specify either 'preset' or 'smiles'"
                )
            mol, xyz_atoms, mw = generate_molecule(smiles, name)
            if charge_source in ("openff_am1bcc", "openff"):
                charges, method = get_openff_charges(mol, smiles, workdir=workdir)
            elif charge_source in ("zero", "zeros"):
                charges = [0.0] * len(xyz_atoms)
                method = "zero"
            elif charge_source in ("external", "external_file"):
                charges = _load_external_charges(comp_cfg, expected=len(xyz_atoms))
                method = "external"
            else:
                raise ValueError(
                    f"component {name!r}: unsupported charge_source "
                    f"{charge_source!r}. Use 'openff_am1bcc', 'zero', "
                    "or 'external'."
                )

        xyz_path = workdir / f"{name}.xyz"
        write_xyz(xyz_path, xyz_atoms, comment=f"{name} ({smiles or preset})")
        net = float(sum(charges))
        print(
            f"    atoms={len(xyz_atoms):3d}  MW={mw:7.3f}  "
            f"net_charge={net:+.4f}  ({method})"
        )
        out.append(
            {
                "cfg": comp_cfg,
                "name": name,
                "xyz_path": xyz_path,
                "xyz_atoms": xyz_atoms,
                "charges": charges,
                "mw": mw,
                "charge_method": method,
            }
        )
    return out


def _load_xyz_coords(xyz: Path) -> np.ndarray:
    """Read Packmol output XYZ into an ``(N, 3)`` coordinate array."""
    return np.loadtxt(xyz, skiprows=2, usecols=(1, 2, 3))


def _assemble_atom_records(
    per_component: list[dict[str, Any]],
    resolved: list[ResolvedComponent],
    system_coords: np.ndarray,
) -> tuple[list[AtomRecord], int]:
    """Build ordered atom records with molecule ids and charges."""
    expected_atoms = sum(
        r.count * len(p["xyz_atoms"]) for p, r in zip(per_component, resolved)
    )
    if system_coords.shape[0] != expected_atoms:
        raise RuntimeError(
            f"Packmol output has {system_coords.shape[0]} atoms, but "
            f"bookkeeping expected {expected_atoms}."
        )

    records: list[AtomRecord] = []
    cursor = 0
    mol_id = 0
    for p, r in zip(per_component, resolved):
        per_mol_atoms = p["xyz_atoms"]
        charges = p["charges"]
        for _ in range(r.count):
            mol_id += 1
            for i, (elem, *_xyz) in enumerate(per_mol_atoms):
                x, y, z = system_coords[cursor]
                records.append(
                    (elem, float(charges[i]), float(x), float(y), float(z), mol_id)
                )
                cursor += 1
    return records, mol_id


def _load_external_charges(
    comp_cfg: dict[str, Any],
    *,
    expected: int,
) -> list[float]:
    """Read per-atom charges from JSON or plain text."""
    path_str: str | None = None
    charge_source = comp_cfg.get("charge_source", "")
    if isinstance(charge_source, str) and ":" in charge_source:
        prefix, rest = charge_source.split(":", 1)
        if prefix.lower() in ("external", "external_file"):
            path_str = rest.strip()
    if path_str is None:
        path_str = comp_cfg.get("charge_file")
    if not path_str:
        raise ValueError(
            f"component {comp_cfg.get('name')!r}: external charges requested "
            "but no path provided."
        )

    path = Path(path_str).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"External charge file not found: {path}")

    if path.suffix.lower() == ".json":
        payload = json.loads(path.read_text())
        if isinstance(payload, dict):
            values = payload.get("charges") or payload.get("partial_charges")
        elif isinstance(payload, list):
            values = payload
        else:
            values = None
        if values is None:
            raise RuntimeError(f"Unsupported JSON charge file: {path}")
        charges = [float(value) for value in values]
    else:
        charges = []
        for raw in path.read_text().splitlines():
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            charges.append(float(line.split()[0]))

    if len(charges) != expected:
        raise RuntimeError(
            f"External charge count mismatch for {comp_cfg.get('name')!r}: "
            f"got {len(charges)}, expected {expected}"
        )
    return charges


__all__ = [
    "BambooSystemResult",
    "build_bamboo_system",
    "load_bamboo_build_report",
]
