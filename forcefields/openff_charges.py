"""OpenFF AM1-BCC partial-charge generator (lightweight wrapper).

We extract only the per-atom partial charges; bonds / angles / dihedrals / LJ
parameters that the OpenFF interchange would normally compute are intentionally
dropped because BAMBOO learns those internally. Long-range electrostatics in
the LAMMPS run is handled by PPPM acting on the charges we write into in.data.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from rdkit import Chem


def get_openff_charges(
    rdkit_mol: Chem.Mol,
    smiles: str,
    *,
    workdir: Path | None = None,
    charge_method: str = "am1bcc",
    fallback: str | None = "gasteiger",
    use_cache: bool = True,
) -> tuple[list[float], str]:
    """Compute partial charges for *rdkit_mol* using the OpenFF toolkit.

    Parameters
    ----------
    rdkit_mol:
        RDKit Mol object with explicit hydrogens (must match the XYZ that will
        be passed to Packmol — charges are returned in this RDKit atom order).
    smiles:
        Original SMILES used only as a cache key.
    workdir:
        Directory in which to keep ``<workdir>/.charge_cache/``. ``None``
        disables caching entirely.
    charge_method:
        Primary OpenFF ``partial_charge_method`` (e.g. ``am1bcc``).
    fallback:
        Optional secondary method tried when the primary raises. Set to
        ``None`` to disable fallback.
    use_cache:
        Cache by canonical SMILES + method. Re-runs hit the cache.

    Returns
    -------
    charges:
        Per-atom partial charges (length == ``rdkit_mol.GetNumAtoms()``) in
        elementary-charge units, in RDKit atom order.
    method_used:
        The method that actually produced the returned charges.
    """
    try:
        from openff.toolkit.topology import Molecule
        from openff.units import unit  # noqa: F401 — sanity check
    except ImportError as exc:
        raise RuntimeError(
            "openff-toolkit is not installed.\n"
            "  pip install openff-toolkit openff-units\n"
            "  (and ambertools available so AM1-BCC can call sqm)"
        ) from exc

    # Cache key uses canonical SMILES (without Hs to be RDKit-canonical-safe).
    canonical = Chem.MolToSmiles(
        Chem.RemoveHs(Chem.Mol(rdkit_mol)),
        canonical=True,
        isomericSmiles=True,
    )

    cache_file: Path | None = None
    if use_cache and workdir is not None:
        cache_dir = Path(workdir) / ".charge_cache"
        cache_dir.mkdir(parents=True, exist_ok=True)
        key_raw = json.dumps(
            {"smiles": canonical, "method": charge_method, "fallback": fallback},
            sort_keys=True,
        )
        key = hashlib.sha256(key_raw.encode("utf-8")).hexdigest()[:16]
        cache_file = cache_dir / f"charges_{key}.json"
        if cache_file.exists():
            payload = json.loads(cache_file.read_text())
            if len(payload.get("charges", [])) == rdkit_mol.GetNumAtoms():
                return list(map(float, payload["charges"])), str(payload["method"])

    off_mol = Molecule.from_rdkit(
        rdkit_mol,
        allow_undefined_stereo=True,
        hydrogens_are_explicit=True,
    )

    methods = [charge_method]
    if fallback and fallback not in methods:
        methods.append(fallback)

    errors: list[str] = []
    method_used: str | None = None
    for method in methods:
        try:
            off_mol.assign_partial_charges(partial_charge_method=method)
            method_used = method
            break
        except Exception as exc:  # noqa: BLE001 — we report every failure together
            errors.append(f"{method}: {exc}")

    if method_used is None:
        raise RuntimeError(
            "Partial-charge assignment failed for all methods:\n  "
            + "\n  ".join(errors)
        )

    charges_raw = off_mol.partial_charges.magnitude.tolist()
    charges = [float(c) for c in charges_raw]

    # Re-center any tiny numerical drift onto the molecule's formal charge so
    # the data file is exactly net-integer-charge (PPPM neutrality check later
    # is on the whole-box scale).
    target = float(sum(atom.GetFormalCharge() for atom in rdkit_mol.GetAtoms()))
    drift = sum(charges) - target
    if abs(drift) > 1e-8 and len(charges) > 0:
        offset = drift / len(charges)
        charges = [c - offset for c in charges]

    if cache_file is not None:
        cache_file.write_text(
            json.dumps(
                {"charges": charges, "method": method_used, "smiles": canonical},
                indent=2,
            )
        )

    return charges, method_used


__all__ = ["get_openff_charges"]
