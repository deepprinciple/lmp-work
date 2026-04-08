"""
LigParGen/BOSS force-field implementation.
"""
from __future__ import annotations

import copy
import json
import subprocess
from pathlib import Path
from typing import Any, Dict

from rdkit import Chem

from ..core.forcefield import ForceField


LIGPARGEN_LAMMPS_SETTINGS: Dict[str, Any] = {
    "pair_style": "lj/cut/coul/long 12.0",
    "pair_modify": "mix geometric tail yes",
    "bond_style": "harmonic",
    "angle_style": "harmonic",
    "dihedral_style": "opls",
    "improper_style": "cvff",
    "special_bonds": {
        "lj": [0.0, 0.0, 0.5],
        "coul": [0.0, 0.0, 0.5],
    },
    "kspace_style": "pppm 1e-5",
}


def get_ligpargen_lammps_settings() -> Dict[str, Any]:
    """Return a copy of the default LigParGen/OPLS LAMMPS settings."""
    return copy.deepcopy(LIGPARGEN_LAMMPS_SETTINGS)


class LigParGen(ForceField):
    """LigParGen/BOSS OPLS-AA small-molecule force field."""

    def __init__(
        self,
        workdir: Path,
        *,
        charge_method: str = "CM1A-LBCC",
        input_mode: str = "generated_pdb",
        residue_name: str = "MOL",
        n_optimizations: int = 0,
        wrapper_script: str | None = None,
        debug: bool = False,
    ):
        super().__init__(workdir)
        self.charge_method = self._normalize_charge_method(charge_method)
        self.input_mode = str(input_mode).strip().lower()
        self.residue_name = str(residue_name).strip() or "MOL"
        self.n_optimizations = int(n_optimizations)
        self.wrapper_script = Path(
            wrapper_script or "/root/lmp-work/.compat/bin/ligpargen_local.sh"
        )
        self.debug = bool(debug)
        self.output_dir = self.workdir / "ligpargen"
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def generate_parameters(
        self,
        mol: Chem.Mol,
        smiles: str,
        name: str,
    ) -> Path:
        """
        Generate a single-molecule LAMMPS data file using LigParGen/BOSS.
        """
        if not self.wrapper_script.exists():
            raise FileNotFoundError(
                "LigParGen wrapper script not found: "
                f"{self.wrapper_script}. "
                "Please provide --forcefield_wrapper_script pointing to your "
                "locally licensed LigParGen/BOSS runtime."
            )

        formal_charge = int(
            sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
        )
        if self.charge_method == "CM1A-LBCC" and formal_charge != 0:
            raise ValueError(
                "CM1A-LBCC is only supported for neutral molecules in LigParGen/BOSS"
            )

        print("=" * 70)
        print("步骤3: 生成 LigParGen/BOSS OPLS 参数")
        print("=" * 70)
        print(f"  charge method : {self.charge_method}")
        print(f"  input mode    : {self.input_mode}")
        print(f"  wrapper       : {self.wrapper_script}")

        input_file = None
        command = [str(self.wrapper_script)]
        if self.input_mode == "smiles":
            command += ["-s", smiles]
        elif self.input_mode in {"generated_pdb", "pdb"}:
            input_file = self.workdir / f"{name}.ligpargen_input.pdb"
            Chem.MolToPDBFile(mol, str(input_file))
            command += ["-i", str(input_file)]
        else:
            raise ValueError(
                "LigParGen input_mode must be one of: smiles, generated_pdb"
            )

        command += [
            "-n",
            name,
            "-p",
            str(self.output_dir),
            "-r",
            self.residue_name,
            "-c",
            str(formal_charge),
            "-o",
            str(self.n_optimizations),
            "-cgen",
            self.charge_method,
        ]
        if self.debug:
            command.append("-debug")

        try:
            completed = subprocess.run(
                command,
                cwd=self.workdir,
                check=True,
                text=True,
                capture_output=True,
            )
        except subprocess.CalledProcessError as exc:
            stderr = exc.stderr.strip() if exc.stderr else ""
            stdout = exc.stdout.strip() if exc.stdout else ""
            detail = stderr or stdout or str(exc)
            raise RuntimeError(
                "LigParGen/BOSS parameterization failed.\n"
                f"Command: {' '.join(command)}\n"
                f"Details: {detail}"
            ) from exc

        output_file = self.output_dir / f"{name}.lammps.lmp"
        if not output_file.exists():
            raise FileNotFoundError(
                "LigParGen finished but the expected LAMMPS data file is missing: "
                f"{output_file}"
            )

        meta = {
            "engine": "ligpargen",
            "charge_method": self.charge_method,
            "input_mode": self.input_mode,
            "residue_name": self.residue_name,
            "n_optimizations": self.n_optimizations,
            "wrapper_script": str(self.wrapper_script),
            "formal_charge": formal_charge,
            "input_file": str(input_file) if input_file else None,
            "stdout_tail": (
                completed.stdout.strip().splitlines()[-20:]
                if completed.stdout
                else []
            ),
        }
        (self.output_dir / "ligpargen_meta.json").write_text(
            json.dumps(meta, indent=2, ensure_ascii=False)
        )

        print(f"  ✓ LigParGen/BOSS 参数已生成: {output_file}")
        return output_file

    def get_lammps_settings(self) -> Dict[str, Any]:
        return get_ligpargen_lammps_settings()

    @property
    def name(self) -> str:
        return "LigParGen-BOSS-OPLSAA"

    @property
    def citation(self) -> str:
        return (
            "Dodda, L. S.; Cabeza de Vaca, I.; Tirado-Rives, J.; Jorgensen, W. L. "
            "LigParGen web server: an automatic OPLS-AA parameter generator for "
            "organic ligands. Nucleic Acids Research 2017, 45(W1), W331-W336. "
            "Localized Bond-Charge Corrected CM1A Charges for Condensed-Phase "
            "Simulations."
        )

    def _normalize_charge_method(self, charge_method: str) -> str:
        token = str(charge_method).strip().upper()
        aliases = {
            "CM1A": "CM1A",
            "CM1A-LBCC": "CM1A-LBCC",
            "LBCC": "CM1A-LBCC",
        }
        if token not in aliases:
            raise ValueError(
                "LigParGen charge_method must be one of: CM1A, CM1A-LBCC"
            )
        return aliases[token]
