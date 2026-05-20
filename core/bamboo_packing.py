"""Multi-component Packmol box builder.

Adapted from ``lmp-work/core/packing.py``, generalised to accept a list of
components each with its own single-molecule XYZ template and copy count. The
box size is sized from the *total* mass of all components and the requested
initial-pack density.
"""
from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from utils.constants import AVOGADRO


@dataclass
class PackmolComponent:
    """One species in the packed box."""

    name: str
    xyz_file: Path
    count: int
    n_atoms_per_molecule: int
    molecular_weight: float  # g/mol


class PackmolBuilder:
    """Build a multi-component liquid box with Packmol."""

    _SUCCESS_TOKEN = "Success!"
    _IMPERFECT_TOKENS = (
        "ENDED WITHOUT PERFECT PACKING",
        "Packmol was not able to find a solution",
    )

    def __init__(self, workdir: Path):
        self.workdir = Path(workdir)
        self.workdir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # public API
    # ------------------------------------------------------------------

    def build_box(
        self,
        components: list[PackmolComponent],
        *,
        density_g_cm3: float,
        density_scale: float = 0.85,
        tolerance: float = 2.0,
        seed: int = 192911,
        max_attempts: int = 5,
        seed_step: int = 17,
        nloop: int | None = 200,
        full_box: bool = False,
        margin: float | None = None,
        allow_imperfect: bool = True,
        accept_forced_output: bool = True,
        box_aspect_ratio: tuple[float, float, float] | None = None,
        box_lengths: tuple[float, float, float] | None = None,
    ) -> tuple[Path, tuple[float, float, float]]:
        """Pack the components into a cuboid liquid box.

        ``density_g_cm3`` is the *target* density used for sizing the box;
        Packmol actually runs at ``density × density_scale`` to leave slack for
        clash resolution. A subsequent NPT equilibration is expected to
        compress the system toward the target density.
        """
        if not components:
            raise ValueError("Need at least one component to pack.")
        if density_g_cm3 <= 0:
            raise ValueError("density_g_cm3 must be > 0")
        if density_scale <= 0:
            raise ValueError("density_scale must be > 0")
        if tolerance <= 0:
            raise ValueError("tolerance must be > 0")
        if max_attempts < 1:
            raise ValueError("max_attempts must be >= 1")

        box_dims = self._compute_box_lengths(
            components=components,
            density_g_cm3=density_g_cm3 * density_scale,
            box_aspect_ratio=box_aspect_ratio,
            box_lengths=box_lengths,
        )
        lx, ly, lz = box_dims

        if margin is None:
            margin = 0.0 if full_box else tolerance
        if margin < 0:
            raise ValueError("margin must be >= 0")
        inner = (lx - 2 * margin, ly - 2 * margin, lz - 2 * margin)
        if not full_box and min(inner) <= 0:
            raise RuntimeError(
                f"Box too small ({lx:.2f}, {ly:.2f}, {lz:.2f} A) "
                f"for the requested margin {margin:.2f} A"
            )

        print("=" * 70)
        print("Packmol: building multi-component liquid box")
        print("=" * 70)
        for c in components:
            print(
                f"  {c.name:<8s}  n={c.count:5d}  "
                f"MW={c.molecular_weight:7.3f} g/mol  "
                f"atoms/mol={c.n_atoms_per_molecule}"
            )
        print(f"  target density (post-NPT) : {density_g_cm3:.3f} g/cm^3")
        print(f"  initial pack density      : {density_g_cm3 * density_scale:.3f} g/cm^3")
        print(f"  box (A)                   : {lx:.3f} x {ly:.3f} x {lz:.3f}")
        print(f"  tolerance / margin (A)    : {tolerance:.2f} / {margin:.2f}")
        print()

        system_xyz = self.workdir / "system.xyz"
        canonical_inp = self.workdir / "packmol.inp"
        canonical_log = self.workdir / "packmol.log"

        attempts_summary: list[str] = []
        accepted: dict | None = None
        chosen_inp: Path | None = None
        chosen_log: Path | None = None

        for i in range(max_attempts):
            attempt_id = i + 1
            seed_i = seed + i * seed_step
            if max_attempts == 1:
                attempt_inp, attempt_log = canonical_inp, canonical_log
            else:
                attempt_inp = self.workdir / f"packmol_try{attempt_id}.inp"
                attempt_log = self.workdir / f"packmol_try{attempt_id}.log"

            self._write_input(
                attempt_inp, system_xyz, components, box_dims, margin, inner,
                tolerance, seed_i, full_box=full_box, nloop=nloop,
            )
            run = self._run(
                attempt_inp, attempt_log, system_xyz,
                accept_forced_output=accept_forced_output,
            )
            attempts_summary.append(
                f"attempt={attempt_id} seed={seed_i} rc={run['returncode']} "
                f"status={run['status']} output={run['output_exists']}"
            )

            if run["usable"] and (run["perfect"] or run["imperfect"] or allow_imperfect):
                accepted, chosen_inp, chosen_log = run, attempt_inp, attempt_log
                if run["perfect"]:
                    break

        if accepted is None:
            raise RuntimeError(
                "Packmol failed for every attempt:\n  - "
                + "\n  - ".join(attempts_summary)
            )

        if max_attempts > 1:
            assert chosen_inp is not None and chosen_log is not None
            shutil.copy2(chosen_inp, canonical_inp)
            shutil.copy2(chosen_log, canonical_log)

        if accepted["perfect"]:
            tag = "perfect"
        elif accepted["imperfect"]:
            tag = "imperfect (accepted as best-effort)"
        else:
            tag = "unknown (output present, accepted)"
        print(f"  Packmol status: {tag}")
        if accepted["used_forced_output"]:
            print("  note: used Packmol's *_FORCED output as fallback")

        self._verify(system_xyz, components)
        return system_xyz, box_dims

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _compute_box_lengths(
        *,
        components: list[PackmolComponent],
        density_g_cm3: float,
        box_aspect_ratio: tuple[float, float, float] | None,
        box_lengths: tuple[float, float, float] | None,
    ) -> tuple[float, float, float]:
        if box_lengths is not None:
            lx, ly, lz = (float(v) for v in box_lengths)
            if min(lx, ly, lz) <= 0:
                raise ValueError("box_lengths must all be > 0")
            return lx, ly, lz

        # mass [g] = (Σ count × MW) / N_A   →   V [cm³] = m / ρ   →   convert to A³
        mass_total_g = sum(c.molecular_weight * c.count for c in components) / AVOGADRO
        volume_a3 = (mass_total_g / density_g_cm3) * 1e24

        if box_aspect_ratio is None:
            box_aspect_ratio = (1.0, 1.0, 1.0)
        rx, ry, rz = (float(v) for v in box_aspect_ratio)
        if min(rx, ry, rz) <= 0:
            raise ValueError("box_aspect_ratio components must all be > 0")
        scale = (volume_a3 / (rx * ry * rz)) ** (1.0 / 3.0)
        return rx * scale, ry * scale, rz * scale

    @staticmethod
    def _write_input(
        inp_file: Path,
        system_xyz: Path,
        components: list[PackmolComponent],
        box_dims: tuple[float, float, float],
        margin: float,
        inner: tuple[float, float, float],
        tolerance: float,
        seed: int,
        *,
        full_box: bool,
        nloop: int | None,
    ) -> None:
        lx, ly, lz = box_dims
        if full_box:
            region = (
                f"  inside box 0.000000 0.000000 0.000000 "
                f"{lx:.6f} {ly:.6f} {lz:.6f}"
            )
        else:
            ix, iy, iz = inner
            region = (
                f"  inside box {margin:.6f} {margin:.6f} {margin:.6f} "
                f"{margin + ix:.6f} {margin + iy:.6f} {margin + iz:.6f}"
            )
        nloop_line = f"nloop {int(nloop)}\n" if nloop is not None else ""

        blocks: list[str] = []
        for c in components:
            blocks.append(
                f"structure {c.xyz_file.absolute()}\n"
                f"  number {c.count}\n"
                f"{region}\n"
                f"end structure\n"
            )

        text = (
            f"tolerance {tolerance:.6f}\n"
            f"seed {seed}\n"
            f"{nloop_line}"
            f"filetype xyz\n"
            f"output {system_xyz.absolute()}\n\n"
            + "\n".join(blocks)
        )
        inp_file.write_text(text)

    def _run(
        self,
        inp_file: Path,
        log_file: Path,
        system_xyz: Path,
        *,
        accept_forced_output: bool,
    ) -> dict:
        if system_xyz.exists():
            system_xyz.unlink()
        forced = system_xyz.parent / f"{system_xyz.name}_FORCED"
        if forced.exists():
            forced.unlink()
        if log_file.exists():
            log_file.unlink()

        try:
            with open(inp_file, "rb") as fin, open(log_file, "wb") as flog:
                result = subprocess.run(
                    ["packmol"],
                    stdin=fin, stdout=flog, stderr=subprocess.STDOUT,
                    cwd=self.workdir,
                )
        except FileNotFoundError as exc:
            raise RuntimeError(
                "`packmol` binary not found on PATH. Install Packmol and ensure "
                "it is reachable as `packmol`."
            ) from exc

        text = log_file.read_text(errors="ignore") if log_file.exists() else ""
        perfect = self._SUCCESS_TOKEN in text
        imperfect = any(t in text for t in self._IMPERFECT_TOKENS)
        if perfect:
            status = "perfect"
        elif imperfect:
            status = "imperfect"
        elif "ERROR" in text.upper():
            status = "error"
        else:
            status = "unknown"

        output_exists = system_xyz.exists()
        used_forced = False
        if not output_exists and accept_forced_output and forced.exists():
            shutil.copy2(forced, system_xyz)
            output_exists = True
            used_forced = True

        usable = output_exists and (result.returncode == 0 or perfect or imperfect)
        return {
            "returncode": result.returncode,
            "status": status,
            "perfect": perfect,
            "imperfect": imperfect,
            "usable": usable,
            "output_exists": output_exists,
            "used_forced_output": used_forced,
        }

    @staticmethod
    def _verify(system_xyz: Path, components: list[PackmolComponent]) -> None:
        expected = sum(c.count * c.n_atoms_per_molecule for c in components)
        with open(system_xyz) as fh:
            n_total = int(fh.readline().strip())
        if n_total != expected:
            raise RuntimeError(
                f"system.xyz atom count mismatch: got {n_total}, expected {expected}"
            )
        n_mol = sum(c.count for c in components)
        print(f"  packed atoms: {n_total} (= {n_mol} molecules)")


__all__ = ["PackmolBuilder", "PackmolComponent"]
