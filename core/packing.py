"""
Packmol liquid-box builder.
"""
import shutil
import subprocess
from pathlib import Path
from ..utils.constants import AVOGADRO


class PackmolBuilder:
    """Builder for generating a liquid simulation box using Packmol."""

    _SUCCESS_TOKEN = "Success!"
    _IMPERFECT_TOKENS = (
        "ENDED WITHOUT PERFECT PACKING",
        "Packmol was not able to find a solution",
    )

    def __init__(self, workdir: Path):
        """
        Args:
            workdir: Working directory.
        """
        self.workdir = Path(workdir)
        self.workdir.mkdir(parents=True, exist_ok=True)

    def build_box(
        self,
        mol_xyz: Path,
        n_molecules: int,
        density: float,
        mol_weight: float,
        tolerance: float = 2.0,
        seed: int = 192911,
        max_attempts: int = 3,
        seed_step: int = 97,
        nloop: int | None = None,
        full_box: bool = False,
        margin: float | None = None,
        allow_imperfect: bool = True,
        accept_forced_output: bool = True,
    ) -> tuple[Path, float]:
        """
        Build a liquid-phase box with Packmol.

        Args:
            mol_xyz: Single-molecule XYZ file.
            n_molecules: Number of molecules.
            density: Target density (g/cm^3).
            mol_weight: Molecular weight (g/mol).
            tolerance: Packmol tolerance (A).
            seed: Packmol random seed.
            max_attempts: Number of attempts with different seeds.
            seed_step: Seed increment between attempts.
            nloop: Optional Packmol nloop value.
            full_box: Use full box "inside box 0 0 0 L L L" when True.
            margin: Margin for "inside cube" mode. Defaults to tolerance.
            allow_imperfect: Accept "ENDED WITHOUT PERFECT PACKING" result.
            accept_forced_output: Accept *_FORCED as fallback output.

        Returns:
            Tuple of (system_xyz_path, box_length).
        """
        if n_molecules <= 0:
            raise ValueError("n_molecules 必须 > 0")
        if density <= 0:
            raise ValueError("density 必须 > 0")
        if mol_weight <= 0:
            raise ValueError("mol_weight 必须 > 0")
        if tolerance <= 0:
            raise ValueError("tolerance 必须 > 0")
        if max_attempts < 1:
            raise ValueError("max_attempts 必须 >= 1")
        if max_attempts > 1 and seed_step == 0:
            raise ValueError("max_attempts>1 时 seed_step 不能为 0")
        if nloop is not None and nloop < 1:
            raise ValueError("nloop 必须 >= 1")

        print("=" * 70)
        print("步骤2: 构建液相盒子（Packmol）")
        print("=" * 70)

        # Compute box size from density
        box_length = self._calculate_box_size(n_molecules, density, mol_weight)

        if margin is None:
            margin = 0.0 if full_box else tolerance
        if margin < 0:
            raise ValueError("margin 必须 >= 0")

        print("  系统参数:")
        print(f"    分子数: {n_molecules}")
        print(f"    目标密度: {density:.3f} g/cm³")
        print(f"    盒子边长: {box_length:.3f} Å")
        print(f"    初始随机种子: {seed}")
        print(f"    尝试次数: {max_attempts}")
        if max_attempts > 1:
            print(f"    seed步长: {seed_step}")
        print(f"    打包模式: {'full box' if full_box else 'inner cube'}")
        if not full_box:
            print(f"    边界余量: {margin:.3f} Å")
        if nloop is not None:
            print(f"    nloop: {nloop}")

        # Prepare Packmol input/output paths
        system_xyz = self.workdir / "system.xyz"
        packmol_inp = self.workdir / "packmol.inp"
        packmol_log = self.workdir / "packmol.log"

        inner_side = box_length - 2.0 * margin

        if not full_box and inner_side <= 0:
            raise RuntimeError(
                f"盒子边长过小({box_length:.3f} Å)，"
                f"无法预留边界余量({margin:.3f} Å)"
            )

        attempt_summaries = []
        accepted = None
        selected_inp = None
        selected_log = None

        for i in range(max_attempts):
            attempt_id = i + 1
            current_seed = seed + i * seed_step

            if max_attempts == 1:
                attempt_inp = packmol_inp
                attempt_log = packmol_log
            else:
                attempt_inp = self.workdir / f"packmol_try{attempt_id}.inp"
                attempt_log = self.workdir / f"packmol_try{attempt_id}.log"

            # Write Packmol input for this attempt
            self._write_packmol_input(
                attempt_inp,
                mol_xyz,
                system_xyz,
                n_molecules,
                box_length,
                margin,
                inner_side,
                tolerance,
                current_seed,
                full_box=full_box,
                nloop=nloop,
            )

            # Run Packmol
            run_result = self._run_packmol(
                attempt_inp,
                attempt_log,
                system_xyz,
                accept_forced_output=accept_forced_output,
            )

            attempt_summaries.append(
                f"attempt={attempt_id}, seed={current_seed}, returncode={run_result['returncode']}, "
                f"status={run_result['status']}, output={run_result['output_exists']}"
            )

            if run_result["usable"] and run_result["perfect"]:
                accepted = run_result
                selected_inp = attempt_inp
                selected_log = attempt_log
                break

            if run_result["usable"] and run_result["imperfect"] and allow_imperfect:
                accepted = run_result
                selected_inp = attempt_inp
                selected_log = attempt_log
                break

            if run_result["usable"] and run_result["status"] == "unknown" and allow_imperfect:
                accepted = run_result
                selected_inp = attempt_inp
                selected_log = attempt_log
                break

            if attempt_id < max_attempts:
                print(
                    f"  [重试 {attempt_id}/{max_attempts}] seed={current_seed} 未获得可用构型 "
                    f"(status={run_result['status']})，继续尝试..."
                )

        if accepted is None:
            summary = "\n    - " + "\n    - ".join(attempt_summaries)
            error_msg = (
                "\n  Packmol失败！所有尝试均未生成可用结构。"
                f"\n  尝试摘要:{summary}"
            )
            last_log = attempt_log if 'attempt_log' in locals() else packmol_log
            if last_log.exists():
                error_msg += f"\n\n最后一次日志: {last_log}\n\n{last_log.read_text(errors='ignore')}"
            raise RuntimeError(error_msg)

        # Keep canonical names for downstream tooling
        if max_attempts > 1:
            shutil.copy2(selected_inp, packmol_inp)
            shutil.copy2(selected_log, packmol_log)

        if accepted["perfect"]:
            print("    ✓ Packmol状态: perfect packing")
        elif accepted["imperfect"]:
            print("    ⚠ Packmol状态: imperfect packing（best-effort）")
        else:
            print("    ⚠ Packmol状态: unknown（已生成输出并继续）")
        if accepted["used_forced_output"]:
            print("    ⚠ 使用 *_FORCED 输出回填为 system.xyz")
        print(f"    ✓ Packmol日志: {packmol_log}")

        # Validate output
        self._verify_output(system_xyz, n_molecules, mol_xyz)

        print(f"    ✓ 输出: {system_xyz}\n")
        return system_xyz, box_length

    @staticmethod
    def _calculate_box_size(n_molecules: int, density: float, mol_weight: float) -> float:
        """
        Compute cubic box length.

        Args:
            n_molecules: Number of molecules.
            density: Density (g/cm^3).
            mol_weight: Molecular weight (g/mol).

        Returns:
            Box length (A).
        """
        mass_total = mol_weight * n_molecules / AVOGADRO  # g
        volume_cm3 = mass_total / density  # cm³
        box_length = (volume_cm3 * 1e24) ** (1.0 / 3.0)  # Å
        return box_length

    @staticmethod
    def _write_packmol_input(
        inp_file: Path,
        mol_xyz: Path,
        output_xyz: Path,
        n_molecules: int,
        box_length: float,
        margin: float,
        inner_side: float,
        tolerance: float,
        seed: int,
        *,
        full_box: bool = False,
        nloop: int | None = None,
    ):
        """Write the Packmol input file."""
        if full_box:
            region_line = (
                f"  inside box 0.000000 0.000000 0.000000 "
                f"{box_length:.6f} {box_length:.6f} {box_length:.6f}"
            )
        else:
            region_line = (
                f"  inside cube {margin:.6f} {margin:.6f} {margin:.6f} {inner_side:.6f}"
            )

        nloop_line = f"nloop {int(nloop)}\n" if nloop is not None else ""

        packmol_input = f"""tolerance {tolerance:.6f}
seed {seed}
{nloop_line}filetype xyz
output {output_xyz.absolute()}

structure {mol_xyz.absolute()}
  number {n_molecules}
{region_line}
end structure
"""
        inp_file.write_text(packmol_input)

    def _run_packmol(
        self,
        inp_file: Path,
        log_file: Path,
        output_xyz: Path,
        *,
        accept_forced_output: bool = True,
    ) -> dict:
        """Run Packmol and return status summary."""
        print("\n  运行Packmol...")

        # Remove stale outputs
        if output_xyz.exists():
            output_xyz.unlink()
        forced_xyz = output_xyz.parent / f"{output_xyz.name}_FORCED"
        if forced_xyz.exists():
            forced_xyz.unlink()
        if log_file.exists():
            log_file.unlink()

        try:
            with open(inp_file, "rb") as fin, open(log_file, "wb") as flog:
                result = subprocess.run(
                    ["packmol"],
                    stdin=fin,
                    stdout=flog,
                    stderr=subprocess.STDOUT,
                    cwd=self.workdir,
                )
        except FileNotFoundError as e:
            raise RuntimeError(
                "未找到 packmol 可执行文件，请确认已安装并在 PATH 中。"
            ) from e

        log_text = log_file.read_text(errors="ignore") if log_file.exists() else ""
        status = self._parse_packmol_status(log_text)

        output_exists = output_xyz.exists()
        used_forced_output = False
        if not output_exists and accept_forced_output and forced_xyz.exists():
            shutil.copy2(forced_xyz, output_xyz)
            output_exists = True
            used_forced_output = True

        # Accept known non-perfect but usable Packmol outputs.
        usable = output_exists and (
            result.returncode == 0 or status["perfect"] or status["imperfect"]
        )

        return {
            "returncode": result.returncode,
            "status": status["status"],
            "perfect": status["perfect"],
            "imperfect": status["imperfect"],
            "usable": usable,
            "output_exists": output_exists,
            "used_forced_output": used_forced_output,
            "log_file": log_file,
        }

    def _parse_packmol_status(self, log_text: str) -> dict:
        """Parse Packmol status from log text."""
        perfect = self._SUCCESS_TOKEN in log_text
        imperfect = any(token in log_text for token in self._IMPERFECT_TOKENS)

        if perfect:
            status = "perfect"
        elif imperfect:
            status = "imperfect"
        elif "ERROR" in log_text.upper():
            status = "error"
        else:
            status = "unknown"

        return {"perfect": perfect, "imperfect": imperfect, "status": status}

    @staticmethod
    def _verify_output(system_xyz: Path, n_molecules: int, mol_xyz: Path):
        """Validate output file integrity."""
        # Read atom count for single molecule
        with open(mol_xyz) as f:
            n_atoms_per_mol = int(f.readline().strip())

        # Read total atom count in the packed system
        with open(system_xyz) as f:
            n_atoms_total = int(f.readline().strip())

        expected_atoms = n_atoms_per_mol * n_molecules
        if n_atoms_total != expected_atoms:
            raise RuntimeError(
                f"原子数不匹配: {n_atoms_total} != {expected_atoms}"
            )

        print(f"    ✓ 原子数验证: {n_atoms_total} = {n_molecules} × {n_atoms_per_mol}")
