"""
OpenFF力场实现
"""
import hashlib
import json
import shutil
from pathlib import Path
from rdkit import Chem
from typing import Dict, Any
from ..core.forcefield import ForceField


class OpenFF(ForceField):
    """OpenFF Sage力场"""

    def __init__(
        self,
        workdir: Path,
        version: str = "openff-2.0.0",
        *,
        strict_stereo: bool = True,
        charge_method: str = "am1bcc",
        charge_fallback: str | None = "gasteiger",
        use_cache: bool = True,
    ):
        """
        Args:
            workdir: 工作目录
            version: OpenFF力场版本
            strict_stereo: True时不允许未定义立体化学
            charge_method: 首选电荷方法
            charge_fallback: 失败时备用电荷方法（None表示不回退）
            use_cache: 是否启用参数缓存
        """
        super().__init__(workdir)
        self.version = version
        self.strict_stereo = strict_stereo
        self.charge_method = charge_method
        self.charge_fallback = charge_fallback
        self.use_cache = use_cache
        self.cache_dir = self.workdir / ".openff_cache"
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def generate_parameters(
        self,
        mol: Chem.Mol,
        smiles: str,
        name: str
    ) -> Path:
        """
        使用OpenFF生成力场参数

        Args:
            mol: RDKit分子对象
            smiles: SMILES字符串
            name: 分子名称

        Returns:
            单分子LAMMPS data文件路径
        """
        print("=" * 70)
        print(f"步骤3: 生成OpenFF力场参数（{self.version}）")
        print("=" * 70)

        try:
            from openff.toolkit.topology import Molecule, Topology
            from openff.toolkit.typing.engines.smirnoff import ForceField as OFFForceField
            from openff.units import unit
        except ImportError:
            raise RuntimeError(
                "OpenFF未安装！\n"
                "安装方法: pip install openff-toolkit openff-interchange\n"
                "或使用conda: conda install -c conda-forge openff-toolkit"
            )

        canonical_smiles = self._canonicalize_smiles(smiles)
        cache_key = self._make_cache_key(canonical_smiles)

        lmp_single = self.workdir / "single_molecule_openff.lmp"
        meta_file = self.workdir / "single_molecule_openff.meta.json"

        if self.use_cache and self._try_restore_from_cache(cache_key, lmp_single, meta_file):
            print(f"  ✓ 命中缓存参数: {lmp_single}")
            return lmp_single

        print("  创建OpenFF分子对象...")
        try:
            off_mol = Molecule.from_rdkit(
                mol,
                allow_undefined_stereo=not self.strict_stereo,
                hydrogens_are_explicit=True,
            )
        except Exception as e:
            raise RuntimeError(
                "OpenFF分子构建失败：立体化学信息不完整。\n"
                "可选方案：\n"
                "1) 在SMILES中补全手性标记；\n"
                "2) 临时设置 strict_stereo=False。"
            ) from e

        print(f"  加载力场 {self.version}.offxml ...")
        forcefield = OFFForceField(f"{self.version}.offxml")

        # 显式电荷策略：先主方法，失败再可选回退，提升跨环境可重复性。
        used_charge_method = self._assign_charges(off_mol)
        print(f"  电荷方法: {used_charge_method}")

        topology = Topology.from_molecules([off_mol])
        interchange = forcefield.create_interchange(topology)

        # 赋予单分子坐标（用于后续导出）
        import numpy as np
        rdkit_positions = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
        interchange.positions = rdkit_positions * unit.angstrom

        # 临时盒子（单分子不需要真实盒子）
        interchange.box = np.eye(3) * 50.0 * unit.angstrom

        # 导出LAMMPS格式
        lmp_base = self.workdir / f"single_molecule_openff_{cache_key}"
        lmp_double = self.workdir / "single_molecule_openff.lmp.lmp"
        lmp_tmp = self.workdir / f"single_molecule_openff_{cache_key}.lmp"

        for stale in (lmp_single, lmp_double, lmp_base, lmp_tmp):
            if stale.exists():
                stale.unlink()

        interchange.to_lammps(str(lmp_base))

        produced = None
        for candidate in (lmp_tmp, lmp_double, lmp_base):
            if candidate.exists():
                produced = candidate
                break

        if produced is None:
            raise RuntimeError("OpenFF导出失败：未生成LAMMPS文件")

        if produced != lmp_single:
            shutil.move(str(produced), str(lmp_single))

        self._validate_lammps_data_sections(lmp_single)

        meta = {
            "name": name,
            "smiles_input": smiles,
            "smiles_canonical": canonical_smiles,
            "forcefield": self.version,
            "charge_method": used_charge_method,
            "strict_stereo": self.strict_stereo,
            "cache_key": cache_key,
        }
        meta_file.write_text(json.dumps(meta, indent=2, ensure_ascii=False))
        self._save_to_cache(cache_key, lmp_single, meta)

        print(f"  ✓ OpenFF参数已生成: {lmp_single}")

        self._print_forcefield_info(mol, used_charge_method)

        return lmp_single

    def get_lammps_settings(self) -> Dict[str, Any]:
        """
        获取OpenFF的LAMMPS设置

        Returns:
            LAMMPS设置字典
        """
        return {
            'pair_style': 'lj/cut/coul/long 12.0',
            'pair_modify': 'mix arithmetic tail yes',  # OpenFF用算术平均
            'bond_style': 'harmonic',
            'angle_style': 'harmonic',
            'dihedral_style': 'fourier',  # OpenFF通常用fourier
            'improper_style': 'cvff',
            'special_bonds': {
                'lj': [0.0, 0.0, 0.5],
                'coul': [0.0, 0.0, 0.8333333333]  # OpenFF: 1-4 Coul = 1/1.2
            },
            'kspace_style': 'pppm 1e-5',
        }

    @property
    def name(self) -> str:
        return f"OpenFF-{self.version}"

    @property
    def citation(self) -> str:
        return (
            "Boothroyd, S., et al. (2023). "
            "Development and Benchmarking of Open Force Field 2.0.0: "
            "The Sage Small Molecule Force Field. "
            "Journal of Chemical Theory and Computation, 19(11), 3251-3275."
        )

    def _canonicalize_smiles(self, smiles: str) -> str:
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            raise ValueError(f"无效SMILES: {smiles}")
        return Chem.MolToSmiles(m, canonical=True, isomericSmiles=True)

    def _make_cache_key(self, canonical_smiles: str) -> str:
        payload = {
            "smiles": canonical_smiles,
            "version": self.version,
            "strict_stereo": self.strict_stereo,
            "charge_method": self.charge_method,
            "charge_fallback": self.charge_fallback,
        }
        raw = json.dumps(payload, sort_keys=True, ensure_ascii=True)
        return hashlib.sha256(raw.encode("utf-8")).hexdigest()[:16]

    def _cache_paths(self, cache_key: str) -> tuple[Path, Path]:
        return (
            self.cache_dir / f"single_molecule_openff_{cache_key}.lmp",
            self.cache_dir / f"single_molecule_openff_{cache_key}.meta.json",
        )

    def _try_restore_from_cache(self, cache_key: str, target_lmp: Path, target_meta: Path) -> bool:
        cache_lmp, cache_meta = self._cache_paths(cache_key)
        if not cache_lmp.exists():
            return False
        if not self._is_valid_lammps_data(cache_lmp):
            return False
        shutil.copy2(cache_lmp, target_lmp)
        if cache_meta.exists():
            shutil.copy2(cache_meta, target_meta)
        return True

    def _save_to_cache(self, cache_key: str, source_lmp: Path, meta: Dict[str, Any]) -> None:
        cache_lmp, cache_meta = self._cache_paths(cache_key)
        shutil.copy2(source_lmp, cache_lmp)
        cache_meta.write_text(json.dumps(meta, indent=2, ensure_ascii=False))

    def _assign_charges(self, off_mol) -> str:
        methods = [self.charge_method]
        if self.charge_fallback and self.charge_fallback not in methods:
            methods.append(self.charge_fallback)

        errors = []
        for method in methods:
            try:
                off_mol.assign_partial_charges(partial_charge_method=method)
                return method
            except Exception as e:
                errors.append(f"{method}: {e}")

        msg = "\n".join(errors)
        raise RuntimeError(
            "OpenFF电荷赋值失败。\n"
            "请检查可用toolkit后端（如 AmberTools/OpenEye）或更换 charge_method。\n"
            f"尝试记录:\n{msg}"
        )

    def _is_valid_lammps_data(self, path: Path) -> bool:
        try:
            self._validate_lammps_data_sections(path)
            return True
        except Exception:
            return False

    @staticmethod
    def _validate_lammps_data_sections(path: Path) -> None:
        text = path.read_text(errors="ignore")
        required = ("Masses", "Pair Coeffs", "Atoms", "Bonds", "Angles", "Dihedrals")
        missing = [s for s in required if s not in text]
        if missing:
            raise RuntimeError(f"LAMMPS data文件缺少必要section: {missing}")

    def _print_forcefield_info(self, mol: Chem.Mol, charge_method: str):
        """打印力场信息"""
        print(f"\n  力场信息:")
        print(f"    力场: OpenFF {self.version}")
        print(f"    原子类型数: {len(set([a.GetAtomicNum() for a in mol.GetAtoms()]))}")
        print(f"    键数: {mol.GetNumBonds()}")
        print(f"    电荷方法: {charge_method}")
        print(f"    立体化学严格模式: {self.strict_stereo}")
        print(f"    参数缓存: {'ON' if self.use_cache else 'OFF'}")
        print(f"    混合规则: arithmetic mean")
        print(f"    1-4缩放: LJ=0.5, Coul=0.833")
        print()
