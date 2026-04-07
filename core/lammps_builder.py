"""
LAMMPS多分子体系构建器
"""
from pathlib import Path
import numpy as np
from typing import Dict, List
from ..utils.io import parse_lammps_data, write_lammps_data, make_lammps_header


class LAMMPSBuilder:
    """LAMMPS多分子体系构建器"""

    def __init__(self, workdir: Path):
        """
        Args:
            workdir: 工作目录
        """
        self.workdir = Path(workdir)

    def build_system(
        self,
        single_data: Path,
        system_xyz: Path,
        n_molecules: int,
        n_atoms_per_mol: int,
        box_lengths: tuple[float, float, float] | float,
        forcefield_name: str = "Unknown"
    ) -> Path:
        """
        将单分子拓扑复制到多分子体系

        Args:
            single_data: 单分子LAMMPS data文件
            system_xyz: Packmol生成的多分子坐标文件
            n_molecules: 分子数
            n_atoms_per_mol: 每个分子的原子数
            box_lengths: 盒子尺寸 (Å)，可为单个边长或 (Lx, Ly, Lz)
            forcefield_name: 力场名称

        Returns:
            多分子体系LAMMPS data文件路径
        """
        if isinstance(box_lengths, (int, float)):
            box_dims = (float(box_lengths), float(box_lengths), float(box_lengths))
        else:
            if len(box_lengths) != 3:
                raise ValueError("box_lengths 必须包含 3 个分量")
            box_dims = tuple(float(v) for v in box_lengths)

        print("=" * 70)
        print("步骤4: 导出LAMMPS data（多分子体系）")
        print("=" * 70)

        # 读取Packmol坐标
        coords = np.loadtxt(system_xyz, skiprows=2, usecols=(1, 2, 3))
        expected_atoms = n_atoms_per_mol * n_molecules

        if coords.shape[0] != expected_atoms:
            raise RuntimeError(f"坐标数不匹配: {coords.shape[0]} != {expected_atoms}")

        print(f"  分子数: {n_molecules}")
        print(f"  每分子原子数: {n_atoms_per_mol}")
        print(f"  总原子数: {expected_atoms}")
        print(f"  盒子尺寸: {box_dims[0]:.3f} x {box_dims[1]:.3f} x {box_dims[2]:.3f} Å")

        # 解析单分子data
        single = parse_lammps_data(single_data)

        # 复制拓扑
        system = self._replicate_topology(single, coords, n_molecules, n_atoms_per_mol)

        # 生成header
        system['header'] = self._make_system_header(
            single, n_molecules, box_dims, forcefield_name
        )

        # 写出
        output_file = self.workdir / "system.data"
        comments = {
            'pair_coeffs': 'lj/cut/coul/long',
            'bond_coeffs': 'harmonic',
            'angle_coeffs': 'harmonic',
            'atoms': 'full'
        }
        write_lammps_data(output_file, system, comments)

        print(f"  ✓ LAMMPS data: {output_file}\n")
        return output_file

    def _replicate_topology(
        self,
        single: Dict[str, List[str]],
        coords: np.ndarray,
        n_mol: int,
        n_atoms_per_mol: int
    ) -> Dict[str, List[str]]:
        """复制单分子拓扑到多分子体系"""

        # 复制原子
        new_atoms = []
        for mol_id in range(n_mol):
            for atom_line in single['atoms']:
                parts = atom_line.split()
                atom_id = int(parts[0]) + mol_id * n_atoms_per_mol
                mol_id_new = mol_id + 1
                atom_type = parts[2]
                charge = parts[3]

                # 使用Packmol坐标
                coord_idx = atom_id - 1
                x, y, z = coords[coord_idx]

                new_atoms.append(
                    f"{atom_id} {mol_id_new} {atom_type} {charge} "
                    f"{x:.8f} {y:.8f} {z:.8f}"
                )

        # 复制bonds/angles/dihedrals
        def replicate_section(section_data: List[str], n_entries_per_mol: int):
            new_section = []
            for mol_id in range(n_mol):
                offset = mol_id * n_atoms_per_mol
                entry_offset = mol_id * n_entries_per_mol

                for entry in section_data:
                    parts = entry.split()
                    new_id = int(parts[0]) + entry_offset
                    entry_type = parts[1]
                    atom_ids = [str(int(x) + offset) for x in parts[2:]]

                    new_section.append(f"{new_id} {entry_type} " + " ".join(atom_ids))

            return new_section

        return {
            'header': [],  # 稍后生成
            'masses': single['masses'],
            'pair_coeffs': single['pair_coeffs'],
            'bond_coeffs': single['bond_coeffs'],
            'angle_coeffs': single['angle_coeffs'],
            'dihedral_coeffs': single['dihedral_coeffs'],
            'improper_coeffs': single['improper_coeffs'],
            'atoms': new_atoms,
            'bonds': replicate_section(single['bonds'], len(single['bonds'])),
            'angles': replicate_section(single['angles'], len(single['angles'])),
            'dihedrals': replicate_section(single['dihedrals'], len(single['dihedrals'])),
            'impropers': replicate_section(single['impropers'], len(single['impropers'])),
        }

    def _make_system_header(
        self,
        single: Dict[str, List[str]],
        n_mol: int,
        box_lengths: tuple[float, float, float],
        forcefield_name: str
    ) -> List[str]:
        """生成多分子体系的header"""

        n_atoms = len(single['atoms']) * n_mol
        n_bonds = len(single['bonds']) * n_mol
        n_angles = len(single['angles']) * n_mol
        n_dihedrals = len(single['dihedrals']) * n_mol
        n_impropers = len(single['impropers']) * n_mol

        n_atom_types = len(single['masses'])
        n_bond_types = len(single['bond_coeffs']) if single['bond_coeffs'] else 0
        n_angle_types = len(single['angle_coeffs']) if single['angle_coeffs'] else 0
        n_dihedral_types = len(single['dihedral_coeffs']) if single['dihedral_coeffs'] else 0
        n_improper_types = len(single['improper_coeffs']) if single['improper_coeffs'] else 0

        title = f"LAMMPS data file - {forcefield_name} force field"

        return make_lammps_header(
            n_atoms, n_bonds, n_angles, n_dihedrals, n_impropers,
            n_atom_types, n_bond_types, n_angle_types, n_dihedral_types, n_improper_types,
            box_lengths,
            title
        )
