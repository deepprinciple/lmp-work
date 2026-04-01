"""
文件I/O工具函数
"""
from pathlib import Path
import numpy as np
from typing import Dict, List, Tuple


def parse_lammps_data(data_file: Path) -> Dict[str, List[str]]:
    """
    解析LAMMPS data文件

    Returns:
        dict: 包含各section的字典
    """
    lines = data_file.read_text().splitlines()

    data = {
        'header': [],
        'masses': [],
        'pair_coeffs': [],
        'bond_coeffs': [],
        'angle_coeffs': [],
        'dihedral_coeffs': [],
        'improper_coeffs': [],
        'atoms': [],
        'bonds': [],
        'angles': [],
        'dihedrals': [],
        'impropers': []
    }

    section = 'header'
    for line in lines:
        stripped = line.strip()

        # 识别section
        if 'Masses' in line:
            section = 'masses'
            continue
        elif 'Pair Coeffs' in line:
            section = 'pair_coeffs'
            continue
        elif 'Bond Coeffs' in line:
            section = 'bond_coeffs'
            continue
        elif 'Angle Coeffs' in line:
            section = 'angle_coeffs'
            continue
        elif 'Dihedral Coeffs' in line:
            section = 'dihedral_coeffs'
            continue
        elif 'Improper Coeffs' in line:
            section = 'improper_coeffs'
            continue
        elif 'Atoms' in line:
            section = 'atoms'
            continue
        elif 'Bonds' in line:
            section = 'bonds'
            continue
        elif 'Angles' in line:
            section = 'angles'
            continue
        elif 'Dihedrals' in line:
            section = 'dihedrals'
            continue
        elif 'Impropers' in line:
            section = 'impropers'
            continue

        if stripped and not stripped.startswith('#'):
            data[section].append(line)

    return data


def write_lammps_data(
    output_file: Path,
    data: Dict[str, List[str]],
    comments: Dict[str, str] = None
) -> None:
    """
    写出LAMMPS data文件

    Args:
        output_file: 输出文件路径
        data: 包含各section的字典
        comments: 可选的section注释
    """
    comments = comments or {}

    with open(output_file, 'w') as f:
        # Header
        for line in data['header']:
            f.write(line + "\n")

        # Coefficients sections
        sections = [
            ('masses', 'Masses'),
            ('pair_coeffs', 'Pair Coeffs'),
            ('bond_coeffs', 'Bond Coeffs'),
            ('angle_coeffs', 'Angle Coeffs'),
            ('dihedral_coeffs', 'Dihedral Coeffs'),
            ('improper_coeffs', 'Improper Coeffs'),
        ]

        for key, title in sections:
            if data[key]:
                comment = comments.get(key, '')
                f.write(f"\n{title}")
                if comment:
                    f.write(f" # {comment}")
                f.write("\n\n")
                for line in data[key]:
                    f.write(line + "\n")

        # Topology sections
        topology_sections = [
            ('atoms', 'Atoms'),
            ('bonds', 'Bonds'),
            ('angles', 'Angles'),
            ('dihedrals', 'Dihedrals'),
            ('impropers', 'Impropers'),
        ]

        for key, title in topology_sections:
            if data[key]:
                comment = comments.get(key, '')
                f.write(f"\n{title}")
                if comment:
                    f.write(f" # {comment}")
                f.write("\n\n")
                for line in data[key]:
                    f.write(line + "\n")


def read_xyz(xyz_file: Path) -> Tuple[int, np.ndarray, List[str]]:
    """
    读取XYZ文件

    Returns:
        (n_atoms, coordinates, elements)
    """
    with open(xyz_file) as f:
        n_atoms = int(f.readline().strip())
        f.readline()  # comment line

        coords = []
        elements = []
        for line in f:
            parts = line.split()
            if len(parts) >= 4:
                elements.append(parts[0])
                coords.append([float(x) for x in parts[1:4]])

    return n_atoms, np.array(coords), elements


def make_lammps_header(
    n_atoms: int,
    n_bonds: int,
    n_angles: int,
    n_dihedrals: int,
    n_impropers: int,
    n_atom_types: int,
    n_bond_types: int,
    n_angle_types: int,
    n_dihedral_types: int,
    n_improper_types: int,
    box_bounds: Tuple[float, float, float],
    title: str = "LAMMPS data file"
) -> List[str]:
    """
    生成LAMMPS data文件的header

    Args:
        box_bounds: (Lx, Ly, Lz) 盒子边长
    """
    Lx, Ly, Lz = box_bounds

    header = [title, ""]

    # Counts
    header.append(f"{n_atoms} atoms")
    if n_bonds > 0:
        header.append(f"{n_bonds} bonds")
    if n_angles > 0:
        header.append(f"{n_angles} angles")
    if n_dihedrals > 0:
        header.append(f"{n_dihedrals} dihedrals")
    if n_impropers > 0:
        header.append(f"{n_impropers} impropers")

    header.append("")

    # Types
    header.append(f"{n_atom_types} atom types")
    if n_bond_types > 0:
        header.append(f"{n_bond_types} bond types")
    if n_angle_types > 0:
        header.append(f"{n_angle_types} angle types")
    if n_dihedral_types > 0:
        header.append(f"{n_dihedral_types} dihedral types")
    if n_improper_types > 0:
        header.append(f"{n_improper_types} improper types")

    header.append("")

    # Box
    header.extend([
        f"0.0 {Lx:.8f} xlo xhi",
        f"0.0 {Ly:.8f} ylo yhi",
        f"0.0 {Lz:.8f} zlo zhi",
        ""
    ])

    return header
