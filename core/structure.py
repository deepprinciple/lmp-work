"""
分子3D结构生成（RDKit）
"""
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors


class MoleculeStructure:
    """分子3D结构生成器"""

    def __init__(self, smiles: str, name: str = "molecule"):
        """
        Args:
            smiles: SMILES字符串
            name: 分子名称
        """
        self.smiles = smiles
        self.name = name
        self.mol = None

    def generate(self, optimize: bool = True, max_iters: int = 500) -> Chem.Mol:
        """
        Generate a 3D structure from SMILES.

        Args:
            optimize: Run MMFF optimization after embedding.
            max_iters: Maximum MMFF optimization iterations.

        Returns:
            RDKit molecule with 3D coordinates.
        """
        mol = Chem.MolFromSmiles(self.smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {self.smiles}")

        mol = Chem.AddHs(mol)

        if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
            raise RuntimeError("Failed to generate a 3D conformer")

        if optimize:
            if not AllChem.MMFFHasAllMoleculeParams(mol):
                raise RuntimeError("MMFF parameters are not available for this molecule")
            status = AllChem.MMFFOptimizeMolecule(mol, maxIters=max_iters)
            if status < 0:
                raise RuntimeError("MMFF optimization failed")

        self.mol = mol
        return self.mol

    def _print_info(self):
        """打印分子信息"""
        print("\n  分子信息:")
        print(f"    SMILES: {Chem.MolToSmiles(Chem.RemoveHs(self.mol))}")
        print(f"    分子式: {Chem.rdMolDescriptors.CalcMolFormula(self.mol)}")
        print(f"    原子数: {self.mol.GetNumAtoms()}")
        print(f"    键数: {self.mol.GetNumBonds()}")
        print(f"    分子量: {Descriptors.MolWt(self.mol):.2f} g/mol")

    def save_xyz(self, output_path: Path) -> Path:
        """
        保存为XYZ文件

        Args:
            output_path: 输出路径

        Returns:
            输出文件路径
        """
        if self.mol is None:
            raise RuntimeError("请先调用 generate()")

        Chem.MolToXYZFile(self.mol, str(output_path))
        print(f"    保存至: {output_path}\n")
        return output_path

    def save_pdb(self, output_path: Path) -> Path:
        """
        保存为PDB文件

        Args:
            output_path: 输出路径

        Returns:
            输出文件路径
        """
        if self.mol is None:
            raise RuntimeError("请先调用 generate()")

        Chem.MolToPDBFile(self.mol, str(output_path))
        return output_path

    def save_mol(self, output_path: Path) -> Path:
        """
        保存为MDL mol文件

        Args:
            output_path: 输出路径

        Returns:
            输出文件路径
        """
        if self.mol is None:
            raise RuntimeError("请先调用 generate()")

        Chem.MolToMolFile(self.mol, str(output_path), includeStereo=False)
        return output_path

    def save_mol2(self, output_path: Path) -> Path:
        """
        保存为MOL2文件（用于LigParGen）

        Args:
            output_path: 输出路径

        Returns:
            输出文件路径
        """
        if self.mol is None:
            raise RuntimeError("请先调用 generate()")

        return self.save_mol(output_path)

    @property
    def n_atoms(self) -> int:
        """原子数"""
        return self.mol.GetNumAtoms() if self.mol else 0

    @property
    def molecular_weight(self) -> float:
        """分子量 (g/mol)"""
        return Descriptors.MolWt(self.mol) if self.mol else 0.0
