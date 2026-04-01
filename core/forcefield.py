"""
力场抽象基类
"""
from abc import ABC, abstractmethod
from pathlib import Path
from rdkit import Chem
from typing import Dict, Any


class ForceField(ABC):
    """力场抽象基类"""

    def __init__(self, workdir: Path):
        """
        Args:
            workdir: 工作目录
        """
        self.workdir = Path(workdir)
        self.workdir.mkdir(parents=True, exist_ok=True)

    @abstractmethod
    def generate_parameters(
        self,
        mol: Chem.Mol,
        smiles: str,
        name: str
    ) -> Path:
        """
        生成力场参数

        Args:
            mol: RDKit分子对象
            smiles: SMILES字符串
            name: 分子名称

        Returns:
            单分子LAMMPS data文件路径
        """
        pass

    @abstractmethod
    def get_lammps_settings(self) -> Dict[str, Any]:
        """
        获取LAMMPS模拟设置

        Returns:
            包含pair_style, special_bonds等设置的字典
        """
        pass

    @property
    @abstractmethod
    def name(self) -> str:
        """力场名称"""
        pass

    @property
    @abstractmethod
    def citation(self) -> str:
        """力场引用文献"""
        pass
