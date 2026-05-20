"""
物理常数和单位转换
"""

# 基本物理常数
BOLTZMANN_J_K = 1.380649e-23      # J/K
AVOGADRO = 6.02214076e23          # 1/mol
ELEMENTARY_CHARGE = 1.602176634e-19  # C
VACUUM_PERMITTIVITY = 8.8541878128e-12  # F/m

# Atomic data for BAMBOO data-file generation.  Atom types are assigned in
# ascending atomic-number order so the data file and pair_coeff elements stay
# synchronized.
ATOMIC_MASS: dict[str, float] = {
    "H":  1.008,
    "Li": 6.941,
    "B":  10.811,
    "C":  12.011,
    "N":  14.007,
    "O":  15.9994,
    "F":  18.998,
    "Na": 22.990,
    "P":  30.974,
    "S":  32.06,
    "Cl": 35.45,
    "K":  39.098,
}

ATOMIC_NUMBER: dict[str, int] = {
    "H": 1,
    "Li": 3,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Na": 11,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "K": 19,
}

# 单位转换
ANGSTROM_TO_METER = 1e-10
ATM_TO_PASCAL = 101325.0
KCAL_MOL_TO_JOULE = 4184.0
FEMTOSECOND_TO_SECOND = 1e-15

# LAMMPS单位系统常数（real units）
class LAMMPSUnits:
    """LAMMPS real单位系统"""
    ENERGY = "kcal/mol"
    DISTANCE = "Angstrom"
    TIME = "femtosecond"
    MASS = "g/mol"
    PRESSURE = "atmosphere"
    TEMPERATURE = "Kelvin"
    CHARGE = "e"  # 元电荷

    # 热流单位转换: kcal/mol/Å^2/fs -> W/m^2
    @staticmethod
    def heat_flux_to_SI():
        return KCAL_MOL_TO_JOULE / (AVOGADRO * ANGSTROM_TO_METER**2 * FEMTOSECOND_TO_SECOND)

# 默认参数
class Defaults:
    """默认模拟参数"""
    TEMPERATURE = 298.0  # K
    PRESSURE = 1.0       # atm
    TIMESTEP = 0.5       # fs
    CUTOFF = 12.0        # Å
    PPPM_ACCURACY = 1e-5
    NEIGHBOR_SKIN = 2.0  # Å
