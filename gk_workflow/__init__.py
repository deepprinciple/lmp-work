"""Compatibility package for the Green-Kubo workflow project.

This package maps the historical ``gk_workflow.*`` import paths onto the
current repository layout, where subpackages such as ``analysis`` and
``core`` live at the project root.
"""

from importlib import import_module
from pathlib import Path

_PACKAGE_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _PACKAGE_DIR.parent

# Let ``import gk_workflow.analysis`` resolve to ``analysis/`` at the project root.
__path__ = [str(_PACKAGE_DIR), str(_PROJECT_ROOT)]

__version__ = "2.0.0"
__all__ = [
    "MoleculeStructure",
    "PackmolBuilder",
    "LAMMPSBuilder",
    "OpenFF",
]


def __getattr__(name: str):
    if name in {"MoleculeStructure", "PackmolBuilder", "LAMMPSBuilder"}:
        module = import_module(".core", __name__)
        return getattr(module, name)
    if name == "OpenFF":
        module = import_module(".forcefields", __name__)
        return getattr(module, name)
    raise AttributeError(name)
