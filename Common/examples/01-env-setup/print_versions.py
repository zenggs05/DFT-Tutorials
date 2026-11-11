#!/usr/bin/env python3
"""Print Python and package versions of the DFT environment.

  python print_versions.py            # table on screen
  python print_versions.py > versions.txt   # keep it with your results

A missing package is reported as "not installed" instead of crashing.
For the full check (environment variables, ulimit, ~/.vaspkit ...) use
Common/scripts/env_check.sh.
"""
import platform
import sys
from importlib import metadata

PACKAGES = ["numpy", "scipy", "matplotlib", "ase", "pymatgen",
            "pymatgen-core", "spglib", "seekpath", "phonopy", "mp-api",
            "sumo", "pyprocar", "py4vasp", "cclib", "goodvibes",
            "lobsterpy", "dpdispatcher", "rdkit"]

# minimum Python required by recent releases (2026-11)
NEEDS = {"pymatgen": (3, 11), "phonopy": (3, 11), "mp-api": (3, 12),
         "lobsterpy": (3, 11)}


def main():
    print("Python   %s (%s)" % (platform.python_version(), sys.executable))
    print("Platform %s" % platform.platform())
    print("-" * 44)
    for name in PACKAGES:
        try:
            ver = metadata.version(name)
        except metadata.PackageNotFoundError:
            ver = "not installed"
        print("%-16s %s" % (name, ver))
    low = [n for n, v in NEEDS.items() if sys.version_info[:2] < v]
    if low:
        print("-" * 44)
        print("WARNING: Python %d.%d is too old for the latest %s"
              % (sys.version_info[0], sys.version_info[1], ", ".join(low)))


if __name__ == "__main__":
    main()
