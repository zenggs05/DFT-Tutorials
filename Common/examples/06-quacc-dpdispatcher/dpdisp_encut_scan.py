#!/usr/bin/env python
"""Batch-submit a VASP ENCUT scan with dpdispatcher and collect the energies.

Directory layout created by this script (inside work/):
    work/encut_300/{INCAR,KPOINTS,POSCAR,POTCAR}  ...  work/encut_600/...
dpdispatcher uploads every task folder, writes one Slurm script per job
(group_size tasks per job), submits, polls until done and downloads OUTCAR/OSZICAR.

Usage:
    python dpdisp_encut_scan.py --template ../06-custodian-vasp-si --dry-run
    python dpdisp_encut_scan.py --template ../06-custodian-vasp-si
    python dpdisp_encut_scan.py --local          # Shell + LocalContext test (Linux/WSL)

The template directory must contain INCAR, KPOINTS, POSCAR and POTCAR.
"""
import argparse
import re
import shutil
from pathlib import Path

from dpdispatcher import Machine, Resources, Submission, Task

ENCUTS = [300, 350, 400, 450, 500, 550, 600]


def prepare(template: Path, work: Path) -> list:
    work.mkdir(exist_ok=True)
    names = []
    for enc in ENCUTS:
        d = work / f"encut_{enc}"
        d.mkdir(exist_ok=True)
        for f in ("KPOINTS", "POSCAR", "POTCAR"):
            src = template / f
            if src.exists():
                shutil.copy(src, d / f)
        incar = (template / "INCAR").read_text()
        incar = re.sub(r"(?im)^\s*(ENCUT|NSW|IBRION|ISIF)\s*=.*$", "", incar)
        incar += f"\nENCUT = {enc}\nNSW = 0\nIBRION = -1\n"
        (d / "INCAR").write_text(incar)
        names.append(d.name)
    return names


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--template", default=".", help="folder with INCAR KPOINTS POSCAR POTCAR")
    p.add_argument("--dry-run", action="store_true", help="upload and write scripts only")
    p.add_argument("--local", action="store_true", help="use Shell + LocalContext (test)")
    args = p.parse_args()

    work = Path("work")
    tasks_dirs = prepare(Path(args.template), work)

    if args.local:
        machine = Machine.load_from_dict({"batch_type": "Shell", "context_type": "LocalContext",
                                          "local_root": "./", "remote_root": "./remote_tmp"})
        resources = Resources.load_from_dict({"number_node": 1, "cpu_per_node": 1,
                                              "group_size": len(tasks_dirs)})
        # fake "VASP" that writes the three files we download (Linux/WSL/macOS only)
        command = ("echo fake > OUTCAR; echo done > vasp.out; "
                   "echo '   1 F= -.10000000E+02 E0= -.10000000E+02  d E =0.0' > OSZICAR")
    else:
        machine = Machine.load_from_json("dpdispatcher_machine.json")
        resources = Resources.load_from_json("dpdispatcher_resources.json")
        command = "srun vasp_std > vasp.out 2>&1"

    inputs = ["INCAR", "KPOINTS", "POSCAR", "POTCAR"]
    if not args.local and not (work / tasks_dirs[0] / "POTCAR").exists():
        raise SystemExit("POTCAR missing in the template folder (it is never shipped in the repo)")
    task_list = [
        Task(command=command, task_work_path=f"{name}/",
             forward_files=[f for f in inputs if (work / name / f).exists()],
             backward_files=["OUTCAR", "OSZICAR", "vasp.out"])
        for name in tasks_dirs
    ]
    sub = Submission(work_base=str(work), machine=machine, resources=resources,
                     task_list=task_list)
    sub.run_submission(dry_run=args.dry_run, check_interval=60)

    if args.dry_run:
        print("dry run: inputs and job scripts uploaded, nothing submitted")
        return
    print(f"{'ENCUT':>6s} {'E0 (eV)':>16s}")
    for name in tasks_dirs:
        osz = work / name / "OSZICAR"
        energy = "n/a"
        if osz.exists():
            m = re.findall(r"E0=\s*([-+0-9.Ee]+)", osz.read_text())
            energy = m[-1] if m else "n/a"
        print(f"{name.split('_')[1]:>6s} {energy:>16s}")


if __name__ == "__main__":
    main()
