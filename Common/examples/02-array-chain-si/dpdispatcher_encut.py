#!/usr/bin/env python3
"""ENCUT scan of Si submitted with dpdispatcher (pip install dpdispatcher).

dpdispatcher writes the Slurm/PBS/LSF script, submits it, polls the queue
and (for SSH contexts) copies files back. Here everything is local: run
this script on the cluster login node inside tmux (it waits until all
jobs are done). --dry-run only prepares and stages the task directories
(dp_encut_*/ and dp_work/<hash>/) without submitting anything.

  python dpdispatcher_encut.py --dry-run     # stage files, do not submit
  python dpdispatcher_encut.py               # submit, poll, collect

Requires: POSCAR, KPOINTS, INCAR.encut.template and POTCAR in this dir.
"""
import argparse
import os
import shutil

from dpdispatcher import Machine, Resources, Submission, Task

ENCUTS = [250, 300, 350, 400, 450, 500, 550]


def prepare():
    tasks = []
    template = open("INCAR.encut.template").read()
    for e in ENCUTS:
        d = "dp_encut_%d" % e
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "INCAR"), "w", newline="\n") as f:
            f.write(template.replace("@ENCUT@", str(e)))
        for fn in ("POSCAR", "KPOINTS", "POTCAR"):
            shutil.copy(fn, d)
        tasks.append(Task(command="mpirun -np 8 vasp_std",
                          task_work_path=d + "/",
                          forward_files=["INCAR", "POSCAR", "KPOINTS",
                                         "POTCAR"],
                          backward_files=["OUTCAR", "OSZICAR"],
                          outlog="vasp.out", errlog="vasp.err"))
    return tasks


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--queue", default="cpu")
    args = p.parse_args()

    machine = Machine.load_from_dict({
        "batch_type": "Slurm",          # also: PBS, Torque, LSF, Shell
        "context_type": "LocalContext",  # SSHContext for a remote cluster
        "local_root": "./",
        "remote_root": "./dp_work",
        # copy files instead of symlinking them (symlinks need extra
        # rights on Windows and break if the link targets move)
        "remote_profile": {"symlink": False},
    })
    resources = Resources.load_from_dict({
        "number_node": 1,
        "cpu_per_node": 8,
        "gpu_per_node": 0,
        "queue_name": args.queue,
        "group_size": 1,                 # tasks per job script
        "custom_flags": ["#SBATCH --time=00:30:00"],
        "module_list": [],               # e.g. ["vasp/6.5.1"]
        "envs": {"OMP_NUM_THREADS": 1},
    })
    submission = Submission(work_base="./", machine=machine,
                            resources=resources, task_list=prepare())
    submission.run_submission(dry_run=args.dry_run)
    print("done; results in dp_encut_*/ (OUTCAR copied back)")


if __name__ == "__main__":
    main()
