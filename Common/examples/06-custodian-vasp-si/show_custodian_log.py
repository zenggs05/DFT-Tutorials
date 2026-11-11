#!/usr/bin/env python
"""Print a readable summary of custodian.json (which errors were fixed and how).

Usage:
    python show_custodian_log.py [path/to/custodian.json]

Only the standard library is used, so it also works on a laptop that does not
have custodian installed.
"""
import json
import sys
from pathlib import Path


def short(obj):
    """Return the class name of an MSON dict such as {"@class": "VaspErrorHandler"}."""
    if isinstance(obj, dict):
        return obj.get("@class", "?")
    return str(obj)


def main() -> int:
    path = Path(sys.argv[1] if len(sys.argv) > 1 else "custodian.json")
    if not path.is_file():
        print(f"{path} not found (custodian writes it in the run directory)")
        return 1
    runs = json.loads(path.read_text())
    for i, run in enumerate(runs, 1):
        job = run.get("job", {})
        suffix = job.get("suffix", "") if isinstance(job, dict) else ""
        print(f"=== job {i}: {short(job)} {suffix}")
        corrections = run.get("corrections", [])
        if not corrections:
            print("    no corrections needed")
        for k, corr in enumerate(corrections, 1):
            errors = ", ".join(map(str, corr.get("errors", [])))
            print(f"  [{k}] handler = {short(corr.get('handler'))}")
            print(f"      errors  = {errors}")
            for act in corr.get("actions") or []:
                print(f"      action  = {json.dumps(act)}")
            if not corr.get("actions"):
                print("      action  = None (handler could not fix it)")
        for key in ("max_errors", "max_errors_per_job", "nonzero_return_code"):
            if run.get(key):
                print(f"  !! {key} triggered")
        if run.get("validator"):
            print(f"  !! validation failed: {short(run['validator'])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
