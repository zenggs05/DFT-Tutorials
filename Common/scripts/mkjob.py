#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mkjob.py -- generate a batch job script for VASP / Gaussian 16 / CP2K / OpenMX
            on Slurm, PBS (PBS Pro/OpenPBS, optional Torque syntax) or LSF.
            按“软件 + 调度器 + 核数”生成作业脚本：资源请求、module、MPI x OpenMP、
            ulimit、CRLF(dos2unix) 检查、可选 scratch 目录和运行后的成功检查。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

What the generated script contains
----------------------------------
* scheduler header (#SBATCH / #PBS / #BSUB) built from the options below;
* a clearly marked "USER SETTINGS" block (executable, input, output, scratch);
* module purge/load lines (or placeholders) and an optional "source" line;
* ulimit -s unlimited, OMP_NUM_THREADS (+ OMP_PLACES/OMP_PROC_BIND for hybrid
  MPI x OpenMP runs) and OMP_STACKSIZE;
* pre-run checks: input files exist, no Windows CRLF line endings
  (fix: dos2unix), program specific checks (VASP: INCAR/POSCAR/POTCAR/KPOINTS,
  number of POTCAR datasets vs. species; Gaussian: %NProcShared/%CPU and %Mem
  vs. the requested resources; CP2K: input file; OpenMX: DATA.PATH directory);
* optional per-job scratch directory (Gaussian: GAUSS_SCRDIR; other codes:
  run in scratch, copy results back on exit);
* post-run check: VASP "General timing and accounting" + "reached required
  accuracy" for relaxations, Gaussian "Normal termination" vs.
  "Error termination", CP2K "PROGRAM ENDED AT" + "SCF run NOT converged",
  OpenMX "Elapsed.Time." in System.Name.out + "scf_conv=0" in System.Name.md.

Usage examples
--------------
    python mkjob.py vasp -N 2 -n 64 -t 24:00:00 -q cpu --module "vasp/6.5.1"
    python mkjob.py vasp -s pbs -N 1 -n 64 --mem 200G -o job.pbs
    python mkjob.py gaussian -i benzene.gjf -c 16 --mem 40G --scratch '${TMPDIR:-/tmp}'
    python mkjob.py cp2k -i h2o_md.inp -N 2 -n 32 -c 2 --source ~/cp2k/tools/toolchain/install/setup
    python mkjob.py openmx -i Si.dat -n 32 -c 2 -s lsf -q normal
    python mkjob.py cp2k -i opt.inp --launcher mpirun -o -     # print to stdout

Submit afterwards with sbatch job.sh / qsub job.sh / bsub < job.sh.
Partition/queue names, accounts, module names and the MPI launcher are site
specific: always compare with the user guide of your cluster.

Dependencies: Python >= 3.9, standard library only (argparse, string.Template).
"""

import argparse
import datetime
import os
import re
import shlex
import stat
import sys
from string import Template

PROGRAMS = ("vasp", "gaussian", "cp2k", "openmx")
ALIASES = {"g16": "gaussian", "gau": "gaussian", "omx": "openmx", "cp2k.psmp": "cp2k"}
SCHEDULERS = ("slurm", "pbs", "lsf")

# Program defaults: executable, input, MPI ranks per node, threads per rank
DEFAULTS = {
    "vasp": {"exe": "vasp_std", "input": "", "ntasks": 32, "threads": 1},
    "gaussian": {"exe": "g16", "input": "input.gjf", "ntasks": 1, "threads": 16},
    "cp2k": {"exe": "cp2k.psmp", "input": "input.inp", "ntasks": 32, "threads": 1},
    "openmx": {"exe": "openmx", "input": "input.dat", "ntasks": 32, "threads": 1},
}


class JobTemplate(Template):
    """string.Template with '@@' as delimiter, so bash '$' needs no escaping."""
    delimiter = "@@"


def die(msg):
    sys.exit("mkjob.py: error: " + msg)


# --------------------------------------------------------------------------
# argument helpers
# --------------------------------------------------------------------------
def parse_time(text):
    """Return wall time in seconds. Accepts D-HH:MM:SS, HH:MM:SS, HH:MM or minutes."""
    t = text.strip()
    m = re.fullmatch(r"(?:(\d+)-)?(\d+):(\d{1,2})(?::(\d{1,2}))?", t)
    if m:
        days = int(m.group(1) or 0)
        h, mi = int(m.group(2)), int(m.group(3))
        s = int(m.group(4) or 0)
        if mi >= 60 or s >= 60:
            raise argparse.ArgumentTypeError("minutes/seconds must be < 60: %r" % text)
        return ((days * 24 + h) * 60 + mi) * 60 + s
    if re.fullmatch(r"\d+", t):
        return int(t) * 60
    raise argparse.ArgumentTypeError(
        "cannot parse wall time %r (use D-HH:MM:SS, HH:MM:SS, HH:MM or minutes)" % text)


def fmt_hms(sec):
    h, rem = divmod(sec, 3600)
    return "%d:%02d:%02d" % (h, rem // 60, rem % 60)


def fmt_slurm_time(sec):
    d, rem = divmod(sec, 86400)
    h, rem = divmod(rem, 3600)
    core = "%02d:%02d:%02d" % (h, rem // 60, rem % 60)
    return ("%d-%s" % (d, core)) if d else core


def fmt_lsf_time(sec):
    """LSF -W takes [hour:]minute; seconds are rounded up to a full minute."""
    minutes = -(-sec // 60)
    return "%d:%02d" % divmod(minutes, 60)


def parse_mem(text):
    """Return memory in MiB from strings like 40G, 40GB, 4000M, 0.5T."""
    m = re.fullmatch(r"\s*(\d+(?:\.\d+)?)\s*([kKmMgGtT]?)[bB]?\s*", text)
    if not m:
        raise argparse.ArgumentTypeError("cannot parse memory %r (examples: 40G, 4000M)" % text)
    val = float(m.group(1))
    unit = m.group(2).upper() or "M"
    factor = {"K": 1.0 / 1024, "M": 1.0, "G": 1024.0, "T": 1024.0 * 1024}[unit]
    mib = int(round(val * factor))
    if mib <= 0:
        raise argparse.ArgumentTypeError("memory must be positive: %r" % text)
    return mib


def fmt_mem_slurm(mib):
    return ("%dG" % (mib // 1024)) if mib % 1024 == 0 else ("%dM" % mib)


def fmt_mem_pbs(mib):
    return ("%dgb" % (mib // 1024)) if mib % 1024 == 0 else ("%dmb" % mib)


def clean_name(name):
    s = re.sub(r"[^A-Za-z0-9_.-]", "_", name)
    if not s or not s[0].isalpha():
        s = "j" + s          # PBS requires the name to start with a letter
    return s[:40]


# --------------------------------------------------------------------------
# scheduler headers
# --------------------------------------------------------------------------
def header_slurm(a):
    L = ["#SBATCH --job-name=%s" % a.name]
    if a.queue:
        L.append("#SBATCH --partition=%s" % a.queue)
    else:
        L.append("##SBATCH --partition=cpu              # sinfo -s lists the partitions")
    L += ["#SBATCH --nodes=%d" % a.nodes,
          "#SBATCH --ntasks-per-node=%d          # MPI ranks per node" % a.ntasks,
          "#SBATCH --cpus-per-task=%d            # OpenMP threads per rank" % a.threads]
    if a.mem_mib:
        L.append("#SBATCH --mem=%s                    # memory PER NODE" % fmt_mem_slurm(a.mem_mib))
    else:
        L.append("##SBATCH --mem=0                      # 0 = all memory of the node (site dependent)")
    L.append("#SBATCH --time=%s" % fmt_slurm_time(a.time_s))
    if a.gpus:
        L.append("#SBATCH --gres=gpu:%d                 # GPUs per node (syntax may differ per site)" % a.gpus)
    L += ["#SBATCH --output=%x-%j.out", "#SBATCH --error=%x-%j.err"]
    L.append(("#SBATCH --account=%s" % a.account) if a.account else "##SBATCH --account=my_project")
    if a.qos:
        L.append("#SBATCH --qos=%s" % a.qos)
    L.append("#SBATCH --exclusive" if a.exclusive else "##SBATCH --exclusive")
    return L


def header_pbs(a):
    L = ["#PBS -N %s" % a.name]
    L.append(("#PBS -q %s" % a.queue) if a.queue else "##PBS -q workq                       # qstat -Q lists the queues")
    if a.torque:
        L.append("#PBS -l nodes=%d:ppn=%d" % (a.nodes, a.ntasks * a.threads))
        if a.mem_mib:
            L.append("#PBS -l mem=%s                       # total memory of the job (Torque)"
                     % fmt_mem_pbs(a.mem_mib * a.nodes))
        if a.gpus:
            L.append("# GPUs: Torque usually needs nodes=N:ppn=P:gpus=G (site dependent)")
    else:
        chunk = "select=%d:ncpus=%d:mpiprocs=%d:ompthreads=%d" % (
            a.nodes, a.ntasks * a.threads, a.ntasks, a.threads)
        if a.mem_mib:
            chunk += ":mem=%s" % fmt_mem_pbs(a.mem_mib)
        if a.gpus:
            chunk += ":ngpus=%d" % a.gpus
        L.append("#PBS -l %s" % chunk)
        L.append("#PBS -l place=scatter:excl" if a.exclusive else "##PBS -l place=scatter:excl          # one chunk per node, exclusive")
    L += ["#PBS -l walltime=%s" % fmt_hms(a.time_s),
          "#PBS -j oe                            # merge stdout and stderr",
          "#PBS -o %s.pbslog" % a.name]
    L.append(("#PBS -A %s" % a.account) if a.account else "##PBS -A my_project")
    return L


def header_lsf(a):
    total = a.nodes * a.ntasks
    L = ["#BSUB -J %s" % a.name]
    L.append(("#BSUB -q %s" % a.queue) if a.queue else "##BSUB -q normal                      # bqueues lists the queues")
    L += ["#BSUB -n %d                           # TOTAL number of MPI slots" % total,
          '#BSUB -R "span[ptile=%d]"            # slots per host' % a.ntasks]
    if a.threads > 1:
        L.append('##BSUB -R "affinity[core(%d)]"      # cores per rank (only if the site enables affinity)'
                 % a.threads)
    L.append("#BSUB -W %s                       # wall time [hour:]minute" % fmt_lsf_time(a.time_s))
    if a.mem_mib:
        per_slot = max(1, a.mem_mib // max(1, a.ntasks))
        L.append('##BSUB -R "rusage[mem=%d]"         # per slot in MB IF LSF_UNIT_FOR_LIMITS=MB (check site)'
                 % per_slot)
    if a.gpus:
        L.append('#BSUB -gpu "num=%d"                 # GPUs per host (LSF 10.1 syntax)' % a.gpus)
    L += ["#BSUB -o %J.out", "#BSUB -e %J.err"]
    L.append(("#BSUB -P %s" % a.account) if a.account else "##BSUB -P my_project")
    L.append("#BSUB -x" if a.exclusive else "##BSUB -x                             # exclusive hosts")
    return L


# scheduler specific shell variables
SCHED_VARS = {
    "slurm": {"workdir": '${SLURM_SUBMIT_DIR:-$PWD}', "jobid": '${SLURM_JOB_ID:-$$}',
              "submit": "sbatch", "nodes": '${SLURM_JOB_NODELIST:-unknown}'},
    "pbs": {"workdir": '${PBS_O_WORKDIR:-$PWD}', "jobid": '${PBS_JOBID:-$$}',
            "submit": "qsub", "nodes": '$(sort -u "${PBS_NODEFILE:-/dev/null}" 2>/dev/null | tr "\\n" " ")'},
    "lsf": {"workdir": '${LS_SUBCWD:-$PWD}', "jobid": '${LSB_JOBID:-$$}',
            "submit": "bsub <", "nodes": '${LSB_HOSTS:-unknown}'},
}


# --------------------------------------------------------------------------
# launcher line
# --------------------------------------------------------------------------
def mpi_prefix(a):
    """Return the MPI launch prefix (without the program)."""
    if a.scheduler == "slurm":
        if a.launcher == "srun":
            return ('srun --ntasks="${SLURM_NTASKS:-@@{np}}" --cpus-per-task="${OMP_NUM_THREADS}" '
                    "--cpu-bind=cores")
        return 'mpirun -np "${SLURM_NTASKS:-@@{np}}"'
    if a.scheduler == "pbs":
        if a.launcher == "srun":
            die("--launcher srun only works with Slurm")
        return 'mpirun -np "${NP}" -machinefile "${PBS_NODEFILE}"'
    if a.launcher == "srun":
        die("--launcher srun only works with Slurm")
    return 'mpirun -np "${NP}"'


# --------------------------------------------------------------------------
# script body templates (bash); placeholders are @@{name}
# --------------------------------------------------------------------------
COMMON_TOP = r"""
set -o pipefail

# ======================= USER SETTINGS (edit here) ==========================
EXE="@@{exe}"                     # executable (full path if not on PATH)
INPUT="@@{input}"                 # main input file (VASP: unused)
OUTPUT="@@{output}"               # file that receives stdout/stderr
NP=@@{np}                         # total MPI ranks = nodes x ranks per node
NTHREADS=@@{threads}                    # OpenMP threads per rank
SCRATCH_ROOT="@@{scratch}"        # empty = run in the submit directory
# ============================================================================

# ----------------------------- environment ----------------------------------
@@{modules}
# A small default stack is the classic cause of "forrtl: severe (174): SIGSEGV"
ulimit -s unlimited 2>/dev/null || {
    echo "WARNING: 'ulimit -s unlimited' not allowed, stack = $(ulimit -s) kB" >&2
    ulimit -s "$(ulimit -H -s)" 2>/dev/null || true
}
@@{omp}

WORKDIR="@@{workdir}"
JOBID="@@{jobid}"
cd "${WORKDIR}" || { echo "cannot cd to ${WORKDIR}" >&2; exit 1; }
echo "Job ${JOBID} (@@{program}) started $(date '+%F %T') on @@{nodelist}"
echo "MPI ranks: ${NP}   OpenMP threads per rank: ${NTHREADS}   workdir: ${WORKDIR}"

# ----------------------------- pre-run checks -------------------------------
fail=0
need_file() {   # need_file FILE : must exist, be non-empty and use LF line endings
    if [ ! -s "$1" ]; then
        echo "ERROR: input file '$1' is missing or empty in ${PWD}" >&2; fail=1; return
    fi
    if grep -q $'\r' "$1"; then
        echo "ERROR: '$1' has Windows (CRLF) line endings -> run: dos2unix $1" >&2; fail=1
    fi
}
command -v "${EXE%% *}" >/dev/null 2>&1 || \
    echo "WARNING: '${EXE%% *}' not found in PATH (module not loaded?)" >&2
"""

VASP_CHECK = r"""
for f in INCAR POSCAR POTCAR; do need_file "$f"; done
if grep -qiE '^[[:space:]]*KSPACING[[:space:]]*=' INCAR 2>/dev/null; then
    echo "KSPACING found in INCAR: KPOINTS file not required"
else
    need_file KPOINTS
fi
if [ -s POTCAR ] && [ -s POSCAR ]; then
    npot=$(grep -c "End of Dataset" POTCAR || true)
    nspec=$(sed -n '6p' POSCAR | awk '{print NF}')
    if [ "${nspec}" -gt 0 ] && ! sed -n '6p' POSCAR | grep -qE '^[[:space:]]*[0-9]' \
       && [ "${npot}" != "${nspec}" ]; then
        echo "ERROR: POTCAR has ${npot} datasets but POSCAR line 6 lists ${nspec} species" >&2; fail=1
    fi
fi
"""

GAUSSIAN_CHECK = r"""
need_file "${INPUT}"
if [ -s "${INPUT}" ]; then
    nps=$(grep -ioE '^%nprocshared=[0-9]+' "${INPUT}" | head -n 1 | cut -d= -f2)
    [ -n "${nps}" ] && [ "${nps}" -gt "${NTHREADS}" ] && {
        echo "ERROR: %NProcShared=${nps} > ${NTHREADS} cores requested" >&2; fail=1; }
    [ -z "${nps}" ] && ! grep -qiE '^%cpu=' "${INPUT}" && \
        echo "WARNING: no %NProcShared/%CPU in ${INPUT}: Gaussian uses 1 core (or GAUSS_PDEF/GAUSS_CDEF)" >&2
    grep -qiE '^%mem=' "${INPUT}" || \
        echo "WARNING: no %Mem in ${INPUT}: Gaussian uses its small default memory" >&2
    [ -n "$(tail -c 2 "${INPUT}" | tr -d '\n')" ] && \
        echo "WARNING: ${INPUT} should end with a blank line" >&2
fi
"""

CP2K_CHECK = r"""
need_file "${INPUT}"
[ -n "${CP2K_DATA_DIR:-}" ] || echo "NOTE: CP2K_DATA_DIR not set; BASIS/POTENTIAL files must be found via the input paths" >&2
"""

OPENMX_CHECK = r"""
need_file "${INPUT}"
if [ -s "${INPUT}" ]; then
    datapath=$(awk 'tolower($1)=="data.path"{print $2; exit}' "${INPUT}")
    if [ -n "${datapath}" ] && [ ! -d "${datapath}" ]; then
        echo "ERROR: DATA.PATH '${datapath}' is not a directory on this node (DFT_DATA19?)" >&2; fail=1
    fi
fi
"""

SCRATCH_GAUSSIAN = r"""
# ----------------------------- Gaussian scratch -----------------------------
# .rwf/.int/.d2e go to node-local disk; the .chk stays where %Chk points.
if [ -n "${SCRATCH_ROOT}" ]; then
    export GAUSS_SCRDIR="${SCRATCH_ROOT}/${USER}/g16_${JOBID}"
    mkdir -p "${GAUSS_SCRDIR}" || { echo "cannot create ${GAUSS_SCRDIR}" >&2; exit 1; }
    trap 'rm -rf "${GAUSS_SCRDIR}"' EXIT
fi
echo "GAUSS_SCRDIR=${GAUSS_SCRDIR:-<not set>}"
"""

SCRATCH_GENERIC = r"""
# ----------------------------- optional scratch -----------------------------
# Copies the submit directory to scratch, runs there, copies everything back
# on exit (also after a crash). The scratch path must be visible to every MPI
# rank: node-local disks are only safe for single-node jobs.
if [ -n "${SCRATCH_ROOT}" ]; then
    SCR="${SCRATCH_ROOT}/${USER}/@@{name}_${JOBID}"
    mkdir -p "${SCR}" || { echo "cannot create ${SCR}" >&2; exit 1; }
    cp -a "${WORKDIR}/." "${SCR}/"
    copy_back() { cp -a "${SCR}/." "${WORKDIR}/" && rm -rf "${SCR}"; }
    trap copy_back EXIT
    cd "${SCR}" || exit 1
    echo "running in scratch ${SCR}"
fi
"""

RUN = r"""
[ "${fail}" -eq 0 ] || { echo "pre-run checks failed, nothing was started" >&2; exit 1; }

# ----------------------------- run ------------------------------------------
t0=$(date +%s)
@@{runline}
status=$?
t1=$(date +%s)
echo "@@{program} exit code ${status}, elapsed $((t1 - t0)) s"
"""

VASP_POST = r"""
# ----------------------------- post-run check -------------------------------
ok=0
if grep -q "General timing and accounting" OUTCAR 2>/dev/null; then
    echo "CHECK: VASP finished normally (timing block in OUTCAR)"
    nsw=$(awk '$1=="NSW" && $2=="="{print $3; exit}' OUTCAR)
    ibrion=$(awk '$1=="IBRION" && $2=="="{print $3; exit}' OUTCAR)
    case "${ibrion}" in
        1|2|3)
            if [ "${nsw:-0}" -gt 0 ]; then
                if grep -q "reached required accuracy" OUTCAR; then
                    echo "CHECK: ionic relaxation converged (reached required accuracy)"
                else
                    echo "CHECK: relaxation stopped after NSW=${nsw} steps WITHOUT reaching EDIFFG" >&2; ok=2
                fi
            fi ;;
    esac
else
    echo "CHECK: no timing block in OUTCAR -> VASP crashed or was killed" >&2; ok=1
fi
grep -E "ZBRENT|EDDDAV|EDDRMM|VERY BAD NEWS|BRMIX|Sub-Space-Matrix|ZPOTRF|PRICEL|internal error" \
    "${OUTPUT}" 2>/dev/null | sort | uniq -c | head -n 10
[ "${status}" -eq 0 ] && [ "${ok}" -eq 1 ] && status=1
exit "${status}"
"""

GAUSSIAN_POST = r"""
# ----------------------------- post-run check -------------------------------
nnorm=$(grep -c "Normal termination of Gaussian" "${OUTPUT}" 2>/dev/null || true)
if grep -q "Error termination" "${OUTPUT}" 2>/dev/null; then
    echo "CHECK: Error termination:" >&2
    grep -B 3 "Error termination" "${OUTPUT}" | tail -n 8 >&2
    [ "${status}" -eq 0 ] && status=1
elif [ "${nnorm:-0}" -gt 0 ]; then
    echo "CHECK: ${nnorm} x Normal termination of Gaussian (one per Link1 job step)"
else
    echo "CHECK: no termination line -> killed (time/memory limit?)" >&2
    [ "${status}" -eq 0 ] && status=1
fi
grep "imaginary frequencies" "${OUTPUT}" 2>/dev/null | tail -n 1
exit "${status}"
"""

CP2K_POST = r"""
# ----------------------------- post-run check -------------------------------
if grep -q "PROGRAM ENDED AT" "${OUTPUT}" 2>/dev/null; then
    echo "CHECK: CP2K ended normally (PROGRAM ENDED AT)"
else
    echo "CHECK: no 'PROGRAM ENDED AT' -> crashed, aborted or killed" >&2
    grep -A 6 "ABORT" "${OUTPUT}" 2>/dev/null | head -n 12 >&2
    [ "${status}" -eq 0 ] && status=1
fi
nbad=$(grep -c "SCF run NOT converged" "${OUTPUT}" 2>/dev/null || true)
[ "${nbad:-0}" -gt 0 ] && echo "CHECK: 'SCF run NOT converged' appears ${nbad} times" >&2
grep -E "GEOMETRY OPTIMIZATION COMPLETED|MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED" "${OUTPUT}" 2>/dev/null | tail -n 1
grep -E "Total number of message passing processes|Number of threads for this process" "${OUTPUT}" 2>/dev/null | head -n 2
exit "${status}"
"""

OPENMX_POST = r"""
# ----------------------------- post-run check -------------------------------
SYSNAME=$(awk 'tolower($1)=="system.name"{print $2; exit}' "${INPUT}")
if [ -n "${SYSNAME}" ] && grep -q "Elapsed.Time." "${SYSNAME}.out" 2>/dev/null; then
    echo "CHECK: OpenMX finished (${SYSNAME}.out has Elapsed.Time.)"
    grep "Utot\." "${SYSNAME}.out" | tail -n 1
else
    echo "CHECK: ${SYSNAME:-System.Name}.out incomplete -> crashed or killed; see ${OUTPUT}" >&2
    [ "${status}" -eq 0 ] && status=1
fi
if [ -f "${SYSNAME}.md" ]; then
    nbad=$(grep -o "scf_conv=0" "${SYSNAME}.md" | wc -l)
    [ "${nbad}" -gt 0 ] && echo "CHECK: ${nbad} step(s) with scf_conv=0 in ${SYSNAME}.md (SCF not converged)" >&2
fi
exit "${status}"
"""


def build_omp(a):
    if a.program == "gaussian":
        return ("# Gaussian ignores OMP_NUM_THREADS: cores come from %NProcShared/%CPU\n"
                "export OMP_STACKSIZE=@@{stack}")
    lines = []
    if a.scheduler == "slurm":
        lines.append('export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-${NTHREADS}}"')
    else:
        lines.append('export OMP_NUM_THREADS="${NTHREADS}"')
    if a.threads > 1:
        lines += ["export OMP_PLACES=cores        # one thread per physical core",
                  "export OMP_PROC_BIND=close     # keep the threads of a rank together"]
    lines.append("export OMP_STACKSIZE=@@{stack}     # private stack of every OpenMP thread")
    return "\n".join(lines)


def build_modules(a):
    out = []
    if a.purge:
        out.append("module purge")
    for m in a.module:
        out.append("module load %s" % m)
    if not a.module:
        examples = {
            "vasp": "intel-oneapi-mpi intel-oneapi-mkl vasp/6.5.1",
            "gaussian": "gaussian/16-C.02",
            "cp2k": "cp2k/2026.2",
            "openmx": "openmx/4.0",
        }
        out.append("# module load %s     # <- site specific, edit" % examples[a.program])
    for src in a.source:
        out.append('source "%s"' % src)
    if a.program == "gaussian" and not a.module and not a.source:
        out.append('# export g16root=/opt/gaussian; source "${g16root}/g16/bsd/g16.profile"')
    return "\n".join(out)


def build_runline(a):
    if a.program == "gaussian":
        return '${EXE} < "${INPUT}" > "${OUTPUT}" 2>&1'
    pre = mpi_prefix(a)
    if a.program == "vasp":
        return pre + ' ${EXE} > "${OUTPUT}" 2>&1'
    if a.program == "cp2k":
        return pre + ' ${EXE} -i "${INPUT}" -o "${OUTPUT}" > "${OUTPUT%.out}.stdout" 2>&1'
    # OpenMX: -nt = OpenMP threads per rank
    return pre + ' ${EXE} "${INPUT}" -nt "${OMP_NUM_THREADS}" > "${OUTPUT}" 2>&1'


def default_output(a):
    if a.program == "vasp":
        return "vasp.out"
    stem = os.path.splitext(os.path.basename(a.input))[0] or "job"
    return {"gaussian": stem + ".log", "cp2k": stem + ".out", "openmx": stem + ".std"}[a.program]


def build_script(a, argv):
    sv = SCHED_VARS[a.scheduler]
    header = {"slurm": header_slurm, "pbs": header_pbs, "lsf": header_lsf}[a.scheduler](a)
    stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    top = ["#!/bin/bash",
           "# " + "=" * 74,
           "# %s job script for %s, generated by mkjob.py (DFT-Tutorials) %s"
           % (a.program, a.scheduler.upper(), stamp),
           "#   python mkjob.py %s" % " ".join(shlex.quote(x) for x in argv),
           "# Resources: %d node(s) x %d MPI rank(s) x %d thread(s) = %d cores"
           % (a.nodes, a.ntasks, a.threads, a.nodes * a.ntasks * a.threads),
           "# Submit   : %s %s" % (sv["submit"], a.outname if a.output != "-" else "job.sh"),
           "# Check partition/queue, account, modules and launcher against your site guide.",
           "# " + "=" * 74]
    checks = {"vasp": VASP_CHECK, "gaussian": GAUSSIAN_CHECK,
              "cp2k": CP2K_CHECK, "openmx": OPENMX_CHECK}[a.program]
    post = {"vasp": VASP_POST, "gaussian": GAUSSIAN_POST,
            "cp2k": CP2K_POST, "openmx": OPENMX_POST}[a.program]
    scratch = SCRATCH_GAUSSIAN if a.program == "gaussian" else SCRATCH_GENERIC
    body = COMMON_TOP + checks + scratch + RUN + post
    body = body.replace("@@{omp}", build_omp(a)).replace("@@{runline}", build_runline(a))
    values = {
        "exe": a.exe, "input": a.input, "output": a.out_file,
        "np": str(a.nodes * a.ntasks), "threads": str(a.threads),
        "scratch": a.scratch or "", "modules": build_modules(a),
        "stack": a.omp_stacksize, "workdir": sv["workdir"], "jobid": sv["jobid"],
        "nodelist": sv["nodes"], "program": a.program, "name": a.name,
    }
    try:
        text = JobTemplate(body).substitute(values)
    except (KeyError, ValueError) as exc:   # programming error in a template
        die("template substitution failed: %s" % exc)
    return "\n".join(top + header) + "\n" + text


def build_parser():
    p = argparse.ArgumentParser(
        prog="mkjob.py",
        description="Generate a Slurm/PBS/LSF job script for VASP, Gaussian 16, CP2K or OpenMX.",
        epilog=("Examples:\n"
                "  python mkjob.py vasp -N 2 -n 64 -t 24:00:00 -q cpu --module vasp/6.5.1\n"
                "  python mkjob.py gaussian -i mol.gjf -c 16 --mem 40G --scratch '${TMPDIR:-/tmp}'\n"
                "  python mkjob.py cp2k -i md.inp -n 32 -c 2 -s pbs -o job.pbs\n"
                "  python mkjob.py openmx -i Si.dat -n 16 -c 2 -s lsf\n"
                "Gaussian runs on one node: -c sets the cores (= %NProcShared)."),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("program", help="vasp | gaussian (g16) | cp2k | openmx")
    p.add_argument("-s", "--scheduler", choices=SCHEDULERS, default="slurm",
                   help="batch system (default: slurm)")
    p.add_argument("-i", "--input", help="main input file (Gaussian .gjf, CP2K .inp, OpenMX .dat)")
    p.add_argument("--output", dest="out_file",
                   help="program output/log file (default: vasp.out, <input>.log/.out/.std)")
    p.add_argument("-N", "--nodes", type=int, default=1, help="number of nodes (default 1)")
    p.add_argument("-n", "--ntasks-per-node", dest="ntasks", type=int,
                   help="MPI ranks per node (default 32; Gaussian: 1)")
    p.add_argument("-c", "--threads", type=int,
                   help="OpenMP threads per rank = cpus-per-task (default 1; Gaussian: cores, default 16)")
    p.add_argument("--mem", type=parse_mem, dest="mem_mib",
                   help="memory per node, e.g. 180G (Gaussian: also suggests %%Mem)")
    p.add_argument("-t", "--time", type=parse_time, dest="time_s", default=parse_time("24:00:00"),
                   help="wall time D-HH:MM:SS, HH:MM:SS, HH:MM or minutes (default 24:00:00)")
    p.add_argument("-q", "--queue", "--partition", dest="queue", help="partition (Slurm) / queue (PBS, LSF)")
    p.add_argument("-A", "--account", help="account / project to charge")
    p.add_argument("--qos", help="Slurm QOS")
    p.add_argument("-J", "--name", help="job name (default: program or input stem)")
    p.add_argument("--gpus", type=int, default=0, help="GPUs per node (VASP GPU: 1 rank per GPU)")
    p.add_argument("--exclusive", action="store_true", help="request whole nodes exclusively")
    p.add_argument("--torque", action="store_true", help="PBS: Torque syntax nodes=N:ppn=P instead of select=")
    p.add_argument("--module", action="append", default=[], metavar="MOD",
                   help="module to load (repeatable), e.g. --module 'intel-oneapi-mpi vasp/6.5.1'")
    p.add_argument("--purge", action="store_true", help="add 'module purge' before the module loads")
    p.add_argument("--source", action="append", default=[], metavar="FILE",
                   help="file to source, e.g. CP2K toolchain 'install/setup' (repeatable)")
    p.add_argument("--exe", help="executable (default: vasp_std, g16, cp2k.psmp, openmx)")
    p.add_argument("--launcher", choices=("srun", "mpirun"),
                   help="MPI launcher (default: srun for Slurm, mpirun otherwise)")
    p.add_argument("--scratch", metavar="DIR",
                   help="scratch root, e.g. '${TMPDIR:-/tmp}' or /scratch (Gaussian: GAUSS_SCRDIR)")
    p.add_argument("--omp-stacksize", default="512m", help="OMP_STACKSIZE (default 512m)")
    p.add_argument("-o", "--output-script", dest="output", default=None,
                   help="script file to write (default job_<program>_<scheduler>.sh; '-' = stdout)")
    p.add_argument("-f", "--force", action="store_true", help="overwrite an existing script file")
    return p


def finalize(a):
    prog = a.program.lower()
    prog = ALIASES.get(prog, prog)
    if prog not in PROGRAMS:
        die("unknown program %r (choose from %s)" % (a.program, ", ".join(PROGRAMS)))
    a.program = prog
    d = DEFAULTS[prog]
    if a.launcher is None:
        a.launcher = "srun" if a.scheduler == "slurm" else "mpirun"
    if prog == "gaussian":
        cores = a.threads or a.ntasks or d["threads"]
        if a.nodes != 1:
            print("NOTE: Gaussian 16 (non-Linda) runs on one node: --nodes set to 1", file=sys.stderr)
        a.nodes, a.ntasks, a.threads = 1, 1, cores
    else:
        a.ntasks = a.ntasks or d["ntasks"]
        a.threads = a.threads or d["threads"]
    for key in ("nodes", "ntasks", "threads"):
        if getattr(a, key) < 1:
            die("--%s must be >= 1" % key)
    if a.gpus < 0:
        die("--gpus must be >= 0")
    a.exe = a.exe or d["exe"]
    a.input = a.input if a.input is not None else d["input"]
    if prog != "vasp" and not a.input:
        die("%s needs an input file (-i)" % prog)
    if a.input and re.search(r"[\s\"'$`]", a.input):
        die("please use an input file name without spaces or quotes: %r" % a.input)
    a.out_file = a.out_file or default_output(a)
    if a.name is None:
        stem = os.path.splitext(os.path.basename(a.input))[0] if a.input else prog
        a.name = stem or prog
    a.name = clean_name(a.name)
    a.outname = a.output or "job_%s_%s.sh" % (prog, a.scheduler)
    if a.scratch and a.nodes > 1 and prog != "gaussian":
        print("WARNING: --scratch with %d nodes: the directory must be on a file system that "
              "all nodes can see (node-local /tmp is NOT)" % a.nodes, file=sys.stderr)
    if prog == "gaussian" and a.mem_mib:
        pct = int(a.mem_mib * 0.8 / 1024)
        print("Gaussian: put %%NProcShared=%d and about %%Mem=%dGB (80%% of --mem) in %s"
              % (a.threads, max(1, pct), a.input), file=sys.stderr)
    if prog == "vasp" and a.threads > 1:
        print("NOTE: VASP with OpenMP threads needs an OpenMP build; VASP then sets NCORE=1 itself",
              file=sys.stderr)
    return a


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    a = finalize(build_parser().parse_args(argv))
    text = build_script(a, argv)
    if a.output == "-":
        sys.stdout.write(text)
        return 0
    path = a.outname
    if os.path.exists(path) and not a.force:
        die("%s exists (use -f to overwrite or -o to choose another name)" % path)
    with open(path, "w", newline="\n", encoding="utf-8") as fh:
        fh.write(text)
    try:
        os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP)
    except OSError:
        pass
    print("written %s  (%s, %d node(s) x %d rank(s) x %d thread(s))"
          % (path, a.scheduler, a.nodes, a.ntasks, a.threads))
    print("submit with: %s %s" % (SCHED_VARS[a.scheduler]["submit"], path))
    return 0


if __name__ == "__main__":
    sys.exit(main())
