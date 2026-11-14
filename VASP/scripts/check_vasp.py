#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_vasp.py -- health check of a (running or finished) VASP job.
                 VASP 作业健康检查：收敛情况、受力、磁矩、耗时与常见报错诊断。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

What it reports
---------------
* VASP version, task type (IBRION/ISIF/NSW), number of ionic steps and the
  number of electronic steps of every ionic step (SCF steps hitting NELM are
  flagged as unconverged);
* 'reached required accuracy' / NSW exhausted / job finished or still running;
* last energies (TOTEN, energy(sigma->0)), T*S per atom (SIGMA check),
  E-fermi, external pressure (+ Pulay stress) for ISIF >= 3;
* maximum force, taking Selective dynamics flags of POSCAR/CONTCAR into account
  (the force is projected onto the allowed displacement directions, exact also
  for non-orthogonal cells) and compared with EDIFFG;
* total and largest site magnetic moments (collinear and non-collinear);
* elapsed time, average time per ionic step, memory;
* WARNING / ADVICE boxes of VASP 6 and known error strings in OUTCAR and in the
  stdout files (vasp.out, slurm-*.out, *.o<jobid>, ...), each with a fix hint
  in Chinese (ZBRENT, BRMIX, EDDDAV, Sub-Space-Matrix, IBZKPT, Tetrahedron,
  PRICEL, ZPOTRF, time limit, out-of-memory, ...).

Usage examples
--------------
    python check_vasp.py                    # current directory
    python check_vasp.py relax/ scf/ band/  # several jobs
    python check_vasp.py */ --brief         # one line per directory
    python check_vasp.py neb/               # 00..NN image directories are detected
    python check_vasp.py -o vasp.out --json report.json
    python check_vasp.py --strict && echo converged   # exit code 1 if not converged

Exit status: 0 = OK, 1 = errors found (or not converged with --strict),
2 = no VASP output found.

Dependencies: python >= 3.9, standard library only (numpy is used if present).
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import os
import re
import sys
import time
from pathlib import Path

__version__ = "1.0.0"

# ---------------------------------------------------------------------------
# Known error / warning strings (VASP 5.4 - 6.6 and batch systems) -> advice
# (key, severity, [substrings], Chinese advice)
# ---------------------------------------------------------------------------
PATTERNS = [
    ("ZBRENT", "error", ["ZBRENT: fatal error in bracketing", "ZBRENT: fatal internal in",
                         "ZBRENT: fatal error: bracketing interval incorrect"],
     "离子优化线搜索失败：cp CONTCAR POSCAR 后续算；减小 POTIM（如 0.2）；EDIFF 收紧到 1E-6 或更小；"
     "改用 IBRION=1，或 VTST 的 FIRE（IBRION=3, POTIM=0, IOPT=7）。"),
    ("BRMIX", "error", ["BRMIX: very serious problems"],
     "电荷混合严重失败：删除旧 CHGCAR/WAVECAR 重算；检查对称性（ISYM=0 或 SYMPREC=1E-4~1E-6）；"
     "必要时减小 AMIX/BMIX。"),
    ("EDDDAV", "error", ["Error EDDDAV: Call to ZHEGV failed", "EDDDAV: Call to ZHEGV failed"],
     "Davidson 对角化失败：先检查原子是否重叠；换 ALGO=Fast 或 All；ISTART=0 不读 WAVECAR；减小 POTIM。"),
    ("SUBSPACE", "warning", ["Sub-Space-Matrix is not hermitian in DAV"],
     "Davidson 子空间矩阵非厄米：偶尔出现可忽略；频繁出现时换 ALGO=Fast/All 并检查结构与步长。"),
    ("EDDRMM", "error", ["WARNING in EDDRMM: call to ZHEGV failed"],
     "RMM-DIIS 失败：改用 ALGO=Normal（Davidson），或减小 POTIM、检查结构。"),
    ("ZPOTRF", "error", ["LAPACK: Routine ZPOTRF failed", "Routine ZPOTRF ZTRTRI"],
     "波函数正交化失败：多因原子过近或上一步结构异常；检查 CONTCAR，ISTART=0 重启，减小 POTIM，ALGO=Normal。"),
    ("PSSYEVX", "error", ["ERROR in subspace rotation PSSYEVX", "ERROR in subspace rotation PDSYEVX"],
     "ScaLAPACK 子空间旋转失败：尝试 ALGO=Normal，或在 INCAR 中设 LSCALAPACK = .FALSE.。"),
    ("EDDIAG", "error", ["ERROR EDDIAG: Call to routine ZHEEV failed",
                         "ERROR in EDDIAG: call to ZHEEV/ZHEEVX/DSYEV/DSYEVX failed"],
     "子空间对角化失败：换 ALGO（Normal 与 Fast 互换），检查结构。"),
    ("EDWAV", "error", ["EDWAV: internal error, the gradient is not orthogonal"],
     "梯度不正交：试 ALGO=Normal 或 Fast；反复出现时检查编译（降低优化级别）与 MPI/数学库。"),
    ("IBZKPT", "error", ["Reciprocal lattice and k-lattice"],
     "k 网格与晶格对称性不相容（常见于六角/三方体系的偶数 Monkhorst-Pack 网格）：改用 Gamma 中心网格或修改网格。"),
    ("TETRAHEDRON", "error", ["Tetrahedron method fails", "tetrahedron method fails",
                              "subroutine BZINTS", "BZINTS: Tetrahedron",
                              "Routine TETIRR needs special values"],
     "ISMEAR=-5 需要至少 4 个 k 点的 Gamma 中心网格：分子、Gamma 点和 slab 计算改用 ISMEAR=0, SIGMA=0.05。"),
    ("DENTET", "warning", ["DENTET"],
     "四面体法无法精确确定费米能级：增加 k 点，或改用 ISMEAR=0。"),
    ("PRICEL", "error", ["internal error in subroutine PRICEL",
                         "PRICELV: current lattice and primitive lattice are incommensurate"],
     "结构对称性有噪声：用 spglib/vaspkit 规整结构，调 SYMPREC（1E-4 或 1E-6），必要时 ISYM=0。"),
    ("ROTMAT", "error", ["Found some non-integer element in rotation matrix", "SGRCON"],
     "对称操作矩阵非整数（结构对称性噪声）：规整结构，调 SYMPREC，或 ISYM=0。"),
    ("INV_ROT", "error", ["rotation matrix was not found (increase SYMPREC)"],
     "找不到逆对称操作：增大 SYMPREC（如 1E-4），或 ISYM=0。"),
    ("SYMPREC_NOISE", "warning", ["determination of the symmetry of your systems shows a strong"],
     "对称性判断对 SYMPREC 很敏感：规整坐标或调整 SYMPREC。"),
    ("POSMAP", "error", ["POSMAP"],
     "对称等价原子映射失败：多为原子过近或 SYMPREC 不当；检查结构，调整 SYMPREC 或 ISYM=0。"),
    ("RHOSYG", "error", ["RHOSYG"],
     "电荷密度对称化失败：ISYM=0 或调整 SYMPREC。"),
    ("POINT_GROUP", "error", ["group operation missing"],
     "点群操作缺失：规整结构或 ISYM=0。"),
    ("BRAVAIS", "error", ["Inconsistent Bravais lattice"],
     "晶胞与识别出的布拉维格子不一致：用标准化晶胞或调 SYMPREC。"),
    ("TOO_FEW_BANDS", "error", ["TOO FEW BANDS", "number of bands is not sufficient"],
     "能带数不足：增大 NBANDS。"),
    ("HIGH_BAND_OCC", "warning", ["Your highest band is occupied at some k-points"],
     "最高能带仍有占据：增大 NBANDS（光学、金属、高温 MD 尤其注意）。"),
    ("NBANDS_CHANGED", "info", ["The number of bands has been changed"],
     "VASP 为适配并行自动增加了 NBANDS；LOBSTER/LOPTICS/MLFF 等需要一致 NBANDS 时请手动设置为并行组数的整数倍。"),
    ("READ_ERROR", "error", ["Error reading item"],
     "INCAR 某标签读取失败：检查拼写与数值；dos2unix 去掉 Windows 换行，不要用 Tab；"
     "MAGMOM 个数需等于原子数（非共线为 3×NIONS），可用 N*value 简写。"),
    ("TRIPLE_PRODUCT", "error", ["ERROR: the triple product of the basis vectors"],
     "晶格矢量构成左手系（三重积为负）：交换两个晶格矢量以及对应的坐标分量。"),
    ("REAL_OPTLAY", "error", ["REAL_OPTLAY: internal error", "REAL_OPT: internal ERROR"],
     "实空间投影出错：设 LREAL=.FALSE.（或保持 LREAL=Auto 但检查 ENCUT/PREC）。"),
    ("RSPHER", "error", ["ERROR RSPHER"],
     "实空间投影出错：设 LREAL=.FALSE.。"),
    ("BRIONS", "error", ["BRIONS problems: POTIM should be increased"],
     "增大 POTIM，或改用 IBRION=1/2。"),
    ("AMIN", "warning", ["One of the lattice vectors is very long (>50 A), but AMIN"],
     "存在很长的晶格矢量（厚真空层）：设置 AMIN = 0.01。"),
    ("SBESSEL", "error", ["SBESSELITER : nicht konvergent"],
     "球贝塞尔函数迭代不收敛：通常与 LREAL 有关，设 LREAL=.FALSE.。"),
    ("ELF_KPAR", "error", ["ELF: KPAR>1 not implemented"],
     "ELF 计算不支持 KPAR>1：设 KPAR=1。"),
    ("ELF_NCL", "error", ["ELF not implemented for non collinear case"],
     "ELF 不支持非共线计算：改用共线计算。"),
    ("DFPT_NCORE", "error", ["PEAD routines do not work for NCORE", "remove the tag NPAR from the INCAR file"],
     "DFPT/LCALCEPS 不支持 NCORE/NPAR>1：删除 INCAR 中的 NCORE 和 NPAR。"),
    ("HNFORM", "error", ["HNFORM: k-point generating"],
     "k 点生成失败：改变 k 网格，或使用标准化晶胞。"),
    ("KMESH_SYMM", "error", ["Fatal error detecting k-mesh", "Fatal error: unable to match k-point"],
     "k 网格与对称性不匹配：改用 Gamma 中心网格，或 ISYM=0。"),
    ("KSHIFT", "error", ["Could not get correct shifts"],
     "无法得到正确的 k 点平移：改用 Gamma 中心网格。"),
    ("WAVECAR", "error", ["while reading WAVECAR", "while reading plane"],
     "WAVECAR 与当前设置（ENCUT/晶胞/k 点）不匹配：删除 WAVECAR 或设 ISTART=0。"),
    ("SET_CORE_WF", "error", ["internal error in SET_CORE_WF"],
     "POTCAR 有问题（版本不兼容或文件损坏）：重新生成 POTCAR。"),
    ("XC_TABLE", "error", ["supplied exchange-correlation table"],
     "交换关联表与 POTCAR 不匹配：检查 POTCAR 的 LEXCH 与 INCAR 中 GGA/METAGGA。"),
    ("ALGO_TET", "warning", ["ALGO=A and IALGO=5X tend to fail"],
     "ALGO=All 与四面体法（ISMEAR=-5）组合不稳定：换 ALGO=Normal 或 ISMEAR=0。"),
    ("SCF_NOT_CONVERGED", "warning", ["electronic self-consistency was not achieved"],
     "有离子步的电子步未收敛：增大 NELM，换 ALGO=All/Normal，减小 AMIX，或增加 NBANDS。"),
    ("IONS_CLOSE", "error", ["distance between some ions is very small"],
     "原子间距过小：检查结构（重复原子、坐标单位或 Direct/Cartesian 写错）。"),
    ("VERY_BAD_NEWS", "error", ["VERY BAD NEWS"],
     "VASP 内部错误：看该行上下文给出的子程序名。"),
    ("MPI_ABORT", "error", ["BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES",
                            "APPLICATION TERMINATED WITH THE EXIT STRING", "MPI_ABORT was invoked"],
     "MPI 进程异常退出：看前面的报错；常见原因是内存不足或段错误（作业脚本中 ulimit -s unlimited）。"),
    ("SEGFAULT", "error", ["forrtl: severe (174)", "SIGSEGV", "Segmentation fault"],
     "段错误：ulimit -s unlimited；检查内存、编译器与数学库版本是否匹配。"),
    ("OOM", "error", ["oom-kill", "Out Of Memory", "out-of-memory", "OOM Killed"],
     "内存不足：减少每节点 MPI 进程数、合理设置 NCORE/KPAR、LREAL=Auto，或申请更多内存/节点。"),
    ("TIME_LIMIT", "error", ["DUE TO TIME LIMIT", "exceeded limit", "TIME LIMIT EXCEEDED"],
     "作业超时被杀：cp CONTCAR POSCAR（MD 还要保留速度）后续算，或延长 walltime。"),
]

ADVICE_BANNERS = ("W    W    AA    RRRRR", "EEEEEEE  RRRRRR", "ADVICE TO THIS USER RUNNING VASP")
ART_RE = re.compile(r"^[WANRIG!EO#\s]*$")
FLOAT = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?"


def fnum(s):
    try:
        return float(s.replace("D", "E").replace("d", "e"))
    except (ValueError, AttributeError):
        return None


# ---------------------------------------------------------------------------
# POSCAR (for selective dynamics)
# ---------------------------------------------------------------------------
def read_poscar_constraints(path: Path):
    """Return (lattice 3x3 list, flags list[N][3] or None, natoms) or None."""
    try:
        lines = path.read_text(errors="replace").splitlines()
        scale = [float(x) for x in lines[1].split()[:3]]
        lat = [[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)]
        i = 5
        tok = lines[i].split()
        if not all(t.isdigit() for t in tok):
            i += 1
        counts = [int(x) for x in lines[i].split()]
        i += 1
        sel = lines[i].strip()[:1] in ("s", "S")
        if sel:
            i += 1
        i += 1
        nat = sum(counts)
        if len(scale) == 1 and scale[0] > 0:
            lat = [[x * scale[0] for x in v] for v in lat]
        elif len(scale) == 3:
            lat = [[v[k] * scale[k] for k in range(3)] for v in lat]
        flags = None
        if sel:
            flags = []
            for j in range(nat):
                t = lines[i + j].split()
                flags.append([x.upper().startswith("T") for x in t[3:6]] if len(t) >= 6 else [True] * 3)
        return lat, flags, nat
    except (OSError, IndexError, ValueError):
        return None


def project_force(force, lat, flag):
    """Project a Cartesian force onto the subspace of allowed displacements
    spanned by the lattice vectors whose selective-dynamics flag is T."""
    free = [lat[k] for k in range(3) if flag[k]]
    if len(free) == 3:
        return force
    if not free:
        return [0.0, 0.0, 0.0]
    # Gram-Schmidt orthonormal basis of the free directions
    basis = []
    for v in free:
        w = list(v)
        for b in basis:
            d = sum(w[i] * b[i] for i in range(3))
            w = [w[i] - d * b[i] for i in range(3)]
        n = math.sqrt(sum(x * x for x in w))
        if n > 1e-10:
            basis.append([x / n for x in w])
    out = [0.0, 0.0, 0.0]
    for b in basis:
        d = sum(force[i] * b[i] for i in range(3))
        out = [out[i] + d * b[i] for i in range(3)]
    return out


# ---------------------------------------------------------------------------
# OSZICAR
# ---------------------------------------------------------------------------
ELEC_RE = re.compile(r"^\s*([A-Za-z]{2,4})\s*:\s*(\d+)\s+(" + FLOAT + ")")
IONIC_RE = re.compile(r"^\s*(\d+)\s+(F=|T=)")


def parse_oszicar(path: Path) -> dict:
    res = {"ionic": [], "nelec": [], "algos": {}, "current_elec": 0}
    n = 0
    with open(path, errors="replace") as fh:
        for line in fh:
            m = ELEC_RE.match(line)
            if m:
                n += 1
                res["algos"][m.group(1)] = res["algos"].get(m.group(1), 0) + 1
                continue
            m = IONIC_RE.match(line)
            if m:
                d = {"step": int(m.group(1))}
                for key, pat in (("F", r"\bF=\s*(" + FLOAT + ")"), ("E0", r"\bE0=\s*(" + FLOAT + ")"),
                                 ("T", r"\bT=\s*(" + FLOAT + ")"), ("E", r"\bE=\s*(" + FLOAT + ")"),
                                 ("dE", r"d E\s*=\s*(" + FLOAT + ")")):
                    mm = re.search(pat, line)
                    if mm:
                        d[key] = fnum(mm.group(1))
                mm = re.search(r"mag=\s*(.*)$", line)
                if mm:
                    vals = [fnum(x) for x in re.findall(FLOAT, mm.group(1))]
                    vals = [v for v in vals if v is not None]
                    d["mag"] = vals[0] if len(vals) == 1 else vals
                res["ionic"].append(d)
                res["nelec"].append(n)
                n = 0
    res["current_elec"] = n
    return res


# ---------------------------------------------------------------------------
# OUTCAR (streamed: works for multi-GB MD files)
# ---------------------------------------------------------------------------
PARAM_RES = {
    "NIONS": re.compile(r"NIONS\s*=\s*(\d+)"),
    "NKPTS": re.compile(r"NKPTS\s*=\s*(\d+)"),
    "NBANDS": re.compile(r"NBANDS\s*=\s*(\d+)"),
    "NELM": re.compile(r"\bNELM\s*=\s*(\d+)"),
    "EDIFF": re.compile(r"\bEDIFF\s*=\s*(" + FLOAT + ")"),
    "EDIFFG": re.compile(r"\bEDIFFG\s*=\s*(" + FLOAT + ")"),
    "NSW": re.compile(r"\bNSW\s*=\s*(\d+)"),
    "IBRION": re.compile(r"\bIBRION\s*=\s*(-?\d+)"),
    "ISIF": re.compile(r"\bISIF\s*=\s*(\d+)"),
    "ISPIN": re.compile(r"\bISPIN\s*=\s*(\d+)"),
    "ISMEAR": re.compile(r"\bISMEAR\s*=\s*(-?\d+)"),
    "SIGMA": re.compile(r"\bSIGMA\s*=\s*(" + FLOAT + ")"),
    "ENCUT": re.compile(r"\bENCUT\s*=\s*(" + FLOAT + ")"),
    "IMAGES": re.compile(r"\bIMAGES\s*=\s*(\d+)"),
    "LSORBIT": re.compile(r"\bLSORBIT\s*=\s*([TF])"),
}


def scan_patterns(line: str, found: dict, src: str) -> None:
    for key, sev, subs, _adv in PATTERNS:
        for s in subs:
            if s in line:
                ent = found.setdefault(key, {"count": 0, "line": line.strip()[:160], "sources": set()})
                ent["count"] += 1
                ent["sources"].add(src)
                break


def _skip_line(line: str) -> bool:
    st = line.lstrip()
    return not st or st[0].isdigit() or st[0] in "-.+"


def parse_outcar(path: Path) -> dict:
    r = {"params": {}, "toten": [], "e0": [], "found": {}, "boxes": {}, "loop_real": [],
         "finished": False, "reached": False, "iterations": 0, "ionic_iter": 0,
         "forces": None, "magmom_tot": None, "site_mag": None, "efermi": None,
         "pressure": None, "pulay": None, "ts_last": None, "version": None}
    header = True
    nions = None
    with open(path, errors="replace") as fh:
        it = iter(fh)
        for lineno, line in enumerate(it):
            if r["version"] is None and lineno < 20:
                m = re.search(r"(vasp\.\d+\.\d+\.\d+)", line)
                if m:
                    r["version"] = m.group(1)
            if header:
                if "Iteration" in line and "(" in line:
                    header = False
                else:
                    for k, rx in PARAM_RES.items():
                        m = rx.search(line)
                        if m:
                            r["params"][k] = m.group(1)
                    if "NIONS" in line and nions is None and "NIONS" in r["params"]:
                        nions = int(r["params"]["NIONS"])
            if "Iteration" in line and "(" in line:
                m = re.search(r"Iteration\s*(\d+)\s*\(\s*(\d+)\)", line)
                if m:
                    r["ionic_iter"] = int(m.group(1))
                    r["iterations"] += 1
                continue
            if "EENTRO" in line:
                m = re.search(r"EENTRO\s*=\s*(" + FLOAT + ")", line)
                if m:
                    r["ts_last"] = fnum(m.group(1))
                continue
            if "FREE ENERGIE OF THE ION-ELECTRON SYSTEM" in line:
                for _ in range(6):
                    try:
                        nxt = next(it)
                    except StopIteration:
                        break
                    m = re.search(r"TOTEN\s*=\s*(" + FLOAT + ")", nxt)
                    if m:
                        r["toten"].append(fnum(m.group(1)))
                    m = re.search(r"energy\(sigma->0\)\s*=\s*(" + FLOAT + ")", nxt)
                    if m:
                        r["e0"].append(fnum(m.group(1)))
                        break
                continue
            if "TOTAL-FORCE" in line and "POSITION" in line:
                if nions is None:
                    continue
                try:
                    next(it)  # dashes
                    fr = []
                    for _ in range(nions):
                        t = next(it).split()
                        fr.append([fnum(x) for x in t[3:6]])
                    if all(None not in f and len(f) == 3 for f in fr):
                        r["forces"] = fr
                except StopIteration:
                    pass
                continue
            if "magnetization (x)" in line and nions:
                try:
                    rows = []
                    for _ in range(3):
                        next(it)
                    for _ in range(nions):
                        rows.append(fnum(next(it).split()[-1]))
                    r["site_mag"] = {"x": rows}
                except (StopIteration, IndexError):
                    pass
                continue
            if ("magnetization (y)" in line or "magnetization (z)" in line) and nions and r["site_mag"]:
                comp = "y" if "(y)" in line else "z"
                try:
                    rows = []
                    for _ in range(3):
                        next(it)
                    for _ in range(nions):
                        rows.append(fnum(next(it).split()[-1]))
                    r["site_mag"][comp] = rows
                except (StopIteration, IndexError):
                    pass
                continue
            if "number of electron" in line and "magnetization" in line:
                parts = line.split()
                try:
                    idx = parts.index("magnetization") + 1
                    vals = [v for v in (fnum(x) for x in parts[idx:]) if v is not None]
                    if vals:
                        r["magmom_tot"] = vals[0] if len(vals) == 1 else vals
                except ValueError:
                    pass
                continue
            if "E-fermi" in line:
                m = re.search(r"E-fermi\s*:\s*(" + FLOAT + ")", line)
                if m:
                    r["efermi"] = fnum(m.group(1))
                continue
            if "external pressure" in line:
                m = re.search(r"external pressure\s*=\s*(" + FLOAT + ")", line)
                p = re.search(r"Pullay stress\s*=\s*(" + FLOAT + ")", line)
                if m:
                    r["pressure"] = fnum(m.group(1))
                if p:
                    r["pulay"] = fnum(p.group(1))
                continue
            if "reached required accuracy" in line:
                r["reached"] = True
                continue
            if "LOOP+" in line:
                m = re.search(r"real time\s*(" + FLOAT + ")", line)
                if m:
                    r["loop_real"].append(fnum(m.group(1)))
                continue
            if "Elapsed time (sec)" in line:
                r["elapsed"] = fnum(line.split(":")[-1].strip())
                continue
            if "Total CPU time used (sec)" in line:
                r["cpu"] = fnum(line.split(":")[-1].strip())
                continue
            if "Maximum memory used (kb)" in line:
                r["maxmem_kb"] = fnum(line.split(":")[-1].strip().rstrip("."))
                continue
            if "General timing and accounting" in line:
                r["finished"] = True
                continue
            if any(b in line for b in ADVICE_BANNERS):
                kind = "ADVICE" if "ADVICE" in line else ("ERROR" if "EEEEEEE" in line else "WARNING")
                text = []
                try:
                    while True:
                        nxt = next(it)
                        if nxt.lstrip().startswith("-----") or not nxt.strip().startswith("|"):
                            break
                        inner = nxt.strip().strip("|").strip()
                        if not inner or ART_RE.match(inner):
                            continue
                        text.append(inner)
                except StopIteration:
                    pass
                if text:
                    # messages are wrapped over several box lines (e.g. "...BZINTS: Tetrahedron" /
                    # "method fails ..."), so match the joined text once per box
                    scan_patterns(" ".join(text), r["found"], "OUTCAR")
                msg = " ".join(text[:2])[:150] if text else "(empty box)"
                key = (kind, msg)
                r["boxes"][key] = r["boxes"].get(key, 0) + 1
                continue
            if line.lstrip().startswith("WARNING"):
                key = ("WARNING", line.strip()[:150])
                r["boxes"][key] = r["boxes"].get(key, 0) + 1
            if not _skip_line(line):
                scan_patterns(line, r["found"], "OUTCAR")
    return r


def scan_stdout(paths, found: dict) -> None:
    for p in paths:
        try:
            with open(p, errors="replace") as fh:
                box = []          # inner text of the current VASP 6 message box
                for line in fh:
                    st = line.strip()
                    if st.startswith("|"):
                        inner = st.strip("|").strip()
                        if inner and not ART_RE.match(inner):
                            box.append(inner)
                        continue
                    if box:
                        # wrapped box messages are matched on their joined text
                        scan_patterns(" ".join(box), found, p.name)
                        box = []
                    if not _skip_line(line):
                        scan_patterns(line, found, p.name)
                if box:
                    scan_patterns(" ".join(box), found, p.name)
        except OSError:
            continue


def stdout_candidates(d: Path, explicit):
    if explicit:
        return [Path(x) for x in explicit if Path(x).is_file()]
    names = ["vasp.out", "vasp.log", "vasp_out", "stdout", "out", "OUTPUT", "log", "job.out", "run.log"]
    c = [d / n for n in names if (d / n).is_file()]
    for pat in ("slurm-*.out", "*.o[0-9]*", "*.out", "*.log"):
        for f in sorted(glob.glob(str(d / pat)), key=os.path.getmtime)[-3:]:
            fp = Path(f)
            if fp not in c and fp.name not in ("OUTCAR",) and fp.is_file():
                c.append(fp)
    return c


# ---------------------------------------------------------------------------
# analysis of one directory
# ---------------------------------------------------------------------------
def analyse(d: Path, args) -> dict:
    rep = {"dir": str(d), "status": "unknown", "problems": [], "notes": []}
    outcar = d / "OUTCAR"
    oszicar = d / "OSZICAR"
    if not outcar.is_file() and not oszicar.is_file():
        rep["status"] = "no-output"
        return rep
    osz = parse_oszicar(oszicar) if oszicar.is_file() else None
    out = parse_outcar(outcar) if outcar.is_file() else None
    found = out["found"] if out else {}
    scan_stdout(stdout_candidates(d, args.stdout), found)
    p = out["params"] if out else {}

    def ip(k, default=None):
        try:
            return int(p[k])
        except (KeyError, ValueError):
            return default

    def fp(k, default=None):
        v = fnum(p.get(k)) if k in p else None
        return default if v is None else v

    nions, nelm, nsw, ibrion, isif = ip("NIONS"), ip("NELM", 60), ip("NSW", 0), ip("IBRION", -1), ip("ISIF", 2)
    ediffg = fp("EDIFFG")
    rep.update(version=out["version"] if out else None, nions=nions, nsw=nsw, ibrion=ibrion, isif=isif,
               nelm=nelm, ediffg=ediffg, ismear=ip("ISMEAR"), sigma=fp("SIGMA"), encut=fp("ENCUT"),
               nkpts=ip("NKPTS"), nbands=ip("NBANDS"), ispin=ip("ISPIN"))
    is_md = ibrion == 0
    relax = ibrion in (1, 2, 3) and nsw > 0

    # ionic / electronic steps
    if osz:
        nion = len(osz["ionic"])
        rep["n_ionic"] = nion
        rep["nelec"] = osz["nelec"]
        rep["algos"] = osz["algos"]
        unconv = [i + 1 for i, n in enumerate(osz["nelec"]) if nelm and n >= nelm]
        rep["scf_unconverged_steps"] = unconv
        if unconv:
            rep["problems"].append(f"{len(unconv)} ionic step(s) hit NELM={nelm} (SCF not converged): "
                                   f"{unconv[:10]}{' ...' if len(unconv) > 10 else ''}")
        if osz["ionic"]:
            last = osz["ionic"][-1]
            rep["E0_oszicar"] = last.get("E0")
            rep["F_oszicar"] = last.get("F")
            if "mag" in last:
                rep["mag_oszicar"] = last["mag"]
            if is_md:
                temps = [x.get("T") for x in osz["ionic"] if x.get("T") is not None]
                if temps:
                    rep["md_T_last"] = temps[-1]
                    tail = temps[len(temps) // 2:]
                    rep["md_T_mean_2nd_half"] = sum(tail) / len(tail)
        if osz["current_elec"]:
            rep["notes"].append(f"{osz['current_elec']} electronic steps after the last ionic step "
                                "(SCF running or interrupted)")
    else:
        rep["n_ionic"] = out["ionic_iter"] if out else 0

    if out:
        rep["finished"] = out["finished"]
        rep["reached"] = out["reached"]
        if out["e0"]:
            rep["E0"] = out["e0"][-1]
        if out["toten"]:
            rep["TOTEN"] = out["toten"][-1]
        rep["efermi"] = out["efermi"]
        if out["ts_last"] is not None and nions:
            rep["TS_per_atom_meV"] = 1000.0 * out["ts_last"] / nions
            if abs(rep["TS_per_atom_meV"]) > 1.0 and rep.get("ismear") not in (-5,):
                rep["problems"].append(f"|T*S| = {abs(rep['TS_per_atom_meV']):.2f} meV/atom > 1 meV: "
                                       "SIGMA too large, reduce it or compare energy(sigma->0)")
        if out["pressure"] is not None:
            rep["pressure_kB"] = out["pressure"]
            rep["pulay_kB"] = out["pulay"]
            if isif >= 3 and relax and abs(out["pressure"]) > 5:
                rep["notes"].append(f"external pressure {out['pressure']:.2f} kB after a variable-cell run: "
                                    "relax again from CONTCAR / raise ENCUT (Pulay stress)")
        if out.get("elapsed") is not None:
            rep["elapsed_s"] = out["elapsed"]
        if out.get("maxmem_kb"):
            rep["maxmem_MB"] = out["maxmem_kb"] / 1024.0
        lr = [x for x in out["loop_real"] if x is not None]
        if lr:
            rep["time_per_ionic_s"] = sum(lr) / len(lr)
            rep["time_ionic_total_s"] = sum(lr)
        # magnetism
        if out["magmom_tot"] is not None and (ip("ISPIN", 1) == 2 or p.get("LSORBIT") == "T"):
            rep["magmom_total"] = out["magmom_tot"]
        sm = out["site_mag"]
        if sm:
            if "y" in sm and "z" in sm:
                mags = [math.sqrt(x * x + y * y + z * z) for x, y, z in zip(sm["x"], sm["y"], sm["z"])]
            else:
                mags = sm["x"]
            order = sorted(range(len(mags)), key=lambda k: -abs(mags[k]))
            rep["site_mag_top"] = [(k + 1, round(mags[k], 3)) for k in order[:args.top]]
        # forces with selective dynamics
        if out["forces"]:
            cons = None
            for nm in ("CONTCAR", "POSCAR"):
                f = d / nm
                if f.is_file() and f.stat().st_size > 0:
                    cons = read_poscar_constraints(f)
                    if cons:
                        break
            forces = out["forces"]
            if cons and cons[1] is not None and cons[2] == len(forces):
                lat, flags, _ = cons
                eff = [project_force(fv, lat, fl) for fv, fl in zip(forces, flags)]
                rep["n_fixed_dof"] = sum(3 - sum(fl) for fl in flags)
            else:
                eff = forces
            norms = [math.sqrt(sum(x * x for x in v)) for v in eff]
            k = max(range(len(norms)), key=lambda i: norms[i])
            rep["fmax"] = norms[k]
            rep["fmax_atom"] = k + 1
            rep["frms"] = math.sqrt(sum(x * x for x in norms) / len(norms))
            if relax and ediffg is not None and ediffg < 0 and rep["fmax"] > abs(ediffg):
                rep["notes"].append(f"max force {rep['fmax']:.4f} eV/A (atom {k + 1}) > |EDIFFG| = {abs(ediffg):g}")
        # boxes
        rep["boxes"] = [{"kind": k[0], "text": k[1], "count": c} for k, c in out["boxes"].items()]
        try:
            age = time.time() - outcar.stat().st_mtime
            rep["outcar_age_s"] = age
        except OSError:
            age = None
    # pattern hits
    hits = []
    for key, sev, _subs, adv in PATTERNS:
        if key in found:
            e = found[key]
            hits.append({"key": key, "severity": sev, "count": e["count"], "line": e["line"],
                         "sources": sorted(e["sources"]), "advice": adv})
    rep["hits"] = hits

    # status
    errors = [h for h in hits if h["severity"] == "error"]
    finished = rep.get("finished", False)
    if relax:
        if rep.get("reached"):
            rep["status"] = "converged"
        elif finished and rep.get("n_ionic", 0) >= nsw:
            rep["status"] = "NSW-exhausted"
            rep["problems"].append(f"NSW={nsw} ionic steps used without reaching EDIFFG: "
                                   "cp CONTCAR POSCAR and continue")
        elif finished:
            rep["status"] = "stopped"
        else:
            rep["status"] = "running/killed"
    elif is_md:
        rep["status"] = "finished" if finished else "running/killed"
    else:
        if finished:
            rep["status"] = "finished" if not rep.get("scf_unconverged_steps") else "SCF-unconverged"
        else:
            rep["status"] = "running/killed"
    if rep["status"] == "running/killed" and out is not None:
        age = rep.get("outcar_age_s")
        if age is not None and age < args.running_age:
            rep["status"] = "running"
        else:
            rep["status"] = "incomplete"
    if errors:
        rep["status"] += "+errors"
    return rep


# ---------------------------------------------------------------------------
# reporting
# ---------------------------------------------------------------------------
def fmt_time(s):
    if s is None:
        return "-"
    s = float(s)
    if s < 120:
        return f"{s:.1f} s"
    if s < 7200:
        return f"{s / 60:.1f} min"
    return f"{s / 3600:.2f} h"


def print_report(rep: dict, args) -> None:
    print("=" * 78)
    print(f"Directory / 目录: {rep['dir']}")
    if rep["status"] == "no-output":
        print("  no OUTCAR/OSZICAR found / 未找到 OUTCAR 或 OSZICAR")
        return
    print(f"  VASP version   : {rep.get('version') or '?'}")
    task = "MD" if rep.get("ibrion") == 0 else ("relaxation" if rep.get("ibrion") in (1, 2, 3) and rep.get("nsw", 0) > 0
                                                 else f"IBRION={rep.get('ibrion')}")
    print(f"  task           : {task}  (NSW={rep.get('nsw')}, ISIF={rep.get('isif')}, NIONS={rep.get('nions')}, "
          f"NKPTS={rep.get('nkpts')}, NBANDS={rep.get('nbands')}, ENCUT={rep.get('encut')})")
    print(f"  status / 状态  : {rep['status']}")
    nel = rep.get("nelec") or []
    print(f"  ionic steps    : {rep.get('n_ionic', 0)}" + (f"   electronic steps: last={nel[-1]}, max={max(nel)}, "
                                                           f"total={sum(nel)} (NELM={rep.get('nelm')})" if nel else ""))
    if nel and args.verbose:
        print("  SCF steps per ionic step: " + " ".join(str(x) for x in nel))
    if rep.get("algos"):
        print("  SCF algorithm  : " + ", ".join(f"{k}x{v}" for k, v in rep["algos"].items()))
    if "E0" in rep:
        print(f"  energy(sigma->0): {rep['E0']:.8f} eV" + (f"   TOTEN: {rep['TOTEN']:.8f} eV" if "TOTEN" in rep else "")
              + (f"   per atom: {rep['E0'] / rep['nions']:.6f} eV" if rep.get("nions") else ""))
    elif rep.get("E0_oszicar") is not None:
        print(f"  E0 (OSZICAR)   : {rep['E0_oszicar']:.8f} eV")
    if "TS_per_atom_meV" in rep:
        print(f"  T*S per atom   : {rep['TS_per_atom_meV']:.3f} meV  (ISMEAR={rep.get('ismear')}, SIGMA={rep.get('sigma')})")
    if rep.get("efermi") is not None:
        print(f"  E-fermi        : {rep['efermi']:.4f} eV")
    if "fmax" in rep:
        extra = f", {rep['n_fixed_dof']} fixed degrees of freedom excluded" if rep.get("n_fixed_dof") else ""
        print(f"  max force      : {rep['fmax']:.4f} eV/A on atom {rep['fmax_atom']} (RMS {rep['frms']:.4f}{extra})"
              + (f"   EDIFFG={rep['ediffg']:g}" if rep.get("ediffg") is not None else ""))
    if "pressure_kB" in rep:
        print(f"  pressure       : {rep['pressure_kB']:.2f} kB (Pulay stress {rep.get('pulay_kB') or 0:.2f} kB)")
    if "magmom_total" in rep:
        m = rep["magmom_total"]
        ms = f"{m:.4f}" if isinstance(m, float) else " ".join(f"{x:.4f}" for x in m)
        print(f"  total magmom   : {ms} mu_B")
    elif "mag_oszicar" in rep:
        m = rep["mag_oszicar"]
        print(f"  total magmom   : {m if isinstance(m, float) else m} mu_B (OSZICAR)")
    if rep.get("site_mag_top"):
        print("  largest |m_i|  : " + ", ".join(f"#{i}: {m}" for i, m in rep["site_mag_top"]))
    if "md_T_last" in rep:
        print(f"  MD temperature : last {rep['md_T_last']:.1f} K, mean of 2nd half {rep['md_T_mean_2nd_half']:.1f} K")
    t = []
    if "elapsed_s" in rep:
        t.append(f"elapsed {fmt_time(rep['elapsed_s'])}")
    if "time_per_ionic_s" in rep:
        t.append(f"{fmt_time(rep['time_per_ionic_s'])}/ionic step")
    if "maxmem_MB" in rep:
        t.append(f"max memory {rep['maxmem_MB']:.0f} MB (root rank)")
    if "outcar_age_s" in rep and not rep.get("finished"):
        t.append(f"OUTCAR last modified {fmt_time(rep['outcar_age_s'])} ago")
    if t:
        print("  timing         : " + "; ".join(t))
    for n in rep["notes"]:
        print(f"  note           : {n}")
    for pmsg in rep["problems"]:
        print(f"  PROBLEM        : {pmsg}")
    boxes = rep.get("boxes") or []
    if boxes:
        print(f"  VASP WARNING/ADVICE messages ({len(boxes)} distinct):")
        for b in boxes[: args.max_boxes]:
            print(f"    [{b['kind']} x{b['count']}] {b['text']}")
        if len(boxes) > args.max_boxes:
            print(f"    ... {len(boxes) - args.max_boxes} more (use --max-boxes)")
    hits = rep.get("hits") or []
    if hits:
        print("  Known messages / 已知报错与建议:")
        for h in hits:
            print(f"    [{h['severity'].upper()}] {h['key']} x{h['count']} in {', '.join(h['sources'])}")
            print(f"        > {h['line']}")
            print(f"        建议: {h['advice']}")
    elif rep["status"] != "no-output":
        print("  no known error strings found / 未发现已知报错")


def brief_line(rep: dict) -> str:
    e0 = rep.get("E0", rep.get("E0_oszicar"))
    e0s = f"{e0:.6f}" if isinstance(e0, float) else "-"
    fm = f"{rep['fmax']:.4f}" if "fmax" in rep else "-"
    mag = rep.get("magmom_total", rep.get("mag_oszicar"))
    if isinstance(mag, list):
        mag = math.sqrt(sum(x * x for x in mag))
    ms = f"{mag:.3f}" if isinstance(mag, float) else "-"
    errs = ",".join(h["key"] for h in rep.get("hits", []) if h["severity"] == "error") or "-"
    return (f"{rep['dir'][:28]:<28} {rep['status'][:22]:<22} {rep.get('n_ionic', 0):>5} "
            f"{e0s:>16} {fm:>9} {ms:>8} {fmt_time(rep.get('elapsed_s')):>10}  {errs}")


def expand_dirs(dirs):
    out = []
    for d in dirs:
        p = Path(d)
        if p.is_file():
            p = p.parent
        if not p.is_dir():
            print(f"WARNING: '{d}' is not a directory, skipped", file=sys.stderr)
            continue
        if not (p / "OUTCAR").is_file() and not (p / "OSZICAR").is_file():
            imgs = sorted(x for x in p.iterdir() if x.is_dir() and re.fullmatch(r"\d{2}", x.name))
            if imgs:
                out += imgs
                continue
        out.append(p)
    return out


def build_parser():
    ap = argparse.ArgumentParser(
        prog="check_vasp.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=("Health check of VASP jobs: SCF/ionic convergence, forces (selective dynamics aware), "
                     "magnetisation, timing and known error messages with fixes.\n"
                     "VASP 作业健康检查：收敛、受力（考虑 Selective dynamics）、磁矩、耗时与常见报错的中文修复建议。"),
        epilog=("Examples / 示例:\n  python check_vasp.py\n  python check_vasp.py */ --brief\n"
                "  python check_vasp.py relax -o relax/vasp.out --json r.json\n"
                "Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License"))
    ap.add_argument("dirs", nargs="*", default=["."], help="job directories (default .); NEB 00..NN detected")
    ap.add_argument("-o", "--stdout", action="append", help="VASP stdout file(s) to scan (default: auto-detect)")
    ap.add_argument("--brief", action="store_true", help="one line per directory / 每个目录一行")
    ap.add_argument("--json", help="write the full report to this JSON file")
    ap.add_argument("--strict", action="store_true", help="exit 1 unless every job converged/finished cleanly")
    ap.add_argument("--top", type=int, default=5, help="number of largest site moments to show (default 5)")
    ap.add_argument("--max-boxes", type=int, default=8, help="max distinct WARNING/ADVICE boxes shown")
    ap.add_argument("--running-age", type=float, default=900,
                    help="OUTCAR modified within this many seconds => 'running' (default 900)")
    ap.add_argument("-v", "--verbose", action="store_true", help="print SCF step count of every ionic step")
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return ap


def main(argv=None) -> int:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except AttributeError:
        pass
    args = build_parser().parse_args(argv)
    dirs = expand_dirs(args.dirs)
    if not dirs:
        print("ERROR: no valid directory given", file=sys.stderr)
        return 2
    reports = []
    for d in dirs:
        try:
            reports.append(analyse(d, args))
        except OSError as exc:
            reports.append({"dir": str(d), "status": "read-error", "problems": [str(exc)], "notes": []})
    if args.brief:
        print(f"{'directory':<28} {'status':<22} {'ionic':>5} {'E0 (eV)':>16} {'Fmax':>9} {'mag':>8} {'elapsed':>10}  errors")
        for r in reports:
            print(brief_line(r))
    else:
        for r in reports:
            print_report(r, args)
    if args.json:
        def conv(o):
            if isinstance(o, set):
                return sorted(o)
            return str(o)
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump(reports, fh, indent=1, ensure_ascii=False, default=conv)
        print(f"JSON report written to {args.json}")
    if all(r["status"] == "no-output" for r in reports):
        return 2
    bad = any("+errors" in r["status"] for r in reports)
    if args.strict:
        good = ("converged", "finished")
        bad = bad or any(r["status"] not in good for r in reports)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
