#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
opt_monitor.py - follow a Gaussian 16 geometry optimization step by step.

逐步解析 Gaussian 16 几何优化（含 TS 优化、柔性扫描）每一步的能量与四项收敛判据
（Maximum/RMS Force、Maximum/RMS Displacement），打印表格、给出“振荡 / 缓慢下降 / 平坦势能面”
的诊断与建议，并画出各判据与阈值的对比图——L9999（步数用完）时用它判断下一步该怎么做。

Usage examples
--------------
  python opt_monitor.py opt.log                    # table + diagnosis + opt_opt.png
  python opt_monitor.py opt.log -o conv.png --dpi 200
  python opt_monitor.py opt.log --no-plot --csv steps.csv
  python opt_monitor.py scan.log --scan-point 7    # only one point of a relaxed scan
  python opt_monitor.py ts.log --last 30 --show    # interactive window instead of a file
  python opt_monitor.py running.log                # works on unfinished jobs, too

Notes
-----
* One record per Berny step = one "Converged?" table.  The energy of a step is the energy of the
  geometry whose forces were tested (Counterpoise > ONIOM extrapolated > MP2/double hybrid >
  TD total > SCF Done).
* Only job steps whose route contains Opt are analysed.  The extra convergence table that a Freq job
  prints after an Opt (re-test with the exact Hessian) is reported separately, not as a step.
* Gaussian declares convergence when all four criteria are YES, or when the forces are 100 times
  smaller than their thresholds even if the displacements are not ("negligible forces").
* Plot: energy relative to the lowest step (kcal/mol) plus four log-scale panels with dashed
  threshold lines; filled markers = criterion met, open markers = not met.

Dependencies: Python >= 3.9; matplotlib (only for the figure; --no-plot works without it).

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License
"""

import argparse
import csv
import math
import os
import re
import sys

HARTREE_TO_KCAL = 627.509474

RE_DASH = re.compile(r"^ -{10,}\s*$")
RE_STEP = re.compile(r"Step number\s+(\d+)\s+out of a maximum of\s+(\d+)(?:\s+on scan point\s+(\d+)\s+out of\s+(\d+))?")
RE_IRC_STEP = re.compile(r"^\s*Pt\s+\d+\s+Step number")
RE_ITEM = re.compile(r"^\s*(Maximum Force|RMS\s+Force|Maximum Displacement|RMS\s+Displacement)"
                     r"\s+(\S+)\s+(\S+)\s+(YES|NO)\s*$")
RE_PRED = re.compile(r"Predicted change in Energy\s*=\s*(\S+)")
RE_SCF = re.compile(r"SCF Done:\s+E\(([^)]+)\)\s*=\s*(-?\d+\.\d+)")
RE_MP2 = re.compile(r"EUMP2\s*=\s*(-?\d+\.\d+[DdEe][+-]\d+)")
RE_DH = re.compile(r"E2\(([^)]+)\)\s*=\s*\S+\s+E\(([^)]+)\)\s*=\s*(-?\d+\.\d+[DdEe][+-]\d+)")
RE_ONIOM = re.compile(r"ONIOM: extrapolated energy\s*=\s*(-?\d+\.\d+)")
RE_CP = re.compile(r"Counterpoise corrected energy\s*=\s*(-?\d+\.\d+)")
RE_TD = re.compile(r"Total Energy, E\(([^)]*)\)\s*=\s*(-?\d+\.\d+)")
RE_OPT_KW = re.compile(r"(^|[\s#/,])opt(?=$|[\s=(,])", re.IGNORECASE)
RE_IRC_KW = re.compile(r"(^|[\s#/,])irc(?=$|[\s=(,])", re.IGNORECASE)

ITEM_KEYS = {"Maximum Force": "max_force", "RMS Force": "rms_force",
             "Maximum Displacement": "max_disp", "RMS Displacement": "rms_disp"}
CRITERIA = ("max_force", "rms_force", "max_disp", "rms_disp")
CRIT_TITLE = {"max_force": "Maximum Force", "rms_force": "RMS Force",
              "max_disp": "Maximum Displacement", "rms_disp": "RMS Displacement"}
ENERGY_PRIORITY = ("cp", "oniom", "mp2", "td", "scf")

# chart styling (single series per panel; recessive grid; threshold as a dashed reference line)
COLOR_SERIES = "#2a78d6"
COLOR_REF = "#52514e"
COLOR_GRID = "#e4e3df"
COLOR_TEXT = "#0b0b0b"
COLOR_MUTED = "#8a8984"
COLOR_SURFACE = "#fcfcfb"


def dfloat(text):
    return float(text.replace("D", "E").replace("d", "e"))


def to_float(text):
    try:
        return dfloat(text)
    except ValueError:
        return float("nan")      # '********' overflow fields


# --------------------------------------------------------------------------------------
# parsing
# --------------------------------------------------------------------------------------
def parse_log(path):
    """Return (steps, info). Each step is a dict with energy, the four criteria and flags."""
    steps = []
    info = {"normal": 0, "error_link": None, "freq_check": None, "stopped": None, "route": None,
            "negligible": False}
    job = 1
    job_kind = "other"
    prev = ""
    capturing = False
    route_buf = []
    energies = {}
    cur = None              # step being filled (after 'Step number')
    table = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if capturing:
                if RE_DASH.match(line):
                    capturing = False
                    route = "".join(route_buf).strip()
                    if RE_IRC_KW.search(route):
                        job_kind = "irc"
                    elif RE_OPT_KW.search(route):
                        job_kind = "opt"
                        info["route"] = info["route"] or route
                    else:
                        job_kind = "other"
                else:
                    route_buf.append(line[1:] if line.startswith(" ") else line)
                prev = line
                continue
            if line.startswith(" #") and RE_DASH.match(prev):
                capturing = True
                route_buf = [line[1:]]
                prev = line
                continue
            prev = line

            if "orientation:" in line and ("Input" in line or "Standard" in line or "Z-Matrix" in line):
                if energies:            # a new geometry starts: forget the energies of the previous one
                    energies = {}
                continue
            m = RE_SCF.search(line)
            if m:
                energies["scf"] = float(m.group(2))
                continue
            if "EUMP2" in line:
                m = RE_MP2.search(line)
                if m:
                    energies["mp2"] = dfloat(m.group(1))
                continue
            if "E2(" in line:
                m = RE_DH.search(line)
                if m:
                    energies["mp2"] = dfloat(m.group(3))
                continue
            if "ONIOM: extrapolated energy" in line:
                m = RE_ONIOM.search(line)
                if m:
                    energies["oniom"] = float(m.group(1))
                continue
            if "Counterpoise corrected energy" in line:
                m = RE_CP.search(line)
                if m:
                    energies["cp"] = float(m.group(1))
                continue
            if "Total Energy, E(" in line:
                m = RE_TD.search(line)
                if m:
                    energies["td"] = float(m.group(2))
                continue
            if "Proceeding to internal job step number" in line:
                job += 1
                job_kind = "other"
                continue
            if "Normal termination of Gaussian" in line:
                info["normal"] += 1
                continue
            if "Error termination via Lnk1e" in line:
                mm = re.search(r"l(\d+)\.exe", line, re.IGNORECASE)
                info["error_link"] = "L" + mm.group(1) if mm else "L?"
                continue

            m = RE_STEP.search(line)
            if m and not RE_IRC_STEP.match(line):
                if job_kind != "opt":
                    if job_kind == "other":
                        cur = {"freq_check": True}
                    continue
                e = None
                etype = None
                for key in ENERGY_PRIORITY:
                    if key in energies:
                        e, etype = energies[key], key.upper()
                        break
                cur = {"job": job, "step": int(m.group(1)), "maxcycles": int(m.group(2)),
                       "scan_point": int(m.group(3)) if m.group(3) else None,
                       "scan_total": int(m.group(4)) if m.group(4) else None,
                       "energy": e, "energy_type": etype, "converged": False, "predicted": None}
                table = {}
                continue
            m = RE_ITEM.match(line)
            if m and cur is not None:
                key = ITEM_KEYS[re.sub(r"\s+", " ", m.group(1))]
                table[key] = (to_float(m.group(2)), to_float(m.group(3)), m.group(4) == "YES")
                if len(table) == 4:
                    if cur.get("freq_check"):
                        info["freq_check"] = "all YES" if all(v[2] for v in table.values()) else "NOT converged"
                        cur = None
                    else:
                        for k, (val, thr, ok) in table.items():
                            cur[k] = val
                            cur[k + "_thr"] = thr
                            cur[k + "_ok"] = ok
                        steps.append(cur)
                continue
            m = RE_PRED.search(line)
            if m and steps and cur is steps[-1]:
                cur["predicted"] = to_float(m.group(1))
                continue
            if "Optimization completed" in line:
                if "negligible forces" in line:
                    info["negligible"] = True
                if job_kind == "opt" and steps:
                    steps[-1]["converged"] = True
                elif job_kind == "other":
                    info["freq_check"] = "converged"
                continue
            if "Optimization stopped" in line:
                info["stopped"] = "Optimization stopped"
                continue
            if "Number of steps exceeded" in line:
                info["stopped"] = "Number of steps exceeded (L9999)"
                continue
            if "Wrong number of Negative eigenvalues" in line:
                info["stopped"] = line.strip()
                continue
    return steps, info


# --------------------------------------------------------------------------------------
# diagnosis
# --------------------------------------------------------------------------------------
def diagnose(steps, info, is_ts=False):
    """Return a list of (Chinese) diagnosis lines for the last optimization in `steps`."""
    out = []
    if not steps:
        return ["没有找到优化步（Converged? 表）。这不是 Opt 作业，或作业在第一步之前就中断了。"]
    last = steps[-1]
    n = len(steps)
    if last["converged"]:
        out.append("优化已收敛（Stationary point found），共 %d 步。" % n)
        if info.get("negligible"):
            out.append("注意：是按“力可以忽略”判据收敛的（位移未达标）。势能面可能很平坦，"
                       "频率计算中留意低频/小虚频。")
        if info.get("freq_check") == "NOT converged":
            out.append("但后续 Freq 步用精确 Hessian 重新检验时判为未收敛：结构不是严格的驻点，"
                       "建议 Opt=(ReadFC 或 CalcFC) 从该结构续算后再做 Freq。")
        return out

    e = [s["energy"] for s in steps if s["energy"] is not None]
    w = min(10, len(e) - 1)
    status = info.get("stopped") or ("作业以 %s 报错结束" % info["error_link"] if info.get("error_link")
                                     else "作业未结束（仍在运行、被杀死或日志被截断）")
    out.append("优化未收敛：%s；当前共 %d 步（MaxCycles=%s）。" % (status, n, last.get("maxcycles")))
    oks = ["YES" if last.get(k + "_ok") else "NO" for k in CRITERIA]
    out.append("最后一步判据：MaxF %s, RMSF %s, MaxD %s, RMSD %s。" % tuple(oks))
    if w < 3:
        out.append("步数太少，无法判断趋势。")
        return out
    tail = e[-(w + 1):]
    d = [b - a for a, b in zip(tail, tail[1:])]
    ups = sum(1 for x in d if x > 1e-6)
    sign_changes = sum(1 for a, b in zip(d, d[1:]) if a * b < 0 and abs(a) > 1e-7 and abs(b) > 1e-7)
    drop = (tail[0] - tail[-1]) * HARTREE_TO_KCAL
    spread = (max(tail) - min(tail)) * HARTREE_TO_KCAL
    out.append("最近 %d 步：能量共下降 %.3f kcal/mol，波动范围 %.3f kcal/mol，%d 次上升，%d 次升降交替。"
               % (w, drop, spread, ups, sign_changes))
    forces_ok = last.get("max_force_ok") and last.get("rms_force_ok")
    if is_ts:
        out.append("这是 TS 优化：能量先升后降属正常；关键是 Hessian 是否始终只有 1 个负本征值。")
    if sign_changes >= 3 and ups >= 2:
        out.append("诊断：振荡。能量和力在几个结构之间来回跳动。")
        out.append("建议：从能量最低的一步重新开始（log2xyz.py --lowest --gjf），并用 Opt=CalcFC（或 ReadFC/"
                   "Recalc=N）获得准确 Hessian，减小步长 Opt=MaxStep=10，必要时 Int=SuperFine；"
                   "同时确认每一步 SCF 都已收敛（无 'Convergence criterion not met'）。")
    elif drop > 0.3 and ups <= 2:
        out.append("诊断：缓慢但稳定地下降，只是步数不够。")
        out.append("建议：从最后一步续算并增大 MaxCycles（Opt=(Restart,MaxCycles=200) 需要 chk；或 "
                   "log2xyz.py --gjf --maxcycles 200），初始结构很差时先用小基组或 xtb 预优化。")
    elif spread < 0.05:
        if forces_ok:
            out.append("诊断：平坦势能面——力已达标，位移迟迟不达标（如甲基/苯环转动等软模式）。")
            out.append("建议：可用 Opt=CalcFC 或 CalcAll 获得准确曲率；结构多半已足够好，可直接做 Freq 检查"
                       "低频模式，也可以从最后结构以更紧的积分格点重新优化。")
        else:
            out.append("诊断：能量几乎不变但力仍未达标，可能卡在数值噪声或坐标系问题上。")
            out.append("建议：Int=SuperFine、SCF=Tight，Opt=CalcFC；若有接近 180° 的键角，改用 Opt=Cartesian。")
    else:
        out.append("诊断：没有明显规律（可能在势能面上漫游或结构发生了较大重排）。")
        out.append("建议：用 GaussView 或 log2xyz.py --all 查看轨迹，确认结构、电荷与多重度合理，"
                   "再从合理的结构用 Opt=CalcFC 重新优化。")
    return out


# --------------------------------------------------------------------------------------
# output
# --------------------------------------------------------------------------------------
def print_table(steps):
    e_valid = [s["energy"] for s in steps if s["energy"] is not None]
    emin = min(e_valid) if e_valid else None
    print("%5s %5s %18s %10s %9s %10s %10s %10s %10s  %s" %
          ("#", "scan", "energy/Eh", "E-Emin", "dE", "MaxForce", "RMSForce", "MaxDisp", "RMSDisp", "YES?"))
    print("%5s %5s %18s %10s %9s %10s %10s %10s %10s" % ("", "", "", "kcal/mol", "kcal/mol", "", "", "", ""))
    prev = None
    for i, s in enumerate(steps, 1):
        e = s["energy"]
        de = "%9.4f" % ((e - prev) * HARTREE_TO_KCAL) if (e is not None and prev is not None) else "%9s" % "-"
        flags = "".join("Y" if s.get(k + "_ok") else "n" for k in CRITERIA)
        if s["converged"]:
            flags += "  <- converged"
        print("%5d %5s %18s %10s %s %10.6f %10.6f %10.6f %10.6f  %s" % (
            i, s["scan_point"] or "-", "%.8f" % e if e is not None else "-",
            "%.4f" % ((e - emin) * HARTREE_TO_KCAL) if e is not None else "-", de,
            s["max_force"], s["rms_force"], s["max_disp"], s["rms_disp"], flags))
        prev = e if e is not None else prev
    if steps:
        s = steps[-1]
        print("thresholds: MaxF %.6f  RMSF %.6f  MaxD %.6f  RMSD %.6f" %
              tuple(s[k + "_thr"] for k in CRITERIA))


def write_csv(steps, path):
    keys = ["index", "job", "step", "scan_point", "energy", "energy_type"]
    for k in CRITERIA:
        keys += [k, k + "_thr", k + "_ok"]
    keys += ["predicted", "converged"]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for i, s in enumerate(steps, 1):
            row = dict(s)
            row["index"] = i
            w.writerow({k: ("" if row.get(k) is None else row.get(k)) for k in keys})


def plot_steps(steps, info, title, out_path, show=False, dpi=150):
    try:
        import matplotlib
        if not show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec
    except ImportError:
        print("ERROR: matplotlib is not installed (pip install matplotlib) - use --no-plot", file=sys.stderr)
        return False

    x = list(range(1, len(steps) + 1))
    e = [s["energy"] for s in steps]
    valid = [v for v in e if v is not None]
    emin = min(valid) if valid else 0.0

    plt.rcParams.update({"font.size": 9, "axes.edgecolor": COLOR_MUTED, "axes.labelcolor": COLOR_TEXT,
                         "xtick.color": COLOR_MUTED, "ytick.color": COLOR_MUTED})
    fig = plt.figure(figsize=(9.0, 8.0), facecolor=COLOR_SURFACE)
    gs = GridSpec(3, 2, figure=fig, height_ratios=[1.15, 1, 1], hspace=0.55, wspace=0.28)
    axes = [fig.add_subplot(gs[0, :])] + [fig.add_subplot(gs[1 + i // 2, i % 2]) for i in range(4)]

    # boundaries between scan points / opt job steps
    bounds = []
    for i in range(1, len(steps)):
        a, b = steps[i - 1], steps[i]
        if a["job"] != b["job"] or a["scan_point"] != b["scan_point"]:
            bounds.append(i + 0.5)

    def style(ax):
        ax.set_facecolor(COLOR_SURFACE)
        ax.grid(True, color=COLOR_GRID, linewidth=0.8, linestyle="-")
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for b in bounds:
            ax.axvline(b, color=COLOR_GRID, linewidth=1.0, zorder=0)
        ax.set_xlim(0.5, len(steps) + 0.5)

    # energy panel
    ax = axes[0]
    style(ax)
    xe = [xi for xi, v in zip(x, e) if v is not None]
    ye = [(v - emin) * HARTREE_TO_KCAL for v in e if v is not None]
    ax.plot(xe, ye, color=COLOR_SERIES, linewidth=2, marker="o", markersize=4, zorder=3)
    conv = [(xi, (s["energy"] - emin) * HARTREE_TO_KCAL) for xi, s in zip(x, steps)
            if s["converged"] and s["energy"] is not None]
    if conv:
        ax.plot([c[0] for c in conv], [c[1] for c in conv], linestyle="none", marker="o", markersize=9,
                markerfacecolor=COLOR_SURFACE, markeredgecolor=COLOR_SERIES, markeredgewidth=2, zorder=4)
    ax.set_title("Energy relative to lowest step (kcal/mol)" + ("   (rings = converged)" if conv else ""),
                 loc="left", color=COLOR_TEXT, fontsize=10)
    ax.set_xlabel("optimization step")

    # criteria panels
    for ax, key in zip(axes[1:], CRITERIA):
        style(ax)
        vals = [s[key] for s in steps]
        thr = steps[-1][key + "_thr"]
        pos = [(xi, v) for xi, v in zip(x, vals) if v is not None and not math.isnan(v) and v > 0]
        if pos:
            ax.plot([p[0] for p in pos], [p[1] for p in pos], color=COLOR_SERIES, linewidth=2, zorder=3)
            met = [(xi, s[key]) for xi, s in zip(x, steps) if s.get(key + "_ok") and s[key] and s[key] > 0]
            unmet = [(xi, s[key]) for xi, s in zip(x, steps) if not s.get(key + "_ok") and s[key] and s[key] > 0]
            if met:
                ax.plot([p[0] for p in met], [p[1] for p in met], linestyle="none", marker="o", markersize=4.5,
                        color=COLOR_SERIES, zorder=4)
            if unmet:
                ax.plot([p[0] for p in unmet], [p[1] for p in unmet], linestyle="none", marker="o",
                        markersize=4.5, markerfacecolor=COLOR_SURFACE, markeredgecolor=COLOR_SERIES, zorder=4)
            ax.set_yscale("log")
        ax.axhline(thr, color=COLOR_REF, linestyle=(0, (4, 3)), linewidth=1.2, zorder=2)
        ax.annotate("threshold %.6f" % thr, xy=(1.0, thr), xycoords=("axes fraction", "data"),
                    xytext=(-2, 3), textcoords="offset points", ha="right", va="bottom",
                    fontsize=8, color=COLOR_REF)
        ax.set_title(CRIT_TITLE[key] + " (a.u.)", loc="left", color=COLOR_TEXT, fontsize=10)
        ax.set_xlabel("optimization step")

    status = ("converged" if steps[-1]["converged"]
              else info.get("stopped") or ("error " + info["error_link"] if info.get("error_link") else "not converged"))
    fig.suptitle("%s  -  %d steps, %s" % (title, len(steps), status), x=0.06, ha="left",
                 color=COLOR_TEXT, fontsize=11)
    if show:
        plt.show()
    else:
        fig.savefig(out_path, dpi=dpi, facecolor=COLOR_SURFACE, bbox_inches="tight")
        print("figure written to %s" % out_path, file=sys.stderr)
    plt.close(fig)
    return True


def build_parser():
    p = argparse.ArgumentParser(
        description="Parse every optimization step of a Gaussian 16 log (energy + 4 convergence criteria), "
                    "diagnose oscillation vs slow descent, and plot against the thresholds.  "
                    "解析每个优化步的能量与四项收敛判据，诊断振荡/缓慢下降，并画图。",
        epilog="Examples:\n"
               "  python opt_monitor.py opt.log\n"
               "  python opt_monitor.py opt.log -o conv.png --csv steps.csv\n"
               "  python opt_monitor.py scan.log --scan-point 7\n"
               "  python opt_monitor.py ts.log --last 30 --show",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("logs", nargs="+", help="Gaussian .log/.out file(s) / 输出文件")
    p.add_argument("-o", "--output", help="figure file (single log; default <log stem>_opt.png)")
    p.add_argument("--no-plot", action="store_true", help="no figure (matplotlib not needed) / 不画图")
    p.add_argument("--show", action="store_true", help="open an interactive window instead of saving")
    p.add_argument("--dpi", type=int, default=150, help="figure resolution (default 150)")
    p.add_argument("--csv", metavar="FILE", help="write the step table to CSV (single log)")
    p.add_argument("--scan-point", type=int, metavar="K", help="only steps of relaxed-scan point K")
    p.add_argument("--job", type=int, metavar="N", help="only optimization steps of job step N (Link1)")
    p.add_argument("--last", type=int, metavar="N", help="only the last N steps")
    p.add_argument("-q", "--quiet", action="store_true", help="no table, only the diagnosis")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    if len(args.logs) > 1 and (args.output or args.csv):
        print("ERROR: -o/--csv can only be used with a single log", file=sys.stderr)
        return 2
    rc = 0
    for log in args.logs:
        if not os.path.isfile(log):
            print("ERROR: %s not found" % log, file=sys.stderr)
            rc = 1
            continue
        try:
            steps, info = parse_log(log)
        except OSError as exc:
            print("ERROR: cannot read %s: %s" % (log, exc), file=sys.stderr)
            rc = 1
            continue
        if args.job is not None:
            steps = [s for s in steps if s["job"] == args.job]
        if args.scan_point is not None:
            steps = [s for s in steps if s["scan_point"] == args.scan_point]
        if args.last:
            steps = steps[-args.last:]
        print("=" * 78)
        print("%s" % log)
        if info["route"]:
            print("route: %s" % info["route"])
        if not steps:
            print("no optimization steps found (not an Opt job, wrong --job/--scan-point, or stopped before step 1)")
            rc = 1
            continue
        scan_total = steps[-1].get("scan_total")
        if scan_total:
            pts = sorted({s["scan_point"] for s in steps})
            nconv = sum(1 for s in steps if s["converged"])
            print("relaxed scan: points %d-%d of %d present, %d converged" % (pts[0], pts[-1], scan_total, nconv))
        if not args.quiet:
            print_table(steps)
        # diagnose the last optimization (last scan point / last job)
        lastgrp = [s for s in steps if s["job"] == steps[-1]["job"] and s["scan_point"] == steps[-1]["scan_point"]]
        is_ts = bool(info["route"] and re.search(r"(?i)\bts\b|qst[23]", info["route"]))
        print("-" * 78)
        for line in diagnose(lastgrp, info, is_ts):
            print(line)
        if args.csv:
            write_csv(steps, args.csv)
            print("CSV written to %s" % args.csv, file=sys.stderr)
        if not args.no_plot:
            out = args.output or os.path.splitext(log)[0] + "_opt.png"
            if not plot_steps(steps, info, os.path.basename(log), out, show=args.show, dpi=args.dpi):
                rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main())
