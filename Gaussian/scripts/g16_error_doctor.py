#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
g16_error_doctor.py - diagnose failed (or suspicious) Gaussian 16 jobs from their output files.

读取 Gaussian 16 的 .log/.out：判断作业是正常结束、报错结束还是被中断；从
"Error termination via Lnk1e in .../lXXX.exe" 识别出错的 Link，在出错的作业步中匹配常见报错原文
（L502 Convergence failure、L9999 Number of steps exceeded、L123 corrector steps、L301 多重度/基组/色散、
End of file in ZSymb、galloc、Erroneous write、Wrong number of Negative eigenvalues 等），输出中文的
原因解释与处理建议；正常结束的作业也会检查虚频个数、SCF 未收敛警告、<S**2> 自旋污染和波函数不稳定。
可选 --draft：调用同目录的 log2xyz.py 起草一个续算 .gjf（只是草稿，提交前必须人工检查）。

Usage examples
--------------
  python g16_error_doctor.py job.log                 # diagnosis in Chinese
  python g16_error_doctor.py *.log -q                # one status line per file
  python g16_error_doctor.py runs/ -r --json report.json
  python g16_error_doctor.py job.log --context 30    # show 30 lines above the error termination
  python g16_error_doctor.py job.log --draft         # also write job_restart.gjf for review
  python g16_error_doctor.py job.log --draft fix.gjf --force
  python g16_error_doctor.py --list-errors           # table of all recognised messages

Notes
-----
* The failing Link is taken from "Error termination via Lnk1e in <path>/lNNN.exe" (Linux) or
  "...\\lNNN.exe" (Windows); "Error termination request processed by link NNN" is used as fallback.
* Error strings are matched only inside the job step that failed (from its route section to the
  error termination), so harmless warnings of earlier --Link1-- steps are not blamed.  A few messages
  that also occur as warnings in healthy runs (Linear angle in Bend, Bend failed ...) are only used
  when they appear in the last 80 lines before the error.
* Message texts follow Gaussian 16 Rev. C.01/C.02; wording can differ slightly between revisions,
  so always read the lines above the error yourself (--context).
* --draft never edits your input: it writes a new file (refuses to overwrite without --force),
  built from the last geometry in the log via log2xyz.py (must sit next to this script).
  Gen/GenECP basis blocks are not in the log and must be pasted by hand.

Dependencies: Python >= 3.9 standard library; log2xyz.py from the same folder (only for --draft).

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License
"""

import argparse
import collections
import glob
import json
import os
import re
import sys
import time

RE_DASH = re.compile(r"^ -{10,}\s*$")
RE_LNK1E = re.compile(r"Error termination via Lnk1e in\s+(\S*?)[\\/]?l(\d+)\.exe", re.IGNORECASE)
RE_LINKREQ = re.compile(r"Error termination request processed by link\s+(\d+)", re.IGNORECASE)
RE_CHARGE = re.compile(r"Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)")
RE_STEPMAX = re.compile(r"Step number\s+\d+\s+out of a maximum of\s+(\d+)")
RE_S2 = re.compile(r"S\*\*2 before annihilation\s+(-?\d+\.\d+),\s+after\s+(-?\d+\.\d+)")
RE_FREQ = re.compile(r"^\s*Frequencies --\s+(.*)$")
RE_INSTAB = re.compile(r"The wavefunction has an?\s+(.*?instability)", re.IGNORECASE)
RE_OPT_KW = re.compile(r"(^|[\s#/,])opt(?=$|[\s=(,])", re.IGNORECASE)
RE_TS_KW = re.compile(r"opt\s*=?\s*\(?[^\s]*\b(ts|qst2|qst3)\b", re.IGNORECASE)

# ---------------------------------------------------------------------------------------------
# what each Link does (Chinese), used when no known message is found
# ---------------------------------------------------------------------------------------------
LINK_INFO = {
    1: "L1：解析 route。通常是 route 语法错误、关键词拼写错误或关键词互相冲突。",
    101: "L101：读取标题、电荷多重度和坐标。常见于空行缺失、缺电荷多重度行、坐标格式错误。",
    103: "L103：Berny 优化器（内坐标、收敛判断）。常见于线性角导致内坐标失效、TS 优化 Hessian 曲率不对。",
    123: "L123：IRC 积分。通常是 corrector 步失败。",
    202: "L202：坐标重定向、对称性与距离矩阵。常见于原子重叠、Bohr/Angstrom 单位弄错。",
    301: "L301：建立基组、ECP 与色散参数。常见于多重度与电子数矛盾、基组不支持该元素、Gen/GenECP 格式错误。",
    401: "L401：初始猜测。常见于 Guess=Read 读不到 chk，或 chk 与当前体系（原子数、基组）不匹配。",
    502: "L502：SCF 迭代。通常是 SCF 不收敛或 DFT 数值积分不准。",
    508: "L508：二次收敛 SCF（SCF=QC/XQC）。QC 本身也没有收敛。",
    801: "L801：积分变换（冻芯设置）。后 HF、双杂化、TD 计算中的冻芯问题。",
    914: "L914：CIS/TD-DFT 激发态迭代。",
    1002: "L1002：CPHF（频率、NMR、极化率）。CPHF 不收敛或数值问题。",
    9999: "L9999：结束作业、写 archive。通常是优化/扫描在给定步数内没有收敛。",
}

# ---------------------------------------------------------------------------------------------
# known messages: id, regex, typical Link(s), scope ("step" = anywhere in failing step,
# "near" = last 80 lines before the error), title, cause, fixes, draft strategy
# ---------------------------------------------------------------------------------------------
KNOWN = [
    dict(id="qperr", rx=r"QPErr --- A syntax error was detected", links=(1,), scope="step",
         title="route 语法错误（QPErr）",
         cause="route 中括号不配对、关键词或选项拼写错误、选项间缺逗号、使用了全角字符，或 route 没有以 # 开头。",
         fixes=["看日志中 QPErr 上方的标记，它指出解析失败的位置；",
                "用 testrt \"#p ...\" 只检查 route 语法；对照 https://gaussian.com/keywords/ 核对拼写",
                "（def2SVP、wB97XD、M062X 等写法不带连字符）。"],
         draft=None),
    dict(id="itype", rx=r"Illegal IType or MSType generated by parse", links=(1,), scope="step",
         title="关键词组合不合法（Illegal IType or MSType）",
         cause="route 中的关键词组合不被支持（例如对某种方法请求了它不支持的计算）。",
         fixes=["简化 route，只保留方法/基组和任务，再逐个加回其他关键词定位冲突项。"], draft=None),
    dict(id="zsymb", rx=r"End of file in ZSymb", links=(101,), scope="step",
         title="读坐标时遇到文件结尾（End of file in ZSymb）",
         cause="坐标后面没有空行（最常见），或标题/电荷多重度段缺失；--Link1-- 后的步骤缺标题或电荷多重度。",
         fixes=["文件末尾至少留一个空行（cat -A job.gjf | tail 查看）；",
                "按“route / 空行 / 标题 / 空行 / 电荷 多重度 / 坐标 / 空行”检查分段；",
                "--Link1-- 步骤用 Geom=AllCheck 时不要再写标题和电荷多重度。"], draft=None),
    dict(id="mult", rx=r"The combination of multiplicity\s+\d+\s+and\s+\d+\s+electrons is impossible",
         links=(301,), scope="step",
         title="多重度与电子数矛盾",
         cause="电子数为奇数时多重度只能是 2、4、6……，偶数时只能是 1、3、5……；电荷或多重度写错。",
         fixes=["电子数 = 原子序数之和 - 总电荷（用 ECP 时也按全电子计数）；",
                "修正电荷/多重度；开壳层用 UB3LYP 等非限制性方法；gjf_builder.py 生成输入时会自动检查。"],
         draft=None),
    dict(id="basis_range", rx=r"Atomic number out of range for .*basis set", links=(301,), scope="step",
         title="基组不支持该元素（Atomic number out of range）",
         cause="所选基组（如 6-31G 系列）不包含体系中的某些元素（过渡金属、I 等重元素）。",
         fixes=["全体系改用 def2SVP/def2TZVP（重元素自动带 def2-ECP）；",
                "或用 GenECP：重元素 LANL2DZ/SDD，其余元素用 Pople 基组。"], draft=None),
    dict(id="r6ds8", rx=r"R6DS8: Unable to choose the S8 parameter", links=(301,), scope="step",
         title="该泛函没有所选的色散参数（R6DS8）",
         cause="EmpiricalDispersion=GD3BJ 只对内置了 BJ 参数的泛函有效。",
         fixes=["M06-2X 等 Minnesota 泛函改用 EmpiricalDispersion=GD3；",
                "wB97XD、B97D3、B2PLYPD3、APFD 等已含色散，不要再加 EmpiricalDispersion。"], draft=None),
    dict(id="gen_basis", rx=r"End of file reading basis center", links=(301,), scope="step",
         title="Gen 基组段格式错误",
         cause="基组段没有正常结束：元素块缺 ****，或缺结尾空行。",
         fixes=["每个元素块以 **** 结束，分子中每种元素只列一次；基组段后空一行，文件末尾留空行。"],
         draft=None),
    dict(id="ecp_card", rx=r"EOF while reading ECP pointer card", links=(301,), scope="step",
         title="ECP 段缺失或格式错误",
         cause="route 用了 GenECP（或 Pseudo=Read），但 ECP 段缺失或缺结尾空行。",
         fixes=["在基组段之后空一行写 ECP 段（只列需要 ECP 的元素），整个文件以空行结束。"],
         draft=None),
    dict(id="distance", rx=r"Small interatomic distances encountered|Problem with the distance matrix",
         links=(202,), scope="step",
         title="原子距离过近或重叠",
         cause="复制坐标时重复了原子、Bohr 与 Angstrom 单位混用、Z-matrix 变量写错、截取团簇时有周期性镜像原子。",
         fixes=["在 GaussView 中打开输入检查重叠原子（报错行列出了原子编号）；Bohr 坐标需加 Units=Bohr。"],
         draft=None),
    dict(id="intcoord", rx=r"FormBX had a problem|Error in internal coordinate system", links=(103,),
         scope="step",
         title="冗余内坐标失效（FormBX / internal coordinate system）",
         cause="某个键角接近 180°（炔、腈、M-C≡O、X···H···Y 型 TS），内坐标定义奇异。",
         fixes=["先从最后的结构重新开始优化（Gaussian 会按新结构重建内坐标）；",
                "仍失败时改用 Opt=Cartesian；或在线性片段旁加哑原子。"],
         draft="same"),
    dict(id="linear", rx=r"Linear angle in (Bend|Tors)|Bend failed for angle|Tors failed for dihedral",
         links=(103,), scope="near",
         title="线性角导致内坐标失效（Linear angle in Bend/Tors）",
         cause="优化过程中出现接近 180° 的键角，键角/二面角坐标失效。",
         fixes=["从最后的结构重启（Geom=Check Guess=Read Opt，或用 --draft 生成的输入）；",
                "仍失败时用 Opt=Cartesian；扫描时避免让坐标跨过线性构型。"],
         draft="same"),
    dict(id="negeig", rx=r"Wrong number of Negative eigenvalues", links=(103,), scope="step",
         title="TS 优化 Hessian 负本征值个数不对",
         cause="TS 优化要求 Hessian 恰好有 1 个负本征值；估计的初始 Hessian 曲率不对或初猜离 TS 太远。",
         fixes=["Opt=(TS,CalcFC,NoEigenTest)：计算准确的初始 Hessian，并允许负本征值个数暂时不对；",
                "改进初猜：柔性扫描的最高点（scan_irc_profile.py 会报告）或 QST2/QST3；",
                "最后必须用频率确认恰好 1 个虚频，并做 IRC。"],
         draft="negeig"),
    dict(id="corrector", rx=r"Maximum number of corrector steps exce(e)?ded", links=(123,), scope="step",
         title="IRC corrector 步失败",
         cause="IRC 的 corrector 迭代没有收敛（原文拼写为 exceded）；常见于 TS 不准、力常数不对或路径走到很平的区域。",
         fixes=["确认 TS 收敛且恰好 1 个虚频，IRC 用 RCFC 读同一水平的力常数；",
                "加大 IRC=(MaxCycle=100)；改用 IRC=(RCFC,LQA) 或 IRC=(RCFC,GS2,StepSize=5)；",
                "可用 IRC=(CalcFC,ReCalc=5) 定期重算 Hessian；已接近端点时直接取最后一点做 Opt 确认连接关系。"],
         draft=None),
    dict(id="scf", rx=r"Convergence failure -- run terminated", links=(502,), scope="step",
         title="SCF 不收敛（L502 Convergence failure）",
         cause="SCF 在最大迭代次数内没有收敛。常见原因：初猜差、HOMO-LUMO 能隙小、电荷/多重度或结构错误。",
         fixes=["先用 grep -E \"^ (Cycle|E=)\" job.log | tail -n 40 看能量是振荡还是缓慢下降；",
                "先检查电荷、多重度、结构和基组/ECP；振荡时加 MaxCycle 没用；",
                "依次尝试 SCF=XQC → 小基组/闭壳层离子收敛后 Guess=Read → SCF=VShift=500 → SCF=Fermi；",
                "开壳层过渡金属可试 Guess=Mix、Guess=Fragment、SCF=NoIncFock；收敛后用 Stable=Opt 检查。"],
         draft="scf"),
    dict(id="qc", rx=r"^\s*Convergence failure\.\s*$", links=(508,), scope="step",
         title="二次收敛 SCF 也失败（L508）",
         cause="SCF=QC/XQC/YQC 的二次收敛算法也没有收敛。",
         fixes=["去掉 QC 类选项，改用更好的初猜（Guess=Read）、SCF=VShift、SCF=Fermi 等方法。"],
         draft="qc"),
    dict(id="quadrature", rx=r"Inaccurate quadrature in CalDSu", links=(502, 1002), scope="step",
         title="DFT 数值积分不准（Inaccurate quadrature）",
         cause="数值积分得到的电子数与实际偏差过大；常见于基组/ECP 写错、结构异常，或弥散函数与拥挤结构同时出现。",
         fixes=["先检查基组、ECP 与结构；改用 Int=SuperFine（同一项目统一格点）；SCF=NoVarAcc；",
                "考虑去掉不必要的弥散函数。"],
         draft="quadrature"),
    dict(id="frozencore", rx=r"Excessive mixing of frozen core and valence orbitals", links=(801,),
         scope="step",
         title="冻芯轨道与价轨道混合过大",
         cause="后 HF、双杂化或 TD 计算中的冻芯设置不合理（常见于含重元素或 ECP 的体系）。",
         fixes=["改用全电子相关（如 MP2(Full)），或调整冻芯设置。"], draft=None),
    dict(id="optsteps", rx=r"Number of steps exceeded", links=(9999,), scope="step",
         title="优化步数用完（L9999 Number of steps exceeded）",
         cause="几何优化（或柔性扫描的某个点）在 MaxCycles 步内没有收敛。",
         fixes=["先诊断：python opt_monitor.py <log> 画出能量与 4 个判据（振荡 / 缓慢下降 / 平坦势能面）；",
                "缓慢下降：Opt=(Restart,MaxCycles=200)（需要 chk）或从最后结构 Opt=(MaxCycles=200)；",
                "振荡：Opt=CalcFC、Opt=(MaxStep=10)、Int=SuperFine，从能量最低的一步重启；",
                "结构越走越离谱：检查初始结构、电荷与多重度。"],
         draft="optsteps"),
    dict(id="galloc", rx=r"galloc:\s+could not allocate memory", links=(), scope="step",
         title="内存分配失败（galloc）",
         cause="%mem 超过节点当前可用的内存。",
         fixes=["降低 %mem，取调度系统申请内存的 75-85%；检查 Default.Route 或环境变量 GAUSS_MDEF 是否另有设置。"],
         draft=None),
    dict(id="write", rx=r"Erroneous write", links=(), scope="step",
         title="写 scratch 文件失败（Erroneous write，通常是磁盘满）",
         cause="GAUSS_SCRDIR 所在磁盘写满、超出配额，或目录不可写。",
         fixes=["df -h $GAUSS_SCRDIR、quota -s 检查空间与配额；清理以前遗留的 Gau-*.rwf；",
                "把 GAUSS_SCRDIR 换到更大的本地盘，或用 %RWF=/disk1/,50GB,/disk2/,-1 拆分。"],
         draft="same"),
    dict(id="read", rx=r"Erroneous read", links=(), scope="step",
         title="读文件失败（Erroneous read）",
         cause="chk/rwf 文件不完整或损坏（例如上一个作业写文件时被杀），或磁盘问题。",
         fixes=["不要再读这个 chk；从 log 中的最后结构重新开始（--draft），必要时去掉 Guess=Read/Geom=Check。"],
         draft="same"),
    dict(id="ntrerr", rx=r"NtrErr Called from FileIO", links=(), scope="step",
         title="chk/rwf 读写问题（NtrErr Called from FileIO）",
         cause="Geom=Check/Guess=Read 指定的 chk 不存在或名字写错，chk 来自其他平台/版本，或 chk 已损坏。",
         fixes=["核对文件名：用 %oldchk 读、%chk 写；跨平台/跨版本先 formchk 再 unfchk；",
                "C.01 之前的 chk 用 c8616 升级（先备份）。真正的原因在这一行更上方。"], draft=None),
    dict(id="segv", rx=r"segmentation violation|Severe Error Message #2070", links=(), scope="step",
         title="程序异常退出（segmentation violation / #2070）",
         cause="这一信息本身不说明原因；常见根本原因有 %mem 超过物理内存、磁盘写满、输入中的特殊错误、个别修订版的缺陷。",
         fixes=["往上找最后一条具体的错误或警告；降低 %mem 和核数重试；换机器或换修订版验证能否重复。"],
         draft=None),
    dict(id="link1", rx=r"Link 1 not found|l1\.exe.*not found", links=(), scope="step",
         title="找不到 Gaussian 程序（Link 1 not found）",
         cause="没有加载 g16.profile、g16root/GAUSS_EXEDIR 错误；Windows 上路径含中文或空格、杀毒软件隔离。",
         fixes=["作业脚本中 source $g16root/g16/bsd/g16.profile；Windows 把文件放到纯英文、无空格的路径。"],
         draft=None),
]
for _k in KNOWN:
    _k["re"] = re.compile(_k["rx"], re.IGNORECASE)
KNOWN_BY_ID = {k["id"]: k for k in KNOWN}
RE_ANY_KNOWN = re.compile("|".join("(?:%s)" % k["rx"] for k in KNOWN), re.IGNORECASE)
NEAR_LINES = 80


# ---------------------------------------------------------------------------------------------
# parsing
# ---------------------------------------------------------------------------------------------
def parse_log(path, context=12):
    """Single pass over a Gaussian output file; returns a dict with everything the doctor needs."""
    info = {
        "file": path, "routes": [], "step_starts": [], "charge": None, "mult": None,
        "normal": 0, "error_line": None, "error_text": None, "link": None, "context": [],
        "matches": [], "conv_not_met": 0, "scf_done": 0, "scf_done_ok": 0, "maxcycles": None,
        "stationary": 0, "s2": None, "freqs": None, "instability": None, "last_line": "",
        "nlines": 0,
    }
    tail = collections.deque(maxlen=max(context, 1))
    prev = ""
    capturing = False
    route_buf = []
    pending_notmet = False
    freq_block = None
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for n, raw in enumerate(fh, 1):
            line = raw.rstrip("\r\n")
            info["nlines"] = n
            if line.strip():
                info["last_line"] = line.strip()
            if capturing:
                if RE_DASH.match(line):
                    capturing = False
                    info["routes"].append("".join(route_buf).strip())
                else:
                    route_buf.append(line[1:] if line.startswith(" ") else line)
                prev = line
                tail.append(line)
                continue
            if line.startswith(" #") and RE_DASH.match(prev):
                capturing = True
                route_buf = [line[1:]]
                info["step_starts"].append(n)
                prev = line
                tail.append(line)
                continue
            prev = line

            if info["error_line"] is None and "Error termination" in line:
                info["context"] = list(tail)
            if "Error termination" in line:
                m = RE_LNK1E.search(line)
                if m:
                    info["link"] = int(m.group(2))
                    info["error_text"] = line.strip()
                    info["error_line"] = info["error_line"] or n
                else:
                    m = RE_LINKREQ.search(line)
                    if m and info["link"] is None:
                        info["link"] = int(m.group(1))
                    if info["error_line"] is None:
                        info["error_line"] = n
                        info["error_text"] = line.strip()
            tail.append(line)

            if "Normal termination of Gaussian" in line:
                info["normal"] += 1
                continue
            if info["charge"] is None:
                m = RE_CHARGE.search(line)
                if m:
                    info["charge"], info["mult"] = int(m.group(1)), int(m.group(2))
                    continue
            if "Convergence criterion not met" in line:
                info["conv_not_met"] += 1
                pending_notmet = True
                continue
            if "SCF Done:" in line:
                info["scf_done"] += 1
                if not pending_notmet:
                    info["scf_done_ok"] += 1
                pending_notmet = False
                continue
            if "Stationary point found" in line:
                info["stationary"] += 1
                continue
            m = RE_STEPMAX.search(line)
            if m:
                info["maxcycles"] = int(m.group(1))
                continue
            m = RE_S2.search(line)
            if m:
                info["s2"] = (float(m.group(1)), float(m.group(2)))
                continue
            if "Harmonic frequencies (cm**-1)" in line:
                freq_block = []
                info["freqs"] = freq_block
                continue
            m = RE_FREQ.match(line)
            if m and freq_block is not None:
                for t in m.group(1).split():
                    try:
                        freq_block.append(float(t))
                    except ValueError:
                        pass
                continue
            m = RE_INSTAB.search(line)
            if m:
                info["instability"] = line.strip()
                continue
            if "wavefunction is stable under the perturbations considered" in line:
                info["instability"] = None
                continue
            if not RE_ANY_KNOWN.search(line):
                continue
            for k in KNOWN:
                if k["re"].search(line):
                    info["matches"].append((n, k["id"], line.strip()))
    return info


def failing_route(info):
    """Route of the job step that contains the error (auto-generated Freq routes skipped)."""
    routes = info["routes"]
    if not routes:
        return None
    user = [r for r in routes if not r.upper().startswith("#N GEOM=ALLCHECK")]
    if info["error_line"] is not None:
        idx = [i for i, s in enumerate(info["step_starts"]) if s <= info["error_line"]]
        if idx:
            r = routes[idx[-1]]
            if not r.upper().startswith("#N GEOM=ALLCHECK"):
                return r
    return user[-1] if user else routes[-1]


def select_matches(info):
    """Known messages inside the failing job step (or the last NEAR_LINES lines for 'near' ones)."""
    end = info["error_line"] or info["nlines"]
    starts = [s for s in info["step_starts"] if s <= end]
    start = starts[-1] if starts else 1
    if info["error_line"] is None:            # killed/truncated: only look at the tail of the file
        start = max(start, end - 400)
    found = []
    seen = set()
    for n, kid, text in info["matches"]:
        k = KNOWN_BY_ID[kid]
        lo = end - NEAR_LINES if k["scope"] == "near" else start
        if lo <= n <= end + 5 and kid not in seen:
            if kid == "qc" and info["link"] not in (508, None):
                continue
            seen.add(kid)
            found.append((n, kid, text))
    return found


# ---------------------------------------------------------------------------------------------
# diagnosis
# ---------------------------------------------------------------------------------------------
def ideal_s2(mult):
    s = (mult - 1) / 2.0
    return s * (s + 1.0)


def diagnose(info):
    """Return a dict: status, link, problems (list of KNOWN ids + texts), notes (Chinese strings)."""
    route = failing_route(info) or ""
    res = {"file": info["file"], "route": route, "link": info["link"], "error": info["error_text"],
           "problems": [], "notes": [], "status": None}
    last_normal = "Normal termination of Gaussian" in info["last_line"]
    if info["error_line"] is not None:
        res["status"] = "error"
    elif last_normal and info["normal"] > 0:
        res["status"] = "normal"
    else:
        res["status"] = "incomplete"

    if res["status"] != "normal":
        for n, kid, text in select_matches(info):
            res["problems"].append({"id": kid, "line": n, "text": text})
        if res["status"] == "error" and not res["problems"]:
            res["notes"].append(LINK_INFO.get(info["link"], "L%s：没有内置说明的 Link。" % info["link"])
                                + " 没有匹配到已知报错原文，请阅读 Error termination 上方的内容（--context 加大行数），"
                                  "并查 https://liyuanhe211.github.io/Solution_for_Every_Gaussian_Error_Message/ 。")
        if res["status"] == "incomplete":
            age = time.time() - os.path.getmtime(info["file"])
            if age < 300:
                res["notes"].append("文件在 %.0f 秒前还在更新：作业可能仍在运行。" % age)
            else:
                res["notes"].append("没有 Normal/Error termination：作业被调度系统杀掉（超时、超内存、scancel）"
                                    "或日志被截断。用 sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS "
                                    "（PBS: qstat -f / tracejob）查看原因；优化作业可从最后结构续算（--draft）。")
        if info["conv_not_met"] and not any(p["id"] in ("scf", "qc") for p in res["problems"]):
            res["notes"].append("日志中有 %d 处 'Convergence criterion not met'：对应的 SCF Done 能量不可用。"
                                % info["conv_not_met"])
        return res

    # normal termination: sanity checks
    notes = res["notes"]
    notes.append("作业正常结束（%d 次 Normal termination）。" % info["normal"])
    if info["conv_not_met"]:
        notes.append("警告：有 %d 处 'Convergence criterion not met'（SCF 未收敛但作业继续），"
                     "相应能量不可用，检查是否用了 SCF=(MaxCycle=N) 之外的放宽选项。" % info["conv_not_met"])
    if RE_OPT_KW.search(route) and info["stationary"] == 0:
        notes.append("警告：route 含 Opt，但没有 'Stationary point found'，优化可能没有真正收敛。")
    if info["freqs"] is not None:
        nimag = sum(1 for f in info["freqs"] if f < 0)
        is_ts = bool(RE_TS_KW.search(" ".join(info["routes"])))
        want = 1 if is_ts else 0
        msg = "最后一次频率计算：%d 个虚频（%s）。" % (
            nimag, ", ".join("%.1f" % f for f in info["freqs"] if f < 0) or "无")
        if nimag != want:
            msg += (" 期望 %d 个（%s）：" % (want, "TS" if is_ts else "极小点")) + (
                "TS 若没有虚频说明优化到了极小点，有 2 个以上需沿多余虚频方向调整结构。" if is_ts else
                "沿虚频方向位移结构后重新优化；频率必须与优化用完全相同的方法/基组/溶剂/格点；"
                "很小的虚频（约 -10 cm^-1 量级）可能是数值噪声，可用 Opt=Tight、Int=SuperFine 复核。")
        notes.append(msg)
    if info["s2"] and info["mult"]:
        before, after = info["s2"]
        ideal = ideal_s2(info["mult"])
        if (ideal == 0 and before > 0.1) or (ideal > 0 and abs(before - ideal) > 0.1 * ideal):
            notes.append("<S**2> 湮灭前 %.4f（理想值 %.4f，多重度 %d），湮灭后 %.4f：自旋污染明显。"
                         "若是有意构造的 broken-symmetry 态可忽略；否则检查是否收敛到了错误的态（Stable=Opt）"
                         "或考虑 ROKS/其他泛函。" % (before, ideal, info["mult"], after))
    if info["instability"]:
        notes.append("Stable 检查发现不稳定：%s → 用 Stable=Opt 得到稳定的波函数再继续。" % info["instability"])
    return res


# ---------------------------------------------------------------------------------------------
# route editing for --draft
# ---------------------------------------------------------------------------------------------
def split_route(route):
    tokens, buf, depth = [], "", 0
    for ch in route:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth = max(depth - 1, 0)
        if ch in " \t" and depth == 0:
            if buf:
                tokens.append(buf)
            buf = ""
        else:
            buf += ch
    if buf:
        tokens.append(buf)
    return tokens


def parse_token(tok):
    """'Opt=(TS,CalcFC)' -> ('opt', 'Opt', ['TS','CalcFC']); None for method/basis tokens."""
    if "/" in tok or tok.startswith("#"):
        return None
    m = re.match(r"^([A-Za-z][A-Za-z0-9]*)\s*(?:=\s*\((.*)\)|=\s*(.*)|\((.*)\))?$", tok)
    if not m:
        return None
    inner = m.group(2) if m.group(2) is not None else (m.group(3) if m.group(3) is not None else m.group(4))
    opts = [o for o in re.split(r",(?![^(]*\))", inner or "") if o.strip()]
    return m.group(1).lower(), m.group(1), [o.strip() for o in opts]


def build_token(name, opts):
    if not opts:
        return name
    if len(opts) == 1 and "(" not in opts[0]:
        return "%s=%s" % (name, opts[0])
    return "%s=(%s)" % (name, ",".join(opts))


def edit_route(route, strategy, info, warnings):
    """Return (new_route, guess_read) for a draft strategy."""
    tokens = split_route(route)
    out = []
    seen = {}
    for tok in tokens:
        p = parse_token(tok)
        if p:
            key, name, opts = p
            low = [o.lower() for o in opts]
            if key == "geom" and any(w in " ".join(low) for w in ("check", "connect")):
                continue
            if key == "guess" and any(o.startswith("read") for o in low):
                continue
            if key == "opt":
                opts = [o for o in opts if o.lower() != "restart"]
            seen[key if key != "integral" else "int"] = len(out)
            out.append([key, name, opts])
        else:
            out.append([None, tok, None])

    def get(key, default_name):
        if key in seen:
            return out[seen[key]]
        out.append([key, default_name, []])
        seen[key] = len(out) - 1
        return out[-1]

    has_opt = "opt" in seen
    guess = has_opt and info["scf_done_ok"] > 0
    if strategy == "scf":
        scf = get("scf", "SCF")
        if any(o.lower() in ("qc", "xqc", "yqc") for o in scf[2]):
            warnings.append("route 已含 QC/XQC：草稿未改 SCF 设置，请手工换初猜（小基组收敛后 Guess=Read）"
                            "或用 SCF=VShift=500、SCF=Fermi")
        else:
            scf[2].append("XQC")
    elif strategy == "qc":
        scf = get("scf", "SCF")
        scf[2] = [o for o in scf[2] if o.lower() not in ("qc", "xqc", "yqc")] + ["VShift=500"]
    elif strategy == "quadrature":
        integ = get("int", "Int")
        integ[2] = [o for o in integ[2] if not re.match(r"(?i)(grid=)?(ultra|super)?fine$|grid=", o)]
        integ[2].append("SuperFine")
    elif strategy == "optsteps":
        if not has_opt:
            warnings.append("route 中没有 Opt，未修改 MaxCycles")
        else:
            opt = out[seen["opt"]]
            old = info["maxcycles"] or 100
            opt[2] = [o for o in opt[2] if not o.lower().startswith("maxcyc")]
            opt[2].append("MaxCycles=%d" % max(200, 2 * old))
            warnings.append("先用 opt_monitor.py 判断是缓慢下降还是振荡；振荡时请改用能量最低的结构并加 CalcFC")
    elif strategy == "negeig":
        if not has_opt:
            warnings.append("route 中没有 Opt")
        else:
            opt = out[seen["opt"]]
            low = [o.lower() for o in opt[2]]
            if not any(o in ("calcfc", "readfc", "calcall", "rcfc") or o.startswith("recalcfc") for o in low):
                opt[2].append("CalcFC")
            if "noeigentest" not in low:
                opt[2].append("NoEigenTest")
    elif strategy == "same":
        if has_opt:
            warnings.append("route 未改动：从最后结构重新开始；若再次出现线性角问题改用 Opt=Cartesian")
    if guess:
        out.append([None, "Guess=Read", None])
    pieces = [build_token(t[1], t[2]) if t[0] else t[1] for t in out]
    if not pieces or not pieces[0].startswith("#"):
        pieces.insert(0, "#p")
    return " ".join(pieces), guess


def write_draft(info, res, out_path, force):
    strategies = [KNOWN_BY_ID[p["id"]]["draft"] for p in res["problems"]]
    strategies = [s for s in strategies if s]
    if res["status"] == "normal":
        print("  [draft] 作业正常结束，不生成续算输入。")
        return False
    if res["problems"] and not strategies:
        print("  [draft] 这类错误需要先修改输入（或资源设置），不适合自动起草续算输入。")
        return False
    strategy = strategies[0] if strategies else "same"
    if os.path.exists(out_path) and not force:
        print("  [draft] %s 已存在；用 --force 覆盖或换文件名。" % out_path)
        return False
    route = failing_route(info)
    if not route:
        print("  [draft] 日志中没有 route，无法起草。")
        return False
    warnings = []
    new_route, guess = edit_route(route, strategy, info, warnings)
    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, here)
    try:
        import log2xyz  # sibling script of DFT-Tutorials
    except ImportError:
        print("  [draft] 需要同目录下的 log2xyz.py（DFT-Tutorials Gaussian/scripts）。")
        return False
    finally:
        sys.path.pop(0)
    xyz = os.path.splitext(out_path)[0] + "_geom.xyz"
    argv = [info["file"], "--last", "-o", xyz, "--gjf", out_path, "--route=" + new_route]
    if guess:
        argv.append("--guess-read")
    rc = log2xyz.main(argv)
    if rc != 0:
        print("  [draft] log2xyz.py 未能写出续算输入（见上方错误）。")
        return False
    print("  [draft] 已起草 %s（route: %s）。这只是草稿：提交前请逐行检查 Link 0、route、电荷多重度、"
          "坐标及 Gen/ModRedundant 段。" % (out_path, new_route))
    for w in warnings:
        print("  [draft] 注意：" + w)
    return True


# ---------------------------------------------------------------------------------------------
# output
# ---------------------------------------------------------------------------------------------
def print_report(info, res, quiet=False):
    tag = {"normal": "OK   ", "error": "ERROR", "incomplete": "INCOMPLETE"}[res["status"]]
    first = res["problems"][0]["id"] if res["problems"] else ""
    one = "%-10s %s" % (tag, info["file"])
    if res["link"] is not None:
        one += "  L%d" % res["link"]
    if first:
        one += "  " + KNOWN_BY_ID[first]["title"]
    if quiet:
        print(one)
        return
    print("=" * 78)
    print(one)
    if res["route"]:
        print("route: %s" % res["route"])
    if res["error"]:
        print("出错行: %s" % res["error"])
    if res["status"] == "error" and info["context"]:
        print("-" * 30 + " Error termination 之前的 %d 行 " % len(info["context"]) + "-" * 20)
        for ln in info["context"]:
            print("  | " + ln)
        print("-" * 78)
    if res["link"] in LINK_INFO and res["problems"]:
        print(LINK_INFO[res["link"]])
    for p in res["problems"]:
        k = KNOWN_BY_ID[p["id"]]
        print("\n[%s] %s" % ("/".join("L%d" % l for l in k["links"]) or "任意 Link", k["title"]))
        print("  原文（第 %d 行）: %s" % (p["line"], p["text"]))
        print("  原因: " + k["cause"])
        print("  处理:")
        for f in k["fixes"]:
            print("    - " + f)
    for note in res["notes"]:
        print("* " + note)
    if res["status"] != "normal":
        print("更多报错见 Gaussian 第18章（Gaussian/tutorials/18-errors-capstone.md）。")


def list_errors():
    print("%-12s %-10s %s" % ("id", "Link", "message (regex) -> 说明"))
    for k in KNOWN:
        print("%-12s %-10s %s -> %s" % (k["id"], "/".join("L%d" % l for l in k["links"]) or "any",
                                         k["rx"], k["title"]))


def collect(paths, recursive):
    files = []
    for p in paths:
        if os.path.isdir(p):
            pat = "**" if recursive else ""
            for ext in ("*.log", "*.out"):
                files += sorted(glob.glob(os.path.join(p, pat, ext), recursive=recursive))
        else:
            files.append(p)
    return files


def build_parser():
    p = argparse.ArgumentParser(
        description="Diagnose Gaussian 16 jobs: failing Link, known error messages with Chinese causes "
                    "and fixes, post-run sanity checks, optional restart draft.  "
                    "识别出错的 Link 与常见报错，给出中文原因和建议，可起草续算输入。",
        epilog="Examples:\n"
               "  python g16_error_doctor.py job.log\n"
               "  python g16_error_doctor.py *.log -q\n"
               "  python g16_error_doctor.py runs/ -r --json report.json\n"
               "  python g16_error_doctor.py job.log --draft            # writes job_restart.gjf\n"
               "  python g16_error_doctor.py --list-errors",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("logs", nargs="*", help="Gaussian .log/.out files or directories / 输出文件或目录")
    p.add_argument("-r", "--recursive", action="store_true", help="search directories recursively")
    p.add_argument("-c", "--context", type=int, default=12, metavar="N",
                   help="lines shown above 'Error termination' (default 12)")
    p.add_argument("-q", "--quiet", action="store_true", help="one status line per file")
    p.add_argument("--json", metavar="FILE", help="also write the diagnosis as JSON ('-' = stdout)")
    p.add_argument("--draft", nargs="?", const="", metavar="GJF",
                   help="draft a restart input (single log; default <log stem>_restart.gjf) - REVIEW IT")
    p.add_argument("--force", action="store_true", help="allow --draft to overwrite an existing file")
    p.add_argument("--list-errors", action="store_true", help="print the recognised messages and exit")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.list_errors:
        list_errors()
        return 0
    if not args.logs:
        build_parser().print_usage(sys.stderr)
        print("ERROR: give at least one log file or directory", file=sys.stderr)
        return 2
    files = collect(args.logs, args.recursive)
    if not files:
        print("ERROR: no .log/.out files found", file=sys.stderr)
        return 2
    if args.draft is not None and len(files) != 1:
        print("ERROR: --draft works with exactly one log file", file=sys.stderr)
        return 2
    rc = 0
    report = []
    for f in files:
        if not os.path.isfile(f):
            print("ERROR: %s not found" % f, file=sys.stderr)
            rc = 1
            continue
        try:
            info = parse_log(f, context=args.context)
        except OSError as exc:
            print("ERROR: cannot read %s: %s" % (f, exc), file=sys.stderr)
            rc = 1
            continue
        if info["nlines"] == 0 or (not info["routes"] and info["normal"] == 0 and info["error_line"] is None
                                   and not info["matches"]):
            print("ERROR: %s does not look like a Gaussian output (no route section)" % f, file=sys.stderr)
            rc = 1
            continue
        res = diagnose(info)
        print_report(info, res, quiet=args.quiet)
        report.append(res)
        if res["status"] != "normal":
            rc = 1
        if args.draft is not None:
            out = args.draft or os.path.splitext(f)[0] + "_restart.gjf"
            if os.path.abspath(out) == os.path.abspath(f):
                print("ERROR: --draft would overwrite the log", file=sys.stderr)
                return 2
            write_draft(info, res, out, args.force)
    if args.json:
        text = json.dumps(report, ensure_ascii=False, indent=2)
        if args.json == "-":
            print(text)
        else:
            with open(args.json, "w", encoding="utf-8") as fh:
                fh.write(text + "\n")
            print("JSON written to %s" % args.json, file=sys.stderr)
    return rc


if __name__ == "__main__":
    sys.exit(main())
