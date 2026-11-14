#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
omx_traj.py -- convert and plot OpenMX optimisation / MD trajectories (System.Name.md).
               把 .md 轨迹转为 extxyz / ASE Trajectory，绘制能量、最大力、温度随步数的变化，
               并可用 dpdata 导出 DeePMD 等机器学习势训练数据。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

Input format (OpenMX 4.0, written by iterout.c in APPEND mode, one frame per step)
-------------------------------------------------------------------------------
    N
    time= t (fs) Energy= U (Hartree) [Temperature= T (Given Temp.= T0)] Cell_Vectors= 9 x Ang  scf_conv=1
    Sym  x y z  fx fy fz  vx vy vz  net_charge  mag_moment  theta phi
    ...
* x y z in Angstrom; f = -dE/dR in Hartree/Bohr; v in m/s (written only for
  NVE / NVT_VS / NVT_VS2 / NVT_VS4 / NVT_NH / NVT_Langevin; zero for
  optimisation and for the NPT integrators); theta/phi are non-collinear spin
  angles in degrees.
* NVT_VS4 prints "Temperature[1]= ... Temperature[2]= ..." per atom group.
* Older files with fewer columns (only coordinates) are read as far as possible.
Because the file is appended, re-running a job with the same System.Name
puts a second trajectory after the first one; the script warns when the time
goes backwards and --last-run keeps only the frames after the last restart.

Optional:  <prefix>.ene (MD only; col 2 time, 15 Utot, 16 Utot+Ukc,
19 temperature, 22 Nose-Hoover Hamiltonian, Hartree) is plotted as the
conserved quantity when it exists.

Usage examples
--------------
    python omx_traj.py h2o_rf.md                       # summary + h2o_rf_traj.png
    python omx_traj.py h2o_rf --fmax-line 0.0003 --force-unit au
    python omx_traj.py H2O_NVE --extxyz H2O_NVE.extxyz --skip 100 --stride 2
    python omx_traj.py Si64_NVT_NH --traj si64.traj        # needs ASE
    python omx_traj.py Si64_NVT_NH --dpdata deepmd_si64 --skip 200 --stride 5
    python omx_traj.py Si64_NVT_NH --csv si64_steps.csv --no-plot

Units written: energy eV, forces eV/Angstrom (Hartree -> eV 27.211386245988,
Bohr -> Angstrom 0.529177210903).  MD.Opt.criterion of OpenMX is a maximum
force in Hartree/Bohr; use --force-unit au to compare directly.

Dependencies: Python >= 3.9, numpy, matplotlib (plot); optional ase (--traj),
dpdata (--dpdata, reads <prefix>.dat + <prefix>.md with fmt="openmx/md").
"""

import argparse
import os
import re
import sys

HA2EV = 27.211386245988
BOHR2ANG = 0.529177210903
HAB2EVA = HA2EV / BOHR2ANG


def die(msg):
    sys.exit("omx_traj.py: error: " + msg)


def num_after(tokens, key):
    if key in tokens:
        i = tokens.index(key)
        if i + 1 < len(tokens):
            try:
                return float(tokens[i + 1])
            except ValueError:
                return None
    return None


def read_md(fname):
    """Return list of frame dicts."""
    import numpy as np

    frames = []
    try:
        f = open(fname)
    except OSError as exc:
        die("cannot read %s (%s)" % (fname, exc))
    with f:
        lineno = 0
        while True:
            line = f.readline()
            lineno += 1
            if not line:
                break
            if not line.strip():
                continue
            try:
                nat = int(line.split()[0])
            except ValueError:
                die("%s line %d: expected the number of atoms, got: %s" % (fname, lineno, line.strip()))
            head = f.readline()
            lineno += 1
            tok = head.split()
            fr = {"time": num_after(tok, "time="), "energy_ha": num_after(tok, "Energy="),
                  "temp": num_after(tok, "Temperature="), "given_temp": None,
                  "group_temps": [float(v) for v in re.findall(r"Temperature\[\d+\]=\s*([-+0-9.eE]+)", head)],
                  "cell": None, "scf_conv": None}
            m = re.search(r"Given Temp\.=\s*([-+0-9.eE]+)", head)
            if m:
                fr["given_temp"] = float(m.group(1))
            if "Cell_Vectors=" in head:
                try:
                    fr["cell"] = np.array(head.split("Cell_Vectors=")[1].split()[:9], float).reshape(3, 3)
                except ValueError:
                    die("%s line %d: cannot read Cell_Vectors" % (fname, lineno))
            m = re.search(r"scf_conv=\s*(\d+)", head)
            if m:
                fr["scf_conv"] = int(m.group(1))
            sym, pos, frc, vel, chg, mag = [], [], [], [], [], []
            for _ in range(nat):
                q = f.readline().split()
                lineno += 1
                if len(q) < 4:
                    die("%s line %d: incomplete atom line (truncated frame?)" % (fname, lineno))
                try:
                    sym.append(q[0])
                    pos.append([float(v) for v in q[1:4]])
                    if len(q) >= 7:
                        frc.append([float(v) for v in q[4:7]])
                    if len(q) >= 10:
                        vel.append([float(v) for v in q[7:10]])
                    if len(q) >= 12:
                        chg.append(float(q[10]))
                        mag.append(float(q[11]))
                except ValueError:
                    die("%s line %d: non-numeric value" % (fname, lineno))
            fr["symbols"] = sym
            fr["pos"] = np.array(pos)
            fr["forces"] = np.array(frc) * HAB2EVA if len(frc) == nat else None
            fr["vel"] = np.array(vel) if len(vel) == nat else None
            fr["charge"] = np.array(chg) if len(chg) == nat else None
            fr["magmom"] = np.array(mag) if len(mag) == nat else None
            frames.append(fr)
    if not frames:
        die("no frame found in %s" % fname)
    return frames


def parse_index_list(text, nat):
    out = set()
    for tok in re.split(r"[,\s]+", text.strip()):
        if not tok:
            continue
        m = re.fullmatch(r"(\d+)-(\d+)", tok)
        try:
            out.update(range(int(m.group(1)), int(m.group(2)) + 1) if m else [int(tok)])
        except ValueError:
            die("bad atom list item: %s" % tok)
    bad = [i for i in out if i < 1 or i > nat]
    if bad:
        die("atom index out of range 1..%d: %s" % (nat, sorted(bad)))
    return sorted(out)


def write_extxyz(fname, frames, pbc, write_vel):
    with open(fname, "w") as f:
        for fr in frames:
            nat = len(fr["symbols"])
            props = "species:S:1:pos:R:3"
            if fr["forces"] is not None:
                props += ":forces:R:3"
            vel = write_vel and fr["vel"] is not None
            if vel:
                props += ":vel_m_per_s:R:3"
            info = []
            if fr["cell"] is not None:
                info.append('Lattice="%s"' % " ".join("%.8f" % v for v in fr["cell"].ravel()))
            info.append('Properties=%s' % props)
            if fr["energy_ha"] is not None:
                info.append("energy=%.8f" % (fr["energy_ha"] * HA2EV))
            if fr["time"] is not None:
                info.append("time_fs=%.4f" % fr["time"])
            if fr["temp"] is not None:
                info.append("temperature_K=%.3f" % fr["temp"])
            info.append('pbc="%s"' % ("T T T" if (pbc and fr["cell"] is not None) else "F F F"))
            f.write("%d\n%s\n" % (nat, " ".join(info)))
            for i in range(nat):
                row = "%-3s %15.8f %15.8f %15.8f" % ((fr["symbols"][i],) + tuple(fr["pos"][i]))
                if fr["forces"] is not None:
                    row += " %15.8f %15.8f %15.8f" % tuple(fr["forces"][i])
                if vel:
                    row += " %12.4f %12.4f %12.4f" % tuple(fr["vel"][i])
                f.write(row + "\n")


def write_ase_traj(fname, frames, pbc):
    try:
        from ase import Atoms
        from ase import units
        from ase.calculators.singlepoint import SinglePointCalculator
        from ase.io.trajectory import Trajectory
    except ImportError:
        die("--traj needs ASE: pip install ase")
    with Trajectory(fname, "w") as traj:
        for fr in frames:
            symbols = [re.match(r"[A-Z][a-z]?", s).group(0) if re.match(r"[A-Z][a-z]?", s) else "X"
                       for s in fr["symbols"]]
            atoms = Atoms(symbols=symbols, positions=fr["pos"],
                          cell=fr["cell"] if fr["cell"] is not None else None,
                          pbc=bool(pbc and fr["cell"] is not None))
            if fr["vel"] is not None and abs(fr["vel"]).max() > 0:
                # m/s -> Angstrom/fs (1e-5) -> ASE velocity unit
                atoms.set_velocities(fr["vel"] * 1e-5 / units.fs)
            kw = {}
            if fr["energy_ha"] is not None:
                kw["energy"] = fr["energy_ha"] * HA2EV
            if fr["forces"] is not None:
                kw["forces"] = fr["forces"]
            if kw:
                atoms.calc = SinglePointCalculator(atoms, **kw)
            traj.write(atoms)


def read_ene(fname):
    import numpy as np

    try:
        d = np.loadtxt(fname, ndmin=2)
    except (OSError, ValueError) as exc:
        print("warning: cannot read %s (%s)" % (fname, exc))
        return None
    if d.shape[1] < 22:
        print("warning: %s has %d columns (< 22); ignored" % (fname, d.shape[1]))
        return None
    return d


def main():
    ap = argparse.ArgumentParser(
        prog="omx_traj.py",
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("md", help="System.Name.md or System.Name (prefix)")
    ap.add_argument("--skip", type=int, default=0, help="drop the first N frames (equilibration)")
    ap.add_argument("--stride", type=int, default=1, help="keep every N-th frame for the exports")
    ap.add_argument("--last-run", action="store_true",
                    help="keep only frames after the last time reset (appended reruns)")
    ap.add_argument("--exclude", default="", help='atoms ignored for the max force, e.g. fixed ones "1-8"')
    ap.add_argument("--force-unit", choices=["ev", "au"], default="ev",
                    help="unit of the printed/plotted max force: eV/Ang (default) or Hartree/Bohr")
    ap.add_argument("--fmax-line", type=float, default=None,
                    help="draw a horizontal criterion line in the force panel (in --force-unit)")
    ap.add_argument("--extxyz", default="", help="write extended XYZ (energy eV, forces eV/Ang)")
    ap.add_argument("--velocities", action="store_true", help="add velocities (m/s) to the extxyz file")
    ap.add_argument("--traj", default="", help="write an ASE .traj file (needs ase)")
    ap.add_argument("--nopbc", action="store_true", help="mark structures as non-periodic (Cluster runs)")
    ap.add_argument("--dpdata", default="", help="export to this directory with dpdata (openmx/md reader)")
    ap.add_argument("--dpdata-fmt", default="deepmd/npy", help="dpdata output format (default deepmd/npy)")
    ap.add_argument("--csv", default="", help="per-step table: step, time, E, dE, Fmax, T")
    ap.add_argument("--ene", default="", help="<prefix>.ene to plot (default: auto-detect)")
    ap.add_argument("-o", "--out", default="", help="figure name (default <prefix>_traj.png)")
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--show", action="store_true", help="interactive window instead of Agg")
    args = ap.parse_args()

    import numpy as np

    mdfile = args.md if args.md.endswith(".md") else args.md + ".md"
    if not os.path.exists(mdfile):
        die("file not found: %s" % mdfile)
    prefix = mdfile[:-3]
    if args.stride < 1 or args.skip < 0:
        die("--stride must be >= 1 and --skip >= 0")

    frames = read_md(mdfile)
    n_total = len(frames)
    times = [fr["time"] for fr in frames]
    resets = [i for i in range(1, len(times))
              if times[i] is not None and times[i - 1] is not None and times[i] < times[i - 1] - 1e-9]
    if resets:
        print("warning: time goes backwards at frame(s) %s - the file contains appended runs"
              % ", ".join(str(i + 1) for i in resets[:5]))
        if args.last_run:
            frames = frames[resets[-1]:]
            print("--last-run: keeping %d frames of the last run" % len(frames))
    nat = len(frames[0]["symbols"])
    if any(len(fr["symbols"]) != nat for fr in frames):
        die("the number of atoms changes between frames")
    frames = frames[args.skip:]
    if not frames:
        die("no frame left after --skip")

    excl = set(i - 1 for i in parse_index_list(args.exclude, nat)) if args.exclude else set()
    keep = np.array([i not in excl for i in range(nat)])
    fconv = 1.0 if args.force_unit == "ev" else 1.0 / HAB2EVA
    funit = "eV/Ang" if args.force_unit == "ev" else "Hartree/Bohr"

    step = np.arange(1, len(frames) + 1) + args.skip
    t = np.array([fr["time"] if fr["time"] is not None else np.nan for fr in frames])
    e = np.array([fr["energy_ha"] * HA2EV if fr["energy_ha"] is not None else np.nan for fr in frames])
    temp = np.array([fr["temp"] if fr["temp"] is not None else np.nan for fr in frames])
    gtemps = [fr["group_temps"] for fr in frames]
    fmax = np.array([np.max(np.linalg.norm(fr["forces"][keep], axis=1)) * fconv
                     if fr["forces"] is not None and keep.any() else np.nan for fr in frames])
    unconv = [i for i, fr in enumerate(frames) if fr["scf_conv"] == 0]
    is_md = np.isfinite(temp).any() or any(gtemps)

    print("%s: %d frames (after skip), %d atoms, %s" % (mdfile, len(frames), nat,
                                                       "MD" if is_md else "optimisation / single points"))
    if frames[0]["cell"] is not None and frames[-1]["cell"] is not None:
        v0, v1 = abs(np.linalg.det(frames[0]["cell"])), abs(np.linalg.det(frames[-1]["cell"]))
        print("cell volume: first %.3f  last %.3f Ang^3" % (v0, v1))
    if np.isfinite(e).any():
        print("energy (eV): first %.6f  last %.6f  min %.6f (frame %d)  last-first = %.3f meV/atom"
              % (e[0], e[-1], np.nanmin(e), int(np.nanargmin(e)) + 1 + args.skip,
                 (e[-1] - e[0]) / nat * 1000.0))
    if np.isfinite(fmax).any():
        print("max |F| (%s): first %.5f  last %.5f%s" % (funit, fmax[0], fmax[-1],
                                                        " (excluding %d atoms)" % len(excl) if excl else ""))
    if is_md and np.isfinite(temp).any():
        print("temperature (K): mean %.2f  std %.2f  (target %s)"
              % (np.nanmean(temp), np.nanstd(temp),
                 "%.1f" % frames[-1]["given_temp"] if frames[-1]["given_temp"] is not None else "n/a"))
    if unconv:
        print("WARNING: %d frame(s) with scf_conv=0 (SCF not converged), e.g. frame %d"
              % (len(unconv), unconv[0] + 1 + args.skip))

    ene = None
    enefile = args.ene or (prefix + ".ene" if os.path.exists(prefix + ".ene") else "")
    if enefile and is_md:
        ene = read_ene(enefile)
        if ene is not None:
            econs = ene[:, 15] * HA2EV
            drift = (econs[-1] - econs[0]) / nat * 1000.0
            print("%s: Utot+Ukc drift (last-first) = %.3f meV/atom over %d steps; span %.3f meV/atom"
                  % (enefile, drift, len(econs), (econs.max() - econs.min()) / nat * 1000.0))

    if args.csv:
        with open(args.csv, "w") as f:
            f.write("step,time_fs,energy_eV,dE_meV_per_atom,fmax_%s,temperature_K\n" % args.force_unit)
            for i in range(len(frames)):
                f.write("%d,%.4f,%.8f,%.4f,%.6f,%.3f\n"
                        % (step[i], t[i], e[i], (e[i] - e[0]) / nat * 1000.0, fmax[i], temp[i]))
        print("table written to %s" % args.csv)

    sub = frames[::args.stride]
    if args.extxyz:
        write_extxyz(args.extxyz, sub, not args.nopbc, args.velocities)
        print("%d frames written to %s" % (len(sub), args.extxyz))
    if args.traj:
        write_ase_traj(args.traj, sub, not args.nopbc)
        print("%d frames written to %s" % (len(sub), args.traj))
    if args.dpdata:
        try:
            import dpdata
        except ImportError:
            die("--dpdata needs dpdata: pip install dpdata")
        if not os.path.exists(prefix + ".dat"):
            die("dpdata's openmx/md reader needs %s.dat next to %s" % (prefix, mdfile))
        ls = dpdata.LabeledSystem(prefix, fmt="openmx/md")
        if len(ls) != n_total:
            print("warning: dpdata read %d frames, this script %d (appended runs?)" % (len(ls), n_total))
        if resets and args.last_run:
            print("note: --last-run is not applied to the dpdata export (dpdata reads the whole file)")
        part = ls[args.skip::args.stride]
        part.to(args.dpdata_fmt, args.dpdata)
        print("dpdata: %d frames written to %s (%s)" % (len(part), args.dpdata, args.dpdata_fmt))

    if args.no_plot:
        return
    import matplotlib
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    panels = ["E"]
    if np.isfinite(fmax).any():
        panels.append("F")
    if is_md:
        panels.append("T")
    fig, axes = plt.subplots(len(panels), 1, figsize=(6.0, 2.2 * len(panels) + 0.4), sharex=True)
    axes = np.atleast_1d(axes)
    xs = t if (is_md and np.isfinite(t).all()) else step
    xlabel = "time (fs)" if xs is t else "step"
    for ax, p in zip(axes, panels):
        if p == "E":
            ax.plot(xs, (e - e[0]) / nat * 1000.0, color="#1f4e79", lw=1.0, label=r"$U_\mathrm{tot}$")
            if ene is not None:
                xe = ene[:, 1] if xs is t else np.arange(1, len(ene) + 1)
                ax.plot(xe, (ene[:, 15] * HA2EV - e[0]) / nat * 1000.0, color="#c0504d", lw=1.0,
                        label=r"$U_\mathrm{tot}+U_\mathrm{kin}$ (.ene)")
                ax.legend(frameon=False, fontsize=8)
            ax.set_ylabel(r"$\Delta E$ (meV/atom)")
        elif p == "F":
            ax.plot(xs, fmax, color="#2e7d32", lw=1.0, marker="o" if len(xs) < 60 else None, ms=3)
            if not is_md and np.nanmin(fmax) > 0:
                ax.set_yscale("log")
            if args.fmax_line is not None:
                ax.axhline(args.fmax_line, color="0.4", ls="--", lw=0.8)
            ax.set_ylabel(r"max $|F|$ (%s)" % funit)
        elif p == "T":
            if np.isfinite(temp).any():
                ax.plot(xs, temp, color="#c0504d", lw=0.9)
            for g in range(max((len(v) for v in gtemps), default=0)):
                ax.plot(xs, [v[g] if len(v) > g else np.nan for v in gtemps], lw=0.9, label="group %d" % (g + 1))
            if any(gtemps):
                ax.legend(frameon=False, fontsize=8)
            gt = frames[-1]["given_temp"]
            if gt is not None:
                ax.axhline(gt, color="0.4", ls="--", lw=0.8)
            ax.set_ylabel("T (K)")
    axes[-1].set_xlabel(xlabel)
    fig.tight_layout()
    if args.show:
        plt.show()
    else:
        out = args.out or prefix + "_traj.png"
        fig.savefig(out, dpi=300)
        print("figure written to %s" % out)


if __name__ == "__main__":
    main()
