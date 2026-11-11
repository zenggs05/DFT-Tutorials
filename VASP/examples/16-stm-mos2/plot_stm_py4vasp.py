#!/usr/bin/env python3
"""Simulated STM images with py4vasp (VASP >= 6.5.0, HDF5 build, LPARDH5 = .TRUE.).

Run in the directory of the LPARD calculation that wrote vaspout.h5.
Saving PNG files needs the 'kaleido' package; otherwise an HTML file is written.

Examples
--------
    python plot_stm_py4vasp.py --mode constant_height --tip-height 3.0 --supercell 6
    python plot_stm_py4vasp.py --mode constant_current --current 1.0 --supercell 6
"""
import argparse

import py4vasp


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--path", default=".", help="directory containing vaspout.h5")
    p.add_argument("--mode", default="constant_height", choices=["constant_height", "constant_current"])
    p.add_argument("--spin", default="total", choices=["total", "up", "down"])
    p.add_argument("--tip-height", type=float, default=3.0, help="Angstrom, constant-height mode")
    p.add_argument("--current", type=float, default=1.0, help="nA, constant-current mode")
    p.add_argument("--supercell", type=int, default=6)
    p.add_argument("-o", "--output", default=None)
    args = p.parse_args()

    calc = py4vasp.Calculation.from_path(args.path)
    print(calc.partial_density.selections())
    selection = "%s(%s)" % (args.mode, args.spin)
    if args.mode == "constant_height":
        graph = calc.partial_density.to_stm(selection, tip_height=args.tip_height, supercell=args.supercell)
    else:
        graph = calc.partial_density.to_stm(selection, current=args.current, supercell=args.supercell)

    fig = graph.to_plotly()
    out = args.output or ("stm_%s.png" % args.mode)
    try:
        fig.write_image(out)
    except Exception as exc:          # kaleido missing or no display backend
        out = out.rsplit(".", 1)[0] + ".html"
        fig.write_html(out)
        print("PNG export failed (%s); wrote interactive HTML instead" % exc)
    print("written", out)


if __name__ == "__main__":
    main()
