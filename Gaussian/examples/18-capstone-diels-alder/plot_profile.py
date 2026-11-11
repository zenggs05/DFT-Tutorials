#!/usr/bin/env python3
"""Draw a simple free-energy profile (horizontal bars joined by dashed lines).

Values are typed on the command line as LABEL=VALUE pairs (kcal/mol, relative to the first
point), e.g. the numbers that goodvibes_run.sh printed:

    python plot_profile.py R=0 TS=<dG_TS> P=<dG_P> -o da_profile.png
    python plot_profile.py R=0 TS1=18.0 INT=3.0 TS2=12.5 P=-20.0 --ylabel "dG (kcal/mol)"

Labels starting with "TS" are drawn in red and annotated above the bar; the others below.
Requires Python >= 3.9 and matplotlib.
"""
import argparse
import sys


def parse_points(items):
    points = []
    for item in items:
        if "=" not in item:
            sys.exit(f"bad point '{item}': use LABEL=VALUE")
        label, value = item.split("=", 1)
        try:
            points.append((label, float(value)))
        except ValueError:
            sys.exit(f"bad value in '{item}': replace placeholders such as <dG_TS> by numbers")
    if len(points) < 2:
        sys.exit("give at least two points")
    return points


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("points", nargs="+", help="LABEL=VALUE pairs in reaction order")
    ap.add_argument("-o", "--output", default="profile.png", help="output image (png, pdf, svg)")
    ap.add_argument("--ylabel", default=r"$\Delta G$ (kcal/mol)", help="y-axis label")
    ap.add_argument("--title", default="", help="optional title")
    args = ap.parse_args()
    points = parse_points(args.points)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    width, gap = 0.6, 0.4          # bar length and spacing (x units)
    fig, ax = plt.subplots(figsize=(1.4 * len(points) + 1.5, 3.6))
    span = max(v for _, v in points) - min(v for _, v in points) or 1.0
    for i, (label, g) in enumerate(points):
        x0 = i * (width + gap)
        is_ts = label.upper().startswith("TS")
        color = "tab:red" if is_ts else "black"
        ax.plot([x0, x0 + width], [g, g], color=color, lw=2.5, solid_capstyle="butt")
        dy = 0.04 * span
        ax.text(x0 + width / 2, g + dy if is_ts else g - dy, f"{label}\n{g:.1f}" if is_ts else f"{g:.1f}\n{label}",
                ha="center", va="bottom" if is_ts else "top", fontsize=9, color=color)
        if i > 0:
            xp = (i - 1) * (width + gap) + width
            ax.plot([xp, x0], [points[i - 1][1], g], color="grey", lw=1, ls="--")
    ax.set_xticks([])
    ax.set_ylabel(args.ylabel)
    ax.set_ylim(min(v for _, v in points) - 0.2 * span, max(v for _, v in points) + 0.2 * span)
    for side in ("top", "right", "bottom"):
        ax.spines[side].set_visible(False)
    if args.title:
        ax.set_title(args.title)
    fig.tight_layout()
    fig.savefig(args.output, dpi=300)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
