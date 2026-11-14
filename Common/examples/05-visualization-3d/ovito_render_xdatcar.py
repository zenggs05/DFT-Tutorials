#!/usr/bin/env python3
"""ovito_render_xdatcar.py - render a VASP XDATCAR (or extxyz/POSCAR) trajectory with OVITO's Python module.

Produces, without any GUI (works on a cluster login node / in a batch job):
  * <prefix>_frame0.png      a still image of the first frame
  * <prefix>.mp4 (or .gif)   an animation of the trajectory (every --every-nth frame)
  * <prefix>_rdf.dat         RDF averaged over the rendered frames (optional, --rdf-cutoff)

Usage:
  python ovito_render_xdatcar.py XDATCAR
  python ovito_render_xdatcar.py XDATCAR --every-nth 5 --movie traj.gif --size 800 600
  python ovito_render_xdatcar.py XDATCAR --rdf-cutoff 6.0 --no-movie

Installation (OVITO 3.16.1+ requires Python >= 3.11; OVITO advises NOT to pip-install
into an Anaconda/conda environment - use a plain venv instead):
  python3.12 -m venv ~/venv-ovito && source ~/venv-ovito/bin/activate
  pip install ovito
"""
import argparse
import sys

try:
    from ovito.io import import_file
    from ovito.modifiers import CoordinationAnalysisModifier
    from ovito.vis import TachyonRenderer, Viewport
except ImportError:  # pragma: no cover
    sys.exit("OVITO Python module not found: pip install ovito (Python >= 3.11, plain venv)")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("traj", help="XDATCAR, OUTCAR-derived extxyz, POSCAR, ... (format detected by OVITO)")
    p.add_argument("--prefix", default="traj", help="prefix of the output files")
    p.add_argument("--size", type=int, nargs=2, default=(800, 600), metavar=("W", "H"))
    p.add_argument("--camera", type=float, nargs=3, default=(2.0, 1.0, -1.0),
                   help="camera viewing direction (x y z)")
    p.add_argument("--ortho", action="store_true", help="orthographic instead of perspective view")
    p.add_argument("--every-nth", type=int, default=1, help="render every n-th frame in the movie")
    p.add_argument("--fps", type=int, default=10)
    p.add_argument("--movie", default=None, help="movie file name (.mp4/.avi/.mov/.gif); default <prefix>.mp4")
    p.add_argument("--no-movie", action="store_true", help="only the still image (and RDF)")
    p.add_argument("--rdf-cutoff", type=float, default=0.0, help="compute an averaged RDF up to this cutoff (A)")
    p.add_argument("--rdf-bins", type=int, default=200)
    return p.parse_args()


def main():
    a = parse_args()
    pipeline = import_file(a.traj)
    nframes = pipeline.source.num_frames
    print(f"{a.traj}: {nframes} frame(s)")

    data0 = pipeline.compute(0)
    lengths = [round(float(sum(x * x for x in data0.cell[:, i]) ** 0.5), 3) for i in range(3)]
    print("particles:", data0.particles.count, " cell lengths (A):", lengths)

    pipeline.add_to_scene()
    pipeline.source.data.cell.vis.line_width = 0.05   # thin cell edges in rendered images

    vp = Viewport(type=Viewport.Type.Ortho if a.ortho else Viewport.Type.Perspective,
                  camera_dir=tuple(a.camera))
    vp.zoom_all(size=tuple(a.size))
    renderer = TachyonRenderer()          # software ray tracer: no GPU / display needed

    still = f"{a.prefix}_frame0.png"
    vp.render_image(size=tuple(a.size), frame=0, filename=still, background=(1, 1, 1), renderer=renderer)
    print("written", still)

    if not a.no_movie and nframes > 1:
        movie = a.movie or f"{a.prefix}.mp4"
        vp.render_anim(filename=movie, size=tuple(a.size), fps=a.fps, every_nth=a.every_nth,
                       background=(1, 1, 1), renderer=renderer)
        print("written", movie)

    if a.rdf_cutoff > 0:
        pipeline.modifiers.append(CoordinationAnalysisModifier(cutoff=a.rdf_cutoff,
                                                               number_of_bins=a.rdf_bins))
        acc, r, n = None, None, 0
        for frame in range(0, nframes, a.every_nth):
            data = pipeline.compute(frame)
            xy = data.tables["coordination-rdf"].xy()
            r = xy[:, 0]
            acc = xy[:, 1].copy() if acc is None else acc + xy[:, 1]
            n += 1
        out = f"{a.prefix}_rdf.dat"
        with open(out, "w") as fh:
            fh.write(f"# total g(r) averaged over {n} frames, cutoff {a.rdf_cutoff} A\n# r(A)  g(r)\n")
            for ri, gi in zip(r, acc / n):
                fh.write(f"{ri:10.4f} {gi:12.6f}\n")
        print("written", out)


if __name__ == "__main__":
    main()
