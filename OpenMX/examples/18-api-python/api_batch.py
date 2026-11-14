#!/usr/bin/env python3
"""Run several OpenMX inputs side by side inside ONE MPI job via the OpenMX 4.0 Python API.

    mpirun -np 8 python3 api_batch.py --nt 1 a.dat b.dat c.dat d.dat

MPI_COMM_WORLD is split into one sub-communicator per input file (here 4 groups of
2 ranks) and every group calls pyopenmx_api.run(comm, input, nthreads) once, which is
the pattern of the official example_openmx_api.py.  Useful for many small, independent
calculations (structure screening, finite-displacement sets) that would otherwise
queue as separate jobs.

Requirements (see README.md):
  * pyopenmx_api*.so built with "make clean; make pyapi" in openmx4.0/source and placed
    next to this script or on PYTHONPATH;
  * mpi4py linked against the SAME MPI library that was used to build OpenMX;
  * number of MPI ranks >= number of input files.
Each input file keeps its own System.Name, so use different names for different inputs.
"""
import argparse
import os
import sys


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("inputs", nargs="+", help="OpenMX input files (.dat)")
    p.add_argument("--nt", type=int, default=1, help="OpenMP threads per MPI rank passed to the API")
    args = p.parse_args()

    from mpi4py import MPI          # imported here so that -h works without MPI
    import pyopenmx_api

    world = MPI.COMM_WORLD
    rank, size = world.Get_rank(), world.Get_size()
    ngroup = len(args.inputs)

    if ngroup > size:
        if rank == 0:
            sys.stderr.write("error: {} inputs but only {} MPI ranks; use at least one rank per input\n"
                             .format(ngroup, size))
        return 1
    missing = [f for f in args.inputs if not os.path.isfile(f)]
    if missing:
        if rank == 0:
            sys.stderr.write("error: input file(s) not found: {}\n".format(", ".join(missing)))
        return 1

    color = rank * ngroup // size            # contiguous blocks of ranks per input
    sub = world.Split(color, rank)
    my_input = args.inputs[color]
    if sub.Get_rank() == 0:
        print("group {:3d}: {:3d} ranks x {} threads -> {}".format(color, sub.Get_size(), args.nt, my_input),
              flush=True)

    ierr = pyopenmx_api.run(sub.py2f(), my_input, args.nt)

    # collect the return codes on rank 0
    codes = world.gather((color, my_input, ierr) if sub.Get_rank() == 0 else None, root=0)
    world.Barrier()
    sub.Free()
    status = 0
    if rank == 0:
        for item in sorted(c for c in codes if c is not None):
            print("input {:<30s} ierr = {}".format(item[1], item[2]))
            if item[2] != 0:
                status = 1
    return world.bcast(status, root=0)


if __name__ == "__main__":
    sys.exit(main())
