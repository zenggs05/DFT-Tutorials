"""Check that charge and multiplicity of a Gaussian input are consistent.
Usage: python check_charge_mult.py job1.gjf [job2.gjf ...]
Handles Cartesian and Z-matrix atom lines, labels such as C1, O(Fragment=1),
C-CT-0.1 and atomic numbers. Ghost atoms (Bq) carry no electrons."""
import re, sys

SYMBOLS = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn "
           "Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd "
           "In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu "
           "Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn").split()
Z = {s.upper(): i + 1 for i, s in enumerate(SYMBOLS)}

def atomic_number(label):
    if "BQ" in label.upper():                      # ghost atom, e.g. O-Bq
        return 0
    head = re.split(r"[-(]", label, maxsplit=1)[0]  # strip MM type / params
    if head.isdigit():
        return int(head)
    m = re.match(r"([A-Za-z]{1,2})", head)
    for cand in (m.group(1), m.group(1)[:1]) if m else ():
        if cand.upper() in Z:
            return Z[cand.upper()]
    return 0                                       # X (dummy) or unknown

def sections(path):
    lines = open(path, encoding="ascii", errors="replace").read().splitlines()
    lines = [l for l in lines if not l.lstrip().startswith("%")]
    blocks, cur = [], []
    for l in lines:
        if l.strip() == "":
            if cur:
                blocks.append(cur)
                cur = []
        else:
            cur.append(l)
    if cur:
        blocks.append(cur)
    return blocks          # [route, title, molecule, extra sections...]

def main(path):
    blocks = sections(path)
    if len(blocks) < 3:                            # e.g. Geom=AllCheck input
        print(f"{path}: no molecule specification found (Geom=AllCheck?) -> skipped")
        return 0
    mol = blocks[2]                                # molecule specification
    nums = [int(x) for x in re.split(r"[\s,]+", mol[0].strip())]
    charge, mult = nums[0], nums[1]                # total charge / multiplicity
    nelec = -charge
    for line in mol[1:]:
        f = line.split()
        if f[0].lower().rstrip(":") in ("variables", "constants"):
            break                                  # Z-matrix variable section
        if len(f) == 2 or ("=" in f[0] and "(" not in f[0]):
            continue                               # "R1 0.97" or "R1=0.97"
        nelec += atomic_number(f[0].split(",")[0])
    ok = nelec % 2 != mult % 2
    print(f"{path}: charge={charge} mult={mult} electrons={nelec} ->",
          "OK" if ok else "IMPOSSIBLE (electron count and multiplicity parity disagree)")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(max(main(p) for p in sys.argv[1:]))
