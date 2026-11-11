#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gjf_builder.py - generate Gaussian 16 input files (.gjf) from xyz, Gaussian output or existing gjf files.

用任务模板（sp / opt / freq / optfreq / ts / qst2 / qst3 / irc / scan / td / tdopt / nmr / cp / oniom）
从 xyz、Gaussian log 或已有 gjf 批量生成输入文件：自动写 %chk/%mem/%nprocshared、溶剂、色散、
--Link1-- 多步作业，并可从本地基组库自动插入 GenECP 块；生成前检查电荷、多重度与电子数是否匹配。

Usage examples
--------------
  python gjf_builder.py water.xyz -t optfreq                        # B3LYP/def2SVP Opt Freq
  python gjf_builder.py water.xyz -t optfreq -m B3LYP -b "6-31G(d)" --disp GD3BJ --solvent Water
  python gjf_builder.py conf.xyz --all-frames -t opt --mem 16GB --nproc 16     # one gjf per frame
  python gjf_builder.py opt.log -t sp -m M062X -b def2TZVP --disp GD3       # last geometry of a log
  python gjf_builder.py pd.xyz -t opt -b "6-31G(d)" --heavy-basis LANL2DZ   # GenECP (Pd: LANL2DZ+ECP)
  python gjf_builder.py ir.xyz -t sp -b def2TZVP --heavy-basis def2-TZVP --basis-lib ~/basis
  python gjf_builder.py r.xyz -t optfreq --link1 "sp:wB97XD/def2TZVP"       # Opt Freq + high-level SP
  python gjf_builder.py ts_guess.xyz -t ts --link1 irc
  python gjf_builder.py react.xyz -t qst2 --product prod.xyz
  python gjf_builder.py butane.xyz -t scan --scan "D 1 2 3 4 S 36 10.0"
  python gjf_builder.py dye.xyz -t td --nstates 20 -m CAMB3LYP --solvent Ethanol
  python gjf_builder.py dimer.xyz -t cp --fragments auto -b "aug-cc-pVDZ"
  python gjf_builder.py enzyme_model.xyz -t oniom --high 1-12 --low-method UFF   # ONIOM stub

Basis-set library (--basis-lib DIR or environment variable GJF_BASIS_LIB):
  put files downloaded from https://www.basissetexchange.org in "Gaussian" format into DIR, e.g.
  DIR/def2-TZVP.gbs.  --heavy-basis def2-TZVP then takes the orbital basis and (if present) the ECP of
  every heavy element from that file.  --heavy-basis can also be a file path, or one of the built-in
  Gaussian keywords LANL2DZ, LANL2MB, SDD (written as keywords inside the GenECP block).

Notes
-----
* Only Cartesian coordinates are read from gjf files (Z-matrix input is rejected).  Atom labels such
  as C(Fragment=1) and freeze flags (-1/0) in a gjf are kept.
* Charge/multiplicity are taken from the source file when available; otherwise charge 0 and the
  lowest multiplicity consistent with the electron count are used.
* Generated files are ASCII with LF line endings and end with blank lines, as Gaussian requires.
* Checks are guides, not guarantees: read the generated input before submitting it.

Dependencies: Python >= 3.9 standard library only.

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License
"""

import argparse
import os
import re
import sys

SYMBOLS = (
    "X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
    "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce "
    "Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn "
    "Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl "
    "Mc Lv Ts Og"
).split()
Z_OF = {s: i for i, s in enumerate(SYMBOLS) if i > 0}

# approximate single-bond covalent radii (Angstrom, Cordero et al. 2008) for bond detection only
_COV = (
    "0.31 0.28 1.28 0.96 0.84 0.76 0.71 0.66 0.57 0.58 1.66 1.41 1.21 1.11 1.07 1.05 1.02 1.06 2.03 1.76 "
    "1.70 1.60 1.53 1.39 1.39 1.32 1.26 1.24 1.32 1.22 1.22 1.20 1.19 1.20 1.20 1.16 2.20 1.95 1.90 1.75 "
    "1.64 1.54 1.47 1.46 1.42 1.39 1.45 1.44 1.42 1.39 1.39 1.38 1.39 1.40 2.44 2.15 2.07 2.04 2.03 2.01 "
    "1.99 1.98 1.98 1.96 1.94 1.92 1.92 1.89 1.90 1.87 1.87 1.75 1.70 1.62 1.51 1.44 1.41 1.36 1.36 1.32 "
    "1.45 1.46 1.48 1.40 1.50 1.50"
)
COV_RADII = [0.0] + [float(x) for x in _COV.split()]

TASKS = ("sp", "opt", "freq", "optfreq", "ts", "qst2", "qst3", "irc", "scan", "td", "tdopt",
         "nmr", "cp", "oniom")
TASK_ALIASES = {"opt+freq": "optfreq", "opt_freq": "optfreq", "energy": "sp", "nmr-giao": "nmr"}

NO_BASIS_METHODS = {"PM6", "PM7", "PM3", "PM3MM", "AM1", "PDDG", "MNDO", "MINDO3", "CNDO", "INDO",
                    "ZINDO", "DFTB", "DFTBA", "UFF", "DREIDING", "AMBER", "EXTERNAL", "HUCKEL"}
BUILTIN_DISP = {"WB97XD", "B97D", "B97D3", "B2PLYPD", "B2PLYPD3", "MPW2PLYPD", "PW6B95D3", "APFD"}
MINNESOTA_PREFIX = ("M05", "M06", "M08", "M11", "MN12", "MN15", "N12", "SOGGA11")
POPLE_RE = re.compile(r"^(3-21|4-31|6-21|6-31|6-311)", re.IGNORECASE)
KEYWORD_ECP = {"LANL2DZ": 11, "LANL2MB": 11, "SDD": 19}   # keyword -> first Z that carries an ECP
SOLVENT_MODELS = {"smd": "SMD", "pcm": "IEFPCM", "iefpcm": "IEFPCM", "cpcm": "CPCM"}


class BuildError(Exception):
    """Fatal problem in the user input."""


def warn(msg):
    print("WARNING: " + msg, file=sys.stderr)


def note(msg):
    print("NOTE: " + msg, file=sys.stderr)


# --------------------------------------------------------------------------------------
# atoms and readers
# --------------------------------------------------------------------------------------
def element_from_label(label):
    """Atomic number from a Gaussian atom label: C, C1, Cl2, C(Fragment=1), C-CT--0.1, 6, Bq."""
    base = label.split("(")[0].split("-")[0].strip()
    if base.isdigit():
        z = int(base)
        if z >= len(SYMBOLS):
            raise BuildError("atomic number %s out of range in label '%s'" % (base, label))
        return z
    m = re.match(r"([A-Za-z]+)", base)
    if not m:
        raise BuildError("cannot read an element from atom label '%s'" % label)
    letters = m.group(1)
    if letters.lower() == "bq":
        return 0
    if len(letters) >= 2 and letters[:2].capitalize() in Z_OF and (letters[1].islower() or letters.isupper()):
        return Z_OF[letters[:2].capitalize()]
    if letters[0].upper() in Z_OF:
        return Z_OF[letters[0].upper()]
    raise BuildError("unknown element in atom label '%s'" % label)


def make_atom(z, xyz, label=None, flag=None, frag=None):
    return {"z": z, "xyz": [float(v) for v in xyz], "label": label, "flag": flag, "frag": frag}


def read_xyz(path):
    """Return a list of frames: (atoms, comment)."""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        lines = [l.rstrip("\r\n") for l in fh]
    frames = []
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        try:
            nat = int(lines[i].split()[0])
        except ValueError:
            raise BuildError("%s line %d: expected the number of atoms" % (path, i + 1))
        if i + 2 + nat > len(lines):
            raise BuildError("%s: frame at line %d is truncated" % (path, i + 1))
        comment = lines[i + 1]
        atoms = []
        for k in range(nat):
            tok = lines[i + 2 + k].split()
            if len(tok) < 4:
                raise BuildError("%s line %d: expected 'element x y z'" % (path, i + 3 + k))
            atoms.append(make_atom(element_from_label(tok[0]), tok[1:4]))
        frames.append((atoms, comment))
        i += 2 + nat
    if not frames:
        raise BuildError("%s contains no xyz frame" % path)
    return frames


def _is_float(t):
    try:
        float(t)
        return True
    except ValueError:
        return False


def read_gjf(path):
    """Read the first job step of a gjf/com file (Cartesian coordinates only)."""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        text = fh.read().replace("\r\n", "\n").replace("\r", "\n")
    step = re.split(r"(?im)^\s*--link1--\s*$", text)[0]
    lines = [l.split("!")[0].rstrip() for l in step.split("\n")]
    i = 0
    n = len(lines)
    while i < n and (not lines[i].strip() or lines[i].lstrip().startswith("%")):
        i += 1
    if i >= n or not lines[i].lstrip().startswith("#"):
        raise BuildError("%s: route line starting with '#' not found" % path)
    route = []
    while i < n and lines[i].strip():
        route.append(lines[i].strip())
        i += 1
    route = " ".join(route)
    while i < n and not lines[i].strip():
        i += 1
    if re.search(r"(?i)geom\s*=\s*\(?[^ ]*allcheck", route):
        raise BuildError("%s uses Geom=AllCheck (no geometry in the file)" % path)
    title = []
    while i < n and lines[i].strip():
        title.append(lines[i].strip())
        i += 1
    while i < n and not lines[i].strip():
        i += 1
    if i >= n:
        raise BuildError("%s: charge/multiplicity line not found" % path)
    cm = lines[i].replace(",", " ").split()
    try:
        cm = [int(x) for x in cm]
    except ValueError:
        raise BuildError("%s: cannot read the charge/multiplicity line '%s'" % (path, lines[i]))
    if len(cm) < 2 or len(cm) % 2:
        raise BuildError("%s: charge/multiplicity line must contain pairs of integers" % path)
    i += 1
    atoms = []
    while i < n and lines[i].strip():
        tok = lines[i].replace(",", " ").split()
        if len(tok) >= 5 and re.fullmatch(r"[+-]?\d+", tok[1]) and all(_is_float(t) for t in tok[2:5]):
            flag, xyz = int(tok[1]), tok[2:5]
        elif len(tok) >= 4 and all(_is_float(t) for t in tok[1:4]):
            flag, xyz = None, tok[1:4]
        else:
            raise BuildError("%s: '%s' is not a Cartesian atom line (Z-matrix input is not supported; "
                             "convert it with GaussView or newzmat first)" % (path, lines[i].strip()))
        label = tok[0]
        m = re.search(r"(?i)fragment\s*=\s*(\d+)", label)
        atoms.append(make_atom(element_from_label(label), xyz, label=None, flag=flag,
                               frag=int(m.group(1)) if m else None))
        i += 1
    if not atoms:
        raise BuildError("%s: no atoms found" % path)
    return {"atoms": atoms, "charge": cm[0], "mult": cm[1], "cm": cm, "route": route, "title": " ".join(title)}


RE_DASH = re.compile(r"^ -{10,}\s*$")
RE_CHARGE = re.compile(r"Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)")


def read_log(path):
    """Return (frames, charge, mult) from a Gaussian output; frames are orientation blocks."""
    frames = []
    charge = mult = None
    prefer = None
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        it = iter(fh)
        for line in it:
            if "orientation:" in line and ("Input" in line or "Standard" in line or "Z-Matrix" in line):
                kind = "standard" if "Standard" in line else "input"
                atoms, dashes = [], 0
                for row in it:
                    if RE_DASH.match(row):
                        dashes += 1
                        if dashes == 3:
                            break
                        continue
                    if dashes == 2:
                        tok = row.split()
                        try:
                            z = int(tok[1])
                            xyz = [float(t) for t in tok[-3:]]
                        except (ValueError, IndexError):
                            continue
                        if z >= 0:
                            atoms.append(make_atom(z, xyz))
                if atoms:
                    prefer = prefer or kind
                    frames.append((kind, atoms))
            elif charge is None:
                m = RE_CHARGE.search(line)
                if m:
                    charge, mult = int(m.group(1)), int(m.group(2))
    if not frames:
        raise BuildError("%s: no geometry (orientation block) found" % path)
    # keep one orientation per step: Input orientation if the file has it
    kind = "input" if any(k == "input" for k, _ in frames) else "standard"
    return [(a, "") for k, a in frames if k == kind], charge, mult


def load_structures(path):
    """Return (list of (atoms, comment), info dict) for xyz / gjf / com / log / out files."""
    if not os.path.isfile(path):
        raise BuildError("%s not found" % path)
    ext = os.path.splitext(path)[1].lower()
    if ext == ".xyz":
        return read_xyz(path), {}
    if ext in (".gjf", ".com", ".gau", ".inp"):
        g = read_gjf(path)
        return [(g["atoms"], g["title"])], g
    if ext in (".log", ".out"):
        frames, charge, mult = read_log(path)
        info = {}
        if charge is not None:
            info = {"charge": charge, "mult": mult}
        return frames, info
    raise BuildError("%s: unknown file type (use .xyz, .gjf/.com or .log/.out)" % path)


def parse_atom_list(text, natoms, what="atom list"):
    """'1-3,7,9-10' -> sorted 0-based indices."""
    out = set()
    for part in re.split(r"[,\s]+", text.strip()):
        if not part:
            continue
        m = re.fullmatch(r"(\d+)(?:-(\d+))?", part)
        if not m:
            raise BuildError("bad %s '%s' (use e.g. 1-3,7)" % (what, text))
        a = int(m.group(1))
        b = int(m.group(2) or a)
        if a < 1 or b > natoms or a > b:
            raise BuildError("%s '%s' out of range 1..%d" % (what, part, natoms))
        out.update(range(a - 1, b))
    return sorted(out)


# --------------------------------------------------------------------------------------
# chemistry checks
# --------------------------------------------------------------------------------------
def electron_count(atoms, charge):
    return sum(a["z"] for a in atoms) - charge


def check_charge_mult(atoms, charge, mult, label="molecule"):
    """Raise BuildError if multiplicity and electron count are inconsistent."""
    if mult < 1:
        raise BuildError("%s: multiplicity must be >= 1 (got %d)" % (label, mult))
    ne = electron_count(atoms, charge)
    if ne < 0:
        raise BuildError("%s: charge %+d removes more electrons than exist" % (label, charge))
    unpaired = mult - 1
    if unpaired > ne or (ne - unpaired) % 2:
        options = "1, 3, 5" if ne % 2 == 0 else "2, 4, 6"
        raise BuildError("%s: %d electrons (charge %+d) and multiplicity %d are impossible; "
                         "with this charge use multiplicity %s ... (or change the charge)"
                         % (label, ne, charge, mult, options))
    return ne


def default_mult(atoms, charge):
    return 1 if electron_count(atoms, charge) % 2 == 0 else 2


def method_variants(method):
    """Method name without U/R/RO prefix and hyphens, e.g. UM06-2X -> {UM062X, M062X}."""
    m = method.upper().replace("-", "")
    out = {m}
    if m.startswith("RO"):
        out.add(m[2:])
    if m[:1] in ("U", "R"):
        out.add(m[1:])
    return out


def check_dispersion(method, disp, force=False):
    """Return the EmpiricalDispersion keyword (or None) after sanity checks."""
    if not disp or disp.lower() == "none":
        return None
    d = disp.upper()
    if d not in ("GD3BJ", "GD3", "GD2", "PFD"):
        raise BuildError("unknown dispersion '%s' (use GD3BJ, GD3, GD2 or none)" % disp)
    names = method_variants(method)
    if names & BUILTIN_DISP:
        warn("%s already contains dispersion; EmpiricalDispersion=%s not added" % (method, d))
        return None
    if d == "GD3BJ" and any(n.startswith(MINNESOTA_PREFIX) for n in names):
        msg = ("%s has no D3(BJ) parameters in Gaussian (L301 'R6DS8: Unable to choose the S8 parameter'); "
               "use --disp GD3" % method)
        if not force:
            raise BuildError(msg + " (or --no-check to override)")
        warn(msg)
    return "EmpiricalDispersion=" + d


# --------------------------------------------------------------------------------------
# basis-set library (Basis Set Exchange "Gaussian" format)
# --------------------------------------------------------------------------------------
RE_ELEM_HDR = re.compile(r"^\s*([A-Za-z]{1,2})\s+0\s*$")


def parse_gbs(path):
    """Parse a Gaussian-format basis file -> ({elem: basis_lines}, {elem: ecp_lines})."""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        raw = [l.rstrip("\r\n") for l in fh]
    lines = [l for l in raw if not l.lstrip().startswith("!")]
    basis, ecp = {}, {}
    i, n = 0, len(lines)
    while i < n and not RE_ELEM_HDR.match(lines[i]):
        i += 1
    # orbital basis blocks: 'El 0' ... '****'
    while i < n:
        m = RE_ELEM_HDR.match(lines[i])
        if not m:
            if lines[i].strip():
                break
            i += 1
            continue
        el = m.group(1).capitalize()
        block = [lines[i].strip().split()[0].capitalize() + "     0"]
        i += 1
        while i < n and lines[i].strip() != "****":
            block.append(lines[i].rstrip())
            i += 1
        if i >= n:
            raise BuildError("%s: basis block of %s does not end with ****" % (path, el))
        block.append("****")
        basis[el] = block
        i += 1
        # a blank line after '****' separates the orbital basis from the ECP section
        if i < n and not lines[i].strip():
            j = i
            while j < n and not lines[j].strip():
                j += 1
            if j < n and RE_ELEM_HDR.match(lines[j]) and j + 1 < n and "****" not in lines[j + 1] \
                    and re.search(r"(?i)-ecp|\s\d+\s+\d+\s*$", lines[j + 1]):
                i = j
                break
            i = j
    # ECP blocks: 'EL 0' / 'EL-ECP lmax ncore' / potentials ... until the next header or EOF
    cur = None
    while i < n:
        m = RE_ELEM_HDR.match(lines[i])
        if m:
            cur = m.group(1).capitalize()
            ecp[cur] = [cur + "     0"]
        elif cur and lines[i].strip():
            ecp[cur].append(lines[i].rstrip())
        i += 1
    return basis, ecp


def find_basis_file(name, libdir):
    if os.path.isfile(name):
        return name
    if not libdir:
        return None
    if not os.path.isdir(libdir):
        raise BuildError("basis library directory %s not found" % libdir)
    norm = re.sub(r"[^a-z0-9]", "", name.lower())
    for fn in sorted(os.listdir(libdir)):
        stem = fn.split(".")[0]
        if re.sub(r"[^a-z0-9]", "", stem.lower()) == norm and fn.lower().endswith((".gbs", ".bas", ".gen", ".txt")):
            return os.path.join(libdir, fn)
    return None


def build_gen_block(atoms, light_basis, heavy_basis, heavy_elements, libdir):
    """Return (basis_keyword, basis_lines, ecp_lines)."""
    present = sorted({a["z"] for a in atoms if a["z"] > 0})
    if heavy_elements:
        heavy = set()
        for s in re.split(r"[,\s]+", heavy_elements.strip()):
            if s:
                if s.capitalize() not in Z_OF:
                    raise BuildError("unknown element '%s' in --heavy-elements" % s)
                heavy.add(Z_OF[s.capitalize()])
        missing = heavy - set(present)
        if missing:
            warn("--heavy-elements %s not present in the molecule, ignored"
                 % ",".join(SYMBOLS[z] for z in sorted(missing)))
        heavy &= set(present)
    else:
        heavy = {z for z in present if z > 36}
    if not heavy:
        raise BuildError("--heavy-basis given but no heavy element in the molecule "
                         "(default: elements beyond Kr; use --heavy-elements)")
    light = [z for z in present if z not in heavy]
    basis_lines, ecp_lines = [], []
    if light:
        basis_lines += [" ".join(SYMBOLS[z] for z in light) + " 0", light_basis, "****"]
    path = find_basis_file(heavy_basis, libdir)
    if path:
        bas, ecp = parse_gbs(path)
        for z in sorted(heavy):
            el = SYMBOLS[z]
            if el not in bas:
                raise BuildError("element %s not found in basis file %s" % (el, path))
            basis_lines += bas[el]
            if el in ecp:
                ecp_lines += ecp[el]
        note("heavy-element basis taken from %s" % path)
    elif heavy_basis.upper() in KEYWORD_ECP:
        kw = heavy_basis.upper()
        for z in sorted(heavy):
            basis_lines += ["%s 0" % SYMBOLS[z], kw, "****"]
            if z >= KEYWORD_ECP[kw]:
                ecp_lines += ["%s 0" % SYMBOLS[z], kw]
    else:
        raise BuildError("heavy basis '%s' is neither a file, a file in the basis library (--basis-lib / "
                         "GJF_BASIS_LIB) nor one of the keywords %s.  Download it in Gaussian format from "
                         "https://www.basissetexchange.org" % (heavy_basis, ", ".join(sorted(KEYWORD_ECP))))
    keyword = "GenECP" if ecp_lines else "Gen"
    return keyword, basis_lines, ecp_lines


# --------------------------------------------------------------------------------------
# geometry helpers (fragments, ONIOM link atoms)
# --------------------------------------------------------------------------------------
def bonded(a, b, scale=1.2):
    ra = COV_RADII[a["z"]] if 0 < a["z"] < len(COV_RADII) else 1.5
    rb = COV_RADII[b["z"]] if 0 < b["z"] < len(COV_RADII) else 1.5
    d2 = sum((p - q) ** 2 for p, q in zip(a["xyz"], b["xyz"]))
    return d2 < (scale * (ra + rb)) ** 2


def connected_fragments(atoms):
    n = len(atoms)
    label = [-1] * n
    nfrag = 0
    for s in range(n):
        if label[s] >= 0:
            continue
        stack = [s]
        label[s] = nfrag
        while stack:
            i = stack.pop()
            for j in range(n):
                if label[j] < 0 and bonded(atoms[i], atoms[j]):
                    label[j] = nfrag
                    stack.append(j)
        nfrag += 1
    return [l + 1 for l in label]


# --------------------------------------------------------------------------------------
# input assembly
# --------------------------------------------------------------------------------------
def clean_title(text):
    t = re.sub(r"[@#!_\\\-\x00-\x1f]", " ", text)
    t = re.sub(r"\s+", " ", t).strip()
    return t or "Gaussian job"


def fmt_atom(a, with_flag=False, frag=None, oniom=None):
    sym = SYMBOLS[a["z"]] if a["z"] > 0 else "Bq"
    label = "%s(Fragment=%d)" % (sym, frag) if frag else sym
    xyz = "%14.8f %14.8f %14.8f" % tuple(a["xyz"])
    if oniom:
        return "%-16s %2d %s %s" % (label, a["flag"] or 0, xyz, oniom)
    if with_flag:
        return "%-16s %2d %s" % (label, a["flag"] or 0, xyz)
    return "%-16s %s" % (label, xyz)


def scrf_keyword(solvent, model):
    if not solvent or solvent.lower() in ("none", "gas"):
        return None
    m = SOLVENT_MODELS[model.lower()]
    return "SCRF=(%s,Solvent=%s)" % (m, solvent)


def opt_keyword(opts, base=None):
    items = list(base or [])
    if opts.get("modred"):
        items.append("ModRedundant")
    if opts.get("calcfc") and "CalcFC" not in items:
        items.append("CalcFC")
    if opts.get("tight"):
        items.append("Tight")
    if opts.get("maxcycles"):
        items = [i for i in items if not i.lower().startswith("maxcycles")]
        items.append("MaxCycles=%d" % opts["maxcycles"])
    if not items:
        return "Opt"
    if len(items) == 1:
        return "Opt=" + items[0]
    return "Opt=(%s)" % ",".join(items)


def td_keyword(args, root=None):
    items = ["NStates=%d" % args.nstates]
    if root:
        items.append("Root=%d" % root)
    if args.states == "triplets":
        items.append("Triplets")
    elif args.states == "50-50":
        items.append("50-50")
    if args.tda:
        return "TDA=(%s)" % ",".join(items)
    return "TD=(%s)" % ",".join(items)


def task_keywords(task, args, nfrag=0, allow_modred=True):
    o = {"modred": allow_modred and (bool(args.modred) or task == "scan"), "calcfc": args.calcfc,
         "tight": args.tight, "maxcycles": args.maxcycles}
    freq = [] if args.no_freq else ["Freq"]
    if task == "sp":
        return ["SP"]
    if task == "opt":
        return [opt_keyword(o)]
    if task == "freq":
        return ["Freq"]
    if task == "optfreq":
        return [opt_keyword(o), "Freq"]
    if task == "ts":
        o["calcfc"] = True
        o["maxcycles"] = args.maxcycles or 100
        return [opt_keyword(o, ["TS", "CalcFC", "NoEigenTest"])] + freq
    if task in ("qst2", "qst3"):
        o["maxcycles"] = args.maxcycles or 100
        return [opt_keyword(o, [task.upper()])] + freq
    if task == "irc":
        items = ["CalcFC", "MaxPoints=%d" % args.irc_points, "StepSize=%d" % args.irc_step]
        if args.irc_direction != "both":
            items.append(args.irc_direction.capitalize())
        if args.irc_integrator != "default":
            items.append(args.irc_integrator.upper())
        return ["IRC=(%s)" % ",".join(items)]
    if task == "scan":
        return [opt_keyword(o)]
    if task == "td":
        return [td_keyword(args)]
    if task == "tdopt":
        return [opt_keyword(o), "Freq", td_keyword(args, root=args.root or 1)]
    if task == "nmr":
        return ["NMR=%s" % args.nmr_method.upper()]
    if task == "cp":
        return ["Counterpoise=%d" % nfrag]
    if task == "oniom":
        job = {"sp": [], "opt": [opt_keyword(o)], "optfreq": [opt_keyword(o), "Freq"]}[args.oniom_job]
        return job
    raise BuildError("unknown task %s" % task)


def method_spec(method, basis):
    if method.upper().replace("-", "") in NO_BASIS_METHODS or not basis:
        return method
    return "%s/%s" % (method, basis)


def make_link0(args, chk, oldchk=None):
    lines = []
    if chk:
        lines.append("%%chk=%s" % chk)
    if oldchk:
        lines.append("%%oldchk=%s" % oldchk)
    if args.mem and args.mem.lower() != "none":
        lines.append("%%mem=%s" % args.mem)
    nproc = 8 if args.nproc is None else args.nproc
    if args.cpu:
        lines.append("%%cpu=%s" % args.cpu)
    elif nproc > 0:
        lines.append("%%nprocshared=%d" % nproc)
    return lines


def parse_link1(spec, args):
    """'task[:method/basis[:extra words]]' -> dict."""
    parts = spec.split(":", 2)
    task = TASK_ALIASES.get(parts[0].strip().lower(), parts[0].strip().lower())
    if task not in ("sp", "opt", "freq", "optfreq", "td", "nmr", "irc", "tdopt"):
        raise BuildError("--link1 task '%s' not supported (sp, opt, freq, optfreq, td, tdopt, nmr, irc)" % task)
    ms = parts[1].strip() if len(parts) > 1 and parts[1].strip() else None
    extra = parts[2].strip() if len(parts) > 2 else ""
    gas = False
    words = []
    for w in extra.split():
        if w.lower() == "gas":
            gas = True
        else:
            words.append(w)
    return {"task": task, "method_spec": ms, "extra": " ".join(words), "gas": gas}


def build_input(atoms, charge, mult, args, task, chk, title, extra_geoms=None, gjf_info=None):
    """Return the text of one .gjf file."""
    natoms = len(atoms)
    for a in atoms:
        a.setdefault("flag", None)

    # ---- freeze flags
    if args.freeze:
        idx = set(parse_atom_list(args.freeze, natoms, "--freeze"))
        for i, a in enumerate(atoms):
            a["flag"] = -1 if i in idx else 0
    with_flag = any(a.get("flag") is not None for a in atoms)

    # ---- charge / multiplicity
    if not args.no_check:
        check_charge_mult(atoms, charge, mult)

    # ---- fragments for counterpoise
    frags = None
    nfrag = 0
    cm_line = "%d %d" % (charge, mult)
    if task == "cp":
        if args.fragments and args.fragments.lower() != "auto":
            frags = [0] * natoms
            for k, part in enumerate(args.fragments.split("/"), start=1):
                for i in parse_atom_list(part, natoms, "--fragments"):
                    frags[i] = k
            if 0 in frags:
                raise BuildError("atoms %s are in no fragment (--fragments)"
                                 % ",".join(str(i + 1) for i, f in enumerate(frags) if f == 0))
        elif not args.fragments and all(a.get("frag") for a in atoms):
            frags = [a["frag"] for a in atoms]
        else:
            frags = connected_fragments(atoms)
            note("fragments detected from bonding: %s" % ", ".join(
                "%d: atoms %s" % (k, compress([i + 1 for i, f in enumerate(frags) if f == k]))
                for k in sorted(set(frags))))
        nfrag = max(frags)
        if nfrag < 2:
            raise BuildError("counterpoise needs at least 2 fragments (found %d)" % nfrag)
        if args.frag_cm:
            vals = [int(x) for x in args.frag_cm.replace(",", " ").split()]
            if len(vals) != 2 * nfrag:
                raise BuildError("--frag-cm needs %d integers (charge mult for each of %d fragments)"
                                 % (2 * nfrag, nfrag))
        else:
            vals = []
            for k in range(1, nfrag + 1):
                fa = [a for a, f in zip(atoms, frags) if f == k]
                vals += [0, default_mult(fa, 0)]
            if charge != 0 or mult != 1:
                warn("charged/open-shell complex: give fragment charges/multiplicities with --frag-cm")
        for k in range(1, nfrag + 1):
            fa = [a for a, f in zip(atoms, frags) if f == k]
            if not args.no_check:
                check_charge_mult(fa, vals[2 * k - 2], vals[2 * k - 1], "fragment %d" % k)
        if sum(vals[0::2]) != charge:
            warn("fragment charges (%d) do not add up to the total charge (%d)" % (sum(vals[0::2]), charge))
        cm_line = " ".join(str(v) for v in [charge, mult] + vals)

    # ---- ONIOM layers
    oniom_cols = None
    oniom_method = None
    if task == "oniom":
        if not args.high:
            raise BuildError("ONIOM needs --high ATOMS (high-layer atoms, e.g. 1-12)")
        layer = ["L"] * natoms
        for i in parse_atom_list(args.high, natoms, "--high"):
            layer[i] = "H"
        if args.medium:
            for i in parse_atom_list(args.medium, natoms, "--medium"):
                if layer[i] == "H":
                    raise BuildError("atom %d is in both --high and --medium" % (i + 1))
                layer[i] = "M"
        rank = {"H": 3, "M": 2, "L": 1}
        oniom_cols = []
        nlink = 0
        for i, a in enumerate(atoms):
            col = layer[i]
            partners = [j for j in range(natoms) if rank[layer[j]] > rank[layer[i]] and bonded(a, atoms[j])]
            if partners:
                j = min(partners, key=lambda j: -rank[layer[j]])
                col += " H %d" % (j + 1)
                nlink += 1
                if len(partners) > 1:
                    warn("atom %d is bonded to several atoms of higher layers; check its link atom" % (i + 1))
            oniom_cols.append(col)
        hm = method_spec(args.method, args.basis)
        if args.medium:
            oniom_method = "ONIOM(%s:%s:%s)" % (hm, args.medium_method, args.low_method)
            cm_line = " ".join(["%d %d" % (charge, mult)] * 6)
        else:
            oniom_method = "ONIOM(%s:%s)" % (hm, args.low_method)
            cm_line = " ".join(["%d %d" % (charge, mult)] * 3)
        warn("ONIOM stub: %d link atom(s) assigned from bonding (H caps). Check layers, link atoms, MM atom "
             "types and the charge/multiplicity of every layer (all set to %d %d) in GaussView before "
             "running. UFF types are assigned automatically by Gaussian; Amber needs explicit types."
             % (nlink, charge, mult))
        with_flag = True

    # ---- basis / GenECP
    gen_basis, gen_ecp = [], []
    basis_kw = args.basis
    if args.heavy_basis:
        basis_kw, gen_basis, gen_ecp = build_gen_block(atoms, args.basis, args.heavy_basis,
                                                       args.heavy_elements, args.basis_lib)
    elif not args.no_check and POPLE_RE.match(args.basis or "") and any(a["z"] > 36 for a in atoms):
        warn("Pople basis %s does not cover elements beyond Kr (L301 'Atomic number out of range'); "
             "use def2 basis sets or --heavy-basis LANL2DZ/SDD" % args.basis)
    extra_kw = []
    pople = POPLE_RE.match(args.basis or "")
    if gen_basis and pople and pople.group(1) != "6-311":
        extra_kw.append("6D")   # Pople 6-31G* family is defined with Cartesian d; keep that inside Gen

    disp = check_dispersion(args.method, args.disp, force=args.no_check)
    scrf = scrf_keyword(args.solvent, args.solvent_model)

    if task == "oniom":
        model = oniom_method
    else:
        model = method_spec(args.method, basis_kw)
    route = [args.route_prefix, model] + extra_kw
    if disp:
        route.append(disp)
    route += task_keywords(task, args, nfrag)
    if scrf:
        route.append(scrf)
    if args.oldchk:
        route.append("Guess=Read")
    if args.extra:
        route.append(args.extra)

    out = make_link0(args, chk, args.oldchk)
    out += [" ".join(route), "", clean_title(title), "", cm_line]
    for i, a in enumerate(atoms):
        out.append(fmt_atom(a, with_flag=with_flag, frag=frags[i] if frags else None,
                            oniom=oniom_cols[i] if oniom_cols else None))
    out.append("")

    # ---- ModRedundant / scan section
    modred = list(args.modred or [])
    if task == "scan":
        if not args.scan:
            raise BuildError("task scan needs at least one --scan 'B 1 2 S 10 0.1' (or A/D)")
        modred += args.scan
    elif args.scan:
        warn("--scan ignored for task %s" % task)
    if modred and "ModRedundant" not in " ".join(route):
        warn("--modred ignored: task %s has no optimization" % task)
        modred = []
    for line in modred:
        check_modred_line(line, natoms)
    if modred:
        out += [l.strip() for l in modred] + [""]

    # ---- QST2 / QST3 additional geometries
    if task in ("qst2", "qst3"):
        names = ["product"] + (["TS guess"] if task == "qst3" else [])
        for (geo, gtitle), name in zip(extra_geoms or [], names):
            if [a["z"] for a in geo] != [a["z"] for a in atoms]:
                raise BuildError("%s geometry must contain the same atoms in the same order as the reactant"
                                 % name)
            out += [clean_title(gtitle or name), "", "%d %d" % (charge, mult)]
            out += [fmt_atom(a) for a in geo] + [""]

    # ---- Gen / GenECP
    if gen_basis:
        out += gen_basis + [""]
        if gen_ecp:
            out += gen_ecp + [""]

    # ---- --Link1-- steps
    for spec in args.link1 or []:
        st = parse_link1(spec, args)
        ms = st["method_spec"] or model if task != "oniom" else st["method_spec"]
        if not ms:
            raise BuildError("--link1 after an ONIOM job needs an explicit method/basis")
        uses_gen = bool(gen_basis) and (st["method_spec"] is None or re.search(r"(?i)/gen(ecp)?\b", ms))
        if st["method_spec"] and re.search(r"(?i)/gen(ecp)?\b", ms) and not gen_basis:
            raise BuildError("--link1 '%s' uses Gen/GenECP but no --heavy-basis was given" % spec)
        sroute = [args.route_prefix, ms]
        if uses_gen and extra_kw:
            sroute += extra_kw
        if st["method_spec"] is None and disp:
            sroute.append(disp)
        sargs_task = task_keywords(st["task"], args, nfrag, allow_modred=False)
        if st["task"] == "td":
            sargs_task = [td_keyword(args)]
        sroute += sargs_task
        sroute += ["Geom=Check", "Guess=Read"]
        if scrf and not st["gas"]:
            sroute.append(scrf)
        if st["extra"]:
            sroute.append(st["extra"])
        out += ["--Link1--"] + make_link0(args, chk)
        out += [" ".join(sroute), "", clean_title("%s step %s" % (title, st["task"])), "", "%d %d" % (charge, mult), ""]
        if uses_gen:
            out += gen_basis + [""]
            if gen_ecp:
                out += gen_ecp + [""]
        if not chk:
            raise BuildError("--link1 needs a checkpoint file (do not use --no-chk)")
    out.append("")
    text = "\n".join(out) + "\n"
    try:
        text.encode("ascii")
    except UnicodeEncodeError:
        raise BuildError("generated input contains non-ASCII characters (check --title / --extra)")
    return text


def check_modred_line(line, natoms):
    """Validate the atom indices of a ModRedundant line such as 'B 1 2 S 10 0.1' or 'D 1 2 3 4 F'."""
    tok = line.split()
    if not tok:
        raise BuildError("empty ModRedundant line")
    nidx = {"X": 1, "B": 2, "A": 3, "D": 4, "L": 4}
    if tok[0].upper() in nidx:
        idx = tok[1:1 + nidx[tok[0].upper()]]
    else:                          # type inferred from the number of leading integers
        idx = []
        for t in tok:
            if re.fullmatch(r"\d+|\*", t):
                idx.append(t)
            else:
                break
    if not idx:
        raise BuildError("cannot read atom numbers from ModRedundant line '%s'" % line)
    for t in idx:
        if t == "*":
            continue
        if not t.isdigit():
            raise BuildError("ModRedundant line '%s': '%s' is not an atom number" % (line, t))
        if not 1 <= int(t) <= natoms:
            raise BuildError("ModRedundant line '%s' refers to atom %s (molecule has %d atoms)"
                             % (line, t, natoms))


def compress(nums):
    """[1,2,3,5] -> '1-3,5'"""
    out = []
    for n in nums:
        if out and n == out[-1][1] + 1:
            out[-1][1] = n
        else:
            out.append([n, n])
    return ",".join(str(a) if a == b else "%d-%d" % (a, b) for a, b in out)


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------
def build_parser():
    p = argparse.ArgumentParser(
        description="Generate Gaussian 16 .gjf inputs from xyz / log / gjf files with task templates.  "
                    "用任务模板从 xyz/log/gjf 批量生成 Gaussian 输入（含电荷多重度检查、GenECP、Link1）。",
        epilog=__doc__[__doc__.index("Usage examples"):__doc__.index("Notes")],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("inputs", nargs="+", help="xyz, gjf/com or log/out files / 输入结构文件")
    p.add_argument("-t", "--task", default="optfreq",
                   help="task template: %s (default optfreq) / 任务类型" % ", ".join(TASKS))
    p.add_argument("-m", "--method", default="B3LYP", help="method/functional keyword (default B3LYP)")
    p.add_argument("-b", "--basis", default="def2SVP",
                   help="basis keyword (default def2SVP); with --heavy-basis this is the light-atom basis")
    p.add_argument("--disp", default="none", help="empirical dispersion: GD3BJ, GD3, GD2 or none (default none)")
    p.add_argument("--solvent", help="solvent name, e.g. Water, Toluene (default gas phase) / 溶剂")
    p.add_argument("--solvent-model", default="smd", choices=sorted(SOLVENT_MODELS),
                   help="implicit solvent model (default smd; pcm = IEFPCM)")
    p.add_argument("-c", "--charge", type=int, help="total charge (default: from source file or 0)")
    p.add_argument("-s", "--mult", type=int, help="spin multiplicity (default: from source or lowest)")
    p.add_argument("--route-prefix", default="#p", choices=["#", "#p", "#n", "#t", "#P", "#N", "#T"],
                   help="route prefix (default #p, verbose output)")
    p.add_argument("--extra", help="extra route keywords, e.g. 'Int=SuperFine NoSymm' / 额外关键词")
    p.add_argument("--title", help="title line (default: file name and task)")
    p.add_argument("--no-check", action="store_true",
                   help="skip charge/multiplicity and dispersion checks / 跳过检查")

    g = p.add_argument_group("Link 0 / resources")
    g.add_argument("--mem", default="8GB", help="%%mem (default 8GB; 'none' to omit)")
    g.add_argument("--nproc", type=int, help="%%nprocshared (default 8; 0 to omit)")
    g.add_argument("--cpu", help="use %%cpu=LIST (e.g. 0-15) instead of %%nprocshared")
    g.add_argument("--chk", help="checkpoint name (default: output file stem)")
    g.add_argument("--no-chk", action="store_true", help="do not write %%chk")
    g.add_argument("--oldchk", help="read the guess from this chk (%%oldchk + Guess=Read)")

    g = p.add_argument_group("geometry / optimization")
    g.add_argument("--frame", type=int, default=-1,
                   help="frame of a multi-frame xyz/log (1-based, negative from the end; default last)")
    g.add_argument("--all-frames", action="store_true", help="one input per frame (conformers) / 每帧一个输入")
    g.add_argument("--freeze", metavar="ATOMS", help="freeze atoms in Opt, e.g. 1-10,15 (flag -1)")
    g.add_argument("--modred", action="append", metavar="LINE",
                   help="ModRedundant line, repeatable, e.g. 'B 1 2 F' (adds Opt=ModRedundant)")
    g.add_argument("--scan", action="append", metavar="LINE",
                   help="scan coordinate for task scan, repeatable, e.g. 'D 1 2 3 4 S 36 10.0'")
    g.add_argument("--maxcycles", type=int, help="Opt MaxCycles")
    g.add_argument("--calcfc", action="store_true", help="Opt=CalcFC")
    g.add_argument("--tight", action="store_true", help="Opt=Tight")
    g.add_argument("--no-freq", action="store_true", help="ts/qst2/qst3: do not add Freq")

    g = p.add_argument_group("TS / IRC / QST")
    g.add_argument("--product", help="product structure for qst2/qst3 (same atom order)")
    g.add_argument("--ts-guess", help="TS guess structure for qst3")
    g.add_argument("--irc-points", type=int, default=30, help="IRC MaxPoints (default 30)")
    g.add_argument("--irc-step", type=int, default=10, help="IRC StepSize in 0.01 amu^1/2 Bohr (default 10)")
    g.add_argument("--irc-direction", choices=["both", "forward", "reverse"], default="both")
    g.add_argument("--irc-integrator", choices=["default", "lqa", "gs2", "hpc"], default="default",
                   help="IRC integrator (default: Gaussian's default HPC)")

    g = p.add_argument_group("TD-DFT / NMR")
    g.add_argument("--nstates", type=int, default=10, help="TD NStates (default 10)")
    g.add_argument("--root", type=int, help="TD Root (tdopt default 1)")
    g.add_argument("--states", choices=["singlets", "triplets", "50-50"], default="singlets")
    g.add_argument("--tda", action="store_true", help="use TDA instead of TD")
    g.add_argument("--nmr-method", choices=["giao", "csgt"], default="giao")

    g = p.add_argument_group("Counterpoise / ONIOM")
    g.add_argument("--fragments", help="'1-3/4-6' (fragments separated by /) or 'auto' (bonding); "
                                       "default: Fragment= labels of a gjf, else auto")
    g.add_argument("--frag-cm", help="charge and multiplicity of each fragment, e.g. '0 1 0 1'")
    g.add_argument("--high", metavar="ATOMS", help="ONIOM high-layer atoms, e.g. 1-12")
    g.add_argument("--medium", metavar="ATOMS", help="ONIOM medium-layer atoms (3-layer ONIOM)")
    g.add_argument("--low-method", default="UFF", help="ONIOM low-level method (default UFF)")
    g.add_argument("--medium-method", default="PM6", help="ONIOM medium-level method (default PM6)")
    g.add_argument("--oniom-job", choices=["sp", "opt", "optfreq"], default="opt")

    g = p.add_argument_group("Gen / GenECP")
    g.add_argument("--heavy-basis", help="basis (+ECP) for heavy elements: LANL2DZ, LANL2MB, SDD, a file, "
                                         "or a name found in --basis-lib (e.g. def2-TZVP)")
    g.add_argument("--heavy-elements", help="elements that get --heavy-basis, e.g. 'Pd,I' (default: beyond Kr)")
    g.add_argument("--basis-lib", default=os.environ.get("GJF_BASIS_LIB"),
                   help="directory with Gaussian-format basis files (default $GJF_BASIS_LIB)")

    g = p.add_argument_group("multi-step / output")
    g.add_argument("--link1", action="append", metavar="SPEC",
                   help="append a --Link1-- step 'task[:method/basis[:extra words]]', repeatable; "
                        "e.g. 'sp:wB97XD/def2TZVP', 'td:CAMB3LYP/def2TZVP', 'sp:M062X/def2TZVP:gas'")
    g.add_argument("-o", "--output", help="output file (single input / frame only)")
    g.add_argument("--outdir", default=".", help="output directory (default .)")
    g.add_argument("--suffix", help="file name suffix (default _TASK)")
    g.add_argument("--stdout", action="store_true", help="print the input instead of writing a file")
    g.add_argument("-f", "--force", action="store_true", help="overwrite existing files")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    task = TASK_ALIASES.get(args.task.lower(), args.task.lower())
    if task not in TASKS:
        print("ERROR: unknown task '%s'; choose from %s" % (args.task, ", ".join(TASKS)), file=sys.stderr)
        return 2
    if args.cpu and args.nproc is not None:
        print("ERROR: give either --nproc or --cpu, not both", file=sys.stderr)
        return 2
    if args.mem and args.mem.lower() != "none" and not re.fullmatch(r"(?i)\d+(KB|MB|GB|TB|KW|MW|GW|TW)?", args.mem):
        print("ERROR: --mem must look like 8GB, 8000MB or 1GW (integer + unit)", file=sys.stderr)
        return 2
    if task in ("qst2", "qst3") and not args.product:
        print("ERROR: task %s needs --product FILE" % task, file=sys.stderr)
        return 2
    if task in ("qst2", "qst3") and args.modred:
        print("ERROR: --modred with QST2/QST3 is not supported by this builder (edit the input by hand)",
              file=sys.stderr)
        return 2
    if task == "qst3" and not args.ts_guess:
        print("ERROR: task qst3 needs --ts-guess FILE", file=sys.stderr)
        return 2

    jobs = []   # (atoms, charge, mult, stem, title_hint, gjf_info)
    try:
        for path in args.inputs:
            frames, info = load_structures(path)
            if args.all_frames:
                chosen = list(enumerate(frames))
            else:
                k = args.frame - 1 if args.frame > 0 else len(frames) + args.frame
                if not 0 <= k < len(frames):
                    raise BuildError("%s: --frame %d out of range (1..%d)" % (path, args.frame, len(frames)))
                chosen = [(k, frames[k])]
            stem = os.path.splitext(os.path.basename(path))[0]
            for k, (atoms, comment) in chosen:
                charge = args.charge if args.charge is not None else info.get("charge", 0)
                if args.mult is not None:
                    mult = args.mult
                elif "mult" in info and args.charge is None:
                    mult = info["mult"]
                else:
                    mult = default_mult(atoms, charge)
                    if mult == 2:
                        note("%s: odd number of electrons, multiplicity set to 2 (use -s to change)" % stem)
                fstem = "%s_f%03d" % (stem, k + 1) if args.all_frames and len(frames) > 1 else stem
                jobs.append((atoms, charge, mult, fstem, path, info))
        if args.output and len(jobs) > 1:
            raise BuildError("-o/--output can only be used for a single structure; use --outdir")
        extra = []
        if task in ("qst2", "qst3"):
            pf, _ = load_structures(args.product)
            extra.append((pf[-1][0], "product"))
            if task == "qst3":
                tf, _ = load_structures(args.ts_guess)
                extra.append((tf[-1][0], "TS guess"))

        suffix = args.suffix if args.suffix is not None else "_" + task
        written = 0
        for atoms, charge, mult, fstem, path, info in jobs:
            out_path = args.output or os.path.join(args.outdir, fstem + suffix + ".gjf")
            chk = None if args.no_chk else (args.chk or os.path.splitext(os.path.basename(out_path))[0] + ".chk")
            title = args.title or "%s %s %s" % (fstem, task, method_spec(args.method, args.basis))
            text = build_input([dict(a) for a in atoms], charge, mult, args, task, chk, title,
                               extra_geoms=extra, gjf_info=info)
            if args.stdout:
                sys.stdout.write(text)
                continue
            if os.path.exists(out_path) and not args.force:
                raise BuildError("%s exists (use -f to overwrite)" % out_path)
            if os.path.abspath(out_path) == os.path.abspath(path):
                raise BuildError("output %s would overwrite the source file" % out_path)
            os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
            with open(out_path, "w", newline="\n", encoding="ascii") as fh:
                fh.write(text)
            written += 1
            print("written %s  (%d atoms, charge %d, mult %d)" % (out_path, len(atoms), charge, mult),
                  file=sys.stderr)
        if not args.stdout:
            print("%d input file(s) generated" % written, file=sys.stderr)
    except BuildError as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        return 1
    except OSError as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
