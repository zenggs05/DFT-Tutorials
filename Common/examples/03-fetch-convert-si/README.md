# 示例 03-fetch-convert-si：下载 Si 结构、检查对称性并转换为四款软件的输入

## 目的

用金刚石结构硅（Si）走一遍“拿到结构 → 检查对称性 → 转成各软件输入”的完整流程：

1. 用 Materials Project **新版 API**（mp-api 或 pymatgen 自带的 `MPRester`）下载 `mp-149`（旧版 REST API 已于 2025-09 下线，旧 key 失效）；
2. 用 spglib 扫描不同 `symprec` 下识别出的空间群，得到原胞和标准惯用胞；
3. 用 ASE 把同一个结构写成 VASP `POSCAR`、xyz/extxyz、Gaussian `.gjf`（周期性体系带 `TV` 行）、CP2K `&CELL/&COORD` 片段和 OpenMX `Atoms.*` 块。

对应教程：[第3章 结构准备](../../tutorials/03-structure-modeling.md) 第 3～5 节。

## 文件说明

| 文件 | 说明 |
|---|---|
| [Si.cif](Si.cif) | 惯用立方胞（8 原子），a = 5.431 Å（室温实验值），故意写成 P1 并逐个列出原子，任何读取器得到的原子都一样 |
| [fetch_mp_structure.py](fetch_mp_structure.py) | 从 Materials Project 下载结构，写出 `<mpid>.cif` 和 `POSCAR_<mpid>`；`--formula` 列出同一化学式的全部条目（按 E_hull 排序）；`--backend mp-api/pymatgen` 选择客户端 |
| [check_symmetry.py](check_symmetry.py) | spglib `symprec` 扫描、原胞/惯用胞、Wyckoff 位置；装了 pymatgen 时再用 `SpacegroupAnalyzer` 给出 `get_primitive_standard_structure()` 结果作对照；`--rattle` 模拟弛豫后带数值噪声的结构 |
| [cif_to_all.py](cif_to_all.py) | 一个结构文件 → `POSCAR`、`.xyz`、`.extxyz`、`.gjf`、`_cp2k.inc`、`_openmx.dat`；`--primitive` 先约化为 spglib 标准原胞 |

## 运行步骤

```bash
conda activate dft                         # Python >= 3.11 environment from Common chapter 1
pip install ase spglib pymatgen mp-api     # pymatgen / mp-api only needed for step 1 and part [4]

# 1) download (needs a NEW 32-character key from the MP dashboard)
export MP_API_KEY=xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
python fetch_mp_structure.py --formula Si           # list candidates, most stable first
python fetch_mp_structure.py mp-149                 # primitive cell -> mp-149.cif, POSCAR_mp-149
python fetch_mp_structure.py mp-149 --conventional --prefix Si_conv

# 2) symmetry
python check_symmetry.py Si.cif
python check_symmetry.py Si.cif --rattle 0.01       # what a slightly noisy CONTCAR looks like

# 3) conversion
python cif_to_all.py Si.cif                          # 8-atom cell -> converted/
python cif_to_all.py Si.cif --primitive -o prim      # 2-atom primitive cell -> prim/
```

没有网络或 API key 时可以跳过第 1 步，直接用仓库里的 `Si.cif` 做第 2、3 步。

## 预期结果 / 如何检查

以下数值由几何关系直接推出，可用来核对你的输出：

- `check_symmetry.py Si.cif`：`symprec` 从 1e-5 到 1e-1 都应识别为 `Fd-3m (227)`；原胞 2 个原子，|a| = |b| = |c| = 5.431/√2 ≈ 3.8403 Å；惯用胞 8 个原子，5.431 Å；只有一组不等价位置（Wyckoff 8a，多重度 8）。
- 加 `--rattle 0.01` 后，`symprec = 1e-5` 时会降到 `P1`，要放宽到与位移相当或更大（1e-2～1e-1 Å）才恢复 `Fd-3m`，这正是第 4 节讲的 `symprec` 取值问题。
- `cif_to_all.py Si.cif --primitive -o prim`：屏幕打印 `Primitive cell: 2 atoms`，体积应为 160.19/4 ≈ 40.05 Å³。

```bash
head -8 converted/POSCAR                 # line 6 "Si", line 7 "8", then "Direct"
grep -c "^TV" converted/Si.gjf           # 3 translation vectors for a 3D-periodic cell
grep -A4 "&CELL" converted/Si_cp2k.inc   # A/B/C in angstrom + PERIODIC XYZ
grep -A3 "<Atoms.UnitVectors" converted/Si_openmx.dat
file converted/*                         # all files must be "ASCII text" (no CRLF)
```

- `fetch_mp_structure.py mp-149`：`Space group : Fd-3m (227)`，`Sites : 2`；MP 的晶格常数来自 PBE（或 r2SCAN）弛豫，比实验值大约 1%，不要和 `Si.cif` 的 5.431 Å 混用。

## 注意事项

- 旧 key（16 位）或旧地址 `materialsproject.org/rest/v2` 会认证失败；pymatgen 的 `MPRester` 在 key 长度不是 32 时直接报 `ValueError: Invalid or old API key ...`。新 key 在 <https://next-gen.materialsproject.org/dashboard> 获取。
- `.xyz` 不保存晶胞，周期体系请用 `.extxyz`、`POSCAR` 或 CIF 传递。
- OpenMX 片段里每个原子的 up/down 电子数按 `DFT_DATA19` 的价电子数一半写入；脚本内置表里没有的元素会写成 `?.?`，用 `--valence Ti=12 O=6` 补上，并核对 `<Definition.of.Atomic.Species` 中的 PAO/VPS 名称。
- Gaussian 周期性计算（PBC）只支持部分方法，`.gjf` 中的 `PBEPBE/6-31G(d)` 只是占位，按 [Gaussian 第3章](../../../Gaussian/tutorials/03-method-basis.md) 选择方法和基组。
- 板块通用转换脚本 [convert_any.py](../../scripts/convert_any.py) 额外支持超胞、加真空、按层固定，用 `python ../../scripts/convert_any.py -h` 查看用法；本目录的 `cif_to_all.py` 不依赖它。

返回：[第3章 结构准备](../../tutorials/03-structure-modeling.md) | [Common 目录](../../README.md)
