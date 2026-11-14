# 04-pymatgen-io：pymatgen 输入集、VASP 输出解析、相图与 Gaussian/CP2K 输入

## 目的

- 用 `MPRelaxSet / MPStaticSet / MPNonSCFSet / MPHSEBSSet / MP24RelaxSet / MatPESStaticSet` 一键生成 VASP 输入，并用 `user_incar_settings` 覆盖默认值；
- 用 `Vasprun / Outcar / Oszicar / Locpot` 解析输出，得到带隙、费米能、磁矩、平面平均势，画能带/DOS；
- 通过新版 Materials Project API 构建 Li–Fe–P–O 相图，计算 `E_above_hull` 与形成能；
- 用 `pymatgen.io.gaussian`、`pymatgen.io.cp2k` 生成并修改 Gaussian、CP2K 输入；
- 用 `Chgcar` 做体数据加减（AECCAR0+AECCAR2、差分电荷及其平面积分曲线）；
- 用 `OptimadeRester` 通过 OPTIMADE 同时检索 MP、OQMD、Alexandria 等多个数据库。

## 文件说明

| 文件 | 说明 | 需要 |
|---|---|---|
| `pmg_write_inputset.py` | Si 的 relax → static → band(line) / DOS(uniform) → HSE 能带 → MP24 / MatPES 输入集；默认只写 `POTCAR.spec` | pymatgen |
| `pmg_parse_vasprun.py` | 解析一个 VASP 目录：能量、收敛、带隙、VBM/CBM、磁矩；`--band`、`--dos`、`--locpot` 作图 | pymatgen、matplotlib |
| `pmg_phase_diagram_LiFePO4.py` | MP 新 API 下载 Li–Fe–P–O 条目并缓存为 JSON，计算 LiFePO4 的 E_hull、分解产物，画 Li–Fe–O 三元子相图 | pymatgen、`MP_API_KEY` |
| `pmg_gaussian_input.py` | 写苯的 opt+freq `.gjf`；`--parse` 读取 Gaussian log | pymatgen |
| `pmg_cp2k_input.py` | `StaticSet` 生成水分子 CP2K 输入，再用 `Cp2kInput.update()` 修改；`--parse` 读取 CP2K 输出 | pymatgen |
| `pmg_chgcar_ops.py` | `sum`：AECCAR0+AECCAR2 生成 Bader 参考密度 `CHGCAR_sum`；`diff`：差分电荷 `CHGCAR_diff` 与沿 z 的平面积分曲线 `diff_z.dat` | pymatgen、numpy；可选 matplotlib |
| `pmg_optimade_query.py` | 按元素/化学式/原始 OPTIMADE filter 检索多个数据库，结构写成 CIF | pymatgen、requests（联网，无需 key） |

## 运行步骤

```bash
pip install "pymatgen>=2026.3"               # 2026.3.23 起核心代码在 pymatgen-core，导入方式不变
python -c "import pymatgen.core as pc; print(pc.__version__)"

python pmg_write_inputset.py                 # 生成 01_relax ... 07_matpes_static
python pmg_parse_vasprun.py path/to/run --dos
python pmg_parse_vasprun.py path/to/band_run --band   # 目录中需有 line-mode KPOINTS

export MP_API_KEY=xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx   # 32 位新 key
python pmg_phase_diagram_LiFePO4.py          # 第一次联网下载，之后可加 --offline
python pmg_phase_diagram_LiFePO4.py --offline --my-vasprun myrun/vasprun.xml   # 加入自己的 MP 兼容计算

python pmg_gaussian_input.py                 # 写 benzene_opt_freq.gjf
python pmg_cp2k_input.py                     # 写 h2o_static.inp
cp2k.psmp -c h2o_static.inp                  # 只做语法检查

python pmg_chgcar_ops.py sum AECCAR0 AECCAR2           # -> CHGCAR_sum，再 bader CHGCAR -ref CHGCAR_sum
python pmg_chgcar_ops.py diff AB/CHGCAR A/CHGCAR B/CHGCAR   # -> CHGCAR_diff, diff_z.dat
python pmg_optimade_query.py --providers mp oqmd --formula TiO2 --max 5
```

如果要在 `--potcar` 模式写出真正的 POTCAR，先配置一次赝势库：

```bash
pmg config -p /path/to/vasp_pp  ~/pmg_psp    # /path/to/vasp_pp 下有 potpaw_PBE/ 等解压后的子目录
pmg config --add PMG_VASP_PSP_DIR ~/pmg_psp
pmg config --add PMG_DEFAULT_FUNCTIONAL PBE_54          # 按你手头的赝势版本
```

## 预期结果 / 如何检查

- `pmg_write_inputset.py` 末尾打印各目录的关键 INCAR 标签：`03_band` 中 `ICHARG = 11`，`05_hse_band` 中 `LHFCALC = True`，`06_mp24_relax` 中 `METAGGA = R2scan`、`KSPACING = 0.22`。
- `pmg_parse_vasprun.py` 对截断的 vasprun.xml 会打印 `WARNING: XML is malformed. Parsing has stopped but partial data is available.` 而不是直接崩溃。
- `pmg_chgcar_ops.py diff` 打印的 `net charge of d_rho` 应接近 0；若提示 `FFT grids differ`，说明三次计算的网格不一致，需显式设置相同的 `NGXF/NGYF/NGZF` 重算。
- `pmg_optimade_query.py` 按数据库分组打印条目数与每个结构的化学式、原子数、体积，CIF 写在 `optimade_cif/<provider>/` 下；各数据库的结构来源和泛函不同，不能直接混用能量。
- 相图脚本打印 LiFePO4 各多形体的 E_hull（最稳定相应为 0 或接近 0）与分解产物，并生成 `LiFeO_pd.png`。

## 注意事项

- MP 默认（`MPRelaxSet`）开自旋极化（`ISPIN = 2`，`MAGMOM` 按 YAML 中的元素表给初值，表中没有的元素用 0.6 μB）；Si 这类非磁体系可在 `user_incar_settings` 中设 `ISPIN: 1`。
- `user_incar_settings` 中值为 `None` 表示删除该标签，例如 `{"LREAL": None}`。
- 自己的能量必须与 MP 条目设置一致（同一泛函、U 值、POTCAR 与能量修正）才能放进同一张相图。`--my-vasprun` 会用 `MaterialsProject2020Compatibility` 做与 MP 条目相同的修正，设置不兼容时直接拒绝；`--my-energy` 不做任何修正，只用于粗略测试。
- `pmg_optimade_query.py` 必须给出 `--elements`、`--nelements`、`--formula` 或 `--filter` 之一，空查询会被拒绝（否则相当于下载整个数据库）。
- OPTIMADE 服务由各数据库自行维护，偶尔会超时或暂时下线；某个数据库没有返回结果或报网络错误时，换个时间或减少 `--providers` 再试。
- pymatgen 的 CP2K 输入集默认打开 UKS、写 `&DFT/&PRINT/&PDOS`（CP2K 2026.2 重构了该部分），务必检查并修改后再用。

返回教程：[第4章 ASE 与 pymatgen 实战](../../tutorials/04-ase-pymatgen.md)
