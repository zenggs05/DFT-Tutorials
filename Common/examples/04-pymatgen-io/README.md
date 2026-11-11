# 04-pymatgen-io：pymatgen 输入集、VASP 输出解析、相图与 Gaussian/CP2K 输入

## 目的

- 用 `MPRelaxSet / MPStaticSet / MPNonSCFSet / MPHSEBSSet / MP24RelaxSet / MatPESStaticSet` 一键生成 VASP 输入，并用 `user_incar_settings` 覆盖默认值；
- 用 `Vasprun / Outcar / Oszicar / Locpot` 解析输出，得到带隙、费米能、磁矩、平面平均势，画能带/DOS；
- 通过新版 Materials Project API 构建 Li–Fe–P–O 相图，计算 `E_above_hull` 与形成能；
- 用 `pymatgen.io.gaussian`、`pymatgen.io.cp2k` 生成并修改 Gaussian、CP2K 输入。

## 文件说明

| 文件 | 说明 | 需要 |
|---|---|---|
| `pmg_write_inputset.py` | Si 的 relax → static → band(line) / DOS(uniform) → HSE 能带 → MP24 / MatPES 输入集；默认只写 `POTCAR.spec` | pymatgen |
| `pmg_parse_vasprun.py` | 解析一个 VASP 目录：能量、收敛、带隙、VBM/CBM、磁矩；`--band`、`--dos`、`--locpot` 作图 | pymatgen、matplotlib |
| `pmg_phase_diagram_LiFePO4.py` | MP 新 API 下载 Li–Fe–P–O 条目并缓存为 JSON，计算 LiFePO4 的 E_hull、分解产物，画 Li–Fe–O 三元子相图 | pymatgen、`MP_API_KEY` |
| `pmg_gaussian_input.py` | 写苯的 opt+freq `.gjf`；`--parse` 读取 Gaussian log | pymatgen |
| `pmg_cp2k_input.py` | `StaticSet` 生成水分子 CP2K 输入，再用 `Cp2kInput.update()` 修改；`--parse` 读取 CP2K 输出 | pymatgen |

## 运行步骤

```bash
pip install "pymatgen>=2025"                 # 2026 起核心代码在 pymatgen-core，导入方式不变
python -c "import pymatgen.core as pc; print(pc.__version__)"

python pmg_write_inputset.py                 # 生成 01_relax ... 07_matpes_static
python pmg_parse_vasprun.py path/to/run --dos
python pmg_parse_vasprun.py path/to/band_run --band   # 目录中需有 line-mode KPOINTS

export MP_API_KEY=xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx   # 32 位新 key
python pmg_phase_diagram_LiFePO4.py          # 第一次联网下载，之后可加 --offline

python pmg_gaussian_input.py                 # 写 benzene_opt_freq.gjf
python pmg_cp2k_input.py                     # 写 h2o_static.inp
cp2k.psmp -c h2o_static.inp                  # 只做语法检查
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
- 相图脚本打印 LiFePO4 各多形体的 E_hull（最稳定相应为 0 或接近 0）与分解产物，并生成 `LiFeO_pd.png`。

## 注意事项

- MP 默认（`MPRelaxSet`）开自旋极化（`ISPIN = 2`，`MAGMOM` 每原子 0.6）；Si 这类非磁体系可在 `user_incar_settings` 中设 `ISPIN: 1`。
- `user_incar_settings` 中值为 `None` 表示删除该标签，例如 `{"LREAL": None}`。
- 自己的能量必须与 MP 条目设置一致（同一泛函、U 值、POTCAR 与能量修正）才能放进同一张相图。
- pymatgen 的 CP2K 输入集默认打开 UKS、写 `&DFT/&PRINT/&PDOS`（CP2K 2026.2 重构了该部分），务必检查并修改后再用。

返回教程：第4章 ASE 与 pymatgen 实战
