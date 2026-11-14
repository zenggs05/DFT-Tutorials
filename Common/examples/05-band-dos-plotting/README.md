# 05-band-dos-plotting：GaAs 能带与态密度——从 VASP 输入到出版级作图

## 目的

用一个两原子的 GaAs 原胞走完“SCF → 能带（ICHARG=11 + line-mode）→ DOS（ISMEAR=-5）→ 作图”的完整流程，并用四种工具画同一组数据，便于对照：

- **sumo**（命令行，出版级 PDF）：`sumo_plot.sh`；
- **PyProcar 6**（投影能带、PDOS、能带+DOS 组合、费米面）：`pyprocar_fatband.py`；
- **py4vasp**（读 `vaspout.h5`，交互式图与 CSV 导出）：`py4vasp_demo.ipynb`；
- **matplotlib + pymatgen**（完全自定义，配套样式文件）：`plot_band_dos_matplotlib.py` + `dft_paper.mplstyle`。

另附 HSE06 零权重 k 点能带的输入（`INCAR_hse_band` + `KPOINTS_hse_zero_weight`）。

## 文件说明

| 文件 | 作用 |
|---|---|
| `POSCAR` | GaAs 闪锌矿原胞，a = 5.653 Å（实验晶格常数） |
| `INCAR_scf` / `KPOINTS_scf` | 自洽计算，8×8×8 Γ 网格，写出 CHGCAR |
| `INCAR_band` / `KPOINTS_band_linemode` | 非自洽能带，`ICHARG = 11`，路径 L–Γ–X–U\|K–Γ，`LORBIT = 11` |
| `INCAR_dos` / `KPOINTS_dos` | 非自洽 DOS，16×16×16 网格，`ISMEAR = -5`、`NEDOS = 3001`、`LORBIT = 11` |
| `INCAR_hse_band` / `KPOINTS_hse_zero_weight` | HSE06 自洽能带：6×6×6 网格的 16 个不可约 k 点（带权重）+ 21 个零权重路径点 L–Γ–X |
| `job_band_dos.slurm` | Slurm 脚本：依次在 `scf/`、`band/`、`dos/`（可选 `hse_band/`）中运行 |
| `sumo_plot.sh` | sumo-bandplot / sumo-dosplot / sumo-bandstats 一键作图 |
| `pyprocar_fatband.py` | PyProcar 6 的 plain/parametric 能带、分元素 PDOS、能带+DOS 组合图、可选三维费米面 |
| `py4vasp_demo.ipynb` | py4vasp 读取 `vaspout.h5` 画能带、投影能带、PDOS，导出 CSV |
| `plot_band_dos_matplotlib.py` | pymatgen 读取 vasprun.xml，matplotlib 画能带+PDOS（支持自旋极化、分元素 fat band） |
| `dft_paper.mplstyle` | 期刊单栏尺寸、嵌入字体、Okabe–Ito 色盲友好配色的 matplotlib 样式文件 |

## 运行步骤

```bash
# 0) 准备 POTCAR（不随仓库提供）：PBE 的 Ga_d 与 As，顺序与 POSCAR 一致
cat $VASP_PP_PATH/potpaw_PBE/Ga_d/POTCAR $VASP_PP_PATH/potpaw_PBE/As/POTCAR > POTCAR
grep TITEL POTCAR                       # 应依次看到 Ga_d、As

# 1) 计算（约 16 核数分钟；加 RUN_HSE=1 追加 HSE06 能带，耗时多一个数量级）
sbatch job_band_dos.slurm
# RUN_HSE=1 sbatch --export=ALL job_band_dos.slurm

# 2) 参考能量：取 SCF 或 DOS 计算的 E-fermi，不要用能带计算的
grep E-fermi scf/OUTCAR dos/OUTCAR

# 3) 作图（任选）
bash sumo_plot.sh
python pyprocar_fatband.py --fermi <E-fermi> --mode parametric --atoms 1 --orbitals 1 2 3
python plot_band_dos_matplotlib.py --band band --dos dos --project Ga,As --formats pdf png
jupyter lab py4vasp_demo.ipynb          # 需要 VASP 以 -DVASP_HDF5 编译
```

## 预期结果 / 如何检查

- `band/OUTCAR`、`dos/OUTCAR` 末尾有 `General timing and accounting informations for this job`；
- `band/vasprun.xml` 能被 pymatgen 完整读取（作业中断时 XML 不完整会报 ParseError）：`python -c "from pymatgen.io.vasp import Vasprun; print(Vasprun('band/vasprun.xml', parse_potcar_file=False).get_band_structure('band/KPOINTS', line_mode=True).get_band_gap())"`；
- PBE 下 GaAs 是 Γ 点直接带隙，但数值明显低于实验值 1.42 eV（低温约 1.52 eV），甚至可能接近 0，这是半局域泛函的已知缺陷，不是作图出错；HSE06 会明显改善。带隙的具体数值以你自己的计算为准；
- 四种工具画出的能带形状应一致；**零点不同**：sumo 默认取 VBM，VASPKIT/py4vasp/本例 matplotlib 脚本默认取 E-fermi，对比时注意。

## 注意事项

- 能带与 DOS 两步计算的 `ENCUT`、`PREC`、POTCAR 必须与 SCF 相同；`LCHARG = .FALSE.`，避免覆盖 SCF 的 CHGCAR。
- `ISMEAR = -5` 需要 Γ 中心且足够密的网格；不可约 k 点少于 4 个时 VASP 会报错停止，此时改用 `ISMEAR = 0` 和较小的 `SIGMA`。
- `KPOINTS_hse_zero_weight` 中的不可约 k 点是按本例 POSCAR 的晶胞取向、6×6×6 Γ 网格和 Td 点群加时间反演算出的；**换了结构或网格就必须重新生成**（VASPKIT 251、`sumo-kgen --hybrid` 或先跑一次 PBE 取 IBZKPT），不要手工改。带权重部分已按对称性核对（16 个点的权重之和为 216 = 6³）。这个文件的第一行不是 VASPKIT 251 的格式，因此 HSE 结果请用 sumo-bandplot 或 py4vasp 作图；VASPKIT 252–257 只支持由 251 生成的 KPOINTS。VASP ≥ 6.3 也可以改用 `KPOINTS_OPT`，见 [VASP 第6章](../../../VASP/tutorials/06-hybrid-metagga.md)。
- 程序输出的具体数值依赖版本与赝势，本目录不提供“标准答案”。

返回教程：[第5章 后处理与可视化](../../tutorials/05-postprocess-visualization.md)
