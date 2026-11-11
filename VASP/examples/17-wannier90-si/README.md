# 示例 17-1：Si 的最大局域 Wannier 函数（VASP2WANNIER90）

## 目的

- 用 `LWANNIER90 = .TRUE.` 生成 Wannier90 的输入（`.win/.amn/.mmn/.eig`），再用独立的 `wannier90.x` 构造 8 个 sp3 型最大局域 Wannier 函数；
- 学会设置解纠缠（disentanglement）的外窗口和冻结窗口；
- 输出 Wannier 插值能带和紧束缚哈密顿量 `wannier90_hr.dat`（供 WannierTools 等使用）。

对应教程：[第17章 高级电子结构](../../tutorials/17-advanced-electronic.md)（§3）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `POSCAR` | 金刚石 Si 原胞，a = 5.431 Å（实验值） |
| `INCAR_scf` / `KPOINTS_scf` | 自洽（有对称性，Γ 中心 8×8×8） |
| `INCAR_wannier` / `KPOINTS_wan` | 非自洽 `ICHARG = 11`、`ISYM = -1`、`NBANDS = 16`、`NUM_WANN = 8`、`LWANNIER90 = .TRUE.`；Γ 中心 6×6×6 完整网格 |
| `wannier90.win` | `num_wann = 8`、`Si:sp3` 投影、解纠缠窗口、`bands_plot`、`write_hr`、k 路径 L–Γ–X–U\|K–Γ |
| `job.slurm` | 自洽 → 按 E-fermi 改写窗口 → 非自洽写 Wannier90 文件 → `wannier90.x` |
| `POTCAR` | **不提供**：`cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR` |

## 前提

- VASP 编译时加入 `-DVASP2WANNIER90` 并链接 Wannier90 库（VASP ≥ 6.2 支持 Wannier90 2.x/3.x；Wannier90 3.x 的库需按串行方式编译）。若可执行程序没有接口，运行后不会出现 `wannier90.amn` 等文件；
- 独立的 `wannier90.x`，建议与链接进 VASP 的库版本相同（例如都用 3.1.0）。

## 运行步骤

```bash
cd VASP/examples/17-wannier90-si
cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR
export W90_EXE=/path/to/wannier90.x
sbatch job.slurm
```

手动执行时注意两点：

1. `dis_win_max` / `dis_froz_max` 是**绝对能量**（与 VASP 本征值同一零点），先 `grep E-fermi scf/OUTCAR`，再按 E_F+11 eV / E_F+2 eV 修改；
2. VASP 会把晶胞、原子、`mp_grid` 和 k 点块追加到运行目录里的 `wannier90.win`，因此保留一份原始文件。

```bash
cd wannier
mpirun -np 16 vasp_std > vasp.out
wannier90.x wannier90
```

## 预期结果 / 如何检查

```bash
grep -A 12 "Final State" wannier90.wout         # 8 个 WF 的中心与 spread，以及总 spread
grep "<-- DIS" wannier90.wout | tail -3          # 解纠缠收敛过程
grep "<-- CONV" wannier90.wout | tail -3         # 局域化收敛过程
gnuplot -persist wannier90_band.gnu              # 插值能带
```

- 8 个 Wannier 中心应位于 Si–Si 键附近（sp3 杂化轨道指向近邻），每个 spread 数值相近；
- 冻结窗口以下（价带和导带底）的插值能带应与直接 DFT 能带（[第5章](../../tutorials/05-band-dos.md) 的 Si 能带）重合，冻结窗口以上可能出现偏差，这是解纠缠的正常结果；
- 横坐标的定义可能与 VASP/py4vasp 相差 2π 因子，对比时以高对称点位置对齐；
- `wannier90_hr.dat` 可直接用于 WannierTools、postw90 等后续计算。

## 注意事项

- 报 `dis_windows: More states in the frozen window than target WFs`：冻结窗口太高，降低 `dis_froz_max`；
- 报窗口内态数少于 `num_wann`：`dis_win_max` 太低或 `NBANDS` 不够；
- `ISYM = -1` 使 k 点数变为完整网格（6×6×6 = 216 个），非自洽计算的成本随之增加；
- SOC 计算需用 `vasp_ncl`，`num_wann` 加倍（`spinors = true` 由 VASP 写入），并保持 `ISYM = -1`；
- 也可以不用独立文件，而在 INCAR 中用 `WANNIER90_WIN = "..."` 写入 win 内容（VASP ≥ 6.2，总长度不超过 8192 字符，见 VASP Known issues）。
