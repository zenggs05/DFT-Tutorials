# 示例 06-2：Si 的 RFC5 晶胞优化与金刚石的 EvsLC 状态方程

## 目的

- 用 `MD.Type RFC5` 从扭曲的晶胞同时优化 Si 的晶格矢量和内坐标；
- 用 `MD.Type EvsLC` 扫描金刚石晶格常数（固定 `scf.Ngrid`），用三阶 Birch–Murnaghan 方程拟合 V0、a0、B0、B0'；
- 提供可复用的拟合脚本 `fit_evslc.py`（第4章基组测试和石墨 c 轴扫描也用它）。

对应教程：[第6章 结构优化、晶胞优化、状态方程与 NEB](../../tutorials/06-optimization-eos-neb.md)（§4、§5）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `Si_RFC5_cellopt.dat` | `System.Name si_rfc5`；初始晶格三个矢量长度和夹角都偏离 fcc 原胞，原子 2 位于 (0.27, 0.24, 0.26)；RFC5，判据 1.0e-4 Hartree/Bohr |
| `Cdia_EvsLC.dat` | `System.Name cdia_evslc`；从 a = 3.45 Å 起，每步 +0.5%，共 15 步；`scf.Ngrid 24 24 24`，12×12×12 k 点 |
| `fit_evslc.py` | 读取 `.EvsLC`（9 列晶格矢量 + 总能），`--mode volume` 做 BM3 拟合，`--mode c` 做 E(c) 多项式拟合；只依赖 numpy，`--plot` 需要 matplotlib |
| `job.slurm` | Slurm 作业脚本 |

## 运行步骤

```bash
cd OpenMX/examples/06-cellopt-eos-si-c
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" *.dat   # once

# variable-cell optimization
mpirun -np 8 openmx Si_RFC5_cellopt.dat -nt 1 > si_rfc5.std
grep -A25 "History of cell optimization" si_rfc5.out

# E-V curve of diamond
mpirun -np 8 openmx Cdia_EvsLC.dat -nt 1 > cdia_evslc.std
python fit_evslc.py cdia_evslc.EvsLC --lattice fcc --natoms 2 --plot cdia_eos.png
python fit_evslc.py -h                       # all options
```

## 预期结果 / 如何检查

**RFC5**：历史表各列为 `MD_iter`、`SD_scaling`、`|Maximum force|`、`Maximum step`、`Utot`、`Enpy`、`Volume`。

- 收敛后三个晶格矢量等长、夹角 60°，原子 2 回到 (¼,¼,¼)；
- 立方晶格常数 $a=(4V)^{1/3}$ 应接近 PBE 文献值约 5.47 Å（实验 5.431 Å）；
- 读取最终晶格：`grep -A4 "<Atoms.UnitVectors" Si_RFC5_cellopt.dat#`，或 `.md` 最后一帧注释行中的 `Cell_Vectors=`。

**EvsLC**：`cdia_evslc.EvsLC` 每行 10 列。拟合输出 V0、a0、E0、B0、B0'、RMS 残差和逐点偏差。

- V0 应位于扫描范围中部，RMS 残差远小于 1 meV/cell；
- 量级参考：实验 a = 3.567 Å、B0 ≈ 443 GPa；PBE 文献值约 3.57–3.58 Å、约 430 GPa（Haas et al., PRB 79, 085104 (2009)）；
- 若数据点呈锯齿，检查 `scf.Ngrid` 是否生效：`grep "Num. of grids" cdia_evslc.out` 应显示 24, 24, 24。

`fit_evslc.py` 已用合成的 Birch–Murnaghan 数据验证：能准确还原输入的 V0、B0、B0'。也可以用 ASE 交叉验证（代码见教程 §5.2）。

## 注意事项

- 晶胞优化只支持共线计算（含 DFT+U）；4.0 手册写明不支持 DFT-D2/D3，需要 vdW 时见示例 [06-vdw-d3-graphite](../06-vdw-d3-graphite/README.md)。
- `RFC5` + `MD.Fixed.Cell.Vectors` 的组合需要 OpenMX ≥ 4.0；3.9 中该约束只对 `OptC1/OptC2/OptC5` 有效。
- 加压（焓）优化：在 `Si_RFC5_cellopt.dat` 末尾取消注释 `MD.applied.pressure 10.0`（GPa，正值为压缩）。
- EvsLC 只均匀缩放晶格，不弛豫内坐标；内坐标可变的晶体请按教程 §5.4 逐点优化。
- 输入文件只含 ASCII 字符和 LF 换行。
