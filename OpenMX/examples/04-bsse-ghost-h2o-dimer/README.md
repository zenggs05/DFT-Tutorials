# 示例 04-2：用空原子 E 做水二聚体的 counterpoise（BSSE）修正

## 目的

- 学会在 `Definition.of.Atomic.Species` 中用 `E` 赝势定义空原子（ghost atom）；
- 用 Boys–Bernardi counterpoise（CP）方法计算水二聚体的结合能，并给出基组重叠误差（BSSE）的大小。

对应教程：[第4章 基组 PAO 与赝势 VPS 的选择](../../tutorials/04-pao-vps.md)（§8）。

## 文件说明

| 文件 | System.Name | 内容 | 对应能量 |
|---|---|---|---|
| `H2O_dimer.dat` | `h2o_dimer` | 完整二聚体 | $E_{AB}^{AB}$ |
| `H2O_dimer_BSSE_ghost.dat` | `h2o_A_ghostB` | 水 A（原子 1–3）+ B 位置的空原子 | $E_{A}^{AB}$ |
| `H2O_dimer_BSSE_ghost2.dat` | `h2o_B_ghostA` | A 位置的空原子 + 水 B（原子 4–6） | $E_{B}^{AB}$ |
| `H2O_monoA.dat` | `h2o_monoA` | 只有水 A | $E_{A}^{A}$ |
| `H2O_monoB.dat` | `h2o_monoB` | 只有水 B | $E_{B}^{B}$ |
| `run_bsse.sh` | — | 依次运行 5 个输入并计算 ΔE_raw、ΔE_CP、BSSE（kcal/mol） | — |
| `job.slurm` | — | Slurm 作业脚本 | — |

几何：两个刚性水分子（r(OH) = 0.9572 Å，∠HOH = 104.52°），供体的 O–H 指向受体 O，R(O···O) = 2.91 Å，受体分子平面垂直于供体平面。这是手工搭建的教学结构，不是 PBE 优化结构。所有输入共用 12 Å 立方盒子、250 Ry 截断能、Standard 基组（`O6.0-s2p2d1`、`H6.0-s2p1`）。

## 运行步骤

```bash
cd OpenMX/examples/04-bsse-ghost-h2o-dimer
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" *.dat   # once

MPI_RUN="mpirun -np 4" NT=1 bash run_bsse.sh    # interactive
sbatch job.slurm run_bsse.sh                    # or on a Slurm cluster
bash run_bsse.sh collect                        # re-print the analysis later
```

## 预期结果 / 如何检查

`run_bsse.sh` 输出格式如下（**示意输出，数值仅供参考**）：

```text
E_AB(AB)  = -34.xxxxxxxxxx Ha
E_A(AB)   = -17.xxxxxxxxxx Ha   E_A(A) = -17.xxxxxxxxxx Ha   BSSE(A) = 0.xxx kcal/mol
E_B(AB)   = -17.xxxxxxxxxx Ha   E_B(B) = -17.xxxxxxxxxx Ha   BSSE(B) = 0.xxx kcal/mol
dE_raw    =   -x.xxx kcal/mol  (no counterpoise)
dE_CP     =   -x.xxx kcal/mol  (counterpoise corrected)
BSSE      =    0.xxx kcal/mol  (dE_CP - dE_raw)
```

检查要点：

1. 每个单体在二聚体基组中的能量应低于（或等于）在自身基组中的能量，即 `BSSE(A)`、`BSSE(B)` ≥ 0；
2. ΔE_CP 的绝对值小于 ΔE_raw，BSSE 为正；Ozaki 的 OpenMX 实践讲义给出的经验是，数据库基组下分子体系的 BSSE 约为 0.5 kcal/mol 量级；
3. 所有 `.md` 文件中都是 `scf_conv=1`；
4. 在 `h2o_A_ghostB.out` 的 Mulliken 布居中，总电子数为 8，空原子上有少量布居，这是 A 借用 B 基函数的直接体现；
5. 量级参考：S22 基准（CCSD(T)）给出的水二聚体结合能约 −5.0 kcal/mol。本例几何未优化，不要求逐位吻合。

## 注意事项

- 空原子的 PAO 必须与对应真实原子完全相同，初始电子数写 `0.0 0.0`。
- `H2O_monoA/B.dat` 中原子相对积分网格的位置与二聚体计算不同，会带来很小的网格误差。要做更精细的比较，可从标准输出 `h2o_dimer.std` 中找到 `Grid_Origin`（`grep Grid_Origin h2o_dimer.std`），在其余输入中加入 `scf.fixed.grid  x y z` 固定网格原点（第5章 §4.4）。CP 修正的三个计算原子位置完全相同，不受此影响。
- 若要比较不同基组下 BSSE 的变化，记得同时修改空原子 `EO`、`EH` 的 PAO。
- 盒子大小只需保证相邻镜像之间基函数不重叠（分子尺寸 + 2 × PAO 半径）。12 Å 对 6.0 bohr 的 PAO 足够。
