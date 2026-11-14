# 示例 18-3：C60 的 Cluster 对角化（ELPA1/ELPA2）、DC 对比与并行基准

配套教程：[第18章 并行、大体系 O(N)、作业脚本与常见错误排查](../../tutorials/18-parallel-on-errors.md)（§3、§7）

## 目的

- 用一个规模适中的分子（C60，60 原子，`C6.0-s2p2d1` 共 780 个基函数）练习 Cluster 计算的 MPI × OpenMP 混合并行；
- 比较 `scf.eigen.lib elpa1` 与 `elpa2` 的耗时（两者结果应一致到 SCF 收敛精度）；
- 把 O(N) DC 方法（`orderN.HoppingRanges 7.0`）的总能与常规对角化比较，直观认识截断误差；
- 作为 `bench_np_nt.sh` 的测试输入，找到本机最快的“进程数 × 线程数”组合。

## 文件说明

| 文件 | System.Name | 说明 |
|---|---|---|
| `C60_cluster_elpa.dat` | `c60_elpa1` | `scf.EigenvalueSolver Cluster`，`scf.eigen.lib elpa1`（Cluster 的默认值） |
| `C60_DC.dat` | `c60_dc` | 同一结构与参数，改用 `DC`、`orderN.HoppingRanges 7.0` |
| `run_compare.sh` | — | 依次运行 ELPA1、ELPA2（自动生成 `C60_cluster_elpa2.dat`）、DC，打印 Utot、SCF 步数、耗时和相对 ELPA1 的能量差（mHartree） |

结构为理想截角二十面体（所有 C–C 键长 1.44 Å，未优化），放在 18 Å 立方盒子中；GGA-PBE，`scf.energycutoff 200`，`scf.criterion 1.0e-8`，`level.of.fileout 0`。这个结构只用于基准测试，真实 C60 有 1.40/1.45 Å 两种键长，做物性计算前请先优化（[第6章](../../tutorials/06-optimization-eos-neb.md)）。

## 运行步骤

```bash
cd OpenMX/examples/18-c60-cluster-elpa
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" *.dat

# ELPA1 vs ELPA2 vs DC
NP=8 NT=1 bash run_compare.sh

# MPI x OpenMP scan on one 16-core node (short test: set scf.maxIter 10 in a copy first)
sed 's/^scf.maxIter .*/scf.maxIter                 10/' C60_cluster_elpa.dat > c60_short.dat
bash ../18-job-scripts/bench_np_nt.sh c60_short.dat "16x1 8x2 4x4 2x8"
```

## 预期结果 / 如何检查

```bash
grep "normally finished" c60_*.std
grep "Utot\." c60_elpa1.out c60_elpa2.out c60_dc.out
grep -A3 "Total Computational Time" c60_elpa1.out
```

- `c60_elpa1` 与 `c60_elpa2` 的 `Utot` 应在 1e-8 Hartree 量级内一致；手册（第 28.1.2 节）指出两者在 Cluster 计算中速度几乎相同，你的机器上可能略有差别。
- `c60_dc` 的 `Utot` 与常规对角化相差 mHartree 量级。作为参照，手册表 5 用最小基组算得 C60 常规对角化 −343.90129 Hartree、DC(7.0 Å) −343.90004 Hartree（差 1.25 mHartree）；本例的基组不同，绝对值不可直接比较，只比较差值的量级。
- 对 60 个原子的体系，DC 不会比常规对角化快——手册估计 DC 的优势要在约 500 原子以上才出现。本例的目的是看误差，而不是看速度。
- 基准扫描中，同一核数下不同组合的耗时差异可能达到数十个百分点；所有组合的 `Utot` 必须一致。

## 注意事项

- Cluster 计算的 MPI 并行先分配到自旋，再分配到本征态；非自旋极化的小分子用太多进程收益有限，可以多用 OpenMP 线程。
- 在共享节点上做基准时，其它用户的作业会干扰计时，最好独占节点（Slurm `--exclusive`）。
- 想看更真实的大分子基准，可以跑官方的 `-runtestL`（含 Mn12、GEL、SialicAcid 等），见教程 §6.7。

返回：[第18章](../../tutorials/18-parallel-on-errors.md) | [OpenMX 目录](../../README.md)
