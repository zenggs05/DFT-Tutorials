# 示例 14-3：主客体复合物 18-冠-6·K⁺ 的 CP 相互作用能与形变能

配套教程：[第14章 弱相互作用与 BSSE（Counterpoise）](../../tutorials/14-weak-bsse.md)（§6）

## 目的

- 处理带电碎片的 Counterpoise 输入（电荷多重度行 `1 1 0 1 1 1`）；
- 掌握主客体结合能的完整流程：复合物优化 → 自由主体优化 → CP 单点 → 形变能 → 结合（自由）能；
- 比较 B3LYP-D3(BJ) 与 ωB97X-D 两种泛函的结果。

## 文件说明

| 文件 | 说明 |
|---|---|
| `host_guest_opt.gjf` | 3 步，B3LYP-D3(BJ)/def2-SVP：① 复合物 Opt Freq（电荷 +1）；② 自由 18-冠-6 Opt Freq（从复合物构型出发）；③ K⁺ 单点 |
| `host_guest_cp_sp.gjf` | 2 步 CP 单点（def2-TZVP）：① B3LYP-D3(BJ)；② ωB97X-D（已含色散，不加 `EmpiricalDispersion`） |
| `run_g16_slurm.sh` | Slurm 作业脚本（16 核、20 GB） |

起始结构用 RDKit 的 MMFF 力场预优化（K···O 2.7–2.9 Å），只作为 DFT 优化的起点。`host_guest_cp_sp.gjf` 中暂时放的也是这套 MMFF 坐标，**正式计算前请替换为 `host_guest_opt.log` 第 1 步的优化坐标**，并保留每个原子的 `(Fragment=n)`：冠醚 42 个原子为碎片 1，K⁺（最后一个原子）为碎片 2。

## 运行步骤

```bash
cd Gaussian/examples/14-cp-crown-k
sbatch run_g16_slurm.sh host_guest_opt.gjf
# after it finishes: extract the optimized complex (step 1) and paste it into host_guest_cp_sp.gjf
python ../../scripts/log2xyz.py -h
sbatch run_g16_slurm.sh host_guest_cp_sp.gjf
```

43 个原子的 def2-TZVP CP 单点每步要做 5 次 SCF，16 核上一般需要数小时；优化作业所需时间取决于起始结构离 DFT 极小点有多远。

## 预期结果 / 如何检查

```bash
grep -c "Normal termination" host_guest_opt.log host_guest_cp_sp.log   # expect 3 and 2
grep "imaginary frequencies" host_guest_opt.log                         # no output = no imaginary modes
python ../14-cp-water-dimer/cp_summary.py host_guest_cp_sp.log
grep "Thermal correction to Gibbs Free Energy" host_guest_opt.log
```

| 检查项 | 预期 |
|---|---|
| 复合物结构 | K⁺ 位于冠醚空腔中心附近，6 个 K–O 距离相近 |
| BSSE | 正值；def2-TZVP 下明显小于双 ζ 基组 |
| ΔE_int(CP) | 较大的负值（以离子-偶极作用为主），两种泛函结果应同号、量级相同 |
| 自由冠醚 | 无虚频；构象可能只是邻近的局部极小（见下） |

形变能需在与 CP 单点相同的水平（def2-TZVP）下计算：对“复合物中的冠醚构型”和“自由冠醚优化构型”各做一次单点，差值即为冠醚形变能；K⁺ 为单原子，形变能为 0。结合自由能的组合方法见教程 §6.4（包括 1 atm→1 M 标准态校正 −1.89 kcal/mol 和溶剂化处理）。

## 注意事项

- 自由 18-冠-6 的全局最低构象与结合 K⁺ 的构象不同，从复合物构型出发的优化通常只得到局部极小，会低估形变能。严格做法是先做构象搜索（如 CREST），见 [第18章](../../tutorials/18-errors-capstone.md)。
- `Counterpoise` 不能与 `SCRF` 合用；溶液中的结合能用“气相 CP 结合能 + SMD 溶剂化能差”的组合（[第9章](../../tutorials/09-solvation.md)）。
- def2 系列对 K 为全电子基组；若换成更重的碱金属（Rb、Cs），def2 会自动带赝势。
- 可对复合物做 IGMH 分析（教程 §9），把冠醚和 K⁺ 各定义为一个片段，直观显示 K–O 作用区域。

返回教程：[第14章](../../tutorials/14-weak-bsse.md)
