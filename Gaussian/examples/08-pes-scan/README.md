# 示例 08：刚性扫描、柔性扫描、二维扫描与 GIC 扫描

## 目的

用小分子演示 Gaussian 16 的各种势能面扫描写法：

- `Scan` + Z-matrix 的刚性扫描（H₂O₂ 二面角）；
- `Opt=ModRedundant` 的柔性扫描（正丁烷二面角；S<sub>N</sub>2 反应的 C···Cl 键长）；
- 两个键长同时扫描的二维势能面；
- `Opt=AddGIC` 沿“键长差”复合坐标扫描。

对应教程：[第8章 势能面扫描（刚性与柔性扫描）](../../tutorials/08-scan.md)。S<sub>N</sub>2 扫描得到的最高点结构可以直接交给 [第7章](../../tutorials/07-ts-irc.md) 做 TS 优化。

## 文件说明

| 文件 | 扫描类型 | 扫描坐标 | 点数 | 理论水平 |
|---|---|---|---|---|
| `h2o2_rigid_zmatrix_scan.gjf` | 刚性（`Scan`） | Z-matrix 变量 `DHOOH 0.0 36 10.0` | 37 | B3LYP/6-311+G(d,p) |
| `butane_dihedral_relaxed_scan.gjf` | 柔性 | `D 1 2 3 4 S 36 10.0`，从 anti（180°）出发 | 37 | B3LYP-D3(BJ)/6-31G(d) |
| `sn2_bond_scan_ts_guess.gjf` | 柔性 | `B 1 2 S 14 -0.1`：C···Cl 3.2 → 1.8 Å | 15 | B3LYP-D3(BJ)/6-31+G(d) |
| `two_dimensional_scan.gjf` | 二维柔性 | `B 1 2 S 6 -0.2`（C···Cl 3.0 → 1.8 Å）× `B 1 3 S 6 0.2`（C–Br 2.0 → 3.2 Å） | 49 | B3LYP-D3(BJ)/6-31+G(d) |
| `gic_bond_difference_scan.gjf` | GIC 柔性 | `RC = (R(1,3) − R(1,2)) × 0.529177`，从约 −1.21 Å 起 12 步 × 0.2 | 13 | B3LYP-D3(BJ)/6-31+G(d) |
| `run_scan.sh` | Slurm 脚本 | 依次运行全部输入，或只运行命令行指定的输入 | — | — |

S<sub>N</sub>2 体系的原子编号：1 = C，2 = Cl，3 = Br，4–6 = H（与第 7 章示例一致）。正丁烷：1–4 为碳链 C1–C4，5–14 为 H。

几何来源：正丁烷坐标由 RDKit 嵌入后经 MMFF 约束优化得到（C1–C2–C3–C4 = 180°）；H₂O₂ 的 Z-matrix 参数和 S<sub>N</sub>2 复合物坐标是手工设定的近似值。刚性扫描之前，最好先在同一水平下优化 H₂O₂，再把优化得到的 ROH、ROO、AOOH 填入输入。

## 运行步骤

```bash
cd Gaussian/examples/08-pes-scan

# 单个运行
g16 < butane_dihedral_relaxed_scan.gjf > butane_dihedral_relaxed_scan.log

# Slurm：全部运行，或只运行其中一个
sbatch run_scan.sh
sbatch run_scan.sh sn2_bond_scan_ts_guess
```

计算量提示：H₂O₂ 刚性扫描最快；正丁烷和一维 S<sub>N</sub>2 扫描需要几十次约束优化；二维扫描需要 49 次约束优化，耗时最长。

## 预期结果 / 如何检查

```bash
# 柔性扫描：收敛的扫描点数应等于总点数（37 / 15 / 49 / 13）
grep -c "Stationary point found" butane_dihedral_relaxed_scan.log

# 柔性扫描汇总表
grep -A 30 "Summary of Optimized Potential Surface Scan" butane_dihedral_relaxed_scan.log

# 刚性扫描汇总表
grep -A 45 "Summary of the potential surface scan" h2o2_rigid_zmatrix_scan.log

# 每个收敛点的能量（相对最低点，kcal/mol）
awk '/SCF Done/{e=$5} /Stationary point found/{n++; printf "%3d %18.8f\n", n, e}' butane_dihedral_relaxed_scan.log > scan_E.dat
awk 'NR==FNR{if(min==""||$2<min)min=$2; next}{printf "%3d %10.3f\n", $1, ($2-min)*627.509}' scan_E.dat scan_E.dat
```

定性预期（文献中的近似范围，仅用于判断曲线形状是否合理）：

- 正丁烷：anti 最低；gauche（约 ±65°）高约 0.6–0.9 kcal/mol；syn（0°）能垒最高，约 4.5–6 kcal/mol；
- H₂O₂：平衡二面角约 110–120°，trans 势垒约 1 kcal/mol，cis 势垒约 7 kcal/mol；
- S<sub>N</sub>2 一维扫描：能量先上升，在 C···Cl ≈ 2.3–2.4 Å 附近达到最高，随后 Br⁻ 离去、能量下降；
- GIC 扫描：输出 `Initial Parameters` 中 RC 的初值应约为 −1.21，此后每个点增加 0.2。

作图：在 GaussView 中打开 `.log`，Results → Scan；或使用仓库脚本 [`scan_irc_profile.py`](../../scripts/scan_irc_profile.py)（`python ../../scripts/scan_irc_profile.py -h` 查看用法）。本仓库不附带输出文件，所有数值以你自己的计算为准。

## 注意事项

- 所有扫描都加了 `NoSymm`：扫描变量一旦破坏对称性，没有 `NoSymm` 的作业会报错。
- ModRedundant 和 GIC 中的原子编号从 1 开始，哑原子不计，与 GaussView 显示的编号一致。修改结构后要重新核对。
- GIC 坐标若是表达式（generic GIC），程序内部以 Bohr 为单位，所以示例中乘以 0.529177 换算为 Å。第一次使用请先把 `NSteps` 改为 2 试跑，确认 RC 的数值和步长符合预期。GIC 需要 G16 Rev. B.01 及以上版本。
- 一定要检查**每个**扫描点是否收敛。没有收敛的点能量不可信。
- 键角扫描不要穿过 180°，否则会触发 `Linear angle in Bend` 等 L103 报错。

[← 返回第8章](../../tutorials/08-scan.md)
