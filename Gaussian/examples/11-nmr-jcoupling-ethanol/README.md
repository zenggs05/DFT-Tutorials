# 示例 11-nmr-jcoupling-ethanol：乙醇的自旋-自旋耦合常数

配套教程：[第11章 NMR 化学位移与其他分子性质](../../tutorials/11-nmr-props.md)（第 7 节）。

## 目的

用 `NMR=(SpinSpin,Mixed)` 计算乙醇（anti 构象）的屏蔽常数和全部核对之间的 J 耦合常数，练习读取 K/J 矩阵、对甲基做旋转平均，并理解同位素对 J 的影响。

## 文件说明

| 文件 | 内容 |
|---|---|
| `ethanol_nmr_spinspin_mixed.gjf` | 第 1 步 B3LYP-D3(BJ)/6-311+G(d,p) Opt+Freq；第 2 步 B3LYP/6-311+G(d,p) `NMR=(SpinSpin,Mixed)` |
| `run.slurm` | Slurm 作业脚本 |

原子编号：C1（甲基碳）、C2（亚甲基碳）、O3、H4–H6（甲基）、H7–H8（亚甲基）、H9（羟基）。

## 运行步骤

```bash
sbatch run.slurm
# 或
g16 < ethanol_nmr_spinspin_mixed.gjf > ethanol_nmr_spinspin_mixed.log
```

## 预期结果 / 如何检查

```bash
grep -c "Normal termination" ethanol_nmr_spinspin_mixed.log
grep "NImag" ethanol_nmr_spinspin_mixed.log
# 最后一次出现的 J 矩阵即最终结果（已包含 Mixed 第一步的 FC 项）
awk '/Total nuclear spin-spin coupling J \(Hz\)/{n=NR} {line[NR]=$0} END{for(i=n;i<=n+15;i++) print line[i]}' \
    ethanol_nmr_spinspin_mixed.log
```

矩阵是下三角形式，数值使用 Fortran 的 `D` 指数（如 `0.7123D+01` = 7.123 Hz），超过 5 列时分块打印。

需要做的处理：

- 甲基-亚甲基的 $^3J_{\rm HH}$：取 H4/H5/H6 与 H7/H8 之间 6 个耦合的平均值，再与实验值（约 7 Hz 量级的三重峰/四重峰裂分）比较；
- 同碳（geminal）H7-H8 的 $^2J_{\rm HH}$ 在实验中因两者化学等价而不可观察；
- 羟基 H9 的耦合在多数溶剂中因质子交换而观察不到。

## 13C 耦合怎么算

默认同位素下碳为 12C（核自旋为 0），因此涉及碳的 J 为 0，而同位素无关的 K 矩阵仍有数值。要得到 $^1J_{\rm CH}$ 等：

1. 用 `python ../../scripts/log2xyz.py -h` 查看如何提取优化后的最终坐标；
2. 新建一个不带 `Geom=Check` 的输入，把碳写成 `C(Iso=13)`，例如 `C(Iso=13)  -1.0129  0.1038  0.1169`（坐标替换为优化后的值）；
3. route 使用同样的 `B3LYP/6-311+G(d,p) NMR=(SpinSpin,Mixed)`。

## 注意事项

- SpinSpin 的计算量约为频率计算的 2 倍，`Mixed` 的第一步使用了更大的修改基组；大分子请用 `NMR=(SpinSpin,Mixed,ReadAtoms)` 只算关心的原子（语法见教程第 7.4 节）。
- 乙醇还有 gauche 构象，严格比较需要 Boltzmann 平均。

返回：[第11章 NMR 化学位移与其他分子性质](../../tutorials/11-nmr-props.md)
