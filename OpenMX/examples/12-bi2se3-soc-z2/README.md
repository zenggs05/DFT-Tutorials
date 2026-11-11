# 算例 12-bi2se3-soc-z2：Bi2Se3 的能带反转与 Z2 拓扑不变量

## 目的

- 比较 Bi2Se3 不含 / 含 SOC 的能带，观察 Γ 点附近的能隙变化与价带"驼峰"形状（能带反转的典型特征）；
- 用 `HS.fileout on` 生成 scfout，再用 OpenMX 自带的 `Z2FH`（Fukui–Hatsugai 方法）计算三维 Z2 不变量 (ν0; ν1ν2ν3)；
- （可选）用 `scf.SO.factor` 连续调节 SOC 强度，观察能隙闭合—重开的拓扑相变。

对应教程：[第12章 自旋轨道耦合与拓扑性质](../../tutorials/12-soc-topology.md) §5、§6。

## 文件说明

| 文件 | 说明 |
|---|---|
| `Bi2Se3_SOC_band.dat` | 菱方原胞（5 原子），Bi8.0-s3p2d2f1 / Se7.0-s3p2d2（Standard），GGA-PBE，nc + SOC，k 网格 8³，路径 Z–F–Γ–L；含被注释掉的 `<scf.SO.factor` 块 |
| `Bi2Se3_Z2.dat` | 同一结构的 SCF，`HS.fileout on`，按官方示例使用 `scf.Generation.Kpoint Regular`，不算能带 |
| `Z2FH.in` | Z2FH 的三行输入：半区网格数 5、计算全部 6 个时间反演不变面（1）、从第 1 个面开始（1） |
| `job.slurm` | SOC / 无 SOC 能带 → `soc_splitting.py` 分析 → Z2 SCF → Z2FH |

结构参数（a = 4.138 Å、c = 28.64 Å、u(Bi) = 0.4008、v(Se2) = 0.2117）为文献中常用的近似实验值。运行前修改两个 `.dat` 中的 `DATA.PATH`，并编译 Z2FH：

```bash
cd ~/openmx4.0/source && make Z2FH      # copied to ../work/Z2FH
```

## 运行步骤

```bash
sbatch job.slurm
# or by hand:
mpirun -np 16 openmx Bi2Se3_SOC_band.dat -nt 1 > bi2se3_soc.std 2>&1
python ../12-gaas-soc-band/soc_splitting.py bi2se3_soc.Band --at G Z F L
mpirun -np 16 openmx Bi2Se3_Z2.dat -nt 1 > bi2se3_z2.std 2>&1
mpirun -np 4 ~/openmx4.0/work/Z2FH bi2se3_z2.scfout < Z2FH.in > Z2FH.out
cat Z2.dat
```

## 预期结果 / 如何检查

- 有 SOC 时 Bi2Se3 为体相绝缘体，带隙在 Γ 附近，GGA 约 0.3 eV 量级（实验约 0.3 eV）；SOC 下 Γ 附近的最高价带常呈"M 形/驼峰"，这是能带反转后 Bi-p 与 Se-p 成分在 Γ 处互换的结果；
- 无 SOC 时 Γ 点的能隙与能带形状明显不同（正常次序）；
- `Z2.dat` 应给出强拓扑绝缘体 (ν0; ν1ν2ν3) = (1; 000)（这是 Bi2Se3 公认的拓扑分类，见 Zhang et al., Nat. Phys. 5, 438 (2009)）。手册中 Z2.dat 的格式为：
  - 第 1 行：6 个面的 Z2（x0, xπ, y0, yπ, z0, zπ）；
  - 第 2 行：(nu0, nu1, nu2, nu3)；
- `LCNum*.dat` 与 `LCNum*.pl` 是各面上整数场 n(k) 的数据与 gnuplot 脚本（`gnuplot LCNum1.pl`）；n(k) 本身依赖规范，在不同机器上可能不同，但 Z2 应可重复。

能隙数值依赖结构、基组与 k 点，以你的计算为准。

## 注意事项

- Z2FH 只适用于**非磁**、时间反演对称的非共线计算；磁性体系请改用 Chern 数（calB）或 Wannier 方法。
- Z2 结果对内坐标敏感：官方 work 目录中 `Bi2Se3_trivial-Z2.dat`、`Bi2Se3_weak-Z2.dat` 只改变了原子内坐标就分别得到平庸相与弱拓扑相。请使用可靠的实验或优化结构，并检查带隙确实打开（金属态下 Z2 无定义）。
- 把 `Z2FH.in` 第 1 行从 5 增大到 8 或 10，确认 Z2 不变；
- SOC 强度扫描：取消 `Bi2Se3_SOC_band.dat` 中 `<scf.SO.factor` 块的注释，把 Bi 的 p 因子取 0.0、0.25、0.5、0.75、1.0（其余为 1.0），比较 Γ 点能隙——在某个临界值附近能隙闭合再重新打开。**每个 species 都必须写到 f 为止**。

[返回第12章](../../tutorials/12-soc-topology.md)
