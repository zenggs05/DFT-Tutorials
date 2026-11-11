# 算例 10-benzene-mo：苯分子的 HOMO / LUMO cube 文件

## 目的

- 练习 `MO.fileout`、`num.HOMOs`、`num.LUMOs`、`MO.Nkpoint` 的用法；
- 认识 OpenMX 轨道 cube 的命名规则（`homo{自旋}_{序号}`）和轨道 cube 第 1 行写入的本征值；
- 观察苯的二重简并 π 轨道（HOMO e1g、LUMO e2u），并在 VESTA 中画出正负相位的等值面。

对应教程：[第10章 电荷分析、分子轨道与 cube 可视化](../../tutorials/10-charge-mo-cube.md) §5。

## 文件说明

| 文件 | 说明 |
|---|---|
| `Benzene_MO.dat` | OpenMX 输入：GGA-PBE，C6.0-s2p2d1 / H6.0-s2p1，Cluster 求解器，输出 2 个 HOMO 和 2 个 LUMO |
| `job.slurm` | 最简 Slurm 脚本，运行后列出轨道文件、HOMO 序号与本征值 |

运行前修改 `DATA.PATH`（指向包含 `PAO/` 与 `VPS/` 的目录）。

## 运行步骤

```bash
sbatch job.slurm
# or
mpirun -np 4 openmx Benzene_MO.dat -nt 1 > benzene.std 2>&1
```

## 预期结果 / 如何检查

```bash
ls benzene.homo0_0.cube benzene.homo0_1.cube benzene.lumo0_0.cube benzene.lumo0_1.cube
grep "HOMO =" benzene.out                      # expected: HOMO = 15 (30 valence electrons)
grep -A 22 "Eigenvalues (Hartree)" benzene.out # eigenvalues 14/15 and 16/17 nearly degenerate
head -n 2 benzene.homo0_0.cube                 # eigenvalue written in the comment lines
```

- 苯有 6×4 + 6×1 = 30 个价电子，非自旋极化时占据 15 个轨道，故 `HOMO = 15`；
- 第 14、15 个本征值（HOMO−1、HOMO）近似相等，第 16、17 个（LUMO、LUMO+1）也近似相等；
- 只生成 `homo0_*`/`lumo0_*`（自旋 0），因为这是非自旋极化计算；
- HOMO–LUMO 能隙 = (ε16 − ε15) × 27.2114 eV。GGA 的 Kohn–Sham 能隙系统性小于实验能隙，只用于定性比较。

在 VESTA 中打开 `benzene.homo0_0.cube`，Isosurface level 取 0.02–0.05，显示正、负两种颜色：应看到垂直于分子平面、上下反号的 π 型轨道。两条简并轨道的形状可能相对教科书图形转过一个角度——它们的任意线性组合都是本征态。

## 注意事项

- `MO.Nkpoint` 在 Cluster 计算中必须为 1，否则 OpenMX 报 `MO.Nkpoint should be 1 for the cluster calculation.` 并退出。
- 想看 LCAO 系数（例如验证 HOMO 只由 C 的 pz 组成），把 `level.of.fileout` 改为 2 重算，在 `benzene.out` 中查找 "LCAO coefficients"。
- 几何未优化；需要定量结果时先按第6章做结构优化。

[返回第10章](../../tutorials/10-charge-mo-cube.md)
