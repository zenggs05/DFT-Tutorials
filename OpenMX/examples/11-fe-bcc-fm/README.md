# 算例 11-fe-bcc-fm：bcc Fe 铁磁共线自旋极化计算（与非磁对照）

## 目的

- 练习 `scf.SpinPolarization on` 与坐标行第 6、7 列（上/下自旋初始电子数）的写法；
- 从 `.out` 读出总自旋磁矩、Mulliken 与 Voronoi 原子磁矩；
- 用同一组参数比较铁磁（FM）与非磁（NM）两种状态的总能量。

对应教程：[第11章 自旋极化、非共线磁性与交换耦合](../../tutorials/11-spin-noncollinear.md) §3。

## 文件说明

| 文件 | 说明 |
|---|---|
| `Fe_bcc_FM.dat` | bcc Fe 原胞（a = 2.866 Å，实验值），Fe5.5H-s3p2d1 / Fe_PBE19H（16 价电子），初始 10.0 / 6.0，k 网格 16³ |
| `job.slurm` | 先算 FM，再用 sed 生成并运行 NM 对照（`scf.SpinPolarization off`、8.0 / 8.0），最后打印能量与磁矩 |

运行前修改 `DATA.PATH`。

## 运行步骤

```bash
sbatch job.slurm
# or by hand:
mpirun -np 8 openmx Fe_bcc_FM.dat -nt 1 > fe_fm.std 2>&1
```

## 预期结果 / 如何检查

```bash
grep "Utot\." fe_fm.out fe_nm.out          # E(FM) should be clearly lower than E(NM)
grep "Total spin moment" fe_fm.out         # total spin moment per cell (1 atom)
grep -A 6 "Mulliken populations" fe_fm.out # Up / Down / Sum / Diff per atom
grep "by Voronoi charges" fe_fm.out        # Voronoi spin moment
```

- 实验饱和磁矩约 2.2 μB/Fe；GGA 通常给出 2.2–2.3 μB 左右。Mulliken 的 `Diff` 在多组 d 基函数时可能略大，Voronoi 值可作对照；
- `Sum` 列应为 16.000（Fe_PBE19H 的价电子数），`Sum of MulP` 的 total 应等于 ideal(neutral)；
- FM 比 NM 低数百 meV/atom 的量级是合理的；具体数值以你的计算为准；
- `fe_fm.sden.cube` 是自旋密度，可在 VESTA 中以 ±0.005 左右的等值面查看。

## 注意事项

- 输出的磁矩采用 OpenMX 约定（与电子的物理磁矩符号相反），对大小与 FM/AFM 判断没有影响。
- k 点与电子温度会影响磁矩：建议把 `scf.Kgrid` 从 12³ 加密到 20³ 检查收敛，电子温度不要设得过高。
- 若磁矩塌缩为 0，检查初始占据是否对称、混合权重是否过大。
- 块（`<...` 与 `...>` 之间）内不要写注释行；本输入把说明都放在块外。

[返回第11章](../../tutorials/11-spin-noncollinear.md)
