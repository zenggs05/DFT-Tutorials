# 算例 09-pt111-slab：Pt(111) p(2×2) 四层 slab 的弛豫（固定底部两层 + 偶极校正）

## 目的

- 认识一个标准金属表面 slab 的全部要素：表面超胞、原子层数、真空层、Selective dynamics 固定底层、Γ 中心 n×n×1 k 网格；
- 只弛豫表面两层原子（`ISIF = 2`），并在 z 方向加偶极校正；
- 为后续功函数（[09-pt111-workfunction-dipole](../09-pt111-workfunction-dipole/README.md)）和吸附计算（[第10章](../../tutorials/10-adsorption-catalysis.md)）准备弛豫好的表面。

对应教程：[第9章 表面模型、表面能、功函数与偶极校正](../../tutorials/09-surface-workfunction.md) §3–§5。

## 文件说明

| 文件 | 说明 |
|---|---|
| `POSCAR` | Pt(111) p(2×2)，4 层共 16 原子，a(体相) = 3.97 Å，上下各约 10 Å 真空（c = 26.88 Å），由 ASE `fcc111` 生成；底部两层 `F F F`，上部两层 `T T T` |
| `KPOINTS` | Γ 中心 8×8×1 |
| `INCAR` | `ISIF=2` 弛豫，MP smearing，`LDIPOL/IDIPOL=3/DIPOL=0.5 0.5 0.5` |
| `job.slurm` | 最简 Slurm 脚本 |

POTCAR（VASP 推荐的 Pt 势就是标准 `Pt`）：

```bash
cat $VASP_PP_PATH/potpaw_PBE/Pt/POTCAR > POTCAR
```

> [!IMPORTANT]
> POSCAR 中的 3.97 Å 是一个接近 PBE 的近似晶格常数（实验约 3.92 Å）。正式计算前请用 [第4章](../../tutorials/04-relax-eos.md) 的方法，以**相同的 ENCUT、POTCAR 和 smearing** 优化体相 Pt，再用你自己的晶格常数重新生成 slab（见教程 §3 的 ASE / pymatgen 脚本）。

## 运行步骤

```bash
sbatch job.slurm
```

## 预期结果 / 如何检查

```bash
grep "reached required accuracy" OUTCAR
grep "F=" OSZICAR | tail -3
grep -i "dipol" OUTCAR | tail -5         # dipole moment along z, should be small for clean Pt(111)
python -c "from ase.io import read; a=read('CONTCAR'); import numpy as np; \
z=np.sort(a.positions[:,2]); print(np.round(np.diff(z[::4]),3))"   # interlayer spacings
```

- 弛豫后表面第一层与第二层的间距相对体相（3.97/√3 ≈ 2.292 Å）略有变化（金属表面通常是百分之一量级的弛豫，Pt(111) 的顶层间距略有扩张，具体以你的计算为准）。
- 固定的底部两层原子坐标应与 POSCAR 完全相同。
- 干净 Pt(111) 的偶极矩很小；如果打印出很大的偶极矩，检查 DIPOL 是否落在真空中、或 slab 是否跨越了 c 方向的周期边界。

## 注意事项

- 层数与真空要做收敛测试：4 层对 Pt(111) 的功函数和吸附能只是"够用的起点"，表面能请看 [09-cu100-surface-energy](../09-cu100-surface-energy/README.md) 中的层数收敛。
- slab 的原子不能跨越 z 方向的周期边界（所有 z 分数坐标都应在 0–1 之间且远离 0 和 1），否则偶极校正和功函数都会出错。
- `ISMEAR = -5` 不适合 slab 弛豫（受力不准，且 k 点少时报错），金属用 1 或 2，最终能量取 `energy(sigma->0)`。

[返回第9章](../../tutorials/09-surface-workfunction.md)
