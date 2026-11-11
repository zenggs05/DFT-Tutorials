# 示例 04-2：单层 MoS2 只优化面内晶格（LATTICE_CONSTRAINTS）

## 目的

- 演示二维材料的变胞优化：放开 a、b，固定真空方向 c；
- 学会使用 `LATTICE_CONSTRAINTS = .TRUE. .TRUE. .FALSE.` 并检查它是否真的生效。

对应教程：[第4章 结构优化与状态方程](../../tutorials/04-relax-eos.md)（§6）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `INCAR` | `ISIF=3` + `LATTICE_CONSTRAINTS = .TRUE. .TRUE. .FALSE.`，`ENCUT=450`，`EDIFFG=-0.01` |
| `POSCAR` | 1H-MoS2，120° 六角晶胞，初始 a = 3.20 Å，c = 20 Å（约 17 Å 真空），Mo 在中间层 |
| `KPOINTS` | Γ 中心 12×12×1（二维体系 z 方向取 1，六角晶格必须用 Γ 中心） |
| `job.slurm` | 单次运行的 Slurm 脚本，结束后打印晶格与应力 |
| `POTCAR` | **不提供**，按 `Mo_pv`、`S` 顺序拼接 |

## 运行步骤

```bash
cd VASP/examples/04-relax-mos2-2d
cat $VASP_PP_PATH/potpaw_PBE/{Mo_pv,S}/POTCAR > POTCAR
grep TITEL POTCAR                 # 顺序必须与 POSCAR 第 6 行（Mo S）一致
sbatch job.slurm
```

## 预期结果 / 如何检查

```bash
grep "reached required accuracy" OUTCAR
sed -n 3,5p CONTCAR               # 第三个晶格矢量必须仍为 0 0 20
grep "in kB" OUTCAR | tail -1     # XX、YY、XY 分量应接近 0；ZZ 不参与优化，可以不为 0
python -c "from ase.io import read; a=read('CONTCAR'); print('a = %.4f A, c = %.4f A' % (a.cell.lengths()[0], a.cell.lengths()[2]))"
```

- 面内晶格常数应从 3.20 Å 收缩；PBE 一般给出约 3.18–3.19 Å（近似值；实验约 3.16 Å）。
- c 轴保持 20.000 Å。**如果 c 轴变了**，说明你的 VASP 版本不支持该组合（见注意事项）。

## 注意事项

- `LATTICE_CONSTRAINTS` 与 `IBRION=1/2` 组合从 **VASP 6.4.3** 起才被读取；VASP Wiki 提示 6.5.0 之前存在相关 bug，与 `ISIF=4/5` 组合需要 **6.5.1+**。推荐 6.5.1 或更新版本。
- 旧版本会静默忽略该标签，此时请改用教程 §6.5 的 E–a 扫描方法（固定 c，缩放 a，`ISIF=2`）。
- 单层的面内晶格基本不受 vdW 影响；若要与双层/异质结比较能量，请统一加 vdW 修正（见 [第8章](../../tutorials/08-vdw-2d.md)）。
- 不要对二维 slab 直接 `ISIF=3` 而不加约束，真空层会被压缩。
