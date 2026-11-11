# 示例 16-3：金刚石 Si 的电子局域函数（ELF）

## 目的

- 用 `LELF = .TRUE.` 写出 ELFCAR，掌握 ELF 计算的并行限制（`KPAR = 1`、`NPAR = 1`）；
- 在 VESTA 中画出 Si–Si 共价键的 ELF 等值面和 (110) 切片。

对应教程：[第16章 电荷与成键分析](../../tutorials/16-charge-bonding.md)（§6）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `INCAR` | `PREC = Accurate`、`LELF = .TRUE.`、`KPAR = 1`、`NPAR = 1` |
| `POSCAR` | 金刚石 Si 原胞，a = 5.431 Å（实验值） |
| `KPOINTS` | Γ 中心 8×8×8 |
| `job.slurm` | Slurm 作业脚本 |
| `POTCAR` | **不提供**：`cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR` |

## 运行步骤

```bash
cd VASP/examples/16-elf-si
cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR
sbatch job.slurm            # 或：mpirun -np 8 vasp_std > vasp.out
```

## 预期结果 / 如何检查

```bash
grep "General timing" OUTCAR      # 计算正常结束
ls -lh ELFCAR                     # 文件存在
head -12 ELFCAR                   # 结构块 + 网格维数
```

- 在 VESTA 中打开 ELFCAR，等值面设为约 0.8：每个 Si–Si 键中点附近出现一个局域区域（共价键）；
- Utilities → 2D Data Display，切 (1 1 0) 面：键中点 ELF 接近 1，原子之间的空旷区很低；
- 想做原胞之外的大图，在 VESTA 中把 Boundary 设为 0–2 即可，不需要重新计算。

## 注意事项

- 用 `KPAR > 1` 会报 `ELF: KPAR>1 not implemented`；用 `vasp_ncl` 会得到 `WARNING: ELF not implemented for non collinear case`；
- ELF 在平面波网格上计算，原子核附近的数值受 PAW 赝化影响，只解读价层区域；
- ELFCAR 的数值不乘体积；用 ASE 的 `VaspChargeDensity` 读取时它会被除以体积，需要乘回来。
