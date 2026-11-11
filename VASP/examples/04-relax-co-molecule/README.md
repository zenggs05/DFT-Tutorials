# 示例 04-3：CO 分子在真空盒子中的结构优化（vasp_gam）

## 目的

- 演示孤立分子的标准设置：大盒子、只用 Γ 点、`vasp_gam`、`ISMEAR=0`、`ISIF=2`；
- 得到 CO 的平衡键长和总能，后者可直接用于第10章的吸附能计算。

对应教程：[第4章 结构优化与状态方程](../../tutorials/04-relax-eos.md)（§7）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `INCAR` | `ENCUT=520`、`ISMEAR=0`、`SIGMA=0.01`、`ISPIN=1`、`IBRION=2`、`ISIF=2` |
| `POSCAR` | 15 Å 立方盒子，CO 沿 z 方向放在中心，初始键长 1.16 Å（笛卡尔坐标） |
| `KPOINTS` | 只含 Γ 点（`Gamma 1 1 1`），与 `vasp_gam` 配套 |
| `job.slurm` | 使用 `vasp_gam` 的 Slurm 脚本 |
| `POTCAR` | **不提供**，按 `C`、`O` 顺序拼接 |

## 运行步骤

```bash
cd VASP/examples/04-relax-co-molecule
cat $VASP_PP_PATH/potpaw_PBE/{C,O}/POTCAR > POTCAR
grep TITEL POTCAR
sbatch job.slurm                  # 或：mpirun -np 8 vasp_gam > vasp.out
```

## 预期结果 / 如何检查

```bash
grep "reached required accuracy" OUTCAR
grep "F=" OSZICAR | tail -1                       # 最后一步的总能（E0 列即 sigma->0 能量）
python -c "from ase.io import read; a=read('CONTCAR'); print('d(C-O) = %.4f A' % a.get_distance(0, 1, mic=True))"
```

- 键长应收缩到约 1.13–1.14 Å（PBE 近似值；实验值 1.128 Å）。
- 一般 5–10 个离子步内收敛。

## 注意事项

- `vasp_gam` 只能用于只有 Γ 点的 KPOINTS；体相计算请用 `vasp_std`。若集群没有编译 `vasp_gam`，用 `vasp_std` 也能得到相同结果，只是更慢。
- 不要用 `ISMEAR=-5`（四面体方法需要至少 4 个 k 点），也不要优化盒子（`ISIF` 必须为 2）。
- O2 等开壳层分子需要 `ISPIN=2` 并设置 `MAGMOM`（如 `MAGMOM = 2*1`）或 `NUPDOWN = 2`。
- 用于吸附能时，气相分子与表面计算的 `ENCUT`、`PREC`、泛函、vdW 修正必须一致；建议测试 12/15/18 Å 三种盒子，能量差 < 1–2 meV 即可。
