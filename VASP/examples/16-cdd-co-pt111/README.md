# 示例 16-2：CO/Pt(111) 的差分电荷密度（CDD）

## 目的

- 按"同一晶胞、同一坐标、同一网格"的原则完成三次单点：CO/Pt(111)、Pt slab、CO；
- 计算 Δρ = ρ(CO/Pt) − ρ(Pt) − ρ(CO)，写出 VESTA 可读的 `CHGCAR_cdd`；
- 计算沿 z 的平面平均 Δλ(z) 与电荷位移曲线 ΔQ(z)。

对应教程：[第16章 电荷与成键分析](../../tutorials/16-charge-bonding.md)（§5）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `POSCAR` | 与示例 16-1 相同的 CO/Pt(111) p(2×2) 初始结构；正式计算请换成[第10章](../../tutorials/10-adsorption-catalysis.md)优化后的 CONTCAR |
| `INCAR` | 三次计算共用的单点设置（`LCHARG = .TRUE.`、`ISPIN = 1`） |
| `KPOINTS` | Γ 中心 6×6×1，三次计算共用 |
| `split_fragments.py` | 从 POSCAR 中按元素或编号切出片段，晶胞和坐标按原文本复制，不做任何舍入 |
| `run_cdd.sh` | `prepare` / `run` / `analyze` 三步；`run` 会把 AB 的 NGXF/NGYF/NGZF 写进片段的 INCAR |
| `cdd_planar.py` | 读 CHGCAR（只用第一块）、检查网格与晶胞、求差、输出 `CHGCAR_cdd`、`cdd_planar.dat`、`cdd_planar.png`（numpy，matplotlib 可选） |
| `job.slurm` | Slurm 脚本，依次执行三步 |
| `POTCAR` | **不提供**；`run_cdd.sh prepare` 在设置了 `VASP_PP_PATH` 时自动为三个目录拼接 |

## 运行步骤

```bash
cd VASP/examples/16-cdd-co-pt111
export VASP_PP_PATH=/path/to/your/potentials      # 目录下应有 potpaw_PBE/

# 一次性提交（推荐）
sbatch job.slurm

# 或者分步执行
bash run_cdd.sh prepare      # 生成 AB/ slab/ co/
NP=32 bash run_cdd.sh run    # 依次运行三个单点
bash run_cdd.sh analyze      # 求差与平面平均
```

也可以用 vaspkit 314（依次输入 `AB/CHGCAR slab/CHGCAR co/CHGCAR`）或仓库脚本 [`chgdiff.py`](../../scripts/chgdiff.py)（`python ../../scripts/chgdiff.py -h` 查看用法）得到同样的结果。

## 预期结果 / 如何检查

```bash
grep -H "dimension x,y,z NGXF" AB/OUTCAR slab/OUTCAR co/OUTCAR   # 三行数字必须相同
head -4 cdd_planar.dat
```

- `cdd_planar.py` 打印的三个电子数分别约为 170、160、10（标准 Pt/C/O POTCAR），`integral of d_rho` 应接近 0；
- `cdd_planar.dat` 第 4 列 ΔQ(z) 在 slab 下方的真空区为 0，越过 CO 后回到约 0；
- 在 VESTA 中打开 `CHGCAR_cdd`，等值面取约 ±0.002 e/Å³：C–Pt 之间和 CO 的 π 区域会出现正负交替的电荷重新分布（σ 给予与 π 反馈），具体形状以自己的计算为准。

## 注意事项

- **片段不能重新优化**，否则原子核附近会出现大块伪影；
- 三次计算的 `ENCUT`、`PREC`、k 点、`ISMEAR/SIGMA`、POTCAR 版本必须一致；
- 本例为非磁性体系，`ISPIN = 1`；磁性片段若用 `ISPIN = 2`，CHGCAR 第二块是磁化密度，`cdd_planar.py` 只使用第一块；
- 大网格的 CHGCAR 可达数百 MB，不要提交到仓库。
