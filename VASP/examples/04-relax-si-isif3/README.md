# 示例 04-1：Si 体相变胞优化（ISIF=3）与 Pulay 应力消除

## 目的

- 演示 `IBRION = 2` + `ISIF = 3` 的体相完全优化；
- 演示"优化 → CONTCAR 拷成 POSCAR → 再优化"的循环，直到体积不再变化，从而消除 Pulay 应力；
- 学会用 `grep` 检查 "reached required accuracy"、`external pressure` 和体积。

对应教程：[第4章 结构优化与状态方程](../../tutorials/04-relax-eos.md)（§4）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `INCAR` | `PREC=Accurate`、`ENCUT=400`（≥1.3×ENMAX）、`IBRION=2`、`ISIF=3`、`EDIFFG=-0.005` |
| `POSCAR` | 金刚石 Si 的 fcc 原胞，初始晶格常数故意取 5.55 Å（比实验值大约 2%） |
| `KPOINTS` | Γ 中心 8×8×8 网格（fcc 晶格推荐 Γ 中心） |
| `job.slurm` | Slurm 作业脚本：最多 4 轮优化，每轮结束后 `cp CONTCAR POSCAR`，某轮离子步数 ≤ 3 且收敛即停止 |
| `POTCAR` | **不提供**，需自行生成（见下） |

## 运行步骤

```bash
cd VASP/examples/04-relax-si-isif3
# 1) 生成 POTCAR（按你本地赝势库的目录名调整 potpaw_PBE）
cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR
grep -E "TITEL|ENMAX" POTCAR

# 2) 修改 job.slurm 中的分区、模块名、核数后提交
sbatch job.slurm

# 3) 查看每一轮的结果
cat relax_rounds.log
```

不想用循环时，也可以手动运行一次：`mpirun -np 4 vasp_std > vasp.out`，然后 `cp CONTCAR POSCAR` 再运行一次。

## 预期结果 / 如何检查

```bash
grep "reached required accuracy" round*/OUTCAR     # 每一轮都应出现
grep "external pressure" OUTCAR | tail -1           # 最后一轮 |P| 应在约 1 kB 以内
head -5 CONTCAR                                     # 优化后的晶格矢量
python -c "from ase.io import read; a=read('CONTCAR'); print('a = %.4f A' % (a.cell.lengths()[0]*2**0.5))"
```

- `relax_rounds.log` 中，第 1 轮体积明显减小，后续轮次体积变化应 < 0.2%，离子步数迅速减少到 1–3 步。
- PBE 下 Si 的晶格常数通常在约 5.46–5.47 Å（近似值；实验值 5.431 Å），具体数值取决于 ENCUT、k 点和 POTCAR 版本。

## 注意事项

- 仓库不包含 POTCAR（VASP 许可证限制），也不要把 POTCAR 提交到公开仓库。
- `ISIF=3` 必须使用足够高的 `ENCUT`；把 `ENCUT` 降到约 250 eV 再跑一遍，可以直观看到 Pulay 应力导致的体积偏小（教程练习 1）。
- 金刚石 Si 的原子位置被对称性固定，受力恒为 0，收敛实际由应力决定；有内坐标的体系同样适用这套流程。
- 若作业中途被杀，`POSCAR.orig` 保存了最初的结构，`round*/` 中保存了每轮的 OUTCAR/CONTCAR。
