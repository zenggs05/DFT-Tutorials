# 示例 05：几何优化（Gaussian 16）

配套教程：[第5章 几何优化](../../tutorials/05-opt.md)

## 目的

用小分子和一个表面簇模型练习 Gaussian 16 几何优化的常见写法：

- `Opt Freq` 联合作业，读懂收敛表；
- `Opt=(CalcFC,MaxCycles=N)` 和 `Opt=ReadFC`（读入低水平频率的 Hessian）；
- 坐标前加冻结码 `-1` 冻结原子；
- `Opt=ModRedundant` 冻结键长（`B i j F`）和二面角（`D i j k l F`）；
- 用 `Opt=Restart` 或 `Geom=Check Guess=Read` 续算。

## 文件说明

| 文件 | 体系 | 演示内容 | 依赖 |
|---|---|---|---|
| `acetaldehyde_opt_freq.gjf` | 乙醛 | 标准 `Opt Freq` | — |
| `opt_maxcycles_calcfc.gjf` | 甘氨酸 | `Opt=(CalcFC,MaxCycles=200) Freq` | — |
| `opt_readfc_from_lowlevel.gjf` | 甘氨酸 | 第 1 步 def2-SVP `Opt Freq`，第 2 步 def2-TZVP `Opt=ReadFC` | — |
| `opt_frozen_atoms_surface_cluster.gjf` | Cu10(111) 两层簇 + 顶位 CO | Cu 冻结码 `-1`，只优化 CO；`NoSymm` | — |
| `opt_modredundant_freeze_bond.gjf` | 水二聚体 | `B 1 4 F`：O···O 距离固定在 3.10 Å | — |
| `opt_modredundant_freeze_dihedral.gjf` | 正丁烷（gauche） | `D 1 2 3 4 F`：C–C–C–C 二面角固定在 65° | — |
| `opt_restart_from_chk.gjf` | 甘氨酸 | `Opt=(Restart,MaxCycles=300)`，只有 Link0 和 route | 需要 `opt_maxcycles_calcfc.chk` |
| `opt_continue_geom_check.gjf` | 甘氨酸 | `%oldchk` + `Geom=Check Guess=Read` 从最后一步继续 | 需要 `opt_maxcycles_calcfc.chk` |
| `run_g16.sh` | — | 最简 Slurm 作业脚本（8 核、16 GB） | — |

起始结构说明：乙醛、甘氨酸、丁烷由 RDKit（ETKDG + MMFF）生成，丁烷的二面角随后设为 65°；水二聚体和 Cu10 簇按典型几何参数手工搭建（Cu–Cu 最近邻 2.556 Å，对应 fcc Cu 晶格常数 3.615 Å；Cu–C 1.85 Å，C–O 1.15 Å）。

## 运行步骤

```bash
cd Gaussian/examples/05-opt

# independent jobs
for f in acetaldehyde_opt_freq opt_maxcycles_calcfc opt_readfc_from_lowlevel \
         opt_modredundant_freeze_bond opt_modredundant_freeze_dihedral \
         opt_frozen_atoms_surface_cluster; do
    g16 < $f.gjf > $f.log
done

# restart demos: they read opt_maxcycles_calcfc.chk
cp opt_maxcycles_calcfc.chk opt_maxcycles_calcfc.chk.bak     # Opt=Restart writes into this chk
g16 < opt_continue_geom_check.gjf > opt_continue_geom_check.log
g16 < opt_restart_from_chk.gjf   > opt_restart_from_chk.log
```

Slurm：`sbatch run_g16.sh acetaldehyde_opt_freq.gjf`（先按本集群修改脚本中加载 Gaussian 的两行）。

**怎样真正练习续算：** `Opt=Restart` 是为中断或失败的优化设计的，对已经正常收敛的 `opt_maxcycles_calcfc` 做续算没有实际意义。想看到续算的效果，先把 `opt_maxcycles_calcfc.gjf` 中的 `MaxCycles=200` 临时改成 `MaxCycles=4` 运行，让它以 L9999 失败，再运行两个续算输入。

## 预期结果 / 如何检查

```bash
# Opt Freq jobs print "Normal termination" twice
grep -c "Normal termination" *.log

# convergence: last table of each job, and the final status line
for f in *.log; do
    echo "== $f"
    grep -E "^ (Maximum|RMS) +(Force|Displacement)" "$f" | tail -n 4
    grep -E "Stationary point found|negligible forces|Number of steps exceeded" "$f" | tail -n 1
done

# imaginary frequencies (should be none for the minima)
grep "imaginary frequencies" *.log

# frozen coordinates keep their input values
grep -A 30 "Optimized Parameters" opt_modredundant_freeze_bond.log | grep "R(1,4)"
grep -A 60 "Optimized Parameters" opt_modredundant_freeze_dihedral.log | grep "D(1,2,3,4)"

# ReadFC really read the Hessian in step 2
grep "D2E/DX2" opt_readfc_from_lowlevel.log | head
```

判断标准：

| 作业 | 应该看到 |
|---|---|
| 乙醛、甘氨酸 | `Stationary point found`；Freq 步最后一张收敛表全部为 `YES`；没有 `imaginary frequencies` 行 |
| `opt_readfc_from_lowlevel` | 第 2 步 `Initial Parameters` 表的 `Derivative Info` 列不再是 `estimate D2E/DX2`；第 2 步的步数明显少于第 1 步 |
| 水二聚体 | `R(1,4)` 最终仍为 3.10 Å，`-DE/DX` 不为零（约束上存在残余力，属正常）；其余坐标已弛豫 |
| 丁烷 | `D(1,2,3,4)` 最终仍为 65°；C–C 键长、C–H 键长已弛豫 |
| Cu10 + CO | 所有 Cu 坐标不变，只有 C、O 移动；CO 在顶位上，Cu–C 与 C–O 距离被优化 |
| 续算 | `opt_continue_geom_check` 从上一作业的最后结构开始，第 1 步 SCF 能量与上一作业最后一步接近 |

如果甘氨酸优化后出现虚频，说明起始构象不是极小点附近的构象，按教程第 10 节处理（沿虚频位移后重新优化，或先做构象搜索）。

## 注意事项

- 受限优化（冻结键长、二面角或原子）得到的不是驻点，不要对其做热化学分析。
- `Opt=Restart` 会继续写入原来的 chk，请先备份；如果原优化是因为结构跑偏而失败，应改用 `Geom=Check Guess=Read`，或从 log 中取能量最低的结构重新开始（[log2xyz.py](../../scripts/log2xyz.py)，`python ../../scripts/log2xyz.py -h`）。
- ModRedundant 中的原子编号就是输入中原子的顺序。在 GaussView 中重新排序原子后，要同步修改约束编号。
- 监控优化过程可以用 GaussView（`Results → Optimization`）或 [opt_monitor.py](../../scripts/opt_monitor.py)（`python ../../scripts/opt_monitor.py -h`）。
- Cu10 簇模型只用于演示冻结原子的写法。研究表面吸附时，簇的大小、自旋态和泛函都需要单独验证，周期性模型见本仓库的 VASP/CP2K 部分。
