# 示例 03：ENCUT、k 点与展宽（SIGMA）收敛测试（Si 与金属 Al）

配套教程：[第3章 收敛测试与计算精度控制](../../tutorials/03-convergence.md)

## 目的

- 对金刚石 Si 做 ENCUT 扫描（能量、压强/Pulay 应力、原子受力三种判据）；
- 对 Si 做 k 点扫描（N×N×N Γ 网格或 KSPACING 两种方式）；
- 对金属 fcc Al 做 ISMEAR/SIGMA 与 k 网格的联合测试，检查熵项 `|T*S|` < 1 meV/atom；
- 提供可直接改造的批量脚本：本地运行、Slurm 提交、汇总与作图。

## 文件说明

| 文件 | 说明 |
|---|---|
| `common.sh` | 公共函数：`natoms`、`outcar_row`（从 OUTCAR 提取 E0/TOTEN/EENTRO/压强/最大力/NKPTS/耗时）、`run_or_submit`、`check_potcar` |
| `plot_conv.py` | 读取 `.dat` 结果表，计算相对参考点的偏差，判断收敛并画图（需要 numpy、matplotlib） |
| `encut_test/INCAR.template` | ENCUT 扫描模板（`__ENCUT__` 占位；EDIFF=1E-7、ISIF=2 输出应力） |
| `encut_test/POSCAR` | 理想金刚石 Si，a = 5.43 Å |
| `encut_test/POSCAR_displaced` | 第二个 Si 原子沿 a1 方向位移约 0.08 Å，用于力的收敛 |
| `encut_test/KPOINTS` | Γ 中心 8×8×8（扫描 ENCUT 时固定） |
| `encut_test/run_encut.sh` | 建目录 → 运行/提交 → 汇总为 `encut_<POSCAR>.dat` |
| `kpoint_test/INCAR` | k 点扫描用 INCAR（ENCUT = 320 eV） |
| `kpoint_test/POSCAR` | 同上 Si |
| `kpoint_test/run_kpts.sh` | `SCAN=mesh`（写 KPOINTS）或 `SCAN=kspacing`（INCAR 追加 KSPACING） |
| `sigma_check_Al/INCAR.template` | Al 展宽测试模板（`__ISMEAR__`、`__SIGMA__` 占位） |
| `sigma_check_Al/POSCAR` | fcc Al 原胞，a = 4.05 Å（实验值约 4.05 Å） |
| `sigma_check_Al/run_sigma.sh` | ISMEAR × SIGMA × k 网格扫描，外加 ISMEAR=−5 + 24³ 网格参考 |

> [!IMPORTANT]
> 本目录不含 POTCAR。每个测试目录运行前都要先放入对应的 POTCAR（Si 或 Al）。

## 运行步骤

```bash
export VASP_PP_PATH=$HOME/vasp/potentials       # 含 potpaw_PBE/Si、potpaw_PBE/Al
export OMP_NUM_THREADS=1

# 1) ENCUT
cd encut_test
cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR
NP=4 bash run_encut.sh                          # 理想结构
NP=4 bash run_encut.sh POSCAR_displaced         # 力的收敛
python ../plot_conv.py encut_POSCAR.dat
python ../plot_conv.py encut_POSCAR_displaced.dat --thr-f 5
cd ..

# 2) k 点
cd kpoint_test
cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR
bash run_kpts.sh                                # N = 2 4 6 ... 16
SCAN=kspacing bash run_kpts.sh                  # KSPACING = 0.5 ... 0.15
python ../plot_conv.py kpts_mesh.dat
python ../plot_conv.py kpts_kspacing.dat --ref min
cd ..

# 3) Al 展宽
cd sigma_check_Al
cat $VASP_PP_PATH/potpaw_PBE/Al/POTCAR > POTCAR
bash run_sigma.sh
```

在 Slurm 集群上：

```bash
MODE=slurm bash run_encut.sh                    # 每个 ENCUT 一个作业（模板默认 ../01-setup-compile/job_cpu.slurm）
JOB_TEMPLATE=$HOME/my_job.slurm MODE=slurm bash run_kpts.sh
# 全部作业完成后：
MODE=collect bash run_encut.sh
MODE=collect bash run_kpts.sh
```

常用环境变量：`ENCUTS`、`MESHES`、`SPACINGS`、`ISMEARS`、`SIGMAS`、`REF_N`、`NP`、`VASP_CMD`、`MODE`、`JOB_TEMPLATE`（含义见各脚本开头的注释）。

## 预期结果 / 如何检查

| 测试 | 产物 | 判断方法 |
|---|---|---|
| ENCUT | `encut_POSCAR.dat`、`encut_POSCAR.png` | `plot_conv.py` 给出 ΔE < 1 meV/atom、ΔP < 1 kB 的最小 ENCUT；常见情况是压强比能量收敛得慢（Pulay 应力） |
| 力 | `encut_POSCAR_displaced.dat` | 最大力偏差 < 5 meV/Å |
| k 点 | `kpts_mesh.dat`、`kpts_kspacing.dat` | ΔE < 1 meV/atom；NKPTS 列与教程第 5.1 节的表对照（如 8³ → 29） |
| 展宽 | `sigma_Al.dat` | 脚本末尾标 `OK` 的组合满足 \|`T*S`\| < 1 meV/atom 且 \|ΔE0\| < 1 meV/atom；同一 SIGMA 下 MP（ISMEAR=1）的 \|`T*S`\| 通常小于 Gaussian |

每个子目录都应在 OUTCAR 末尾出现 `General timing and accounting`；未正常结束的目录会在汇总时打印 WARNING，结果列显示为 NaN。

## 注意事项

- 扫描时其余参数全部固定，且每个目录 `ISTART = 0`（模板已设置），避免读入旧 WAVECAR；
- 本示例的 Si ENCUT = 320 eV 只是 k 点测试的起点；请用你自己的 ENCUT 测试结果替换 `kpoint_test/INCAR` 中的值；
- 生产计算的 EDIFF 可以放宽到 1E-6，但测试阶段建议 1E-7；
- 改造成自己的体系时，替换 POSCAR/POTCAR，按需要修改 INCAR（自旋、DFT+U、泛函），并相应调整扫描范围（例如含 O 的体系 ENCUT 从 400 eV 起扫）；
- 输入文件均为纯 ASCII、LF 换行；生成的 `*.dat`、`*.png` 与各计算子目录不要提交到仓库。
