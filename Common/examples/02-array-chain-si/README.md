# 示例 02-array-chain-si：用作业数组做 ENCUT 扫描、用依赖链串联弛豫→静态→能带（VASP，Si）

## 目的

以两原子 Si 原胞为例，演示批量计算的三种组织方式：

1. **作业数组**（Slurm `--array`、PBS `-J`、LSF `-J "name[1-7]"`）：同一脚本按数组下标取不同的 ENCUT，7 个计算并行排队，最多同时运行 4 个；
2. **作业依赖链**（`--dependency=afterok`）：弛豫 → 静态自洽（写 CHGCAR）→ 非自洽能带，前一步失败时后续作业自动取消；
3. **dpdispatcher**：用 Python 描述任务和资源，由它生成 Slurm 脚本、提交、轮询并回收文件。

对应教程：[第2章 HPC 作业调度：Slurm/PBS/LSF 模板与并行设置](../../tutorials/02-job-scheduler.md)，第 7–9 节。

## 文件说明

| 文件 | 说明 |
|---|---|
| [POSCAR](POSCAR) | Si 金刚石结构原胞，a = 5.431 Å（实验值，弛豫步会修正） |
| [KPOINTS](KPOINTS) | Γ 中心 8×8×8 网格，用于 ENCUT 扫描、弛豫和静态计算 |
| [KPOINTS.band](KPOINTS.band) | fcc 布里渊区 L–Γ–X–U\|K–Γ 线模式路径，每段 40 个点 |
| [INCAR.encut.template](INCAR.encut.template) | ENCUT 扫描模板，占位符 `@ENCUT@` 由脚本替换 |
| [INCAR.relax](INCAR.relax) | 第 1 步：`IBRION = 2`、`ISIF = 3` 同时弛豫原子与晶胞 |
| [INCAR.static](INCAR.static) | 第 2 步：`ISMEAR = -5` 静态自洽，`LCHARG = .TRUE.` 写出 CHGCAR |
| [INCAR.band](INCAR.band) | 第 3 步：`ICHARG = 11` 读入 CHGCAR 做非自洽能带 |
| [slurm_array_encut.sh](slurm_array_encut.sh) | Slurm 作业数组：`--array=0-6%4`，任务 i 在 `encut_<E>/` 中计算 |
| [pbs_array_encut.sh](pbs_array_encut.sh) | 同一扫描的 PBS Pro 写法（`#PBS -J 0-6`、`$PBS_ARRAY_INDEX`），注释中给出 Torque 写法 |
| [lsf_array_encut.sh](lsf_array_encut.sh) | 同一扫描的 LSF 写法（`#BSUB -J "si-encut[1-7]%4"`、`$LSB_JOBINDEX`） |
| [collect_encut.sh](collect_encut.sh) | 汇总 `encut_*/OUTCAR`：每原子能量及其与最大 ENCUT 的差（meV/atom） |
| [job_step.sh](job_step.sh) | 依赖链中的单步作业，由环境变量 `STEP=relax/static/band` 选择；未正常结束时以非零码退出 |
| [chain_relax_static_band.sh](chain_relax_static_band.sh) | 在登录节点运行，用 `sbatch --parsable` 和 `--dependency=afterok` 依次提交三步 |
| [dpdispatcher_encut.py](dpdispatcher_encut.py) | 用 dpdispatcher 提交同样的 ENCUT 扫描；`--dry-run` 只准备目录不提交 |

## 运行步骤

先准备 POTCAR（受许可证保护，本仓库不提供），例如用自己的 PBE 赝势库：

```bash
cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR      # 或: vaspkit -task 103
grep TITEL POTCAR                                     # 确认是 PAW_PBE Si
```

再按需修改各脚本头部的 `--partition`、核数和 `module load`（默认每个计算 8 核、30 分钟以内）。

**A. 作业数组做 ENCUT 扫描**

```bash
sbatch slurm_array_encut.sh            # PBS: qsub pbs_array_encut.sh   LSF: bsub < lsf_array_encut.sh
squeue -u $USER -r                     # -r: 每个数组任务单独一行
sacct -j <数组作业号> --format=JobID,State,Elapsed,ExitCode
bash collect_encut.sh                  # 全部完成后汇总
```

只重算其中几个值：`sbatch --array=2,5 slurm_array_encut.sh`。

**B. 依赖链：弛豫 → 静态 → 能带**

```bash
bash chain_relax_static_band.sh        # 在登录节点运行，只调用 sbatch
squeue -u $USER -o '%.10i %.12j %.2t %.10M %R'
```

**C. dpdispatcher（可选）**

```bash
pip install dpdispatcher               # 或使用第1章的 dft conda 环境
python dpdispatcher_encut.py --dry-run # 只准备 dp_encut_*/ 和 dp_work/，不提交
python dpdispatcher_encut.py           # 提交并等待（建议在 tmux 中运行）
```

## 预期结果 / 如何检查

- 数组作业：`squeue -r` 中作业号形如 `123456_0`、`123456_1`……，同时处于 `R` 状态的不超过 4 个；完成后出现 `encut_250/` 到 `encut_550/` 七个目录。
- `collect_encut.sh` 输出三列：ENCUT、每原子能量、与 550 eV 结果的差。下面是格式示意（示意输出，数值仅供参考）：

  ```text
     ENCUT      E/atom (eV)  dE (meV/atom)
       250        -5.4xxxxx          xx.xx
       ...
       550        -5.4xxxxx           0.00
  ```

  判断标准：取 |dE| < 1 meV/atom 的最小 ENCUT，具体数值以你自己的扫描结果为准；Si 的 PBE POTCAR 中 ENMAX 约为 245 eV，弛豫晶胞时还要在收敛值之上留出余量（INCAR.relax 用 400 eV，约 1.6 × ENMAX）。

- 依赖链：提交后 `squeue` 中 si-static、si-band 的 `REASON` 列显示 `(Dependency)`；弛豫完成后依次运行。`1-relax/OUTCAR` 中应有 `reached required accuracy`，`2-static/` 中有非空 CHGCAR，`3-band/` 中有 EIGENVAL 和 vasprun.xml，可用 [VASP 第5章](../../../VASP/tutorials/05-band-dos.md) 的方法作图。
- 若某一步失败（例如弛豫在 `NSW` 步内未收敛，`job_step.sh` 以退出码 3 结束），后续作业因 `--kill-on-invalid-dep=yes` 被自动取消，`sacct` 中显示 `CANCELLED`。

## 注意事项

- 这是演示调度功能的小算例，ENCUT 和 k 点并未针对生产计算优化；ENCUT 扫描请保持其余参数完全一致。
- `ISIF = 3` 的晶胞弛豫受 Pulay 应力影响，生产计算建议在收敛后用新结构再弛豫一次（见 [VASP 第4章](../../../VASP/tutorials/04-relax-eos.md)）。
- 能带步的 `ISMEAR = 0`：线模式 k 点不能用四面体方法。
- `dpdispatcher_encut.py` 中 `mpirun -np 8` 要与 `cpu_per_node` 一致；在远程集群上使用时把 `context_type` 改为 `SSHContext` 并填写 `remote_profile`（主机名、用户名、密钥），具体字段见 dpdispatcher 文档。
- 所有脚本都应以 LF 换行保存；Windows 上编辑后先 `dos2unix *.sh INCAR*`。
