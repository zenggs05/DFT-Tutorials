# 示例 01：VASP 编译配置、作业脚本与安装验证（Si）

配套教程：[第1章 环境准备、VASP 编译安装（CPU/GPU）与作业提交](../../tutorials/01-setup-compile.md)

## 目的

1. 给出编译 VASP 6.6.x 时需要对官方模板做的**改动清单**（不包含官方 makefile.include 原文）；
2. 给出 VTST 补丁的操作步骤与只读检查脚本；
3. 提供 Slurm（CPU/GPU）与 PBS 作业脚本模板；
4. 用 2 原子金刚石 Si 快速验证安装：并行一致性、vasp_gam 与 vasp_ncl 的一致性。

## 文件说明

| 文件 | 说明 |
|---|---|
| `makefile.include.oneapi.changes` | 基于 `arch/makefile.include.oneapi`（ifx + MKL + Intel MPI）的改动清单与 HDF5/可选库追加片段 |
| `makefile.include.nvhpc_acc.changes` | 基于 NVHPC OpenACC GPU 模板的改动清单（GPU 架构、CUDA 版本、QD、FFTW/MKL、HDF5） |
| `vtst_patch_steps.md` | VTST 补丁逐步操作清单（复制文件、改 `.objects`、`main.F`、`src/makefile`） |
| `check_vtst_patch.sh` | 只读检查脚本：确认 VTST 的 4 处修改都已完成 |
| `job_cpu.slurm` | Slurm 单节点纯 MPI 作业模板 |
| `job_gpu.slurm` | Slurm NVIDIA GPU 作业模板（1 rank/GPU + OpenMP） |
| `job_pbs.sh` | PBS/Torque 作业模板（附 PBS Pro 写法） |
| `test_Si/INCAR` | Si 单点能（PREC=Accurate、ENCUT=320、ISMEAR=0） |
| `test_Si/INCAR_ncl` | 同上，加 `LNONCOLLINEAR = .TRUE.`、`MAGMOM = 6*0`，用于 vasp_ncl 测试 |
| `test_Si/POSCAR` | 金刚石 Si 原胞，a = 5.43 Å（接近实验值 5.431 Å） |
| `test_Si/KPOINTS` | Γ 中心 6×6×6 |
| `test_Si/KPOINTS_gamma_only` | 只含 Γ 点，用于 vasp_gam 与 vasp_std 对比 |

> [!IMPORTANT]
> 本目录**不含 POTCAR**（受 VASP 许可证保护）。请用自己的赝势库生成，见下文。

## 运行步骤

### 1. 编译（概要）

```bash
cd ~/software/vasp.6.6.1
cp arch/makefile.include.oneapi makefile.include
# 按 makefile.include.oneapi.changes 修改并追加 HDF5 等内容
make veryclean && make DEPS=1 -j 16 all
```

如需 VTST：

```bash
bash check_vtst_patch.sh ~/software/vasp.6.6.1     # 全部 PASS 后再编译
```

### 2. 准备 Si 算例

```bash
cd test_Si
export VASP_PP_PATH=$HOME/vasp/potentials          # 目录下应有 potpaw_PBE/Si/POTCAR
cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR
grep TITEL POTCAR                                   # 应为 PAW_PBE Si ...
```

### 3. 并行一致性测试（vasp_std）

```bash
export OMP_NUM_THREADS=1
ulimit -s unlimited
for np in 1 4; do
    mpirun -np $np vasp_std > out.np$np 2>&1
    cp OUTCAR OUTCAR.np$np
    grep "energy(sigma->0)" OUTCAR | tail -1
done
```

### 4. vasp_gam 与 vasp_std 一致性测试（只用 Γ 点）

```bash
cp KPOINTS KPOINTS.666 && cp KPOINTS_gamma_only KPOINTS
mpirun -np 4 vasp_std > out.std_G 2>&1 && grep "energy(sigma->0)" OUTCAR | tail -1
mpirun -np 4 vasp_gam > out.gam_G 2>&1 && grep "energy(sigma->0)" OUTCAR | tail -1
cp KPOINTS.666 KPOINTS
```

### 5. vasp_ncl 一致性测试

```bash
cp INCAR INCAR.std && cp INCAR_ncl INCAR
mpirun -np 4 vasp_ncl > out.ncl 2>&1 && grep "energy(sigma->0)" OUTCAR | tail -1
cp INCAR.std INCAR
```

### 6. 在集群上提交

```bash
cp ../job_cpu.slurm . && sbatch job_cpu.slurm        # Slurm
# 或 cp ../job_pbs.sh . && qsub job_pbs.sh           # PBS
```

## 预期结果 / 如何检查

| 检查项 | 命令 | 判断标准 |
|---|---|---|
| 正常结束 | `grep "General timing and accounting" OUTCAR` | 有输出 |
| 电子步收敛 | `cat OSZICAR` | 最后一行有 `F=`、`E0=`；电子步数远小于 NELM=100 |
| 1 核与 4 核一致 | 比较 `OUTCAR.np1` 与 `OUTCAR.np4` 的 `energy(sigma->0)` | 差值约 $10^{-5}$ eV 量级或更小 |
| vasp_gam 与 vasp_std（Γ 点） | 第 4 步两次输出 | 应在 $10^{-4}$ eV 量级内一致 |
| vasp_ncl 与 vasp_std | 第 5 步与第 3 步 | 非磁性体系应在约 1 meV 内一致（经验判断） |
| 能量量级 | 第 3 步 | 2 原子 Si 约 −10.8 eV（经验量级，仅供参考） |

> [!NOTE]
> 只有 Γ 点的 Si 原胞在物理上严重欠收敛，第 4 步只用于检验两个程序是否给出相同结果，不要把那里的能量当作 Si 的总能。

## 注意事项

- 作业脚本中的分区名、模块名、可执行文件路径都需要按你的集群修改；
- GPU 脚本要求 INCAR 中设置 `KPAR` = GPU 数、`NSIM` = 16–32，每块 GPU 只用 1 个 MPI rank；
- `*.changes` 文件只是改动说明，不能直接当作 makefile.include 使用；最终的 makefile.include 含官方模板内容，**请保存在私有位置**，不要提交到公开仓库；
- 所有输入文件均为 ASCII、LF 换行。
