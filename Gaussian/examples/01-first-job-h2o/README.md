# 示例 01：第一个 Gaussian 计算——水分子单点与几何优化

配套教程：[第1章 环境准备与第一个 Gaussian 计算](../../tutorials/01-setup-first-job.md)

## 目的

- 用最小的输入文件（水分子，HF/6-31G(d)）跑通 Gaussian 16，熟悉 Link 0 / route / title / 电荷多重度 / 坐标 / 空行这几部分。
- 运行第一个几何优化（B3LYP/6-31G(d) Opt），观察每一步的能量和收敛判据。
- 学会判断作业是否 Normal termination，失败时定位出错的 Link。

## 文件说明

| 文件 | 内容 |
|---|---|
| `h2o_hf_sp.gjf` | `#p HF/6-31G(d) SP`：水分子单点能，最小示例 |
| `h2o_b3lyp_opt.gjf` | `#p B3LYP/6-31G(d) Opt`：水分子几何优化，初始结构故意给得较粗糙 |
| `run_g16_local.sh` | 本机（Linux/WSL）依次运行 gjf，自动去除 CRLF、检查末尾空行，并报告 OK/FAILED、最后的 SCF 能量和出错的 Link |
| `run_g16_slurm.sh` | 最小 Slurm 提交脚本（每个作业独立 scratch，结束时清理） |

## 运行步骤

```bash
# 0) make sure Gaussian is set up (see chapter 1, section 4)
export g16root=/path/to/parent/of/g16
source $g16root/g16/bsd/g16.profile
export GAUSS_SCRDIR=/tmp/$USER/g16scr && mkdir -p $GAUSS_SCRDIR

cd Gaussian/examples/01-first-job-h2o

# 1a) run directly
g16 < h2o_hf_sp.gjf > h2o_hf_sp.log
g16 < h2o_b3lyp_opt.gjf > h2o_b3lyp_opt.log

# 1b) or use the helper script (runs every *.gjf and checks termination)
bash run_g16_local.sh

# 1c) or submit to a Slurm cluster
sbatch run_g16_slurm.sh h2o_b3lyp_opt.gjf
```

两个作业在普通电脑上都只需几秒。

## 预期结果 / 如何检查

```bash
tail -n 1 h2o_hf_sp.log                     # Normal termination of Gaussian 16 ...
grep "SCF Done" h2o_hf_sp.log               # E(RHF) about -76.01 Hartree
grep -m1 "basis functions" h2o_hf_sp.log    # 19 basis functions for H2O/6-31G(d)
grep -A1 "Dipole moment" h2o_hf_sp.log

grep -c "SCF Done" h2o_b3lyp_opt.log        # number of optimization steps
grep "Stationary point found" h2o_b3lyp_opt.log
grep -A8 "Optimized Parameters" h2o_b3lyp_opt.log
```

判断标准：

- 两个 log 的最后一行都是 `Normal termination of Gaussian 16 at ...`。
- HF/6-31G(d) 单点能在 −76.01 Hartree 附近；水分子用 6-31G(d)（6D）共有 19 个基函数（O 15 个 + 每个 H 2 个）。
- 优化出现 `-- Stationary point found.`，O–H 约 0.97 Å、H–O–H 约 103–104°。实验值为 0.9572 Å 和 104.52°（NIST CCCBDB）；B3LYP/6-31G(d) 与实验有零点几个百分点的偏差，属于正常。以上计算值为近似值，请以你自己的输出为准。
- 用 GaussView 打开 `h2o_b3lyp_opt.log`，Results → Optimization 可以看到能量随步数下降的曲线。

也可以用仓库脚本批量汇总（参数用 `-h` 查看）：

```bash
python ../../scripts/g16_log_summary.py -h
```

## 注意事项

- 文件末尾必须有空行，否则 L101 报 `End of file in ZSymb.`。在 Windows 上编辑过的文件先 `dos2unix`。
- `%nprocshared` 和 `%mem` 按自己机器的情况调整。在集群上，`%nprocshared` 要等于 `--cpus-per-task`，`%mem` 约为 `--mem` 的 75–85%。
- 登录节点上只能跑这种几秒钟的测试；正式计算请用 `run_g16_slurm.sh` 或第17章的模板。
- `Opt` 只保证找到驻点，正式工作中应写 `Opt Freq` 并确认没有虚频（第6章）。
- 生成的 `.chk` 可以用 `formchk h2o_b3lyp_opt.chk` 转成 `.fchk`，在 GaussView 中查看轨道。
