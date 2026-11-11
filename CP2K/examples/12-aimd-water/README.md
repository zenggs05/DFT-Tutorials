# 示例 12：32 个水分子的 AIMD（NVT → 续算 → NVE / NPT → REFTRAJ）

配套教程：第12章 从头算分子动力学（AIMD）

## 目的

以 32 个水分子（96 原子、立方盒子 L = 9.8654 Å，对应 0.997 g/cm³）为例，走完一条完整的 Born–Oppenheimer AIMD 工作流：

1. `water32_nvt_csvr.inp`：NVT 平衡，CSVR 恒温器 `TIMECON 50 fs`，2 ps；
2. `md_restart.inp`：用 `&EXT_RESTART` 和 `RESTART.wfn` 续算 NVT 生产段，`TIMECON 200 fs`，5 ps；
3. `water32_nve.inp`：从平衡态出发跑 NVE，检查能量守恒，并每 1 fs 输出速度和体系偶极矩（用于 VDOS/IR）；
4. `water32_npt_i.inp`：从平衡态出发跑各向同性 NPT（1 atm），检验密度；
5. `reftraj_eval.inp`：用 REFTRAJ 对生产轨迹中的 10 帧换用 TZV2P 基组重算能量和力（2026.2 的 `EVAL ENERGY_FORCES` 写法）。

电子结构：PBE-D3(BJ)，`DZVP-MOLOPT-PBE-GTH`（BASIS_MOLOPT_UZH）+ `GTH-PBE`（POTENTIAL_UZH），`CUTOFF 400 / REL_CUTOFF 60`（NPT 用 600 Ry），OT + `OUTER_SCF`，`EPS_SCF 1E-6`，ASPC 外推（阶数 3），步长 0.5 fs。

## 文件说明

| 文件 | 说明 |
|---|---|
| `water32.xyz` | 初始结构：O 原子按 fcc 排布、分子随机取向（O–O 最近 3.49 Å，H…H ≥ 1.8 Å），**不是**平衡液体结构 |
| `gen_water_box.py` | 生成 `water32.xyz` 的脚本（只用标准库），`python gen_water_box.py --seed 2026` 可复现 |
| `water32_nvt_csvr.inp` | 第 1 步：NVT 平衡 |
| `md_restart.inp` | 第 2 步：NVT 生产段（续算） |
| `water32_nve.inp` | 第 3a 步：NVE + 偶极矩输出 |
| `water32_npt_i.inp` | 第 3b 步：NPT_I（`STRESS_TENSOR ANALYTICAL`、`CELL_REF`） |
| `reftraj_eval.inp` | 第 4 步：REFTRAJ 重算 |
| `job_slurm.sh` | Slurm 作业脚本（参数为输入文件名），用 `--dependency` 串联各步 |

## 运行步骤

```bash
cd CP2K/examples/12-aimd-water
export CP2K_DATA_DIR=/path/to/cp2k/data
cp2k.psmp -c water32_nvt_csvr.inp            # syntax check

# Slurm: chain the steps
jid=$(sbatch --parsable job_slurm.sh water32_nvt_csvr.inp)
pid=$(sbatch --parsable --dependency=afterok:$jid job_slurm.sh md_restart.inp)
sbatch --dependency=afterok:$jid job_slurm.sh water32_nve.inp
sbatch --dependency=afterok:$jid job_slurm.sh water32_npt_i.inp
sbatch --dependency=afterok:$pid job_slurm.sh reftraj_eval.inp

# local test run (shorten STEPS first, e.g. to 40)
mpirun -np 8 cp2k.psmp -i water32_nvt_csvr.inp -o water32_nvt_csvr.out
```

轨迹分析使用本仓库脚本（用法见 `-h`）：

```bash
python ../../scripts/md_analyze.py -h   # .ener drift / temperature, RDF, MSD, VACF, IR from dipoles
```

## 预期结果 / 如何检查

主要输出文件（以第 1 步为例，项目名 `water32_nvt`）：

| 文件 | 内容 |
|---|---|
| `water32_nvt-1.ener` | 每步：`Step Nr.`、`Time[fs]`、`Kin.[a.u.]`、`Temp[K]`、`Pot.[a.u.]`、`Cons Qty[a.u.]`、`UsedTime[s]` |
| `water32_nvt-pos-1.xyz` | 坐标轨迹（Å），每 10 步一帧 |
| `water32_nvt-vel-1.xyz`、`water32_nvt-frc-1.xyz` | 速度（bohr/au_time）、力（Hartree/bohr） |
| `water32_nvt-1.restart`（及 `.bak-*`） | 运动重启文件（ASCII，本身就是完整的 CP2K 输入） |
| `water32_nvt-1_1000.restart` 等 | `RESTART_HISTORY` 每 1000 步保留一份 |
| `water32_nvt-RESTART.wfn` | 波函数重启文件（二进制） |

快速检查（跳过前 1000 步平衡段；`$4` 为温度，`$6` 为守恒量）：

```bash
# mean temperature and its fluctuation
awk '!/^#/ && $1>=1000 {n++; s+=$4; ss+=$4*$4} END {m=s/n; printf "T = %.1f +- %.1f K (%d steps)\n", m, sqrt(ss/n-m*m), n}' water32_nvt-1.ener
# drift of the conserved quantity in the NVE run, per atom and per ps (96 atoms)
awk '!/^#/ {if (!f) {c0=$6; t0=$2; f=1}; c=$6; t=$2} END {printf "drift = %.2e Ha/atom/ps\n", (c-c0)/96/((t-t0)/1000)}' water32_nve-1.ener
# SCF iterations per MD step should stay small (ASPC working)
grep -c "SCF run converged" water32_nve.out
grep "MD| Energy drift per atom" water32_nve.out | tail -n 3
```

- **温度**：平衡后的平均温度应接近 300 K；96 个原子的正则涨落约为 $T\sqrt{2/(3N)}\approx 25$ K，这是正常现象；
- **能量守恒（NVE）**：守恒量应无明显的系统漂移；常用的经验判据是漂移小于约 1 meV/atom/ps（≈ 3.7×10⁻⁵ Ha/atom/ps），具体要求取决于研究目标；
- **NPT 密度**：由 `water32_npt-1.cell` 最后一列（体积，Å³）换算，$\rho = 32\times18.015/(N_A V)$；小盒子的体积涨落很大，至少要数 ps 才能得到有意义的平均；PBE-D3 水的密度与实验 0.997 g/cm³ 存在偏差，结果以你的计算为准；
- **RDF**：O–O 径向分布函数的第一峰位置约 2.8 Å（实验值约 2.8 Å）；纯 PBE 在 300 K 下普遍过度结构化（峰过高），这是泛函本身的已知问题；
- **REFTRAJ**：`water32_reftraj-1.ener` 的 `Pot.` 列给出 10 帧的新能量，`water32_reftraj-frc-1.xyz` 为新力。

本仓库不附带参考输出，所有数值以你的计算为准。

## 注意事项

- 32 个水分子的 2 ps AIMD 在 32 核上大约需要数小时（取决于硬件），本地测试请先把 `STEPS` 改成几十步；
- `water32.xyz` 没有氢键网络，前 1–2 ps 温度会有明显过冲，第 1 步用强耦合（50 fs）就是为此；统计分析请丢弃平衡段；
- 续算时不要改变原子顺序、元素或基组/赝势的 q 值，否则 `.wfn` 不能用；`.restart` 与 `.wfn` 都来自同一次运行；
- 2026.2 删除了 REFTRAJ 的 `EVAL_ENERGY_FORCES`/`EVAL_FORCES`，旧输入要改成 `EVAL ENERGY_FORCES`（`EVAL` 至少自 2024.1 起已经存在）；
- 不要在 MD 中使用 `IGNORE_CONVERGENCE_FAILURE`：未收敛步的力是错的，会悄悄破坏整条轨迹。
