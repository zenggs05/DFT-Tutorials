# 示例 10-uvvis-emission-coumarin：香豆素的吸收、荧光与磷光（TD-DFT）

配套教程：[第10章 激发态 TD-DFT](../../tutorials/10-tddft.md)（第 4–7 节）。

## 目的

以香豆素（coumarin，C9H6O2）为例，走完一套完整的“吸收—发射”计算流程：

1. SMD(乙醇) 中优化 S0 并做频率，再做 20 个单重态的垂直激发（UV-Vis）；
2. 优化 S1 并做解析 TD 频率，得到垂直发射能、绝热激发能和 0-0 能量；
3. 用 cLR 态特定溶剂化重新计算垂直发射能（官方 SCRF 页的方案）；
4. 用 UKS 优化 T1，计算垂直磷光能和 S-T 能隙。

## 文件说明

| 文件 | 内容 | 依赖 |
|---|---|---|
| `coumarin_td_camb3lyp_smd.gjf` | 第 1 步 S0 Opt+Freq（CAM-B3LYP/def2-SVP, SMD 乙醇）；第 2 步 TD 单点（CAM-B3LYP/def2-TZVP，NStates=20） | 无 |
| `s1_opt_freq_emission.gjf` | S1 优化 + 频率（`TD=(NStates=6,Root=1) Opt Freq NoSymm`） | `coumarin_s0.chk` |
| `s1_emission_state_specific.gjf` | 两步：S1 的 cLR 平衡溶剂化（`NonEquilibrium=Save`）→ S1 结构上的非平衡 S0（`NonEquilibrium=Read`），IEFPCM | `coumarin_s1.chk` |
| `t1_uks_opt.gjf` | 两步：UKS 三重态 Opt+Freq → T1 结构上的 S0 单点 | `coumarin_s0.chk` |
| `run.slurm` | 按依赖顺序依次运行以上 4 个输入，并对所有 chk 执行 formchk | — |

初始结构由 RDKit（ETKDG + MMFF）生成，只作为优化起点。

## 运行步骤

```bash
# 集群（Slurm）：先按集群修改 run.slurm 中的 module / g16root / scratch 路径
sbatch run.slurm

# 或者在交互节点上手动按顺序运行
g16 < coumarin_td_camb3lyp_smd.gjf  > coumarin_td_camb3lyp_smd.log
g16 < s1_opt_freq_emission.gjf      > s1_opt_freq_emission.log
g16 < s1_emission_state_specific.gjf > s1_emission_state_specific.log
g16 < t1_uks_opt.gjf                > t1_uks_opt.log
```

`s1_opt_freq_emission.gjf` 的 TD 解析频率最耗时（8 核上通常需要数小时，视硬件而定）。

## 预期结果 / 如何检查

```bash
# 1) 所有作业正常结束
grep -c "Normal termination" *.log
# 2) S0、S1、T1 优化后都没有虚频（归档段 NImag=0）
grep "NImag" coumarin_td_camb3lyp_smd.log s1_opt_freq_emission.log t1_uks_opt.log
# 3) 吸收：激发能、波长、振子强度
grep "Excited State" coumarin_td_camb3lyp_smd.log | awk '{print $3, $4, $5, $7, $9}'
# 4) 发射：S1 优化最后一步的第 1 激发态 = LR 垂直发射能
grep "Excited State   1:" s1_opt_freq_emission.log | tail -1
grep "Total Energy, E(TD-HF/TD-DFT)" s1_opt_freq_emission.log | tail -1
# 5) cLR 发射：两步能量差
grep -A3 "PCM Corrected Linear-Response" s1_emission_state_specific.log
grep "SCF Done" s1_emission_state_specific.log | tail -1
# 6) T1：<S**2> 应接近 2.0
grep "S\*\*2" t1_uks_opt.log | tail -3
```

需要自己计算的量（单位换算 1 Hartree = 27.2114 eV，λ(nm) = 1239.84 / E(eV)）：

- 绝热激发能 = S1 优化最后的 `Total Energy, E(TD-HF/TD-DFT)` − S0 优化最后的 `SCF Done`（同为 CAM-B3LYP/def2-SVP + SMD）；
- 0-0 能量 = 绝热激发能 + [ZPE(S1) − ZPE(S0)]，ZPE 取各自的 `Zero-point correction`；
- Stokes 位移 = 垂直吸收能 − 垂直发射能（建议换算为 cm⁻¹）；
- 垂直磷光能 ≈ T1 优化最后的 `SCF Done` − `t1_uks_opt.log` 第 2 步的 `SCF Done`。

光谱展宽可使用 GaussView（Results → UV-VIS）、Multiwfn 主功能 11，或仓库脚本（`python ../../scripts/spectrum_uvvis_ecd.py -h`）。

## 注意事项

- 吸收单点用 def2-TZVP，而优化用 def2-SVP：比较绝热量时只能用同一基组的能量，不要把 def2-TZVP 的激发能和 def2-SVP 的总能量混用。
- `s1_emission_state_specific.gjf` 按官方示例使用 IEFPCM，而 S1 结构来自 SMD 优化；在方法部分说明这一点即可。
- 如果 S1 优化中第 1 激发态的主要组态或振子强度突然改变，可能发生了根翻转，请参考教程第 6.2 节。
- 所有输入假定 `%mem=8GB`、`%nprocshared=8`，与 `run.slurm` 中的资源申请一致；修改时两边要同步。

返回：[第10章 激发态 TD-DFT](../../tutorials/10-tddft.md)
