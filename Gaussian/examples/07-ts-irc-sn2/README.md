# 示例 07：Cl⁻ + CH₃Br S<sub>N</sub>2 反应的 TS 优化、IRC 与端点验证

## 目的

用一个只有 6 个原子的阴离子体系，完整走一遍过渡态研究的标准流程：

1. `Opt=(TS,CalcFC,NoEigenTest) Freq` 优化过渡态 [Cl···CH₃···Br]⁻，并确认恰好 1 个虚频；
2. `IRC=(RCFC,...)` 从 TS 出发双向积分反应路径；
3. 对 IRC 两端（离子–偶极复合物）做 `Opt Freq`，确认 0 个虚频；
4. 计算分离反应物 CH₃Br 和 Cl⁻，得到能垒的另一个参考态。

对应教程：[第7章 过渡态搜索与 IRC 验证](../../tutorials/07-ts-irc.md)。

## 文件说明

| 文件 | 说明 |
|---|---|
| `sn2_cl_ch3br_ts_calcfc.gjf` | 手工搭建的 C<sub>3v</sub> TS 初猜（C–Cl 2.32 Å、C–Br 2.47 Å、CH₃ 平面），`Opt=(TS,CalcFC,NoEigenTest,MaxCycles=100) Freq`，写出 `sn2_ts.chk` |
| `sn2_irc_both_directions.gjf` | 用 `%oldchk=sn2_ts.chk` 读取 TS 结构和力常数，`IRC=(RCFC,Phase=(1,3),MaxPoints=80,StepSize=10,MaxCycle=50)`。`Phase=(1,3)` 把“C(1)–Br(3) 变长”定义为正向，即正向走向产物 |
| `irc_endpoint_opt_freq.gjf` | 两个 `--Link1--` 作业：反应物侧复合物 Cl⁻···CH₃Br 和产物侧复合物 ClCH₃···Br⁻ 的 `Opt Freq`。**坐标是手工搭建的近似结构，正式使用时请替换为你自己 IRC 的端点** |
| `sn2_separated_reactants.gjf` | 分离反应物：CH₃Br 的 `Opt Freq` 和 Cl⁻ 的 `Freq`（单原子只有平动和电子贡献） |
| `run_ts_irc.sh` | Slurm 脚本：依次运行以上 4 个作业，在 IRC 之前检查 TS 是否恰好 1 个虚频 |

原子编号（ModRedundant、`Phase` 中会用到）：1 = C，2 = Cl，3 = Br，4–6 = H。

理论水平：B3LYP-D3(BJ)/6-31+G(d)。体系是阴离子，所以加了弥散函数；这是教学用的设置，正式研究请做泛函和基组测试。

## 运行步骤

```bash
cd Gaussian/examples/07-ts-irc-sn2

# 方式 1：本地或交互节点逐个运行
g16 < sn2_cl_ch3br_ts_calcfc.gjf  > sn2_cl_ch3br_ts_calcfc.log
grep "imaginary frequencies" sn2_cl_ch3br_ts_calcfc.log      # 应恰好 1 个
g16 < sn2_irc_both_directions.gjf > sn2_irc_both_directions.log
g16 < irc_endpoint_opt_freq.gjf   > irc_endpoint_opt_freq.log
g16 < sn2_separated_reactants.gjf > sn2_separated_reactants.log

# 方式 2：Slurm 集群（先按集群情况修改脚本中的 module/g16root）
sbatch run_ts_irc.sh
```

IRC 作业依赖 TS 作业生成的 `sn2_ts.chk`，因此两者必须在同一目录中运行，并且 TS 作业要先完成。

## 预期结果 / 如何检查

| 检查项 | 命令 | 期望 |
|---|---|---|
| TS 作业正常结束 | `grep -c "Normal termination" sn2_cl_ch3br_ts_calcfc.log` | `2`（Opt 和 Freq 各一次） |
| TS 虚频个数 | `grep "imaginary frequencies" sn2_cl_ch3br_ts_calcfc.log` | `1 imaginary frequencies` |
| 虚频模式 | GaussView → Results → Vibrations 播放负频率 | CH₃ 沿 Cl–C–Br 轴伞形翻转，Cl 靠近、Br 离开 |
| IRC 汇总表 | `grep -A 200 "Summary of reaction path following" sn2_irc_both_directions.log` | TS 位于表中间（能量和坐标均为 0），两侧能量单调下降 |
| IRC 两端 | GaussView → Results → IRC/Path，查看两端结构 | 正向（Path 1）C–Br 变长、C–Cl 变短；反向相反 |
| 端点 | `grep "imaginary frequencies" irc_endpoint_opt_freq.log` | 无输出（两个复合物都是 0 个虚频） |
| 自由能 | `grep "Sum of electronic and thermal Free Energies" *.log` | 用于计算 ΔG‡ 和 ΔG<sub>r</sub> |

所有能量结果请以你自己的计算为准。本仓库没有附带输出文件，也不提供参考数值。

作图与汇总可用仓库脚本：[`scan_irc_profile.py`](../../scripts/scan_irc_profile.py)（IRC 能量曲线）、[`g16_log_summary.py`](../../scripts/g16_log_summary.py)（终止状态、虚频、热化学汇总）、[`log2xyz.py`](../../scripts/log2xyz.py)（导出 IRC 路径或端点结构），均用 `-h` 查看用法，例如 `python ../../scripts/scan_irc_profile.py -h`。

## 注意事项

- 气相 S<sub>N</sub>2 的势能面是“双井”形状：两侧离子–偶极复合物的能量低于分离反应物。计算 ΔG‡ 时需要说明参考态是复合物还是分离反应物（见教程第 8 节）。
- 双分子反应换算到 1 mol/L 标准态时，ΔG‡(1 M) = ΔG‡(1 atm) − 1.89 kcal/mol（298.15 K）。
- 如果 IRC 报 `Max corrector steps exceded`（L123），先把 `MaxCycle` 加大，再尝试在 `IRC=(...)` 中加入 `LQA`，或改用 `GS2` 并减小 `StepSize`。
- `%nprocshared=8`、`%mem=8GB` 与 `run_ts_irc.sh` 中的 `--cpus-per-task=8`、`--mem=10G` 配套，修改时要同步调整。
- 本示例的输入文件按照 Gaussian 16 官方文档的语法编写，但本仓库没有附带运行结果。如果在你的版本上遇到问题，请对照教程的“常见问题与报错”一节排查。

[← 返回第7章](../../tutorials/07-ts-irc.md)
