# 示例 06：频率计算、红外拉曼与热化学（Gaussian 16）

配套教程：[第6章 频率计算、红外拉曼与热化学](../../tutorials/06-freq-thermo.md)

## 目的

用几个小分子练习 Gaussian 16 频率计算和热化学分析的完整流程：

- `Opt Freq=Raman`：谐振频率、IR 强度和 Raman 活性；
- `Temperature=`/`Pressure=`：改变热化学的温度和压力；
- `Freq=(ReadFC,ReadIsotopes)` 与 `H(Iso=2)`：两种同位素替换写法；
- `Freq=Anharmonic`：VPT2 非谐频率；
- `Freq=ReadFC` + `Pressure=24.4654`：不重算 Hessian，直接得到 1 mol/L 标准态的自由能；
- 低水平 `Opt Freq` + 高水平单点（`--Link1--`），再用 GoodVibes 做 qRRHO 和 1 M 校正。

## 文件说明

| 文件 | 体系 | 演示内容 | 依赖 |
|---|---|---|---|
| `h2o_opt_freq_raman.gjf` | H2O | `Opt Freq=Raman`，B3LYP-D3(BJ)/def2-TZVP | — |
| `freq_temperature_pressure.gjf` | 甲醇 | `Opt Freq Temperature=373.15 Pressure=2.0` | — |
| `d2o_freq_readisotopes.gjf` | D2O | `Freq=(ReadFC,ReadIsotopes) Geom=AllCheck`，从 H2O 的 chk 读力常数 | 需要 `h2o_opt_freq_raman.chk` |
| `d2o_freq_iso_inline.gjf` | D2O | 坐标中写 `H(Iso=2)`，从头 `Opt Freq` | — |
| `formaldehyde_freq_anharmonic.gjf` | 甲醛 | `Opt=VeryTight Int=SuperFine Freq=Anharmonic` | — |
| `freq_readfc_rethermo.gjf` | H2O | `Freq=ReadFC Geom=AllCheck Pressure=24.4654`，即 298.15 K 下 1 mol/L 标准态 | 需要 `h2o_opt_freq_raman.chk` |
| `highlevel_sp_plus_thermo.gjf` | 乙醇 | 第 1 步 B3LYP-D3(BJ)/def2-SVP `Opt Freq`，第 2 步 ωB97X-D/def2-TZVP 单点 | — |
| `goodvibes_run.sh` | — | GoodVibes 命令示例（qRRHO、1 M、关闭缩放、单点组合、Truhlar 准谐振） | 需要上面各作业的 log 和 GoodVibes |
| `run_g16.sh` | — | 最简 Slurm 作业脚本（8 核、16 GB） | — |

起始结构：H2O、D2O、甲醛按常用实验近似几何参数手工搭建；甲醇和乙醇由 RDKit（ETKDG + MMFF）生成。

## 运行步骤

```bash
cd Gaussian/examples/06-freq-thermo

# 1) independent jobs
for f in h2o_opt_freq_raman freq_temperature_pressure d2o_freq_iso_inline \
         formaldehyde_freq_anharmonic highlevel_sp_plus_thermo; do
    g16 < $f.gjf > $f.log
done

# 2) jobs that read h2o_opt_freq_raman.chk (run after step 1)
g16 < d2o_freq_readisotopes.gjf > d2o_freq_readisotopes.log
g16 < freq_readfc_rethermo.gjf  > freq_readfc_rethermo.log

# 3) quasi-harmonic thermochemistry (pip install goodvibes)
bash goodvibes_run.sh
```

Slurm：`sbatch run_g16.sh h2o_opt_freq_raman.gjf`，依赖 chk 的作业要等前一个作业结束后再提交（或用 `sbatch --dependency=afterok:<jobid>`）。`formaldehyde_freq_anharmonic` 要沿每个简正模式做数值微分（甲醛约十几次 Hessian 计算），比其他作业慢得多。

## 预期结果 / 如何检查

```bash
# termination (Opt Freq jobs: 2 lines; ReadFC jobs: 1 line; Link1 job: 3 lines)
grep -c "Normal termination" *.log

# minima: no imaginary frequencies
grep "imaginary frequencies" *.log

# harmonic frequencies, IR and Raman
grep -E "Frequencies --|IR Inten|Raman Activ" h2o_opt_freq_raman.log

# isotopes: both D2O jobs should give the same frequencies
grep "Frequencies --" d2o_freq_readisotopes.log d2o_freq_iso_inline.log

# thermochemistry block, temperature and pressure actually used
grep -E "Temperature .* Kelvin|Rotational symmetry number" *.log
grep -E "Zero-point correction|Thermal correction to Gibbs|Sum of electronic and thermal Free" *.log

# anharmonic results
grep -A 12 "Fundamental Bands" formaldehyde_freq_anharmonic.log
```

判断标准：

| 作业 | 应该看到 |
|---|---|
| H2O | 3 个正的频率，次序为弯曲（A1）、对称伸缩（A1）、反对称伸缩（B2）；有 `Raman Activ` 行；`Rotational symmetry number  2.` |
| 甲醇 | `Temperature   373.150 Kelvin.  Pressure   2.00000 Atm.` |
| 两个 D2O 作业 | 频率相同（数值误差范围内），并且明显低于 H2O；ReadIsotopes 作业几乎不耗时 |
| `freq_readfc_rethermo` | 与 H2O 作业相比，`Thermal correction to Gibbs Free Energy` 大约增加 0.00302 Eh（约 1.89 kcal/mol） |
| 甲醛 | 有 `Anharmonic Zero Point Energy` 和 `Fundamental Bands` 表，`E(anharm)` 一般低于 `E(harm)` |
| 乙醇 Link1 | 第 1 步给出热校正，第 2 步给出 ωB97X-D 的 `SCF Done`；GoodVibes 的 `--spc link` 输出中带 `_SPC` 后缀的列 |

参考：气相实验基频 H2O 约 1595/3657/3756 cm⁻¹，D2O 约 1178/2671/2788 cm⁻¹，甲醛约 2782、1746、1500、1167、2843、1249 cm⁻¹（NIST CCCBDB 实验振动数据）。谐振计算值一般比这些数值高几个百分点。

本目录不提供计算输出，具体数值请以自己的结果为准。

## 注意事项

- 频率必须与优化使用同一理论水平。所有示例都用 `Opt Freq` 写在同一 route 中，或者用 `Freq=ReadFC` 读取同一水平的力常数。
- `Geom=AllCheck` 的输入在 route 之后**不能**再写 title 和电荷多重度；ReadIsotopes 的输入节紧跟在 route 后的空行之后。
- 用 `ReadFC` 时基组从 chk 读取，不要再写 Gen 基组输入节；上一作业已经用了 `Freq=Raman` 时，ReadFC 作业中不要再写 Raman。
- `Pressure=24.4654` 只在 298.15 K 下对应 1 mol/L；换温度时按 P = cRT 重新计算。GoodVibes 的 `--conc 1.0` 与之等价，两者**不要同时使用**。
- GoodVibes v4 会按识别出的理论水平自动套用频率缩放因子；与 Gaussian 原始输出对照时加 `-v 1.0`。
- 本仓库的 [thermo_qrrho.py](../../scripts/thermo_qrrho.py)（`python ../../scripts/thermo_qrrho.py -h`）和 [spectrum_ir_raman.py](../../scripts/spectrum_ir_raman.py)（`python ../../scripts/spectrum_ir_raman.py -h`）分别可以做教学版 qRRHO 热化学和 IR/Raman 光谱绘制。
