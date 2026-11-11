# 示例 09：SMD 隐式溶剂——溶液中优化、溶剂化自由能、自定义溶剂、pKa 与氧化还原电位

## 目的

用苯酚、苯氧负离子和四硫富瓦烯（TTF）演示 Gaussian 16 中隐式溶剂模型的典型用法：

1. 在 SMD(水) 中做 Opt Freq，得到溶液中的自由能；
2. 用“气相优化 + SMD 单点”计算溶剂化自由能 ΔG_solv；
3. 用 `Solvent=Generic,Read` 自定义溶剂，并用内置溶剂做对照，验证写法；
4. 用苯酚/苯氧负离子计算 pKa（直接法）；
5. 用 TTF/TTF•⁺ 在乙腈中计算氧化电位（低水平 Opt Freq + 高水平单点）。

对应教程：[第9章 隐式溶剂模型 PCM/SMD](../../tutorials/09-solvation.md)。

## 文件说明

| 文件 | 作业 | 理论水平 |
|---|---|---|
| `phenol_smd_water_opt_freq.gjf` | 苯酚，SMD(水)，Opt Freq | M06-2X/6-311+G(d,p) |
| `phenoxide_smd_water_opt_freq.gjf` | 苯氧负离子（−1 1），SMD(水)，Opt Freq | M06-2X/6-311+G(d,p) |
| `solvation_free_energy_gas_and_smd.gjf` | Link1：① 气相 Opt Freq；② 同一结构 SMD(水) 单点（`Geom=Check Guess=Read`） | M05-2X/6-31G(d)（SMD 参数化所用水平之一） |
| `custom_solvent_generic_read.gjf` | Link1：① 内置 `Solvent=Ethanol`；② `Solvent=Generic,Read` 手动输入乙醇的 7 个 SMD 参数 | M05-2X/6-31G(d) |
| `redox_potential_cycle_neutral_and_cation.gjf` | Link1 四步：① TTF Opt Freq；② TTF•⁺（1 2）Opt Freq；③ TTF 单点；④ TTF•⁺ 单点，全部在 SMD(乙腈) 中 | ①② B3LYP-D3(BJ)/def2-SVP；③④ M06-2X/def2-TZVP |
| `run_solvation.sh` | Slurm 脚本：依次运行全部输入，最后汇总 SCF 能量和自由能 | — |

初始坐标由 RDKit 嵌入 + MMFF 优化得到，只作为 DFT 优化的起点。

## 运行步骤

```bash
cd Gaussian/examples/09-solvation-smd

g16 < phenol_smd_water_opt_freq.gjf        > phenol_smd_water_opt_freq.log
g16 < phenoxide_smd_water_opt_freq.gjf     > phenoxide_smd_water_opt_freq.log
g16 < solvation_free_energy_gas_and_smd.gjf > solvation_free_energy_gas_and_smd.log
g16 < custom_solvent_generic_read.gjf      > custom_solvent_generic_read.log
g16 < redox_potential_cycle_neutral_and_cation.gjf > redox_potential_cycle_neutral_and_cation.log

# 或者在 Slurm 集群上
sbatch run_solvation.sh
```

## 预期结果 / 如何检查

**通用检查**

```bash
grep -c "Normal termination" *.log                 # Opt Freq 为 2，Link1 作业按步数累计
grep "imaginary frequencies" *.log                 # 极小点不应出现这一行
grep "SMD-CDS" *.log                               # 确认 SMD 的非静电项被计算
```

**ΔG_solv（苯酚，水）**

```bash
awk '/SCF Done/{e=$5} /Normal termination/{n++; E[n]=e}
     END{printf "dG_solv = %.2f kcal/mol\n", (E[n]-E[1])*627.509}' solvation_free_energy_gas_and_smd.log
```

苯酚的实验水合自由能约为 −6.6 kcal/mol，SMD 对中性分子的典型误差约为 1 kcal/mol。

**自定义溶剂**：`grep "SCF Done" custom_solvent_generic_read.log` 得到的两个能量应当一致（差别在 10⁻⁶ Hartree 量级以内）。如果不一致，说明参数输入或格式有误。

**pKa（直接法）**：

```bash
grep "Sum of electronic and thermal Free Energies" phenol_smd_water_opt_freq.log phenoxide_smd_water_opt_freq.log
```

然后按教程第 8 节的公式计算：ΔG = [G(A⁻) − G(HA)] × 627.509 + (−270.3) kcal/mol，pKa = ΔG / 1.364。苯酚的实验 pKa 约为 10.0。直接法误差常达几个 pKa 单位，这是方法本身的局限，不一定说明计算出错。

**氧化电位（TTF，乙腈）**：

- 从步骤 ①② 的日志中取 `Thermal correction to Gibbs Free Energy`，从 ③④ 取 `SCF Done`；
- G(X) = E<sub>③/④</sub> + 热校正<sub>①/②</sub>，ΔG<sub>ox</sub> = G(TTF•⁺) − G(TTF)；
- 换算为 eV（× 27.2114）后，E°(vs SHE) = ΔG<sub>ox</sub>(eV) − 4.28 V（或 4.44 V，需注明采用哪个值）；更推荐在同一水平下计算二茂铁，报告相对 Fc⁺/Fc 的电位；
- 检查 TTF•⁺ 的自旋污染：`grep "S\*\*2" redox_potential_cycle_neutral_and_cation.log`，`<S**2>` 应接近 0.75。

仓库脚本 `python ../../scripts/g16_log_summary.py -h`、`python ../../scripts/thermo_qrrho.py -h` 可以批量提取上述数据。本仓库不附带输出文件，所有数值以你自己的计算为准。

## 注意事项

- Gaussian 的热化学对应 1 atm。溶液中的自由能换算到 1 mol/L 时，每个溶质加 1.89 kcal/mol（298.15 K）。pKa 直接法中 HA 和 A⁻ 的校正相互抵消，H⁺ 的值已包含该校正；氧化还原过程两边溶质个数相同，校正也相互抵消。
- `custom_solvent_generic_read.gjf` 第二步的分子说明只有一行 `0 1`（结构取自 chk），后面空一行，然后是 PCM 参数段，段末再空一行。空行写错是最常见的读入错误。
- SMD 自定义溶剂的 7 个参数全部必需，并且都必须写成实数。`SurfaceTensionAtInterface` 的单位是 cal mol⁻¹ Å⁻²（1 dyn/cm = 1.4393 cal mol⁻¹ Å⁻²）。
- `%nprocshared=8`、`%mem=8GB` 与 `run_solvation.sh` 中的申请（8 核、10G）配套。

[← 返回第9章](../../tutorials/09-solvation.md)
