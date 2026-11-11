# 示例 03：方法、泛函、色散校正与混合基组（Gen/GenECP）

配套教程：[第3章 方法、泛函与基组选择（含 Gen/GenECP 与赝势）](../../tutorials/03-method-basis.md)

## 目的

- 写出带正确色散校正的 DFT 输入：B3LYP + `GD3BJ`，M06-2X + `GD3`（零阻尼）。
- 演示 high//low 两步法：B3LYP-D3(BJ)/def2-SVP 优化 + 频率，再用双杂化 DSD-PBEP86/def2-TZVP 做单点。
- 演示含过渡金属体系的两种基组设置：
  - `GenECP` 混合基组：C/H/P/Cl 用 6-31G(d)，Pd 用 LANL2DZ 基组 + ECP；
  - 从 Basis Set Exchange（BSE）下载 Gaussian 格式的基组与 ECP 并粘贴进输入，再与内置 `Def2TZVP` 关键词的结果对比，验证粘贴是否正确。

## 文件说明

| 文件 | 体系 | route 要点 | 估计耗时（8 核，仅供参考） |
|---|---|---|---|
| `benzene_b3lyp_gd3bj_def2svp.gjf` | 苯 | `B3LYP/def2SVP EmpiricalDispersion=GD3BJ Opt Freq` | 几分钟 |
| `dsdpbep86_def2tzvp_sp.gjf` | 苯（读上一步的 chk） | `%oldchk` + `DSDPBEP86/def2TZVP Geom=AllCheck SP` | 几分钟 |
| `ethylene_m062x_gd3_def2tzvp_sp.gjf` | 乙烯 | `M062X/def2TZVP EmpiricalDispersion=GD3 SP` | 不到 1 分钟 |
| `pd_complex_genecp_lanl2dz.gjf` | trans-PdCl2(PMe3)2 | `B3LYP/GenECP 6D EmpiricalDispersion=GD3BJ Opt Freq` | 半小时到数小时 |
| `ir_complex_gen_def2tzvp_from_bse.gjf` | trans-IrCl(CO)(PH3)2（Vaska 配合物的简化模型） | `PBE1PBE/GenECP EmpiricalDispersion=GD3BJ SP`；轻原子用内置 `Def2TZVP`，Ir 的基组和 ECP 从 BSE 粘贴 | 几分钟 |
| `ir_complex_def2tzvp_builtin.gjf` | 同上 | `PBE1PBE/Def2TZVP EmpiricalDispersion=GD3BJ SP`（内置关键词，ECP 自动加上） | 几分钟 |

结构来源：苯和乙烯取自 G2 数据集的几何（通过 ASE 获得）；两个金属配合物是按典型键长（Pd–Cl 2.31、Pd–P 2.33、Ir–Cl 2.38、Ir–C 1.82、C–O 1.15、Ir–P 2.33 Å）手工搭建的近似结构，只作为计算起点。

Ir 的基组和 ECP 数据下载自 BSE（def2-TZVP，Gaussian94 格式，BSE 0.12，数据来源 Turbomole 7.3），去掉了文件头的注释和空行后原样粘贴。

## 运行步骤

```bash
cd Gaussian/examples/03-method-basis

# 1) high//low pair: run the def2-SVP opt+freq FIRST, then the double hybrid
g16 < benzene_b3lyp_gd3bj_def2svp.gjf > benzene_b3lyp_gd3bj_def2svp.log
g16 < dsdpbep86_def2tzvp_sp.gjf      > dsdpbep86_def2tzvp_sp.log

# 2) Minnesota functional with zero-damping D3
g16 < ethylene_m062x_gd3_def2tzvp_sp.gjf > ethylene_m062x_gd3_def2tzvp_sp.log

# 3) GenECP mixed basis (longest job; submit to a cluster)
g16 < pd_complex_genecp_lanl2dz.gjf > pd_complex_genecp_lanl2dz.log

# 4) BSE-pasted basis vs. built-in keyword: the two should agree
g16 < ir_complex_gen_def2tzvp_from_bse.gjf > ir_complex_gen_def2tzvp_from_bse.log
g16 < ir_complex_def2tzvp_builtin.gjf      > ir_complex_def2tzvp_builtin.log
```

`dsdpbep86_def2tzvp_sp.gjf` 通过 `%oldchk=benzene_b3lyp_gd3bj_def2svp.chk` 读取苯的优化结构，所以**必须先运行苯的优化**，并且两个文件在同一目录下。

在集群上请用 Slurm 提交，并保持 `%nprocshared` 等于 `--cpus-per-task`、`%mem` 约为 `--mem` 的 80%（参考 `../01-first-job-h2o/run_g16_slurm.sh` 或第17章）。

## 预期结果 / 如何检查

```bash
# termination
for f in *.log; do echo "$f: $(tail -n 1 $f)"; done

# number of basis functions (compare with the table below)
grep -m1 "basis functions" *.log

# dispersion contribution printed by Gaussian
grep -m1 "Dispersion energy" benzene_b3lyp_gd3bj_def2svp.log ethylene_m062x_gd3_def2tzvp_sp.log

# double hybrid: the final energy is NOT the "SCF Done" line
grep "SCF Done" dsdpbep86_def2tzvp_sp.log
grep "E2(" dsdpbep86_def2tzvp_sp.log

# BSE-pasted vs built-in: same NBasis and (nearly) identical energy expected
grep -m1 "basis functions" ir_complex_*.log
grep "SCF Done" ir_complex_*.log
```

| 检查项 | 预期 |
|---|---|
| 苯 def2-SVP 基函数数 | 114（C [3s2p1d] = 14 个，H [2s1p] = 5 个；按标准收缩方案计算，以输出为准） |
| 苯 def2-TZVP 基函数数 | 222（C 31 个，H 6 个） |
| 乙烯 def2-TZVP 基函数数 | 86 |
| 苯的频率 | 无虚频（`grep "Frequencies --"` 第一行全为正值） |
| 双杂化输出 | `SCF Done` 给出的只是杂化泛函部分的能量；总能量在含 `E2(` 的那一行的 `E(...)=` 之后，用 Fortran D 格式（如 `-0.232D+03`）打印 |
| Ir：BSE 粘贴版 vs 内置版 | 两个 log 的 `basis functions` 数相同，SCF 能量应一致或只在数值噪声范围内有差别。如果基函数数不同，说明粘贴有误或 5D/6D 设置不一致 |
| Pd 配合物 | Opt 收敛（`Stationary point found`），频率无虚频；Pd–P、Pd–Cl 键长在合理范围（约 2.3 Å 上下，近似值） |
| ECP 是否生效（电子数） | `grep -m1 "alpha electrons"`：Pd 配合物全电子 164 个，LANL2DZ 的 ECP 替换 28 个芯电子后剩 136 个，即 alpha、beta 各 68 个；Ir 配合物全电子 144 个，def2-ECP 替换 60 个后剩 84 个，即各 42 个 |

## 注意事项

- **色散校正要与泛函匹配**：M06-2X 只能用 `GD3`，加 `GD3BJ` 会在 L301 报 `R6DS8: Unable to choose the S8 parameter`。`wB97XD`、`B97D3`、`B2PLYPD3`、`DSDPBEP86` 已自带色散，不要再加 `EmpiricalDispersion`。
- **Gen/GenECP 格式**：每个元素块以 `****` 结束；基组段与 ECP 段之间、ECP 段之后各有一个空行；分子中出现的每种元素都必须定义一次基组；Gen 段里不能列出分子中不存在的元素（否则作业中止）。
- 只要有任何原子使用 ECP，就必须用 `GenECP`（等价于 `Gen Pseudo=Read`），不能只写 `Gen`。
- `pd_complex_genecp_lanl2dz.gjf` 中的 `6D` 让 6-31G(d) 使用笛卡尔 d 函数，与内置关键词 `6-31G(d)` 的定义一致（Gen 输入默认是 5D）。是否加 `6D` 会改变基函数数和能量，**同一项目中所有要比较的计算必须保持一致**。
- 从 BSE 粘贴时删掉文件开头的空行。空行会被 Gaussian 当作基组段的结束，常见报错为 `End of file reading basis center` 或 `EOF while reading ECP pointer card`。
- 双杂化泛函的成本与 MP2 相当，内存和磁盘需求都明显大于普通杂化泛函，一般只用于单点。
