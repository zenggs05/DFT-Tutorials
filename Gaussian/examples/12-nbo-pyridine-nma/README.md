# 示例 12-nbo-pyridine-nma：Mulliken 基组依赖性、NBO、Wiberg 键级与 E(2)

配套教程：[第12章 布居分析、NBO 与静电势/RESP 电荷](../../tutorials/12-population-charges.md)（第 4–6 节）。

## 目的

1. 在吡啶上运行完整的 Mulliken 分析和内置 NBO 3.1 分析（NPA 电荷、Lewis 结构、杂化轨道）；
2. 在同一吡啶结构上用 5 个基组比较 Mulliken 与 Hirshfeld/CM5 电荷，直观感受 Mulliken 的基组依赖性；
3. 对反式 N-甲基乙酰胺（NMA）用 `Pop=NBORead` 加 `$NBO BNDIDX $END` 得到 Wiberg 键级，并从 E(2) 表中找出酰胺共振 $n_N\to\pi^*_{C=O}$。

## 文件说明

| 文件 | 内容 | 依赖 |
|---|---|---|
| `pyridine_pop_nbo.gjf` | 吡啶：B3LYP-D3(BJ)/def2-SVP Opt+Freq → B3LYP/def2-TZVP `Pop=(Full,NBO)` | 无 |
| `pyridine_mulliken_basis.gjf` | 5 个 Link1 单点：STO-3G、6-31G(d)、6-311++G(d,p)、def2-TZVP、aug-cc-pVTZ，均为 `Pop=Hirshfeld` | `pyridine_opt.chk` |
| `nboread_bndidx_e2.gjf` | NMA：Opt+Freq → B3LYP/def2-TZVP `Pop=NBORead`，`$NBO BNDIDX $END` | 无 |
| `run.slurm` | 按依赖顺序运行以上 3 个输入 | — |

原子编号：吡啶中 N 为 4 号；NMA 中 C1（乙酰基甲基）、C2（羰基碳）、O3、N4、C5（N-甲基）、H9（N-H）。

## 运行步骤

```bash
sbatch run.slurm
# 或按顺序手动运行
g16 < pyridine_pop_nbo.gjf        > pyridine_pop_nbo.log
g16 < pyridine_mulliken_basis.gjf > pyridine_mulliken_basis.log
g16 < nboread_bndidx_e2.gjf       > nboread_bndidx_e2.log
```

## 预期结果 / 如何检查

```bash
# 吡啶 NPA 电荷与 NBO 摘要
grep -A20 "Summary of Natural Population Analysis" pyridine_pop_nbo.log | head -22
grep -A40 "Natural Bond Orbitals (Summary)" pyridine_pop_nbo.log | grep "LP"

# 5 个基组下 N4 的 Mulliken 与 CM5 电荷
grep -A12 "Mulliken charges:" pyridine_mulliken_basis.log | grep -E "^ +4 +N"
grep -A12 "Hirshfeld charges, spin densities" pyridine_mulliken_basis.log | grep -E "^ +4 +N"

# NMA 的 Wiberg 键级和 E(2)
grep -A20 "Wiberg bond index matrix in the NAO basis" nboread_bndidx_e2.log | head -24
grep -A80 "SECOND ORDER PERTURBATION THEORY ANALYSIS" nboread_bndidx_e2.log | grep "LP ( 1) N   4"
```

预期现象（定性）：

- 吡啶 N 带负的 NPA 电荷，NBO 列表中有一个位于环平面内的 sp² 型 `LP ( 1) N 4`；
- N4 的 Mulliken 电荷随基组明显变化，CM5 电荷变化小得多；
- NMA 的 C2-N4 Wiberg 键级明显大于 1，`LP(1) N4 → BD*(2) C2-O3` 的 E(2) 是表中最大的项之一（通常为数十 kcal/mol）。

## 注意事项

- `pyridine_mulliken_basis.gjf` 需要 `pyridine_pop_nbo.gjf` 第 1 步生成的 `pyridine_opt.chk`。
- `$NBO` 行必须放在 `0 1` 之后的空行之后，并以空行结束；只有 `Pop=NBORead`（或 NBO7Read）才会读取它。
- E(2) 只适合在同一水平下做相对比较，不能把各项相加当作“总共轭能”。
- 输出中的轨道编号、原子编号与本 README 一致的前提是没有改动输入中的原子顺序。

返回：[第12章 布居分析、NBO 与静电势/RESP 电荷](../../tutorials/12-population-charges.md)
