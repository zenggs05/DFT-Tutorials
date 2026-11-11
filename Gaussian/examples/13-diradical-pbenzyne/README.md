# 示例 13-2：对苯炔单重态双自由基——Stable=Opt、Guess=Mix 与 UNO

配套教程：[第13章 开壳层体系、波函数稳定性与 Broken-Symmetry](../../tutorials/13-openshell-stability.md)（§7、§8）

## 目的

- 从闭壳层 RB3LYP 出发，用 `Stable=Opt` 发现 RHF→UHF 不稳定性并自动得到 broken-symmetry（BS）UKS 解；
- 用 `Guess=Mix` 走另一条路线，验证两者得到同一个解；
- 在 BS 势能面上优化，并在最终结构上重新做稳定性检查；
- 打印并保存自然轨道（UNO），用占据数判断双自由基特征；
- 计算三重态，得到单重态-三重态能隙，并做近似自旋投影（AP）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `benzyne_stable_opt_guess_mix.gjf` | 6 步 `--Link1--` 作业，UB3LYP/6-31G(d)，全部 `NoSymm`（见下表） |
| `run_g16_slurm.sh` | Slurm 作业脚本（8 核、10 GB） |

| 步骤 | route 要点 | 写入的 chk |
|---|---|---|
| 1 | `RB3LYP Stable=Opt` | `pbz_rks.chk`（若不稳定，里面存的是 UKS 解） |
| 2 | `UB3LYP Guess=Mix Stable=Opt Geom=Check` | `pbz_mix.chk` |
| 3 | `UB3LYP Opt Freq Geom=Check Guess=Read`（读第 1 步） | `pbz_bs_opt.chk` |
| 4 | `UB3LYP Stable Geom=AllCheck Guess=Read Pop=NaturalOrbitals` | `pbz_bs_stab.chk` |
| 5 | `Guess=(Save,Only,NaturalOrbitals) Geom=AllCheck ChkBasis` | `pbz_uno.chk`（含自然轨道） |
| 6 | `UB3LYP Opt Freq Geom=Check`，多重度 3 | `pbz_triplet.chk` |

起始结构是按理想六边形（C–C 1.39 Å，C–H 1.08 Å）搭建的 D2h 结构，去掉了 1、4 位的氢，不是文献坐标。

## 运行步骤

```bash
cd Gaussian/examples/13-diradical-pbenzyne
g16 < benzyne_stable_opt_guess_mix.gjf > benzyne_stable_opt_guess_mix.log
# or on a Slurm cluster
sbatch run_g16_slurm.sh benzyne_stable_opt_guess_mix.gjf
```

8 核上通常十几分钟到一小时（主要耗时在两次 `Opt Freq`）。完成后把 UNO 转成 fchk 以便查看：

```bash
formchk pbz_uno.chk pbz_uno.fchk
```

## 预期结果 / 如何检查

```bash
grep -c "Normal termination" benzyne_stable_opt_guess_mix.log          # expect 6
grep -i "instability\|stable under\|already stable" benzyne_stable_opt_guess_mix.log
grep "SCF Done" benzyne_stable_opt_guess_mix.log | head -4
grep "S\*\*2 before annihilation" benzyne_stable_opt_guess_mix.log | tail -3
grep "Sum of electronic and zero-point Energies" benzyne_stable_opt_guess_mix.log
grep -c "imaginary frequencies" benzyne_stable_opt_guess_mix.log        # 0 if no imaginary modes
```

定性预期：

| 检查项 | 预期 |
|---|---|
| 第 1 步 | 先出现 `The wavefunction has an RHF -> UHF instability.`，随后 `SCF Done: E(UB3LYP)`，能量低于最初的 RB3LYP |
| 第 1、2 步 | 最终能量一致（差别在 10⁻⁶ Hartree 量级），⟨S²⟩ 一致 |
| BS 单重态 ⟨S²⟩ | 明显大于 0、小于 1，体现部分双自由基特征 |
| 第 4 步 | `stable under the perturbations considered` |
| 自然轨道 | HONO/LUNO 占据数成对出现（约 $2-n$ 与 $n$），$n$ 明显大于 0 |
| 第 6 步 | 三重态 ⟨S²⟩ 接近 2.0 |

单重态-三重态能隙用第 3 步与第 6 步的 `Sum of electronic and zero-point Energies` 之差；再按教程 §2.5 的 AP 公式修正：

$$
E_{\rm S}^{\rm AP} = E_{\rm BS} + \frac{\langle S^2\rangle_{\rm BS}}{\langle S^2\rangle_{\rm T}-\langle S^2\rangle_{\rm BS}}(E_{\rm BS}-E_{\rm T})
$$

与实验约 3.8 kcal/mol（单重态更稳定；Wenthold, Squires, Lineberger, *J. Am. Chem. Soc.* 1998, 120, 5279）比较。

## 注意事项

- 第 3 步之后的所有单重态步骤都显式写 `UB3LYP`，否则会退回闭壳层。
- 第 5 步的写法来自官方 Population 关键词页面（把自然轨道写入 chk 的推荐方法）；如果只需要占据数，看第 4 步输出即可。
- 若第 1、2 步得到不同能量的解，说明存在多个 UKS 解，选能量更低且稳定的那个继续后续步骤（修改第 3 步的 `%oldchk`）。
- 要把轨道交给 CASSCF/NEVPT2，参考教程 §10 的 MOKIT 流程；MOKIT 文档建议读入稳定的 UHF（而非 UDFT）轨道，需另做一个 UHF + `Stable=Opt` 计算。

返回教程：[第13章](../../tutorials/13-openshell-stability.md)
