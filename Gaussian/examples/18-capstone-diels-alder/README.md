# 示例 18b：完整研究流程——丁二烯 + 乙烯 Diels–Alder 反应的自由能剖面

## 目的

用一个经典的周环反应（1,3-丁二烯 + 乙烯 → 环己烯）走完第18章的完整研究流程：

1. CREST（GFN2-xTB/ALPB 甲苯）搜索丁二烯的构象，确定最低能构象（s-trans），把各构象转成 DFT 输入；
2. B3LYP-D3(BJ)/def2-SVP + SMD(甲苯) 优化 + 频率：丁二烯、乙烯、环己烯（0 个虚频）；
3. 同一水平下搜索过渡态（恰好 1 个虚频），再做双向 IRC 确认它连接反应物和产物；
4. ωB97X-D/def2-TZVP + SMD(甲苯) 高精度单点；
5. GoodVibes：Grimme 准谐振（qRRHO）熵校正 + 1 mol/L 标准态 + 高精度单点能，得到 ΔG‡ 和 ΔG<sub>r</sub>；
6. 画自由能剖面图。

整套记号为 ωB97X-D/def2-TZVP(SMD)//B3LYP-D3(BJ)/def2-SVP(SMD)。

对应教程：[第18章 常见报错排查大全与完整研究流程案例](../../tutorials/18-errors-capstone.md)（第 8–9 节）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `butadiene_strans.xyz` | 丁二烯起始结构（RDKit/MMFF 生成），CREST 的输入 |
| `crest_conformers.sh` | 运行 `crest ... --gfn2 --alpb toluene --imtdgc`，并把 `crest_conformers.xyz` 中每个构象写成 `conf_NN_opt.gjf`（与本流程相同的 DFT 水平） |
| `reactant_opt.gjf` | 丁二烯（s-trans）`Opt Freq`，写 `reactant_opt.chk` |
| `reactant2_opt.gjf` | 乙烯 `Opt Freq` |
| `ts_opt.gjf` | 过渡态：手工搭建的 C<sub>s</sub> 初猜（两根正在形成的 C–C 键均为 2.27 Å），`Opt=(TS,CalcFC,NoEigenTest,MaxCycles=100) Freq` |
| `irc.gjf` | `%oldchk=ts_opt.chk`，`IRC=(RCFC,Phase=(1,6),MaxPoints=50,StepSize=10,MaxCycle=50)`；`Phase=(1,6)` 规定 C1···C6 距离变长为正向（走向反应物） |
| `product_opt.gjf` | 环己烯 `Opt Freq`，原子顺序与 TS 相同（C1–C6，然后是 H），便于与 IRC 端点对比 |
| `sp_highlevel.gjf` | 高精度单点模板（写的是 TS：`%oldchk=ts_opt.chk`，`Geom=AllCheck Guess=Read`），ωB97X-D/def2-TZVP + SMD(甲苯) |
| `make_sp_inputs.sh` | 由模板生成 4 个单点输入 `<root>_TZ.gjf`（GoodVibes 用 `--spc TZ` 按文件名配对 `X.log` 与 `X_TZ.log`） |
| `goodvibes_run.sh` | 检查所有作业正常结束、虚频个数正确（极小点 0、TS 1），然后运行 GoodVibes：`--spc TZ --temp 298.15 --conc 1.0 -f 100 --qs grimme`，并按 `da_profile.yaml` 输出反应剖面 |
| `da_profile.yaml` | GoodVibes ≥ 4.4 的 reaction-profile 文档：R = 丁二烯 + 乙烯，TS，P = 环己烯 |
| `plot_profile.py` | 不依赖 GoodVibes 的简单剖面图脚本：`python plot_profile.py R=0 TS=<数值> P=<数值> -o da_profile.png` |
| `run_workflow_slurm.sh` | 一个 Slurm 作业中按顺序运行全部 Gaussian 计算，每一步检查失败即停止（TS 不是恰好 1 个虚频时不做 IRC） |

原子编号（TS、IRC、产物一致）：C1、C4 为丁二烯末端碳，C2、C3 为丁二烯中间碳，C5、C6 为乙烯碳；新形成的键为 C1–C6 和 C4–C5。

所有 `.gjf` 都**没有**写 `%nprocshared`/`%mem`：在集群上由作业脚本设置 `GAUSS_CDEF`/`GAUSS_MDEF`（见第17章）；在工作站上手动运行时，先执行例如 `export GAUSS_CDEF=0-7 GAUSS_MDEF=12GB`，否则 Gaussian 只用 1 个核和 800 MB 内存。

## 运行步骤

```bash
cd Gaussian/examples/18-capstone-diels-alder

# 0) optional: conformer search of the diene (needs crest + xtb)
NP=8 bash crest_conformers.sh
#    -> compare the DFT energies of conf_NN_opt; the lowest one replaces reactant_opt.gjf

# 1)-5) all Gaussian steps in one Slurm job (skips steps that already finished)
sbatch run_workflow_slurm.sh

#    ... or by hand, in this order:
# g16 < reactant_opt.gjf  > reactant_opt.log
# g16 < reactant2_opt.gjf > reactant2_opt.log
# g16 < product_opt.gjf   > product_opt.log
# g16 < ts_opt.gjf        > ts_opt.log      # check: exactly 1 imaginary frequency
# g16 < irc.gjf           > irc.log
# bash make_sp_inputs.sh
# for r in reactant_opt reactant2_opt ts_opt product_opt; do g16 < ${r}_TZ.gjf > ${r}_TZ.log; done

# 6) thermochemistry and profile (pip install goodvibes; >= 4.4 for the YAML profile)
bash goodvibes_run.sh

# 7) figure (type in the relative qh-G values printed by GoodVibes)
python plot_profile.py R=0 TS=<dG_TS> P=<dG_P> -o da_profile.png
```

## 预期结果 / 如何检查

| 检查项 | 命令 | 期望 |
|---|---|---|
| 极小点 | `grep "imaginary frequencies" reactant_opt.log reactant2_opt.log product_opt.log` | 无输出 |
| TS 虚频 | `grep "imaginary frequencies" ts_opt.log` | `1 imaginary frequencies` |
| TS 虚频模式 | GaussView 打开 `ts_opt.log` → Results → Vibrations，播放负频率 | 两根 C···C 键同步伸缩（协同、近同步的 [4+2] 成键） |
| TS 结构 | GaussView 测量 C1–C6、C4–C5 | 两者相等（C<sub>s</sub> 对称），约 2.2–2.3 Å |
| IRC 两端 | `grep -A 60 "Summary of reaction path following" irc.log`；GaussView → Results → IRC/Path | 一端 C1–C6、C4–C5 都变长（s-cis 丁二烯 + 乙烯），另一端成键（环己烯） |
| 单点 | `grep "SCF Done" *_TZ.log` | 4 个能量 |
| 自由能 | `goodvibes_profile.txt` | 列出 ΔG(T) 与 Δqh-G(T)（kcal/mol），TS 相对反应物为正，产物为负 |

参考：该反应在气相中的实验活化能约为 27.5 kcal/mol（Rowley & Steiner, *Discuss. Faraday Soc.* **1951**, 10, 198，DOI: [10.1039/DF9511000198](https://doi.org/10.1039/DF9511000198)）。注意实验 E<sub>a</sub> 与计算的 ΔG‡ 不是同一个量：对气相双分子反应 E<sub>a</sub> = ΔH‡ + 2RT，而 ΔG‡ 还包含较大的 −TΔS‡ 项和标准态的选择；本例又是在 SMD(甲苯) 中计算的。比较时请按教程第 8.11 节换算。本仓库没有附带计算结果，所有数值请以你自己的计算为准。

## 注意事项

- **构象**：丁二烯的最低能构象是 s-trans，而参与反应的是 s-cis（gauche）构象。计算能垒必须以最低能构象（s-trans）为参考，否则会低估能垒。IRC 反应物端得到的是 s-cis 丁二烯与乙烯的弱复合物，不能直接作为参考态。
- **一致性**：所有物种的优化、频率和单点必须用完全相同的方法、基组、溶剂模型和格点；不要混用不同水平的热校正。
- **标准态**：`--conc 1.0` 对每个物种加上 RT ln(24.46)（298.15 K 时 +1.89 kcal/mol）。双分子反应 A + B → TS 中，反应物有两个物种，TS 只有一个，所以 ΔG‡(1 M) = ΔG‡(1 atm) − 1.89 kcal/mol。
- **频率缩放**：GoodVibes 会根据识别到的方法和基组自动选择频率缩放因子；想关闭缩放用 `-v 1.0`，并在文中说明。
- **TS 初猜失败**时，可以用第8章的方法做柔性扫描：同时缩短 C1–C6 与 C4–C5（`Geom=GIC` 定义两键之和），取能量最高点作为初猜；或者用 QST2（反应物与产物原子编号一致）。
- 仓库脚本可替代或辅助：[`thermo_qrrho.py`](../../scripts/thermo_qrrho.py)（教学版 qRRHO + 1 M 校正 + Boltzmann 加权）、[`scan_irc_profile.py`](../../scripts/scan_irc_profile.py)（IRC 能量曲线）、[`log2xyz.py`](../../scripts/log2xyz.py)（导出 IRC 端点结构）、[`g16_log_summary.py`](../../scripts/g16_log_summary.py)（汇总所有 log），均用 `-h` 查看用法。

[← 返回第18章](../../tutorials/18-errors-capstone.md)
