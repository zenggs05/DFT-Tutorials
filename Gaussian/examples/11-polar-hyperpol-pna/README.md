# 示例 11-polar-hyperpol-pna：p-硝基苯胺的极化率与第一超极化率

配套教程：[第11章 NMR 化学位移与其他分子性质](../../tutorials/11-nmr-props.md)（第 9 节）。

## 目的

计算推-拉型非线性光学分子 p-硝基苯胺（pNA）的偶极矩、静态与 1064 nm 下的极化率 α 和第一超极化率 β，练习 `Polar CPHF=RdFreq` 的输入格式、输出读取与单位换算。

## 文件说明

| 文件 | 内容 |
|---|---|
| `pna_polar_hyperpolarizability.gjf` | 第 1 步 B3LYP-D3(BJ)/6-311+G(d,p) Opt+Freq（气相）；第 2 步 CAM-B3LYP/aug-cc-pVDZ `Polar CPHF=RdFreq`，频率 `1064nm` |
| `run.slurm` | Slurm 作业脚本 |

## 运行步骤

```bash
sbatch run.slurm
# 或
g16 < pna_polar_hyperpolarizability.gjf > pna_polar_hyperpolarizability.log
```

## 预期结果 / 如何检查

```bash
grep -c "Normal termination" pna_polar_hyperpolarizability.log     # 应为 3
grep "NImag" pna_polar_hyperpolarizability.log
grep -A1 "Dipole moment (field-independent basis, Debye)" pna_polar_hyperpolarizability.log | tail -2
grep -n "Alpha(" pna_polar_hyperpolarizability.log      # 静态与 1064 nm 的极化率块
grep -n "Beta("  pna_polar_hyperpolarizability.log      # 静态与 1064 nm 的超极化率块
grep "Isotropic polarizability" pna_polar_hyperpolarizability.log
```

- 静态块标签为 `Alpha(0;0)`、`Beta(0;0,0)`，动态块的标签中带有频率 w，具体包括哪些 β（如电光效应、二次谐波）以输出为准；
- β 主要沿分子长轴（氨基—硝基方向）；$\beta_{\rm tot}=\sqrt{\beta_x^2+\beta_y^2+\beta_z^2}$，$\beta_i=\beta_{iii}+\sum_{j\neq i}\beta_{ijj}$（静态）；
- 单位：1 a.u.(α) ≈ 0.14818 Å³；1 a.u.(β) ≈ 8.6392×10⁻³³ esu；
- 输出注明 “Values do not include the 1/n! factor of 1/2”，与文献比较前确认约定（可能相差 2 倍）。

## 注意事项

- 超极化率对弥散函数非常敏感，不要用 6-31G(d) 之类无弥散的基组得出定量结论。
- 常规杂化泛函（如 B3LYP）往往高估推-拉分子的 β，可以用本例比较 B3LYP、CAM-B3LYP 和 ωB97X-D。
- 频率行 `1064nm` 必须位于电荷多重度行之后的空行之后，并以空行结束。
- 溶液中的 NLO 性质可加 `SCRF`，但溶剂对 β 的影响较大，需要说明所用模型。

返回：[第11章 NMR 化学位移与其他分子性质](../../tutorials/11-nmr-props.md)
