# 示例 12-hirshfeld-cm5-acetate：醋酸根的 Hirshfeld/CM5、MK 与 NPA 电荷对比

配套教程：[第12章 布居分析、NBO 与静电势/RESP 电荷](../../tutorials/12-population-charges.md)（第 8 节）。

## 目的

在 SMD(水) 中优化醋酸根阴离子（CH3COO⁻），并在**同一个溶液相密度**上比较 Hirshfeld、CM5、MK 和 NPA 电荷。演示官方推荐的做法：追加 `Geom=AllCheck Density=Check` 作业步，读取已保存的密度做不同的布居分析，而不重复 SCF。

## 文件说明

| 文件 | 内容 |
|---|---|
| `acetate_hirshfeld_cm5.gjf` | 第 1 步 B3LYP-D3(BJ)/6-311+G(d,p) SMD(水) Opt+Freq，`Pop=Hirshfeld`；第 2 步 `Density=Check Pop=MK`；第 3 步 `Density=Check Pop=NPA` |
| `run.slurm` | Slurm 作业脚本 |

原子编号：C1（甲基碳）、C2（羧基碳）、O3、O4、H5–H7。电荷多重度为 `-1 1`。

## 运行步骤

```bash
sbatch run.slurm
# 或
g16 < acetate_hirshfeld_cm5.gjf > acetate_hirshfeld_cm5.log
```

## 预期结果 / 如何检查

```bash
grep -c "Normal termination" acetate_hirshfeld_cm5.log       # 应为 4（Opt、Freq、MK、NPA）
grep "NImag" acetate_hirshfeld_cm5.log                         # NImag=0
# Hirshfeld (Q-H) 与 CM5 (Q-CM5)：取最后一次输出
grep -A10 "Hirshfeld charges, spin densities, dipoles, and CM5 charges" acetate_hirshfeld_cm5.log | tail -11
# MK 电荷
grep -A9 "ESP charges:" acetate_hirshfeld_cm5.log | tail -9
# NPA 电荷
grep -A12 "Summary of Natural Population Analysis" acetate_hirshfeld_cm5.log | tail -9
```

把四种电荷整理成表（每种方法的总电荷都应为 −1）。定性上：

- 两个 O 在所有方法中都应近似等价（羧酸根的共振）；
- Hirshfeld 电荷绝对值最小，CM5 与 MK 的量级接近，NPA 的绝对值通常最大；
- 甲基 H 的 MK 电荷可能因格点离得远而不太稳定，这是 ESP 拟合电荷的常见现象。

## 注意事项

- 第 2、3 步使用 `Geom=AllCheck`，所以没有标题行和分子说明，route 后只有一个空行。
- `Density=Check` 隐含 `Guess=Only ChkBasis`：不会重新做 SCF，基组也从 chk 读取，因此 route 中不写基组。
- Opt+Freq 中布居分析会出现多次（优化首步、末步和 Freq 步），务必读取最后一次。
- 想把 CM5 电荷存为后续 MM/ONIOM 计算使用的电荷，可在相应步骤中写 `Pop=(Hirshfeld,SaveCM5)`（C.01 起），参见第 15 章。

返回：[第12章 布居分析、NBO 与静电势/RESP 电荷](../../tutorials/12-population-charges.md)
