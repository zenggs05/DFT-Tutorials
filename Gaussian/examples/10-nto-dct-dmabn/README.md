# 示例 10-nto-dct-dmabn：DMABN 激发态的 NTO、DCT 与空穴-电子分析

配套教程：[第10章 激发态 TD-DFT](../../tutorials/10-tddft.md)（第 9 节）。

## 目的

4-(二甲氨基)苯甲腈（DMABN）是研究电荷转移激发的经典推-拉分子。本例在一个 Link1 输入里完成：

1. 基态优化与频率（SMD 乙腈）；
2. TD-DFT 计算 6 个单重态，并用 `IOp(9/40=4)` 打印足够多的组态系数，供 Multiwfn 使用；
3. 分别生成 S1、S2 的自然跃迁轨道（NTO），各存一个 chk；
4. 对 S2 计算弛豫密度并给出 `Pop=DCT` 电荷转移指标。

## 文件说明

| 文件 | 内容 |
|---|---|
| `nto_analysis_link1.gjf` | 5 个 Link1 步骤：Opt+Freq → TD → NTO(S1) → NTO(S2) → DCT(S2) |
| `run.slurm` | Slurm 作业脚本，结束后对所有 chk 执行 formchk |

生成的 chk：`dmabn_opt.chk`、`dmabn_td.chk`（基态轨道 + TD 结果）、`dmabn_nto_s1.chk`、`dmabn_nto_s2.chk`（轨道已替换为 NTO）、`dmabn_dct_s2.chk`。

## 运行步骤

```bash
sbatch run.slurm
# 或
g16 < nto_analysis_link1.gjf > nto_analysis_link1.log
for c in dmabn_td dmabn_nto_s1 dmabn_nto_s2; do formchk $c.chk $c.fchk; done
```

## 预期结果 / 如何检查

```bash
grep -c "Normal termination" nto_analysis_link1.log      # 应为 6（Opt、Freq、TD、NTO、NTO、DCT）
grep "Excited State" nto_analysis_link1.log | head -6 | awk '{print $3, $5, $7, $9}'
grep -i -A12 "DCT" nto_analysis_link1.log | head -40     # 先定位 DCT 输出块的格式
```

NTO 可视化：

```bash
cubegen 0 MO=HOMO dmabn_nto_s2.fchk nto_s2_hole.cub 0 h
cubegen 0 MO=LUMO dmabn_nto_s2.fchk nto_s2_particle.cub 0 h
```

也可以直接在 GaussView 中打开 `dmabn_nto_s2.fchk` 查看轨道；最高“占据”NTO 是空穴，最低“空”NTO 是粒子，其“能量”位置存放的是 NTO 本征值 λ。

Multiwfn 空穴-电子分析：`Multiwfn dmabn_td.fchk` → 主功能 18 → 子功能 1 → 输入含 TD 结果的输出文件路径 → 选择激发态编号，查看 D 指数、Sr、t 指数并导出空穴/电子 cube。本例的 log 中第 2 步和第 5 步都含有 TD 输出，为避免歧义，建议把 TD 那一段（log 中第 3 个 `Entering Link 1` 开始，前两段分别是 Opt 和自动追加的 Freq，到其后第一个 `Normal termination` 为止）另存为 `dmabn_td_only.log` 再交给 Multiwfn，或者把 TD 步骤单独作为一个作业运行。

判断依据：空穴集中在二甲氨基-苯环、粒子集中在苯环-氰基，且 D 指数、DCT 距离明显大于另一个态，即该态更具 CT 特征。

## 注意事项

- NTO 步骤会把 chk 中的正则轨道替换掉，所以一定通过 `%oldchk`（读）+ `%chk`（新文件）写到副本里，不要写回 `dmabn_td.chk`。
- Multiwfn 需要的是 **TD 作业的 fchk**（`dmabn_td.fchk`），不是 NTO 的 fchk。
- 气相与乙腈中 DMABN 低能态的顺序和 CT 程度不同，可以把 `SCRF` 去掉重算比较。
- `Pop=DCT` 的输出格式随版本可能略有差异，写解析脚本前先 grep 确认。

返回：[第10章 激发态 TD-DFT](../../tutorials/10-tddft.md)
