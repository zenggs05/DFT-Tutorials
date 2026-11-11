# 示例 15-3：三层 ONIOM(B3LYP/6-31G(d):PM6:UFF)——4-叔丁基环己酮

配套教程：[第15章 ONIOM 多层方法与 QM/MM](../../tutorials/15-oniom.md)（§3.4、§7）

## 目的

- 写出三层 ONIOM 的层标记 `H`/`M`/`L` 和两种边界上的链接原子；
- 理解三层 ONIOM 的 6 对电荷多重度；
- 认识真实体系、中间体系、模型体系分别是什么分子。

## 文件说明

| 文件 | 说明 |
|---|---|
| `oniom3_high_medium_low.gjf` | `ONIOM(B3LYP/6-31G(d):PM6:UFF) Opt`，29 个原子 |
| `run_g16_slurm.sh` | Slurm 作业脚本（8 核、10 GB） |

分层（起始结构为 RDKit/MMFF 预优化的椅式构象，叔丁基在平伏键）：

| 层 | 原子序号 | 内容 | 链接原子 |
|---|---|---|---|
| High | 1–8 | C=O、两个 α-CH2 | — |
| Medium | 9–16 | 两个 β-CH2、γ-CH | 原子 9 `M H 3`、原子 12 `M H 6`（连到 α-碳） |
| Low | 17–29 | 叔丁基 | 原子 17 `L H 15`（季碳连到 γ-碳） |

加上链接原子后：模型体系（高层计算）是丙酮，中间体系（中层计算）是环己酮，真实体系是 4-叔丁基环己酮。电荷多重度行写 6 对 `0 1`，顺序为 real-low、int-medium、int-low、model-high、model-medium、model-low。

## 运行步骤

```bash
cd Gaussian/examples/15-oniom3-cyclohexanone
g16 < oniom3_high_medium_low.gjf > oniom3_high_medium_low.log
# or
sbatch run_g16_slurm.sh oniom3_high_medium_low.gjf
```

体系很小，几分钟到二十分钟即可完成。想先看各子计算的输入，可把 route 改为 `ONIOM(B3LYP/6-31G(d):PM6:UFF)=OnlyInputFiles` 单独运行一次（只生成输入，不计算）。

## 预期结果 / 如何检查

```bash
grep -c "Normal termination" oniom3_high_medium_low.log        # expect 1
grep "ONIOM: gridpoint" oniom3_high_medium_low.log | tail -5
grep "ONIOM: extrapolated energy" oniom3_high_medium_low.log | tail -1
grep -c "Optimization completed" oniom3_high_medium_low.log    # expect 1
```

| 检查项 | 预期 |
|---|---|
| gridpoint 行 | 每个几何步有 5 个子计算：低层算真实体系和中间体系、中层算中间体系和模型体系、高层算模型体系 |
| 外推能量 | $E^{\rm low}(R)+E^{\rm medium}(IM)+E^{\rm high}(SM)-E^{\rm low}(IM)-E^{\rm medium}(SM)$ |
| 优化结构 | 保持椅式、叔丁基平伏；C=O 键长与全 QM 计算接近 |

## 注意事项

- 三层 ONIOM 的电荷多重度只写一对时，所有子计算都用这一对；本例写全 6 对是为了演示顺序。带电或开壳层体系必须按顺序为每个子体系写对。
- 中层原子 9 的 `M H 3` 表示：在模型体系（高层和中层都会计算它）中，原子 9 被一个连到原子 3 的 H 取代。低层原子 17 的 `L H 15` 表示：在中间体系（中层和低层计算它）中，原子 17 被一个连到原子 15 的 H 取代；模型体系不含原子 15，所以这个链接 H 不会出现在模型体系里。
- 练习：给 route 加 `Freq` 计算频率，或把高层换成 `M062X/def2SVP` 观察外推能量的变化；比较时保持分层不变。

返回教程：[第15章](../../tutorials/15-oniom.md)
