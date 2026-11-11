# 示例 14-1：水二聚体的 Counterpoise（BSSE）校正与基组收敛

配套教程：[第14章 弱相互作用与 BSSE（Counterpoise）](../../tutorials/14-weak-bsse.md)（§3、§4、§7）

## 目的

- 写出正确的 `Counterpoise=2` 输入：`(Fragment=n)` 标记与 `0 1 0 1 0 1` 电荷多重度行；
- 读懂 `Counterpoise corrected energy`、`BSSE energy`、`sum of fragments` 和两行 `complexation energy`；
- 比较 def2-SVP、def2-TZVP、aug-cc-pVTZ 三个基组下 BSSE 的大小；
- 用自由单体能量求形变能与结合能。

## 文件说明

| 文件 | 说明 |
|---|---|
| `water_dimer_counterpoise.gjf` | 4 步作业：第 1–3 步为 B3LYP-D3(BJ) 下三个基组的 CP 单点（同一结构）；第 4 步为 def2-TZVP 下自由水分子的 Opt Freq |
| `cp_summary.py` | 逐个 Link1 步骤提取 CP 汇总块，换算为 kcal/mol；`--mono` 给出自由单体能量后计算形变能与结合能（仅用 Python 标准库） |
| `run_g16_slurm.sh` | Slurm 作业脚本（8 核、10 GB） |

二聚体结构为手工搭建的 Cs 型氢键构型：O···O 2.91 Å，O–H 0.96 Å，H–O–H 104.5°，O–H···O 近似线性，受体分子平面与供体平面垂直。它是合理的近似结构，但不是文献或基准集坐标。碎片 1 为氢键供体（原子 1–3），碎片 2 为受体（原子 4–6）。

## 运行步骤

```bash
cd Gaussian/examples/14-cp-water-dimer
g16 < water_dimer_counterpoise.gjf > water_dimer_counterpoise.log
# or
sbatch run_g16_slurm.sh water_dimer_counterpoise.gjf
```

整个作业在 8 核上通常只需几分钟（aug-cc-pVTZ 一步最慢）。

汇总结果：

```bash
python cp_summary.py water_dimer_counterpoise.log
EMONO=$(grep "SCF Done" water_dimer_counterpoise.log | tail -1 | awk '{print $5}')
python cp_summary.py water_dimer_counterpoise.log --mono $EMONO $EMONO
```

第二条命令把第 4 步（def2-TZVP 自由水）的能量作为两个单体的参考，所以只有第 2 步（同为 def2-TZVP）给出的形变能和结合能有意义。

## 预期结果 / 如何检查

```bash
grep -c "Normal termination" water_dimer_counterpoise.log            # expect 4
grep -A4 "Counterpoise corrected energy" water_dimer_counterpoise.log
grep "Counterpoise: doing" water_dimer_counterpoise.log | head -4     # order of the fragment calculations
grep -c "SCF Done" water_dimer_counterpoise.log                        # 5 per CP step + monomer opt steps
```

定性预期：

| 检查项 | 预期 |
|---|---|
| BSSE | 三步都为正值，def2-SVP 最大，aug-cc-pVTZ 最小 |
| ΔE_int(raw) 与 ΔE_int(CP) | 小基组差别明显，大基组下两者接近 |
| ΔE_int(CP) | 负值，量级为几 kcal/mol；高精度参考值约 −5.0 kcal/mol（CCSD(T)/CBS，S22 基准集，对应其参考几何；近似值） |
| 形变能（第 2 步） | 很小（本例结构未优化，单体键长键角与优化值略有差别） |

`cp_summary.py` 算出的 raw/corrected 相互作用能应与 Gaussian 输出的两行 `complexation energy` 一致，可借此确认自己理解了各项的定义。

## 注意事项

- CP 作业不产生分子轨道，所以 CP 步骤没有写 `%chk`；需要波函数（IGMH 等）时另做普通单点。
- `Counterpoise` 不能与 `SCRF`、`ONIOM` 同时使用。
- 教程 §7 给出了用 `O-Bq`/`H-Bq` 鬼原子手动计算 $E_A^{\,AB}(AB)$ 的写法，可以与本例的 `BSSE energy` 互相核对。
- 正式研究应先优化二聚体结构（普通优化或 `Opt Counterpoise=2`），再把优化坐标连同 `(Fragment=n)` 完整写入 CP 单点输入。

返回教程：[第14章](../../tutorials/14-weak-bsse.md)
