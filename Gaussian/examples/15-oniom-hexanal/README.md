# 示例 15-1：二层 ONIOM——己醛的 ONIOM(B3LYP:UFF)、ONIOM(B3LYP:PM6) 与全 QM 参考

配套教程：[第15章 ONIOM 多层方法与 QM/MM](../../tutorials/15-oniom.md)（§3、§4、§5）

## 目的

- 写出二层 ONIOM 的原子行：冻结标记、层标记（`H`/`L`）和链接原子（`L H 7`）；
- 理解电荷多重度行 `0 1 0 1 0 1` 的三对含义；
- 读懂 `ONIOM: gridpoint` 与 `ONIOM: extrapolated energy` 输出；
- 对比 MO:MM（UFF）、MO:MO（PM6）两种低层与全 QM 参考，评估分层误差。

## 文件说明

| 文件 | route | 说明 |
|---|---|---|
| `oniom2_b3lyp_uff_organic.gjf` | `ONIOM(B3LYP/6-31G(d):UFF) Opt Freq` | 醛基 + α、β-CH2 为高层，丙基为低层（UFF） |
| `oniom2_b3lyp_pm6_momo.gjf` | `ONIOM(B3LYP/def2SVP:PM6) Opt Freq` | 相同分层，低层换成 PM6（MO:MO） |
| `hexanal_full_qm_reference.gjf` | `B3LYP/6-31G(d) Opt Freq` | 全分子 DFT 参考，原子顺序与 ONIOM 输入相同 |
| `run_g16_slurm.sh` | — | Slurm 作业脚本（8 核、10 GB） |

分层：原子 1–9（O、醛基 C 和 H、两个 CH2）为高层；原子 10–19 为低层。原子 10（C4）的行末写 `L H 7`：它与高层原子 7（C3）成键，在模型体系中被 H 取代。模型体系因此是丙醛 CH3CH2CHO。起始结构为 RDKit/MMFF 预优化的全反式构象。

## 运行步骤

```bash
cd Gaussian/examples/15-oniom-hexanal
for f in oniom2_b3lyp_uff_organic oniom2_b3lyp_pm6_momo hexanal_full_qm_reference; do
    g16 < $f.gjf > $f.log
done
# or on Slurm
for f in *.gjf; do sbatch run_g16_slurm.sh $f; done
```

三个作业都很小，8 核上一般几分钟到二十分钟。

## 预期结果 / 如何检查

```bash
grep -c "Normal termination" *.log                                   # 1 each
grep "ONIOM: gridpoint" oniom2_b3lyp_uff_organic.log | tail -3
grep "ONIOM: extrapolated energy" oniom2_b3lyp_uff_organic.log | tail -1
grep "imaginary frequencies" *.log                                    # no output = no imaginary modes
grep -c "Optimization completed" *.log
```

| 检查项 | 预期 |
|---|---|
| gridpoint 1 / 3（UFF） | 数值很小的 MM 能量（相对力场零点），单独没有绝对意义 |
| gridpoint 2（B3LYP，丙醛） | 普通 DFT 能量，与该步 `SCF Done` 一致 |
| `extrapolated energy` | $E_3 + E_2 - E_1$，即 ONIOM 总能量，比较能量时用它 |
| MO:MO 作业 | 三个 gridpoint 都是量子化学能量 |
| 频率 | 三个作业都应无虚频 |

几何对比：从三个 log 中取出最终结构（[../../scripts/log2xyz.py](../../scripts/log2xyz.py)，`python ../../scripts/log2xyz.py -h` 查看用法），比较 C=O、C1–C2 键长和 O=C1–C2–C3 二面角。醛基附近的几何差别应很小；差别越大，说明分层（或低层方法）越不合适。注意 PM6 作业的高层用的是 def2-SVP，与另两个作业的 6-31G(d) 不同，比较时要考虑这一点，也可以自行把它改成 6-31G(d) 再比较。

## 注意事项

- 不能用 `SCF Done` 比较 ONIOM 结构之间的能量，必须用 `ONIOM: extrapolated energy`。
- 链接原子写在**低层**原子的行上，`bonded-to` 填与它相连的**高层**原子序号；增删原子后要同步修改。
- UFF 会自动分配原子类型，这里没有写类型；改用 Amber 时所有原子都必须显式写类型（见示例 [15-qmmm-phenol-water](../15-qmmm-phenol-water/README.md)）。
- 想看三个子计算的独立输入，在 route 中把 `ONIOM(...)` 改成 `ONIOM(...)=OnlyInputFiles`，Gaussian 会只生成输入不计算。

返回教程：[第15章](../../tutorials/15-oniom.md)
