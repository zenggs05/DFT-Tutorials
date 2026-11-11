# 示例 14-2：苯二聚体的 CP 校正优化（Opt Counterpoise=2）

配套教程：[第14章 弱相互作用与 BSSE（Counterpoise）](../../tutorials/14-weak-bsse.md)（§5）

## 目的

- 学会写 `Opt Counterpoise=2` 做 CP 校正的几何优化；
- 对比 CP 校正优化与普通优化得到的 π-π 堆积距离；
- 体会色散校正对 π-π 堆积体系的决定性作用。

## 文件说明

| 文件 | 说明 |
|---|---|
| `benzene_dimer_cp_opt.gjf` | 3 步作业，B3LYP-D3(BJ)/def2-SVP：① `Opt Counterpoise=2 NoSymm`；② 普通 `Opt NoSymm`；③ 自由苯 `Opt Freq` |
| `run_g16_slurm.sh` | Slurm 作业脚本（16 核、20 GB，对应 `%nprocshared=16`、`%mem=16GB`） |

起始结构：两个理想苯环（C–C 1.39 Å，C–H 1.09 Å）平行放置，层间距 3.50 Å，沿 x 方向错位 1.70 Å（平行错位型），是手工给出的起始猜测。碎片 1 为原子 1–12，碎片 2 为原子 13–24。

## 运行步骤

```bash
cd Gaussian/examples/14-cp-benzene-dimer
sbatch run_g16_slurm.sh benzene_dimer_cp_opt.gjf
# or locally
g16 < benzene_dimer_cp_opt.gjf > benzene_dimer_cp_opt.log
```

CP 校正优化每一步要做 5 次 SCF，耗时明显高于普通优化；16 核上整个作业一般需要一到数小时。可以在运行中监控：

```bash
grep "Counterpoise corrected energy" benzene_dimer_cp_opt.log | tail -3
python ../../scripts/opt_monitor.py -h     # plot convergence of each opt step
```

## 预期结果 / 如何检查

```bash
grep -c "Normal termination" benzene_dimer_cp_opt.log        # expect 3
grep -c "Optimization completed" benzene_dimer_cp_opt.log    # expect 3
grep "Counterpoise corrected energy" benzene_dimer_cp_opt.log | tail -1
grep -A4 "Counterpoise corrected energy" benzene_dimer_cp_opt.log | tail -5
```

取出第 1、2 步的最终结构（[../../scripts/log2xyz.py](../../scripts/log2xyz.py)，`python ../../scripts/log2xyz.py -h` 查看用法），在 GaussView 中定义两个苯环的质心并测量层间距和错位距离。

定性预期：

| 检查项 | 预期 |
|---|---|
| 两步优化后的构型 | 仍为平行错位型（若滑向 T 型，说明势能面很平坦，记录并比较能量） |
| 层间距 | CP 校正优化的结果一般大于普通优化 |
| 最后一次 `complexation energy (corrected)` | 负值，约几 kcal/mol；高精度参考值约 −2.6 ~ −2.7 kcal/mol（S22 系列基准，近似值），双 ζ 基组只宜定性比较 |
| 第 3 步 | 自由苯无虚频，用于形变能 |

## 注意事项

- 本例刻意使用 def2-SVP，使 BSSE 对结构的影响足够明显；定量的相互作用能请在优化结构上用 def2-TZVP 或更大基组做 CP 单点（教程 §8）。
- 去掉 `EmpiricalDispersion=GD3BJ` 重新优化（教程练习 2），可以看到 B3LYP 本身几乎不能把两个苯环维系在一起。
- 若想在第 1 步结构上继续算高水平 CP 单点，请把优化坐标连同 `(Fragment=n)` 完整写出，不要依赖 `Geom=Check` 携带碎片信息（官方文档未说明是否保存）。
- 弱络合物势能面平坦，若频率计算出现很小的虚频，可尝试 `Opt=Tight` 和 `Int=SuperFine` 后重新优化（[第6章](../../tutorials/06-freq-thermo.md)）。

返回教程：[第14章](../../tutorials/14-weak-bsse.md)
