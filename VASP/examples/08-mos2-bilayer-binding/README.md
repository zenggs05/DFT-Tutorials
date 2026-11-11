# 算例 08-mos2-bilayer-binding：双层 MoS2 层间结合能曲线（刚性层扫描）

## 目的

- 构建 2H（AA'）堆垛的双层 MoS2，固定单层几何，只改变层间距 d（上下两层 Mo 原子的竖直距离），扫描结合能曲线；
- 计算层间结合能 $E_b(d) = [E_{\rm bilayer}(d) - 2E_{\rm mono}]/A$（单位 meV/Å²，负值表示结合），找到平衡层间距；
- 体会不同 vdW 方法（`IVDW = 11/12/13`、vdW-DF）对层间距和结合能的影响。

对应教程：[第8章](../../tutorials/08-vdw-2d.md) §8。

## 文件说明

| 文件 | 说明 |
|---|---|
| `make_bilayer.py` | 生成双层（`--d`）或单层参考（`--mono`）POSCAR；参数 `--a`、`--t`（单层 S–S 高度）、`--c` |
| `INCAR` | PBE-D3(BJ) 单点能（`IVDW = 12`，`ENCUT = 450`） |
| `KPOINTS` | Γ 中心 12×12×1，单层与双层共用 |
| `run.sh` | `prepare`：建立 `mono/` 与 `d_5.6 … d_10.0/` 并提交；`collect`：汇总 |
| `analyze_binding.py` | 计算 E_b(d)，在最小值附近三次多项式拟合得到 d_eq，输出 `binding_curve.png` |
| `job.slurm` | 最简 Slurm 脚本 |

POTCAR：

```bash
cat $VASP_PP_PATH/potpaw_PBE/{Mo_pv,S}/POTCAR > POTCAR
```

## 运行步骤

先用 [08-mos2-monolayer](../08-mos2-monolayer/README.md) 优化单层，读出 a 和 S–S 高度 t，再：

```bash
export A=3.19 T=3.13 C=30.0     # replace by your relaxed monolayer values
bash run.sh prepare             # 1 monolayer + 12 bilayer single points
# wait ...
bash run.sh collect
```

## 预期结果 / 如何检查

- 结合能曲线在 d ≈ 6.1–6.3 Å 附近出现极小值（体相 2H-MoS2 的实验层间距 c/2 约 6.15 Å），远处趋于 0；
- 极小值深度在 −20 meV/Å² 量级（Björkman et al., PRL 108, 235502 (2012) 指出多数层状材料的层间结合能约为 20 meV/Å²；D3/vdW-DF 结果会有百分之几十以内的差别，仅作量级参考）；
- 纯 PBE（删掉 `IVDW`）几乎没有结合：曲线很浅或单调——可以作为对照组。

示意输出（数值仅供参考）：

```text
 d(Mo-Mo) A   E_b meV/A^2  E_b meV/f.u.
       6.00        -1x.xx        -1xx.x
       6.20        -2x.xx        -1xx.x
cubic fit around the minimum: d_eq ~ 6.2 A, E_b ~ -2x meV/A^2
```

## 注意事项

- 这里是"刚性层"扫描：层内几何不弛豫。得到 d_eq 后，建议在该距离附近做一次完整弛豫（`ISIF = 2`，或 `ISIF = 3` + `LATTICE_CONSTRAINTS`）确认。
- 单层参考必须放在**同一个晶胞**（相同 a、c）中，使用相同的 ENCUT、k 网格、`IVDW`、smearing，才能让误差相互抵消。
- 层间结合能（binding）与剥离能（exfoliation，从体相剥离一层）定义不同，比较文献数值时要先确认定义与单位（meV/Å²、J/m²、meV/atom）。1 meV/Å² ≈ 0.016 J/m²。
- 想比较不同堆垛（AA、AB、AA'……）时，只需修改 `make_bilayer.py` 中上层原子的面内位置。

[返回第8章](../../tutorials/08-vdw-2d.md)
