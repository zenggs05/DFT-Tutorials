# 示例 14-2：石墨烯 pz 双带 MLWF（解纠缠与紧束缚参数）

配套教程：[第14章 Wannier 函数：MLWF 与 CWF](../../tutorials/14-wannier.md)

## 目的

为单层石墨烯的每个 C 原子构造一个 pz 型 MLWF，得到描述 π/π* 能带与 Dirac 锥的最小紧束缚模型，并从 `.HWR` 中读出在位能和最近邻跳跃积分 $t_1$。π* 带与其他态纠缠，所以需要外窗口 + 内（冻结）窗口的解纠缠。

## 文件说明

| 文件 | 说明 |
|---|---|
| `Graphene_pz_MLWF.dat` | SCF（15×15×1）+ 2 个 pz MLWF（`Wannier.Kgrid 12 12 1`）+ 插值能带 + cube 输出 |
| `hwr_hopping.py` | 读取 `graphene_pz.HWR`，打印在位能（相对 E_F）和最大的若干跳跃积分（eV），仅用标准库 |
| `job.slurm` | Slurm 作业脚本：运行并调用分析脚本 |

能带比较和中心/spread 提取复用示例 14-1 的脚本：`../14-si-mlwf-sp3/band_compare.py`、`../14-si-mlwf-sp3/wf_centers.py`。

结构：$a=2.46$ Å，120° 约定的六角晶格，$c=15$ Å；两个 C 位于 (1/3, 2/3, 1/2) 与 (2/3, 1/3, 1/2)。这个约定下 K = (1/3, 1/3, 0)，M = (1/2, 0, 0)。基组 `C6.0-s2p2d1` + `C_PBE19`；投影子物种 `projC` 为 `C6.0-s1p1`。

## 运行步骤

```bash
cd OpenMX/examples/14-graphene-pz-mlwf
sed -i "s#^DATA.PATH .*#DATA.PATH                   $HOME/openmx4.0/DFT_DATA19#" Graphene_pz_MLWF.dat
export OMP_NUM_THREADS=1
mpirun -np 8 openmx Graphene_pz_MLWF.dat -nt 1 > graphene_pz.std 2>&1
# or
sbatch job.slurm
```

建议先用同一输入（把 `Wannier.Func.Calc` 改为 `off`）只算能带，确认 σ 带顶（Γ 点）和 Γ 点 π* 的能量，再决定内外窗口。

## 预期结果 / 如何检查

```bash
grep -A 8 "Selected bands within the Outer and Inner Windows" graphene_pz.std
grep "CONVERGENCE ACHIEVED" graphene_pz.std | tail -n 2
python3 ../14-si-mlwf-sp3/wf_centers.py graphene_pz.std
python3 hwr_hopping.py graphene_pz.HWR --top 12
python3 ../14-si-mlwf-sp3/band_compare.py graphene_pz.Band graphene_pz.Wannier_Band \
        --emin -2.0 --emax 1.5 --plot graphene_pz_bands.png --ylim -10 10
```

判断标准（定性；仓库不附参考输出）：

- 每个 k 点外窗口内至少 2 条带，内窗口内至多 2 条带；
- 两个 WF 的中心在 C 原子上，spread 相等；
- 最大的跳跃积分是最近邻 $t_1$（R = (0,0,0) 的 $H_{12}$ 以及两个等价的近邻晶胞项），为负值，量级 −2.5 ~ −3 eV（紧束缚文献常用 $|t|\approx$ 2.7–2.8 eV，见 Castro Neto 等, Rev. Mod. Phys. 81, 109 (2009)）；
- 插值能带在 K 点形成 Dirac 锥，内窗口内与原始能带重合。

## 注意事项

- 输入中的窗口是根据石墨烯能带的一般特征给出的起始值，不是计算结果；换基组或泛函后要重新确认。
- 外窗口上限必须足够高，把 Γ 附近的 π* 包进来，否则 Γ 点窗口内没有 pz 成分的态，初始投影接近奇异。
- 二维体系：`Wannier.Kgrid` 的第三个数为 1，`Wannier.Function.Plot.SuperCells` 的第三个数为 0。
- CWF 的 AO 模式只能整壳层选择（`1p` = px、py、pz），所以这种只取 pz 的模型用 MLWF 投影子更方便。
- `.HWR` 中每个 R 块的整数是 Wigner–Seitz 简并度，`hwr_hopping.py` 原样打印；R = (0,0,0) 的简并度为 1。

[返回第14章](../../tutorials/14-wannier.md)
