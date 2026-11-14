# 算例 12-gaas-soc-band：GaAs 能带的自旋轨道劈裂（Γ 与 L 点）

## 目的

- 练习 SOC 计算的最小设置：`scf.SpinPolarization nc` + `scf.SpinOrbit.Coupling on`，DFT_DATA19 的 j 依赖赝势可直接使用；
- 用同一组 j 依赖赝势、只切换 `scf.SpinOrbit.Coupling on/off`，比较全相对论与标量相对论能带；
- 定量读出价带顶在 Γ 点的 SOC 劈裂（Γ8–Γ7，split-off 能量 Δ0）和 L 点的劈裂，与手册和实验对照。

对应教程：[第12章 自旋轨道耦合与拓扑性质](../../tutorials/12-soc-topology.md) §3。

## 文件说明

| 文件 | 说明 |
|---|---|
| `GaAs_SOC_band.dat` | 闪锌矿 GaAs 原胞（a = 5.653 Å，实验值），Ga7.0-s3p2d2 / As7.0-s3p2d2（Standard），GGA-PBE，nc + SOC，k 网格 8³，能带路径 L–Γ–X–U\|K–Γ |
| `GaAs_noSOC_band.dat` | 同上，仅 `scf.SpinOrbit.Coupling off`（OpenMX 自动按 j 简并度平均赝势 = 标量相对论） |
| `soc_splitting.py` | 读取 `.Band`，在指定标签的 k 点按简并度分组列出 E_F 附近能级，并给出最高两组价带、最低两组导带之间的劈裂（numpy） |
| `job.slurm` | 依次运行两个计算并调用 `soc_splitting.py`；若 PATH 中有 `bandgnu13` 则顺便生成 gnuplot 数据 |

运行前修改两个 `.dat` 中的 `DATA.PATH`。

## 运行步骤

```bash
sbatch job.slurm
# or by hand:
mpirun -np 8 openmx GaAs_SOC_band.dat   -nt 1 > gaas_soc.std   2>&1
mpirun -np 8 openmx GaAs_noSOC_band.dat -nt 1 > gaas_nosoc.std 2>&1
python soc_splitting.py gaas_soc.Band   --at G L
python soc_splitting.py gaas_nosoc.Band --at G L
```

画能带可用 `bandgnu13 gaas_soc.Band` 后 `gnuplot gaas_soc.GNUBAND`，或仓库脚本 [`omx_band_plot.py`](../../scripts/omx_band_plot.py)（`-h` 查看用法）。

## 预期结果 / 如何检查

- **无 SOC**：Γ 点价带顶是 6 重简并（非共线计算中每个空间态带两个自旋分量，Γ15 的 3 重轨道简并 × 2）；
- **有 SOC**：价带顶分裂为 4 重简并的 Γ8 与其下方 2 重简并的 Γ7。`soc_splitting.py` 在 G 处打印的 "splitting of the two highest valence groups" 就是 Δ0。参考值：实验约 0.34 eV；OpenMX 4.0 手册表 6（LDA、s2p2d2 基组）给出 Γ15v 处 0.344 eV、L3v 处 0.213 eV；
- 在 L 点，最高价带（无 SOC 时 4 重）分裂为两组 2 重简并能级，劈裂约为 Γ 点的 2/3 量级；
- SOC 使 Γ8 上移约 Δ0/3，因此带隙比无 SOC 时小约 0.1 eV；GGA 本身严重低估 GaAs 带隙（实验 1.52 eV），不要用来比较绝对带隙；
- 闪锌矿结构没有反演中心，SOC 计算中远离 Γ 的一般 k 点上（如 Γ–K 方向）能带会出现很小的自旋劈裂（Dresselhaus 型），而 Γ–L 等高对称线上可能仍简并。

具体数值依赖泛函、基组与 k 点，以你的计算为准。

## 注意事项

- 能带路径坐标以原胞倒格矢为基（本例原胞基矢顺序为 a/2(1,1,0)、a/2(1,0,1)、a/2(0,1,1)）；若改了基矢顺序，必须重新计算 L、X、U、K 的分数坐标，或改用 `<Band.KPath.UnitCell`。
- SOC 只能在非共线框架下打开；共线计算中打开 SOC 会报 `Spin-orbit coupling is not supported for collinear DFT calculations.`
- 两个计算除 `scf.SpinOrbit.Coupling` 外参数完全相同，劈裂才可比较。
- patch 4.0.1 修正了官方 work 目录中的 `GaAs.dat`（使其能复现手册 38.1 节的结果）；用官方文件对照时请先打补丁。

[返回第12章](../../tutorials/12-soc-topology.md)
