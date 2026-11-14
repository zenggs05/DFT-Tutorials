# 示例 17-1：碳原子链的 NEGF 输运（三步法、Step 3 重算、I–V 曲线）

配套教程：[第17章 NEGF 电子输运](../../tutorials/17-negf-transport.md)

## 目的

用一维单原子碳链（C–C 1.5 Å）完整走一遍 OpenMX 的 NEGF 三步法：

1. **Step 1**：用电极（lead）晶胞做普通能带计算，生成 `lead-chain.hks`；
2. **Step 2**：对 L0|C0|R0 扩展中心区做自洽 NEGF 计算，生成 `negf-chain.tranb`；
3. **Step 3**：Step 2 结束后自动计算透射谱、电流、电导、本征通道和电流密度。

在此基础上，再演示只重算 Step 3（`NEGF.tran.SCF.skip on`）、有限偏压计算，以及用脚本扫描偏压得到 I–V 曲线。理想无缺陷链的透射谱在费米能级附近应接近整数 2（两个简并的 π 通道），这是检验整套设置是否正确的最简单判据。

## 文件说明

| 文件 | System.Name | 说明 |
|---|---|---|
| `Lead-Chain.dat` | `lead-chain` | Step 1：4 原子电极胞（6.0 Å），`scf.Kgrid 160 1 1`，`NEGF.output_hks on`，同时输出电极能带 |
| `NEGF-Chain.dat` | `negf-chain` | Step 2+3，0 V：C0 为 8 个原子，L0/R0 各 4 个原子 |
| `NEGF-Chain_step3.dat` | `negf-chain` | 只重算 Step 3：读取 `negf-chain.tranb`，用更细的能量网格（−2…2 eV，2001 点）并在 3 个能量处算本征通道 |
| `NEGF-Chain_bias0.5.dat` | `negf-chain-0.5` | Step 2+3，偏压 0.5 V |
| `run_iv.sh` | `negf-chain-V*` | 扫描偏压 0–1.0 V，汇总电流到 `iv.csv`；`RESTART=1` 时以前一偏压的结果为初猜 |
| `plot_tran.py` | — | k 平均透射谱作图，并计算零偏压 Landauer 电导和 Seebeck 系数 |
| `job_slurm.sh` | — | 按顺序运行 Step 1 → 0 V → Step 3 重算 → 0.5 V，遇到失败即停止 |

基组为 `C6.0-s2p2d1` + `C_PBE19`（GGA-PBE），电子温度 300 K，截断能 220 Ry。Step 1 和 Step 2 的这几项设置保持一致。

**为什么电极胞取 4 个原子（6 Å）而不是 3 个（4.5 Å）？** NEGF 要求电极第 i 个胞的基函数只与相邻的 i±1 个胞重叠。C6.0 的截断半径是 3.175 Å，两个原子的轨道在 6.35 Å 以内会有重叠。4.5 Å 的胞里，第 i 与第 i+2 个胞最近的原子相距 2×4.5−3.0=6.0 Å，小于 6.35 Å，不满足条件；6 Å 胞对应的距离是 7.5 Å，满足。（官方示例用的是更短程的 C5.0，所以 4.5 Å 也够。）

## 运行步骤

```bash
cd OpenMX/examples/17-negf-carbon-chain
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" *.dat

# Step 1 -> Step 2+3 (0 V) -> Step 3 only -> 0.5 V
mpirun -np 8 openmx Lead-Chain.dat -nt 1 > lead-chain.std 2>&1
mpirun -np 8 openmx NEGF-Chain.dat -nt 1 > negf-chain.std 2>&1
mpirun -np 8 openmx NEGF-Chain_step3.dat -nt 1 > negf-chain-step3.std 2>&1
mpirun -np 8 openmx NEGF-Chain_bias0.5.dat -nt 1 > negf-chain-0.5.std 2>&1
#   or on a cluster: sbatch job_slurm.sh

# I-V curve (0.0 ... 1.0 V)
NP=8 bash run_iv.sh                 # RESTART=1 NP=8 bash run_iv.sh  (reuse previous bias)

# transmission plot, conductance and Seebeck coefficient at 300 K
python3 plot_tran.py negf-chain.tran0_0 --mu 0.0 --plot negf-chain_T.png
```

## 预期结果 / 如何检查

```bash
grep "normally finished" *.std
grep -A8 "<TRAN_Calc_GridBound>" negf-chain.std | tail -n 6     # region labels 12/1/13 of every atom
grep "Intrinsic chemical potential" -A2 negf-chain.std
grep "Chemical potential of left lead\|Chemical potential of right lead" negf-chain-0.5.out
grep "average conductance" negf-chain.conductance
grep "average current" negf-chain-0.5.current
cat iv.csv
```

- `lead-chain.Band`（用 `bandgnu13` 或 [`omx_band_plot.py`](../../scripts/omx_band_plot.py) 画图）：两条简并的 π 带穿过化学势，链是金属性的。
- `<TRAN_Calc_GridBound>` 的标签：L0 的 4 个原子是 `12`（与 L1 有重叠），C0 的 8 个是 `1`，R0 的 4 个是 `13`。如果出现别的组合，或者 C0 原子被标成 12/13，说明几何不满足 NEGF 的重叠条件，需要加长 L0/R0 或 C0。
- `negf-chain.tran0_0`：第 4 列是能量（eV，以左电极化学势为零），第 6、8 列是自旋向上/向下的透射。E=0 附近 T 应接近 2，并呈现平台/台阶结构，这是理想弹道输运的特征。`.conductance` 中的平均电导（单位 $G_0=e^2/h$，每个自旋通道各列一次）应接近 2 和 2。
- 0.5 V：`.out` 中左右电极的化学势相差 0.5 eV（4.0 中偏压对称地加在两侧，各为 ±V/2）；`negf-chain-0.5.current` 给出以安培为单位的电流。理想链的电流约等于 $(e^2/h)\times(T_\uparrow+T_\downarrow)\times V$，0.5 V 时为几十 μA 量级，这是粗估。
- `iv.csv`：理想链的 I–V 在低偏压下近似线性。

上面的描述都是定性判据，具体数值以你的计算为准。

## 注意事项

- Step 2 与 Step 1 的 `scf.XcType`、`scf.ElectronicTemperature`、`scf.energycutoff`、基组和 b/c 方向的几何必须一致；L0/R0 的原子坐标必须是 Step 1 晶胞的整体平移。
- NEGF 输入**不能**使用 FRAC 坐标；不写 `<Atoms.UnitVectors`，bc 面晶胞由 `.hks` 文件提供。
- 只重算 Step 3 时，`System.Name` 必须和 Step 2 相同，因为程序要读取 `System.Name.tranb`。
- NEGF 的 `scf.criterion` 作用于 `NormRD`，而不是 `dUele`。
- `RESTART=1` 用到的 `scf.restart.filename` 在手册中只在 `c2n` 场景下出现，不过 4.0 源码在读取重启文件时普遍使用这个关键词。重启文件中的密度按 MPI 进程分块存储，所以各偏压之间必须保持相同的进程数。请先在本例上确认可用，再用到大体系。
- 练习：把 C0 中间一根键拉长到 2.0 Å（把原子 5–8 及 R0 整体沿 x 平移 0.5 Å），看透射谱如何从平台变成共振峰。
