# 示例 17-2：7-扶手椅型石墨烯纳米带（7-AGNR）的 NEGF 透射：纯净 vs N 掺杂

配套教程：[第17章 NEGF 电子输运](../../tutorials/17-negf-transport.md)

## 目的

计算氢钝化的 7-AGNR 的零偏压透射谱，比较纯净纳米带和中心区一个 C 被 N 取代后的情况：

- 纯净纳米带：电极与中心区完全相同，$T(E)$ 应是整数台阶，费米能级附近有一段 $T=0$ 的能隙（7-AGNR 是半导体）；
- N 掺杂：取代式 N 是施主，会在导带一侧引入散射，透射台阶上出现凹陷，本征通道可以显示散射发生在哪里。

本例同时示范怎样为二维材料纳米带选取电极胞长度（满足“只与相邻胞重叠”的条件），以及怎样输出 NEGF 的态密度。

## 文件说明

| 文件 | System.Name | 说明 |
|---|---|---|
| `Lead-7AGNR.dat` | `lead-7agnr` | Step 1：电极胞含 2 个周期（8.52 Å，28 C + 8 H），`scf.Kgrid 100 1 1`，输出电极能带 |
| `NEGF-7AGNR.dat` | `negf-7agnr` | Step 2+3：L0、C0、R0 各 2 个周期（各 36 个原子），零偏压，本征通道取在 E=±1.2 eV |
| `NEGF-7AGNR-N.dat` | `negf-7agnr-n` | 同上，但中心区第 12 号原子（x=10.65 Å，y=11.23 Å，内部碳位）换成 N（`N6.0-s2p2` + `N_PBE19`，2.5 2.5） |
| `job_slurm.sh` | — | 依次运行 Step 1 → 纯净 → N 掺杂 |

几何由 ASE `graphene_nanoribbon(3.5, 1, type='armchair', saturated=True)` 生成：C–C 1.42 Å，C–H 1.09 Å，未弛豫。纳米带的周期方向转到 x（输运方向，即 a 轴），宽度方向为 y，法向为 z；b=20 Å、c=12 Å 中留有足够真空。为节省 NEGF 的计算量，基组用 DFT_DATA19 推荐表的 Quick 档（`C6.0-s2p2`、`H5.0-s2`）。正式研究请换成 Standard 档，并做收敛测试。

**电极长度的检查**：C6.0 的截断半径为 3.175 Å。1 个周期（4.26 Å）的胞中，第 i 与第 i+2 个胞最近的原子相距仅 5.68 Å < 6.35 Å，不满足条件。取 2 个周期（8.52 Å）时，这一距离为 9.94 Å，满足。

## 运行步骤

```bash
cd OpenMX/examples/17-negf-7agnr
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" *.dat
mpirun -np 32 openmx Lead-7AGNR.dat  -nt 1 > lead-7agnr.std 2>&1
mpirun -np 32 openmx NEGF-7AGNR.dat  -nt 1 > negf-7agnr.std 2>&1
mpirun -np 32 openmx NEGF-7AGNR-N.dat -nt 1 > negf-7agnr-n.std 2>&1
#   or: sbatch job_slurm.sh

python3 ../17-negf-carbon-chain/plot_tran.py negf-7agnr.tran0_0   --plot T_pristine.png --csv T_pristine.csv
python3 ../17-negf-carbon-chain/plot_tran.py negf-7agnr-n.tran0_0 --plot T_Ndoped.png  --csv T_Ndoped.csv
```

## 预期结果 / 如何检查

```bash
grep "normally finished" *.std
grep -A12 "<TRAN_Calc_GridBound>" negf-7agnr.std | tail -n 16    # 36 x '12', 36 x '1', 36 x '13'
grep "average conductance" negf-7agnr.conductance negf-7agnr-n.conductance
ls negf-7agnr.tranec0_*_r.cube negf-7agnr.Dos.val
```

- `lead-7agnr.Band`：Γ 点是直接带隙，半导体性。半局域 DFT（LDA/GGA）给出的 7-AGNR 带隙在 1–2 eV 范围内（参见 Son, Cohen & Louie, *PRL* **97**, 216803 (2006)，<https://doi.org/10.1103/PhysRevLett.97.216803>；GW 修正后会显著增大，见 Yang *et al.*, *PRL* **99**, 186801 (2007)，<https://doi.org/10.1103/PhysRevLett.99.186801>）。Quick 基组与未弛豫的几何会让数值有偏差，这里只要求定性一致。
- 纯净纳米带：`negf-7agnr.tran0_0` 在带隙内 $T\approx0$，带边以外是 1、2… 的整数台阶。透射台阶的边缘应与电极能带的带边对应（注意 tran 文件的能量以左电极化学势为零点）。
- N 掺杂：导带一侧的台阶出现凹陷或 Fano 型线形；比较 E=+1.2 eV 处两种结构的本征通道 cube 文件（`*.tranec0_1_0_0_r.cube`，文件名中的数字依次为 k、能量、自旋、通道序号），可以看到 N 附近的散射。
- `Dos.fileout on` 生成 `negf-7agnr.Dos.val/.Dos.vec`，用 DosMain 处理的方法与普通计算相同（第9章）；能量零点同样是左电极化学势。

## 注意事项

- 半导体电极的化学势落在带隙中，具体位置受电子温度影响。比较不同计算时，要对齐各自左电极的化学势。
- 中心区有掺杂原子时，NEGF 的 SCF 比纯净体系难收敛。可以把 `scf.Max.Mixing.Weight` 降到 0.02、加大 `NEGF.SCF.Iter.Band`，或改用 `NEGF.Poisson.Solver FFT`。
- 想得到更光滑的透射谱，或要在本征通道里加入更多能量点时，可以在输入中加 `NEGF.tran.SCF.skip on`，只重算 Step 3（不要改 System.Name）。
