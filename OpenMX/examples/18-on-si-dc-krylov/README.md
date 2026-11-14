# 示例 18-2：Si 超胞的 O(N) 计算——DC、DC-LNO、Krylov 的精度测试与 512 原子生产计算

配套教程：[第18章 并行、大体系 O(N)、作业脚本与常见错误排查](../../tutorials/18-parallel-on-errors.md)（§6）

## 目的

1. 在 64 原子 Si 超胞上，把三种 O(N) 求解器（`DC`、`DC-LNO`、`Krylov`）的总能与常规对角化（`Band`，4×4×4 k 点）比较，扫描 `orderN.HoppingRanges`，确定满足精度要求的截断球半径；
2. 用选定的半径对 512 原子超胞做 O(N) 自洽计算（`Si512_DC.dat`、`Si512_DC-LNO.dat`）；
3. 演示手册第 29.2 节的“O(N) 自洽 + 一次常规对角化”流程：读取 DC-LNO 的自洽密度，在 Γ 点对角化一次，得到 O(N) 方法给不出的本征值和 DOS（`Si512_onestep_band.dat`）；
4. 用 `make_si_supercell.py` 生成更大的超胞（如 4096 原子），在集群上做 O(N) 的扩展性测试。

## 文件说明

| 文件 | System.Name | 说明 |
|---|---|---|
| `Si64_Band_ref.dat` | `si64_band` | 参考：常规对角化，`scf.Kgrid 4 4 4`，`scf.eigen.lib elpa2` |
| `Si64_DC.dat` | `si64_dc` | DC，`orderN.HoppingRanges 7.0` |
| `Si64_DC-LNO.dat` | `si64_dclno` | DC-LNO，`orderN.HoppingRanges 8.0`，`orderN.LNO.Buffer 0.2`，`orderN.LNO.Occ.Cutoff 0.1`，`scf.dclno.threading on`；文件头注释给出了 `<LNOs.Num` 的写法 |
| `Si64_Krylov.dat` | `si64_krylov` | Krylov，`orderN.KrylovH.order 400`，其余三个开关取默认值 `on` |
| `Si512_DC.dat` | `si512_dc` | 512 原子 DC |
| `Si512_DC-LNO.dat` | `si512_dclno` | 512 原子 DC-LNO（第一步） |
| `Si512_onestep_band.dat` | `si512_dclno`（**与上一步相同**） | 第二步：`scf.restart on`、`scf.maxIter 1`、`Band` + `scf.Kgrid 1 1 1`，输出 DOS |
| `scan_hopping.sh` | — | 自动运行参考计算，再对三种求解器扫描 R = 7…11 Å，输出 `hopping_scan.csv`（ΔE 以 meV/atom 计） |
| `make_si_supercell.py` | — | 生成任意 n×n×n 的 Si 超胞输入（可选求解器、半径、随机位移），`python3 make_si_supercell.py -h` 查看用法 |

共同设置：GGA-PBE，`Si7.0-s2p1d1`（Quick 基组）+ `Si_PBE19`，晶格常数 5.43 Å 的理想金刚石结构，`scf.energycutoff 200`，电子温度 300 K，RMM-DIISK 混合，`scf.criterion 1.0e-8`，`level.of.fileout 0`。所有 Si64 输入的晶胞和截断能相同，因此实空间网格相同，总能可以直接相减。

## 运行步骤

```bash
cd OpenMX/examples/18-on-si-dc-krylov
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" *.dat

# 1. accuracy scan on Si64 (about 16 runs; use one node)
NP=16 NT=1 bash scan_hopping.sh
cat hopping_scan.csv

# 2. edit orderN.HoppingRanges in Si512_*.dat according to the scan, then
sbatch ../18-job-scripts/job_slurm.sh Si512_DC-LNO.dat
#   or interactively: mpirun -np 64 openmx Si512_DC-LNO.dat -nt 2 > si512_dclno.std 2>&1

# 3. one diagonalization on top of the O(N) density (same directory, same System.Name)
mpirun -np 64 openmx Si512_onestep_band.dat -nt 1 > si512_onestep.std 2>&1
DosMain si512_dclno.Dos.val si512_dclno.Dos.vec      # interactive; see chapter 9

# 4. bigger systems for scaling tests
python3 make_si_supercell.py 8 --solver DC-LNO --rcut 8.0 -o Si4096_DC-LNO.dat
```

## 预期结果 / 如何检查

```bash
grep "normally finished" *.std
grep -m1 -i "SNAN" si64_dc_R7.std          # FNAN/SNAN = atoms in each truncated cluster
grep "Utot\." si64_band.out si64_*_R*.out
awk '/Total Computational Time/{print FILENAME, $NF}' *.out
grep -A40 "SCF history at MD" si512_dclno.out | tail -n 5
```

- `hopping_scan.csv`：对 Si 这样有能隙的共价体系，随着 R 增大，|ΔE| 应近似指数下降（手册图 23 的规律）。具体数值以你的计算为准；一般把 |ΔE| ≤ 1–5 meV/atom 的最小 R 作为生产值。
- DC 与 DC-LNO 在相同 R 下精度相近（手册图 25），DC-LNO 的耗时应明显更短；Krylov 的精度还受 `orderN.KrylovH.order` 控制。
- 截断球中的原子数（FNAN+SNAN）随 R 增加大致按 $R^3$ 增长，单步 SCF 耗时也随之增加。
- `Si512_onestep_band.dat` 只做 1 步 SCF，应很快结束；生成 `si512_dclno.Dos.val/.Dos.vec`。注意它会覆盖 `si512_dclno.out`，需要保留第一步的 `.out` 时请先备份。

## 注意事项

- O(N) 方法**不输出波函数**（`MO.fileout` 会被关掉），本征值、能带、DOS、cube 轨道要用第 3 步的“一次对角化”获得。
- O(N) 方法依赖截断球内原子数，**改变晶胞（EOS、晶胞优化）时 FNAN+SNAN 会跳变，能量曲线出现台阶**，此时用 `<orderN.FNAN+SNAN` 块固定每个原子的团簇大小（教程 §6.6）。
- 512 原子的 DC/DC-LNO 计算对内存要求较高，在笔记本上不要尝试；作业脚本中务必保留 `ulimit -s unlimited`。
- 理想晶格的受力为零，只适合做单点能和扩展性测试；要做 MD，用 `make_si_supercell.py --rattle 0.03` 加随机位移，并用 `<LNOs.Num` 固定 LNO 数。
- 第一步和第二步必须在同一目录、用同一个 `System.Name`，并保留 `si512_dclno_rst/`。

返回：[第18章](../../tutorials/18-parallel-on-errors.md) | [OpenMX 目录](../../README.md)
