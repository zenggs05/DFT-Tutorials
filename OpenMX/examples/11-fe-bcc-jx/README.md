# 算例 11-fe-bcc-jx：用 jx 计算 bcc Fe 的交换耦合 J_ij 与平均场居里温度

## 目的

- 完成"共线 SCF（`HS.fileout on`）→ jx 配置文件 → `jx` → 分壳与 T_C 估算"的完整流程；
- 学会由原胞基矢推导 `<ijpairs.cellid` 中的格矢指标 (l1, l2, l3)；
- 用"同一壳层内 J 应相等"检查 k 点与极点数是否足够。

对应教程：[第11章 自旋极化、非共线磁性与交换耦合](../../tutorials/11-spin-noncollinear.md) §11、§12。

## 文件说明

| 文件 | 说明 |
|---|---|
| `Fe_bcc_jx.dat` | bcc Fe 原胞，Fe6.0H-s2p2d2 / Fe_PBE19H（手册建议 jx 使用较小基组），GGA-PBE 共线自旋极化，k 网格 16³，`HS.fileout on` |
| `Fe_bcc_jx.config` | jx 配置：`Flag.PeriodicSum off`、`Num.Poles 60`、`Num.Kgrid 16 16 16`，列出前三个近邻壳层共 26 对 |
| `jx_tc_mf.py` | 解析 jx 输出，按距离分壳，给出截断到所列壳层的平均场 T_C（numpy） |
| `job.slurm` | SCF → jx → `jx_tc_mf.py` 一条龙；用环境变量 `JX` 指定 jx 可执行文件路径 |

运行前修改 `DATA.PATH`，并编译 jx（只需一次）：

```bash
cd ~/openmx4.0/source && make jx      # the executable is copied to ../work/jx
```

## 运行步骤

```bash
sbatch job.slurm
# or step by step:
mpirun -np 8 openmx Fe_bcc_jx.dat -nt 1 > fe_jx.std 2>&1
mpirun -np 16 ~/openmx4.0/work/jx fe_jx.scfout Fe_bcc_jx.config > jx.log
python jx_tc_mf.py jx.log Fe_bcc_jx.dat
```

## 预期结果 / 如何检查

- `fe_jx.scfout` 生成（二进制，几 MB 量级）；jx 启动时会打印读取的 scfout 格式版本与 `Previous eigenvalue solver = Band`、`atomnum = 1`；
- `jx.log` 中每行 `i j c1 c2 c3 J[meV] J[mRy] ...`；同一壳层的 8（或 6、12）个 J 应在数值噪声范围内相等，`jx_tc_mf.py` 的 `spread` 列给出其离散度；
- 第一、二近邻 J 为正（铁磁耦合）且为 meV 到十几 meV 量级是合理的；手册图 43 给出了 bcc Fe 的 J 随距离的完整曲线；
- `jx_tc_mf.py` 给出的 T_C 只包含所列 3 个壳层，是粗略估计。手册用数百对 (i, j, R) 得到 bcc Fe 的平均场 T_C 为 1321 K（实验 1040 K）。

具体数值以你的计算为准。

## 注意事项

- 4.0 的 jx 只支持**共线**计算；需要 DMI、各向异性交换或 SOC 时用 TB2J-OpenMX（`pip install TB2J-OpenMX`，`openmx2J.py --prefix fe_jx --elements Fe --kmesh ...`）。
- `Num.Kgrid` 应不小于 SCF 的 k 网格；把 `Num.Poles` 改为 100 检查 J 的变化是否 < 0.05 meV。
- jx 只对 k 点做 MPI 并行，进程数不要超过 k 点总数。
- 内存不够时调小 `Bunch.ij.pairs`。
- 想加入更多壳层：按 R = l1·a1 + l2·a2 + l3·a3 计算距离，补充行并同步修改 `Num.ij.pairs`。

[返回第11章](../../tutorials/11-spin-noncollinear.md)
