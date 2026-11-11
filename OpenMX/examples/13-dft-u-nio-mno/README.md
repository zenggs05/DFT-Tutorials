# 示例 13：NiO / MnO 的 DFT+U（Dudarev、U+J 与 cFLL、U 扫描）

配套教程：[第13章 DFT+U（Hubbard U 与 Hund J）](../../tutorials/13-dft-u.md)

## 目的

用两个 AFM-II 岩盐氧化物比较 GGA 与 GGA+U：

| 输入 | 体系 | 方法 | 演示内容 |
|---|---|---|---|
| `NiO_AFM_GGA.dat` | NiO，4 原子菱方胞 | GGA-PBE，无 +U | 基准：带隙偏小、Ni 磁矩偏小 |
| `NiO_AFM_U_Dudarev.dat` | 同上 | `scf.DFTU.Type 1`，U_eff = 5 eV（Ni 1d） | Dudarev 简化形式；也是 `scan_U.sh` 的模板 |
| `NiO_AFM_UJ_cFLL.dat` | 同上 | `scf.DFTU.Type 2`，U = 6 eV、J = 1 eV，`scf.dc.Type cFLL` | 一般形式与双计数项 |
| `MnO_AFM_U.dat` | MnO，4 原子菱方胞 | `scf.DFTU.Type 1`，U_eff = 4 eV（Mn 1d） | d⁵ 高自旋：磁矩接近 5 μB |

## 文件说明

| 文件 | 说明 |
|---|---|
| `NiO_AFM_GGA.dat` | GGA 基准，`RMM-DIISK` 混合 |
| `NiO_AFM_U_Dudarev.dat` | GGA+U（Dudarev），`RMM-DIISH` 混合；坐标行第 8 列为轨道极化开关 |
| `NiO_AFM_UJ_cFLL.dat` | GGA+U+J（一般形式，cFLL） |
| `MnO_AFM_U.dat` | MnO GGA+U |
| `scan_U.sh` | 以 Dudarev 输入为模板扫描 U，调用 `gap_from_out.py` 输出 `scan_U_summary.csv` |
| `gap_from_out.py` | 从 `.out` 提取 U、`Utot`、`Uhub`、SCF 网格上的 VBM/CBM/带隙与 Mulliken 磁矩（仅用 Python 标准库） |
| `job_slurm.sh` | Slurm 作业脚本，依次运行四个输入并汇总 |

所有输入使用 DFT_DATA19：`Ni6.0S-s3p2d2` + `Ni_PBE19S`（16 个价电子）、`Mn6.0-s3p2d2` + `Mn_PBE19`（15）、`O6.0-s2p2d1` + `O_PBE19`（6）。

## 运行步骤

1. 修改 `DATA.PATH`，指向同时包含 `PAO/` 与 `VPS/` 的目录（输入中默认写的是官方 deb 包的安装位置）：

   ```bash
   cd OpenMX/examples/13-dft-u-nio-mno
   export OMX_DATA=$HOME/openmx4.0/DFT_DATA19       # your path
   sed -i "s#^DATA.PATH .*#DATA.PATH                   ${OMX_DATA}#" *.dat
   grep "^DATA.PATH" *.dat
   ```

2. 依次计算（每个 4 原子体系在 8 核上通常几分钟到几十分钟，取决于机器）：

   ```bash
   export OMP_NUM_THREADS=1
   mpirun -np 8 openmx NiO_AFM_GGA.dat       -nt 1 | tee nio_gga.std
   mpirun -np 8 openmx NiO_AFM_U_Dudarev.dat -nt 1 | tee nio_dudarev_u5.std
   mpirun -np 8 openmx NiO_AFM_UJ_cFLL.dat   -nt 1 | tee nio_cfll_u6_j1.std
   mpirun -np 8 openmx MnO_AFM_U.dat         -nt 1 | tee mno_u4.std
   # or on a cluster
   sbatch job_slurm.sh
   ```

3. U 扫描（可选）：

   ```bash
   NP=8 bash scan_U.sh "0 2 4 5 6 7 8"
   RESTART=1 NP=8 bash scan_U.sh "0 2 4 5 6 7 8"   # chained restarts, fewer SCF steps
   ```

4. 用 OpenMX 3.9 时，删除各输入中的 `scf.Generation.Kpoint` 一行。

## 预期结果 / 如何检查

```bash
grep -A3 "DFT+U Type and DC" nio_dudarev_u5.std nio_cfll_u6_j1.std   # Type 1(Simplified) / 2(General) + cFLL
grep -E "Uhub\.|Utot\." *.out                                         # Uhub must be non-zero with +U
grep -n "Occupation Number in LDA+U" nio_dudarev_u5.out               # occupation-matrix eigenvalues
python3 gap_from_out.py nio_gga.out nio_dudarev_u5.out nio_cfll_u6_j1.out mno_u4.out
```

判断标准（定性；本仓库不附参考输出，具体数值取决于版本、编译与并行设置）：

- 每个 `.std` 最后的 SCF 步满足 `scf.criterion`，没有在 `scf.maxIter` 处被迫停止；
- NiO、MnO 中两个金属原子的 Mulliken `Diff` 大小相等、符号相反，两个 O 的 `Diff` 为 0，`Total spin moment` ≈ 0；
- NiO：GGA+U 的带隙明显大于 GGA（文献中 PBE 通常不到 1 eV，PBE+U 在 U_eff ≈ 5–6 eV 时约 3 eV；实验约 4.3 eV），Ni 磁矩随 U 增大；
- MnO：Mn 磁矩接近（略小于）5 μB，几乎不随 U 变化，带隙随 U 增大；
- `scan_U.sh` 中 U = 0 的 `Utot` 应与 `NiO_AFM_GGA.dat` 一致（同一基组、同一参数）。

`gap_from_out.py` 只用 SCF k 网格上的本征值，给出的是带隙估计；正式结果请用能带（第 8 章）或 DOS（第 9 章）确认。PDOS 的做法见教程第 11 节。

## 注意事项

- `<Hubbard.U.values` / `<Hund.J.values` 的项数和顺序必须与 PAO 定义一致；换基组（例如改用 Standard 的 `Ni6.0S-s3p2d1`）时要同步删去 `2d` 项；
- 块内不要插入注释行或空行，结束标记要紧跟最后一行数据；
- 不要比较不同 U 下的 `Utot`；
- k 点（7×7×7）、截断能（250 Ry）只是演示值，正式计算请先做收敛测试（第 5 章）；
- 这里的 U 值只用于教学，U 与投影定义绑定，不能直接照搬其他程序的数值。
