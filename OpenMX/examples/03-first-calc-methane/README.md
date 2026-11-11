# 示例 03：第一个 OpenMX 计算——CH4、H2O、CO 分子（First OpenMX Calculation）

配套教程：[第3章 输入文件结构与第一个计算（甲烷分子）](../../tutorials/03-input-first-calc.md)

## 目的

1. 用一个最小基组的甲烷算例（`Methane_manual_check.dat`）与 OpenMX 4.0 手册第 5 节的参考总能量对比，确认安装与数据库正确；
2. 用 Standard 基组计算 CH4、H2O、CO 三个分子，练习阅读 `.std`、`.out`，提取总能量、HOMO/LUMO、Mulliken 电荷与偶极矩；
3. 用 `MO.fileout` 输出 CO 的前线分子轨道 cube 文件并可视化。

## 文件说明

| 文件 | 说明 |
|---|---|
| `Methane_manual_check.dat` | 与官方 `work/Methane.dat` 物理设置相同（`H5.0-s1`/`C5.0-s1p1`、GGA-PBE、120 Ry、10 Å 盒子、C 在原点），`System.Name = metref`，只用于安装验证 |
| `Methane.dat` | CH4 单点能：Standard 基组 `H6.0-s2p1`/`C6.0-s2p2d1`，200 Ry，C–H = 1.089 Å，位于 10 Å 盒子中心，`System.Name = met` |
| `H2O.dat` | H2O 单点能：实验几何（O–H 0.9572 Å，∠HOH 104.52°），12 Å 盒子，`System.Name = h2o` |
| `CO.dat` | CO 单点能 + 分子轨道输出（`MO.fileout on`、`num.HOMOs 2`、`num.LUMOs 2`），C–O 1.128 Å，`System.Name = co` |
| `job_omx_slurm.sh` | 依次运行上述 4 个输入；可 `sbatch` 提交，也可在工作站上直接 `bash` 运行 |
| `check_results.sh` | 汇总当前目录下所有 `.out`：是否正常结束、SCF 步数、Utot（Ha/eV）、HOMO/LUMO 与能隙（eV）、偶极矩、Mulliken 电子数；对 `metref` 还会与手册参考值比较 |

## 运行步骤

```bash
# 0) copy the example to a scratch directory
mkdir -p ~/omx_runs/ch03 && cd ~/omx_runs/ch03
cp /path/to/DFT-Tutorials/OpenMX/examples/03-first-calc-methane/* .

# 1) point DATA.PATH to your DFT_DATA19 (absolute path containing PAO/ and VPS/)
export OMX_DATA=$HOME/software/openmx4.0/DFT_DATA19      # deb package: usually /opt/openmx4.0/DFT_DATA19
sed -i "s#^DATA.PATH.*#DATA.PATH                     ${OMX_DATA}#" *.dat
grep -H "^DATA.PATH" *.dat

# 2a) run one input interactively
export OMP_NUM_THREADS=1
mpirun -np 4 openmx Methane.dat -nt 1 | tee met.std

# 2b) or run everything (workstation)
NP=4 bash job_omx_slurm.sh
# 2c) or submit to Slurm (edit the module lines first)
sbatch job_omx_slurm.sh

# 3) summary
bash check_results.sh
```

`job_omx_slurm.sh` 在设置了 `OMX_DATA` 时会自动改写每个输入的 `DATA.PATH`，并在运行前检查该目录下是否有 `PAO/` 与 `VPS/`。若 `openmx` 不在 `PATH` 中，用 `OPENMX=/full/path/to/openmx bash job_omx_slurm.sh`。

## 预期结果 / 如何检查

- 每个 `.std` 的最后一行应为 `The calculation was normally finished.`（`check_results.sh` 的 `end` 列为 `ok`）；
- `nSCF` 应明显小于 `scf.maxIter`（100）；这几个闭壳层小分子通常在几十步以内收敛；
- `MulP/ideal` 两个数应相等：CH4 与 H2O 为 8，CO 为 10；
- **`metref`**：`Utot` 应与手册第 5 节参考值 −8.033540953061 Hartree 一致（脚本以 1e-6 Ha 为判据，正常安装通常好得多）；手册中该算例 16 步收敛；
- **CH4**：`HOMO = 4`，HOMO 为三重简并（本征值第 2–4 个几乎相同），偶极矩接近 0；
- **H2O**：`HOMO = 4`，偶极矩应在约 1.8–1.9 D（实验值约 1.85 D；具体数值随基组与截断能略有变化）；
- **CO**：`HOMO = 5`，偶极矩很小（实验约 0.11 D）；目录中会多出 `co.homo0_0.cube`、`co.homo0_1.cube`、`co.lumo0_0.cube`、`co.lumo0_1.cube`，可以用 VESTA 或 OpenMX Viewer 打开查看 σ 与 π 型轨道；
- 也可以手动检查：

```bash
grep "Utot\." *.out
grep -A 12 "Eigenvalues (Hartree) for SCF KS-eq." met.out
sed -n '/Mulliken populations/,/Decomposed Mulliken/p' h2o.out
grep -A 8 "Dipole moment (Debye)" co.out
awk '/<coordinates.forces/,/coordinates.forces>/' met.out
```

本目录不附带任何计算输出；以上除 `metref` 的参考值来自官方手册外，均为定性预期，请以你自己的计算结果为准。

## 注意事项

- 输入文件中的 `DATA.PATH` 是占位符 `/path/to/openmx4.0/DFT_DATA19`，**不修改直接运行会因找不到 PAO/VPS 而退出**；
- `System.CurrrentDirectory` 的官方拼写是三个 r，不要“改正”成两个 r；
- 四个输入的 `System.Name` 各不相同（`metref`、`met`、`h2o`、`co`），可以放在同一目录中运行而不互相覆盖；
- 小分子用 1–4 个 MPI 进程即可，进程更多反而更慢；
- 若用 Windows 编辑器修改过 `.dat`，请确认仍是 ASCII + LF（`file *.dat`、`dos2unix *.dat`）；
- 最小基组（`Methane_manual_check.dat`）只能用来验证安装，不能用于研究。
