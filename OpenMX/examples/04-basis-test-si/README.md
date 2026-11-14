# 示例 04-1：Si 的基组收敛测试（Quick / Standard / Precise）

## 目的

- 在相同结构、截断能和 k 点下，比较 DFT_DATA19 推荐表中 Si 的三档 PAO 基组：`Si7.0-s2p1d1`（Quick，10 个基函数/原子）、`Si7.0-s2p2d1`（Standard，13 个）、`Si7.0-s3p3d2`（Precise，22 个）；
- 比较总能、晶格常数（EvsLC 扫描 + Birch–Murnaghan 拟合）和能带，体会"基组也需要收敛测试"。

对应教程：[第4章 基组 PAO 与赝势 VPS 的选择](../../tutorials/04-pao-vps.md)（§9）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `Si_basis_quick.dat` | Quick 基组，SCF + L–Γ–X–W–K–Γ 能带 |
| `Si_basis_standard.dat` | Standard 基组（与另两个文件只差基组和 `System.Name` 两行） |
| `Si_basis_precise.dat` | Precise 基组 |
| `run_basis_test.sh` | `scf`：依次计算并汇总 Utot；`eos`：生成并运行 EvsLC 扫描输入；`collect`：只汇总 |
| `job.slurm` | Slurm 作业脚本，可提交单个输入，也可在作业内运行驱动脚本 |

## 运行步骤

```bash
cd OpenMX/examples/04-basis-test-si
# 1) point DATA.PATH to your DFT_DATA19 (the directory that contains PAO/ and VPS/)
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" *.dat

# 2) SCF + band for the three basis sets (interactive node)
MPI_RUN="mpirun -np 4" NT=1 bash run_basis_test.sh scf
#    or on a cluster:
sbatch job.slurm run_basis_test.sh scf

# 3) lattice constant for each basis (EvsLC, 12 points from a = 5.24 Ang)
MPI_RUN="mpirun -np 4" bash run_basis_test.sh eos
python ../06-cellopt-eos-si-c/fit_evslc.py si_standard_eos.EvsLC --lattice fcc --natoms 2
```

## 预期结果 / 如何检查

1. `run_basis_test.sh scf` 最后打印一张表：每档基组的 Utot（Hartree/原胞）以及相对 Precise 的能量差（meV/atom）。应满足 Quick > Standard > Precise（基组越大能量越低）。数据库 Si 页面给出的参照是：E–V 曲线极小处 `Si7.0-s2p2d1` 与 `Si7.0-s3p3d2f1` 相差 0.025 eV/atom。
2. 表下面的 `Num. of grids` 三行必须相同，否则说明比较的不是同一个积分网格。
3. 用 `grep -o "scf_conv=[01]" si_*.md` 确认每次计算都是 `scf_conv=1`。
4. EvsLC 拟合：Standard 与 Precise 的 a0 应相差约 0.01 Å 量级；PBE 文献值约 5.47 Å，实验 5.431 Å。
5. 能带：用 bandgnu13 或本仓库脚本 [`omx_band_plot.py`](../../scripts/omx_band_plot.py)（`python ../../scripts/omx_band_plot.py -h` 查看用法）画 `si_*.Band`。三档基组的价带应几乎重合，导带对基组更敏感。Si 的 PBE 带隙为间接带隙、约 0.5–0.6 eV 量级（实验约 1.17 eV），被低估属正常。

## 注意事项

- 这里所有计算共用 `scf.energycutoff 220`、`scf.Kgrid 8 8 8`。截断能和 k 点本身的收敛测试在 [第5章示例](../05-scf-conv-si/README.md) 中完成。
- `eos` 模式用 GNU sed 改写输入：把晶格从 a = 5.43 Å 缩到 5.24 Å，设 `MD.Type EvsLC`、`MD.EvsLC.Step 0.7`、`MD.maxIter 12`，并加入 `scf.Ngrid 32 32 32`。固定网格是为了避免 E–V 曲线因网格点数跳变而出现锯齿（第6章 §5）。
- 结果文件 `si_*_eos.EvsLC` 的前 9 列是三个晶格矢量（Å），最后一列是总能（Hartree）。
- 输入文件只含 ASCII 字符；在 Windows 下编辑后请确认仍是 LF 换行（`dos2unix *.dat`）。
