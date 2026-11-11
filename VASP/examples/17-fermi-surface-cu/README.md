# 示例 17-4：fcc Cu 的费米面（IFermi / vaspkit + FermiSurfer）

## 目的

- 在稠密均匀 k 网格上做非自洽计算，得到画费米面所需的本征值；
- 用 IFermi 插值并绘制三维费米面、计算费米面面积和费米速度；
- 用 vaspkit 26x 生成 FermiSurfer / XCrySDen 格式文件作为对照。

对应教程：[第17章 高级电子结构](../../tutorials/17-advanced-electronic.md)（§6）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `POSCAR` | fcc Cu 原胞，a = 3.615 Å（实验值） |
| `INCAR_scf` / `KPOINTS_scf` | 自洽，`ISMEAR = 1`、`SIGMA = 0.1`，Γ 中心 15×15×15 |
| `INCAR_fs` / `KPOINTS_fs` | 非自洽 `ICHARG = 11`，Γ 中心 25×25×25，`LORBIT = 11` |
| `job.slurm` | 自洽 → 非自洽 → `ifermi info` / `ifermi plot` |
| `POTCAR` | **不提供**：`cat $VASP_PP_PATH/potpaw_PBE/Cu/POTCAR > POTCAR` |

## 运行步骤

```bash
cd VASP/examples/17-fermi-surface-cu
cat $VASP_PP_PATH/potpaw_PBE/Cu/POTCAR > POTCAR
sbatch job.slurm
```

IFermi（推荐在单独的 conda 环境中安装，依赖 pymatgen 与 BoltzTraP2）：

```bash
conda create -n ifermi python=3.11 pip cmake numpy && conda activate ifermi
conda install -c conda-forge pymatgen boltztrap2 pyfftw
pip install ifermi

cd fs
ifermi info                                   # 费米面个数、面积、维度
ifermi info --property velocity               # 平均费米速度
ifermi plot --output cu_fs.html               # 交互式三维图（plotly）
ifermi plot --interpolation-factor 10 --property velocity --output cu_fs_velocity.html
ifermi plot --slice 1 1 0 0 --output cu_fs_slice110.png   # (110) 截面；保存图片需 plotly 的导出组件或改用 --type matplotlib
```

vaspkit（在 `fs/` 目录中）：

```bash
vaspkit      # 26 -> 262：写 FERMISURFACE.bxsf（XCrySDen）；26 -> 263：FermiSurfer 格式
```

## 预期结果 / 如何检查

- `ifermi info` 应报告一个穿过费米能级的能带形成的费米面（Cu 的 4s 带，标准 POTCAR 下为第 6 条带）；
- 三维图中是一个近似球形、在 L 点方向形成"颈"（neck）并与相邻布里渊区连通的费米面，这是贵金属费米面的经典特征（de Haas–van Alphen 实验早已证实）；
- 网格太粗时费米面会出现锯齿或颈部断开：增大 `KPOINTS_fs` 或 `--interpolation-factor`。

## 注意事项

- vaspkit 的教程要求在**原胞**上计算费米面，并用 261 生成 KPOINTS；本例的 25×25×25 网格可替换为 vaspkit 261 生成的网格；
- 投影自旋（`--property spin`）需要非共线计算且 `ISYM = -1` 的完整网格；
- 半导体需要用 `--mu` 把能量移到导带或价带中，否则得不到费米面；
- EIGENVAL、PROCAR、vasprun.xml 在稠密网格下体积较大，不要提交到仓库。
