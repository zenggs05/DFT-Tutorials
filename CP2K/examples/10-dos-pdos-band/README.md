# 示例 10：DOS、PDOS 与能带（H2O、金红石 TiO2、石墨烯、Si）

配套教程：[第10章 电子结构分析：DOS、PDOS 与能带](../../tutorials/10-dos-pdos-band.md)

## 目的

用四个小体系演示 CP2K 中态密度（DOS）、投影态密度（PDOS）和能带结构的完整流程，并对比 2026.2 的新 `&DFT/&PRINT/&DOS` 语法与 2026.1 及以前的旧 `&DFT/&PRINT/&PDOS` 语法：

| 输入 | 体系 | SCF 方式 | 演示内容 | 适用版本 |
|---|---|---|---|---|
| `h2o_pdos_2026.inp` | 孤立 H2O（12 Å 盒子，MT 求解器） | OT | `&DOS` + `NLUMO`（SCF 后用 OT 本征求解器补算空轨道）、`&PDOS COMPONENTS`、`&CURVE` 展宽、`ENERGY_ZERO HOCO` | ≥ 2026.2 |
| `tio2_pdos_legacy_2025.inp` | 金红石 TiO2 2×2×3 超胞（72 原子，Γ 点） | 对角化 + `ADDED_MOS` | 旧语法 `&PDOS NLUMO`、`&LDOS LIST`；文件末尾附 2026.2 等价写法 | ≤ 2026.1 |
| `graphene_dos_kpoints.inp` | 石墨烯原胞，真空 20 Å | 对角化 + Fermi–Dirac | k 点 DOS、`MP_GRID` 加密总 DOS、k 点 PDOS（仅曲线输出）、`ENERGY_ZERO FERMI` | ≥ 2026.2 |
| `si_band_structure.inp` | 体相 Si 原胞 | 对角化 + k 点 | `&BAND_STRUCTURE` + `&KPOINT_SET`（L–Γ–X–U \| K–Γ），同时输出以 VBM 为零点的 DOS | 能带部分各版本通用；`&DOS` 部分 ≥ 2026.2 |

## 文件说明

| 文件 | 说明 |
|---|---|
| `h2o_pdos_2026.inp` | H2O 分子 PDOS（新语法） |
| `tio2_pdos_legacy_2025.inp` | 金红石 PDOS（旧语法），注释中给出新语法替换块 |
| `graphene_dos_kpoints.inp` | 石墨烯 k 点 DOS/PDOS |
| `si_band_structure.inp` | Si 能带与 DOS |
| `job_slurm.sh` | Slurm 作业脚本，参数为输入文件名；运行前会删除旧的 `.bs` 文件（CP2K 以追加方式写能带文件） |

所有输入使用 `BASIS_MOLOPT_UZH` + `POTENTIAL_UZH`（PBE），基组与赝势的 `-qN` 成对匹配。

## 运行步骤

```bash
cd CP2K/examples/10-dos-pdos-band
export CP2K_DATA_DIR=/path/to/cp2k/data   # 若已内置数据路径可省略
cp2k.psmp --version | head -n 3           # 确认版本：2026.2 才能跑新语法输入

# 先做语法检查（版本不对时会直接报 unknown subsection / unknown keyword）
cp2k.psmp -c h2o_pdos_2026.inp

# 本地运行（4 个 MPI 进程）
mpirun -np 4 cp2k.psmp -i h2o_pdos_2026.inp       -o h2o_pdos_2026.out
mpirun -np 4 cp2k.psmp -i graphene_dos_kpoints.inp -o graphene_dos_kpoints.out
rm -f si_band.bs && mpirun -np 4 cp2k.psmp -i si_band_structure.inp -o si_band_structure.out

# 旧版本（<= 2026.1）上运行旧语法示例
mpirun -np 16 cp2k.psmp -i tio2_pdos_legacy_2025.inp -o tio2_pdos_legacy_2025.out

# 集群
sbatch job_slurm.sh si_band_structure.inp
```

作图使用本仓库脚本（用法见各自的 `-h`）：

```bash
python ../../scripts/pdos_plot.py -h    # 读旧版 .pdos 与 2026.2 曲线文件，按元素/角动量求和作图
python ../../scripts/band_plot.py -h    # 读 .bs 文件，按 KPOINT_SET 生成横轴和高对称点标签
```

## 预期结果 / 如何检查

通用检查：

```bash
grep "SCF run converged" *.out
grep "PROGRAM ENDED AT"  *.out
ls -1 *.dos *.pdos *.bs
```

- **H2O（2026.2）**：生成 `h2o_pdos-1.dos`（离散能级表）、`h2o_pdos-curve-1.dos`（展宽总 DOS）、`h2o_pdos-k1-1.pdos` / `h2o_pdos-k2-1.pdos`（逐态投影；`kN` 的编号与输出中 `ATOMIC KIND INFORMATION` 的顺序一致，用文件第一行 `# Projected DOS for atomic kind O ...` 确认是哪种元素）以及对应的 `*_curve-1.pdos`（展宽 PDOS，含 `total` 列）。文件名以你的实际输出为准（`ls` 查看）。曲线文件表头应含 `# Energy zero: HOCO`。在 O 的 PDOS 中，HOMO（第 4 个态，1b1）几乎完全来自 `px`，因为分子位于 yz 平面，1b1 轨道垂直于分子平面。
- **TiO2（旧语法）**：生成 `tio2_pdos-k1-1.pdos`、`tio2_pdos-k2-1.pdos`（分别对应 Ti 和 O，以文件第一行为准）和 `tio2_pdos-list1-1.pdos`（原子 1）。价带顶以 O 2p 为主，导带底以 Ti 3d 为主。旧格式只有离散态，需要自己展宽（见教程第 6 节的 Python 片段）。
- **石墨烯**：总 DOS 在费米能级附近应接近“V”形并在 0 eV 处趋于零（狄拉克点）；`pz` 分量主导 −2 到 +2 eV 范围。输出中会出现警告 `State-resolved k-point PDOS output is not implemented yet`，这是正常的：k 点 PDOS 只有展宽曲线文件。
- **Si**：`si_band.bs` 中有两个 `# Set` 块（91 个和 31 个 k 点）。以 VBM 为零点作图后，价带顶在 Γ，导带底位于 Γ–X 线上靠近 X 处，为间接带隙。PBE 带隙远小于实验值（1.17 eV，低温），且 DZVP 基组对导带的描述不如 TZV2P 完备，具体数值请以你的计算为准。

本仓库不附带参考输出文件；数值取决于 CP2K 版本、基组和并行设置。

## 注意事项

- **版本**：`&DOS` 下的 `&PDOS`/`&CURVE`/`&BROADEN`/`NLUMO` 只在 2026.2 及以后可用；2026.2 不再识别 `&PRINT/&PDOS`，会报 `unknown subsection PDOS of section PRINT`。反过来，新语法拿到旧版本上也会报错。
- **DELTA_E 没有单位**：2026.2 中 `DELTA_E` 是 Hartree 原子单位的纯数（默认 0.001），不要写 `[eV]`。`WIDTH` 有单位，可以写 `WIDTH [eV] 0.1`。
- **`&SMEAR` 一定要写 `METHOD`**：2026.2 中 `METHOD` 的默认值是 `GAUSSIAN`，2026.1 及以前是 `ENERGY_WINDOW`，不写会得到你不想要的展宽方式。
- **能带文件追加写入**：`&BAND_STRUCTURE FILE_NAME` 指定的文件以追加方式打开，重复运行前先删除旧文件，否则 `.bs` 中会有多套数据。
- **k 点与对称性**：分析 `COMPONENTS` 分量时不要打开 `&KPOINTS SYMMETRY`；对称约化会旋转轨道分量（官方文档的警告）。
- **基组**：DZVP 对价带足够，但空带（导带）对基组更敏感；定量比较带隙或导带形状时用 TZV2P 复查。
