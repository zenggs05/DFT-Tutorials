# 示例 10-3：VASPsol 隐式水溶剂与线性化 Poisson–Boltzmann 电解质

## 目的

- 演示 VASPsol 的标准两步法：真空单点（写 WAVECAR）→ `ISTART = 1` + `LSOL = .TRUE.` 的溶剂单点；
- 得到 CO/Pt(111) 体系的溶剂化能 $E_\mathrm{sol}-E_\mathrm{vac}$；
- （可选）加入 `LAMBDA_D_K` 的线性化 Poisson–Boltzmann 电解质，查看 `FERMI_SHIFT`，为恒电势方法做准备。

对应教程：[第10章 吸附、催化与自由能](../../tutorials/10-adsorption-catalysis.md)（§11）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `POSCAR` | CO/Pt(111) top 位初始结构（与示例 10-1 的 `3_co_top/POSCAR` 相同）；生产计算请换成优化后的 CONTCAR |
| `KPOINTS` | Γ 中心 6×6×1 |
| `INCAR_vac` | 第 1 步：真空单点，`LWAVE = .TRUE.` |
| `INCAR_sol` | 第 2 步：`ISTART = 1`、`LSOL = .TRUE.`、`EB_K = 78.4`、`TAU = 0.000525`（VASPsol 对水的默认值） |
| `INCAR_lpb` | 第 3 步（可选）：再加 `LAMBDA_D_K = 3.0`（Debye 长度，约对应 1 M 1:1 电解质） |
| `run_vaspsol.sh` | 依次在 `vac/`、`sol/`（、`lpb/`）中运行并打印溶剂化能 |

## 运行步骤

**前提**：VASP 必须打过 VASPsol 补丁（[henniggroup/VASPsol](https://github.com/henniggroup/VASPsol) 的 README；编译要点见[第1章](../../tutorials/01-setup-compile.md)与第10章 §11.1）。

```bash
cd VASP/examples/10-vaspsol-implicit-solvent
cp ../10-co-pt111-adsorption/3_co_top/CONTCAR POSCAR     # 推荐：用优化好的结构
export VASP_PP_PATH=/path/to/your/paw/library
bash ../10-co-pt111-adsorption/make_potcar.sh .

# 在作业脚本或交互式节点中运行（run_vaspsol.sh 内部调用 mpirun）
bash run_vaspsol.sh            # vac -> sol
bash run_vaspsol.sh lpb        # vac -> sol -> lpb
```

用 Slurm 时，把 `bash run_vaspsol.sh` 写进作业脚本（可参考 `../10-co-pt111-adsorption/job.slurm`，把其中的 `mpirun` 行替换掉）。

## 预期结果 / 如何检查

- 脚本最后打印 `E(solvent) - E(vacuum) = ... eV`，对金属表面上的小分子通常是负的、量级在零点几 eV 以内（依体系而定，此处仅为量级描述）；
- `sol/` 的 SCF 步数明显多于 `vac/`；
- 脚本会列出 `sol/` 中与溶剂化有关的输出行；若两次能量几乎相同且找不到任何溶剂化相关输出，说明可执行文件没有包含 VASPsol；
- `lpb/vasp.out` 中应能找到 `FERMI_SHIFT`，把它加到 E_F 上才得到以电解质深处电势为零的费米能级。

## 注意事项

- `NC_K`、`SIGMA_K`、`TAU` 的默认值是针对水拟合的；换溶剂只改 `EB_K` 不够。
- 空腔能（TAU 项）对 FFT 网格敏感：比较 ENCUT = 450 与 500 eV 的溶剂化能，确认收敛。
- 计算吸附能的溶剂化校正时，slab、吸附体系都要做"真空 + 溶剂"两步；气相分子的参考态按实验条件选择（见教程 §11.3）。
- 原 VASPsol 仓库已停止开发，新版本为 VASPsol++（[VASPsol/VASPsol](https://github.com/VASPsol/VASPsol)），关键词可能不同。
- 本示例不加 `LDIPOL`；如需与偶极校正同时使用，请先做对照测试。
