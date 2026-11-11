# 示例 09：CP2K 杂化泛函与 ADMM（PBE0、HSE06、截断库仑）

配套教程：[第9章 杂化泛函与 ADMM](../../tutorials/09-hybrid-admm.md)

## 目的

- 孤立水分子的 PBE0 + ADMM 几何优化（库仑算符，无截断）；
- 32 水分子周期盒子：先 PBE 得到波函数，再以其为初猜做 PBE0-TC-LRC + ADMM（截断库仑 + 长程修正，`SCREEN_ON_INITIAL_P T`）；
- 64 原子 Si 超胞的 HSE06 + ADMM（`SHORTRANGE` 算符；2026.2 起 `CUTOFF_RADIUS` 生效），同一文件内包含 PBE 预计算阶段。

## 文件说明

| 文件 | 说明 | 版本 |
|---|---|---|
| `h2o_pbe0_admm.inp` | H2O，PBE0 + ADMM（`BASIS_PROJECTION` + `MO_DIAG`），GEO_OPT | 2023.1+ |
| `water32_pbe_step1.inp` | 32 H2O（9.8528 Å 立方盒子），PBE 单点，写出 `water32_pbe-RESTART.wfn` | 2023.1+ |
| `water32_pbe0_admm_truncated.inp` | 32 H2O，PBE0-TC-LRC + ADMM，`CUTOFF_RADIUS 4.5`，从第一步重启 | 2023.1+ |
| `water32.xyz` | 32 个水分子，fcc 格点 + 随机取向生成，**未平衡**（最近 O–O 3.48 Å） | — |
| `si_hse06_admm.inp` | Si64，`@SET STAGE PBE` 或 `HSE06`；HSE06 阶段从 PBE 波函数重启，`ADMM_TYPE ADMM2` | HSE06 截断半径需 2026.2 |
| `si64.xyz` | 2×2×2 常规胞金刚石 Si，a = 5.431 Å（实验值），晶胞 10.862 Å | — |
| `job_slurm.sh` | Slurm 脚本，默认 16 MPI × 2 OpenMP（杂化泛函需要更多每进程内存） | — |

所有输入使用 UZH 的杂化专用基组和赝势：`*-MOLOPT-PBE0-GTH-qN`（别名 `-HYB-`）+ `GTH-PBE0-qN`，ADMM 辅助基组 `admm-dzp-qN`（`BASIS_ADMM_UZH`）。

## 运行步骤

```bash
# 0) environment; make sure BASIS_ADMM_UZH and t_c_g.dat are in the data directory
ls "$CP2K_DATA_DIR"/BASIS_ADMM_UZH "$CP2K_DATA_DIR"/t_c_g.dat
cp2k.psmp -c water32_pbe0_admm_truncated.inp

# 1) isolated water (minutes)
mpirun -np 4 cp2k.psmp -i h2o_pbe0_admm.inp -o h2o_pbe0_admm.out

# 2) water box: PBE first, then PBE0 restarted from its wavefunction
mpirun -np 32 cp2k.psmp -i water32_pbe_step1.inp -o water32_pbe_step1.out
export OMP_NUM_THREADS=2
mpirun -np 16 cp2k.psmp -i water32_pbe0_admm_truncated.inp -o water32_pbe0.out
#    or: sbatch --export=ALL,INP=water32_pbe0_admm_truncated.inp job_slurm.sh

# 3) Si64: PBE stage, then HSE06 stage
sed 's/^@SET STAGE .*/@SET STAGE PBE/' si_hse06_admm.inp > si64_pbe.inp
mpirun -np 32 cp2k.psmp -i si64_pbe.inp -o si64_pbe.out          # writes si64_PBE-RESTART.wfn
mpirun -np 16 cp2k.psmp -i si_hse06_admm.inp -o si64_hse06.out
```

`MAX_MEMORY`（每个 MPI 进程，MiB）按节点内存调整：示例中 water32 为 2500、Si64 为 3000，即 16 个进程各需约 2.5–3 GB 的 HFX 内存，另加程序本身的占用。

## 预期结果 / 如何检查

```bash
grep "Hartree-Fock Exchange energy:" *.out | awk '{print $1, $NF}'
grep "ADMM|" h2o_pbe0_admm.out | head
grep "HFX_MEM_INFO| Number of sph. ERI's calculated on the fly" water32_pbe0.out | head -1
grep -i "larger than half the minimal cell dimension" *.out          # should be empty
grep -i -A3 "cutoff radius" si64_hse06.out | head                    # see the note below
grep "GEOMETRY OPTIMIZATION COMPLETED" h2o_pbe0_admm.out
grep -i "SCF run NOT converged" *.out                                # should be empty
```

- 每个杂化计算的能量分解中都应出现非零的 `Hartree-Fock Exchange energy:`。
- water32 PBE0：从 PBE 重启后 SCF 通常十几步内收敛；`calculated on the fly` 为 0 或很小说明 `MAX_MEMORY` 够用。第一步 SCF 明显慢于后续各步（ERI 计算），属于正常现象。
- Si64 HSE06（2026.2）：可能出现提示截断半径小于 erfc 判据所需半径的警告，这是预期的（5 Å 小于自动值约 14 Å），截断半径需要按教程第 7.2 节做收敛测试；在 2026.1 及以前，`CUTOFF_RADIUS` 对 `SHORTRANGE` 不起作用，会出现 “larger than half the minimal cell dimension” 警告。
- Γ 点超胞的 HOMO–LUMO 差不是 Si 的间接带隙；带隙请用第10章的 DOS/能带方法计算（实验值约 1.17 eV）。

## 注意事项

- PBE 预计算与杂化计算必须使用同一主基组、赝势和原子顺序，否则波函数重启无效。
- `SCREEN_ON_INITIAL_P T` 只能在从重启文件起步时使用；`si_hse06_admm.inp` 的 HSE06 阶段依赖 `si64_PBE-RESTART.wfn`。
- `water32.xyz` 只用于演示输入，没有经过平衡；研究工作请先平衡结构。
- 截断库仑要求 `CUTOFF_RADIUS` 不超过晶胞最小晶面间距的一半（本例 L/2 = 4.93 Å 与 5.43 Å）。
- 2026.2 新增的 `&HF/&ACE` 需要对角化 SCF，不能直接加到这些使用 OT 的示例中，见教程第 9 节。
- k 点杂化（RI-HFXk）的输入片段见教程第 8 节，本目录未提供可直接运行的 k 点杂化示例。
