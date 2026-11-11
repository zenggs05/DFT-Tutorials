# 示例 08：CP2K 色散校正、自旋极化与 DFT+U

配套教程：[第8章 色散校正、自旋极化与 DFT+U](../../tutorials/08-dispersion-spin-dftu.md)

## 目的

- 用苯二聚体比较 PBE 与 PBE-D3(BJ) 的分子间相互作用能，并用 `RUN_TYPE BSSE` 做 counterpoise 校正；
- 用 O2 演示 `UKS` + `MULTIPLICITY`，比较三重态与闭壳层解；
- 用 NiO 演示反铁磁（两个 KIND + 相反的 `MAGNETIZATION`）和 DFT+U，分别给出 2026.2 的 k 点版本和适用于旧版本的 Γ 点超胞版本；
- 用石墨演示 PBE-D4（或 D3(BJ)）对层间距的影响（CELL_OPT）。

## 文件说明

| 文件 | 说明 | 版本要求 |
|---|---|---|
| `benzene_dimer_d3bj.inp` | 苯二聚体/单体单点；`@SET SYS` 取 dimer、a 或 b，`@SET VDW` 取 D3BJ 或 NONE | 2023.1+ |
| `benzene_dimer_bsse.inp` | `RUN_TYPE BSSE`，两个 `&FRAGMENT`（原子 1–12、13–24） | 2023.1+ |
| `benzene_dimer.xyz`、`benzene_a.xyz`、`benzene_b.xyz` | 平行错位构型（层间 3.50 Å，错位 1.60 Å），已置于 20 Å 盒子中心；单体与二聚体坐标完全一致 | — |
| `o2_triplet_uks.inp` | O2 的 UKS 几何优化，`@SET MULT 3`（可改为 1） | 2023.1+ |
| `nio_afm_dftu_kp.inp` | NiO AFM-II，4 原子菱方胞，k 点 4×4×4，PBE+U（Mulliken），`@SET UEFF 6.0` | **2026.2+** |
| `nio_afm_dftu_gamma.inp` | NiO AFM-II，64 原子立方超胞，Γ 点 + OT，含 `U_RAMPING` | 2023.1+ |
| `nio_afm_64.xyz` | 64 原子坐标，第一列为 KIND 名 `Ni1`/`Ni2`/`O`（晶胞 8.34 Å，写在输入中） | — |
| `graphite_d4_cell_opt.inp` | AB 石墨 CELL_OPT；`@SET VDW` 取 D4 或 D3BJ，初始 c = 7.20 Å | D4 需链接 libdftd4 |
| `job_slurm.sh` | 通用 Slurm 作业脚本 | — |

## 运行步骤

```bash
# environment, then check that the data files and optional libraries are available
cp2k.psmp -v | grep -i -o -E "libdftd4|libxc|spglib" | sort -u
ls "$CP2K_DATA_DIR"/dftd3.dat

# 1) benzene dimer: 3 systems x 2 dispersion settings (edit @SET lines with sed;
#    "-E VAR=value" on the command line is overridden by @SET inside the file)
for s in dimer a b; do for v in D3BJ NONE; do
  sed -e "s/^@SET SYS .*/@SET SYS $s/" -e "s/^@SET VDW .*/@SET VDW $v/" \
      benzene_dimer_d3bj.inp > bz_${s}_${v}.inp
  mpirun -np 8 cp2k.psmp -i bz_${s}_${v}.inp -o bz_${s}_${v}.out
done; done
mpirun -np 8 cp2k.psmp -i benzene_dimer_bsse.inp -o benzene_dimer_bsse.out

# 2) O2 triplet and closed-shell UKS solution
mpirun -np 4 cp2k.psmp -i o2_triplet_uks.inp -o o2_mult3.out
sed 's/^@SET MULT .*/@SET MULT 1/' o2_triplet_uks.inp > o2_mult1.inp
mpirun -np 4 cp2k.psmp -i o2_mult1.inp -o o2_mult1.out

# 3) NiO: k-points (2026.2) or Gamma supercell (any version)
sbatch --export=ALL,INP=nio_afm_dftu_kp.inp job_slurm.sh
sbatch --export=ALL,INP=nio_afm_dftu_gamma.inp job_slurm.sh

# 4) graphite (use VDW D3BJ if your build has no libdftd4)
sbatch --export=ALL,INP=graphite_d4_cell_opt.inp job_slurm.sh
```

## 预期结果 / 如何检查

```bash
grep "Dispersion energy:" bz_dimer_D3BJ.out
grep "BSSE-free interaction energy" benzene_dimer_bsse.out          # hartree
grep "Ideal and single determinant S\*\*2" o2_mult3.out | tail -1
grep "ENERGY| Total FORCE_EVAL" o2_mult3.out o2_mult1.out | awk '{print $1, $NF}'
grep "DFT+U energy:" nio_afm_dftu_kp.out | tail -1
grep -A8 "Mulliken Population Analysis" nio_afm_dftu_kp.out | tail -8
grep "CELL| Vector c" graphite_cellopt_D4.out | tail -1
```

- 苯二聚体：PBE 的相互作用能很弱（可能为正），PBE-D3(BJ) 为明显的负值。S22 基准（CCSD(T)/CBS）约 −2.7 kcal/mol，但本例几何不是 S22 参考几何，只做定性比较。教程第 3.3 节给出了用 Python 计算 $E_\mathrm{int}$ 的命令。
- O2：三重态能量低于 `MULT 1` 解；三重态 $\langle S^2\rangle$ 接近 2.0；两个 O 的 Mulliken `Spin moment` 各约 1；键长比实验值 1.2075 Å 略长（PBE 的典型行为）。
- NiO：Ni1、Ni2 的 `Spin moment` 大小相近、符号相反，总自旋约为 0（AFM 解）。实验磁矩约 1.6–1.9 μB（文献范围），Mulliken 值依赖基组，只作半定量对比。
- 石墨：c 明显小于初始 7.20 Å，接近实验 6.711 Å（具体数值以计算为准）；纯 PBE 时层间几乎不结合。

## 注意事项

- `&PAIR_POTENTIAL/TYPE` 的默认值在 2025.2（`DFTD3`）与 2026.2（`DFTD3(BJ)`）不同，所有示例都显式写了 `TYPE`。
- `U_MINUS_J` 的默认单位是 hartree，示例中一律写 `[eV]`。
- `nio_afm_dftu_kp.inp` 需要 2026.2（k 点 DFT+U 是 2026.2 的新功能）；旧版本请用 `nio_afm_dftu_gamma.inp`。
- U = 6.0 eV 只是演示值，CP2K（Mulliken 投影）中的 U 与 VASP 等程序的数值不能直接互换，正式工作请按教程第 8.5 节标定。
- NiO 示例使用实验晶格常数 4.17 Å，未做晶胞优化。
- `nio_afm_64.xyz` 的第一列是 KIND 名（`Ni1`/`Ni2`），部分可视化软件无法识别，查看结构时可临时把它们替换为 `Ni`。
- 石墨和 NiO 的 k 点计算使用对角化（2026.2 及以前 OT 不支持 k 点）。
