# 示例 16-4：单层 MoS2 的部分电荷（PARCHG）与 STM 模拟

## 目的

- 完成"自洽（写 WAVECAR）→ `LPARD` 后处理"的两步流程；
- 分别得到占据态（E_F−1 eV … E_F，对应负偏压）和空态（导带底附近，对应正偏压）的部分电荷；
- 用 Tersoff–Hamann 近似生成恒高、恒流 STM 图像：任意版本可用 `stm_from_parchg.py`，VASP ≥ 6.5 还可以用 `LPARDH5` + py4vasp。

对应教程：[第16章 电荷与成键分析](../../tutorials/16-charge-bonding.md)（§7）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `POSCAR` | 2H-MoS2 单层，a = 3.19 Å（接近 PBE 值，实验约 3.16 Å），S–S 层间高度 3.13 Å，c = 23 Å（约 20 Å 真空）；未经优化，正式计算请先按 [第8章](../../tutorials/08-vdw-2d.md) 优化 |
| `KPOINTS` | Γ 中心 12×12×1 |
| `INCAR_scf` | 自洽：`NBANDS = 24`（留足空带）、`LWAVE = .TRUE.` |
| `INCAR_parchg_occ` | `LPARD`、`NBMOD = -3`、`EINT = -1.0 0.0` |
| `INCAR_parchg_empty` | `EINT = 1.0 3.0`（相对 E_F，**需按自己的导带底调整**） |
| `stm_from_parchg.py` | 从 PARCHG 生成恒高 / 恒流 STM 图（numpy + matplotlib） |
| `plot_stm_py4vasp.py` | VASP ≥ 6.5 且 `LPARDH5 = .TRUE.` 时用 py4vasp 画 STM |
| `job.slurm` | 依次运行自洽、两个 LPARD 步骤和作图 |
| `POTCAR` | **不提供**：`cat $VASP_PP_PATH/potpaw_PBE/{Mo_pv,S}/POTCAR > POTCAR` |

## 运行步骤

```bash
cd VASP/examples/16-stm-mos2
cat $VASP_PP_PATH/potpaw_PBE/{Mo_pv,S}/POTCAR > POTCAR
grep -E "TITEL|ZVAL" POTCAR          # Mo_pv: ZVAL = 12；S: ZVAL = 6 -> NELECT = 24

sbatch job.slurm
```

手动执行时的关键命令：

```bash
# 自洽
mkdir scf && cp INCAR_scf scf/INCAR && cp POSCAR KPOINTS POTCAR scf/ && cd scf
mpirun -np 16 vasp_std > vasp.out && cd ..

# 先确定窗口：E-fermi、价带顶、导带底
grep E-fermi scf/OUTCAR | tail -1
grep -E "val. band max|cond. band min|fundamental gap" scf/OUTCAR | tail -3   # VASP >= 6.4.3

# 部分电荷（复制 WAVECAR，INCAR 中 ENCUT/NBANDS/ISMEAR 与自洽一致）
mkdir parchg_occ && cp INCAR_parchg_occ parchg_occ/INCAR && cp POSCAR KPOINTS POTCAR scf/WAVECAR parchg_occ/
cd parchg_occ && mpirun -np 4 vasp_std > vasp.out
python3 ../stm_from_parchg.py PARCHG --mode height --height 3.0 --repeat 6 -o stm_height.png
```

VASP ≥ 6.5.0（HDF5 编译）：在 `INCAR_parchg_*` 中取消 `LPARDH5 = .TRUE.` 的注释，此时不再写 PARCHG，而是写入 `vaspout.h5`：

```bash
pip install py4vasp kaleido
python3 ../plot_stm_py4vasp.py --mode constant_height --tip-height 3.0 --supercell 6
python3 ../plot_stm_py4vasp.py --mode constant_current --current 1.0 --supercell 6
```

## 预期结果 / 如何检查

- `parchg_occ/PARCHG` 存在且不全为 0：`awk 'NR>15 && NF==5' PARCHG | head`；
- 若空态的 PARCHG 很小或全为 0，说明 `EINT` 窗口落在带隙内：根据 `cond. band min` 重新设置（`NBMOD = -3` 时窗口相对 E_F）；
- `stm_from_parchg.py` 会打印恒高平面上的密度范围或恒流模式下的起伏（corrugation）；
- 单层 MoS2 的模拟图像呈六角对称排列的亮点，亮点位于顶层 S 还是 Mo 上方取决于偏压窗口和针尖高度，请结合能带的轨道投影（`LORBIT = 11`）解释，不要只凭图像判断。

## 注意事项

- 只设 `LPARD = .TRUE.` 而不设 `EINT`/`IBAND` 时，VASP 只写一个不含 augmentation 的 CHGCAR，不会产生 PARCHG；
- LPARD 步骤的 `ENCUT`、`NBANDS`、`ISPIN`、`ISYM`、k 点必须与生成 WAVECAR 的自洽计算相同；
- `KPUSE`/`LSEPK` 会把 k 点权重设为 1，需要自洽与后处理都用 `ISYM = -1`；
- 非共线（SOC）计算不支持 LPARD；
- STM 图像的亮暗与偏压、针尖高度、平滑参数有关，只适合定性比较。
