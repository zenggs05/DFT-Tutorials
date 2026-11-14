# 示例 05：Si 的能带结构与态密度（ICHARG=11 两步法 + KPOINTS_OPT 一步法）

## 目的

- 用 SCF → `ICHARG=11` 的传统两步法计算 Si 的 PBE 能带和密网格 DOS/PDOS；
- 用 `KPOINTS_OPT`（VASP ≥ 6.3）在一次运行中得到同样的能带，与两步法对比；
- 用 py4vasp、sumo、VASPKIT 或自写脚本作图并读取带隙。

对应教程：[第5章 能带结构与态密度](../../tutorials/05-band-dos.md)。

## 文件说明

| 路径 | 说明 |
|---|---|
| `1-scf/INCAR, POSCAR, KPOINTS` | 自洽：Γ 中心 8×8×8，`LCHARG=.TRUE.`；POSCAR 用 a = 5.47 Å（请换成[第4章](../../tutorials/04-relax-eos.md)优化得到的值） |
| `2-band-icharg11/INCAR, KPOINTS` | 非自洽能带：`ICHARG=11`、`NBANDS=16`、`LORBIT=11`；Line-mode 路径 L–Γ–X–U\|K–Γ，每段 40 点 |
| `3-dos/INCAR, KPOINTS` | 密网格 DOS：Γ 中心 16×16×16、`ISMEAR=-5`、`LORBIT=11`、`NEDOS=3001` |
| `4-band-kpoints-opt/INCAR, KPOINTS, KPOINTS_OPT` | 一步法：KPOINTS 为均匀网格，KPOINTS_OPT 为同一条路径 |
| `job.slurm` | 依次运行 4 个步骤，自动复制 POSCAR、POTCAR、CHGCAR |
| `plot_with_py4vasp.ipynb` | py4vasp 作图与数据导出（需要 HDF5 版 VASP 生成的 vaspout.h5） |
| `POTCAR` | **不提供**，在本目录生成一次即可 |

## 运行步骤

```bash
cd VASP/examples/05-band-dos-si
cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR
sbatch job.slurm
```

手动逐步运行（便于理解每一步）：

```bash
cd 1-scf && cp ../POTCAR . && mpirun -np 8 vasp_std > vasp.out && cd ..
cd 2-band-icharg11 && cp ../POTCAR ../1-scf/{POSCAR,CHGCAR} . && mpirun -np 8 vasp_std > vasp.out && cd ..
cd 3-dos && cp ../POTCAR ../1-scf/{POSCAR,CHGCAR} . && mpirun -np 8 vasp_std > vasp.out && cd ..
cd 4-band-kpoints-opt && cp ../POTCAR ../1-scf/POSCAR . && mpirun -np 8 vasp_std > vasp.out && cd ..
```

## 预期结果 / 如何检查

```bash
grep "E-fermi" 1-scf/OUTCAR | tail -1                      # 作图用的费米能
grep -A12 "Band structure" 1-scf/OUTCAR                    # VASP >= 6.4.3 的带隙报告
grep -A7 "k-point     1 :" 1-scf/OUTCAR | tail -6          # SCF 的 Gamma 点本征值
grep -A7 "k-point    40 :" 2-band-icharg11/OUTCAR | tail -6 # 路径上的 Gamma 点，应与上一行一致（< 1 meV）
sed -n 6p 3-dos/DOSCAR                                     # EMAX EMIN NEDOS E-fermi
ls 4-band-kpoints-opt/PROCAR_OPT 4-band-kpoints-opt/vaspout.h5
```

作图：

```bash
jupyter notebook plot_with_py4vasp.ipynb                     # py4vasp
cd 2-band-icharg11 && sumo-bandplot && cd ..                 # sumo
cd 3-dos && sumo-dosplot && cd ..
python ../../scripts/plot_band.py -h                         # 本仓库脚本（查看用法）
python ../../scripts/plot_dos.py -h
```

- 能带应呈现 Si 典型的间接带隙：VBM 在 Γ，CBM 在 Γ–X 线上靠近 X 处。PBE 带隙约 0.5–0.6 eV（文献量级；实验约 1.17 eV），被低估属正常。
- 两步法与 KPOINTS_OPT 的能带应几乎重合。
- DOS 中 E_F 处的积分态密度应等于 8（Si 原胞的价电子数）。

## 注意事项

- 两步法中 `2-band-icharg11`、`3-dos` 的 `ENCUT`、`PREC`、POTCAR、POSCAR 必须与 `1-scf` 完全一致；磁性或 DFT+U 体系还要保证 `ISPIN`、`MAGMOM`、`LDAU*`、`LMAXMIX` 一致。
- 能带作图请使用 `1-scf`（或 `3-dos`）的 E-fermi，不要用路径计算的 E-fermi。
- `KPOINTS_OPT` 不能用 `vasp_gam`；`PROCAR_OPT` 在 6.5.0/6.5.1 中有 bug（6.6.0 修复），自旋极化 + PDOS 请用 6.6.1。
- 高对称点坐标对应本 POSCAR 的晶胞约定（Setyawan–Curtarolo fcc 原胞）；换了晶胞请重新生成路径（教程 §5）。
- KPOINTS 中的标签只能用 ASCII（如 `GAMMA`），不要写 Γ。
