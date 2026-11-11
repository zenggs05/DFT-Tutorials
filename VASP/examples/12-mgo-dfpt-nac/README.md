# 示例 12-2：MgO 的 DFPT 力常数与非解析项修正（LO–TO 劈裂）

## 目的

- 用 VASP DFPT（`IBRION = 8`）在 2×2×2 超胞中一次计算全部力常数，再交给 phonopy；
- 用 `LEPSILON = .TRUE.` 计算 ε∞ 与 Born 有效电荷，`phonopy-vasp-born` 生成 BORN 文件；
- 比较有、无 NAC 的声子谱，观察极性材料 Γ 点的 LO–TO 劈裂。

对应教程：[第12章 声子：phonopy 有限位移与 DFPT](../../tutorials/12-phonon.md)（§8–§9）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `POSCAR-unitcell` | 岩盐 MgO 惯用胞（8 原子，Mg 在前），a = 4.25 Å（PBE 近似值）；请先优化 |
| `INCAR_relax` | 紧致优化（`ISIF = 3`、`EDIFF = 1E-8`、`EDIFFG = -1E-3`、`ENCUT = 520`） |
| `INCAR_dfpt`、`KPOINTS_dfpt` | 64 原子超胞 DFPT：`IBRION = 8`、`EDIFF = 1E-8`；**不设 NCORE/NPAR、不设 IVDW**；Γ 中心 3×3×3 |
| `INCAR_born`、`KPOINTS_born` | 单胞 `LEPSILON = .TRUE.`、`LRPA = .FALSE.`；Γ 中心 8×8×8 |
| `band.conf` | `FORCE_CONSTANTS = READ`、`NAC = .TRUE.`、fcc `PRIMITIVE_AXES`、Γ–X–W–K–Γ–L |
| `setup_mgo.sh` | 不带参数：建立 `dfpt/` 与 `born/`；`post`：生成 FORCE_CONSTANTS、BORN 并画声子谱 |
| `POTCAR` | **不提供**；例如 `cat $VASP_PP_PATH/potpaw_PBE/{Mg_pv,O}/POTCAR > POTCAR`（元素顺序 Mg、O） |

## 运行步骤

```bash
cd VASP/examples/12-mgo-dfpt-nac
cat $VASP_PP_PATH/potpaw_PBE/{Mg_pv,O}/POTCAR > POTCAR
grep TITEL POTCAR

# 0) 优化单胞（略，同示例 12-1；完成后覆盖 POSCAR-unitcell）

# 1) 准备 DFPT 与 Born 计算
bash setup_mgo.sh
(cd dfpt && sbatch ../../12-si-phonopy-fd/job.slurm)
(cd born && sbatch ../../12-si-phonopy-fd/job.slurm)

# 2) 后处理
bash setup_mgo.sh post
#   = phonopy-init --fc dfpt/vasprun.xml
#   + (cd born && phonopy-vasp-born --pa "0 1/2 1/2 1/2 0 1/2 1/2 1/2 0") > BORN
#   + phonopy --config band.conf -p -s

# 3) 对比：不加 NAC
phonopy --config band.conf --nonac -p -s        # phonopy >= 4；旧版去掉 conf 中的 NAC 行
```

## 预期结果 / 如何检查

```bash
grep -c 'name="hessian"' dfpt/vasprun.xml           # >= 1
grep -A 4 "MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)" born/OUTCAR
grep -A 9 "BORN EFFECTIVE CHARGES" born/OUTCAR | head -10
cat BORN
```

- VASP Wiki 的 MgO 示例给出 ε∞ ≈ 3.13、Z\*(Mg) ≈ +1.97、Z\*(O) ≈ −1.97，可作为量级参考；
- 所有原子的 Z\* 之和应接近 0；
- phonopy 日志出现 `NAC parameters were read from "BORN"`；
- 有 NAC 时 Γ 点光学支劈裂为 TO 与 LO（LO 明显更高，实验 TO/LO 约 12/21 THz 量级，近似值），无 NAC 时三支光学模在 Γ 点简并；
- 声子谱无虚频。

## 注意事项

- DFPT（IBRION=7/8）和 LEPSILON 不支持 `NCORE`/`NPAR` > 1，常见报错为 `PEAD routines do not work for NCORE` 或提示删除 NPAR；
- DFPT 不支持 vdW 修正与杂化泛函，这类情况改用有限位移（示例 12-1 的流程）与 LCALCEPS；
- `phonopy-vasp-born` 的 `--pa` 必须与 band.conf 的 `PRIMITIVE_AXES` 一致，且只接受数值矩阵（phonopy 4.10 实测不接受 `auto`/`F`）；
- 旧版 phonopy（2.x/3.x）的后处理写成 `phonopy -c POSCAR-unitcell -p -s band.conf`（`setup_mgo.sh post` 会自动选择）；
- VASP 6.6.1 修复了 OpenACC GPU 版 LCALCEPS 的 Born 电荷错误；phonopy 4.8.0 修复了从 vaspout.h5 读取 Born 电荷时的转置错误（默认读 vasprun.xml 不受影响）。
