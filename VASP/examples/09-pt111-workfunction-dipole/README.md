# 算例 09-pt111-workfunction-dipole：Pt(111) 功函数与偶极校正

## 目的

- 用 `LVHAR = .TRUE.` 输出静电势（离子势 + Hartree 势），沿表面法向做平面平均，得到真空能级 V_vac，功函数 Φ = V_vac − E_F；
- 对比 O 只吸附在一侧的**不对称** slab 在有/无偶极校正时真空区电势的差别：无校正时真空区存在人为电场（电势有斜率），有校正时两侧真空各自平坦、中间出现一个台阶，台阶高度就是两侧功函数之差；
- VASP ≥ 6.4.3：用 `LVACPOTAV` 让 VASP 直接在 OUTCAR 中给出两侧真空能级。

对应教程：[第9章](../../tutorials/09-surface-workfunction.md) §7–§8。

## 文件说明

| 文件 | 说明 |
|---|---|
| `INCAR` | 静态计算：`LDIPOL/IDIPOL=3/DIPOL`、`AMIN=0.01`、`LVHAR`、`LVACPOTAV`（≥ 6.4.3） |
| `KPOINTS` | 与 [09-pt111-slab](../09-pt111-slab/README.md) 相同的 Γ 中心 8×8×1 |
| `POSCAR_O_fcc` | 在未弛豫 Pt(111) p(2×2) 顶面 fcc 空位放 1 个 O（覆盖度 0.25 ML），底面干净——人为构造的不对称 slab |
| `plot_locpot_planar.py` | 读取 LOCPOT 做平面平均，自动识别真空区和偶极台阶，输出两侧 V_vac、Φ、真空区斜率；生成 `planar_avg.dat/png` |
| `run.sh` | `prepare`：建立 `clean/`、`O_nodip/`、`O_dip/` 三个计算并提交；`collect`：汇总 |
| `job.slurm` | 最简 Slurm 脚本 |

POTCAR（两种体系分开准备）：

```bash
cat $VASP_PP_PATH/potpaw_PBE/Pt/POTCAR > POTCAR_Pt
cat $VASP_PP_PATH/potpaw_PBE/{Pt,O}/POTCAR > POTCAR_PtO
```

## 运行步骤

```bash
# 0) relax the clean slab first in ../09-pt111-slab (its CONTCAR is used here)
bash run.sh prepare
# wait ...
bash run.sh collect
# single directory, with figure:
cd O_dip && python ../plot_locpot_planar.py LOCPOT --outcar OUTCAR
```

## 预期结果 / 如何检查

```bash
grep "E-fermi" O_dip/OUTCAR | tail -1
grep upper O_dip/OUTCAR | tail -1     # LVACPOTAV (>= 6.4.3): vacuum level upper / lower side
grep -i dipol O_dip/OUTCAR | tail -3  # dipole moment along z
```

- `clean/`：两侧功函数几乎相等（Δ Φ ≈ 0），真空区电势平坦。PBE 计算的 Pt(111) 功函数一般在 5.6–5.9 eV 范围（近似值；实验约 5.9–6.1 eV；系统的 PBE 数据库见 Tran et al., Surf. Sci. 687, 48 (2019)）。
- `O_nodip/`：脚本报告 `dipole step: not found`，并给出 `WARNING: potential is not flat`——真空区残留电场，此时"功函数"没有确定值。
- `O_dip/`：两侧真空各自平坦，O 覆盖一侧（upper）的功函数明显高于干净一侧（O 吸电子形成指向表面内侧的偶极，使该侧功函数升高，量级约零点几到 1 eV，示意）。
- 若 `LVACPOTAV` 打开，OUTCAR 中会有一行 `vacuum level on the upper side and lower side of the slab ...`，应与脚本给出的 V_vac 基本一致（两者对真空区的取法略有不同，差别通常在 0.01 eV 量级）。

## 注意事项

- `DIPOL` 必须放在 slab 的质心附近（只有法向分量有意义）。如果 slab 被放在晶胞底部、真空在上方，DIPOL 的 z 分量要相应修改，否则偶极层会落在原子之间。
- 偶极校正会让 SCF 收敛变慢（电荷在两侧之间来回振荡）。本例已设置 `AMIN = 0.01`、`NELM = 200`、`EDIFF = 1E-6`；仍不收敛时，先不加 `LDIPOL` 收敛波函数，再读入 WAVECAR 打开偶极校正。
- 功函数对真空厚度、层数、k 点都要做收敛；用 LVHAR（不含交换关联势）可以在较小真空下就得到平坦的真空能级。
- `POSCAR_O_fcc` 是未弛豫的演示结构；做真实研究时请先按 [第10章](../../tutorials/10-adsorption-catalysis.md) 优化吸附结构（同样要加偶极校正）。

[返回第9章](../../tutorials/09-surface-workfunction.md)
