# 示例 12-3：石墨烯声子谱（二维材料的有限位移计算）

## 目的

- 演示二维材料声子计算的特殊设置：真空层、面内超胞（`DIM = 5 5 1`）、z 方向 1 个 k 点、只优化面内晶格；
- 观察 ZA（面外弯曲）支在 Γ 附近的二次色散，以及可能出现的小虚频"口袋"；
- 复用示例 12-1 的 `run_disp.sh` 完成全部步骤。

对应教程：[第12章 声子：phonopy 有限位移与 DFPT](../../tutorials/12-phonon.md)（§12）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `POSCAR` | 石墨烯原胞（2 个 C），a = 2.467 Å（PBE 近似值），c = 20 Å 真空 |
| `INCAR_relax`、`KPOINTS_relax` | 面内优化：`ISIF = 3` + `LATTICE_CONSTRAINTS = .TRUE. .TRUE. .FALSE.`（VASP ≥ 6.4）；Γ 中心 24×24×1 |
| `INCAR`、`KPOINTS` | 50 原子超胞单点力；Γ 中心 6×6×1（包含折叠后的 K 点） |
| `band.conf` | Γ–M–K–Γ，`DIM = 5 5 1` |
| `POTCAR` | **不提供**：`cat $VASP_PP_PATH/potpaw_PBE/C/POTCAR > POTCAR` |

## 运行步骤

```bash
cd VASP/examples/12-graphene-phonon
cat $VASP_PP_PATH/potpaw_PBE/C/POTCAR > POTCAR

# 1) 面内晶格优化
mkdir -p relax
cp INCAR_relax relax/INCAR; cp KPOINTS_relax relax/KPOINTS; cp POSCAR POTCAR relax/
(cd relax && sbatch ../../12-si-phonopy-fd/job.slurm)
cp relax/CONTCAR POSCAR
sed -n 3,5p POSCAR                 # 检查 c 轴仍为 20 A

# 2) 位移超胞（只需 1 个位移）、计算、收集
DIM="5 5 1" bash ../12-si-phonopy-fd/run_disp.sh prepare
(cd disp-001 && sbatch ../../12-si-phonopy-fd/job.slurm)
DIM="5 5 1" bash ../12-si-phonopy-fd/run_disp.sh collect

# 3) 声子谱
phonopy --config band.conf -p -s    # phonopy 2.x/3.x: phonopy -p -s band.conf
```

## 预期结果 / 如何检查

- 6 支声子：3 支声学（ZA、TA、LA）+ 3 支光学（ZO、TO、LO）；
- ZA 支在 Γ 附近呈二次色散（ω ∝ q²），TA、LA 为线性；
- Γ 点最高光学模（E2g）约 1580 cm⁻¹ ≈ 47 THz 量级（对应 Raman G 峰，近似比较）；
- K 点附近最高光学支有 Kohn 异常，对 k 点和 SIGMA 敏感；
- 若 ZA 支只在 Γ 极近处出现很小的虚频，通常是数值问题而非真实不稳定。

## 注意事项

- VASP < 6.4 没有 `LATTICE_CONSTRAINTS`：删去该行，用 `ISIF = 2` 做 E–a 扫描确定晶格常数。
- 减小 ZA 小虚频的办法：更大的超胞（6×6×1、8×8×1）、`EDIFF = 1E-8`、`PREC = Accurate`、对力常数施加对称化与求和规则（新版 phonopy 默认的 symfc 投影有帮助）。
- 不要在真空方向扩胞；KPOINTS 的 z 方向保持 1，六角晶格必须用 Γ 中心网格。
- 二维材料的极性 NAC 处理与体相不同，本例（石墨烯非极性）不涉及。
