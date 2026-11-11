# 示例 13-1：O2 与 CH2 的自旋多重度比较（UKS、RKS、ROKS 与 Guess=Mix）

配套教程：[第13章 开壳层体系、波函数稳定性与 Broken-Symmetry](../../tutorials/13-openshell-stability.md)（§4、§5）

## 目的

- 学会为开壳层分子写正确的电荷与多重度，比较三重态和单重态的能量；
- 读懂 `<S**2>`（湮灭前/后）和 `Stable` 的输出；
- 用 `Guess=Mix` 得到 O2 的 broken-symmetry（BS）单重态，并用求和法 $E_{\rm S}\approx 2E_{\rm BS}-E_{\rm T}$ 估算 a¹Δg 激发能；
- 计算 CH2 的单重态-三重态能隙（含零点能），并比较 UKS 与 ROKS。

## 文件说明

| 文件 | 体系 | 步骤 | 演示内容 |
|---|---|---|---|
| `o2_triplet_vs_singlet.gjf` | O2 | 3 | ① UKS 三重态 Opt Freq；② RKS 闭壳层单重态 + `Stable`；③ UKS `Guess=Mix SCF=NoSymm Stable=Opt` BS 单重态 |
| `ch2_carbene_singlet_triplet.gjf` | CH2 | 4 | ① 三重态 ³B1 Opt Freq；② 单重态 ¹A1 Opt Freq；③ 单重态 `Stable`；④ 三重态 ROB3LYP 单点 |
| `run_g16_slurm.sh` | — | — | Slurm 作业脚本（8 核、10 GB，对应输入中的 `%nprocshared=8`、`%mem=8GB`） |

方法：B3LYP/def2-TZVP。起始结构按典型键长、键角手工给出（O–O 1.208 Å；CH2 三重态 H–C–H 134°，单重态 102°），由各步 `Opt` 优化。

## 运行步骤

本地：

```bash
cd Gaussian/examples/13-spin-states-o2-ch2
g16 < o2_triplet_vs_singlet.gjf > o2_triplet_vs_singlet.log
g16 < ch2_carbene_singlet_triplet.gjf > ch2_carbene_singlet_triplet.log
```

Slurm：

```bash
sbatch run_g16_slurm.sh o2_triplet_vs_singlet.gjf
sbatch run_g16_slurm.sh ch2_carbene_singlet_triplet.gjf
```

脚本里的 `module load` 一行需要按集群实际情况修改；scratch 默认放在 `/tmp/$USER`，可用环境变量 `SCR_BASE` 指定其他位置。两个作业都很小，8 核上一般几分钟内完成。

## 预期结果 / 如何检查

```bash
# 1) all job steps finished: 3 for O2, 4 for CH2
grep -c "Normal termination" o2_triplet_vs_singlet.log ch2_carbene_singlet_triplet.log

# 2) energies and <S**2> of every step
grep "SCF Done" o2_triplet_vs_singlet.log | tail -3
grep "S\*\*2 before annihilation" o2_triplet_vs_singlet.log | tail -2

# 3) stability messages
grep -i "instability\|stable under\|already stable" *.log

# 4) zero-point corrected energies for the CH2 gap
grep "Sum of electronic and zero-point Energies" ch2_carbene_singlet_triplet.log
```

应当看到（定性判据，具体数值取决于版本与机器）：

| 检查项 | 预期 |
|---|---|
| O2 第 1 步 ⟨S²⟩ | 非常接近 2.00，频率无虚频 |
| O2 第 2 步 | 稳定性分析报告 `RHF -> UHF instability`（闭壳层 (π*)² 不是最低单重态解） |
| O2 第 3 步 ⟨S²⟩ | 约 1.0（BS 态 = 单重态与三重态 $M_S=0$ 分量各半） |
| O2 能量顺序 | $E_{\rm T} < E_{\rm BS} < E_{\rm RKS}$ |
| CH2 | 三重态能量低于单重态；单重态 `Stable` 应报告稳定；ROKS 的 ⟨S²⟩ 严格为 2.0000，能量略高于 UKS |

O2 的 a¹Δg 激发能估算：$\Delta E \approx 2(E_{\rm BS}-E_{\rm T})$，与实验约 0.98 eV（NIST Chemistry WebBook，近似值）比较。CH2 的能隙用两个 `Sum of electronic and zero-point Energies` 之差，与实验 $T_0\approx 9.0$ kcal/mol（Jensen & Bunker, *J. Chem. Phys.* 1988, 89, 1327）比较。1 Hartree = 627.5095 kcal/mol = 27.2114 eV。

批量汇总也可用 [../../scripts/g16_log_summary.py](../../scripts/g16_log_summary.py)（`python ../../scripts/g16_log_summary.py -h` 查看用法）。

## 注意事项

- 开壳层单重态必须显式写 `UB3LYP`；多重度为 1 时 Gaussian 默认按 RKS 计算。
- 每步用 `%oldchk` 读、`%chk` 写，避免后一步覆盖前一步的 chk。
- `Geom=AllCheck` 的步骤在 route 之后只有一个空行，不写标题和电荷多重度。
- 输入中的标题和说明为英文 ASCII，中文解释都放在本 README 和教程中。
- 文件末尾必须有空行，否则 Gaussian 可能报 `End of file` 类错误。

返回教程：[第13章](../../tutorials/13-openshell-stability.md)
