# 算例 14-2：GaAs 的静态介电张量与 Born 有效电荷（LEPSILON / DFPT）

## 目的

1. 用密度泛函微扰理论（`LEPSILON = .TRUE.`）计算 GaAs 的电子（钳位离子）静态介电常数 ε∞，比较"不含 / 含局域场效应"的差别，并读出 Born 有效电荷 Z*；
2. 再加上 `IBRION = 8`（DFPT 声子），得到离子贡献 ε_ion 与静态介电常数 ε0 = ε∞ + ε_ion，并用 Lyddane–Sachs–Teller（LST）关系做自洽检查。

这两组数据也正是极性材料声子计算中非解析修正（NAC，phonopy 的 BORN 文件）所需要的输入，见 [第12章](../../tutorials/12-phonon.md)。

对应教程：[第14章 介电与光学性质](../../tutorials/14-optics-dielectric.md)。

## 文件说明

| 文件 | 作用 |
|---|---|
| `POSCAR` | 闪锌矿 GaAs 原胞，采用实验晶格常数 a = 5.653 Å（PBE 在自身平衡晶格下带隙接近 0，用实验晶格可得到有限带隙，便于演示） |
| `KPOINTS` | Γ 中心 12×12×12（线性响应 / PEAD 例程要求 Γ 中心或奇数网格） |
| `INCAR` | Part A：`LEPSILON = .TRUE.`，`IBRION = -1` |
| `INCAR.ionic` | Part B：`LEPSILON = .TRUE.` + `IBRION = 8` |
| `job.slurm` | 在 `A_lepsilon/`、`B_ionic/` 两个子目录中依次运行，并 grep 关键结果 |

POTCAR 由作业脚本按 `cat $VASP_PP_PATH/potpaw_PBE/{Ga_d,As}/POTCAR > POTCAR` 生成（仓库不提供 POTCAR）。

## 运行步骤

```bash
cd VASP/examples/14-dielectric-gaas-lepsilon
export VASP_PP_PATH=/path/to/your/paw_library
sbatch job.slurm
```

> [!WARNING]
> 不要在这两个 INCAR 里设置 `NCORE` 或 `NPAR`（保持默认 NCORE = 1），否则线性响应 / PEAD 例程会报错退出。可以用 `KPAR` 做 k 点并行。

## 预期结果 / 如何检查

```bash
cd A_lepsilon
grep -A 5 "MACROSCOPIC STATIC DIELECTRIC TENSOR" OUTCAR     # 依次列出不同近似下的 eps_inf
grep -A 10 "BORN EFFECTIVE CHARGES" OUTCAR                  # 每个离子一个 3x3 张量
grep -A 6 "PIEZOELECTRIC TENSOR" OUTCAR                      # 闪锌矿有非零的 e14
cd ../B_ionic
grep -A 5 "MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONTRIBUTION" OUTCAR
grep "THz" OUTCAR                                           # Gamma 点频率（TO 模、3 个约为 0 的声学模）
```

判断标准：

- 立方对称：介电张量和 Born 电荷张量都应是对角且三个对角元相等；Ga 与 As 的 Z* 大小相等、符号相反（电荷中性的声学求和规则），偏差较大说明 k 点不够。
- 同一 OUTCAR 中会有多组 ε∞：不含局域场（excluding local field effects）、只含 Hartree 局域场（RPA）、含 Hartree + 交换关联局域场（in DFT）。Hartree 局域场一般使 ε∞ 减小，再加入交换关联项一般又使其增大（见官方 `LRPA` 说明）；与实验比较时通常取"including local field effects in DFT"一组。
- 量级参考（实验，近似值）：ε∞ ≈ 10.9，ε0 ≈ 12.9，Z* ≈ 2.1。PBE 低估带隙，因此 ε∞ 往往明显偏大；这是泛函误差而不是输入错误。
- LST 关系：ε0/ε∞ = (ω_LO/ω_TO)²，可与实验 TO/LO 频率比较；`IBRION = 8` 只给出 Γ 点不含长程场的 TO 频率，LO 频率需要非解析修正。

## 注意事项

- `LEPSILON` 不需要增加 `NBANDS`，增加反而更慢。
- 对含 vdW 修正（`IVDW`）的体系，DFPT 声子（`IBRION = 7/8`）不可用，可改用 `IBRION = 6` + `LEPSILON` 或 `LCALCEPS`。
- 非共线（SOC）体系：VASP < 6.5.0 中 `LEPSILON` + 对称性的结果有误（known issue 20，旧版本用 `ISYM = -1`）；GPU 上非共线线性响应也有已知问题（known issue 114，可用 CPU、`OMP_NUM_THREADS=1` 规避）。
- 生成 phonopy 的 BORN 文件：`phonopy-vasp-born` 读取 Part A（或 Part B）目录中的 vasprun.xml，详见 [第12章](../../tutorials/12-phonon.md)。

返回教程：[第14章 介电与光学性质](../../tutorials/14-optics-dielectric.md)
