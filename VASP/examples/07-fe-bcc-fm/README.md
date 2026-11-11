# 算例 07-fe-bcc-fm：bcc Fe 铁磁自旋极化计算（FM / AFM / NM 对比）

## 目的

- 练习 `ISPIN = 2`、`MAGMOM`、`LORBIT = 11` 的基本用法；
- 学会从 OSZICAR / OUTCAR 读出总磁矩和各原子磁矩；
- 用同一个 2 原子立方胞比较铁磁（FM）、反铁磁（AFM）和非磁（NM）三种状态的能量，体会"初始磁矩决定最终磁态"。

对应教程：[第7章 磁性、DFT+U 与自旋轨道耦合](../../tutorials/07-magnetism-u-soc.md) §3。

## 文件说明

| 文件 | 说明 |
|---|---|
| `INCAR` | 铁磁静态计算（PBE，`ISPIN=2`，`MAGMOM = 2*3.0`，`LORBIT=11`） |
| `POSCAR` | bcc Fe 的 2 原子惯用立方胞，a = 2.83 Å（接近 PBE 平衡值；实验约 2.87 Å） |
| `KPOINTS` | Γ 中心 12×12×12 网格 |
| `job.slurm` | 最简 Slurm 脚本（默认 `vasp_std`） |
| `run_fm_afm_nm.sh` | 自动建立 FM / AFM / NM 三个子目录、提交并汇总能量与磁矩 |

POTCAR 不随仓库提供，请自行生成（VASP 推荐 `Fe_pv`）：

```bash
cat $VASP_PP_PATH/potpaw_PBE/Fe_pv/POTCAR > POTCAR
grep TITEL POTCAR        # should show PAW_PBE Fe_pv
```

## 运行步骤

单个铁磁计算：

```bash
sbatch job.slurm
```

三种磁态对比：

```bash
bash run_fm_afm_nm.sh            # create FM/ AFM/ NM/ and submit
# ... after the jobs finish ...
bash run_fm_afm_nm.sh collect    # print E0 and total magnetisation
```

## 预期结果 / 如何检查

```bash
grep "F=" OSZICAR | tail -1                 # last line contains "mag=" (total moment, muB per cell)
grep -A6 "magnetization (x)" OUTCAR | tail -7   # per-atom s/p/d/tot moments (LORBIT=11)
grep "energy(sigma->0)" OUTCAR | tail -1
grep "EENTRO" OUTCAR | tail -1              # T*S, should be < ~1 meV/atom
```

- FM：每个 Fe 的 `tot` 磁矩约 2.2–2.3 μB（PBE 的典型量级，示意，具体以你的计算为准；实验约 2.22 μB/Fe），OSZICAR 中 `mag=` 约为其 2 倍。
- AFM（`MAGMOM = 3 -3`）：两个原子磁矩符号相反，总磁矩≈0；bcc Fe 的 AFM 态能量明显高于 FM，且局域磁矩通常比 FM 小，甚至可能塌缩成非磁解。
- NM（`ISPIN=1`）：能量最高。三者能量顺序应为 E(FM) < E(AFM)、E(FM) < E(NM)。

> [!NOTE]
> `magnetization (x)` 表中的原子磁矩是在 PAW 球内积分得到的，球外（间隙区）的磁化不计入，所以各原子 `tot` 之和一般不严格等于 OSZICAR 中的 `mag=`。

## 注意事项

- 磁性金属一定要检查 smearing 熵项 `T*S`，比较能量时使用 `energy(sigma->0)`。
- 本例晶格常数固定为 2.83 Å；要做 FM 与 NM 的平衡体积对比，请按 [第4章](../../tutorials/04-relax-eos.md) 分别对每种磁态做 E–V 扫描。
- 若 AFM 计算"翻转"成 FM 或 NM，可尝试：更大的初始磁矩、`ISYM = 0`、或者先从非磁 CHGCAR 起步（见教程 §3.5）。
- VASP < 6.6.0 中，若 `MAGMOM` 含负值，默认 `NBANDS` 可能偏少（known issue 51）；遇到奇怪的崩溃时请手动设置 `NBANDS`。

[返回第7章](../../tutorials/07-magnetism-u-soc.md)
