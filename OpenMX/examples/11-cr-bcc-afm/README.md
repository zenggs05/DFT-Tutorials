# 算例 11-cr-bcc-afm：bcc Cr 的反铁磁 / 铁磁 / 非磁对比

## 目的

- 学会用相反的初始占据（第 6、7 列互换）设定反铁磁（AFM）初态；
- 在完全相同的数值参数下比较 AFM、FM、NM 三种状态的能量与磁矩；
- 理解"初始磁构型决定落入哪个自洽解"。

对应教程：[第11章 自旋极化、非共线磁性与交换耦合](../../tutorials/11-spin-noncollinear.md) §4。

## 文件说明

| 文件 | 说明 |
|---|---|
| `Cr_bcc_AFM.dat` | 2 原子简单立方胞（a = 2.88 Å），Cr6.0-s3p2d1 / Cr_PBE19（14 价电子），原子 1 为 8.5/5.5、原子 2 为 5.5/8.5 |
| `run_fm_afm_nm.sh` | 用 sed 生成 `Cr_bcc_FM.dat`（两原子都 8.5/5.5）与 `Cr_bcc_NM.dat`（自旋关闭，7.0/7.0），依次运行并汇总；可直接 `sbatch` |

运行前修改 `DATA.PATH`。

## 运行步骤

```bash
sbatch run_fm_afm_nm.sh            # or: bash run_fm_afm_nm.sh
bash run_fm_afm_nm.sh collect      # print the summary again later
```

## 预期结果 / 如何检查

```bash
grep "Utot\." cr_afm.out cr_fm.out cr_nm.out
grep -A 5 "Mulliken populations" cr_afm.out
grep -A 8 "Voronoi charges" cr_afm.out
```

- AFM：两个 Cr 的 `Diff`（局域磁矩）大小相近、符号相反，`Total spin moment` ≈ 0；
- 能量顺序预期 E(AFM) < E(NM)；FM 初态往往不稳定，可能塌缩为非磁或很小的磁矩；
- GGA 一般高估 Cr 的局域磁矩（实验的自旋密度波振幅约 0.6 μB 量级），这是泛函的已知问题；
- 在 VESTA 中查看 `cr_afm.sden.cube`，两个 Cr 周围应是符号相反的自旋密度。

具体数值以你的计算为准。

## 注意事项

- 如果 AFM 计算翻转成 FM 或 NM：减小 `scf.Init.Mixing.Weight` / `scf.Max.Mixing.Weight`，或先在较高电子温度下收敛，再 `scf.restart on` 降温续算。
- 三种状态必须使用相同的晶胞、基组、截断能、k 点与电子温度，能量才可比较。
- 另一种写法是为同一元素定义两个 species（如 `Cr1`、`Cr2`），便于按子格统计 PDOS 或给不同子格不同的 U，见教程 §4.1。

[返回第11章](../../tutorials/11-spin-noncollinear.md)
