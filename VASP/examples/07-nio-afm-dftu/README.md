# 算例 07-nio-afm-dftu：NiO 反铁磁 PBE+U 与交换耦合常数 J1/J2

## 目的

- 用 Dudarev 方案（`LDAUTYPE = 2`）对 Ni 3d 加 U，计算 AFM-II 型 NiO 的带隙与局域磁矩；
- 理解 `LDAUL / LDAUU / LDAUJ` "每个元素一个值"的写法，以及 `LMAXMIX = 4`、`LDAUPRINT` 的作用；
- 用 FM、AFM-I、AFM-II 三种共线构型的能量差映射 Heisenberg 模型，得到近邻 J1 与次近邻 J2，并用平均场粗估 Néel 温度。

对应教程：[第7章](../../tutorials/07-magnetism-u-soc.md) §5（DFT+U）与 §7（交换耦合 J）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `INCAR` | AFM-II 静态 PBE+U（U_eff = 6.2 eV，`MAGMOM = 2 -2 0 0`） |
| `POSCAR` | AFM-II 菱方胞（2 个 NiO 化学式单元，自旋沿 [111] 方向逐层交替），a(立方) = 4.17 Å（实验值） |
| `POSCAR_AFM1` | AFM-I 四方胞（2 个化学式单元，自旋沿 [001] 逐层交替） |
| `KPOINTS` / `KPOINTS_AFM1` | 两种晶胞的 Γ 中心网格，k 点间距相近（约 0.25 Å⁻¹，含 2π） |
| `run_exchange.sh` | 建立 FM_R / AFM2 / AFM1 / FM_T 四个子目录并提交；`collect` 子命令汇总 |
| `exchange_J.py` | 由三种构型能量解出 J1、J2，给出平均场 T_N 估计 |
| `job.slurm` | 最简 Slurm 脚本 |

POTCAR（VASP 推荐 `Ni_pv`；顺序必须与 POSCAR 第 6 行 `Ni O` 一致）：

```bash
cat $VASP_PP_PATH/potpaw_PBE/{Ni_pv,O}/POTCAR > POTCAR
grep TITEL POTCAR
```

## 运行步骤

1. 单个 AFM-II 计算：

   ```bash
   sbatch job.slurm
   ```

2. 交换耦合常数（4 个计算，耗时约为单个计算的 4–6 倍）：

   ```bash
   bash run_exchange.sh            # prepare FM_R AFM2 AFM1 FM_T and submit
   # ... wait until all four jobs have finished ...
   bash run_exchange.sh collect    # = python exchange_J.py --fm ... --afm1 ... --afm2 ...
   ```

## 预期结果 / 如何检查

```bash
grep -A8 "magnetization (x)" OUTCAR | tail -9     # Ni: about +-1.7 muB, O: ~0
grep "F=" OSZICAR | tail -1                      # mag= should be ~0 for AFM-II
grep -B2 -A12 "occupancies and eigenvectors" OUTCAR | tail -30   # LDAUPRINT=1 output
```

- AFM-II 中两个 Ni 磁矩大小相同、符号相反，总磁矩≈0。PBE+U（U_eff≈6 eV）给出的 Ni 磁矩通常在 1.6–1.8 μB 量级（示意，具体以你的计算为准；实验约 1.6–1.9 μB）。
- 带隙：PBE 只有不到 1 eV，PBE+U(6.2) 一般在 3 eV 左右（示意；实验光学带隙约 4 eV）。可用 `BANDGAP` 标签（VASP ≥ 6.4.3）或 [第5章](../../tutorials/05-band-dos.md) 的 DOS 方法读取。
- `exchange_J.py` 的输出中，AFM-II 能量应最低，J2 为负（强反铁磁超交换，主导项），|J1| 远小于 |J2|。若 `FM (other cell) - FM` 超过约 1 meV/Ni，说明两种晶胞的 k 点或 ENCUT 不够收敛，J 的数值不可信。

## 注意事项

- **不同 U 值之间的总能量不能直接比较**；本例四个计算必须使用完全相同的 INCAR（只改 MAGMOM）。
- 6.2 eV 是 Materials Project（pymatgen `MPRelaxSet`）对氧化物中 Ni 使用的经验值，并非"正确答案"；需要第一性原理的 U 时，见算例 [07-linear-response-u](../07-linear-response-u/README.md)。
- 若要用 `ICHARG = 11` 画 DFT+U 能带，自洽步骤必须已经设置 `LMAXMIX = 4`，否则非自洽能带与自洽结果不一致。
- DFT+U 容易收敛到亚稳态（不同的 d 轨道占据）。比较不同初始 `MAGMOM` 或从 U=0 的波函数逐步增大 U（U ramping），取能量最低的解，并用 `LDAUPRINT = 1` 输出的占据矩阵核对。
- 平均场 T_N 只是数量级估计（通常偏高）；定量的居里/奈尔温度请用经典蒙特卡洛（例如 VAMPIRE 等程序）。

[返回第7章](../../tutorials/07-magnetism-u-soc.md)
