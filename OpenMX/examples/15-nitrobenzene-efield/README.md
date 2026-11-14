# 示例 15-3：硝基苯在外加均匀电场中的响应（scf.Electric.Field + diff_gcube）

配套教程：[第15章 Berry 相宏观极化与外加电场](../../tutorials/15-berry-efield.md)

## 目的

对孤立分子施加锯齿势表示的均匀外电场（`scf.Electric.Field`，单位 GV/m），完成三件事：

1. 读出不同场强下的偶极矩，用有限差分求静态极化率 $\alpha_{xx}=\partial\mu_x/\partial E_x$；
2. 用 `diff_gcube` 计算场致电子密度变化 $\Delta\rho=\rho(E)-\rho(0)$，重现 OpenMX 手册第 30 章图 33 的定性图像（−NO2 的 O 原子与对位 C、H 之间的大幅电荷转移）；
3. 用偶极变化的方向实际确认 `scf.Electric.Field` 的符号约定。

## 文件说明

| 文件 | 说明 |
|---|---|
| `Nitro_Benzene_EF.dat` | 零场模板：理想化平面几何，分子 C2 轴沿 a（= x），盒子 24 × 20 × 14 Å，GGA-PBE，`cluster` 求解器，`level.of.fileout 1` |
| `run_fields.sh` | 由模板生成 `Nitro_Benzene_EF_{m10,m2,p2,p10}.dat`（场强 −10、−2、+2、+10 GV/m，沿 a）；`RUN=1` 时依次运行 |
| `dipole_alpha.py` | 读各 `.out` 的 `Dipole moment (Debye)` 块，只用 \|F\| ≤ 3 GV/m 的点做线性拟合求 α（a.u. 与 Å³），并做符号检查 |
| `job_slurm.sh` | Slurm 作业：全部场强 → `dipole_alpha.py` → `diff_gcube` 生成 `drho_p10.cube`、`drho_m10.cube` |

## 运行步骤

```bash
cd OpenMX/examples/15-nitrobenzene-efield
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" Nitro_Benzene_EF.dat

# diff_gcube is a single C file; the manual (sec. 66) says it is in work/
gcc $(find ~/openmx4.0 -name diff_gcube.c | head -1) -lm -o diff_gcube

sbatch --export=ALL,DIFF_GCUBE=$PWD/diff_gcube job_slurm.sh

# 或在工作站上：
RUN=1 NP=8 bash run_fields.sh
python3 dipole_alpha.py 0:nben_f0.out -2:nben_m2.out 2:nben_p2.out -10:nben_m10.out 10:nben_p10.out
./diff_gcube nben_p10.tden.cube nben_f0.tden.cube drho_p10.cube
```

## 预期结果 / 如何检查

```bash
grep -A2 "Applied External Electric Field" Nitro_Benzene_EF_p10.std   # direction (1,0,0), magnitude 10
grep -A12 "Dipole moment (Debye)" nben_f0.out | grep -E "Absolute|Total"
grep "Utot" nben_*.out
```

- 加场的计算在标准输出开头会打印 `<Applied External Electric Field>`，以及方向（单位矢量）和大小（10^9 V/m）。没有这一段说明关键词没有生效（拼写错误或三个分量都是 0）。
- 零场时偶极沿分子轴（x），硝基苯的偶极较大（实验约 4 D 量级）；+2 与 −2 GV/m 下的 Dx 应关于零场值近似对称变化。
- `dipole_alpha.py` 给出的 α_xx 必须为正；若为负，说明在这一套约定下，正的关键词值对应 −x 方向的物理电场，脚本会给出提示。±10 GV/m 的点通常已明显偏离线性，脚本把它们单独列出。
- 在 VESTA 中打开 `drho_p10.cube` 和 `drho_m10.cube`，两者的增减区域应大致互换；电荷主要在 O 原子与对位 C–H 之间转移。
- 本目录的输入未经作者实际运行，以上都是定性判据。

## 注意事项

- 锯齿势在晶胞中有一个不连续面。分子放在盒子中央，沿场方向每侧留出约 8 Å 以上的真空，让不连续面落在真空中。可以用 `gcube2oned`（第16章）对 `v0.cube` 做沿 a 的平面平均来检查跳变位置；若加场与零场的 `v0` 之差是一条线性斜坡，就说明外场项包含在 `v0` 中。
- 三个分量分别对应 a、b、c 晶格矢量方向（见 4.0 源码 `Input_std.c` 的处理）；本例为正交盒子，所以等同于 x、y、z。
- 场强不宜过大：10 GV/m = 1 V/Å，已经足以让 HOMO/LUMO 能级明显移动；求线性响应用 1–3 GV/m。
- 几何是理想化的，定量比较前请先在零场下弛豫（第6章），并保持所有场强使用同一几何。
- `diff_gcube` 要求两个 cube 的网格完全一致：只改场强、不改结构和晶胞时，这一条件自动满足。若屏幕上出现 `Found a difference in ...`，说明网格不一致；程序不会停下，但写出的差分无效。
- `dipole_alpha.py` 的参数可以以负号开头（如 `-2:nben_m2.out`），脚本已做处理。

返回教程：[第15章](../../tutorials/15-berry-efield.md)
