# 算例 10-co-cu111-cdd：CO/Cu(111) 吸附的电荷密度差 ρ(AB) − ρ(A) − ρ(B)

## 目的

- 掌握 OpenMX 计算吸附体系电荷密度差的标准流程（手册第 68 章）：AB 复合体系 → 读取 `Grid_Origin` → 片段 A、B 用 `scf.fixed.grid` 对齐网格 → `add_gcube` / `diff_gcube`；
- 用平面平均与累计积分估计 CO 与 Cu 表面之间的电荷转移；
- 学会检查三个 cube 文件的网格是否完全一致。

对应教程：[第10章 电荷分析、分子轨道与 cube 可视化](../../tutorials/10-charge-mo-cube.md) §7。

## 文件说明

| 文件 | 说明 |
|---|---|
| `CO_on_Cu111_AB.dat` | 复合体系：Cu(111) 2×2 表面晶胞、3 层（12 个 Cu）+ 顶位竖直吸附的 CO |
| `CO_on_Cu111_A.dat` | 片段 A：只有 CO，晶胞与几何和 AB 完全相同；末尾留有 `#scf.fixed.grid` 占位行 |
| `CO_on_Cu111_B.dat` | 片段 B：只有 Cu 板，同上 |
| `run_cdd.sh` | Slurm 作业脚本：依次计算 AB、A、B，自动填入 `scf.fixed.grid`，检查网格，做差分与平面平均 |
| `cdd_planar.py` | 沿 c 轴的平面平均 Δλ(z) 与累计电荷 ΔQ(z)（numpy；画图需 matplotlib） |

模型参数：Cu 晶格常数约 3.61 Å（实验值），Cu–C 1.90 Å、C–O 1.15 Å 均为**未优化的近似值**，仅用于演示流程；基组为 DFT_DATA19 推荐表中的 Quick（Cu6.0S-s2p1d1 / Cu_PBE19S，11 个价电子）与 Standard（C、O）。

## 运行步骤

1. 修改三个 `.dat` 中的 `DATA.PATH`；
2. 编译 OpenMX 自带的两个小工具（只需一次）：

   ```bash
   cd ~/openmx4.0/source
   gcc add_gcube.c  -lm -o add_gcube
   gcc diff_gcube.c -lm -o diff_gcube
   cp add_gcube diff_gcube ~/bin/          # any directory in PATH
   ```

3. 提交作业（16 核约数十分钟量级，视机器而定）：

   ```bash
   sbatch run_cdd.sh
   ```

   手动执行时，关键两步是：

   ```bash
   grep "Grid_Origin" co_cu_AB.std          # x y z in Bohr
   # put "scf.fixed.grid  x y z" into the A and B inputs, then run them
   ```

## 预期结果 / 如何检查

```bash
grep "scf.fixed.grid" CO_on_Cu111_A_run.dat CO_on_Cu111_B_run.dat
for s in AB A B; do awk 'NR==3{print $2,$3,$4} NR>=4&&NR<=6' co_cu_${s}.tden.cube; echo; done
tail -n 3 cdd_z.dat
```

- 三个 tden.cube 的原点与网格矢量必须逐字相同（`run_cdd.sh` 会自动比较）；
- `add_gcube`/`diff_gcube` 打印 `Found a difference in the number of atoms` 是正常的（片段原子数不同）；若出现 `number of grid` 或 `coordinate of the origin` 的差异提示，结果无效；
- `cdd_z.dat` 最后一行给出的全晶胞积分应接近 0（电荷守恒）；
- 用 `python cdd_planar.py co_cu_dAB.cube --zref 13.1` 读取 Cu 顶层（z ≈ 12.17 Å）与 C（z ≈ 14.07 Å）之间平面上的 ΔQ，作为跨越界面的净电子转移估计；
- 在 VESTA 中打开 `co_cu_dAB.cube`，等值面取 ±0.001 ～ ±0.003 e/Bohr³。按 Blyholder 图像，常见特征是 CO 轴向（σ）区域亏损、π 区域积累、表面 Cu 附近电荷重新分布；以你的计算结果为准。

## 注意事项

- A、B 的几何必须与 AB 中逐位相同（不要单独优化片段），三个计算的晶胞、`scf.energycutoff`、`scf.Kgrid`、电子温度都必须相同。
- 这个模型没有偶极修正；只看密度差无妨，若要算功函数或吸附能随覆盖度变化，请参考第16章的 ESM 方法。
- 若要做自旋极化体系，可同样对 `sden.cube` 求差。
- 3 层板太薄，只适合演示；研究工作中至少 4–5 层并固定底层、优化顶层与吸附物。

[返回第10章](../../tutorials/10-charge-mo-cube.md)
