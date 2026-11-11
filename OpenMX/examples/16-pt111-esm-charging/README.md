# 示例 16-2：带电 Pt(111) 电极（ESM on3/on4）与固定底层弛豫

配套教程：[第16章 表面与界面：ESM、功函数、电荷掺杂](../../tutorials/16-surface-esm-workfunction.md)

## 目的

| 输入 | 演示内容 |
|---|---|
| `Pt111_ESM_on3.dat` | ESM `on3`（真空\|板\|理想金属）：给 Pt(111) 板加 `scf.system.charge`，过剩电荷由右侧理想金属里的镜像电荷补偿，相当于一个“电极 + 对电极”的电容器 |
| `run_charge_scan.sh` | 扫描 q = −0.02…+0.02 e/胞，汇总 q、面电荷密度 σ、化学势 μ 与 Utot，用于估算微分电容 |
| `Pt111_ESM_on4.dat` | ESM `on4`：两侧都是理想金属，用 `ESM.potential.diff 1.0`（eV）在两块金属之间加 1 V 电势差，形成均匀外电场 |
| `Pt111_relax_fixbottom.dat` | ESM `on1` 下的结构优化：`<MD.Fixed.XYZ` 固定底部两层（原子 1–2），只弛豫顶部两层 |

## 文件说明

- 结构：4 层 Pt(111)，1×1 表面胞，采用实验晶格常数 3.92 Å，体相截断。a 轴沿 x，长 40 Å。板占据 x = 12.00–18.79 Å，到 x=0 的 ESM(I) 为 12 Å，到 x=40 Å 的 ESM(II) 为 21.2 Å。
- 基组与赝势：`Pt7.0-s3p2d2f1` + `Pt_PBE19`（DFT_DATA19 推荐表中的 Standard 档，16 个价电子，坐标行写 8.0 8.0）。
- 面积：$A=\frac{\sqrt3}{2}a^2=6.654$ Å$^2$，于是 σ(μC/cm²) = 1602.18 × q / A。q=0.01 时 σ≈2.41 μC/cm²。

## 运行步骤

```bash
cd OpenMX/examples/16-pt111-esm-charging
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" *.dat

# 1) charge scan (5 SCF runs; NP/NT/OPENMX can be overridden)
NP=16 NT=1 bash run_charge_scan.sh
#    only generate the inputs:   DRYRUN=1 bash run_charge_scan.sh
#    only collect the results:   COLLECT_ONLY=1 bash run_charge_scan.sh

# 2) uniform field between two ideal metals
mpirun -np 16 openmx Pt111_ESM_on4.dat -nt 1 > pt111_on4_dv1.std 2>&1

# 3) relaxation with the bottom two layers fixed
mpirun -np 16 openmx Pt111_relax_fixbottom.dat -nt 1 > pt111_relax.std 2>&1
```

## 预期结果 / 如何检查

```bash
grep "Boundary condition" pt111_on3_qp0.000.std pt111_on4_dv1.std   # Vacuum|Vacuum|Metal ; Metal|Vacuum|Metal + field
cat charge_scan.csv
# planar-averaged excess charge of q=+0.01 relative to the neutral slab
diff_gcube pt111_on3_qp0.010.tden.cube pt111_on3_qp0.000.tden.cube dq.cube   # tool from source/, see chapter 10
gcube2oned dq.cube 1 > dq_1d.txt
grep -A6 "MD.Opt" pt111_relax.out | head; tail -n 3 pt111_relax.std
```

- `charge_scan.csv`：μ 应随 q 近似线性变化。在 on3 下，板与理想金属之间是真空电容器，斜率 $\mathrm d\mu/\mathrm d\sigma$ 大致对应 $d_\mathrm{eff}/\varepsilon_0$，其中 $d_\mathrm{eff}$ 是板表面到 ESM(II) 的有效距离。只有一个点或趋势不单调时，先检查各点 SCF 是否都收敛。
- `dq_1d.txt`：过剩电荷集中在朝向理想金属（x 较大）的那一侧表面。这与手册第 49.2 节 Al–Si(111) 示例的现象一致：被掺入的电荷受镜像电荷吸引，偏聚到靠近金属的一侧。
- `pt111_on4_dv1`：用 `gcube2oned pt111_on4_dv1.v0.cube 1` 画出平面平均势，真空区应能看到线性斜坡，板内部的电场则被金属屏蔽。
- `pt111_relax`：原子 1–2 的坐标保持不变，顶层间距通常会有百分之几的弛豫。最终结构见 `pt111_relax.xyz`，如需续算，用 `pt111_relax.dat#`。

以上只是定性预期，具体数值取决于层数、k 点和基组，请以自己的计算为准。

## 注意事项

- `ESM.potential.diff` 只在 `on4` 下合法；在其他模式下设成非零值，程序会报错退出（`<ESM:Warning> Non-zero ESM.ElectricField is not valid except for on4.`）。
- `on1` 要求体系电中性，带电计算请改用 `on2`/`on3`。
- 不同电荷态的总能量不能直接相减来比较，要考虑电子化学势的参考（巨正则自由能），详见教程第 7 节。
