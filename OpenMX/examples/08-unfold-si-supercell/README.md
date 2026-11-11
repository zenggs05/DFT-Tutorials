# 算例 8-3：Si 超胞的能带展开（Band Unfolding）

## 目的

用 8 原子的 Si 惯用立方胞做 SCF，并用 OpenMX 的能带展开功能（`Unfolding.Electronic.Band`）把超胞能带"展开"回 fcc 原胞的布里渊区：

1. 完美 Si：展开后的谱权重应完全重现原胞能带——这是检查 `<Unfolding.ReferenceVectors`、`<Unfolding.Map` 和 k 点写法是否正确的最好方法；
2. Si₇Ge₁（一个 Si 被 Ge 替代）：观察杂质引起的谱权重减弱与展宽。

对应教程：[第8章 能带结构（含自动 k 路径与能带展开）](../../tutorials/08-band.md)。

## 文件说明

| 文件 | 作用 |
|---|---|
| `Si_supercell_unfold.dat` | 完美 Si 的 8 原子立方胞（a = 5.43 Å），参考胞为 fcc 原胞；原子 1–4、5–8 分属两套 fcc 子格（组 1、2）；路径 L–Γ–X–W–K–Γ；`System.Name si8_unfold` |
| `Si7Ge1_unfold.dat` | 原子 1 换成 Ge（`Ge7.0-s3p2d2` + `Ge_PBE19`），映射中给 Ge 单独的组号 3；结构未弛豫，仅作演示；`System.Name si7ge1_unfold` |
| `plot_unfold.py` | 读取 `*.unfold_totup`（或 `_totdn`、非共线的 `_tot`）与 `*.unfold_plotexample` 中的高对称点位置，画谱权重散点图，可叠加原胞的 `BANDDAT1`，也可按轨道列求和 |
| `job.slurm` | 依次运行两个输入并作图（若 `../08-band-cdia-si/si.BANDDAT1` 存在则自动叠加原胞能带） |

## 运行步骤

```bash
cd OpenMX/examples/08-unfold-si-supercell
# edit DATA.PATH in both inputs first
mpirun -np 8 openmx Si_supercell_unfold.dat -nt 1 | tee si8_unfold.std
mpirun -np 8 openmx Si7Ge1_unfold.dat       -nt 1 | tee si7ge1_unfold.std
ls si8_unfold.unfold_*            # unfold_totup, unfold_orbup, unfold_plotexample

# quick look with gnuplot
gnuplot -p -e "load 'si8_unfold.unfold_plotexample'"

# primitive-cell reference bands with the same path (from ../08-band-cdia-si)
python3 plot_unfold.py si8_unfold    --band ../08-band-cdia-si/si.BANDDAT1 --out si8_unfold.png
python3 plot_unfold.py si7ge1_unfold --band ../08-band-cdia-si/si.BANDDAT1 --out si7ge1_unfold.png

# orbital-resolved weight: columns = (orbital serial number in si7ge1_unfold.out) + 2
grep -n -i "orbital weights" -A 30 si7ge1_unfold.out | head -40
python3 plot_unfold.py si7ge1_unfold --orb-cols 3,4,5,6 --out ge_orbitals.png
```

可选：强度图（intensity map）

```bash
cd /path/to/openmx4.0/source && gcc intensity_map.c -lm -o intensity_map
./intensity_map si7ge1_unfold.unfold_totup -c 3 -k 0.1 -e 0.1 -l -13 -u 6 > si7ge1_intmap.txt
```

## 预期结果 / 如何检查

- 完美 Si：每个散点的权重接近 0 或 1，权重大的点精确落在原胞能带（灰线）上；被折叠进来的"假"能带权重接近 0。若点与灰线整体错位，先检查两者是否用了相同路径（本例都是 L–Γ–X–W–K–Γ）以及 `Unfolding.kpoint` 是否按参考胞倒格矢书写。
- Si₇Ge₁：大部分能带仍清晰（Ge 与 Si 同族、电子结构相近），但在部分能量和 k 点处权重变小、出现弱的"影子"带，费米能附近的带隙略有变化（数值依赖于未弛豫的结构，仅作定性观察）。
- `plot_unfold.py` 打印（示意输出，数值仅供参考）：

  ```text
  si8_unfold.unfold_totup: xxxx points, total weight range 0.000 - 1.0xx
  figure written to si8_unfold.png
  ```

## 注意事项

- `<Unfolding.ReferenceVectors` 的单位与 `Atoms.UnitVectors.Unit` 相同；不写时默认用 `Atoms.UnitVectors`（即不展开，只做轨道分解分析）。
- `<Unfolding.kpoint` 每行格式为 `标签 k1 k2 k3`，以**参考胞**的倒格矢为单位；`Unfolding.Nkpoint` 等于行数；路径按列出顺序依次连接。
- `<Unfolding.Map` 第一列必须按 1, 2, 3… 顺序列出全部原子；第二列是组号（正整数，不必从 1 开始）。同一组原子在参考胞中对应同一个位置；杂质原子给新的组号。
- 能带展开只支持 `Band` 和 `Cluster` 求解器；共线、非共线（含 SOC）与 +U 都支持，无需额外关键词。4.0 改进了展开代码的内存与速度。
- 输出能量以费米能为零点（eV），横轴为路径长度（Bohr⁻¹）。

返回教程：[第8章 能带结构（含自动 k 路径与能带展开）](../../tutorials/08-band.md)
