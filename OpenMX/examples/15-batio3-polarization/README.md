# 示例 15-2：四方 BaTiO3 的自发极化（中心对称参考 vs 极化结构）

配套教程：[第15章 Berry 相宏观极化与外加电场](../../tutorials/15-berry-efield.md)

## 目的

自发极化是一个**差值**：$\Delta\mathbf P=\mathbf P(\lambda=1)-\mathbf P(\lambda=0)$，其中 λ = 0 是中心对称参考结构，λ = 1 是铁电结构。本例在**同一个四方晶胞**里把离子从中心对称位置（P4/mmm）线性移动到极化位置（P4mm），每个 λ 用 polB 计算 $P_z$，再把各点连成连续曲线，从而确定正确的分支，避开极化量子 $eR/\Omega$ 带来的多值问题。

## 文件说明

| 文件 | 说明 |
|---|---|
| `BaTiO3_cubic_ref.dat` | λ = 0：四方晶胞（a = 3.992 Å，c = 4.036 Å）中的中心对称位置，`HS.fileout on`；之所以不用真正的立方胞，是为了让路径上只有离子位移，没有晶胞变化 |
| `BaTiO3_tetra_relax.dat` | 固定晶胞、固定 Ba，从一个小的偏心初猜出发用 `MD.Type RF` 弛豫内坐标 |
| `BaTiO3_tetra_polB.dat` | λ = 1：极化结构（初猜坐标），`HS.fileout on`；如果已经跑过弛豫，`make_path.py --from-out` 会改用弛豫后的坐标 |
| `make_path.py` | 读两个端点（FRAC 坐标），按最小像插值，生成 `BaTiO3_lam0.00.dat` … `BaTiO3_lam1.00.dat` |
| `polB.in` | `6 6 8` / `0 0 1`：只算沿 c（极化轴）的 Berry 相，c 方向的 k 点取得更密 |
| `polar_path.py` | 解析各 λ 的 polB 输出，按极化量子展开分支，输出 ΔP（μC/cm²） |
| `job_slurm.sh` | `STAGE=relax`：弛豫；`STAGE=path`：生成路径 → SCF → polB → 分析 |

## 运行步骤

```bash
cd OpenMX/examples/15-batio3-polarization
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" BaTiO3_*.dat

# 1) 弛豫极化结构的内坐标（可跳过，跳过则直接用初猜坐标）
sbatch --export=ALL,STAGE=relax job_slurm.sh
grep -A8 "Fractional coordinates of the final structure" bto_relax.out

# 2) lambda 路径 + polB
sbatch --export=ALL,STAGE=path,POLB=$HOME/openmx4.0/work/polB job_slurm.sh

# 手动分析（job 脚本最后会自动执行）
python3 polar_path.py 0.00:polB_bto_lam0.00.out 0.25:polB_bto_lam0.25.out \
    0.50:polB_bto_lam0.50.out 0.75:polB_bto_lam0.75.out 1.00:polB_bto_lam1.00.out
```

## 预期结果 / 如何检查

- 弛豫后 Ti 和 O 应沿 c 方向相对 Ba 反向位移，`bto_relax.out` 末尾的 Total energy 应低于 λ = 0 结构的能量（`grep "Utot" bto_relax.out bto_lam0.00.out`）。如果弛豫回到了中心对称位置（Ti 的 z 回到 0.5），说明在当前泛函和晶胞下没有铁电不稳定性，此时 ΔP ≈ 0 才是正确结果。
- λ = 0 的 $P_z$ 必须是 f = 1 量子（≈100.5 μC/cm²）的整数倍：本例自旋简并、价电荷全为偶数，真正的量子是 2ec/Ω，中心对称结构只能取它的 0 或 1/2（教程 5.4 节）；对称性正确的话，`polar_path.py` 的 unwrapped 列从 λ = 0 开始应随 λ 平滑、近似线性地增大。
- 极化量子：沿 c 为 $e c/\Omega=e/a^2\approx100.5$ μC/cm²（f = 1），应与 polB 输出末尾 `eRiz/Omega` 行 R3 列一致（`grep -A5 arbitrariness polB_bto_lam0.00.out`）。任何相邻两点之差都应远小于这个值。
- `grep 'may be metallic' polB_bto_lam*.out` 应没有输出；若有，说明该点被 polB 判为金属，结果无效。
- 参考：四方 BaTiO3 室温实验自发极化约 26 μC/cm²（Wieder, 1955）。LDA 在实验晶胞下的计算值通常与之同量级；具体数值取决于弛豫后的位移，**不要把初猜坐标算出的数当作结果**。

## 注意事项

- 端点必须使用完全相同的晶胞和原子顺序，`make_path.py` 会检查原子顺序。
- `scf.Kgrid 8 8 8` 和 `polB.in` 中的 `6 6 8` 都需要做收敛测试（ΔP 变化 < 0.5 μC/cm²）。
- 只算 c 方向时，polB 输出中 x、y 分量的电子部分没有意义，只读 Pz。
- 换成 PbTiO3：把 Ba 换成 Pb（赝势与基组以 DFT_DATA19 中实际存在的 `Pb_CA19.vps` 和推荐基组表为准），晶格换成 PbTiO3 的四方晶格常数，其余流程不变。PbTiO3 的四方畸变和位移都明显更大，路径上建议多取几个 λ 点。
- 本目录的输入未经作者实际运行，文中数值均为参考值或判据。

返回教程：[第15章](../../tutorials/15-berry-efield.md)
