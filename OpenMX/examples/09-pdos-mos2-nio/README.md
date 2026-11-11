# 算例 9-2：单层 MoS₂ 与反铁磁 NiO 的投影态密度

## 目的

- **单层 2H-MoS₂**：二维体系的 DOS/PDOS（k 网格沿真空方向取 1），分析三棱柱配位下 Mo 4d 的 $a_1'$（d3z²−r²）、$e'$（dx²−y²、dxy）、$e''$（dxz、dyz）成分，以及 S 3p 与 Mo 4d 的杂化；体会半芯态在 `Dos.Erange` 之外对积分 DOS 的影响。
- **AFM-II NiO**：用坐标行上/下初始电子数设定反铁磁初态（Ni1 +2 μB、Ni2 −2 μB），得到自旋分辨 PDOS，区分 Ni 的 $e_g$（d1 + d2）与 $t_{2g}$（d3 + d4 + d5），观察 GGA 对带隙的低估（与第 13 章 DFT+U 对比）。

对应教程：[第9章 态密度 DOS/PDOS、COHP 与光学性质](../../tutorials/09-dos-pdos-cohp-optics.md)。

## 文件说明

| 文件 | 作用 |
|---|---|
| `MoS2_monolayer_pdos.dat` | a = 3.16 Å（近似实验值，未弛豫）、c = 20 Å；`Mo7.0-s3p2d2` / `S7.0-s2p2d1f1`；`scf.Kgrid 12 12 1`，`Dos.Kgrid 30 30 1`，`Dos.Erange -20 10`；`System.Name mos2` |
| `NiO_AFM_pdos.dat` | 岩盐 NiO（a = 4.17 Å）的 AFM-II 菱方胞（2 Ni + 2 O），`Ni6.0S-s3p2d1` + `Ni_PBE19S`（16 价电子），`scf.SpinPolarization on`，`Dos.Kgrid 14 14 14`；`System.Name nio_afm` |
| `job.slurm` | 依次运行两个体系，调用 [../09-dos-pdos-cdia/run_dosmain.sh](../09-dos-pdos-cdia/run_dosmain.sh) 和 [../09-dos-pdos-cdia/sum_pdos.py](../09-dos-pdos-cdia/sum_pdos.py) |

## 运行步骤

```bash
cd OpenMX/examples/09-pdos-mos2-nio
# edit DATA.PATH in both inputs first
mpirun -np 8  openmx MoS2_monolayer_pdos.dat -nt 1 | tee mos2.std
mpirun -np 16 openmx NiO_AFM_pdos.dat        -nt 1 | tee nio_afm.std

bash ../09-dos-pdos-cdia/run_dosmain.sh mos2    tetra 0.1
bash ../09-dos-pdos-cdia/run_dosmain.sh nio_afm tetra 0.1

python3 ../09-dos-pdos-cdia/sum_pdos.py mos2    --dat MoS2_monolayer_pdos.dat --resolve-m --plot --emin -8 --emax 5
python3 ../09-dos-pdos-cdia/sum_pdos.py nio_afm --dat NiO_AFM_pdos.dat --by atom --resolve-m --plot --emin -9 --emax 5

grep -i "spin moment" nio_afm.out | tail -n 3      # total moment, should be close to 0 for AFM
```

## 预期结果 / 如何检查

- MoS₂：费米能两侧都以 Mo d 为主，价带顶附近 d3z²−r² 与 (dx²−y², dxy) 贡献明显，S p 与之杂化；带隙为 PBE 量级（文献 PBE 约 1.6–1.7 eV，K 点直接带隙；本例未弛豫，仅作参考）。0 eV 处的积分 DOS 会明显小于 26（Mo 14 + 2×6），因为 Mo 4s/4p 半芯态在 −20 eV 以下的窗口之外。
- NiO：两个 Ni 的 PDOS 上下自旋互为镜像；Ni1 的 $e_g$ 上自旋占据、下自旋空（Ni2 相反），$t_{2g}$ 两个自旋均占据；价带顶 O 2p 与 Ni 3d 强烈杂化。GGA 带隙很小（文献 PBE 约 0.5–1 eV 量级，实验约 4 eV），Ni 磁矩约 1.2–1.4 μB（实验约 1.6–1.9 μB）——这正是需要 DFT+U 的原因。
- `sum_pdos.py` 对 NiO 的输出（示意输出，数值仅供参考）：

  ```text
  method: Tetrahedron, xx files, spin-polarized: True
  Ni1 d                        -> nio_afm.pdos_Ni1_d.dat             N(E<E_F) = x.xxx + x.xxx
  Ni2 d                        -> nio_afm.pdos_Ni2_d.dat             N(E<E_F) = x.xxx + x.xxx
  O3 p                         -> nio_afm.pdos_O3_p.dat              N(E<E_F) = x.xxx + x.xxx
  ```

## 注意事项

- NiO 的晶格矢量写在立方坐标系中，因此 PDOS 的 d1 + d2 对应 $e_g$、d3 + d4 + d5 对应 $t_{2g}$；如果旋转了晶胞，这种对应关系就不成立。
- AFM 初态的两列之和必须等于 VPS 价电子数（`Ni_PBE19S` 为 16，O 为 6）。同一元素也可以定义两个物种（如 `Ni1`、`Ni2`）来区分，配合 `sum_pdos.py --by species`。
- 磁性体系 SCF 较难收敛：本例用 `RMM-DIISK`、较小混合权重和较长历史；若收敛到铁磁或非磁解，检查 `.out` 中各 Ni 的 Mulliken 磁矩。OpenMX 自旋磁矩的符号约定见第 11 章。
- 二维体系 DOS 是"每个超胞"的值，k 网格沿 c 方向只取 1。

返回教程：[第9章 态密度 DOS/PDOS、COHP 与光学性质](../../tutorials/09-dos-pdos-cohp-optics.md)
