# 示例 15-1：NaCl 中 Na 的 Born 有效电荷（Berry 相 polB + 有限位移）

配套教程：[第15章 Berry 相宏观极化与外加电场](../../tutorials/15-berry-efield.md)

## 目的

用 `HS.fileout on` 跑 SCF 生成 `System.Name.scfout`，再用 OpenMX 自带的后处理程序 `polB` 按 Berry 相方法计算宏观电偶极（Debye）和极化（μC/cm²）。把 Na 沿 x 方向移动 −0.05、0、+0.05 Å，对极化做中心差分，得到 Born 有效电荷 $Z^*_{xx}(\mathrm{Na})=\dfrac{\Omega}{|e|}\dfrac{\Delta P_x}{\Delta u_x}=\dfrac{\Delta D_x}{|e|\,\Delta u_x}$。

## 文件说明

| 文件 | 说明 |
|---|---|
| `NaCl_polB.dat` | 2 原子 fcc 原胞，LDA（`*_CA19`），a0 = 5.64 Å，`HS.fileout on`，Na 在原点（未位移的模板） |
| `make_disp.sh` | 从模板生成 `NaCl_dx_m0.05.dat`、`NaCl_dx_p0.00.dat`、`NaCl_dx_p0.05.dat`（只改 Na 的 x 坐标和 `System.Name`） |
| `polB.in` | polB 的批处理输入：第 1 行是倒格矢 a、b、c 方向的 k 网格分割数，第 2 行是要计算的方向开关（1 = 计算） |
| `born_charge.py` | 读取多个 polB 输出，线性拟合 $D_\alpha(u)$，给出 $Z^*_{\alpha x}$ 及其芯电荷/电子分解，并用极化块交叉核对 |
| `job_slurm.sh` | Slurm 作业：生成位移输入 → SCF → polB → `born_charge.py` |

## 运行步骤

```bash
cd OpenMX/examples/15-nacl-polb-born
# 1) DATA.PATH 改成你的 DFT_DATA19 绝对路径
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" NaCl_polB.dat

# 2) 编译 polB（与 openmx 用同一个 makefile；make all 也会编译它）
cd ~/openmx4.0/source && make polB && cd -
#    手册说可执行文件生成在 work 目录；找不到时：find ~/openmx4.0 -name polB -type f

# 3) 集群：
sbatch --export=ALL,POLB=$HOME/openmx4.0/work/polB job_slurm.sh

# 3') 工作站手动执行：
bash make_disp.sh
for u in m0.05 p0.00 p0.05; do
  mpirun -np 8 openmx NaCl_dx_$u.dat -nt 1 > NaCl_dx_$u.std 2>&1
  mpirun -np 8 polB nacl_dx_$u.scfout < polB.in > polB_NaCl_dx_$u.out
done
python3 born_charge.py --atom Na --dir x \
  -0.05:polB_NaCl_dx_m0.05.out 0.00:polB_NaCl_dx_p0.00.out 0.05:polB_NaCl_dx_p0.05.out
```

## 预期结果 / 如何检查

```bash
grep "The calculation was normally finished." NaCl_dx_*.std      # 3 hits
ls nacl_dx_*.scfout                                               # 3 files
grep -A6 "Electric dipole" polB_NaCl_dx_p0.05.out                 # Dx/Dy/Dz rows
grep -A6 "Electric polarization" polB_NaCl_dx_p0.05.out           # Px/Py/Pz rows
grep "Cell_Volume" polB_NaCl_dx_p0.00.out
grep -A5 "arbitrariness of eRi/Omega" polB_NaCl_dx_p0.00.out  # polarization quanta
grep "may be metallic" polB_NaCl_dx_*.out                        # must print nothing
```

- 未位移的 `p0.00` 结构是立方对称的：Total 列的 Dx、Dy、Dz 应相等（它们本身的绝对值没有物理意义，见教程 2.3 节）。
- `born_charge.py` 输出的 Z*_xx 中，芯电荷部分应恰好等于 Na_CA19 的价电子数 9（包含 2s2p 半芯态），电子部分为负值、抵消掉大部分。
- 参考值：OpenMX 4.0 手册第 44 章用官方 `work/NaCl.dat`（由手册输出的晶格矢量可知 a0 = 5.63 Å）得到 Z* = 1.05，表中列出的另一计算值为 1.09，实验值为 1.12。本示例的基组和晶格常数与手册不同，结果应落在 1.0–1.1 附近；偏差明显更大时，先检查 SCF 是否收敛、polB 的 k 网格是否足够。
- 非对角元 Z*_yx、Z*_zx 应接近 0。若 `born_charge.py` 打印 `WARNING ... one-sided slopes differ`，说明三个点不在一条直线上：要么位移太大，要么某个点的 Berry 相跳到了另一个分支。

## 注意事项

- 块（`<...>` 与 `...>`）中的行不能加注释；`make_disp.sh` 靠行首的 `1  Na` 定位要移动的原子。
- polB 的结果只对绝缘体有意义；NaCl 带隙很大，`scf.ElectronicTemperature 300` 不会造成分数占据。
- `polB.in` 必须恰好两行，第一行三个数都 ≥ 1；若第一行含 0，批处理的 polB 会不停打印 `The number of grids should be one or more`（见教程 3.2 节）。
- `born_charge.py` 的参数可以以负号开头（如 `-0.05:file`），脚本已做处理，无需加 `--`。
- 收敛测试：把 `polB.in` 第 1 行改为 `6 6 6`、`12 12 12`，Z* 变化应小于 0.01；再用 `DISPS="-0.02 0.00 0.02"` 复算，检查线性。
- 声学求和规则：再对 Cl 做同样的位移计算（修改 `make_disp.sh` 中匹配的原子），应有 $Z^*(\mathrm{Na})+Z^*(\mathrm{Cl})\approx 0$。
- 本目录的输入未经作者实际运行，文中所有数值都是参考值或判据，请以你自己的输出为准。

返回教程：[第15章](../../tutorials/15-berry-efield.md)
