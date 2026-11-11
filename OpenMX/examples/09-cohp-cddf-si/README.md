# 算例 9-3：硅的 COHP/COOP 成键分析与光学介电函数（CDDF）

## 目的

- **COHP/COOP（OpenMX ≥ 4.0）**：在 DOS 计算的同时用 `COHP.calc.flag` 输出 `si_cohp.cohp`，再用后处理程序 `cohp` 得到最近邻（两组等价的键）与次近邻 Si–Si 的能量分辨 COHP/COOP 和积分值 ICOHP，判断成键/反键。
- **光学性质（CDDF）**：用 Kubo–Greenwood 公式计算 Si 的光学电导率、介电函数 ε(ω)、吸收和反射谱，练习 `CDDF.Kgrid`、`CDDF.FWHM` 等关键词，并与实验 E₁/E₂ 峰对比。

对应教程：[第9章 态密度 DOS/PDOS、COHP 与光学性质](../../tutorials/09-dos-pdos-cohp-optics.md)。

## 文件说明

| 文件 | 作用 |
|---|---|
| `Si_COHP.dat` | Si 原胞（a = 5.43 Å，`Si7.0-s2p2d1`）；`Dos.fileout on`、`Dos.Kgrid 16 16 16`；3 个原子对：(1,2,000)、(1,2,−100) 为最近邻，(1,1,100) 为次近邻；`System.Name si_cohp` |
| `cohp_answers.txt` | `cohp` 的回答：四面体法（1）、能量网格 3000、原子对 1 |
| `Si_CDDF.dat` | 同一结构；`CDDF.start on`、`CDDF.Kgrid 30 30 30`、`CDDF.FWHM 0.2`、最大能量 10 eV（+1 eV）；`System.Name si_cddf` |
| `plot_cohp_cddf.py` | `cohp` 子命令画 −COHP / −ICOHP 并打印 ICOHP(E_F)；`optics` 子命令画 ε₁、ε₂ 与吸收谱并给出 ε₂ 峰位 |
| `job.slurm` | 完整流程（COHP 三个原子对循环 + CDDF） |

## 运行步骤

1. 修改两个 `.dat` 中的 `DATA.PATH`；编译 `cohp`：

   ```bash
   cd /path/to/openmx4.0/source && make cohp && cp cohp ~/bin/
   ```

2. COHP：

   ```bash
   cd OpenMX/examples/09-cohp-cddf-si
   mpirun -np 8 openmx Si_COHP.dat -nt 1 | tee si_cohp.std
   cohp si_cohp.cohp < cohp_answers.txt                       # pair 1
   for P in 2 3; do printf '1\n3000\n%d\n' $P | cohp si_cohp.cohp; done
   python3 plot_cohp_cddf.py cohp si_cohp_i1_s0_tetra.cohp si_cohp_i2_s0_tetra.cohp \
           si_cohp_i3_s0_tetra.cohp --out si_cohp.png
   ```

3. 光学：

   ```bash
   mpirun -np 32 openmx Si_CDDF.dat -nt 1 | tee si_cddf.std
   grep -n "Optical calculations" si_cddf.std
   head -n 4 si_cddf.df_im
   python3 plot_cohp_cddf.py optics si_cddf --out si_optics.png
   ```

## 预期结果 / 如何检查

- COHP 输出文件：`si_cohp_i1_s0_tetra.cohp`（全轨道求和）、`si_cohp_i1_o0_s0_tetra.cohp` 等（按原子 1 的轨道分解）、`si_cohp_i1_o0_o0_s0_tetra.cohp` 等（完全轨道分辨），以及对应的 `.coop`。列：E（eV）、COHP、积分 COHP。
- `plot_cohp_cddf.py cohp` 输出（示意输出，数值仅供参考）：

  ```text
  si_cohp_i1_s0_tetra.cohp                      ICOHP(E_F) =  -x.xxxx   (bonding if negative)
  si_cohp_i2_s0_tetra.cohp                      ICOHP(E_F) =  -x.xxxx   (bonding if negative)
  si_cohp_i3_s0_tetra.cohp                      ICOHP(E_F) =  -0.xxxx   (bonding if negative)
  ```

  判据：对 1 与对 2 是等价的最近邻键，ICOHP 应相同；费米能以下几乎全是成键态（−COHP > 0），导带出现反键态；次近邻的 |ICOHP| 明显更小。
- CDDF 输出：`si_cddf.cd_re/.cd_im`（电导率，S/m）、`.df_re/.df_im`（介电函数）、`.absorption`、`.extinction`、`.transmission`、`.reflection`、`.refractive_index`；每个文件的列为能量、xx、xy、…、zz、(xx+yy+zz)/3，单位写在文件头。
- `plot_cohp_cddf.py optics` 输出（示意输出，数值仅供参考）：

  ```text
  energy grid: 0.0000 - 10.0000 eV, 2001 points
  eps1 at the lowest energy (0.000 eV) = 1x.xxx  (electronic, no local-field effects)
  eps2 peaks near: 3.xx eV, 4.xx eV
  ```

  实验 ε₂ 主峰 E₁ ≈ 3.4 eV、E₂ ≈ 4.3 eV（Aspnes & Studna, PRB 27, 985, 1983）；PBE 低估带隙，计算峰位整体红移约 0.5 eV。30×30×30 的 k 网格尚未完全收敛（手册测试到 100×100×100），正式计算至少比较 50×50×50。

## 注意事项

- COHP 与 DOS 同时计算，`Dos.fileout` 必须为 `on`；`<COHP.pairs` 每行 6 个整数：序号、原子 1、原子 2、原子 2 所在晶胞的平移 l1 l2 l3。符号约定：COHP 负为成键，COOP 正为成键。
- COHP 绝对值依赖基组，适合比较趋势；不要与平面波 + LOBSTER 的数值直接比较。
- CDDF 在 4.0 中只支持共线计算；nc/SOC 下可能段错误或报 `Parameter 10 was incorrect on entry to ZHEEVX`（参见 klsk007/openmx-soc-cddf-fix）。4.0 删除了 `CDDF.material_type`，新增 `CDDF.shift_ChemP`（平移化学势，不是剪刀算符）。
- 安装检查：`mpirun -np 56 ./openmx -runtestCDDF -nt 2`，查看 `runtestCDDF.result`（差值 ≤ 0.1 视为正常）。

返回教程：[第9章 态密度 DOS/PDOS、COHP 与光学性质](../../tutorials/09-dos-pdos-cohp-optics.md)
