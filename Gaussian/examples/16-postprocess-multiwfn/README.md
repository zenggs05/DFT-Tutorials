# 示例 16：formchk、cubegen 与 Multiwfn 后处理（苯酚与水二聚体）

## 目的

用两个小体系演示第16章的完整后处理链：

1. 苯酚：B3LYP-D3(BJ)/def2-SVP 优化 + 频率，再在 B3LYP/def2-TZVP 下做单点，保留波函数；
2. `formchk` → 从 fchk 读取 HOMO/LUMO 能量 → `cubegen` 生成 HOMO/LUMO、电子密度和同格点 ESP 的 cube；
3. Multiwfn 命令流：苯酚分子表面 ESP 定量分析（主功能 12）、TD-DFT UV-Vis 曲线导出（主功能 11）；
4. 水二聚体的 IGMH 弱相互作用分析（主功能 20 → 11）；
5. 用 `Output=WFX` 写出 AIM 风格的波函数文件。

对应教程：[第16章 检查点文件、cubegen 与 Multiwfn 后处理](../../tutorials/16-postprocess-multiwfn.md)。

## 文件说明

| 文件 | 说明 |
|---|---|
| `phenol_opt_sp.gjf` | 两步 `--Link1--`：① `B3LYP EmpiricalDispersion=GD3BJ def2SVP Opt Freq`；② `B3LYP/def2TZVP Geom=Check Guess=Read` 单点。两步都写 `phenol.chk`，所以最后 chk 中是 def2-TZVP 的波函数。初始坐标由 RDKit/MMFF 生成，只是优化起点 |
| `phenol_td.gjf` | `%oldchk=phenol.chk` 读取结构，`CAM-B3LYP/def2TZVP TD(NStates=20) IOp(9/40=4)`，写 `phenol_td.chk`。`IOp(9/40=4)` 让 Gaussian 打印所有绝对值大于 10⁻⁴ 的组态系数，供 Multiwfn 激发态分析使用 |
| `water_dimer.gjf` | 水二聚体：def2-SVP 优化 + 频率，再做 def2-TZVP 单点。原子 1–3 为氢键供体水，4–6 为受体水（IGMH 的片段定义依赖这个顺序） |
| `wfx_output_example.gjf` | 水的单点，`Output=WFX`，在坐标后的输入段中给出文件名 `water.wfx` |
| `formchk_cubegen_homo_lumo.sh` | chk → fchk；从 fchk 读出 HOMO/LUMO 能量和能隙；生成 `*_homo.cub`、`*_lumo.cub`（默认 `npts=-3`） |
| `density_esp_cube.sh` | 生成密度 cube，再用 `npts=-1` + 模板 cube 生成**同格点**的 ESP cube，并检查两者格点是否一致 |
| `multiwfn_esp_surface_cmds.txt` | Multiwfn 命令流：主功能 12 → 开始分析 → 导出极值点（txt、pdb）和表面顶点（vtx.pdb） |
| `multiwfn_igmh_dimer_cmds.txt` | Multiwfn 命令流：主功能 20 → IGMH → 两个片段 1-3 / 4-6 → 中等格点 → 导出 cube |
| `multiwfn_uvvis_plot_cmds.txt` | Multiwfn 命令流：主功能 11 → UV-Vis → 导出谱线和展宽曲线数据 |
| `run_multiwfn_examples.sh` | 依次运行上面三个命令流，并把 Multiwfn 固定名称的输出文件改名 |
| `run_g16_jobs_slurm.sh` | Slurm 脚本：按顺序运行 4 个 Gaussian 作业并生成 fchk |

### 命令流逐行说明

命令流文件中不能写注释（Multiwfn 会把整行当作输入），所以说明放在这里：

`multiwfn_esp_surface_cmds.txt`

| 行 | 输入 | 含义 |
|---|---|---|
| 1 | `12` | 主功能 12：分子表面定量分析 |
| 2 | `0` | 用默认设置开始（ρ = 0.001 a.u. 等值面，映射函数为 ESP） |
| 3 | `1` | 后处理菜单：极值点导出为 `surfanalysis.txt` |
| 4 | `2` | 极值点导出为 `surfanalysis.pdb`（B 因子列为 ESP 值） |
| 5 | `6` | 所有表面顶点导出为 `vtx.pdb` |

`multiwfn_igmh_dimer_cmds.txt`

| 行 | 输入 | 含义 |
|---|---|---|
| 1 | `20` | 主功能 20：弱相互作用可视化 |
| 2 | `11` | IGMH 分析 |
| 3 | `2` | 定义 2 个片段 |
| 4 | `1-3` | 片段 1 的原子 |
| 5 | `4-6` | 片段 2 的原子 |
| 6 | `2` | 中等质量格点 |
| 7 | `3` | 把 cube 文件导出到当前目录（作图需要 `sl2r.cub` 和 `dg_inter.cub`） |

`multiwfn_uvvis_plot_cmds.txt`

| 行 | 输入 | 含义 |
|---|---|---|
| 1 | `11` | 主功能 11：绘制光谱 |
| 2 | `3` | UV-Vis |
| 3 | `2` | 导出 `spectrum_line.txt`（离散线）和 `spectrum_curve.txt`（展宽曲线），使用默认展宽参数 |

三个命令流结束时都没有退回主菜单，Multiwfn 读到文件结尾后会打印 `end-of-file during read` 之类的信息并退出。Multiwfn 手册说明这不是错误。

> [!IMPORTANT]
> 这些命令流是按 Multiwfn 在线手册（2026 年 11 月版本）中的菜单编号编写的，本仓库没有附带运行结果。Multiwfn 更新频繁，**第一次使用前请先交互运行一遍**，确认每一步的提示与上表一致；如果版本差异导致菜单编号变化，按屏幕提示修改命令流即可。

## 运行步骤

```bash
cd Gaussian/examples/16-postprocess-multiwfn

# 1) Gaussian jobs (phenol_td must run after phenol_opt_sp: it reads phenol.chk)
sbatch run_g16_jobs_slurm.sh                 # or run the four inputs one by one:
# g16 < phenol_opt_sp.gjf > phenol_opt_sp.log
# g16 < phenol_td.gjf     > phenol_td.log
# g16 < water_dimer.gjf   > water_dimer.log
# g16 < wfx_output_example.gjf > wfx_output_example.log

# 2) fchk, HOMO/LUMO energies and orbital cubes
bash formchk_cubegen_homo_lumo.sh phenol.chk -3

# 3) density + ESP on the same grid (4 cores for the ESP)
bash density_esp_cube.sh phenol.fchk 4

# 4) Multiwfn analyses (ESP surface, IGMH, UV-Vis)
bash run_multiwfn_examples.sh
```

在 Windows 上：G16W 运行 `.gjf`，用 G16W 的 Utilities 菜单（或命令提示符）运行 formchk/cubegen；Multiwfn 命令流可以在命令提示符中运行，例如 `Multiwfn phenol.fchk < multiwfn_esp_surface_cmds.txt > out.txt`。PowerShell 中写作 `Get-Content multiwfn_esp_surface_cmds.txt | .\Multiwfn.exe phenol.fchk > out.txt`。

## 预期结果 / 如何检查

| 检查项 | 命令 | 期望 |
|---|---|---|
| 苯酚两步都正常结束 | `grep -c "Normal termination" phenol_opt_sp.log` | `3`（Opt、Freq、单点各一次） |
| 苯酚是极小点 | `grep "imaginary frequencies" phenol_opt_sp.log` | 无输出 |
| fchk 中的波函数来自单点 | `sed -n 2p phenol.fchk` | 第 2 行显示 B3LYP 和 def2TZVP |
| HOMO/LUMO | `bash formchk_cubegen_homo_lumo.sh` 的输出 | 打印 HOMO、LUMO 能量（eV）和能隙；HOMO 应为苯环 π 轨道，在 GaussView 中检查 |
| ESP 与密度同格点 | `density_esp_cube.sh` 的输出 | `OK: ... share the same grid` |
| ESP 表面极值 | `grep "Global surface" phenol_esp_surface.out` | 全局极小位于 O 孤对电子一侧，全局极大位于羟基 H 外侧 |
| TD 激发态 | `grep "Excited State" phenol_td.log \| head` | 20 个激发态，列出 eV、nm 和振子强度 f |
| UV-Vis 数据 | `head phenol_td_spectrum_curve.txt` | 两列（能量或波长、摩尔吸光系数） |
| IGMH | `ls igmh_cubes` | 含 `sl2r.cub` 和 `dg_inter.cub`；在 VMD 中用 Multiwfn 自带的 `IGM_inter.vmd` 作图，O–H···O 区域应为蓝色小圆盘 |
| wfx | `ls -l water.wfx` | 文件存在，可被 Multiwfn 直接读取 |

本仓库没有附带输出文件和参考数值，所有数值请以你自己的计算为准。也可以用仓库脚本 `batch_cube_multiwfn.sh` 对一个目录下的所有 chk 批量生成 cube 并调用 Multiwfn（`bash ../../scripts/batch_cube_multiwfn.sh -h` 查看用法），用 [`spectrum_uvvis_ecd.py`](../../scripts/spectrum_uvvis_ecd.py) 在不装 Multiwfn 的情况下画 UV-Vis（`python ../../scripts/spectrum_uvvis_ecd.py -h`）。

## 注意事项

- `%nprocshared=8`、`%mem=8GB` 与 `run_g16_jobs_slurm.sh` 中的 `--cpus-per-task=8`、`--mem=10G` 配套；修改时要同步调整（[第17章](../../tutorials/17-hpc-jobs.md)）。
- 激发态分析（Multiwfn 主功能 18）必须用 `phenol_td.fchk` 配合 `phenol_td.log`，不能用基态的 `phenol.fchk`，因为两者的泛函不同。
- Multiwfn 的输出文件名是固定的（`surfanalysis.txt`、`spectrum_curve.txt`、`sl2r.cub` 等），批处理时每次运行后都要改名，`run_multiwfn_examples.sh` 已经这样做了。
- Linux 下运行 Multiwfn 前确认 `ulimit -s unlimited` 和 `OMP_STACKSIZE` 已设置，否则大体系可能崩溃。
- 想让 Multiwfn 调用 cubegen 加速 ESP 计算，在 `settings.ini` 中设置 `cubegenpath`。
- `water_dimer.gjf` 中的初始结构是手工搭建的近似 C<sub>s</sub> 结构，优化后应保持一个线性的 O–H···O 氢键。

[← 返回第16章](../../tutorials/16-postprocess-multiwfn.md)
