# 算例 10-h2o-cube-charge：水分子的 cube 输出与四种原子电荷

## 目的

- 用 `level.of.fileout 2` 输出全部实空间数据（tden、dden、v0、vhart、vxc0 的 cube 文件）；
- 在同一次计算中得到 Mulliken、Voronoi（`Voronoi.charge on`）和 NPA（`NBO.switch on1`）电荷，再用后处理程序 `esp` 得到 ESP 拟合电荷；
- 用 Python 对 cube 积分，核对价电子数。

对应教程：[第10章 电荷分析、分子轨道与 cube 可视化](../../tutorials/10-charge-mo-cube.md) §3、§4。

## 文件说明

| 文件 | 说明 |
|---|---|
| `H2O_cube.dat` | OpenMX 输入：GGA-PBE，O6.0-s2p2d1 / H6.0-s2p1（DFT_DATA19 Standard），12 Å 立方盒，Cluster 求解器 |
| `cube_integrate.py` | 读取 cube 并积分（numpy），打印网格、体元与积分值 |
| `job.slurm` | 最简 Slurm 脚本（4 MPI 进程），算完后 grep 关键结果 |

> [!IMPORTANT]
> 运行前把 `H2O_cube.dat` 中的 `DATA.PATH` 改成你机器上包含 `PAO/` 与 `VPS/` 两个子目录的那一层目录（默认写的是官方 deb 包的安装位置 `/opt/openmx4.0/DFT_DATA19`）。

## 运行步骤

```bash
# 1) SCF (a few seconds to a minute)
sbatch job.slurm
# or interactively:
mpirun -np 4 openmx H2O_cube.dat -nt 1 > h2o.std 2>&1

# 2) check cube files
ls -lh h2o.*.cube
python cube_integrate.py h2o.tden.cube h2o.dden.cube

# 3) ESP fitting charges (compile once)
cd ~/openmx4.0/source && make esp && cd -     # esp is copied to ~/openmx4.0/work/
~/openmx4.0/work/esp h2o -c 0 -s 1.4 2.0 | tee h2o.esp
```

## 预期结果 / 如何检查

```bash
grep "Utot\." h2o.out                               # total energy (Hartree)
grep -A 12 "Dipole moment (Debye)" h2o.out          # dipole
grep -A 8 "Mulliken populations" h2o.out | head -12 # Mulliken populations
grep -A 14 "Voronoi charges" h2o.out                # Voronoi charges and volumes
grep -B 2 -A 6 "Total  :" h2o.std                   # NPA (natural populations)
grep "Fitting Effective Charge" h2o.esp             # ESP charges
```

- `cube_integrate.py` 对 `h2o.tden.cube` 的积分应接近 8（O 6 + H 1 + H 1 个价电子），对 `h2o.dden.cube` 的积分应接近 0；
- 偶极矩应与实验气相值（约 1.85 D）同一量级；
- Mulliken 表中的 Sum 是**电子布居**，净电荷 = 价电子数 − Sum；Mulliken、ESP、NPA 都应给出 O 带负电、H 带正电；Voronoi 在分子中可能给出与直觉相反的符号（H 的多面体伸入真空），这是方法本身的特性；
- ESP 电荷之和应为 0（`-c 0` 约束总电荷），输出的偶极应与 `.out` 中的 DFT 偶极接近。

具体数值依赖基组与网格，以你的计算为准。

## 注意事项

- 这是一个 Cluster 计算，盒子必须足够大；改小盒子后 ESP 电荷与偶极会变化。
- `level.of.fileout 2` 会生成较大的 `.grid` 文件，只是为了演示；日常使用 1 即可。
- `NBO.switch on1` 只对 Cluster 有效；大体系请改用 Krylov 求解器 + `NBO.switch on2`。
- 若想把 scfout 交给其它后处理程序（jx、polB 等），把 `HS.fileout` 改为 `on`。
- 想看 Bader 电荷，可以运行外部程序 `bader h2o.tden.cube`；注意 OpenMX 的密度只含价电子。

[返回第10章](../../tutorials/10-charge-mo-cube.md)
