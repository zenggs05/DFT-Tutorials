# 示例 06-4：NEB——Si 中间隙 H 在相邻四面体位之间的扩散

## 目的

- 完整走一遍 OpenMX 的 NEB 三步流程：优化初态、优化终态、优化最小能量路径；
- 学会用 `<NEB.Atoms.SpeciesAndCoordinates` 给出终态、设置 image 数和弹簧常数、选择 MPI 进程数；
- 从 `.neb.ene` 读出势垒，检查每个 image 的 SCF 是否收敛。

对应教程：[第6章 结构优化、晶胞优化、状态方程与 NEB](../../tutorials/06-optimization-eos-neb.md)（§7）。

## 体系

8 原子立方 Si 胞（a = 5.43 Å）加 1 个间隙 H。H 从四面体间隙位 T1 = (½,½,½)a 扩散到相邻的 T2 = (¾,¾,¾)a，直线路径经过六角间隙位 (⅝,⅝,⅝)a；T1 与 T2 通过六角位处的反演中心互相对称。中性 H 使总电子数为 33，所以做自旋极化计算（`scf.SpinPolarization on`，H 初始占据 `1.0 0.0`）。

> [!NOTE]
> 这是教学用的小胞：H 与周期镜像只相距 5.43 Å，势垒不能直接与文献比较。实际研究需要 64/216 原子超胞，并考虑 H 的电荷态。

## 文件说明

| 文件 | 说明 |
|---|---|
| `H_Si8_T1_opt.dat` | 初态优化（`System.Name h_si8_t1`），RF，判据 1.0e-4 Hartree/Bohr |
| `H_Si8_T2_opt.dat` | 终态优化（`System.Name h_si8_t2`），与上一个文件只差 H 坐标和名称 |
| `H_diffusion_NEB.dat` | NEB 模板（`System.Name h_si8_neb`）：5 个 image，弹簧常数 0.1 Hartree/Bohr²；内含未弛豫的理想端点坐标，可直接运行用于测试 |
| `build_neb.py` | 从两个 `*.dat#` 读取弛豫后的端点坐标，生成 `H_diffusion_NEB_relaxed.dat`；检查坐标单位、原子顺序和跨周期边界的原子（只依赖 numpy） |
| `job.slurm` | Slurm 作业脚本，默认 10 个 MPI 进程（image 数 5 的整数倍） |

## 运行步骤

```bash
cd OpenMX/examples/06-neb-h-in-si
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" *.dat   # once

# steps 1 and 2: relax both end points with identical settings
mpirun -np 8 openmx H_Si8_T1_opt.dat -nt 1 > h_si8_t1.std
mpirun -np 8 openmx H_Si8_T2_opt.dat -nt 1 > h_si8_t2.std
grep "Utot\." h_si8_t1.out | tail -1 ; grep "Utot\." h_si8_t2.out | tail -1

# step 3: build the NEB input from the relaxed end points and run it
python build_neb.py                          # -> H_diffusion_NEB_relaxed.dat
mpirun -np 10 openmx H_diffusion_NEB_relaxed.dat -nt 1 > h_si8_neb.std
#   or: sbatch job.slurm H_diffusion_NEB_relaxed.dat
```

## 预期结果 / 如何检查

```bash
tail -n 5 h_si8_neb.neb.opt                  # |Maximum force| must drop below MD.Opt.criterion
cat h_si8_neb.neb.ene                        # image, Etot (Ha), distance to neighbour, distance from image 0 (Bohr)
awk '!/^#/ && NF>=4 { if (n==0) e0=$2; n++; if (n==1 || $2>emax) emax=$2;
      printf "%3d  s=%8.4f A  dE=%8.4f eV\n", $1, $4*0.529177, ($2-e0)*27.211386 }
     END { printf "forward barrier = %.3f eV\n", (emax-e0)*27.211386 }' h_si8_neb.neb.ene
for f in h_si8_neb_*.DFTSCF; do echo "== $f"; tail -n 3 "$f"; done
```

- 两个端点的总能应几乎相等（对称性），差值反映了网格位置等数值误差；
- 能量曲线应关于中点对称，最高点在六角位附近；
- image 间距（`.neb.ene` 第 3 列）应大致相等；
- 标准输出只显示 image 1 的信息，**其余 image 的 SCF 是否收敛要看各自的 `.DFTSCF` 文件或 `h_si8_neb_#.out`**；
- `h_si8_neb.neb.xyz` 可在 OpenMX Viewer、VESTA 或 OVITO 中做动画。

## 注意事项

- NEB 的坐标必须用 `Ang` 或 `AU`，不能用 `FRAC`（手册51.8节）；`build_neb.py` 会检查这一点。
- 三步计算的晶胞、截断能、基组、赝势、电子温度、k 点必须完全一致。
- MPI 进程数取 image 数的整数倍；内存不足时加 `MD.NEB.Parallel.Number 1`（或 2、3）。
- 未收敛时用 `H_diffusion_NEB_relaxed.dat#`（每步更新，含所有 image 坐标）续算。
- 想直接测试流程、跳过端点优化，也可以运行未弛豫的 `H_diffusion_NEB.dat`，但得到的不是真正的势垒。
- 输入文件只含 ASCII 字符和 LF 换行。
