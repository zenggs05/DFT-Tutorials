# 示例 16-1：Al(111) 功函数（ESM on1 + 空原子层）

配套教程：[第16章 表面与界面：ESM、功函数、电荷掺杂](../../tutorials/16-surface-esm-workfunction.md)

## 目的

用 ESM `on1`（真空|真空|真空）边界计算 6 层 Al(111) 板模型的功函数，并比较“不加空原子层”与“每侧加 2 层空原子（empty atom，`E.vps`）”的差别。这样可以体会两件事：用局域 PAO 基组描述真空一侧的波函数尾巴时要多加小心；功函数 $\Phi=\bar V_\mathrm{vac}-\mu$ 是怎么从 `v0.cube` 的平面平均势里读出来的。

## 文件说明

| 文件 | 说明 |
|---|---|
| `Al111_ESM_on1.dat` | 6 层 Al(111)，1×1 表面胞，法向沿 a 轴（=x），胞长 50 Å，ESM `on1`，不加空原子 |
| `Al111_ESM_on1_2E.dat` | 同上，每侧外加 2 层空原子（物种 `E`：Al7.0-s2p2d1 基组 + `E` 赝势，电子数 0），共 10 个格点 |
| `make_slab_omx.py` | 用 ASE 生成符合 OpenMX ESM 约定（a 轴=法向）的坐标块和晶胞块，可选空原子层、`ESM.direction x/y/z` |
| `workfunc.py` | 读取 `gcube2oned` 的输出和 `.out` 中的化学势，分别给出左右两侧真空的功函数，可选作图 |
| `job_slurm.sh` | Slurm 作业脚本：跑 OpenMX，检查是否正常结束，若 `gcube2oned` 可用则自动分析 |

结构用的是实验晶格常数 $a_0=4.05$ Å，体相截断、未弛豫：面内 $a=2.8638$ Å，层间距 2.3383 Å。

## 运行步骤

```bash
cd OpenMX/examples/16-al111-esm-workfunction
# 1) 把 DATA.PATH 改为你的 DFT_DATA19 绝对路径（目录下应有 PAO/ 与 VPS/）
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" *.dat

# 2) 编译 gcube2oned（4.0 源码包里它在 source/ 目录，手册写的是 work/）
gcc ~/openmx4.0/source/gcube2oned.c -lm -o gcube2oned

# 3) 运行（工作站示例：8 MPI x 1 线程；集群上用 sbatch job_slurm.sh <input>）
for f in Al111_ESM_on1.dat Al111_ESM_on1_2E.dat; do
  n=$(awk 'tolower($1)=="system.name"{print $2}' $f)
  mpirun -np 8 openmx $f -nt 1 > $n.std 2>&1
  grep -q "The calculation was normally finished." $n.std || echo "$f FAILED"
done

# 4) 平面平均势与功函数
for n in al111_esm_on1 al111_esm_on1_2e; do
  ./gcube2oned $n.v0.cube 1 > ${n}_1d_pot.txt
  python3 workfunc.py ${n}_1d_pot.txt $n.out --plot ${n}_wf.png
done
```

## 预期结果 / 如何检查

```bash
grep "Boundary condition" al111_esm_on1_2e.std              # Vacuum|Vacuum|Vacuum
grep "Required cutoff energy\|Num. of grids" al111_esm_on1_2e.out
grep -A3 "SCF history" al111_esm_on1_2e.out | head; tail -n 5 al111_esm_on1_2e.std
grep "Chemical potential (Hartree)" al111_esm_on1*.out
```

- 两侧真空的势应当平坦（`workfunc.py` 打印的 std 小于 1 meV 量级），而且左右两侧的功函数相同，因为这是对称板。
- 作为参考，OpenMX 4.0 手册第 50 章给出的数值是：用 7 层空原子时 Al(111) 计算值 4.19 eV，实验值 4.26±0.03 eV；手册的测试表明，加 2 层空原子后结果已基本收敛。本示例的层数、基组和网格都和手册不同，结果会有零点几 eV 以内的差别，属正常。
- 不加空原子的 `Al111_ESM_on1` 通常给出偏离的功函数，两者之差就是“真空侧基组不完备”带来的误差。

## 注意事项

- 块（`<...>` 与 `...>` 之间）里的行不能加注释，OpenMX 是按列位置读取的。
- ESM 要求第一条晶格矢量（a 轴）就是表面法向，b、c 与之正交；`ESM.direction` 只用来说明 a 轴沿哪个笛卡尔方向。
- 4.0 源码中 `ESM.buffer.range` 只有在 `ESM.wall.switch on` 时才会读取，原子与边界距离的检查代码也已停用，所以程序不会替你检查原子离边界够不够远。要自己保证原子（包括空原子）与 x=0、x=|a| 之间留有足够距离，本例为 14.5 Å。
- 收敛测试建议依次改变：空原子层数（0/1/2/3）、Al 层数（6/8/10）、`scf.Kgrid`（1×12×12 → 1×18×18）、`scf.energycutoff`（200/250/300 Ry），观察功函数变化是否小于 0.02 eV。
