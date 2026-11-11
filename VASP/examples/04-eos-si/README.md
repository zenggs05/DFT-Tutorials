# 示例 04-4：Si 的 E–V 扫描与 Birch–Murnaghan 状态方程拟合

## 目的

- 在 7 个固定体积（V/V_ref = 0.94–1.06）下计算 Si 的总能；
- 用三阶 Birch–Murnaghan 方程拟合得到平衡体积 V0、体模量 B0 与 B0'，并换算立方晶格常数 a0；
- 与 [示例 04-1](../04-relax-si-isif3/README.md) 中 `ISIF=3` 得到的体积互相验证。

对应教程：[第4章 结构优化与状态方程](../../tutorials/04-relax-eos.md)（§8）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `INCAR` | 固定体积单点：`NSW=0`、`IBRION=-1`、`ISMEAR=-5`、`EDIFF=1E-7`；注释中给出非立方体系改用 `ISIF=4` 的写法 |
| `POSCAR` | 参考结构（a = 5.47 Å 的 fcc 原胞），第 2 行缩放因子会被脚本按体积缩放 |
| `KPOINTS` | Γ 中心 8×8×8 |
| `run_eos.sh` | `prep` 生成 `V_*` 目录；`submit` 逐个 `sbatch`；`run` 在当前分配中串行运行；`collect` 汇总为 `eos.dat` |
| `job.slurm` | 每个 `V_*` 目录中提交的单点作业脚本 |
| `fit_bm.py` | 只依赖 numpy 的三阶 BM 拟合（`--plot` 时需要 matplotlib） |
| `POTCAR` | **不提供**，需自行生成 |

## 运行步骤

```bash
cd VASP/examples/04-eos-si
cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR

bash run_eos.sh prep          # 生成 V_0.94 ... V_1.06
bash run_eos.sh submit        # 每个目录提交一个作业（或在交互节点上：MPI_RUN="mpirun -np 4" bash run_eos.sh run）
# 全部结束后：
bash run_eos.sh collect       # 生成 eos.dat（体积 A^3/cell，能量 eV/cell）
python fit_bm.py eos.dat --lattice fcc --natoms 2 --plot
```

想改体积点，可以用环境变量覆盖：`SCALES="0.90 0.93 0.96 0.99 1.02 1.05 1.08 1.11" bash run_eos.sh prep`。

## 预期结果 / 如何检查

`fit_bm.py` 会打印 V0、E0、B0（GPa）、B0'、拟合残差，以及每个点相对拟合曲线的偏差（meV）。判断标准：

- V0 落在扫描区间中部（脚本会在靠近边缘或超出范围时给出 WARNING）；
- RMS 残差远小于 1 meV；
- `--lattice fcc` 给出的 a0 应与 `ISIF=3` 重新优化后的结果相差 < 0.3%；
- B0 与 PBE 文献量级一致（约 88–90 GPa；实验约 98–99 GPa，GGA 通常低估）。

也可以用 ASE 或 pymatgen 交叉验证（代码见教程 §8.4）。用合成数据测试时，`fit_bm.py` 与 `ase.eos.EquationOfState(eos="birchmurnaghan")` 给出一致的 V0 和 B0。

## 注意事项

- `ISMEAR=-5` 只适合半导体/绝缘体的单点能；金属请改用 `ISMEAR=1`、`SIGMA=0.1~0.2` 并测试 SIGMA 的影响，能量取 `energy(sigma->0)`。
- 非立方晶体（如六角 ZnO、MoS2 体相）每个体积点都要 `ISIF=4` 优化原子与形状，并提高 `ENCUT`。
- `collect` 只收集 OUTCAR 中出现 `Total CPU time used` 的完成算例；未完成的目录会被跳过并给出警告。
- 体积范围 ±6% 对 Si 合适；软材料（碱金属、分子晶体）可以取更大范围，硬材料（金刚石）可以取更小范围，但不宜小于 ±3%，否则 B0' 不稳定。
