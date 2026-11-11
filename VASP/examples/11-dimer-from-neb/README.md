# 示例 11-3：从 NEB 结果出发用 Dimer 精修鞍点（VTST，ICHAIN = 2）

## 目的

- 用 `neb2dim.pl` 把一个（粗）收敛的 NEB 转成 Dimer 输入（POSCAR + MODECAR）；
- 用 Dimer 方法把鞍点收敛到更严格的受力标准，或在更高精度设置下重新收敛鞍点而不必重跑整条 NEB；
- 学会读 DIMCAR、判断收敛并续算。

对应教程：[第11章 过渡态搜索：NEB / CI-NEB / Dimer](../../tutorials/11-neb-dimer.md)（§10）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `INCAR` | Dimer 设置：`ICHAIN = 2`、`IBRION = 3`、`POTIM = 0`、`IOPT = 2`（CG，VTST 推荐）、`EDIFF = 1E-7`、`EDIFFG = -0.01`；`DdR`、`DRotMax`、`DFNMin`、`DFNMax`；**不含 IMAGES** |
| `KPOINTS` | 与示例 11-1 相同（Γ 中心 4×4×1） |
| POSCAR、MODECAR | 由 `neb2dim.pl` 生成，不随仓库提供 |

### MODECAR 是什么

MODECAR 给出 dimer 的初始方向：**每个原子一行、三个笛卡尔分量**（原子顺序与 POSCAR 相同），整体是一个单位向量，应在反应相关原子上有较大的分量。通常不需要手写：

- `neb2dim.pl`：利用 NEB 路径在最高点附近的切线自动生成；
- `modemake.pl POSCAR1 POSCAR2`：由两个构型之差生成（例如 IS 与 TS 猜测结构）。

## 运行步骤

以示例 11-1 的 NEB 结果为起点：

```bash
cd VASP/examples/11-pt-adatom-neb          # 已完成一次 NEB 运行的根目录
vfin.pl neb_final                           # 构型写回各图像的 POSCAR
nebresults.pl                               # 生成 exts.dat（neb2dim.pl 需要）
neb2dim.pl                                  # -> saddle_dimer/（或 neb2dim.pl 2 指定图像）
cd saddle_dimer
cp ../../11-dimer-from-neb/INCAR ../KPOINTS ../POTCAR .
ls                                          # 应有 INCAR KPOINTS POTCAR POSCAR MODECAR
mpirun -np 32 vasp_std > vasp.out           # 单个构型，普通进程数即可
```

监控与续算：

```bash
tail -5 DIMCAR                    # Step  Force  Torque  Energy  Curvature  Angle
grep "energy(sigma->0)" OUTCAR | tail -1     # 鞍点能量以 OUTCAR 为准
# 未收敛时续算：
cp CENTCAR POSCAR && cp NEWMODECAR MODECAR && mpirun -np 32 vasp_std > vasp.out
```

## 预期结果 / 如何检查

- Force 逐步下降到 < 0.01 eV/Å，Curvature 保持为负，Torque 随旋转系统性下降，Angle 逐渐变小；
- 收敛后的中心构型应与 CI-NEB 的爬坡图像非常接近（本例为吸附原子位于 bridge 位）；
- 能量（OUTCAR）与 CI-NEB 鞍点能量相差应很小（主要来自收敛标准不同）；
- 再对收敛构型做示例 11-4 的频率分析，确认只有一个虚频。

## 注意事项

- Dimer 必须设置 `IOPT`（并配合 `IBRION = 3`、`POTIM = 0`），否则作业报错退出；
- Torque 不下降、曲率一直为正：降低 `EDIFF`（1E-7～1E-8），或把 `DdR` 增大到 0.01，或增大 `DRotMax`；
- DIMCAR 的 Energy 列是中心点能量（非外推到零展宽的能量），不要直接当作鞍点能量；
- 换更高 ENCUT 或更密 k 点复核鞍点时，只需在 saddle_dimer 中修改 INCAR/KPOINTS 重新运行 Dimer；
- Dimer 也可以不依赖 NEB，从极小点附近沿随机或指定方向出发寻找未知鞍点，但可能收敛到与目标反应无关的鞍点。
