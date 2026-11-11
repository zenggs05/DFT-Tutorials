# 示例 02：输入文件结构——笛卡尔坐标、Z-matrix、Link1、片段标签与连接表

配套教程：[第2章 输入文件结构详解](../../tutorials/02-input-anatomy.md)

## 目的

用 5 个小输入演示 Gaussian 输入文件的各种写法：

- 完整的 Link 0 / route / title / 电荷多重度 / 笛卡尔坐标；
- Z-matrix 与 `Variables:` 段；
- `--Link1--` 多步作业，用 `Geom=Check Guess=Read` 从 chk 读结构和初猜；
- 原子的 `(Fragment=n)` 标签与多片段电荷多重度行；
- GaussView 风格的 `Geom=Connectivity` 连接表。

另附一个电子数奇偶校验小脚本 `check_charge_mult.py`。

## 文件说明

| 文件 | 体系 / 方法 | 知识点 |
|---|---|---|
| `ethanol_cartesian.gjf` | 乙醇，B3LYP-D3(BJ)/6-31G(d) Opt Freq | 笛卡尔坐标的标准输入 |
| `h2o2_zmatrix.gjf` | H2O2，B3LYP/6-31G(d) Opt | Z-matrix、变量共用（两个 O–H 共用 R1） |
| `link1_opt_then_sp.gjf` | 甲醛，B3LYP-D3(BJ)/def2-SVP Opt Freq → ωB97X-D/def2-TZVP 单点 | `--Link1--`、`Geom=Check`、`Guess=Read` |
| `fragment_labels_example.gjf` | H2O···NH3 氢键复合物，B3LYP-D3(BJ)/6-31+G(d,p) | `(Fragment=n)`、`0 1 0 1 0 1`、`Counterpoise=2` |
| `geom_connectivity_example.gjf` | 乙酸，B3LYP-D3(BJ)/6-31G(d) Opt Freq | `Geom=Connectivity` + 连接表（键级 1.0/2.0） |
| `check_charge_mult.py` | — | 检查电荷、多重度与电子数奇偶性是否相容（只依赖 Python 标准库） |

初始结构：乙醇、H2O2、甲醛、乙酸取自 G2 数据集的几何（经 ASE 的 `ase.build.molecule` 获得），H2O···NH3 复合物是手工搭建的近似结构（N···H 约 2.0 Å），均只作为优化或单点的起点。

## 运行步骤

```bash
cd Gaussian/examples/02-input-anatomy

# 0) check charge/multiplicity parity of every input (no Gaussian needed)
python check_charge_mult.py *.gjf

# 1) check the route syntax of one input (needs Gaussian's testrt)
testrt "#p B3LYP/6-31G(d) EmpiricalDispersion=GD3BJ Opt Freq Geom=Connectivity"

# 2) run the jobs (each takes seconds to a few minutes on 4 cores)
for f in ethanol_cartesian h2o2_zmatrix link1_opt_then_sp fragment_labels_example geom_connectivity_example; do
    g16 < $f.gjf > $f.log
done
# or reuse the helper from example 01:
bash ../01-first-job-h2o/run_g16_local.sh *.gjf
```

在集群上请改用 Slurm 提交（参考 `../01-first-job-h2o/run_g16_slurm.sh` 或第17章的模板）。

## 预期结果 / 如何检查

```bash
# every job must end normally
for f in *.log; do echo "$f: $(grep -c 'Normal termination' $f) normal, $(grep -c 'Error termination' $f) error"; done
```

| 文件 | Normal termination 行数 | 检查点 |
|---|---|---|
| `ethanol_cartesian.log` | 2（Opt + 自动生成的 Freq 步） | `grep "Stationary point found"`；`grep "Frequencies --"` 第一行不应有负值 |
| `h2o2_zmatrix.log` | 1 | 优化后 O–O 约 1.45–1.47 Å，H–O–O–H 二面角在 110–120° 范围（近似值，请以你的输出为准） |
| `link1_opt_then_sp.log` | 3（Opt、Freq、单点） | 最后一个 `SCF Done` 是 `E(RwB97XD)`；单点步中没有坐标输入，结构来自 chk |
| `fragment_labels_example.log` | 1 | `grep "Counterpoise corrected energy"`、`grep "BSSE energy"`（读法见第14章） |
| `geom_connectivity_example.log` | 2 | 结果应与不写连接表时一致（普通 DFT 计算中连接表基本不起作用） |

```bash
grep "SCF Done" link1_opt_then_sp.log | tail -n 1
grep -e "Counterpoise corrected energy" -e "BSSE energy" fragment_labels_example.log
```

## 注意事项

- 所有输入都以空行结尾，`--Link1--` 前也有空行。修改时保持这一点。
- `link1_opt_then_sp.gjf` 两步的 `%chk` 必须同名，第 2 步才能 `Geom=Check`。如果想保留第 1 步的 chk 不被覆盖，可以把第 2 步改成 `%oldchk=h2co_link1.chk` 加 `%chk=h2co_sp_tz.chk`，并改用 `Geom=AllCheck`（这时要删掉第 2 步的 title 和 `0 1` 行）。
- 删除 `geom_connectivity_example.gjf` 的连接表时，必须同时删掉 route 中的 `Geom=Connectivity`。
- `fragment_labels_example.gjf` 的几何不是优化结构，只用来演示输入格式。真正的结合能计算要先优化复合物（见第14章）。
- `%mem`、`%nprocshared` 按自己的机器和作业申请调整。
