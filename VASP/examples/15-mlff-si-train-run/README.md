# 算例 15-3：VASP 机器学习力场全流程（train → refit → run）

## 目的

以 64 原子晶体 Si 为例，完整走一遍 VASP 内置机器学习力场（MLFF，VASP ≥ 6.3；`ML_MODE` 需 ≥ 6.4.0）的标准流程：

1. **train**：NpT Langevin MD 从 200 K 升温到 1000 K，VASP 根据贝叶斯误差估计自动决定哪些步做 DFT、哪些步用力场，边跑边学，输出 `ML_ABN`（训练集）和 `ML_FFN`（力场）；
2. **refit**：`cp ML_ABN ML_AB`，用 `ML_MODE = refit` 重新拟合得到支持"快速预测模式"的力场（官方给出的提速约 20–100 倍）；
3. **run**：`cp ML_FFN ML_FF`，`ML_MODE = run` 做 800 K、200 ps 的纯力场 NpT 生产运行，并用溢出因子（spilling factor）监控外推风险。

对应教程：[第15章 从头算分子动力学与机器学习力场](../../tutorials/15-md-mlff.md)。

## 文件说明

| 文件 | 作用 |
|---|---|
| `POSCAR` | 64 原子晶体 Si（2×2×2 惯用胞，a = 5.468 Å，接近 PBE 平衡值） |
| `KPOINTS` | Γ 中心 2×2×2（64 原子半导体只用 Γ 点会使训练数据不够准确） |
| `INCAR_train` | 第 1 步：`ML_LMLFF = .TRUE.`、`ML_MODE = train`，NpT，`TEBEG = 200`、`TEEND = 1000`，`ENCUT = 350` |
| `INCAR_refit` | 第 2 步：`ML_MODE = refit` |
| `INCAR_run` | 第 3 步：`ML_MODE = run`、`ML_OUTBLOCK = 10`、`ML_ESTBLOCK = 50` |
| `run_mlff.sh` | Slurm 作业脚本：在 `1_train/`、`2_refit/`、`3_run/` 中依次完成三步，并自动复制 ML_ABN→ML_AB、ML_FFN→ML_FF |
| `parse_ml_logfile.py` | 解析 ML_LOGFILE：统计 STATUS 各状态步数、DFT 次数，输出最后的 ERR/BEEF/LCONF，导出 `ERR.dat`、`BEEF.dat` 等，`--plot` 画误差曲线 |

## 运行步骤

```bash
cd VASP/examples/15-mlff-si-train-run
export VASP_PP_PATH=/path/to/your/paw_library
sbatch run_mlff.sh
```

分步手动执行（推荐第一次这样做，便于在每一步之后检查）：

```bash
mkdir 1_train && cp INCAR_train 1_train/INCAR && cp POSCAR KPOINTS POTCAR 1_train/
cd 1_train && mpirun -np 32 vasp_std > vasp.out && cd ..
python3 parse_ml_logfile.py -f 1_train/ML_LOGFILE --plot      # 检查训练过程

mkdir 2_refit && cp INCAR_refit 2_refit/INCAR && cp POSCAR KPOINTS POTCAR 2_refit/
cp 1_train/ML_ABN 2_refit/ML_AB
cd 2_refit && mpirun -np 32 vasp_std > vasp.out && cd ..

mkdir 3_run && cp INCAR_run 3_run/INCAR && cp KPOINTS POTCAR 3_run/
cp 1_train/CONTCAR 3_run/POSCAR && cp 2_refit/ML_FFN 3_run/ML_FF
cd 3_run && mpirun -np 32 vasp_std > vasp.out && cd ..
```

## 预期结果 / 如何检查

训练（`1_train/ML_LOGFILE`）：

```bash
grep "Total memory consumption" ML_LOGFILE                           # 先确认内存估计不超过节点内存
grep STATUS ML_LOGFILE | grep -E "learning|critical" | grep -v "#" | wc -l
grep ERR ML_LOGFILE | grep -v "#" | tail -3                          # 训练集 RMSE: 能量/力/应力
grep BEEF ML_LOGFILE | grep -v "#" | awk '{print $2, $4, $6}' | tail # 步数、最大力贝叶斯误差、当前阈值
```

- 健康的训练：开始时几乎每步都在学习（`critical`/`learning`），之后 DFT 计算越来越稀疏；若 50 步之后仍然每步都学习，或很早就完全不再学习，都说明设置有问题。
- 升温训练中 ERR 会随温度缓慢上升，这是正常的；突然跳升往往意味着晶胞发生了异常形变或相变。
- 官方经验值：300–1000 K 下能量误差 < 1 meV/atom、力误差 30–100 meV/Å 量级通常可以接受。

refit（`2_refit/ML_LOGFILE` 的 ERR 行）一般比训练阶段的误差略小或相当。生产运行（`3_run`）：

```bash
cd 3_run
python3 ../parse_ml_logfile.py -f ML_LOGFILE         # 若有 SFF 行，会给出溢出因子统计
grep "T=" OSZICAR | tail -3
python3 ../../15-aimd-si-liquid-nvt/analyze_md.py --potim 2.0 --nblock 10 --skip 1000 --pair Si-Si
```

`analyze_md.py` 中 `--nblock 10` 对应 `ML_OUTBLOCK = 10`（每 10 步输出一帧）。

## 注意事项

- 训练和续训之间**不能**修改 INCAR 中的 DFT 参数和 POTCAR；续训时必须先 `cp ML_ABN ML_AB`（直接拿旧的 ML_AB 会丢掉上一轮的数据）。
- MLFF 只在训练覆盖过的相和温度范围内可靠：本例只训练了晶体，不能直接用于液态 Si 或表面。
- `ML_MODE = run` 只有 MPI 并行（不能用 OpenMP/GPU 加速力场部分）；需要以 scaLAPACK 编译。
- 版本差异：6.3 用 `ML_ISTART = 0/1/2/3` 代替 `ML_MODE`；6.5.0 及以前用 `ML_IERR` 代替 `ML_ESTBLOCK`；6.6 起溢出因子默认始终计算，超过 0.9 时计算会终止，`ML_OUTBLOCK` 默认 10。
- 仓库中的 [../../scripts/md_analysis.py](../../scripts/md_analysis.py) 也能解析 ML_LOGFILE，用 `python ../../scripts/md_analysis.py -h` 查看用法。

返回教程：[第15章 从头算分子动力学与机器学习力场](../../tutorials/15-md-mlff.md)
