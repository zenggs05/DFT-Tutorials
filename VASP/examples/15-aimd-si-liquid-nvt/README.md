# 算例 15-1：液态 Si 的从头算分子动力学（NVT，Nosé–Hoover）

## 目的

用 64 个 Si 原子的立方超胞做 Born–Oppenheimer AIMD：先在 3000 K 把晶体熔化，再在 1800 K 做生产运行，然后从 XDATCAR/OSZICAR 中得到温度曲线、径向分布函数 g(r)、配位数、均方位移 MSD 与自扩散系数 D。练习 MD 的基本 INCAR 设置（`IBRION = 0`、`MDALGO`、`POTIM`、`TEBEG`、`SMASS`、`NBLOCK`）和轨迹分析流程。

对应教程：[第15章 从头算分子动力学与机器学习力场](../../tutorials/15-md-mlff.md)。

## 文件说明

| 文件 | 作用 |
|---|---|
| `POSCAR` | 64 个原子（2×2×2 金刚石惯用胞的原子排布），立方盒子边长 10.511 Å，对应液态 Si 在熔点附近的密度约 2.57 g/cm³（近似值） |
| `KPOINTS` | 只用 Γ 点 → 使用 `vasp_gam` |
| `INCAR.melt` | 第 1 段：3000 K，1500 步 × 2 fs = 3 ps，让晶体完全熔化 |
| `INCAR.prod` | 第 2 段：1800 K，5000 步 × 2 fs = 10 ps，生产运行 |
| `job.slurm` | 依次运行两段并调用 `analyze_md.py` |
| `analyze_md.py` | 读取 XDATCAR（定胞/变胞均可）和 OSZICAR：温度/能量曲线、g(r)、n(r)、MSD、D、可选振动态密度 |

## 运行步骤

```bash
cd VASP/examples/15-aimd-si-liquid-nvt
export VASP_PP_PATH=/path/to/your/paw_library
sbatch job.slurm
```

手动执行时：

```bash
cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR
cp INCAR.melt INCAR && mpirun -np 32 vasp_gam > melt.out
cp CONTCAR POSCAR          # CONTCAR 中包含速度，第 2 段接着跑
cp INCAR.prod INCAR && mpirun -np 32 vasp_gam > prod.out
python3 analyze_md.py --potim 2.0 --skip 500 --pair Si-Si
```

`--skip 500` 表示丢弃前 500 帧（1 ps）作为平衡段。

## 预期结果 / 如何检查

- 运行中监控：`grep "T=" OSZICAR | tail`（温度围绕 TEBEG 涨落）；`grep -E "external pressure|total pressure" OUTCAR | tail`（定容下的压力；MD 中应关注包含离子动能项的总压力，可用来判断密度是否合理）。
- 熔化判据：MSD 随时间线性增长（固体中 MSD 很快饱和）；g(r) 第一峰之后没有晶体那样尖锐的第二、三峰。
- 量级参考（文献值，近似）：液态 Si 的 g(r) 第一峰约 2.4–2.5 Å，平均配位数约 6（晶体为 4）；熔点附近自扩散系数为 10⁻⁴ cm²/s 量级。64 原子、10 ps 的统计误差可达 20–30%。
- `analyze_md.py` 输出格式（示意输出，数值仅供参考）：

  ```text
  <T> = 1801.3 +- 180.2 K (after skip)
  RDF Si-Si: first peak at 2.4x A, first minimum ~3.0x A, coordination n(r_min) = 6.x
  MSD Si: D = x.xxe-04 cm^2/s  (fit 1.80-5.40 ps)
  ```

## 注意事项

- 小体系温度涨落大：N 个原子的相对温度涨落约为 sqrt(2/(3N))，64 原子约 ±10%，这是正常的。
- `MAXMIX = 40` 可以显著减少每个 MD 步的电子步数；但在机器学习力场训练中**不要**使用（见本章 MLFF 部分）。
- 想得到定压下的密度，请改用 NpT（`MDALGO = 3`、`ISIF = 3`，ENCUT 提高约 30%），参见 [算例 15-2](../15-aimd-water-npt/README.md)。
- 振动态密度：`python3 analyze_md.py --potim 2.0 --skip 500 --vdos`（需要 `NBLOCK = 1`）。
- 仓库中还有功能更全的 [../../scripts/md_analysis.py](../../scripts/md_analysis.py)，用 `python ../../scripts/md_analysis.py -h` 查看用法。

返回教程：[第15章 从头算分子动力学与机器学习力场](../../tutorials/15-md-mlff.md)
