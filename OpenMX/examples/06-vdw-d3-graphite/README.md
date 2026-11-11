# 示例 06-3：DFT-D3 范德华修正——石墨层间晶格常数与双层石墨烯层间距

## 目的

- 学会在 OpenMX 中打开 DFT-D3（BJ 阻尼）并设置 `DFTD.IntDirection`；
- 由于 4.0 手册写明晶胞优化不支持 DFT-D，用 E(c) 扫描确定 AB 石墨的层间晶格常数 c，并与不加色散修正的 PBE 对比；
- 在板模型中把 D3 与 RF 原子优化结合，优化双层石墨烯的层间距。

对应教程：[第6章 结构优化、晶胞优化、状态方程与 NEB](../../tutorials/06-optimization-eos-neb.md)（§6）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `Graphite_D3_opt.dat` | `System.Name graphite_d3`；AB 石墨（4 原子），a = 2.464 Å 固定，`MD.Type EvsLC` + `MD.EvsLC.flag 0 0 1` 只扫描 c：从 6.20 Å 起每步 +1%，共 16 步；`scf.Ngrid 24 24 72`，k 点 15×15×5 |
| `scan_c.sh` | 不依赖 EvsLC 的单点扫描（同一组 c 值），`DFTD=on/off` 切换色散修正，输出与 `.EvsLC` 同格式的文件 |
| `Graphene_bilayer_D3_RF.dat` | `System.Name blg_d3`；AB 双层石墨烯，c = 20 Å，初始层间距 3.70 Å，`MD.Type RF`，`DFTD.IntDirection 1 1 0` |
| `job.slurm` | Slurm 作业脚本 |

拟合使用 [06-cellopt-eos-si-c](../06-cellopt-eos-si-c/README.md) 中的 `fit_evslc.py`。

## 运行步骤

```bash
cd OpenMX/examples/06-vdw-d3-graphite
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" *.dat   # once

# 1) c-axis scan with PBE-D3(BJ)
mpirun -np 8 openmx Graphite_D3_opt.dat -nt 1 > graphite_d3.std
python ../06-cellopt-eos-si-c/fit_evslc.py graphite_d3.EvsLC --mode c --natoms 4 --plot graphite_d3_c.png

# 2) same c values as single points, with and without D3
DFTD=on  MPI_RUN="mpirun -np 8" bash scan_c.sh
DFTD=off MPI_RUN="mpirun -np 8" bash scan_c.sh
python ../06-cellopt-eos-si-c/fit_evslc.py graphite_on_scan.EvsLC  --mode c --natoms 4
python ../06-cellopt-eos-si-c/fit_evslc.py graphite_off_scan.EvsLC --mode c --natoms 4
#    on Slurm (the environment is exported by default):
DFTD=off sbatch job.slurm scan_c.sh

# 3) bilayer graphene: atomic relaxation with D3
mpirun -np 8 openmx Graphene_bilayer_D3_RF.dat -nt 1 > blg_d3.std
```

## 预期结果 / 如何检查

```bash
grep -E "UvdW\.|Utot\." graphite_d3.out | tail -4    # dispersion energy is part of Utot
grep -o "scf_conv=[01]" graphite_d3.md | sort | uniq -c
grep -A30 "History of" blg_d3.out
tail -n 6 blg_d3.xyz                                  # final structure (4 atoms)
```

- PBE-D3(BJ) 的 E(c) 曲线应在实验值 c ≈ 6.71 Å 附近有清晰极小（以你的拟合结果为准）；
- 纯 PBE（`DFTD=off`）的能量随 c 几乎单调下降或只有极浅的极小，`fit_evslc.py` 会提示"最低点在扫描边缘"，说明 GGA 本身描述不了层间色散作用；
- 双层石墨烯优化后层间距应收缩到 3.3–3.4 Å 量级（石墨层间距实验值约 3.35 Å），即两层 z 坐标之差；
- `fit_evslc.py --mode c` 打印的 `E(c_max) - E_min` 只是层间结合能的下限（扫描没有到无穷远）。

## 注意事项

- 4.0 手册第19.1节写明变胞优化不支持 DFT-D2/D3（3.9.7 补丁说明曾称已支持，两者不一致），所以这里用扫描代替晶胞优化。若你的版本在 EvsLC 与 DFT-D 同时使用时报错，直接用 `scan_c.sh`。
- 板模型、表面、二维材料中，把真空方向的 `DFTD.IntDirection` 设为 0，避免计入跨真空的镜像色散作用。
- 层间结合能是小量，局域基组的 BSSE 不可忽略；需要定量结合能时，用空原子做 counterpoise 修正（第4章 §8）。
- `scan_c.sh` 通过匹配 `Graphite_D3_opt.dat` 中 c 矢量那一行（`0.000000   0.000000   6.200000`）来改写 c，修改输入时请保持这一行的格式。
- 输入文件只含 ASCII 字符和 LF 换行。
