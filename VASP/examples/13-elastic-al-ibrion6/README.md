# 算例 13-1：fcc Al 弹性常数（应力-应变法，IBRION=6 + ISIF=3）

## 目的

用 VASP 内置的有限差分（`IBRION = 6`, `ISIF = 3`）一次性得到 fcc Al 的弹性张量 C_ij，练习三件事：

1. 先用同一套参数做 `ISIF = 3` 充分弛豫，消除残余应力和 Pulay 应力；
2. 从 OUTCAR 读出 `TOTAL ELASTIC MODULI (kBar)`，除以 10 换成 GPa，并把 VASP 的 XX YY ZZ **XY YZ ZX** 顺序重排成 Voigt 的 xx yy zz **yz zx xy**；
3. 判断 Born 稳定性，计算 Voigt-Reuss-Hill 平均的 B、G、E、ν、Pugh 比和德拜温度。

对应教程：[第13章 弹性常数与力学性质](../../tutorials/13-elastic.md)。

## 文件说明

| 文件 | 作用 |
|---|---|
| `POSCAR` | fcc Al 原胞（1 个原子），初始 a = 4.04 Å；原胞的立方轴与笛卡尔轴平行，所以输出的张量就是标准立方取向 |
| `KPOINTS` | Γ 中心 21×21×21 网格（金属，C44 对 k 点敏感，请自行测试 15/21/27） |
| `INCAR.relax` | 第 1 步：`ISIF = 3` 变胞弛豫（ENCUT、k 点、SIGMA 与第 2 步完全相同） |
| `INCAR.elastic` | 第 2 步：`IBRION = 6`、`ISIF = 3`、`NFREE = 2`、`POTIM = 0.015`、`EDIFF = 1E-8` |
| `job.slurm` | Slurm 作业脚本：自动拼 POTCAR → 弛豫两轮 → 弹性计算 → 调用 `voigt_vrh.py` |
| `voigt_vrh.py` | 后处理：读 OUTCAR、换单位、重排 Voigt、Born 判据、VRH 平均、德拜温度（只依赖 numpy） |

> [!IMPORTANT]
> 仓库里不提供 POTCAR。作业脚本会执行 `cat $VASP_PP_PATH/potpaw_PBE/Al/POTCAR > POTCAR`，请先把环境变量 `VASP_PP_PATH` 指向你自己的赝势库（目录下应有 `potpaw_PBE/`）。也可以用 `vaspkit` 的 103 功能生成。

## 运行步骤

```bash
cd VASP/examples/13-elastic-al-ibrion6
export VASP_PP_PATH=/path/to/your/paw_library      # 包含 potpaw_PBE/ 的目录
sbatch job.slurm                                     # 或交互式按下面手动执行
```

手动执行（便于理解每一步）：

```bash
cat $VASP_PP_PATH/potpaw_PBE/Al/POTCAR > POTCAR
cp INCAR.relax INCAR
mpirun -np 16 vasp_std > relax1.out
cp CONTCAR POSCAR                    # 再弛豫一次，直到体积不再变化
mpirun -np 16 vasp_std > relax2.out
grep "volume of cell" OUTCAR | tail -2
cp INCAR.elastic INCAR
mpirun -np 16 vasp_std > elastic.out
python3 voigt_vrh.py -f OUTCAR --poscar POSCAR
```

## 预期结果 / 如何检查

1. 弛豫：`grep "external pressure" OUTCAR | tail -1` 应接近 0（|P| < 1 kB）；两轮弛豫的体积差 < 0.1%。
2. 弹性计算结束后：

   ```bash
   grep -A 10 "TOTAL ELASTIC MODULI" OUTCAR
   grep -A 10 "SYMMETRIZED ELASTIC MODULI" OUTCAR      # 刚性离子（clamped-ion）部分
   ```

   对单原子原胞，离子弛豫贡献为零，两组数应相同；立方晶系应满足 C11=C22=C33、C12=C13=C23、C44=C55=C66，其余为 0（数值噪声 < 0.5 GPa）。
3. 参考量级：室温实验值约为 C11 ≈ 107、C12 ≈ 61、C44 ≈ 28 GPa（Kittel《固体物理导论》表格，近似）；用这组实验值计算得到 B ≈ 76 GPa、G_H ≈ 26 GPa、E ≈ 70 GPa、ν ≈ 0.35，可用来检验脚本。PBE 结果与实验通常相差 5–15%，C44 尤其依赖 k 点与 SIGMA。

`voigt_vrh.py` 的输出格式（示意输出，数值仅供参考）：

```text
# source: TOTAL ELASTIC MODULI (kBar)  (converted to GPa, Voigt order 1..6 = xx yy zz yz zx xy)
[[ 105.  59.  59.   0.   0.   0.]
 ...
Born stability (all eigenvalues > 0): YES
                       Voigt     Reuss      Hill
           B (GPa)     ...
```

## 注意事项

- 弛豫与弹性计算必须使用**完全相同**的 ENCUT、k 点、SIGMA、PREC、LREAL；只改 IBRION/ISIF/NFREE/POTIM。
- 金属体系：把 SIGMA 从 0.2 → 0.1 → 0.05 测试一遍，同时加密 k 点，C44 变化 < 2 GPa 再算收敛。
- 不要在 `IBRION = 6` 计算中使用 Selective dynamics（官方说明只有 `IBRION = 5` 支持）。
- 想看刚性离子张量：`python3 voigt_vrh.py --clamped`。
- 仓库中另有更完整的后处理脚本 [../../scripts/elastic_post.py](../../scripts/elastic_post.py)（晶系判据、声速等），用 `python ../../scripts/elastic_post.py -h` 查看用法。

返回教程：[第13章 弹性常数与力学性质](../../tutorials/13-elastic.md)
