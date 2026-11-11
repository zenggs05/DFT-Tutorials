# 算例 13-2：单层 MoS2 面内弹性常数（能量-应变法，单位 N/m）

## 目的

二维材料的晶胞里有真空层，"每单位体积"的弹性常数没有意义，应当使用"每单位面积"的二维弹性常数（单位 N/m）。本算例用**能量-应变法**：对单层 2H-MoS2 施加 4 种面内应变（xx、yy、双轴、剪切），每种 7 个应变值（±1.5% 以内），固定晶胞、弛豫原子（`ISIF = 2`），再用二次多项式拟合 E(δ) 得到 C11、C22、C12、C66，最后计算二维杨氏模量 Y 和泊松比 ν。

对应教程：[第13章 弹性常数与力学性质](../../tutorials/13-elastic.md)（"实战三：二维材料"一节）。

## 文件说明

| 文件 | 作用 |
|---|---|
| `POSCAR` | 单层 2H-MoS2 原胞（3 原子），a = 3.19 Å（PBE 近似值），c = 20 Å 真空方向沿 z |
| `KPOINTS` | Γ 中心 18×18×1 |
| `INCAR` | 每个应变点共用：`ISIF = 2` 只弛豫原子，`EDIFF = 1E-7`、`EDIFFG = -1E-3` |
| `make_strain.py` | 对 POSCAR 施加指定应变：`python make_strain.py POSCAR 模式 δ 输出文件` |
| `run_strain.sh` | 批量生成 `strain_<模式>/<δ>/` 共 28 个目录，可选直接 `sbatch` |
| `job.slurm` | 每个应变目录里使用的 Slurm 作业脚本 |
| `fit_2d_elastic.py` | 收集各目录 OUTCAR 中最后的 `energy(sigma->0)`，拟合并输出 C_ij、Y、ν 与二维 Born 判据 |

## 运行步骤

1. **先得到平衡晶格常数**。参考结构必须处于零应力状态，否则拟合中会出现大的一次项。做法见 [第4章](../../tutorials/04-relax-eos.md)：`ISIF = 3` + `LATTICE_CONSTRAINTS = .TRUE. .TRUE. .FALSE.`（VASP ≥ 6.4.3）或做 E–a 扫描。把得到的 a 写回本目录的 POSCAR。
2. 生成应变结构并提交：

   ```bash
   cd VASP/examples/13-elastic-mos2-2d
   export VASP_PP_PATH=/path/to/your/paw_library   # 包含 potpaw_PBE/
   bash run_strain.sh            # 只建目录，检查无误后
   bash run_strain.sh submit     # 建目录并逐个 sbatch
   ```

   POTCAR 由脚本按 `cat $VASP_PP_PATH/potpaw_PBE/{Mo_pv,S}/POTCAR > POTCAR` 拼接（仓库不提供 POTCAR）。
3. 全部完成后拟合：

   ```bash
   grep -L "reached required accuracy" strain_*/*/OUTCAR    # 列出未收敛的目录，应为空
   python3 fit_2d_elastic.py --plot
   ```

## 预期结果 / 如何检查

- 输出表中 `c1 (eV)` 列（一次项）应接近 0（|c1| ≲ 0.01 eV）；若明显偏离，说明参考晶格没有弛豫好。
- `max|res|` 是拟合残差，应在 0.01 meV 量级；若出现 meV 级残差，检查 EDIFF、k 点或是否有目录弛豫失败。
- 六角对称要求 C11 ≈ C22，且 C66 ≈ (C11 − C12)/2，这是很好的自洽检查。
- 量级参考：文献中单层 MoS2 的 DFT 结果大约为 C11 ≈ 130–150 N/m、C12 ≈ 30–50 N/m，二维杨氏模量约 120–140 N/m、ν ≈ 0.25（不同泛函和参数差异明显，仅供参考）；AFM 纳米压痕实验值为 180 ± 60 N/m（Bertolazzi et al., ACS Nano 5, 9703 (2011), DOI: [10.1021/nn203879f](https://doi.org/10.1021/nn203879f)）。

脚本输出格式（示意输出，数值仅供参考）：

```text
A0 = 8.8128 A^2 (from POSCAR)

    mode  npts    c1 (eV)   2c2/A0 (N/m)  max|res| (meV)
      xx     7     0.0003         13x.xx           0.004
 ...
C11 =   ... N/m
C12 =   ... N/m
C66 =   ... N/m   (hexagonal check: (C11-C12)/2 = ...)
```

## 注意事项

- 所有应变点必须使用相同的 ENCUT、k 点和 POTCAR；本算例固定晶胞（`ISIF = 2`），所以 Pulay 应力不影响结果，但能量差只有 meV 量级，`EDIFF` 必须足够小。
- 应变不要太大（|δ| ≤ 2%），否则非谐项会污染二次系数；可以用 `--deg 3` 加三次项比较 C 的变化。
- 剪切模式使用工程剪切应变 γ = δ（张量分量 ε_xy = δ/2），与脚本中的能量表达式一致。
- 换算关系：若用 3D 方法（如 `IBRION = 6`）对含真空的晶胞算出 C_ij（GPa），则 C_2D(N/m) = C_3D(GPa) × c(Å) / 10，其中 c 为晶胞 z 方向长度；只有面内分量（11、22、12、66）有意义。

返回教程：[第13章 弹性常数与力学性质](../../tutorials/13-elastic.md)
