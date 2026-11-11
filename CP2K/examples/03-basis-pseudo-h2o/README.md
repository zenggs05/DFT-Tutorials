# 示例 03：基组、赝势与 GAPW 全电子计算（Basis Sets, GTH Pseudopotentials & GAPW for H2O）

## 目的

用水分子（以及一个甲醇的模板示例）练习：

1. 经典 MOLOPT 基组阶梯 SZV → DZVP → TZVP → TZV2P → TZV2PX，以及短程版 MOLOPT-SR，比较总能量、偶极矩与耗时；
2. 新推荐的 **UZH 基组 + UZH 赝势成对使用**，以及 2026.2 新增的带版本号快照文件；
3. **GAPW 全电子计算**（`METHOD GAPW` + `POTENTIAL ALL` + 全电子基组 `SVP-MOLOPT-PBE-ae`、`6-31G*`）；
4. 用 `@INCLUDE` 复用的 `&KIND` 模板，按泛函家族（PBE/SCAN/PBE0）自动拼出匹配的基组与赝势名，并可选加入 ADMM 辅助基组。

对应教程 [第3章 GPW/GAPW 方法、基组与赝势](../../tutorials/03-gpw-basis-pseudo.md)。

## 文件说明

| 文件 | 说明 |
|---|---|
| `h2o_szv_dzvp_tzv2p.inp` | 基组阶梯：`-E BASIS=SZV/DZVP/TZVP/TZV2P/TZV2PX`、`-E FAMILY=MOLOPT/MOLOPT-SR`；打印偶极矩（`&MOMENTS`） |
| `run_basis_scan.sh` | 依次运行上面的基组阶梯，汇总总能量、偶极矩、SCF 步数和耗时 |
| `h2o_uzh_pair.inp` | `BASIS_MOLOPT_UZH` + `POTENTIAL_UZH` 成对使用；`-E Q=DZVP/TZVP/TZV2P`；`-E BFILE=… -E PFILE=…` 切换到版本快照 |
| `h2o_gapw_allelectron.inp` | GAPW 全电子；`-E AEB=SVP-MOLOPT-PBE-ae / TZVPP-MOLOPT-PBE-ae / 6-31Gx` |
| `kind_templates.inc` | H/C/N/O 的 `&KIND` 模板（变量 `Q`、`XCTAG`、`ADMM`） |
| `ch3oh_kind_templates.inp` | 用 `@INCLUDE 'kind_templates.inc'` 计算甲醇（近似构型，未优化） |

## 运行步骤

```bash
cd CP2K/examples/03-basis-pseudo-h2o

# 0) 语法检查
for f in h2o_szv_dzvp_tzv2p h2o_uzh_pair h2o_gapw_allelectron ch3oh_kind_templates; do cp2k.psmp -c $f.inp; done

# 1) 基组阶梯（7 个计算，单核共需几分钟）
bash run_basis_scan.sh
#    并行：CP2K_EXE="mpirun -np 4 cp2k.psmp" bash run_basis_scan.sh

# 2) UZH 成对基组/赝势
rm -f h2o_uzh*.out
cp2k.ssmp -i h2o_uzh_pair.inp -o h2o_uzh.out
cp2k.ssmp -E Q=TZV2P -i h2o_uzh_pair.inp -o h2o_uzh_tzv2p.out
#    CP2K >= 2026.2：使用冻结的数据快照（保证与以后版本可复现）
cp2k.ssmp -E BFILE=BASIS_MOLOPT_UZH_2026.2 -E PFILE=POTENTIAL_UZH_2026.2 \
          -i h2o_uzh_pair.inp -o h2o_uzh_2026.2.out

# 3) GAPW 全电子
rm -f h2o_gapw*.out
cp2k.ssmp -i h2o_gapw_allelectron.inp -o h2o_gapw_svp.out
cp2k.ssmp -E AEB=TZVPP-MOLOPT-PBE-ae -i h2o_gapw_allelectron.inp -o h2o_gapw_tzvpp.out
cp2k.ssmp -E AEB=6-31Gx -i h2o_gapw_allelectron.inp -o h2o_gapw_631gs.out

# 4) KIND 模板
cp2k.ssmp -e ch3oh_kind_templates.inp | grep -A3 "&KIND"     # 查看展开后的基组/赝势名
cp2k.ssmp -i ch3oh_kind_templates.inp -o ch3oh_dzvp.out
```

## 预期结果 / 如何检查

```bash
grep "ENERGY| Total FORCE_EVAL" *.out
grep -A1 "Dipole moment \[Debye\]" h2o_*.out        # 偶极矩（最后一列 Total=）
grep "count for grid" h2o_TZV2P-MOLOPT.out           # 弥散函数多 -> 更多高斯落在粗网格
grep "Total charge density (r-space)" h2o_gapw_svp.out
grep "Total Rho_soft + Rho1_hard - Rho1_soft (r-space)" h2o_gapw_svp.out
```

- **基组阶梯**：总能量随基组增大单调下降（变分原理），SZV → DZVP 降幅最大，TZV2P → TZV2PX 变化很小。偶极矩从 SZV 到 DZVP 变化明显，之后逐渐稳定；可与气相实验值约 1.85 D 比较（PBE 本身也有误差，不要求完全一致）。MOLOPT-SR 的能量略高于同级 MOLOPT，但凝聚相计算更快。（定性预期，具体数值以你的运行结果为准。）
- **不同基组/赝势家族的总能量不能直接比较**：UZH 与经典 GTH 的赝势参数不同，GAPW 全电子的总能量（约 −76 Ha 量级）与赝势计算（约 −17 Ha 量级）更是完全不同的参考点。比较的应是能量差、几何结构、偶极矩等可观测量。
- **GAPW 电荷检查**：上面两行电荷积分应非常接近 −10（10 个电子，CP2K 中电子密度记为负值）；偏离较大时需提高 `LEBEDEV_GRID`/`RADIAL_GRID` 或收紧 `EPSFIT` 等参数。
- **版本快照**：在 2026.2 中，`BASIS_MOLOPT_UZH` 与 `BASIS_MOLOPT_UZH_2026.2` 内容相同，所以 `h2o_uzh.out` 与 `h2o_uzh_2026.2.out` 的能量应一致；在更早的版本中，快照文件不存在，会报 “The specified OLD file ... cannot be opened”。

## 注意事项

- 基组名后缀 `-qN` 必须与赝势 `GTH-xxx-qN` 的价电子数一致；UZH 基组与 UZH 赝势、经典 MOLOPT 基组与 `GTH_POTENTIALS` 应**成对**使用。两个库里都有名为 `GTH-PBE-q6` 的氧赝势，但参数并不相同。
- CP2K 使用球谐高斯函数（d 壳层 5 个分量），而 Gaussian 程序对 6-31G* 默认使用笛卡尔 6d，因此同名基组的结果与 Gaussian 不会逐位一致。
- `kind_templates.inc` 中的 `XCTAG` 必须与 `&XC` 中实际使用的泛函一致（例如 `XCTAG=PBE0` 要配合 PBE0 杂化泛函设置，见 [第9章](../../tutorials/09-hybrid-admm.md)）；本目录的示例只使用 PBE。
- 截断 `CUTOFF`/`REL_CUTOFF` 在这里取得偏保守，以便比较基组时网格误差可以忽略；收敛方法见 [第4章](../../tutorials/04-cutoff-convergence.md)。

返回教程：[第3章 GPW/GAPW 方法、基组与赝势](../../tutorials/03-gpw-basis-pseudo.md)
