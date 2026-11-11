# 算例 08-mos2-monolayer：单层 MoS2——二维材料的 vdW 一致设置、真空收敛与真空能级

## 目的

- 用与双层/异质结**相同的 vdW 设置**（`IVDW = 12`）优化单层 MoS2 的面内晶格（`LATTICE_CONSTRAINTS`，c 轴固定）；
- 做真空层厚度收敛测试（c = 15/20/25/30 Å）；
- 用 `LVHAR = .TRUE.` 输出静电势 LOCPOT，求真空能级，把带边（VBM/CBM）对齐到真空能级（电离能 IP、电子亲和能 EA）。

只做"面内晶格优化"的基础版本见 [第4章的 MoS2 算例](../04-relax-mos2-2d/README.md)；本算例侧重第 8 章的二维材料专题。对应教程：[第8章](../../tutorials/08-vdw-2d.md) §7。

## 文件说明

| 文件 | 说明 |
|---|---|
| `POSCAR` | 1H-MoS2 单层，a = 3.19 Å，c = 23.13 Å（约 20 Å 真空），原子层位于 z = c/2 |
| `KPOINTS` | Γ 中心 12×12×1（真空方向只取 1 个 k 点） |
| `INCAR_relax` | PBE-D3(BJ) + `ISIF=3` + `LATTICE_CONSTRAINTS = .TRUE. .TRUE. .FALSE.` |
| `INCAR_scf` | 静态计算，`LVHAR=.TRUE.`（LOCPOT 只含离子势 + Hartree 势）、`LORBIT=11`；可选 `IDIPOL=3` + `LVACPOTAV`（≥ 6.4.3） |
| `set_vacuum.py` | 保持原子笛卡尔坐标不变、修改 c 长度并把原子层居中 |
| `vacuum_test.sh` | 真空收敛测试：`prepare` / `collect` |
| `job.slurm` | 最简 Slurm 脚本 |

POTCAR（与第4章保持一致，用 `Mo_pv`；VASP 官方推荐列表中 Mo 为 `Mo_sv`，两者都可以，但同一项目中要统一）：

```bash
cat $VASP_PP_PATH/potpaw_PBE/{Mo_pv,S}/POTCAR > POTCAR
```

## 运行步骤

```bash
# 1) in-plane relaxation
mkdir -p relax && cp INCAR_relax relax/INCAR && cp POSCAR KPOINTS POTCAR job.slurm relax/
(cd relax && sbatch job.slurm)

# 2) vacuum convergence on the relaxed geometry (uses relax/CONTCAR)
bash vacuum_test.sh prepare
bash vacuum_test.sh collect

# 3) vacuum level and band-edge alignment for the converged c (e.g. 20 A)
python ../09-pt111-workfunction-dipole/plot_locpot_planar.py c_20/LOCPOT --outcar c_20/OUTCAR
```

## 预期结果 / 如何检查

```bash
head -5 relax/CONTCAR                # a changed, third lattice vector still 0 0 23.13
grep "in kB" relax/OUTCAR | tail -1  # XX, YY, XY ~ 0 (ZZ is not relaxed)
grep -i -A6 "bandgap" c_20/OUTCAR | head -12   # VASP >= 6.4.3 prints VBM/CBM/gap
```

- 面内晶格常数：PBE(-D3) 一般在 3.18–3.19 Å 左右（实验约 3.16 Å；近似值）。单层的 D3 修正对 a 的影响很小。
- 真空测试：c ≥ 20 Å 时总能量变化通常已 < 1 meV/atom；E_F 的绝对数值随 c 变化，这是正常的——**只有相对真空能级的能量才有物理意义**。
- 带隙：单层 MoS2 的 PBE 带隙为 K 点直接带隙，约 1.6–1.7 eV（近似值，不含 SOC；加 SOC 后价带顶在 K 点劈裂约 0.15 eV）。
- 电离能 IP = E_vac − E_VBM，电子亲和能 EA = E_vac − E_CBM。`plot_locpot_planar.py` 会打印真空能级；VBM/CBM 取自 OUTCAR 的带隙报告或 [第5章](../../tutorials/05-band-dos.md) 的能带计算。

## 注意事项

- 单层 MoS2 上下对称，没有净偶极，因此不需要偶极校正；Janus MoSSe 等不对称二维材料必须加 `LDIPOL`/`IDIPOL`，见 [第9章](../../tutorials/09-surface-workfunction.md)。
- 二维材料的应力张量是按"含真空的整个晶胞体积"归一化的，换算面内弹性常数/应力时要乘以 c，见 [第13章](../../tutorials/13-elastic.md)。
- `LATTICE_CONSTRAINTS` 在 VASP ≤ 6.4.2 与 `IBRION=1/2` 搭配时会被静默忽略，c 轴照样会变；请检查 CONTCAR 第 5 行。

[返回第8章](../../tutorials/08-vdw-2d.md)
