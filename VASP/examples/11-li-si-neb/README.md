# 示例 11-2：Li 间隙原子在晶体 Si 中的扩散（体相 CI-NEB）

## 目的

- 演示体相扩散 NEB 与表面 NEB 的区别：所有原子放开、三维 k 点、固定晶胞；
- 计算 Li 沿 Td → Hex → Td 路径的迁移能垒；
- 练习超胞尺寸收敛（64 → 216 原子）。

对应教程：[第11章 过渡态搜索：NEB / CI-NEB / Dimer](../../tutorials/11-neb-dimer.md)（§13）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `POSCAR_IS` | 2×2×2 立方超胞 Si64（a = 5.469 Å，PBE 近似值）+ Li 在四面体间隙位 Td，分数坐标 (0.75, 0.75, 0.75)，4 个 Si 近邻约 2.37 Å |
| `POSCAR_FS` | Li 在相邻 Td 位 (0.875, 0.875, 0.875)，跳跃距离约 2.37 Å；原子顺序与 IS 相同（Si 在前、Li 为第 65 号） |
| `INCAR_relax` | 端点优化，`ENCUT = 500`（`Li_sv`），`ISIF = 2` |
| `INCAR_neb` | CI-NEB：`IMAGES = 3`、`LCLIMB = .TRUE.`、`IOPT = 7`（FIRE）、`EDIFFG = -0.03` |
| `KPOINTS` | Γ 中心 2×2×2 |
| `POTCAR` | **不提供**；Si 用 `Si`，Li 用 `Li_sv` |

插值脚本和 NEB 作业脚本复用 [示例 11-1](../11-pt-adatom-neb/README.md) 中的 `make_images.sh`、`idpp_interpolate.py`、`job.slurm`。

## 运行步骤

```bash
cd VASP/examples/11-li-si-neb
export VASP_PP_PATH=/path/to/your/paw/library

mkdir -p IS FS
cp POSCAR_IS IS/POSCAR; cp POSCAR_FS FS/POSCAR
for d in IS FS; do cp INCAR_relax $d/INCAR; cp KPOINTS $d/; done
POT_Li=Li_sv bash ../10-co-pt111-adsorption/make_potcar.sh IS FS
grep TITEL IS/POTCAR                 # 应依次为 Si、Li_sv
cp IS/POTCAR .
(cd IS && sbatch ../../10-co-pt111-adsorption/job.slurm)
(cd FS && sbatch ../../10-co-pt111-adsorption/job.slurm)

# 端点完成后
bash ../11-pt-adatom-neb/make_images.sh IS FS 3
cp INCAR_neb INCAR
sbatch ../11-pt-adatom-neb/job.slurm           # 96 进程 = 3 图像 x 32
nebresults.pl && cat exts.dat
```

## 预期结果 / 如何检查

- 优化后 IS、FS 能量相同（两个 Td 位等价），`dist.pl IS/CONTCAR FS/CONTCAR` 约 2.4 Å 量级；
- 中间图像 02 初始位于六角间隙位（Hex），也是能量最高点；能垒曲线对称；
- 文献中 DFT 计算的能垒多在约 0.5–0.65 eV（近似范围，依赖超胞与泛函），经典实验扩散激活能约 0.66 eV（近似值）；
- 对 TS（图像 02）做频率分析（只放开 Li 及其近邻 Si），应只有一个虚频，方向沿 ⟨111⟩。

## 注意事项

- 本例路径对称、初始猜测很好，所以直接开启 CI；一般情况仍建议先普通 NEB 再 CI（见示例 11-1）。
- 所有原子都放开时，整体平移不受约束；若观察到整体漂移，可在 IS、FS 中同时固定一个远离 Li 的 Si 原子。
- 体相间隙扩散对超胞尺寸敏感：用 3×3×3 超胞（216 Si）复核能垒（教程练习 4）。
- 若改用不带半芯态的 `Li` POTCAR，可降低 ENCUT 以快速测试，但端点与 NEB 必须一致。
