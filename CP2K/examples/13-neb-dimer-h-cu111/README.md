# 示例 13：H 在 Cu(111) 上 fcc → hcp 扩散的 CI-NEB 与 Dimer 过渡态搜索

配套教程：[第13章 过渡态搜索：NEB 与 Dimer](../../tutorials/13-neb-dimer.md)

## 目的

用一个小而完整的表面反应演示 CP2K 的过渡态搜索流程：

1. 用同一套 `&FORCE_EVAL` 严格优化初态（H 在 fcc 空位）和末态（H 在 hcp 空位）；
2. 用 `RUN_TYPE BAND` + `BAND_TYPE CI-NEB` 求最小能量路径和能垒；
3. 对 TS 猜测结构做只动 H 的部分 Hessian，得到虚频模式（molden 文件）；
4. 用 Dimer（`GEO_OPT TYPE TRANSITION_STATE`）从该模式出发精修鞍点；
5. 用[第14章](../../tutorials/14-vibration-thermo.md)的 `ts_freq_check.inp` 确认只有一个虚频。

体系：3×3 表面晶胞、3 层 Cu（27 个原子，底层 9 个固定）+ 1 个 H，PBE，`DZVP-MOLOPT-PBE-GTH-q11/q1`（`BASIS_MOLOPT_UZH` + `POTENTIAL_UZH`），只算 Γ 点，对角化 + Broyden + Fermi–Dirac 500 K。

## 文件说明

| 文件 | 说明 |
|---|---|
| `cu111_h_force_eval.inc` | 所有步骤共享的 `&FORCE_EVAL`（通过 `@INCLUDE` 引入），坐标文件名由变量 `COORD` 给出 |
| `initial.xyz` / `final.xyz` | 未优化的初态（fcc）/末态（hcp）结构，原子顺序相同：1–9 底层，10–18 中间层，19–27 表层，28 为 H |
| `ts_guess.xyz` | 按几何关系放在桥位的 TS 猜测结构（没有 NEB 结果时使用） |
| `geo_opt_endpoint.inp` | 步骤 1：端点优化（`-E COORD=... -E PROJ=...` 切换初/末态） |
| `h_diffusion_cineb.inp` | 步骤 2：CI-NEB，7 个副本（含端点），`NPROC_REP` 默认 8 |
| `ts_guess_vib.inp` | 步骤 3a：只位移 H 的振动分析，写出 `ts_guess_vib-VIBRATIONS-1.mol` |
| `dimer_ts.inp` | 步骤 3b：Dimer 精修，2026.2 用 `INITIALIZATION_METHOD MOLDEN`（旧版本见文件内注释） |
| `last_frame.sh` | 取 xyz 轨迹最后一帧 |
| `neb_profile.py` | 从 BAND 主输出中提取最后一步的能量剖面、能垒，可选作图 |
| `job_slurm.sh` | Slurm 脚本，用 `STEP=` 选择步骤，并按进程数自动设置 `NPROC_REP` |

晶胞：`ABC 7.70039 7.70039 22.0`，`ALPHA_BETA_GAMMA 90 90 60`（xyz 文件的标题行也写了这一信息）。

## 运行步骤

```bash
cd CP2K/examples/13-neb-dimer-h-cu111
export CP2K_DATA_DIR=/path/to/cp2k/data     # 若已内置可省略

# 0. 语法检查（-c 只检查，不计算）
for f in geo_opt_endpoint.inp h_diffusion_cineb.inp ts_guess_vib.inp dimer_ts.inp; do
  cp2k.psmp -c "$f" > /dev/null && echo "$f OK"
done
# 注：h_diffusion_cineb.inp 需要 initial_opt.xyz/final_opt.xyz 存在后才能完整检查

# 1. 端点优化（同一输入，命令行切换坐标与项目名）
mpirun -np 32 cp2k.psmp -i geo_opt_endpoint.inp -o opt_initial.out
mpirun -np 32 cp2k.psmp -E COORD=final.xyz -E PROJ=opt_final -i geo_opt_endpoint.inp -o opt_final.out
bash last_frame.sh opt_initial-pos-1.xyz initial_opt.xyz
bash last_frame.sh opt_final-pos-1.xyz   final_opt.xyz

# 2. CI-NEB：进程数 = 7 个副本 x NPROC_REP
mpirun -np 56 cp2k.psmp -i h_diffusion_cineb.inp -o h_diffusion_cineb.out
#    只有 28 个核时：
# mpirun -np 28 cp2k.psmp -E NPROC_REP=4 -i h_diffusion_cineb.inp -o h_diffusion_cineb.out
python3 neb_profile.py h_diffusion_cineb.out --png mep.png

# 取出 TS 猜测：2026.2+ 用 FINAL 文件；旧版本取最高副本轨迹的最后一帧（这里假设是副本 4）
grep -B1 -A28 "Ener_loc_max=T" h_cu111_neb-FINAL-1.xyz > ts_from_neb.xyz
# bash last_frame.sh h_cu111_neb-pos-Replica_nr_4-1.xyz ts_from_neb.xyz

# 3a. 部分 Hessian 得到虚频模式
mpirun -np 48 cp2k.psmp -E COORD=ts_from_neb.xyz -i ts_guess_vib.inp -o ts_guess_vib.out
# 3b. Dimer 精修（CP2K 2026.2）
mpirun -np 32 cp2k.psmp -E COORD=ts_from_neb.xyz -i dimer_ts.inp -o dimer_ts.out
bash last_frame.sh h_cu111_dimer-pos-1.xyz ts_dimer.xyz

# 或者用 Slurm 逐步提交
sbatch --export=ALL,STEP=opt_initial job_slurm.sh
```

不做 NEB 也可以直接从 `ts_guess.xyz` 出发运行 3a 和 3b（去掉 `-E COORD=...` 即可）。

## 预期结果 / 如何检查

```bash
grep "GEOMETRY OPTIMIZATION COMPLETED" opt_initial.out opt_final.out dimer_ts.out
grep "ENERGY| Total FORCE_EVAL" opt_initial.out | tail -1      # E(IS)
grep "ENERGY| Total FORCE_EVAL" opt_final.out   | tail -1      # E(FS)
grep -A3 "ENERGIES \[au\]" h_diffusion_cineb.out | tail -4      # 最后一步各副本能量
grep "MAX FORCE" h_diffusion_cineb.out | tail -2               # 2026.2；旧版本在 h_cu111_neb-BAND*.out
grep "BAND TASK COMPLETED" h_diffusion_cineb.out               # 2026.2：收敛横幅
grep -l "ABORT" h_cu111_neb-BAND*.out                          # 应为空
grep "VIB|Frequency" ts_guess_vib.out                          # 3 个模式中 1 个为负
grep "DIMER|" dimer_ts.out | head -5                           # 初始化方式（MOLDEN）
grep "OPT| Local curvature" dimer_ts.out | tail -3             # 应为负值
```

- fcc 与 hcp 端点能量非常接近（两者都是三重空位），能量剖面应只有一个峰，最高点在路径中间（副本 4 附近），对应 H 越过两个表层 Cu 之间的桥位；
- `neb_profile.py` 给出正/逆向能垒和反应能；正式结论请以收敛的 k 点、层数和晶格常数为准；
- Dimer 收敛后，H 位于桥位附近，两侧 Cu–H 距离近似相等，`OPT| Local curvature` 为负；
- 最后用 `../14-vibration-thermo/ts_freq_check.inp` 做频率验证：只有一个虚频，其模式指向 fcc ↔ hcp 方向。

本仓库没有附带参考输出；能量、步数和耗时取决于 CP2K 版本、编译方式和并行设置。

## 注意事项

- 端点与 NEB/Dimer 必须使用同一个 `cu111_h_force_eval.inc` 和相同的固定原子（底层 1–9），否则能量不可比；
- slab 体系必须保留 `ALIGN_FRAMES F`、`ROTATE_FRAMES F`；
- `PLOT_REL_ENERGY`、`<PROJECT>-FINAL-1.xyz`、`INITIALIZATION_METHOD MOLDEN` 需要 CP2K ≥ 2026.2；在 2023.x–2026.1 上运行时，删掉 `dimer_ts.inp` 中 MOLDEN 相关的四行（或改写成 `&DIMER_VECTOR`，见教程 §8.3）；
- 只算 Γ 点的 3×3 金属 slab 只适合教学演示，正式计算请加 `&KPOINTS` 并检查收敛；
- `ts_guess_vib.inp` 固定了全部 Cu，只是为 Dimer 提供方向的“快速近似”，不能代替正式的 TS 频率验证；
- 不要在这些输入中加 `IGNORE_CONVERGENCE_FAILURE`：NEB 某个副本 SCF 不收敛时，应检查该副本的结构和 SCF 设置。
