# 示例 07-convergence-templates：四款软件通用的收敛测试模板

## 目的

把“建目录 → 提交 → 收集数据 → 画图判断”这一套收敛测试流程做成可复用的模板，四款软件各给一个最小算例：

| 软件 | 扫描对象 | 目录 / 文件 | 判据（默认） |
|---|---|---|---|
| VASP | `ENCUT` | `encut_scan/` | ΔE < 1 meV/atom，变胞计算 ENCUT ≥ 1.3×max(ENMAX) |
| VASP | k 点网格或 `KSPACING` | `kpoint_scan.py` | ΔE < 1 meV/atom |
| VASP | `ISMEAR`/`SIGMA` | `vasp_tsigma_check.py` | 每原子 \|T*S\| < 1 meV |
| CP2K | `CUTOFF` → `REL_CUTOFF` | `cp2k_cutoff_scan/` | ΔE < 1 meV/atom，同时看 “count for grid” 的分布 |
| OpenMX | `scf.energycutoff`、`scf.Kgrid` | `openmx_cutoff_scan/` | ΔE < 1 meV/atom |
| Gaussian | 基组 × 积分格点 | `gaussian_basis_scan/` | 关心的**相对能量**变化 < 0.5 kJ/mol（自定） |

所有扫描都输出“参数 能量 原子数”三列的表格，统一用 `plot_convergence.py` 画图并给出收敛点。

## 文件说明

| 文件 | 说明 |
|---|---|
| [encut_scan/INCAR.template](encut_scan/INCAR.template) | Si 单点能 INCAR，`@ENCUT@` 为占位符 |
| [encut_scan/POSCAR](encut_scan/POSCAR) | 金刚石 Si 原胞（2 原子） |
| [encut_scan/KPOINTS](encut_scan/KPOINTS) | 扫描 ENCUT 时固定的 Γ 中心 8×8×8 网格 |
| [encut_scan/run_encut_scan.sh](encut_scan/run_encut_scan.sh) | `prepare / run / submit / collect`；若存在 POTCAR，按 max(ENMAX) 的倍数自动生成扫描列表 |
| [kpoint_scan.py](kpoint_scan.py) | 按 n×n×n（随晶胞自动缩放）或 KSPACING 建目录，自动跳过网格相同的重复点；`collect` 汇总成 `kpoints.dat` |
| [vasp_tsigma_check.py](vasp_tsigma_check.py) | 读取 OUTCAR 中的 `EENTRO` 和 `NIONS`，判断每原子熵项 |
| [cp2k_cutoff_scan/template.inp](cp2k_cutoff_scan/template.inp) | 8 原子 Si 立方胞，`@SET CUTOFF` / `@SET REL_CUTOFF` 两行由脚本改写 |
| [cp2k_cutoff_scan/run_scan.sh](cp2k_cutoff_scan/run_scan.sh) | 两步扫描：先 `cutoff` 后 `rel`；汇总能量、各级网格上的高斯函数个数和网格电子数积分误差 |
| [openmx_cutoff_scan/template.dat](openmx_cutoff_scan/template.dat) | Si 原胞，`@ECUT@`、`@KGRID@` 占位符（OpenMX 4.0 + DFT_DATA19） |
| [openmx_cutoff_scan/run_scan.sh](openmx_cutoff_scan/run_scan.sh) | `ecut` 或 `kgrid` 扫描，从 `.out` 中提取 `Utot.` |
| [gaussian_basis_scan/template.gjf](gaussian_basis_scan/template.gjf) | HCN 与 HNC 两步 Link1 优化，`@BASIS@`、`@GRID@` 占位符 |
| [gaussian_basis_scan/run_scan.sh](gaussian_basis_scan/run_scan.sh) | 基组 × 格点组合，汇总异构化能 ΔE(HNC−HCN) |
| [slurm_array_dirs.sh](slurm_array_dirs.sh) | 通用 Slurm 作业数组：每个 task 进入 `dirs.txt` 中的一个目录运行程序 |
| [plot_convergence.py](plot_convergence.py) | 读入任意扫描表，换算为 meV/atom，打印收敛点并画 \|ΔE\| 对数图 |

## 运行步骤

以 VASP 为例（其余软件用法相同，见各脚本开头的注释）：

```bash
cd encut_scan
# 1) 自己生成 POTCAR（有版权，仓库不提供），例如 cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR
bash run_encut_scan.sh prepare                 # 建 encut_<E>/ 目录
cd ..
ls -d encut_scan/encut_* > dirs.txt            # 路径相对于提交目录
sbatch --array=1-$(wc -l < dirs.txt)%4 --export=ALL,RUN="mpirun vasp_std" slurm_array_dirs.sh
# 2) 全部结束后
bash encut_scan/run_encut_scan.sh collect
python plot_convergence.py encut_scan/encut.dat --xlabel "ENCUT (eV)" --out encut.png

# 3) 用收敛的 ENCUT 扫 k 点（另建 kscan/）
python kpoint_scan.py prepare --template encut_scan --encut 400 --mesh 4 6 8 10 12 14 16
python kpoint_scan.py collect --outdir kscan
python plot_convergence.py kpoints.dat --xlabel "N1*N2*N3" --logx
```

CP2K、OpenMX、Gaussian：

```bash
cd cp2k_cutoff_scan && bash run_scan.sh cutoff run && bash run_scan.sh cutoff collect
CUT=400 bash run_scan.sh rel run && CUT=400 bash run_scan.sh rel collect
python ../plot_convergence.py cutoff.dat rel.dat --unit hartree --labels CUTOFF REL_CUTOFF

cd ../openmx_cutoff_scan && bash run_scan.sh ecut run && bash run_scan.sh ecut collect
python ../plot_convergence.py ecut.dat --unit hartree --xlabel "scf.energycutoff (Ry)"

cd ../gaussian_basis_scan && bash run_scan.sh run && bash run_scan.sh collect
```

作业脚本也可以用本板块的 [mkjob.py](../../scripts/mkjob.py) 生成（用 `python ../../scripts/mkjob.py -h` 查看用法）；批量查看各目录是否正常结束可用 [job_status.py](../../scripts/job_status.py)（`python ../../scripts/job_status.py -h`）。

## 预期结果 / 如何检查

- `plot_convergence.py` 会逐行打印 ΔE（meV/atom），并给出“从这一点往后都在阈值以内”的最小参数值，例如（示意输出，数值仅供参考）：

  ```text
  == encut_scan/encut.dat  (reference x = 600)
             x    dE (meV/atom)
           250           12.345  *
           ...
  converged (|dE| < 1.0 meV/atom from here on): x = 400
  ```

- VASP：每个 OUTCAR 末尾应有 `General timing`；`external pressure` 随 ENCUT 增大应趋于稳定（Pulay 应力减小）。
- CP2K：输出中无 `SCF run NOT converged`；`count for grid` 不应全部集中在最粗或最细的一级网格；`Electronic density on regular grids` 的第二个数（与精确电子数之差）应随截断增大而趋近 0。
- OpenMX：`.out` 中有 `Utot.`；`.md` 末尾的 SCF 收敛标志正常。
- Gaussian：两步都出现 `Normal termination`，再比较 ΔE 随基组/格点的变化。

## 注意事项

- 扫描时**只改一个参数**，其余（PREC、EDIFF、展宽、k 点、赝势、基组）保持一致；SCF 收敛阈值要比目标精度严格两个数量级。
- 判据是 1 meV/atom 这样的**能量差**；你真正关心的量（吸附能、能垒、晶格常数、声子频率）可能需要更严格或单独测试。
- 参考点本身也是有限值：最后两个点应彼此一致，否则继续扩大扫描范围。
- 所有脚本都不会自动删除文件，但 `prepare` 会覆盖同名目录中的输入文件。
- 输入只是教学模板：Si 的 8 原子 CP2K 立方胞只算 Γ 点，适合测 CUTOFF，不适合测 k 点。

返回教程：[第7章 DFT 理论速览、学习路线与资源汇总](../../tutorials/07-learning-path.md)
