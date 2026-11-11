# 通用脚本：NEP/GPUMD 机器学习势与数据集工具（NEP, GPUMD & dataset scripts）

本目录收录仓库作者 Gaoshan Zeng 与其他贡献者整理的 Python/Bash 脚本，覆盖机器学习势的整个流程：用 VASP/ABACUS 的 AIMD、单点和微扰结构生成 NEP（neuroevolution potential）训练集（dpdata、ASE、hiphive），检查 NEP 训练结果、挑选结构做主动学习，最后用 GPUMD 做热输运分析（NEMD、HNEMD、EMD、SHC 谱分解、量子修正）。这里逐个说明每个脚本的用途、依赖、命令行参数（或需要修改的变量）、输入输出和注意事项，方便你在自己的 DFT + NEP + GPUMD 项目里直接用上。

> [!IMPORTANT]
> 脚本按原样收录，多数会**在当前工作目录**读写文件，部分会删除或覆盖已有文件（下文各节的“注意”里有说明）。第一次使用请先在一个只放了输入文件的临时目录里，用小数据集试跑一遍，确认输出正确后再用到正式数据上。

---

## 1. 本目录目标与适用场景

| 你想做的事 | 推荐脚本 | 所在小节 |
|---|---|---|
| 在超胞里随机替换掺杂原子（如 MoS2 中 Mo→W） | `modeling/replace_dopant_symbols.py` | 5.1 节 |
| 生成指定曲率半径的球面原子片段（XYZ/CIF） | `modeling/fabricate_spherical_carbon.py` | 5.2 节 |
| 从 VASP AIMD 抽帧，批量建单点能目录 | `dataset/aimd_OUTCAR_xyz2POSCAR.py` | 6.1 节 |
| 用 dpdata 给结构加晶格和原子微扰 | `dataset/create-perturb.py` | 6.2 节 |
| 一次生成应变、形变、MC rattle、微扰四类结构 | `dataset/create_strain_deform_rattle_perturb.py` | 6.3 节 |
| 只扰动指定元素或指定原子 | `dataset/dpdata-system-perturb.py`（代码片段） | 6.4 节 |
| 从 GPUMD 轨迹均匀抽帧，回算 DFT（主动学习） | `dataset/dumpxyz2POSCAR.py` | 6.5 节 |
| 各帧原子数不同的轨迹逐帧回算 DFT | `dataset/messxyz2POSCAR.py` | 6.6 节 |
| ABACUS SCF 结果转 NEP 训练集 extxyz | `nep/abacus_scfoutput2nepxyz.py`、`nep/singleFrame-abacus2nep-exyz.sh` | 7.1 节、7.2 节 |
| ABACUS MD 轨迹转 NEP 训练集 | `nep/multipleFrames-abacus2nep-exyz.sh` | 7.3 节 |
| 旧版 NEP `*.in` 训练集转 extxyz，或按比例拆分训练/测试集 | `nep/inorexyz2exyz_ratio.py` | 7.4 节 |
| 画 NEP loss、能量/力/位力对角图、电荷、描述符 | `nep/plot_nep_results.py` | 7.5 节 |
| 按误差或数值区间挑出“离群”结构 | `nep/get_nep_points2strucs.py` | 7.6 节 |
| 描述符降维、最远点采样（FPS）选结构 | `nep/select_pick_structure.py` | 7.7 节 |
| NEP 声子谱和 DFT（phonopy）声子谱对比 | `nep/create_phonon_compare.py` | 7.8 节 |
| GPUMD 模型分组（NEMD、SHC） | `gpumd/split_group.py`、`gpumd/split_heterojunction_group.py` | 8.1 节、8.2 节 |
| NEMD + SHC：温度分布、热流、谱导热 G(ω) | `gpumd/plot_compute_shc.py` | 8.3 节 |
| HNEMD + SHC：κ(t)、κ(ω)、λ(ω)、κ(L) | `gpumd/plot_kappa_shc.py`、`gpumd/plot_shc_multiple.py` | 8.4 节、8.6 节 |
| 多次独立 HNEMD 求平均 | `gpumd/plot_kappa_multiple.py` | 8.5 节 |
| 谱热导率的量子修正 | `gpumd/plot_shc_quantum_correction.py` | 8.7 节 |
| EMD（Green-Kubo）多次运行平均 | `gpumd/plot_hac_multiple.py` | 8.8 节 |
| CP2K Hartree 势 cube 求真空能级（功函数） | `analysis/extract_vacuum_level_cube.py` | 9.1 节 |

相关教程：

- AIMD 与 VASP 机器学习力场（MLFF），即训练数据的上游来源：[VASP 第15章](../VASP/tutorials/15-md-mlff.md)
- 机器学习势的理论背景与学习路线：Common 第7章
- ASE、pymatgen 读写结构和解析输出：Common 第4章
- 超胞、表面、缺陷、掺杂建模：Common 第3章
- 作图与后处理通用技巧：Common 第5章
- 高通量与工作流：Common 第6章
- 在集群上提交 GPU 作业：Common 第2章

---

## 2. 目录结构

```text
Scripts/
├── modeling/    建模：掺杂替换、球面片段
├── dataset/     训练集结构生成：AIMD 抽帧、微扰、应变/rattle、轨迹回算
├── nep/         NEP 数据格式转换、训练结果分析、结构筛选、声子对比
├── gpumd/       GPUMD 分组建模与热输运后处理（NEMD/HNEMD/EMD/SHC）
└── analysis/    其他分析：cube 文件平面平均静电势
```

| 子目录 | 脚本 | 一句话功能 | 调用方式 | 主要依赖 |
|---|---|---|---|---|
| modeling | `replace_dopant_symbols.py` | 超胞中按数量随机替换元素，输出 extxyz | 改文件头变量 | ase |
| modeling | `fabricate_spherical_carbon.py` | 生成球面原子片段（XYZ + CIF） | `argparse` 参数 | numpy |
| dataset | `aimd_OUTCAR_xyz2POSCAR.py` | VASP AIMD → aimd.xyz → 抽 N 帧 → 单点目录 | `python3 脚本 N` | ase, numpy, vaspkit |
| dataset | `create-perturb.py` | dpdata 微扰 N 个结构并建立单点目录 | `python3 脚本 N POSCAR vasp`（或 `abacus`） | dpdata, vaspkit/atomkit |
| dataset | `create_strain_deform_rattle_perturb.py` | 应变、形变、MC rattle、微扰四类结构 | `python3 脚本 vasp`（或 `abacus`） | ase, dpdata, hiphive |
| dataset | `dpdata-system-perturb.py` | dpdata `perturb` 的修改片段（按元素扰动） | 不可直接运行 | dpdata |
| dataset | `dumpxyz2POSCAR.py` | 从等原子数轨迹均匀抽帧并建单点目录 | `python3 脚本 xyz N vasp`（或 `abacus`） | ase, vaspkit |
| dataset | `messxyz2POSCAR.py` | 原子数不固定的轨迹逐帧建单点目录 | `python 脚本 xyz` | ase, vaspkit |
| nep | `abacus_scfoutput2nepxyz.py` | ABACUS SCF（log/json）→ extxyz | `python 脚本 目录 输出.xyz` | ase-abacus, numpy |
| nep | `singleFrame-abacus2nep-exyz.sh` | 批量 ABACUS 单点 → NEP-dataset.xyz | `bash 脚本 目录` | bash, awk, dos2unix（ase 可选） |
| nep | `multipleFrames-abacus2nep-exyz.sh` | ABACUS MD（MD_dump）→ NEP-dataset.xyz | `bash 脚本 .`（在 MD 输出目录内运行） | bash, awk, dos2unix |
| nep | `inorexyz2exyz_ratio.py` | 旧 `*.in` → extxyz；按比例随机拆分 | `python 脚本 文件 [比例]` | numpy |
| nep | `plot_nep_results.py` | loss、对角图、分量图、电荷、描述符 PCA | 改文件头变量后直接运行 | numpy, matplotlib（ase, seaborn, scikit-learn 可选） |
| nep | `get_nep_points2strucs.py` | 按 `*_train.out` 数值或误差挑结构 | `python 脚本 模式 xyz out ...` | ase, numpy |
| nep | `select_pick_structure.py` | 描述符降维、FPS 选点、框选结构 | `python 脚本 all/select/pick 方法 ...` | ase, scikit-learn, scipy（umap-learn, NepTrain 可选） |
| nep | `create_phonon_compare.py` | 生成旧版 GPUMD 声子输入并对比 DFT 声子谱 | 改文件头变量，运行两次 | ase, numpy, matplotlib |
| gpumd | `split_group.py` | 单一材料沿 x/y/z 分组，写 model.xyz | 改文件头变量 | ase |
| gpumd | `split_heterojunction_group.py` | 多块材料拼接（带间隙）并分组 | 改文件头变量 | ase |
| gpumd | `plot_compute_shc.py` | NEMD 温度/热流 + SHC G(ω)，保存 Gc.npy | 改文件头变量 | ase, matplotlib |
| gpumd | `plot_kappa_shc.py` | 单次 HNEMD κ(t) + SHC κ(ω)、λ(ω)、κ(L) | 改文件头变量 | ase, scipy, matplotlib |
| gpumd | `plot_kappa_multiple.py` | 多次 HNEMD 运行的 κ(t) 平均 | `python 脚本 运行次数 每次步数` | scipy, matplotlib |
| gpumd | `plot_shc_multiple.py` | 多次 HNEMD-SHC 平均并求 λ(ω)、κ(L) | `python 脚本 运行次数` | matplotlib |
| gpumd | `plot_shc_quantum_correction.py` | κ 统计误差 + 谱热导率量子修正 | `python 脚本 运行次数 每次步数` | numpy ≥ 2.0, scipy |
| gpumd | `plot_hac_multiple.py` | 多次 EMD 的 HAC 与 κ(t) 平均 | 直接运行 | matplotlib |
| analysis | `extract_vacuum_level_cube.py` | cube 平面平均静电势，取最大值作真空能级 | 改文件末尾文件名 | ase, matplotlib |

---

## 3. 环境准备

### 3.1 Python 包

建议单独建一个 conda 环境（Python ≥ 3.10，因为 dpdata 1.1.0 要求 Python ≥ 3.10）：

```bash
conda create -n nep python=3.11 -y
conda activate nep
pip install numpy scipy matplotlib ase           # nearly every script
pip install dpdata hiphive                         # dataset/ perturbation and rattle scripts
pip install scikit-learn seaborn                   # plot_nep_results.py (PCA, charge histogram), select_pick_structure.py
pip install umap-learn                             # optional: select_pick_structure.py with method "umap"
pip install NepTrain                               # optional: select_pick_structure.py computes descriptors itself
```

| 包 | 用到它的脚本 | 说明 |
|---|---|---|
| numpy, matplotlib | 几乎全部 | 绘图脚本大多用 `from pylab import *` |
| scipy | `plot_kappa_shc.py`、`plot_kappa_multiple.py`、`plot_shc_quantum_correction.py`、`select_pick_structure.py` | 脚本按 SciPy 版本自动选 `cumtrapz`（< 1.14）或 `cumulative_trapezoid` |
| ase | 大部分脚本 | 读写 POSCAR、extxyz、cube、vasprun.xml |
| dpdata | `create-perturb.py`、`create_strain_deform_rattle_perturb.py`、`dumpxyz2POSCAR.py`（abacus 模式） | 微扰，以及 POSCAR → ABACUS STRU |
| hiphive | `create_strain_deform_rattle_perturb.py` | `generate_mc_rattled_structures` |
| scikit-learn | `plot_nep_results.py`（描述符 PCA）、`select_pick_structure.py` | PCA、t-SNE、Isomap 等降维 |
| umap-learn | `select_pick_structure.py umap` | 包名是 `umap-learn`，不是 `umap` |
| seaborn | `plot_nep_results.py`（电荷分布图） | 只在存在 `charge_train.out` 时才需要 |
| NepTrain | `select_pick_structure.py` | 没有现成描述符文件时，用 `Nep3Calculator` 计算每个结构的平均描述符 |
| ase-abacus | `abacus_scfoutput2nepxyz.py` | `pip install git+https://gitlab.com/1041176461/ase-abacus.git` |

> [!WARNING]
> `ase-abacus` 是 ASE 的分支（fork），包名同样叫 `ase`，安装后会**替换**环境里原有的 ASE。建议给 ABACUS 数据处理单独建一个环境，免得影响其他脚本。

### 3.2 外部程序

| 程序 | 用途 | 用到的脚本 |
|---|---|---|
| GPUMD（`nep`、`gpumd`，可选 `gnep`） | 训练 NEP、跑 MD 和热输运 | `nep/`、`gpumd/` 的输入输出 |
| VASPKIT | `vaspkit -task 102 -kpr 0.03`：按 $2\pi\times0.03$ Å⁻¹ 的 k 点分辨率生成 KPOINTS，同时自动生成 POTCAR | `dataset/` 中的 VASP 模式 |
| ATOMKIT | VASPKIT 团队的多软件建模工具；脚本用它把 POSCAR 转成带赝势/轨道信息的 STRU，并生成 ABACUS 的 KPT | `dataset/` 中的 ABACUS 模式 |
| dos2unix | 统一换行符 | 两个 `.sh` 脚本最后一步 |
| phonopy | `phonopy-bandplot --gnuplot > phonon.out` 导出 DFT 声子谱 | `create_phonon_compare.py` |

VASPKIT 需要先在 `~/.vaspkit` 里配好赝势库路径，否则自动生成 POTCAR 会失败。ATOMKIT 同理，需要在 `~/.atomkit` 中配置 ABACUS 赝势和数值轨道路径。

### 3.3 如何调用

脚本基本都只操作当前目录的文件，推荐做法是进入数据目录，再用绝对路径调用：

```bash
export DFT_SCRIPTS=/path/to/DFT/Scripts      # adjust to your clone
cd /work/my_nep_project/train_01
python $DFT_SCRIPTS/nep/plot_nep_results.py
```

需要修改文件头变量的脚本（表中注明“改文件头变量”的那些），请先**复制**到工作目录再改，不要直接改仓库里的版本：

```bash
cp $DFT_SCRIPTS/gpumd/split_group.py .
vi split_group.py      # edit cx, cy, cz, dic, group_cyclical
python split_group.py
```

在没有图形界面的集群节点上，用 `pylab` 的脚本会自动退回非交互后端，只保存 PNG。`analysis/extract_vacuum_level_cube.py` 是例外：它显式调用了 `matplotlib.use('TkAgg')`，需要显示环境（见 9.1 节）。

### 3.4 GPU 作业模板（Slurm）

NEP 训练和 GPUMD 模拟都在 GPU 上运行。下面是一个最小的 Slurm 模板（分区名、模块名要按你的集群改；更多写法见 Common 第2章）：

```bash
#!/bin/bash
#SBATCH --job-name=nep-train
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --output=%x-%j.out

module load cuda                       # cluster specific
export PATH=/path/to/GPUMD/src:$PATH   # directory containing nep and gpumd

cd $SLURM_SUBMIT_DIR
nep > nep.log 2>&1                     # reads nep.in, train.xyz (test.xyz) in this directory
# gpumd > gpumd.log 2>&1               # reads run.in, model.xyz, nep.txt
```

---

## 4. 典型工作流

```text
 (1) DFT 数据                    (2) 训练集                 (3) 训练与检查                    (4) 应用
 VASP AIMD ──aimd_OUTCAR_xyz2POSCAR──┐
 微扰/应变 ──create-perturb / create_strain_...──> 单点 DFT ──┐
 ABACUS SCF ──abacus_scfoutput2nepxyz / singleFrame.sh ───────┼─> train.xyz ─> nep ─> plot_nep_results
 ABACUS MD  ──multipleFrames.sh ──────────────────────────────┘   test.xyz          get_nep_points2strucs
                                                                                    select_pick_structure
                                                                                    create_phonon_compare
                                                                                          │
         dumpxyz2POSCAR / messxyz2POSCAR <── gpumd MD 轨迹（主动学习，回算 DFT） <────┘
                                                                                          │
                     split_group / split_heterojunction_group ─> gpumd（NEMD/HNEMD/EMD）─> gpumd/plot_*.py
```

1. **初始数据**：VASP AIMD（见 [VASP 第15章](../VASP/tutorials/15-md-mlff.md)）或对原型结构做微扰、应变，再算单点能。单点 INCAR 的写法见 6.0 节。
2. **转格式**：把 DFT 结果转成 NEP 的 extxyz（必须有 `energy`、`lattice`、`properties=species:S:1:pos:R:3:force:R:3`，`virial` 可选），合并成 `train.xyz`，可以用 `inorexyz2exyz_ratio.py` 随机划出 `test.xyz`。
3. **训练与检查**：`nep` 训练后用 `plot_nep_results.py` 看 loss 和对角图；用 `get_nep_points2strucs.py` 或 `select_pick_structure.py` 删掉离群结构、给数据集去冗余。
4. **主动学习**：用 GPUMD 跑 MD，从轨迹抽帧回算 DFT，加入训练集后重新训练，反复迭代直到收敛。
5. **应用**：用 `split_group.py` 给模型分组，跑 NEMD/HNEMD/EMD，再用 `gpumd/plot_*.py` 分析热导率。

> [!NOTE]
> NEP 的 extxyz 格式：第二行的关键字不区分大小写；位力（virial）单位是 eV，指整个晶胞而不是每个原子；`virial` 和 `stress` 同时存在时以 `virial` 为准；力的列名写 `force:R:3` 或 `forces:R:3` 都可以。详见 [train.xyz 格式文档](https://gpumd.org/nep/input_files/train_test_xyz.html)。

---

## 5. modeling/ 建模脚本

### 5.1 replace_dopant_symbols.py：随机掺杂替换

**用途**：读入原胞 `POSCAR`，扩成 `cx×cy×cz` 超胞，按设定数量把目标元素随机替换为一种或几种掺杂元素，可以一次生成多个随机构型，输出 extxyz。适合给 NEP 训练集准备合金或掺杂构型，也适合准备 SQS 之前的随机初猜。

**依赖**：`ase`、`numpy`。

**用法**：没有命令行参数，修改文件开头的变量后运行 `python replace_dopant_symbols.py`。

| 变量 | 默认值 | 说明 |
|---|---|---|
| `periodic_order` | `0` | `0`：在整个超胞中随机选位；`1`：把替换数平均分到每个复制单元（每个原胞拷贝替换数相同）。取 `1` 时，`num_replacements_list` 的每个子列表只能有一个数，且总替换数要能被 `cx*cy*cz` 整除 |
| `num_structures` | `1` | 生成的随机构型数 |
| `uc = read('POSCAR')` | `POSCAR` | 原胞文件，ASE 能读的格式都行 |
| `cx, cy, cz` | `4, 2, 1` | 扩胞倍数 |
| `target_symbols` | `['Mo']` | 被替换的元素列表 |
| `replacement_symbols_list` | `[['W']]` | 每个目标元素对应的掺杂元素列表，例如 `[['W'], ['Se', 'Te']]` |
| `num_replacements_list` | `[[8]]` | 每种掺杂元素的替换个数，形状与上一项一致，例如 `[[8], [4, 4]]` |
| `extend` | `None` | 可选，按位点种类限制替换。用 ASE 读带对称信息的 CIF 时，`atoms.arrays['spacegroup_kinds']` 会标出 Wyckoff 位点种类；令 `extend = atoms.arrays.get('spacegroup_kinds')` 即可启用 |
| `extend_infos` | `[[2]]` | `extend` 不为 `None` 时生效，表示每个目标元素只在这些 `spacegroup_kinds` 编号的位点上替换 |
| `out_file` | `replace_structures.xyz` | 输出文件 |

**输入/输出**：输入 `POSCAR`；输出 `replace_structures.xyz`，含 `num_structures` 帧。

**注意**：

- 脚本固定了 `random.seed(42)`，每次运行结果相同。要得到不同构型，可以改种子或加大 `num_structures`。
- 替换数大于目标元素原子数时会报 `ValueError` 并说明原因。
- 输出的 extxyz 没有能量和力，只用于建模。要拿去 DFT 计算，可以再交给 6.6 节的 `messxyz2POSCAR.py` 拆成单点目录。

### 5.2 fabricate_spherical_carbon.py：球面原子片段生成

**用途**：生成具有给定曲率半径的“球冠”状原子片段，例如研究曲率对碳材料吸附或催化的影响。做法分四步：先在平面上生成交错排列的六角形点阵（列间距 $\sqrt3\,a$，行间距 $1.5\,a$，$a$ 为 `--bond_length`）；在圆盘内截取点阵，并随机抽样到目标原子数；再把每个点竖直投影到以原点为球心、半径为 $R$ 的球面上（$z=\sqrt{R^2-x^2-y^2}$）；原子数不够时逐步扩大圆盘补点。

**依赖**：`numpy`（只用标准库 + numpy）。

**用法**：

```bash
python fabricate_spherical_carbon.py --radius 40 --num_atoms 51 --element C --bond_length 1.42
python fabricate_spherical_carbon.py --radius 500 --num_atoms 300 --element Si --bond_length 2.35
```

| 参数 | 代码中的默认值 | 说明 |
|---|---|---|
| `--radius` | `40.0` | 球面曲率半径（Å） |
| `--num_atoms` | `51` | 原子数 |
| `--element` | `C` | 元素符号 |
| `--bond_length` | `1.42` | 点阵特征长度 $a$（Å） |

`--help` 里显示的默认值说明文字以代码中的 `default=` 为准。

**输出**：`{元素}{原子数}_spherical_r{半径}.xyz` 和同名 `.cif`；屏幕打印平均球心距离（即实际曲率半径）、曲率误差、在 $a$ 的 ±10% 范围内的近邻距离统计、三个方向的尺寸，原子数 ≤ 1000 时还会打印 ASCII 俯视图。

**注意**：

- `save_xyz_file()` 和 `save_cif_file()` 的默认输出文件名带有作者本机的目录前缀，使用前请改成你自己的路径（或直接改成只有文件名）。
- 抽样用的是 `np.random` 且没有设种子，所以每次结果不同。
- 球心在原点，球冠位于 $z\approx R$ 附近；CIF 的立方晶胞边长取 2.2 倍最大坐标，R 很大时晶胞会非常大，做周期性计算前请先平移结构，再按需要设置真空层。
- 生成的只是几何初猜。用于 DFT 或 MD 前，请先用 VESTA/OVITO 检查近邻距离和配位，再做结构优化。

---

## 6. dataset/ 训练集结构生成

### 6.0 单点能模板（INCAR-scf / INCAR-single / INPUT-scf）

`dataset/` 下的脚本会把一个“单点能模板”复制到每个子目录。不同脚本要求的模板文件名不一样：

| 脚本 | VASP 模板 | ABACUS 模板 |
|---|---|---|
| `aimd_OUTCAR_xyz2POSCAR.py` | `INCAR-single` | 不支持 |
| `messxyz2POSCAR.py` | `INCAR-single` | 不支持 |
| `create-perturb.py` | `INCAR-scf` | `INPUT-scf` |
| `create_strain_deform_rattle_perturb.py` | `INCAR-scf` | `INPUT-scf` |
| `dumpxyz2POSCAR.py` | `INCAR-scf` | `INPUT-scf` |

一个可以参考的 VASP 单点模板如下（示例，ENCUT、ISMEAR、自旋等要按体系调整；整个训练集必须用**同一套**参数）：

```text
SYSTEM = NEP single point
PREC   = Accurate
ENCUT  = 520          # use one value for the whole dataset (>= 1.3 x max ENMAX helps stress)
EDIFF  = 1E-6
NELM   = 120
ISMEAR = 0
SIGMA  = 0.05
IBRION = -1           # no ionic update
NSW    = 0
ISIF   = 2            # stress tensor is computed and written
LREAL  = .FALSE.
ISYM   = 0
LWAVE  = .FALSE.
LCHARG = .FALSE.
```

ABACUS 的 `INPUT-scf` 至少要包含 `calculation scf`、`cal_force 1`、`cal_stress 1`（否则没有位力，脚本会提示 `please add cal_stress in INPUT`），并设置合理的 `scf_nmax`。转换脚本靠电子步数是否达到 `scf_nmax` 判断未收敛，见 7.1 节。

### 6.1 aimd_OUTCAR_xyz2POSCAR.py：VASP AIMD 抽帧并建立单点目录

**来源**：脚本注释说明，读取 vasprun.xml 并写 extxyz 的部分使用了知乎作者“叠加态”的代码（https://zhuanlan.zhihu.com/p/397866915）。

**用途**：递归查找当前目录下所有 `vasprun.xml`（读取失败时改读同目录的 `OUTCAR`），把所有离子步写进 `aimd.xyz`（含能量、力、位力 $\mathbf W=-\boldsymbol\sigma V$，`config_type=aimd`）；再均匀抽出 N 帧写成 `dump.xyz`；每帧建一个单点目录，放入 POSCAR、INCAR（复制自 `INCAR-single`），并调用 VASPKIT 生成 KPOINTS 和 POTCAR。

**依赖**：`ase`、`numpy`；外部命令 `find`、`grep`、`awk`、`vaspkit`。

**用法**：

```bash
cd aimd_runs/                # contains run1/vasprun.xml, run2/vasprun.xml, ... and INCAR-single
python3 $DFT_SCRIPTS/dataset/aimd_OUTCAR_xyz2POSCAR.py 100
```

参数 `100` 是要抽取的帧数。抽取方式：设总帧数为 $M$，取第 $\lfloor M/N\rfloor\cdot(i+1)$ 帧（$i=0,\dots,N-1$，从 1 开始计数），也就是每个等长区间的最后一帧。

**输出**：

```text
aimd.xyz                         all ionic steps with energy/forces/virial (extxyz)
dump.xyz                         N sampled frames
train_folders/aimd-1/aimd-1.xyz
train_folders/aimd-1/POSCAR
train_folders/aimd-1/INCAR       copied from INCAR-single
train_folders/aimd-1/KPOINTS     from vaspkit -task 102 -kpr 0.03
train_folders/aimd-1/POTCAR      from vaspkit
...
```

**注意**：

- 脚本启动时会**删除**当前目录已有的 `aimd.xyz`、`dump.xyz` 和 `train_folders/`。注意 GPUMD 的轨迹也常叫 `dump.xyz`，不要在放着 GPUMD 结果的目录里运行。
- 抽帧按“每帧行数 = 第一帧原子数 + 2”切分，要求所有帧原子数相同（AIMD 满足这一点）。
- 写入训练集前，建议检查 AIMD 中有没有电子步达到 NELM 上限的离子步：`grep -B1 E0 OSZICAR | grep -E 'DAV|RMM'` 的第 2 列就是该离子步最后一次电子迭代的步数。
- 如果只想要 `aimd.xyz`（AIMD 直接作训练数据），跑完后保留它即可，`train_folders/` 可以删掉。

### 6.2 create-perturb.py：dpdata 微扰结构批量生成

**用途**：用 `dpdata.System.perturb` 对一个结构生成 N 个“晶格 + 原子坐标”同时扰动的结构，每个结构一个目录，并准备好 VASP 或 ABACUS 单点输入。

**依赖**：`dpdata`；VASP 模式需要 `vaspkit`，ABACUS 模式需要 `atomkit`。

**用法**：

```bash
python3 $DFT_SCRIPTS/dataset/create-perturb.py 50 POSCAR vasp      # needs INCAR-scf in cwd
python3 $DFT_SCRIPTS/dataset/create-perturb.py 50 POSCAR abacus    # needs INPUT-scf in cwd
```

| 位置参数 | 含义 |
|---|---|
| 1 | 微扰结构数 N |
| 2 | 结构文件名（VASP POSCAR 格式） |
| 3 | `vasp` 或 `abacus` |

扰动参数直接写在 `perturb()` 调用里，要改请编辑脚本：

| 参数 | 脚本中的值 | 含义 |
|---|---|---|
| `cell_pert_fraction` | `0.03` | 晶格扰动幅度：对角分量在 ±3% 内均匀取值，非对角分量在 ±1.5% 内取值 |
| `atom_pert_distance` | `0.1` | 原子位移幅度（Å） |
| `atom_pert_style` | `'uniform'` | `uniform`：在半径 0.1 Å 的球内均匀取点；`normal`：正态分布；`const`：固定距离 |
| `atom_pert_prob` | `1.0` | 被扰动原子的比例 |

**输出**（以 `POSCAR` 为例）：

```text
perturb-POSCAR/POSCAR
perturb-POSCAR/train-POSCAR-1/{POSCAR, INCAR, KPOINTS, POTCAR}            # vasp
perturb-POSCAR/train-POSCAR-1/{POSCAR, STRU, INPUT, KPT}                  # abacus
...
```

**ABACUS 模式细节**：脚本先让 atomkit 把原始 POSCAR 转成 `<POSCAR>.STRU`，从中截取前 `2*ntype+3` 行（`ATOMIC_SPECIES` 和 `NUMERICAL_ORBITAL` 两段，含赝势和轨道文件名），替换掉 dpdata 写出的每个 STRU 的开头部分；最后用 atomkit 按 k 点间距 0.03 生成 KPT。所以 atomkit 必须事先配置好赝势和轨道路径。

**注意**：微扰幅度要和你关心的温区匹配。0.03/0.1 Å 适合近平衡构型；高温或相变研究可以加大扰动，或配合 AIMD 数据。

### 6.3 create_strain_deform_rattle_perturb.py：四类扰动结构一键生成

**用途**：对多个原型结构（默认 `POSCAR1`～`POSCAR5`）一次生成四类训练结构：

| 类别 | 生成方式 | 每个原型的数量 | 输出目录 |
|---|---|---|---|
| strain（应变） | 沿笛卡尔 x、y、z 分别施加正应变 $\varepsilon_i\in U(-0.05,0.05)$（晶格矩阵右乘 $\mathrm{diag}(1+\varepsilon_i)$），原子按分数坐标随晶胞缩放 | 10 | `strain_deform/strain-deform-k/` |
| deform（形变） | 晶格矩阵左乘 $\mathbf I+\mathbf R$，$\mathbf R$ 为 3×3 随机矩阵，元素 $\in U(-0.05,0.05)$，相当于对晶格矢量做随机线性组合 | 10 | 同上（与 strain 交替编号） |
| rattle | 先做 strain，再用 hiphive 的 Monte Carlo rattle（`rattle_std=0.03`、`d_min=1.5`、`n_iter=5`） | 20 | `rattle/<名>_size_<扩胞>_iter_<i>/` |
| perturb | `dpdata.perturb`（`cell_pert_fraction=0.05`、`atom_pert_distance=0.2`、`uniform`） | 20 | `perturb-<i>/perturb-<i>-<j>/` |

脚本注释说明，hiphive 直接写 POSCAR 时作者遇到过问题，所以先写 extxyz 再用 ASE 转成 POSCAR。MC rattle 的典型位移约为 $\sqrt{n_\mathrm{iter}}\times$`rattle_std`。

**依赖**：`ase`、`numpy`、`dpdata`、`hiphive`；外部 `vaspkit` 或 `atomkit`。

**用法**：

```bash
ls     # POSCAR1 POSCAR2 POSCAR3 POSCAR4 POSCAR5 INCAR-scf   (or INPUT-scf for abacus)
python3 $DFT_SCRIPTS/dataset/create_strain_deform_rattle_perturb.py vasp
```

需要修改的变量：

| 变量 | 默认 | 说明 |
|---|---|---|
| `prototype_structures['1'..'5']` | `read('POSCAR1')`…`read('POSCAR5')` | 原型结构，增删时下面三处要一起改 |
| `strain_lim` | `[-0.05, 0.05]` | 应变/形变范围 |
| `n_structures` | 第一段 `10`，rattle 段 `20` | 每个原型的生成数 |
| `rattle_std`, `d_min`, `n_iter` | `0.03`, `1.5`, `5` | MC rattle 参数（Å） |
| `size_vals['1'..'5']` | `[(1,1,1)]` | rattle 前的扩胞倍数，可以给多个 |
| `for i in range(1, 6)` | 5 个原型 | perturb 段写死了原型个数，需要与 `prototype_structures` 一致 |
| `num_perturb` 及 `perturb()` 参数 | `20`、`0.05`、`0.2` | perturb 段参数 |

**输出**：上表中的目录，每个子目录有 POSCAR（ABACUS 模式另有 STRU），以及复制进去的 INCAR/INPUT 和生成的 KPOINTS/KPT。VASP 模式用 `vaspkit -task 102 -kpr 0.04`，ABACUS 模式用 atomkit，k 点间距都是 0.04。rattle 目录名中的括号和空格会被去掉，例如 `1_size_1,1,1_iter_0`。

> [!CAUTION]
> 脚本最后会遍历**当前目录下所有子目录的所有子目录**，把 `INCAR-scf`（或 `INPUT-scf`）复制进去，并在每个目录里运行 vaspkit/atomkit。务必在一个只有 `POSCAR1..5` 和模板文件的干净目录里运行，否则其他目录也会被写入文件。

### 6.4 dpdata-system-perturb.py：dpdata `perturb` 的修改片段

**来源**：摘自 dpdata（https://github.com/deepmodeling/dpdata ，LGPL-3.0）`System.perturb` 方法并做了修改，新增参数 `elem_pert_list`，只扰动指定元素（如 `['O', 'Si']`）或指定原子序号（如 `[2, 34]`）。

**用途**：需要“只让某些原子动”（例如只扰动吸附分子、只扰动掺杂原子）时参考。

**用法**：这个文件只是带缩进的方法体片段，**不能直接运行或导入**。目前上游 dpdata（v1.1.0）的 `perturb` 签名是 `perturb(pert_num, cell_pert_fraction, atom_pert_distance, atom_pert_style="normal", atom_pert_prob=1.0)`，没有 `elem_pert_list` 参数。要用这个功能，需要自己把逻辑整合进一个独立函数（也可以参考 Benrui Tang 在 https://github.com/tang070205/tools 中维护的同名文件，其中已改写为可独立导入的函数）。修改或再分发 dpdata 代码时要遵守 LGPL-3.0。

### 6.5 dumpxyz2POSCAR.py：GPUMD 轨迹抽帧（主动学习）

**用途**：主动学习的一步。用 NEP 跑完 GPUMD MD 后，从轨迹（GPUMD `dump_xyz` 等关键字输出的 extxyz，要求每帧原子数相同）中均匀抽出 N 帧，写到 `activate-learning.xyz`，再为每帧建立 DFT 单点目录。

**依赖**：`ase`；ABACUS 模式还需要 `dpdata` 和 `atomkit`；VASP 模式需要 `vaspkit`。

**用法**：

```bash
python3 $DFT_SCRIPTS/dataset/dumpxyz2POSCAR.py dump.xyz 50 vasp      # needs INCAR-scf
python3 $DFT_SCRIPTS/dataset/dumpxyz2POSCAR.py dump.xyz 50 abacus    # needs INPUT-scf and a POSCAR (for the STRU header)
```

| 位置参数 | 含义 |
|---|---|
| 1 | 轨迹文件（extxyz） |
| 2 | 抽取帧数 N（抽法与 6.1 节相同） |
| 3 | `vasp` 或 `abacus` |

**输出**：`activate-learning.xyz` 和 `train_folders/activate-learning-<i>/`（内含 xyz、POSCAR、INCAR、KPOINTS、POTCAR）。

**注意**：`train_folders/` 已存在时，`os.makedirs` 会报 `FileExistsError`，重跑前要先删掉或改名。回算结果的转换方法见 7.1～7.3 节（ABACUS），VASP 可用 [GPUMD tools 中的 vasp2xyz](https://github.com/brucefan1983/GPUMD/tree/master/tools/Format_Conversion/vasp2xyz)。

### 6.6 messxyz2POSCAR.py：原子数不固定的轨迹逐帧拆分

**用途**：输入一个各帧原子数不同的 extxyz（例如沉积、吸附、反应轨迹，或 5.1 节生成的多个掺杂构型），用 ASE 逐帧读取并拆分，**每一帧**都建一个 VASP 单点目录。

**依赖**：`ase`；外部 `vaspkit`。

**用法**：

```bash
python $DFT_SCRIPTS/dataset/messxyz2POSCAR.py candidates.xyz       # needs INCAR-single in cwd
```

**输出**：`mess_strucs_folders/activate-learning-<i>/{activate-learning-<i>.xyz, POSCAR, INCAR, KPOINTS, POTCAR}`。

**注意**：不做抽样，帧多的时候请先用 7.7 节的 FPS 筛一遍。`mess_strucs_folders/` 已存在时同样要先删掉。

---

## 7. nep/ NEP 数据转换与训练结果分析

### 7.1 abacus_scfoutput2nepxyz.py：ABACUS SCF 结果转 extxyz

**来源**：文件头标注作者 Gaoshan Zeng。GPUMD 官方仓库 [tools/Format_Conversion/abacus2xyz](https://github.com/brucefan1983/GPUMD/tree/master/tools/Format_Conversion/abacus2xyz) 和 [tang070205/tools](https://github.com/tang070205/tools) 中有同名脚本，除文件头外内容基本相同。

**用途**：递归遍历目录，找到 ABACUS SCF 输出，提取能量、晶格、坐标、力和位力，写成 NEP extxyz。支持两种输出：

| 找到的文件 | 读取方式 | 位力换算 |
|---|---|---|
| `running_scf.log` | ase-abacus 的 `read(..., format='abacus-out')` | $\mathbf W=-\boldsymbol\sigma_\mathrm{ASE}\,V$（eV） |
| `abacus.json`（同目录没有 log 时） | 直接解析 JSON 的 `init`/`output` 字段 | $\mathbf W=\boldsymbol\sigma_\mathrm{kbar}\,V/1602.1766208$（eV） |

收敛判断：统计 log 中含 `ALGORITHM` 的行数（即电子步数），等于 `scf_nmax` 时视为未收敛；没有 `Total  Time`（计算未正常结束）也跳过。`scf_nmax` 从同目录的 `INPUT` 读取（ABACUS 会在 `OUT.<suffix>/` 下写一份完整 INPUT），找不到 INPUT 时按 100 处理。

**依赖**：`numpy`；读取 log 需要 ase-abacus（见 3.1 节）。

**用法**：

```bash
python $DFT_SCRIPTS/nep/abacus_scfoutput2nepxyz.py ./scf_runs abacus_scf.xyz
```

| 位置参数 | 含义 |
|---|---|
| 1 | 要递归搜索的根目录 |
| 2 | 输出 extxyz 文件名 |

**输出**：每个结构的注释行含 `energy=`、`Lattice=`、`Virial=`（有应力时）、`config_type=<目录绝对路径>`、`Properties=species:S:1:pos:R:3:forces:R:3`。

**注意**：

- 判断收敛所依赖的日志关键词（`ALGORITHM`、`Total  Time`）随 ABACUS 版本可能变化，换新版本后请先用 `grep -c ALGORITHM running_scf.log` 确认能匹配上。
- 文件末尾有一段被注释掉的 ASE 写法（`write(..., format='extxyz', append=True)`），作为另一种写出方式保留。
- 转换完请检查结构数：`grep -c Lattice abacus_scf.xyz` 应等于收敛的计算数。

### 7.2 singleFrame-abacus2nep-exyz.sh：ABACUS 单点批量转 extxyz

**来源**：脚本头部注明 `Modified by Benrui Tang`（tang070205@proton.me）。GPUMD 仓库 tools/Format_Conversion/abacus2xyz/ 收录了几乎相同的版本。这一系列 Bash 脚本沿用了 GPUMD tools 中 Yanzhou Wang 等人 `outcar2nep-exyz` 系列脚本的写法（GPUMD 收录的 multipleFrames 版本头部注明改编自 `multipleFrames-outcars2nep-exyz.sh`）。

**用途**：在指定目录下递归查找（会跟随软链接）所有 `running_scf.log`，跳过未收敛的计算（电子步数等于 `scf_nmax`，或日志里没有 `FINAL_ETOT_IS`），把能量、晶格、位力、元素、坐标和力写入 `NEPdataset/NEP-dataset.xyz`。

**依赖**：bash、grep、awk、paste、dos2unix；装了 ASE 时，坐标和晶格从 log 同目录的 `STRU.cif` 读取；没有 ASE 时改为解析上一级目录的 `STRU`，以及 log 中的 `taud`/`tauc` 坐标。

**用法**：

```bash
bash $DFT_SCRIPTS/nep/singleFrame-abacus2nep-exyz.sh ./scf_runs
```

| 可改变量（文件头） | 默认 | 说明 |
|---|---|---|
| `isol_ener` | `0` | 每个原子要减去的能量（eV），输出 $E-N\cdot$`isol_ener`；对所有元素用同一个值 |

**输出**：`NEPdataset/NEP-dataset.xyz`（运行前会删除旧的 `NEPdataset/`），`Config_type` 取输出目录路径的最后三级，`Weight=1.0`。

**注意**：

- 有 `TOTAL-STRESS` 时，按 $\mathbf W=\boldsymbol\sigma_\mathrm{kbar}V/1602.1766208$ 写 Virial；没有时打印 `Warning: No virial found`。
- 建议在运行环境中安装 ASE，让脚本走读取 `STRU.cif` 的分支。脚本用 `pip list | grep ase` 判断 ASE 是否可用。

### 7.3 multipleFrames-abacus2nep-exyz.sh：ABACUS MD 轨迹转 extxyz

**来源**：与 7.2 节相同（Modified by Benrui Tang，GPUMD tools 收录）。

**用途**：读取 ABACUS MD 的 `running_md.log`（能量 `etot`、原子数、每步电子步数）和 `MD_dump`（晶格、坐标、力，以及可选的 `VIRIAL (kbar)`），逐帧写入 `NEPdataset/NEP-dataset.xyz`，电子步数等于 `scf_nmax` 的帧会被跳过。

**用法**：脚本读取的是**当前工作目录**下的 `running_md.log`、`MD_dump` 和 `INPUT`，所以请先进入 MD 输出目录（通常是 `OUT.<suffix>/`，ABACUS 会在其中写一份完整 INPUT）再运行。命令行参数不能为空，给 `.` 即可：

```bash
cd md_run/OUT.ABACUS
bash $DFT_SCRIPTS/nep/multipleFrames-abacus2nep-exyz.sh .
```

**输出**：`NEPdataset/NEP-dataset.xyz`，`Config_type=<路径后三级>-<帧号>`。

**注意**：

- 循环从第 2 个 `MDSTEP` 块开始（索引 `i=1`）。转换后请抽查几帧，确认能量和坐标对得上，例如对比 `grep etot running_md.log` 的输出与 xyz 中的 `Energy=`。
- 日志关键词（`STEP OF MOLECULAR DYNAMICS`、`MDSTEP`、`LATTICE_VECTORS`、`INDEX`）随 ABACUS 版本可能变化。

### 7.4 inorexyz2exyz_ratio.py：旧格式转换与按比例拆分

**用途**：

1. 把早期 GPUMD/NEP 使用的 `*.in` 训练集格式（第一行是结构数 N，接着 N 行各以该结构的原子数开头，之后依次是每个结构的能量（及位力）行、9 个数的晶格行和原子行）转换为现在的 extxyz。
2. 转换时或对现有 extxyz 按比例随机抽取一部分（例如划出 test.xyz）。

**依赖**：`numpy`（只用标准库 + numpy）。

**用法**：

```bash
python inorexyz2exyz_ratio.py train.in            # -> train.xyz
python inorexyz2exyz_ratio.py train.in 0.1        # -> train.xyz, train_extract.xyz (10 %), train_surplus.xyz (90 %)
python inorexyz2exyz_ratio.py all.xyz 0.1         # -> all_extract.xyz, all_surplus.xyz
```

**输出**：`<前缀>.xyz`、`<前缀>_extract.xyz`、`<前缀>_surplus.xyz`。写出的注释行格式为 `Energy=... Virial="9 个分量" Lattice="..." Config_type=infile-<k> Weight=1 Properties=species:S:1:pos:R:3:force:R:3 pbc="T T T"`。

> [!CAUTION]
> `.in` 模式有三个副作用：(1) 会把输入的 `.in` 文件**原地改写**，所有小数统一格式化为 8 位；(2) 结束时会**删除当前目录下所有 `*.out` 文件**（包括 `loss.out`、`energy_train.out` 等 NEP 输出）；(3) 输出的 `<前缀>.xyz` 会覆盖同名文件（`train.in` 对应 `train.xyz`）。请把 `.in` 文件复制到一个空目录里再转换。

**补充**：随机抽取使用 `random`，没有设种子，所以每次结果不同；需要可复现时请在脚本开头加 `random.seed(...)`。转换后建议用 ASE 读一遍输出（`ase.io.read('train.xyz', ':')`），并抽查几个结构的能量、位力和原子数。

### 7.5 plot_nep_results.py：NEP 训练与预测结果作图

**用途**：在 NEP 训练目录（或 `prediction 1` 的预测目录）中直接运行，自动识别文件并作图：

| 情形 | 判断依据 | 输出图片 |
|---|---|---|
| 训练（`nep`） | 存在 `loss.out`、`nep.in` | `nep-loss.png` |
| 训练（`gnep`） | 存在 `loss.out`、`gnep.in` | `gnep-loss-learning_rate.png`（loss + 学习率） |
| 能量/力对角图 | `lambda_v 0` 或所有结构都没有位力 | `nep-ef-diagonals.png` |
| 能量/力/位力 | 有 `virial_train.out`、无 `stress_train.out` | `nep-efv-diagonals.png` |
| 能量/力/位力/应力 | 还有 `stress_train.out` | `nep-efvs-diagonals.png` |
| 偶极/极化率模型 | 存在 `dipole*` 或 `polarizability*` 文件 | `nep-dipole-diagonal.png` 或 `nep-polarizability-diagonal.png` |
| qNEP 电荷 | 存在 `charge_train.out` | `nep-charge.png`（需 seaborn） |
| 描述符 | 预测模式且存在 `descriptor.out` | `nep-descriptor.png`（PCA，需 scikit-learn） |
| 分元素力 | `element_force = 1` | `<元素>-force.png` |
| 分量图 | `component = 'force'` 等 | `<分量>-train-components.png`、`<分量>-test-components.png` |

图例中给出 RMSE 和 $R^2$。RMSE < 1 时自动换成 meV/atom、meV/Å、MPa 等单位。`lambda_v` 和 `batch` 从 `nep.in`/`gnep.in` 读取；文件中有 `prediction` 关键字时视为全批次（full batch）。位力文件里目标值为 `-1e+06` 的结构（GPUMD 用 −1e6 标记“无位力参考”，脚本按字符串 `-1e+06` 识别）不参与位力对角图，其行号会存到 `train_no_virial_indices.txt` 和 `test_no_virial_indices.txt`。

**依赖**：`numpy`、`matplotlib`；分元素力、电荷、逐原子描述符需要 `ase`（读 `train.xyz`/`test.xyz` 和 `nep.txt` 首行的元素列表）；电荷图需要 `seaborn`；描述符图需要 `scikit-learn`。

**用法**：先改文件开头的变量，再在 NEP 目录中运行 `python plot_nep_results.py`。

| 变量 | 默认 | 说明 |
|---|---|---|
| `three_six_component` | `0` | `1`：对角图中力/偶极按 3 个分量、位力/应力/极化率按 6 个分量分色 |
| `use_range` | `0` | `0`：坐标范围取数据最小/最大值；`1` 或 `2`：用 `plot_range` 中的范围 |
| `element_force` | `0` | `1`：按元素分别画力对角图（此时不画其他图） |
| `component` | `'0'` | `'force'`、`'virial'`、`'stress'`、`'dipole'`、`'polarizability'`：分量子图 |
| `charge_sign` | `1` | 电荷图正负号看起来反了时改成 `-1` |
| `charge_plot_method` | `'hist'` | `'hist'`（直方图 + KDE）或 `'kde'` |
| `plot_range` | 见脚本 | `use_range` 非 0 时各物理量的坐标范围 |
| `train_colors`, `test_colors` | 见脚本 | 分量着色 |

**注意（版本相关）**：

- loss 图按标准势函数模型的 `loss.out` 列序读取（`generation total L1 L2 E_train F_train V_train E_test F_test V_test`）。qNEP（`charge_mode`）模型的 `loss.out` 多了电荷/BEC 列，此时测试集曲线的列号与此不同，需要自己调整。
- GPUMD ≥ 5.9 的 `loss.out` 带以 `#` 开头的表头，`np.loadtxt` 会自动跳过。x 轴标签 `Generation/100` 对应默认每 100 代输出一次，改过 `output_interval` 的话标签要跟着改。
- 训练模式下，非全批次训练时 `*_train.out` 的行顺序在旧版本中不一定与 `train.xyz` 一致（脚本会提示），所以分元素力和电荷图建议用 `prediction 1` 的输出画。GPUMD ≥ 5.9 已改为按 `train.xyz` 的结构顺序写出所有 `*_train.out`。

### 7.6 get_nep_points2strucs.py：按对角图数值筛选结构

**用途**：根据 `energy_train.out`、`force_train.out`、`virial_train.out` 或 `stress_train.out` 中的数值，把对应结构挑出来写成 `deviate.xyz`，其余写成 `reserve.xyz`，用于删除异常数据或重点检查。文件类型由文件名中是否含 `energy`/`force`/`virial`/`stress` 判断；`*.out` 前一半列是 NEP 预测，后一半列是 DFT 参考值（见 [energy_*.out 文档](https://gpumd.org/nep/output_files/energy_out.html)）。

**依赖**：`ase`、`numpy`。

**用法**（脚本无参数运行时会打印同样的用法说明）：

| 模式 | 命令 | 选中条件 |
|---|---|---|
| `direct` | `python get_nep_points2strucs.py direct train.xyz energy_train.out nep -5.2 -5.0 [-4.8 -4.6 ...]` | NEP（或 `dft`）值落在任一 (起, 止) 区间内；力文件中只要有一个原子的某个分量落入区间，就选中整个结构 |
| `allin` | `python get_nep_points2strucs.py allin train.xyz force_train.out nep -10 10` | 一个结构的所有数值都在区间内 |
| `differ` | `python get_nep_points2strucs.py differ train.xyz energy_train.out -0.02 0.02` | 差值（代码中为 NEP − DFT）小于下限或大于上限 |
| `rmse` | `python get_nep_points2strucs.py rmse train.xyz force_train.out 50` | 每个结构的 RMSE 最大的 50 个 |

**输出**：`deviate.xyz`（选中）、`reserve.xyz`（其余）。

**注意**：`*.out` 的行序必须与 xyz 的结构顺序一致。旧版 GPUMD 非全批次训练时请先做一次 `prediction 1`，用预测输出的 `*.out`（GPUMD ≥ 5.9 训练输出已按 `train.xyz` 顺序）。

### 7.7 select_pick_structure.py：描述符降维与最远点采样

**用途**：把一个或多个 xyz 文件中每个结构的描述符（对原子取平均）合在一起降到二维，然后做三件事之一：

- `all`：画出所有点，按每原子能量着色（`combine` 模式），用来观察构型空间覆盖情况；
- `select`：最远点采样（FPS），按最小距离阈值或目标个数选出代表性结构，用于训练集去冗余；
- `pick`：在降维图上框选矩形区域，把区域内的结构挑出来。

**依赖**：`ase`、`numpy`、`scipy`、`matplotlib`、`scikit-learn`；`umap` 方法需要 `umap-learn`；没有描述符文件时需要 `NepTrain`（`from NepTrain.core.nep import *` 中的 `Nep3Calculator`）和一个 `nep.txt`。

**文件头变量**：

| 变量 | 默认 | 说明 |
|---|---|---|
| `nep_name` | `"nep.txt"` | 计算描述符用的 NEP 势 |
| `all_points_plot` | `'combine'` | `combine`：所有文件的点画在一起；`separate`：按文件分色 |
| `sel_des_method` | `'before'` | `before`：在原始描述符空间做 FPS；`after`：在二维降维空间做 FPS |
| `xyz_file` | `["train.xyz", "test.xyz"]` | 输入结构文件 |
| `des_file` | 全为 `"0"` | 与 `xyz_file` 一一对应的描述符文件；`"0"` 表示先找缓存 `descriptor-<名>.out`，没有就用 NepTrain 计算并缓存。也可以填 nep 预测模式 `output_descriptor 1` 输出的逐结构 `descriptor.out` |

**命令行**：

```bash
python select_pick_structure.py all pca
python select_pick_structure.py select pca 0.05 200        # 0.05 -> min_distance; 200 -> select 200 structures
python select_pick_structure.py pick umap -2 0 1 3 4 6 -1 1  # two boxes: x[-2,0] y[1,3] and x[4,6] y[-1,1]
```

可选的降维方法：`pca umap tsne kpca svd isomap se ica fa mds grp srp`。`select` 后面的数，小于 1 的当作 FPS 最小距离，大于等于 1 的当作要选的结构个数，可以一次给多个。

**输出**：

| 模式 | 文件 |
|---|---|
| `all` | `all-points-<方法>.png` |
| `select` | `select_<值>.xyz`、`abandon_<值>.xyz`、`select-<值>-<方法>.png` |
| `pick` | `pick.xyz`、`retain.xyz`、`retain-pick-<方法>.png` |

**注意**：脚本给 UMAP、MDS 等方法固定了 `random_state=42`，但 t-SNE 没有固定种子，每次运行坐标可能不同，`pick` 框选的坐标要用同一次运行的图来确定。FPS 的最小距离阈值与描述符尺度有关，`before` 和 `after` 两种模式的合适阈值差别很大。

### 7.8 create_phonon_compare.py：NEP 与 DFT 声子谱对比

**用途**：用 GPUMD 的 `compute_phonon`（有限位移法）算 NEP 声子谱，并与 phonopy 导出的 DFT 声子谱画在同一张图上，用来检验 NEP 的二阶力常数精度。

> [!IMPORTANT]
> 本脚本生成的是 **GPUMD 5.0 之前**的声子输入：`basis.in` + 显式 k 点的 `kpoints.in`（第一行为 k 点数，其后为 1/Å 单位的 k 矢量），对应旧语法 `compute_phonon <cutoff> <displacement>`。GPUMD ≥ 5.0 简化了 `compute_phonon`：语法变为 `compute_phonon <displacement>`，run.in 开头要写 `replicate`，`kpoints.in` 改为“高对称点分数坐标 + 名称”，不再需要 `basis.in`，`omega2.out` 第一列改为路径坐标。用新版 GPUMD 时请按[新文档](https://gpumd.org/gpumd/input_parameters/compute_phonon.html)准备输入，或参考 tang070205/tools 中已适配新格式的版本。另外，GPUMD 5.1 和 5.2 的 `compute_phonon` 存在倒格矢相关的 bug（5.3 已修复），这两个版本应避免使用。

**依赖**：`ase`、`numpy`、`matplotlib`。

**文件头变量**：

| 变量 | 默认 | 说明 |
|---|---|---|
| `cx, cy, cz` | `20, 20, 1` | 超胞倍数（示例为二维六方材料） |
| `npoints` | `300` | 路径上的总 k 点数，必须是路径段数的整数倍 |
| `special_points` | `G, K, M` | 高对称点（分数坐标），可以从 `vaspkit 305` 等工具得到 |
| `points_path` | `['GMKG']` | 路径，断开的路径写成 `['GM', 'KG']` |
| `dft_file` | `'phonon.out'` | `phonopy-bandplot --gnuplot > phonon.out` 的输出，单位 THz |
| `uc = read('POSCAR')` | `POSCAR` | 原胞（xyz、cif 也可以） |

**用法（运行两次）**：

```bash
python create_phonon_compare.py     # 1st: writes model.xyz, basis.in, kpoints.in, then stops (no omega2.out yet)
gpumd                               # run.in contains: potential nep.txt / compute_phonon <cutoff> <displacement>
python create_phonon_compare.py     # 2nd: reads omega2.out (+ phonon.out) -> phonon.png
```

第一次运行时，生成输入文件后会因为还没有 `omega2.out` 而报错退出，这是正常现象。频率换算为 $\nu=\operatorname{sign}(\omega^2)\sqrt{|\omega^2|}/2\pi$（THz），虚频画为负值。

**输出**：`model.xyz`、`basis.in`、`kpoints.in`、`phonon.png`（DFT 橙色、NEP 蓝色，纵轴 0～35 THz 写死在脚本里，需按体系改）。

**补充**：脚本中有一段被注释掉的 Quantum ESPRESSO 处理代码，用于把 `matdyn` 的 `*.freq.gp`（cm⁻¹）转成 THz（除以 33.35641）后作为 `phonon.out` 使用。用 VASP + phonopy 时不需要它。

---

## 8. gpumd/ 分组建模与热输运后处理

### 8.0 这些脚本读取的 run.in 字段

热输运绘图脚本会解析 `run.in`，按固定的位置取参数。关键词语法见 GPUMD 文档：[compute](https://gpumd.org/gpumd/input_parameters/compute.html)、[compute_shc](https://gpumd.org/gpumd/input_parameters/compute_shc.html)、[compute_hnemd](https://gpumd.org/gpumd/input_parameters/compute_hnemd.html)、[compute_hac](https://gpumd.org/gpumd/input_parameters/compute_hac.html)、[heat_lan 等热输运系综](https://gpumd.org/gpumd/input_parameters/ensemble_heat.html)。

| run.in 行（GPUMD 语法） | 脚本取用的字段 | 用到的脚本 |
|---|---|---|
| `time_step <dt>` | `dt`（fs） | `plot_compute_shc`、`plot_kappa_shc`、`plot_kappa_multiple`、`plot_shc_quantum_correction` |
| `ensemble nvt_... <T_1> ...` | 第 3 个词作为温度 T（整数） | `plot_kappa_shc`、`plot_shc_multiple`、`plot_shc_quantum_correction` |
| `ensemble heat_lan <T> <T_coup> <delta_T> <source> <sink>` | T、delta_T、热源/热沉组号 | `plot_compute_shc` |
| `compute <grouping> <sample> <output> temperature` | 只能有一个量（整行 5 个词）；输出间隔 = sample × output 步 | `plot_compute_shc` |
| `compute_shc <sample_interval> <Nc> <direction> <num_omega> <max_omega> [group ...]` | Nc、方向、频点数、最大角频率 | 所有 SHC 脚本 |
| `compute_hnemd <output_interval> <Fe_x> <Fe_y> <Fe_z>` | 输出间隔，以及哪个分量 > 0（判断方向） | HNEMD 脚本 |
| `compute_hac <sample> <Nc> <output_interval>` | Nc | `plot_hac_multiple` |

> [!TIP]
> 脚本用 `int()` 解析温度、`max_omega` 等数值，请在 run.in 中写成整数（例如 `compute_shc 2 250 0 1000 400`，不要写 `400.0`）。另外，SHC 和 HNEMD 脚本要用 `thermo.out` 最后 10 行的盒子尺寸求体积，因此平衡阶段要加 `dump_thermo`（脚本运行时也会打印这个提醒）。

一个 HNEMD + SHC 的 run.in 示意片段（参数仅作说明，需按体系调整）：

```text
potential     nep.txt
velocity      300
time_step     1
ensemble      nvt_nhc 300 300 100
dump_thermo   1000
run           1000000

ensemble      nvt_nhc 300 300 100
compute_hnemd 1000 0.00001 0 0
compute_shc   2 250 0 1000 400
run           10000000
```

用到的 GPUMD 输出文件格式：[kappa.out](https://gpumd.org/gpumd/output_files/kappa_out.html)（5 列：κ_μx^in、κ_μx^out、κ_μy^in、κ_μy^out、κ_μz^tot）、[shc.out](https://gpumd.org/gpumd/output_files/shc_out.html)（每次运行有 2Nc−1 行 K(t) 和 num_omega 行 J(ω)）、[hac.out](https://gpumd.org/gpumd/output_files/hac_out.html)（11 列）、[compute.out](https://gpumd.org/gpumd/output_files/compute_out.html)（各组温度，最后两列为热源/热沉恒温器累计交换的能量）、[thermo.out](https://gpumd.org/gpumd/output_files/thermo_out.html)（18 列，盒子分量 ax…cz 在最后 9 列）。in/out 分解（面内/面外模式）主要用于二维材料，三维材料把两者相加即可。

### 8.1 split_group.py：沿一个方向按周期数分组

**用途**：用原胞扩出超胞，沿 x/y/z 按“每组包含多少个原胞周期”分组，写出带 `group` 列的 `model.xyz`，供 NEMD（固定区、热源、热沉、中间区）或 SHC 分组使用。

**依赖**：`ase`、`numpy`。

**文件头变量**：

| 变量 | 默认 | 说明 |
|---|---|---|
| `uc = read('POSCAR')` | `POSCAR` | 原胞，要求为正交胞（取 `cell[0][0]` 等对角元作为周期长度） |
| `cx, cy, cz` | `15, 1, 1` | 扩胞倍数 |
| `struc.set_pbc(...)` | `[True, True, True]` | 周期性，二维材料可改成 `[True, True, False]` |
| `dic` | `'x'` | 分组方向 |
| `group_cyclical` | `[3,3,3,3,3]` | 每组的周期数，之和应等于该方向的扩胞倍数 |

**输出**：`model.xyz`；屏幕打印各组边界坐标和每组原子数，可以用它们核对 `heat_lan` 的组号。

### 8.2 split_heterojunction_group.py：异质结/晶界拼接与分组

**用途**：把多块材料（各自扩胞）沿 `dir` 方向依次拼接，块与块之间留出间隙，然后按每块各自的周期数分组，适用于异质结、晶界和多层结构的 NEMD/SHC 分组。

**文件头变量**：

| 变量 | 默认 | 说明 |
|---|---|---|
| `dir` | `'z'` | 拼接和分组的方向 |
| `merge_interface` | `0` | `1`：把前一块的最后一组和后一块的第一组合并为一个“界面组” |
| `vacuums` | `[1, 1]` | 相邻两块之间的间隙（Å），长度为块数 − 1 |
| `ucs` | 三个示例 | 每项为 `(文件, cx, cy, cz, 每组周期数列表)`，例如 `('Si.vasp', 8, 8, 120, [8,20,4,...])` |

**注意**：最终晶胞在另外两个方向沿用第一块材料的尺寸，各块在这两个方向上的长度要事先匹配好（脚本不处理晶格失配）。拼接后请用 OVITO/VESTA 检查界面距离，必要时先做弛豫。屏幕打印的各组边界可以用来核对分组；每组原子数请以 `model.xyz` 中 `group` 列的统计为准。

### 8.3 plot_compute_shc.py：NEMD + SHC 单次结果

**用途**：用于一次 NEMD（`ensemble heat_lan`）+ `compute_shc` 运行，画四联图 `compute-shc.png`：

- (a) 热源到热沉之间各组的时间平均温度（取后半段平均）；
- (b) 热源、热沉恒温器累计交换的能量随时间的变化；
- (c) 维里-速度关联函数 $K(t)$；
- (d) 谱热导 $G(\omega)$。

同时保存 `Gc.npy`，供 HNEMD 脚本计算声子平均自由程谱。

脚本开头会先按与 8.1 节相同的方式建模并**重写 `model.xyz`**（`pbc=[T, T, F]`），参数要与模拟时一致。

**文件头变量**：`dic`（输运方向）、`uc = read('POSCAR')`、`cx, cy, cz = 52, 30, 1`、`group_cycl = [9,20,1,1,1,20]`（每组周期数）、`Emax = 5`（图 b 纵轴范围，keV）。

**计算关系**（与 GPUMD 文档和 Fan 等人 2019 年的论文一致）：

$$Q=\frac12\left[\frac{E_\mathrm{in}(t_{1/2})-E_\mathrm{in}(t_\mathrm{end})}{t_\mathrm{end}-t_{1/2}}+\frac{E_\mathrm{out}(t_\mathrm{end})-E_\mathrm{out}(t_{1/2})}{t_\mathrm{end}-t_{1/2}}\right],\qquad G(\omega)=\frac{J_q(\omega)}{V\,\Delta T}$$

$Q$ 取后半段（稳态）热源、热沉能量变化率的平均（eV/ps）。$\Delta T$ 是热源组与热沉组的平均温差，$V$ 是“一个原胞周期长度 × 横截面积”，所以 `compute_shc ... group` 指定的 SHC 组应正好为 1 个周期（例如 `group_cycl` 中的 `1`）。$G(\omega)$ 的单位是 GW m⁻² K⁻¹ THz⁻¹，转换系数 $1.6\times10^4$ 已写在脚本中。

**输入**：`POSCAR`、`run.in`、`compute.out`、`shc.out`、`thermo.out`。

**注意**：时间轴按 `time_step 1`（1 fs）换算成 ps；各子图的坐标范围写死在脚本中，需按体系修改。脚本内部还按 $G=Q/(A\,\Delta T)$ 算出了整体热导（变量 `G`，单位 GW m⁻² K⁻¹，$A$ 为横截面积），需要这个数值时可以在末尾加一行 `print(G)`。

### 8.4 plot_kappa_shc.py：HNEMD + SHC 单次结果

**用途**：针对一次 HNEMD + SHC 运行：

- `hnemd.png`：κ_in、κ_out 和总 κ 的累积平均 $\kappa(t)=\frac1t\int_0^t\kappa(t')\,dt'$（z 方向只有总 κ）；
- `shc.png`：(a) $K(t)$；(b) 谱热导率 $\kappa(\omega)=\dfrac{1602.17662\,J_q(\omega)}{V\,T\,F_e}$（W m⁻¹ K⁻¹ THz⁻¹）；(c) 声子平均自由程 $\lambda(\omega)=|\kappa(\omega)/G(\omega)|$（nm，$G(\omega)$ 来自 8.3 节保存的 `Gc.npy`）；(d) 长度依赖热导率 $\kappa(L)=\int\dfrac{\kappa(\omega)}{1+\lambda(\omega)/L}\,\dfrac{d\omega}{2\pi}$。

**依赖**：`ase`、`numpy`、`scipy`、`matplotlib`。

**输入**：`POSCAR`、`run.in`、`kappa.out`、`shc.out`、`thermo.out`、`Gc.npy`。

**注意**：

- 脚本开头也会用 `uc*(10,10,10)`、全周期**重写 `model.xyz`**。如果你的 model.xyz 是用别的方式构建的，请先注释掉这几行。
- 时间轴按“每 1000 步输出一次 kappa.out”换算成 ns（与 `compute_hnemd 1000 ...` 对应）。
- HNEMD 和 NEMD 两次运行的 `compute_shc` 必须用相同的 `num_omega` 和 `max_omega`，`Gc.npy` 才能与 $\kappa(\omega)$ 逐点相除。

### 8.5 plot_kappa_multiple.py：多次独立 HNEMD 运行

**用途**：同一个 run.in 中重复多段“平衡 + HNEMD”（或多次运行，输出追加在同一个 `kappa.out` 里）时，把 `kappa.out` 按次数切开，画每次的累积平均（灰色细线）和多次平均（粗线），并在曲线末端标注最终数值。输出 `hnemd-x-multiple.png`、`hnemd-y-multiple.png` 或 `hnemd-z-multiple.png`。

**用法**：

```bash
python plot_kappa_multiple.py 5 10000000    # 5 independent runs, each "run 10000000" steps
```

| 位置参数 | 含义 |
|---|---|
| 1 | 独立运行次数 |
| 2 | 每次 HNEMD 生产阶段的步数（即 `run` 后面的数），脚本据此求每次的行数 = 步数 / `compute_hnemd` 输出间隔 |

**注意**：驱动力方向按 `compute_hnemd` 中哪个分量 > 0 来判断；各子图的纵轴范围写死在脚本中，需要按体系修改。

### 8.6 plot_shc_multiple.py：多次 HNEMD-SHC 平均

**用途**：把多次运行追加在同一个 `shc.out` 中的数据按次数平均，再像 8.4 节一样画 $K(t)$、$\kappa(\omega)$、$\lambda(\omega)$、$\kappa(L)$，输出 `shc-multiple.png`。

**用法**：`python plot_shc_multiple.py 5`（参数为运行次数）。需要 `run.in`、`shc.out`、`thermo.out`、`Gc.npy`。

**注意**：脚本用 `mask = shc_nu <= 50` 只取 50 THz 以下的频率作图，`max_omega` 设得更大时请按注释修改这个阈值。

### 8.7 plot_shc_quantum_correction.py：谱热导率量子修正

**用途**：(1) 对多次 HNEMD 的总 κ 求平均和标准误差（$\sigma/\sqrt{n}$），并打印 `k = ... ± ...`；(2) 对谱热导率做量子修正：

$$\kappa^\mathrm{q}=\int \kappa(\omega)\,\frac{x^2e^x}{(e^x-1)^2}\,\frac{d\omega}{2\pi},\qquad x=\frac{h\nu}{k_\mathrm BT}$$

打印 `k_quantum_corrected = ... ± ...`，并输出 `kappa_shc.png`：(a) 各次及平均 κ(t) ± 误差；(b) 经典与量子修正后的 κ(ν)。方法见 Wang 等人 2023 年的论文（PRB 107, 054303）。

**依赖**：`numpy ≥ 2.0`（用到 `np.trapezoid`）、`scipy`、`matplotlib`。

**用法**：

```bash
python plot_shc_quantum_correction.py 5 10000000    # runs, steps per HNEMD production run
```

**注意**：(a) 图的时间轴按“每输出点 1 ps”换算成 ns（对应 `time_step 1` 和 `compute_hnemd 1000`）；两张图的坐标范围写死，需要修改。这种修正只是把每个模式的经典热容换成量子谐振子热容，低温下的结果仍需谨慎解读。

### 8.8 plot_hac_multiple.py：EMD（Green-Kubo）多次运行

**用途**：多次独立 EMD 运行的 `hac.out` 追加在一起时，按次数切分，画四联图 `emd-multiple.png`：(a) 归一化 HAC（面内 x+y 的 in、out，双对数坐标）；(b) 每次及平均的 $\kappa^{in}$，按 $(\kappa_x+\kappa_y)/2$ 计算；(c) $\kappa^{out}$；(d) in、out 和总 κ。

**用法**：直接运行 `python plot_hac_multiple.py`。脚本从 `compute_hac` 行读取 Nc，并按 `output_interval = 10` 计算每次运行的行数（Nc/10），所以 run.in 应写成 `compute_hac <sample> <Nc> 10` 的形式，否则要修改脚本中的 `/10`。

**注意**：这里把 x、y 平均，适用于面内各向同性的二维材料（如石墨烯）；各向异性材料请分别处理。坐标范围（0～1000 ps 等）写死在脚本中。

---

## 9. analysis/ 其他分析

### 9.1 extract_vacuum_level_cube.py：平面平均静电势与真空能级

**用途**：读取静电势 cube 文件（脚本默认文件名是 CP2K `&V_HARTREE_CUBE` 的输出 `<项目名>-v_hartree-1_0.cube`），沿 z 方向做 xy 平面平均，乘以 27.2114 把 Hartree 换成 eV，取曲线最大值作为真空能级 $E_\mathrm{vac}$。功函数为 $W=E_\mathrm{vac}-E_\mathrm F$，其中 $E_\mathrm F$ 从同一次计算的输出中读取。

**依赖**：`numpy`、`matplotlib`、`ase`（`ase.io.cube.read_cube_data`）。

**用法**：把文件末尾的 `cube_filename` 改成你的文件名，然后运行 `python extract_vacuum_level_cube.py`。大 cube 文件读取可能要几十秒。

**输出**：

| 文件 | 内容 |
|---|---|
| `Electrostatic_Potential_1D.txt` | 两列：z（Å）、平面平均势（eV），可以直接导入 Origin |
| `Vacuum_Level.txt` | 真空能级数值（eV） |
| `Electrostatic_Potential.png` | 曲线图，红色虚线为 $E_\mathrm{vac}$ |

**注意**：

- 脚本写了 `matplotlib.use('TkAgg')` 并在最后调用 `plt.show()`，在没有显示环境的集群上会报 `Cannot load backend 'TkAgg' which requires the 'tk' interactive framework, as 'headless' is currently running`。在服务器上运行时，请把这一行改成 `matplotlib.use('Agg')`，并删掉 `plt.show()`。
- z 轴取 `cell[2,2]`，要求表面法向沿 z、晶胞 c 轴与 z 平行。
- 取最大值作 $E_\mathrm{vac}$ 的前提是真空区足够厚、势能平台足够平。不对称的 slab 加了偶极修正后，两侧平台高度不同，请看图确认取的是哪一侧。
- CP2K 输出 cube 和读取费米能级的方法见 [CP2K 部分](../CP2K/README.md)，通用作图见 Common 第5章。

---

## 10. GPUMD/NEP 版本说明

| GPUMD 版本 | 与本目录脚本相关的变化 |
|---|---|
| 3.7 | `virial_*.out` 由 2 列改为 12 列（6 个 NEP 分量 + 6 个参考分量），`plot_nep_results.py` 按 12 列处理 |
| 4.0 | `nep` 预测模式可输出逐结构/逐原子描述符（`output_descriptor 1/2`），即 `descriptor.out` |
| 4.6 | 引入 qNEP（带电荷），新增 `charge_*.out` |
| 4.8 | 新增 `gnep` 可执行程序（用解析梯度和 Adam 训练 NEP，输入文件是 `gnep.in`） |
| 4.9.1 | 没有位力参考的结构，目标值统一写成 −1e6 |
| 5.0 | `compute_phonon` 输入简化（见 7.8 节） |
| 5.1–5.2 | `compute_phonon` 有倒格矢相关 bug，5.3 修复 |
| 5.9 | `loss.out` 增加 `#` 表头；所有 `*_train.out` 按 `train.xyz` 顺序输出 |

截至 2026-11，GPUMD 最新版本为 v5.9.1（该版本把 `src/` 的许可证改为 LGPL-3.0-or-later，计算功能没有变化）。

---

## 11. 常见问题与报错

| 报错或现象 | 原因 | 解决 |
|---|---|---|
| `ModuleNotFoundError: No module named 'umap'` | 没装 UMAP | `pip install umap-learn` |
| `ModuleNotFoundError: No module named 'hiphive'` / `'dpdata'` | 没装对应的包 | `pip install hiphive dpdata` |
| `AttributeError: module 'numpy' has no attribute 'trapezoid'` | `plot_shc_quantum_correction.py` 需要 NumPy ≥ 2.0 | `pip install -U "numpy>=2"` |
| 读 `running_scf.log` 时提示 `abacus-out` 格式不存在 | 环境中是官方 ASE，没有 ABACUS 读取器 | 在单独的环境中安装 ase-abacus（见 3.1 节） |
| `FileNotFoundError: ... 'INCAR-single'`（或 `INCAR-scf`、`INPUT-scf`） | 模板文件名与脚本要求的不同 | 按 6.0 节的表格准备对应文件名 |
| `FileExistsError: [Errno 17] File exists: 'train_folders/...'` | 上次运行留下了目录 | 删掉或改名旧目录后重跑 |
| 生成的 POTCAR 不对或报 POTCAR 错误 | `~/.vaspkit` 中赝势路径或类型没配好 | 修改 `~/.vaspkit` 后，删掉旧 POTCAR 再运行 |
| `dos2unix: command not found` | 系统没装 dos2unix | `sudo apt install dos2unix` 或 `conda install -c conda-forge dos2unix` |
| `NameError: name 'Ns' is not defined`（plot_compute_shc） | run.in 中的 `compute` 行不是 5 个词（例如同时算多个量） | 写成 `compute 0 10 100 temperature` 的形式 |
| `ValueError: invalid literal for int() with base 10: '400.0'` | run.in 中 `max_omega` 或温度写成了小数 | 改成整数 |
| `FileNotFoundError: ... 'Gc.npy'` | HNEMD-SHC 脚本需要 NEMD 的谱导热 | 先运行 `plot_compute_shc.py` 生成 `Gc.npy`，复制到 HNEMD 目录 |
| `create_phonon_compare.py` 第一次运行报找不到 `omega2.out` | 两次运行流程中的第一次，正常现象 | 先跑 gpumd 再运行第二次 |
| 对角图中位力点集中在 −1e6 附近 | 部分结构没有位力 | `plot_nep_results.py` 会自动剔除这些行；确认训练集中有位力的结构确实写了 `virial=` |

---

## 12. 来源与致谢

| 脚本 | 代码中可见的来源信息 | 补充说明 |
|---|---|---|
| `nep/multipleFrames-abacus2nep-exyz.sh`、`nep/singleFrame-abacus2nep-exyz.sh` | 头部注明 `Modified by Benrui Tang`（tang070205@proton.me） | GPUMD 仓库 [tools/Format_Conversion/abacus2xyz](https://github.com/brucefan1983/GPUMD/tree/master/tools/Format_Conversion/abacus2xyz) 收录了几乎相同的版本；这类脚本沿用了 GPUMD tools 中 Yanzhou Wang 等人 `outcar2nep-exyz` 系列的写法 |
| `dataset/dpdata-system-perturb.py` | dpdata `System.perturb` 的摘录和修改 | dpdata 采用 LGPL-3.0 许可，再分发修改版时须遵守 |
| `dataset/aimd_OUTCAR_xyz2POSCAR.py` | 注释：部分代码来自知乎“叠加态”的文章 https://zhuanlan.zhihu.com/p/397866915 | 文件头标注作者 Gaoshan Zeng；tang070205/tools 中有同名脚本（版本不同） |
| `dataset/create-perturb.py` | 注释中引用 https://github.com/tang070205/tools ；文件头标注作者 Gaoshan Zeng | tang070205/tools 中有同名脚本（版本不同） |
| `nep/abacus_scfoutput2nepxyz.py` | 文件头标注作者 Gaoshan Zeng | GPUMD tools/abacus2xyz 与 tang070205/tools 中的同名脚本除文件头外基本相同 |
| `nep/plot_nep_results.py` | 文件头标注作者 Gaoshan Zeng | tang070205/tools 中有同名脚本（版本不同） |
| `gpumd/` 下全部 8 个脚本，以及 `nep/get_nep_points2strucs.py`、`nep/inorexyz2exyz_ratio.py`、`nep/select_pick_structure.py`、`dataset/create_strain_deform_rattle_perturb.py`、`modeling/replace_dopant_symbols.py` | 代码中无作者信息 | 与 Benrui Tang 公开的 [tang070205/tools](https://github.com/tang070205/tools) 中的同名脚本内容一致 |
| `nep/create_phonon_compare.py` | 代码中无作者信息 | tang070205/tools 中有同名脚本，上游版本已改为适配 GPUMD ≥ 5.0 的输入格式 |
| `dataset/dumpxyz2POSCAR.py`、`dataset/messxyz2POSCAR.py` | 代码中无作者信息 | tang070205/tools 的 README 中介绍过同名脚本（上游当前提供的是 mdxyz2POSCAR.py） |
| `modeling/fabricate_spherical_carbon.py`、`analysis/extract_vacuum_level_cube.py` | 代码中无作者信息 | 本仓库贡献者 |

感谢以上作者公开代码。在自己的项目中复用或改写这些脚本时，请保留原作者信息；需要再分发时，请先确认原作者的授权。

---

## 参考资料

**官方文档**

- GPUMD 文档首页：https://gpumd.org
- NEP（`nep` 可执行程序）文档：https://gpumd.org/nep/index.html
- `gpumd` 可执行程序文档（输入关键词、输出文件）：https://gpumd.org/gpumd/index.html
- NEP 训练集格式 train.xyz/test.xyz：https://gpumd.org/nep/input_files/train_test_xyz.html
- NEP 输出文件 loss.out、energy/force/virial_*.out、descriptor.out：https://gpumd.org/nep/output_files/index.html
- GPUMD 新版 compute_phonon 与 kpoints.in：https://gpumd.org/gpumd/input_parameters/compute_phonon.html 、https://gpumd.org/gpumd/input_files/kpoints_in.html
- ABACUS 文档：https://abacus.deepmodeling.com/
- VASPKIT 文档（命令行模式 `vaspkit -task 102 -kpr 0.04`）：https://vaspkit.com/tutorials.html
- ATOMKIT：https://vaspkit.com/atomkit.html
- hiphive 文档：https://hiphive.materialsmodeling.org/

**相关 GitHub/GitLab 仓库**

- [brucefan1983/GPUMD](https://github.com/brucefan1983/GPUMD)：GPUMD/NEP 源码，`tools/` 下有各种 DFT 输出转 NEP extxyz 的官方收录脚本（vasp2xyz、abacus2xyz、cp2k2xyz 等）。
- [brucefan1983/GPUMD-Tutorials](https://github.com/brucefan1983/GPUMD-Tutorials)：GPUMD 官方教程合集，有 NEMD/HNEMD/SHC、声子等案例，可以和本目录脚本对照使用。
- [tang070205/tools](https://github.com/tang070205/tools)：Benrui Tang 的 GPUMD/NEP 脚本集，本目录多数脚本与它同源，可以在那里找到更新的版本。
- [deepmodeling/dpdata](https://github.com/deepmodeling/dpdata)：多种 DFT/MD 格式互转，以及结构微扰（LGPL-3.0）。
- [ase-abacus](https://gitlab.com/1041176461/ase-abacus)：带 ABACUS 读写接口的 ASE 分支。
- [hiphive](https://gitlab.com/materials-modeling/hiphive)：力常数拟合工具，这里用到它的 MC rattle 结构生成。
- [aboys-cb/NepTrain](https://github.com/aboys-cb/NepTrain)：NEP 自动训练流程（微扰、主动学习、单点、训练），`select_pick_structure.py` 用它算描述符。
- [aboys-cb/NepTrainKit](https://github.com/aboys-cb/NepTrainKit)：NEP 训练集可视化与筛选的图形界面工具，可以代替手工筛点。
- [bigd4/PyNEP](https://github.com/bigd4/PyNEP)：NEP 的 Python 接口，提供数据集处理和采样工具。
- 更多相关仓库见 [Resources/github-repos.md](../Resources/github-repos.md)。

**关键文献**

- Z. Fan *et al.*, Neuroevolution machine learning potentials: Combining high accuracy and low cost in atomistic simulations and application to heat transport, *Phys. Rev. B* **104**, 104309 (2021). DOI: [10.1103/PhysRevB.104.104309](https://doi.org/10.1103/PhysRevB.104.104309)
- Z. Fan *et al.*, GPUMD: A package for constructing accurate machine-learned potentials and performing highly efficient atomistic simulations, *J. Chem. Phys.* **157**, 114801 (2022). DOI: [10.1063/5.0106617](https://doi.org/10.1063/5.0106617)
- Z. Fan *et al.*, Homogeneous nonequilibrium molecular dynamics method for heat transport and spectral decomposition with many-body potentials, *Phys. Rev. B* **99**, 064308 (2019)（HNEMD、SHC、$\lambda(\omega)$ 与 $\kappa(L)$）。DOI: [10.1103/PhysRevB.99.064308](https://doi.org/10.1103/PhysRevB.99.064308)
- Z. Fan *et al.*, Thermal conductivity decomposition in two-dimensional materials: Application to graphene, *Phys. Rev. B* **95**, 144309 (2017)（in/out 分解）。DOI: [10.1103/PhysRevB.95.144309](https://doi.org/10.1103/PhysRevB.95.144309)
- Y. Wang *et al.*, Quantum-corrected thickness-dependent thermal conductivity in amorphous silicon predicted by machine learning molecular dynamics simulations, *Phys. Rev. B* **107**, 054303 (2023)（量子修正）。DOI: [10.1103/PhysRevB.107.054303](https://doi.org/10.1103/PhysRevB.107.054303)
- Efficient GPU-accelerated training of a neuroevolution potential with analytical gradients, *Comput. Phys. Commun.* **320** (2026)（GNEP）。DOI: [10.1016/j.cpc.2025.109994](https://doi.org/10.1016/j.cpc.2025.109994)
- F. Eriksson, E. Fransson, P. Erhart, The Hiphive Package for the Extraction of High-Order Force Constants by Machine Learning, *Adv. Theory Simul.* **2**, 1800184 (2019). DOI: [10.1002/adts.201800184](https://doi.org/10.1002/adts.201800184)

---

[返回仓库首页](../README.md) | Common 第7章：学习路线与机器学习势 | [VASP 第15章：MD 与 MLFF](../VASP/tutorials/15-md-mlff.md)
