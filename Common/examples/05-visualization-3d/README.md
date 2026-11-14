# 05-visualization-3d：结构、体数据、轨迹与振动的可视化（VESTA、OVITO、nglview、Multiwfn、phonon website）

## 目的

能带和 DOS 之外，DFT 结果还有大量“三维”信息需要看：电荷密度和差分电荷、ELF、分子轨道与 ESP 表面、AIMD 轨迹、声子和分子振动。本目录给出几种常用工具的可直接使用的脚本与操作清单：

- `vesta_isosurface_howto.md`：VESTA 打开 CHGCAR/PARCHG/ELFCAR/LOCPOT/cube，设置等值面、截面，在 VESTA 中做 AB − A − B 差分电荷，导出位图与矢量图；
- `ovito_render_xdatcar.py`：用 OVITO 的 pip 版 Python 模块在无图形界面的环境下把 XDATCAR 渲染成图片和动画，并计算时间平均的 RDF；
- `nglview_demo.ipynb`：在 Jupyter 中用 nglview 和 weas-widget 交互查看结构与轨迹，编辑后写回 POSCAR；
- `run_multiwfn_batch.sh`：把同一个 Multiwfn 命令流批量应用到多个 fchk/molden/wfx/cube 文件，每个文件单独一个输出目录；
- `band_phononwebsite.conf`：phonopy 的 band.conf，输出带本征矢的 band.yaml，可上传到 phonon website 播放声子振动动画。

## 文件说明

| 文件 | 依赖 | 说明 |
|---|---|---|
| `vesta_isosurface_howto.md` | VESTA 3.5.8 或 3.90.x | 操作清单（GUI） |
| `ovito_render_xdatcar.py` | `pip install ovito`（OVITO 3.16.1+ 要求 Python ≥ 3.11，官方建议装在普通 venv 中而不是 Anaconda 环境） | 静态图 `traj_frame0.png`、动画 `traj.mp4`、可选 `traj_rdf.dat` |
| `nglview_demo.ipynb` | nglview、weas-widget、ase、JupyterLab | 没有 DFT 输出时用 ASE 自建的 GaAs 与 CO/Cu(111) 演示 |
| `run_multiwfn_batch.sh` | Multiwfn（集群用 noGUI 版）；读 .chk 时需要 Gaussian 的 formchk | 命令流可直接用 [Gaussian 第16章示例](../../../Gaussian/examples/16-postprocess-multiwfn/README.md) 中经过测试的文件 |
| `band_phononwebsite.conf` | phonopy（4.x 用 `--config`） | `EIGENVECTORS = .TRUE.` + `BAND_CONNECTION = .TRUE.` |

## 运行步骤

```bash
# 1) OVITO：单独建一个 venv（不要装进 conda 的 base 或 dft 环境）
python3.12 -m venv ~/venv-ovito && source ~/venv-ovito/bin/activate
pip install ovito
python ovito_render_xdatcar.py /path/to/aimd/XDATCAR --every-nth 5 --movie water.gif --rdf-cutoff 6.0

# 2) Jupyter 交互查看
conda activate dft
conda install -c conda-forge nglview && pip install weas-widget
jupyter lab nglview_demo.ipynb

# 3) Multiwfn 批处理：苯酚等多个分子的 ESP 表面分析
export NTHREADS=8
bash run_multiwfn_batch.sh \
    ../../../Gaussian/examples/16-postprocess-multiwfn/multiwfn_esp_surface_cmds.txt  mols/*.fchk
grep -H "Global surface minimum" *_mwfn/*.out

# 4) 声子动画：在已有 FORCE_SETS 的 phonopy 目录中
phonopy --config band_phononwebsite.conf      # phonopy 4.x；2.x/3.x 写 phonopy band_phononwebsite.conf
# 浏览器打开 https://henriquemiranda.github.io/phononwebsite/ ，Choose files 选择 band.yaml
```

## 预期结果 / 如何检查

- OVITO：屏幕打印帧数、原子数和晶胞长度，当前目录出现 `traj_frame0.png` 与动画文件；动画中原子应在晶胞内连续运动，若原子“跳”到另一侧，是周期性边界的正常表现（可在 OVITO GUI 中加 Unwrap trajectories 修饰器）。RDF 第一峰位置应与该体系的典型键长一致（例如液态水 O–O 约 2.8 Å，来自实验文献）。
- nglview：单元格出现可旋转的三维结构和晶胞框；轨迹单元格下方有帧滑块。空白输出通常是前端插件未加载（见教程“常见问题”）。
- Multiwfn：每个输入文件得到一个 `<name>_mwfn/` 目录，内含 `.out` 与 Multiwfn 写出的文件；表格中状态为 `OK`。
- phonon website：页面左侧显示声子色散，点击某个模式即可看到对应的原子振动动画。

## 注意事项

- VESTA、OVITO 渲染的图片可直接用于论文，但图注中应写明等值面数值、所用软件与版本，并按各软件要求引用（VESTA：K. Momma, F. Izumi, *J. Appl. Crystallogr.* 44, 1272 (2011)；OVITO：A. Stukowski, *Modelling Simul. Mater. Sci. Eng.* 18, 015012 (2010)）。
- OVITO 基础功能免费，部分高级功能属于 OVITO Pro；pip 版 Python 模块的许可条款以官网为准。
- Multiwfn 读取 Gaussian 结果需要 `.fchk`（`formchk x.chk x.fchk`），不能直接读 `.chk`；大体系崩溃时先确认 `ulimit -s unlimited` 和 `OMP_STACKSIZE` 已生效。
- 振动动画：Gaussian 输出可在 GaussView、Avogadro 2 或 Jmol 中播放；周期体系的声子用 phonopy 计算后交给 phonon website 播放。

返回教程：[第5章 后处理与可视化](../../tutorials/05-postprocess-visualization.md)
