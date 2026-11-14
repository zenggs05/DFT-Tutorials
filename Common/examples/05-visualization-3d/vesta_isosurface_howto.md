# VESTA 操作清单：等值面、截面、差分电荷与出图

本清单按 VESTA 3.5.8（稳定版，2021）的菜单名称编写；3.90.x 开发版（2026 年提供 Apple Silicon、Windows on Arm、Linux GTK3 原生构建）菜单基本相同，个别位置可能调整。下载：<https://jp-minerals.org/vesta/en/download.html>。对应教程：[第5章 §10](../../tutorials/05-postprocess-visualization.md)。

## 1. 能打开哪些文件

| 文件 | 内容 | 产生方式 |
|---|---|---|
| `POSCAR` / `CONTCAR` / `*.vasp` / `*.cif` | 结构 | VASP、数据库 |
| `CHGCAR`、`CHG` | 电荷密度（文件中存的是 ρ×V_cell） | `LCHARG = .TRUE.` |
| `PARCHG` | 部分电荷（选定能带/能量窗口） | `LPARD = .TRUE.` |
| `ELFCAR` | 电子局域函数 | `LELF = .TRUE.` |
| `LOCPOT` | 局域势 / 静电势 | `LVTOT` / `LVHAR` |
| `CHGDIFF.vasp` | 差分电荷 | VASPKIT 314 |
| `*.cube` / `*.cub` | Gaussian cube 体数据 | cubegen、Multiwfn、CP2K、OpenMX |
| `*.xsf` | XCrySDen 格式 | VASPKIT 328/428、OpenMX |

> [!TIP]
> 自旋极化的 CHGCAR 包含两组数据（总密度 ρ↑+ρ↓ 和磁化密度 ρ↑−ρ↓），VESTA 打开时默认显示第一组。只看自旋密度时，先用 VASPKIT 312 或 pymatgen 提取，再打开。

## 2. 等值面（isosurface）

1. File → Open 打开 CHGCAR（大文件打开需要几秒到几分钟）。
2. Properties → Isosurfaces（或左侧 Objects 面板中的 Properties）：
   - `Isosurface level` 设定数值；电荷密度从约 0.01–0.05 开始试，差分电荷常用 0.001–0.005 量级；
   - 差分电荷、波函数等有正负的量：添加第二个等值面，Level 取相反数，选择 “Negative”/改颜色，常用黄色（正、电子积累）与青色（负、电子耗散）；
   - `Opacity` 调到 60–80% 可同时看到原子。
3. 记录所用的等值面数值，写进图注（审稿人经常问）。

> [!NOTE]
> VESTA 显示的是除以晶胞体积后的密度，等值面数值与 CHGCAR 中直接存储的数（ρ×V）不同；不同软件给出的同一体系的等值面数值不能直接混用，具体单位以 VESTA 手册为准。

## 3. 截面（2D slice）

- Edit → Lattice Planes：输入 (hkl) 和到原点的距离，在三维视图中显示切面；
- Utilities → 2D Data Display：打开二维数据窗口，Slice 设置 (hkl) 与位置，可调色标范围、等高线间隔，并导出图片或文本数据。适合画成键平面上的电荷密度或 ELF。

## 4. 在 VESTA 里直接做差分电荷

前提：AB、A、B 三个 CHGCAR **晶胞、原子坐标、FFT 网格（NGXF/NGYF/NGZF）完全一致**。A、B 片段使用 AB 弛豫后的坐标、不再弛豫，三次计算显式写相同的 NGXF/NGYF/NGZF 与 ENCUT。

1. 打开 AB 的 CHGCAR；
2. Edit → Edit Data → Volumetric Data；
3. 点 Import…，选择 A 的 CHGCAR，在弹出的对话框中选择减去（Subtract from current data），确认；
4. 再 Import B 的 CHGCAR，同样选择减去；
5. 回到主窗口，按第 2 节设置正、负两个等值面；
6. File → Save As 可另存为 `.vesta` 工程文件，下次直接打开。

等价的命令行方法：VASPKIT 314（输出 `CHGDIFF.vasp`）、VTST 的 `chgdiff.pl`、本仓库的 [chgdiff.py](../../../VASP/scripts/chgdiff.py)。网格不一致时三种方法都会报错或得到条纹状的假信号。

## 5. 出图

下表中的菜单位置随版本略有差异，找不到时在各菜单中查找同名对话框。

| 操作 | 菜单 |
|---|---|
| 正交投影、取消透视畸变 | View → Overall Appearance（或工具栏的投影按钮），选择 Parallel |
| 白色背景 | Edit → Preferences / Overall Appearance → Background |
| 隐藏晶胞边框/坐标轴 | 左侧 Boundary / Properties → Unit cell |
| 位图（PNG/TIFF） | File → Export Raster Image…，可设放大倍数（scale），出版建议 ≥ 300 dpi |
| 矢量图 | File → Export Vector Image…（等值面较复杂时文件会很大） |

## 6. 常见问题

| 现象 | 原因 | 处理 |
|---|---|---|
| 打开 CHGCAR 非常慢或内存不够 | 网格很大（NGXF 等数百） | 只在本地电脑打开；必要时用 pymatgen/VASPKIT 先做截面或降采样 |
| 差分电荷满屏条纹 | 三个 CHGCAR 网格或坐标不一致 | 显式统一 NGXF/NGYF/NGZF 与坐标后重算 |
| 等值面只出现在晶胞一角 | 体数据跨越周期边界 | Boundary 中把显示范围扩到 −0.1–1.1，或平移原点 |
| 元素显示错误或无法识别 | POSCAR 缺少第 6 行元素符号（VASP 4 旧格式），或文件带 CRLF 换行 | 补上元素行；`dos2unix POSCAR` |

返回：[05-visualization-3d](README.md)
