# VASPKIT 1.5.x 常用功能与命令行速查

适用于 VASPKIT 1.5.0/1.5.1（功能编号以 <https://vaspkit.com/features.html> 为准；旧版中文手册中的编号可能已经变动）。对应教程：[第18章 后处理生态、自动化与报错排查](../../tutorials/18-postprocess-troubleshooting.md)（§3）。

## 1. 安装后的配置

```bash
cd vaspkit.1.5.1
cp -f how_to_set_environment_variable ~/.vaspkit     # 官方安装说明的做法（也可运行包内的 setup.sh）
vi ~/.vaspkit
# 关键条目（格式为 "关键词 值 # 注释"）：
#   PBE_PATH           ~/POTCAR/PBE          # PBE 赝势库目录，103/104 生成 POTCAR 时使用；通常只需改这一项
#   POTCAR_TYPE        PBE
#   RECOMMENDED_POTCAR .TRUE.                # 使用推荐的赝势后缀
#   PYTHON_BIN         ~/miniconda3/bin/python3
#   PLOT_MATPLOTLIB    .TRUE.                # 自动出图（需要 numpy、matplotlib）
```

## 2. 运行方式

| 方式 | 示例 | 说明 |
|---|---|---|
| 交互菜单 | `vaspkit` | 逐级输入编号 |
| 命令行参数 | `vaspkit -task 102 -kpr 0.04` | 官方说明：只有部分功能支持参数模式；`vaspkit -help` 查看 |
| 管道输入 | `echo -e "102\n2\n0.04" \| vaspkit` | 适用于所有功能，按菜单顺序依次输入 |
| 输入文件 | `vaspkit < cmd.in` | `cmd.in` 每行一个输入 |

批处理示例（每个子目录生成 KPOINTS、提取带隙）：

```bash
for d in */; do
  (cd "$d" && vaspkit -task 102 -file POSCAR -kpr 0.04 > /dev/null)
done
for d in */; do
  (cd "$d" && echo -e "911\n" | vaspkit > gap.log)
done
```

## 3. 功能编号速查

| 一级菜单 | 常用子功能 | 用途 | 本仓库相关章节 |
|---|---|---|---|
| 01 输入文件 | 101 INCAR 模板、102 KPOINTS（按分辨率）、103/104 POTCAR、105 cif→POSCAR | 准备输入 | [第2章](../../tutorials/02-input-output-files.md) |
| 02 力学性质 | 200/201 应力-应变/能量-应变弹性常数、203 从 OUTCAR 读弹性常数、205 EOS 拟合 | 弹性 | [第13章](../../tutorials/13-elastic.md) |
| 03 k 路径 | 301/302/303 一维/二维/体相能带路径、304 Wannier90、305 phonopy | 能带 | [第5章](../../tutorials/05-band-dos.md) |
| 04 结构编辑 | 400 重定义晶格、401 超胞、402 固定原子、412 分数/笛卡尔坐标转换 | 建模 | |
| 05 催化 | 501/502 吸附物/气体热校正、503 能带中心、505 NEB 线性插值 | 吸附与自由能 | [第10章](../../tutorials/10-adsorption-catalysis.md) |
| 06 对称性 | 601 查对称性、602 找原胞、603 标准惯用胞 | 结构检查 | |
| 08 高级建模 | 801 真空层、803 切表面、804 异质结、821/822 缺陷 | 表面与缺陷 | [第9章](../../tutorials/09-surface-workfunction.md) |
| 11 态密度 | 111 总 DOS、113 分元素 PDOS、115 指定原子和轨道 | DOS | [第5章](../../tutorials/05-band-dos.md) |
| 21 能带 | 211 能带、213 分元素投影能带 | 能带 | [第5章](../../tutorials/05-band-dos.md) |
| 25 杂化泛函能带 | 251 生成零权重 KPOINTS、252 后处理 | HSE 能带 | [第6章](../../tutorials/06-hybrid-metagga.md) |
| 26 费米面 | 261 生成 KPOINTS、262 XCrySDen、263/264 FermiSurfer | 费米面 | [第17章](../../tutorials/17-advanced-electronic.md) |
| 28 能带反折叠 | 281 生成 KPOINTS、282 EBS | 反折叠 | [第17章](../../tutorials/17-advanced-electronic.md) |
| 31 电荷 | 311 电荷密度、312 自旋密度、314 差分电荷、316 平面平均、325 STM | 电荷分析 | [第16章](../../tutorials/16-charge-bonding.md) |
| 42 势 | 426 平面平均势（功函数）、427 宏观平均、428/429 LOCPOT/ELFCAR 格式转换 | 功函数 | [第9章](../../tutorials/09-surface-workfunction.md) |
| 51 波函数 | 511/512 实空间波函数 | 轨道可视化 | |
| 62/65 磁性与自旋纹理 | 621 磁各向异性、651–653 自旋纹理（SOC） | 磁性 | [第7章](../../tutorials/07-magnetism-u-soc.md) |
| 71 光学 | 711 体相线性光学、716/717 联合态密度 | 光学 | [第14章](../../tutorials/14-optics-dielectric.md) |
| 72 分子动力学 | 722 MSD（FFT）、723 扩散系数、726 RDF、728 振动 DOS | MD | [第15章](../../tutorials/15-md-mlff.md) |
| 91 半导体 | 911 带隙、912 有效质量、914 角分辨有效质量 | 半导体 | [第17章](../../tutorials/17-advanced-electronic.md) |
| 92 二维材料 | 921 原子层居中、922 调整真空层、927 相对真空能级的带边 | 二维 | [第8章](../../tutorials/08-vdw-2d.md) |
| 95 声子 | 1.5.0 起可从 OUTCAR 读取声子色散等（952–964） | 声子 | [第12章](../../tutorials/12-phonon.md) |

## 4. 注意

- 很多功能会在当前目录写出固定文件名（如 `KPATH.in`、`BAND.dat`、`PDOS_*.dat`），批处理时每个计算用单独目录；
- 读取 vasprun.xml 或 PROCAR 的功能对文件完整性敏感，作业中断的计算请先确认正常结束；
- 引用：V. Wang et al., *Comput. Phys. Commun.* 267, 108033 (2021)，DOI: [10.1016/j.cpc.2021.108033](https://doi.org/10.1016/j.cpc.2021.108033)。
