# 05-vaspkit-batch：VASPKIT 1.5.x 非交互批处理与常用后处理输入

## 目的

把 VASPKIT 的交互菜单变成可以写进脚本、对几十个目录批量执行的命令：

- `vaspkit_batch.sh`：对一组计算目录执行同一个 VASPKIT 功能（102、303、251、211、213、252、111/113、911、711、426、501、502），检查输入文件、保存屏幕输出、从 VASPKIT 写出的数据文件中提取结果并汇总成一张表；
- `inputs/cmd_*.in`：每个功能对应的“应答文件”，可直接 `vaspkit < inputs/cmd_211_band.in` 使用，也方便核对提示顺序；
- 三个配套 INCAR：线性光学（`INCAR_optics_LOPTICS`，供 711）、功函数（`INCAR_workfunction`，供 426）、频率与热校正（`INCAR_freq_thermo`，供 501/502）。

`~/.vaspkit` 的模板见 [01-env-setup/vaspkit_rc.example](../01-env-setup/vaspkit_rc.example)。

## 文件说明

| 文件 | 内容 |
|---|---|
| `vaspkit_batch.sh` | 批处理主脚本，用法见文件头注释 |
| `inputs/cmd_102_kpoints_gamma_0.04.in` | 102 → 2（Γ 中心）→ 0.04（分辨率，单位 2π/Å） |
| `inputs/cmd_303_kpath_bulk.in` | 303：体相 k 路径，写出 KPATH.in、PRIMCELL.vasp、HIGH_SYMMETRY_POINTS |
| `inputs/cmd_251_hse_kpoints.in` | 251 → 2 → 0.04（自洽网格）→ 0.06（路径分辨率）：杂化泛函零权重 KPOINTS，需先有 KPATH.in |
| `inputs/cmd_211_band.in`、`cmd_213_band_elements.in`、`cmd_252_hse_band.in` | 普通能带、分元素投影能带（需 PROCAR）、杂化泛函能带后处理（KPOINTS 须由 251 生成） |
| `inputs/cmd_111_tdos.in`、`cmd_113_pdos_elements.in` | 总 DOS、分元素 PDOS |
| `inputs/cmd_911_gap.in` | 带隙 |
| `inputs/cmd_711_optics_eV.in` | 711 → 1（能量单位 eV）：体相线性光学谱 |
| `inputs/cmd_426_planar_avg_z.in` | 426 → 3（z 方向）：LOCPOT 平面平均与真空能级 |
| `inputs/cmd_502_gas_298K_1atm_singlet.in` | 502 → 298.15 K → 1 atm → 自旋多重度 1：气体分子热校正 |
| `inputs/cmd_501_adsorbate_298K.in` | 501 → 298.15 K：吸附物热校正（只含振动） |
| `inputs/cmd_314_chgdiff.in` | 314 → 三个 CHGCAR 路径（AB − A − B）：差分电荷，输出 CHGDIFF.vasp |
| `INCAR_optics_LOPTICS` | LOPTICS 第二步（读 WAVECAR、`ALGO = Exact`、大 NBANDS、`CSHIFT`） |
| `INCAR_workfunction` | `LVHAR = .TRUE.` + 偶极校正的 slab 单点 |
| `INCAR_freq_thermo` | `IBRION = 5`、`NFREE = 2`、`POTIM = 0.015` 有限差分频率 |

## 运行步骤

```bash
# 0) 第一次使用：交互运行一遍，核对提示顺序与 inputs/ 中的应答是否一致
cd some_band_calc && vaspkit           # 输入 211，观察提示

# 1) 单个目录、单个功能
vaspkit < ../inputs/cmd_211_band.in > vaspkit_211.log

# 2) 批量：所有应变目录的带隙
bash vaspkit_batch.sh band strain_*/band

# 3) 批量生成 KPOINTS，分辨率改成 0.03
KPR=0.03 bash vaspkit_batch.sh kpoints relax_*/

# 4) 一批气体分子在 298.15 K、1 atm 的热校正（O2 是三重态，单独设置 MULT=3）
bash vaspkit_batch.sh thermo-gas H2O_freq CO2_freq H2_freq
MULT=3 bash vaspkit_batch.sh thermo-gas O2_freq

# 5) 差分电荷（AB、A、B 三个子目录中的 CHGCAR 网格必须一致）
vaspkit < inputs/cmd_314_chgdiff.in
```

## 预期结果 / 如何检查

脚本对每个目录输出一行 `folder | status | result`。`OK` 表示 VASPKIT 正常结束且日志中没有 “error” 字样；`SKIP` 表示缺少必要的输入文件（例如没有 LOCPOT 时跳过 426）；`CHECK` 表示日志中出现 error，需要打开 `vaspkit_<task>.log` 查看。下面是格式示意：

```text
folder                                   status   result
strain_+1/band                           OK       Band Gap (eV): ...
strain_+2/band                           SKIP     missing: EIGENVAL
```

（示意输出，数值仅供参考；`result` 列直接取自 `BAND_GAP` 等文件中的相应行。）

## 注意事项

- 功能编号与提示顺序以 VASPKIT 1.5.1 为准（[功能列表](https://vaspkit.com/features.html)）。老版中文手册（0.73）中的编号和提示可能已经变化，只能参考。
- VASPKIT 在当前目录写固定文件名（`BAND.dat`、`KLABELS`、`PDOS_*.dat` 等），**每个计算一个目录**，不要在同一目录连续处理不同计算。
- 211/213 从 DOSCAR 读 E-fermi。能带计算（ICHARG = 11）给出的 E-fermi 不可靠：能带计算结束后先把它自己写出的 DOSCAR 改名备份，再把 SCF/DOS 计算的 DOSCAR 复制过来，或在 `~/.vaspkit` 中关掉 `SET_FERMI_ENERGY_ZERO` 后自己平移。
- 252–257（杂化泛函能带后处理）只支持由 251 生成的 KPOINTS：251 把生成参数写在 KPOINTS 第一行（`Parameters to Generate KPOINTS (Don't Edit This Line)`），手写或 sumo 生成的零权重 KPOINTS 不能用 252 处理。
- `~/.vaspkit` 中 `AUTO_PLOT = .TRUE.` 时，部分功能结束前会多问一个作图选项；应答文件里缺这一行可能导致报错或卡住，批处理前先交互核对，或暂时关掉 `AUTO_PLOT`。
- 426 的平面平均数据文件名随版本不同（1.5 教程中为 `POTPAVG.dat`，较早版本为 `PLANAR_AVERAGE.dat`），真空能级见屏幕输出的 `Vacuum-Level (eV)` 一行。
- 501 只算振动贡献，并把低于 50 cm⁻¹ 的频率按 50 cm⁻¹ 处理；502 包括平动、转动、振动和电子贡献，需要正确的自旋多重度。两者都从当前目录的 OUTCAR 读取频率，频率计算要先充分弛豫。
- 屏幕输出的格式会随版本变化，按行号截取（如 `sed -n '40p'`）的写法升级后容易失效，脚本只读数据文件或用关键字匹配。

返回教程：[第5章 后处理与可视化](../../tutorials/05-postprocess-visualization.md)
