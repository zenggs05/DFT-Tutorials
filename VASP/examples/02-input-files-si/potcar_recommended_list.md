# VASP 推荐 PAW 赝势（POTCAR）后缀速查表（PBE，potpaw.64）

本表整理自 VASP 官方 Wiki 的 “Choosing pseudopotentials / Recommended pseudopotentials” 页面（**推荐列表自 2026 年 10 月起生效**），只列出赝势名称、价电子数（ZVAL）和 ENMAX，**不包含任何 POTCAR 文件内容**。POTCAR 本身受 VASP 许可证保护，请从你自己的 VASP Portal 下载，并用 [make_potcar.sh](make_potcar.sh) 按 POSCAR 拼接。配套教程：[第2章 四大输入文件与输出文件全解析](../../tutorials/02-input-output-files.md)。

## 后缀含义

| 后缀 | 含义 | 何时使用 |
|---|---|---|
| 无后缀（如 `Si`、`O`） | 标准价电子划分 | 多数主族元素的默认选择 |
| `_pv` | 半芯 p 态（如 3p）作为价电子 | 过渡金属、碱土金属；磁性、高压、短键体系 |
| `_sv` | 半芯 s 和 p 态都作为价电子 | 碱金属、碱土金属、早期过渡金属 |
| `_d` | 半芯 d 态（如 Ga 3d）作为价电子 | Ga、Ge、In、Sn、Tl、Pb、Bi 等 p 区元素 |
| `_s` | 软赝势（ENMAX 低） | 粗略结构筛选；**不要**用于杂化泛函、HF 或切换泛函的计算 |
| `_h` | 硬赝势（ENMAX 高） | 极短键、需要高精度 HF/杂化泛函的场合 |
| `_GW` | 针对费米能级以上高能散射性质优化 | GW、BSE、光学性质等需要大量空带的计算；整个工作流都应使用 |
| `_2`、`_3` | 镧系元素把 f 电子冻结在芯内，固定二价/三价 | 不关心 4f 物理的镧系化合物 |

## 2026 年 10 月推荐列表的变化（标准 DFT，表中标 `(*)`）

`Be → Be_sv`、`Mg → Mg_pv`、`Fe → Fe_pv`、`Co → Co_pv`、`Ni → Ni_pv`、`Cu → Cu_pv`、`Eu_2 → Eu_3`、`Yb_2 → Yb_3`、`Re → Re_pv`、`Os → Os_pv`。官方说明：GW 推荐与参考计算推荐未变；LDA 推荐做了相同调整（Eu、Yb 无 LDA 势）。

> [!IMPORTANT]
> **复现旧工作时，请使用当时的推荐和完全相同的 POTCAR**（官方页面提供 SHA256 校验值）。同一项研究中所有需要比较能量的计算必须使用同一组 POTCAR。论文中应写明每个元素的 `TITEL` 行（`grep TITEL POTCAR`）。

## 推荐列表

“1.3×ENMAX”一栏是涉及晶格应力（ISIF=3 变胞优化、状态方程、弹性常数）时 ENCUT 的常用下限；多元素体系取各元素中的最大值。GW 一栏为“需要大量空带的计算”的推荐，`-` 表示官方未列出。

| 元素 | 标准 DFT 推荐 | ZVAL | ENMAX (eV) | 1.3×ENMAX (eV) | GW/光学推荐 | GW ENMAX (eV) |
|---|---|---|---|---|---|---|
| H | H | 1 | 250.0 | 325 | H_GW | 300.0 |
| He | He | 2 | 478.9 | 623 | He_GW | 405.8 |
| Li | Li_sv | 3 | 499.0 | 649 | Li_sv_GW | 433.7 |
| Be | Be_sv (\*) | 4 | 308.8 | 401 | Be_sv_GW | 537.5 |
| B | B | 3 | 318.6 | 414 | B_GW | 318.6 |
| C | C | 4 | 400.0 | 520 | C_GW | 414.0 |
| N | N | 5 | 400.0 | 520 | N_GW | 420.9 |
| O | O | 6 | 400.0 | 520 | O_GW | 414.6 |
| F | F | 7 | 400.0 | 520 | F_GW | 487.7 |
| Ne | Ne | 8 | 343.6 | 447 | Ne_GW | 432.3 |
| Na | Na_pv | 7 | 259.6 | 337 | Na_sv_GW | 372.9 |
| Mg | Mg_pv (\*) | 8 | 403.9 | 525 | Mg_sv_GW | 429.9 |
| Al | Al | 3 | 240.3 | 312 | Al_GW | 240.3 |
| Si | Si | 4 | 245.3 | 319 | Si_GW | 245.3 |
| P | P | 5 | 255.0 | 332 | P_GW | 255.0 |
| S | S | 6 | 258.7 | 336 | S_GW | 258.7 |
| Cl | Cl | 7 | 262.5 | 341 | Cl_GW | 262.5 |
| Ar | Ar | 8 | 266.4 | 346 | Ar_GW | 290.6 |
| K | K_sv | 9 | 259.3 | 337 | K_sv_GW | 249.0 |
| Ca | Ca_sv | 10 | 266.6 | 347 | Ca_sv_GW | 281.4 |
| Sc | Sc_sv | 11 | 222.7 | 289 | Sc_sv_GW | 379.0 |
| Ti | Ti_sv | 12 | 274.6 | 357 | Ti_sv_GW | 383.8 |
| V | V_sv | 13 | 263.7 | 343 | V_sv_GW | 382.3 |
| Cr | Cr_pv | 12 | 265.7 | 345 | Cr_sv_GW | 384.9 |
| Mn | Mn_pv | 13 | 269.9 | 351 | Mn_sv_GW | 384.6 |
| Fe | Fe_pv (\*) | 14 | 293.2 | 381 | Fe_sv_GW | 387.8 |
| Co | Co_pv (\*) | 15 | 271.0 | 352 | Co_sv_GW | 387.5 |
| Ni | Ni_pv (\*) | 16 | 368.0 | 478 | Ni_sv_GW | 389.6 |
| Cu | Cu_pv (\*) | 17 | 368.6 | 479 | Cu_sv_GW | 467.3 |
| Zn | Zn | 12 | 276.7 | 360 | Zn_sv_GW | 401.7 |
| Ga | Ga_d | 13 | 282.7 | 367 | Ga_d_GW | 404.6 |
| Ge | Ge_d | 14 | 310.3 | 403 | Ge_d_GW | 375.4 |
| As | As | 5 | 208.7 | 271 | As_GW | 208.7 |
| Se | Se | 6 | 211.6 | 275 | Se_GW | 211.6 |
| Br | Br | 7 | 216.3 | 281 | Br_GW | 216.3 |
| Kr | Kr | 8 | 185.3 | 241 | Kr_GW | 252.2 |
| Rb | Rb_sv | 9 | 220.1 | 286 | Rb_sv_GW | 221.2 |
| Sr | Sr_sv | 10 | 229.4 | 298 | Sr_sv_GW | 224.8 |
| Y | Y_sv | 11 | 202.6 | 263 | Y_sv_GW | 339.8 |
| Zr | Zr_sv | 12 | 229.9 | 299 | Zr_sv_GW | 346.4 |
| Nb | Nb_sv | 13 | 293.2 | 381 | Nb_sv_GW | 353.9 |
| Mo | Mo_sv | 14 | 242.7 | 315 | Mo_sv_GW | 344.9 |
| Tc | Tc_pv | 13 | 263.5 | 343 | Tc_sv_GW | 351.0 |
| Ru | Ru_pv | 14 | 240.0 | 312 | Ru_sv_GW | 348.1 |
| Rh | Rh_pv | 15 | 247.4 | 322 | Rh_sv_GW | 351.2 |
| Pd | Pd | 10 | 250.9 | 326 | Pd_sv_GW | 356.1 |
| Ag | Ag | 11 | 249.8 | 325 | Ag_sv_GW | 354.4 |
| Cd | Cd | 12 | 274.3 | 357 | Cd_sv_GW | 361.8 |
| In | In_d | 13 | 239.2 | 311 | In_d_GW | 278.6 |
| Sn | Sn_d | 14 | 241.1 | 313 | Sn_d_GW | 260.1 |
| Sb | Sb | 5 | 172.1 | 224 | Sb_d_GW | 263.1 |
| Te | Te | 6 | 175.0 | 227 | Te_GW | 175.0 |
| I | I | 7 | 175.6 | 228 | I_GW | 175.6 |
| Xe | Xe | 8 | 153.1 | 199 | Xe_GW | 179.5 |
| Cs | Cs_sv | 9 | 220.3 | 286 | Cs_sv_GW | 198.1 |
| Ba | Ba_sv | 10 | 187.2 | 243 | Ba_sv_GW | 267.0 |
| La | La | 11 | 219.3 | 285 | La_GW | 313.7 |
| Ce | Ce | 12 | 273.0 | 355 | Ce_GW | 304.6 |
| Pr | Pr_3 | 11 | 181.7 | 236 | - | - |
| Nd | Nd_3 | 11 | 182.6 | 237 | - | - |
| Pm | Pm_3 | 11 | 177.0 | 230 | - | - |
| Sm | Sm_3 | 11 | 177.1 | 230 | - | - |
| Eu | Eu_3 (\*) | 9 | 129.1 | 168 | - | - |
| Gd | Gd_3 | 9 | 154.3 | 201 | - | - |
| Tb | Tb_3 | 9 | 155.6 | 202 | - | - |
| Dy | Dy_3 | 9 | 155.7 | 202 | - | - |
| Ho | Ho_3 | 9 | 154.1 | 200 | - | - |
| Er | Er_3 | 9 | 155.0 | 202 | - | - |
| Tm | Tm_3 | 9 | 149.2 | 194 | - | - |
| Yb | Yb_3 (\*) | 9 | 188.4 | 245 | - | - |
| Lu | Lu_3 | 9 | 155.0 | 201 | - | - |
| Hf | Hf_pv | 10 | 220.3 | 286 | Hf_sv_GW | 309.0 |
| Ta | Ta_pv | 11 | 223.7 | 291 | Ta_sv_GW | 286.0 |
| W | W_sv | 14 | 223.1 | 290 | W_sv_GW | 317.1 |
| Re | Re_pv (\*) | 13 | 226.2 | 294 | Re_sv_GW | 317.0 |
| Os | Os_pv (\*) | 14 | 228.0 | 296 | Os_sv_GW | 319.8 |
| Ir | Ir | 9 | 210.9 | 274 | Ir_sv_GW | 319.8 |
| Pt | Pt | 10 | 230.3 | 299 | Pt_sv_GW | 323.7 |
| Au | Au | 11 | 229.9 | 299 | Au_sv_GW | 306.7 |
| Hg | Hg | 12 | 233.2 | 303 | Hg_sv_GW | 312.0 |
| Tl | Tl_d | 13 | 237.1 | 308 | Tl_d_GW | 237.1 |
| Pb | Pb_d | 14 | 237.8 | 309 | Pb_d_GW | 237.8 |
| Bi | Bi_d | 15 | 242.8 | 316 | Bi_d_GW | 261.9 |
| Po | Po_d | 16 | 264.6 | 344 | Po_d_GW | 267.8 |
| At | At | 7 | 161.4 | 210 | At_d_GW | 266.3 |
| Rn | Rn | 8 | 151.5 | 197 | Rn_d_GW | 267.3 |
| Fr | Fr_sv | 9 | 214.5 | 279 | - | - |
| Ra | Ra_sv | 10 | 237.4 | 309 | - | - |
| Ac | Ac | 11 | 172.4 | 224 | - | - |
| Th | Th | 12 | 247.3 | 321 | - | - |
| Pa | Pa | 13 | 252.2 | 328 | - | - |
| U | U | 14 | 252.5 | 328 | - | - |
| Np | Np | 15 | 254.3 | 331 | - | - |
| Pu | Pu | 16 | 254.4 | 331 | - | - |
| Am | Am | 17 | 255.9 | 333 | - | - |
| Cm | Cm | 18 | 258.0 | 335 | - | - |

## 使用提示

1. **ENMAX 取最大值**：ENCUT 默认等于 POTCAR 中所有元素 ENMAX 的最大值，但强烈建议在 INCAR 中手动写出 ENCUT。
2. **半芯态影响电子数**：例如 `Ti_sv` 的 ZVAL = 12，而 `Ti` 只有 4；NELECT、NBANDS、计算量都随之改变，Bader 电荷的“净电荷 = ZVAL − Bader 布居”也要用对应的 ZVAL。
3. **LDA 计算使用 LDA 赝势**：PBE 赝势适用于 PBE 及其他 GGA、meta-GGA 和杂化泛函；只有 LDA 计算才用 LDA 赝势。
4. **不要编辑 POTCAR**：换泛函用 INCAR 中的 `GGA`/`METAGGA`/`XC` 标签，而不是修改 POTCAR 里的 `LEXCH`。

## 数据来源

- Choosing pseudopotentials：https://www.vasp.at/wiki/index.php/Choosing_pseudopotentials
- Recommended pseudopotentials（含 SHA256 与历史推荐）：https://www.vasp.at/wiki/index.php/Recommended_pseudopotentials
- Available pseudopotentials：https://www.vasp.at/wiki/index.php/Available_pseudopotentials
