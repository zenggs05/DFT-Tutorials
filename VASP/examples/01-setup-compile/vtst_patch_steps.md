# VTST 补丁编译步骤（VASP 6.4.3–6.6.1 + vtstcode）

本文件是 [第1章 第9节](../../tutorials/01-setup-compile.md) 的操作清单，依据 VTST 官方安装说明（https://theory.cm.utexas.edu/vtsttools/installation.html ）整理，只描述“做什么”，不包含 VASP 源码。

## 0. 准备

```bash
cd ~/software
git clone https://github.com/henkelmangroup/vtstcode.git     # Apache-2.0
ls vtstcode                                                   # vtstcode6.4.3 vtstcode6.5.1 vtstcode6.6.0 vtstcode6.6.1 ...
cd ~/software/vasp.6.6.1
cp -r src src.orig                                            # 备份，出问题时可以 diff 或还原
```

> 选择**与 VASP 版本一致**的 vtstcode 目录。若没有完全一致的目录，优先使用版本号最接近且不高于你的 VASP 的目录，并预期可能需要手动排错。

## 1. 复制 VTST 文件

```bash
cp -r ../vtstcode/vtstcode6.6.1/* src/
ls src/{neb,dimer,lanczos,dynmat,opt,fire,lbfgs,ml_pyamff}.F
ls -d src/pyamff_fortran
```

`chain.F` 会被 VTST 版本覆盖；其余为新增文件。

## 2. 修改 `src/.objects`

```bash
grep -n "chain.o" src/.objects
```

在 `SOURCE` 列表中 `chain.o` 那一行**之前**插入（每行末尾的反斜杠保持与文件中其他行一致）：

```text
	bfgs.o dynmat.o instanton.o lbfgs.o sd.o cg.o dimer.o bbm.o \
	fire.o lanczos.o neb.o qm.o \
	pyamff_fortran/*.o ml_pyamff.o \
	opt.o \
```

规则：`dynmat.o neb.o dimer.o lanczos.o instanton.o` 在 `chain.o` 之前；`sd.o cg.o qm.o lbfgs.o bfgs.o fire.o ml_pyamff.o pyamff_fortran/*.o` 在 `opt.o` 之前。

## 3. 修改 `src/main.F`（两处）

```bash
grep -n -A1 "CALL CHAIN_FORCE" src/main.F
grep -n "chain_init" src/main.F
```

1. `CALL CHAIN_FORCE(...)`：在参数 `TIFOR,` 与 `LATT_CUR%A,` 之间插入 `TSIF,`（续行符 `&` 的位置保持合法）。
2. `IF (LCHAIN) CALL chain_init( T_INFO, IO)` → `CALL chain_init( T_INFO, IO)`。

## 4. 修改 `src/makefile`

```bash
grep -n "^LIB\|^dependencies" src/makefile
```

- 在 `LIB` 变量末尾追加 `pyamff_fortran`，例如 `LIB= lib parser` → `LIB= lib parser pyamff_fortran`；
- 并行编译时把依赖规则改为 `dependencies: sources libs`。

## 5. 自动检查（只读）

```bash
bash /path/to/DFT/VASP/examples/01-setup-compile/check_vtst_patch.sh ~/software/vasp.6.6.1
```

全部 `[PASS]` 后再编译；任何 `[FAIL]` 都会导致编译或运行出错。

## 6. 编译与验证

```bash
make veryclean
make DEPS=1 -j 16 all 2>&1 | tee make_vtst.log
grep -i -c "error" make_vtst.log
nm bin/vasp_std | grep -i -c dimer     # > 0
nm bin/vasp_std | grep -i -c fire      # > 0
```

建议把带 VTST 的程序装到单独目录（如 `~/software/vasp.6.6.1-vtst/bin`）或改名（`vasp_std_vtst`），避免与原版混淆。功能验证见 [第11章 过渡态搜索](../../tutorials/11-neb-dimer.md)：INCAR 中 `IMAGES` 或 `ICHAIN` 设置后才会进入 VTST 代码；VTST 优化器需要 `IBRION = 3`、`POTIM = 0` 与 `IOPT`。

## 7. 安装 VTST 脚本

```bash
git clone https://github.com/henkelmangroup/vtstscripts.git ~/software/vtstscripts
echo 'export PATH=$HOME/software/vtstscripts:$PATH' >> ~/.bashrc
which nebmake.pl nebresults.pl chgsum.pl      # 需要系统有 perl
```

## 常见错误

| 现象 | 原因 | 处理 |
|---|---|---|
| 编译 main.F 时 CHAIN_FORCE 参数不匹配 | 第 3 步第 1 处未改 | 加入 `TSIF,` |
| 运行 NEB 时 VTST 代码没被初始化 | 第 3 步第 2 处未改 | 去掉 `IF (LCHAIN)` |
| 找不到 PyAMFF 模块/对象 | `.objects` 顺序错或 `LIB` 未加 `pyamff_fortran` | 按第 2、4 步修改 |
| `make -j` 时偶发找不到库 | 依赖规则缺 `libs` | `dependencies: sources libs` |
| 与 VASP 版本不匹配的各种编译错误 | vtstcode 目录版本不对 | 换用对应版本目录，从 `src.orig` 重新开始 |
