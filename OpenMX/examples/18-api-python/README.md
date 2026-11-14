# 示例 18-4：用 OpenMX 4.0 Python API 在一个 MPI 作业中并行跑多个输入

配套教程：[第18章 并行、大体系 O(N)、作业脚本与常见错误排查](../../tutorials/18-parallel-on-errors.md)（§9）

## 目的

OpenMX 4.0 新增了 API（手册第 71 章），可以从 C、Fortran 或 Python 程序中把 OpenMX 当作库调用：`openmx_api(comm, input_file, nthreads)`。本例按官方示例 `example_openmx_api.py` 的思路，把 `MPI_COMM_WORLD` 按输入文件数切分成若干子通信器，每个子通信器各跑一个输入，适合“很多个互不相关的小计算”（结构筛选、有限位移、参数扫描），一次排队即可全部完成。

## 文件说明

| 文件 | 说明 |
|---|---|
| `api_batch.py` | 主程序：`mpirun -np 8 python3 api_batch.py --nt 1 a.dat b.dat c.dat d.dat`。进程按连续块分给各输入；进程数必须不少于输入数；结束时在 rank 0 汇总每个输入的返回码 |
| `job_api_slurm.sh` | Slurm 模板：默认运行当前目录下所有 `*.dat`，结束后列出每个输入的 `Utot` |

输入文件请自备，例如把 `openmx4.0/work/Methane.dat`、`work/Cdia.dat`（官方 API 示例用的两个输入）复制过来，并把其中的 `DATA.PATH` 改成绝对路径。

## 运行步骤

```bash
# 1. build the Python extension once (PIC objects are required, hence make clean)
cd ~/openmx4.0/source
make clean
make pyapi                                   # -> pyopenmx_api*.so
# make clean works in source/ (objects); the copy ../work/openmx made by make install is normally kept.
# If you later need to rebuild openmx itself: make install

# 2. Python environment with mpi4py linked to the SAME MPI as OpenMX
python3 -c "import mpi4py; print(mpi4py.get_config())"

# 3. run
cd /path/to/OpenMX/examples/18-api-python
cp ~/openmx4.0/source/pyopenmx_api*.so .
cp ~/openmx4.0/work/Methane.dat ~/openmx4.0/work/Cdia.dat .
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" Methane.dat Cdia.dat
mpirun -np 4 python3 api_batch.py --nt 1 Methane.dat Cdia.dat
#   or: sbatch job_api_slurm.sh Methane.dat Cdia.dat
```

## 预期结果 / 如何检查

- 屏幕上先打印每个组的分配情况（如 `group 0: 2 ranks x 1 threads -> Methane.dat`），最后打印每个输入的 `ierr`，全部为 0 表示成功。
- 每个输入照常生成自己的 `System.Name.out` 等文件；用 `grep "Utot\." *.out` 检查总能，与单独用 `mpirun openmx` 运行同一输入的结果应一致。
- 各组的标准输出会交错打印在同一个终端/日志中，这是正常现象；判断结果以各自的 `.out` 为准。

## 注意事项

- **MPI 一致性最关键**：编译 OpenMX、启动用的 `mpirun` 和 mpi4py 链接的 MPI 必须是同一套，否则可能在启动后立即中止（手册第 71.4 节）。
- `make pyapi` 前必须先 `make clean`，因为普通编译得到的目标文件不是 `-fPIC` 的；`make clean` 删除 `source/` 中已编译的目标文件（`make install` 复制到 `work/` 的 `openmx` 一般不受影响），之后需要重编 `openmx` 时再 `make install`。
- `python3-config: command not found` 或 `Python.h` 找不到：安装 Python 开发包（如 `python3-dev`），或改用对应版本的 `python3.x-config`。
- 线程数由 API 参数（`--nt`）传入，但实际可用线程数仍受 `OMP_NUM_THREADS` 和调度器分配的核数限制。
- 本脚本让每个进程只调用一次 API（与官方示例相同）。输入数多于进程数时，请分批提交或使用作业数组，而不是在同一进程里反复调用。
- 不同输入必须使用不同的 `System.Name`，否则输出文件互相覆盖。

返回：[第18章](../../tutorials/18-parallel-on-errors.md) | [OpenMX 目录](../../README.md)
