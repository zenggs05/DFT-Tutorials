import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from ase.io.cube import read_cube_data


def get_vacuum_level_and_export(cube_file_path):
    print(f"正在读取 {cube_file_path} ... 这可能需要几十秒钟。")
    data, atoms = read_cube_data(cube_file_path)

    planar_avg_hartree = np.mean(data, axis=(0, 1))

    planar_avg_ev = planar_avg_hartree * 27.2114

    cell_z = atoms.cell[2, 2]
    nz = data.shape[2]
    z_axis = np.linspace(0, cell_z, nz)

    e_vac = np.max(planar_avg_ev)
    print(f"成功! 计算得到的真空能级 (Evac) 为: {e_vac:.6f} eV\n")


    print("正在导出数据文件...")

    potential_data = np.column_stack((z_axis, planar_avg_ev))
    np.savetxt('Electrostatic_Potential_1D.txt', potential_data,
               fmt='%.6f', delimiter='\t',
               header='Z_axis(Angstrom)\tPotential(eV)', comments='')
    print(" ➔ 静电势曲线数据已保存至: Electrostatic_Potential_1D.txt (可直接拖入 Origin 绘图)")

    with open('Vacuum_Level.txt', 'w') as f:
        f.write("Vacuum Level (Evac) in eV:\n")
        f.write(f"{e_vac:.6f}\n")
    print(" ➔ 真空能级数值已单独保存至: Vacuum_Level.txt\n")

    plt.figure(figsize=(8, 5))
    plt.plot(z_axis, planar_avg_ev, color='blue', linewidth=2, label='Planar Average Potential')
    plt.axhline(y=e_vac, color='red', linestyle='--', label=f'Vacuum Level: {e_vac:.4f} eV')

    plt.xlabel('Z axis ($\AA$)', fontsize=12)
    plt.ylabel('Electrostatic Potential (eV)', fontsize=12)
    plt.title('Macroscopic Average Electrostatic Potential', fontsize=14)
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.tight_layout()
    plt.savefig('Electrostatic_Potential.png', dpi=300)
    print(" ➔ 预览曲线图已保存至: Electrostatic_Potential.png")

    plt.show()


if __name__ == "__main__":
    cube_filename = "Ni9_LDH_work-function-cube-v_hartree-1_0.cube"
    get_vacuum_level_and_export(cube_filename)