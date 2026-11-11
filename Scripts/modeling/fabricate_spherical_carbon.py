#!/usr/bin/env python3
"""
可定制碳原子球面片段坐标生成器
用户可以自定义曲率半径、原子数量和原子种类
输出为XYZ和CIF格式文件
"""

import numpy as np
import math
import argparse
from itertools import product


def get_user_input():
    """
    获取用户输入的参数
    """
    parser = argparse.ArgumentParser(description='生成球面片段碳原子坐标')

    parser.add_argument('--radius', type=float, default=40.0,
                        help='球面曲率半径（埃），默认800.0')
    parser.add_argument('--num_atoms', type=int, default=51,
                        help='原子数量，默认500')
    parser.add_argument('--element', type=str, default='C',
                        help='原子种类（元素符号），默认C')
    parser.add_argument('--bond_length', type=float, default=1.42,
                        help='原子间键长（埃），默认1.42（适合碳原子）')

    return parser.parse_args()


def generate_hexagonal_layer(center, radius, bond_length=1.42):
    """
    生成六边形网格的一层原子

    参数:
    center: 中心点坐标 (x, y, z)
    radius: 层半径
    bond_length: 原子间键长

    返回:
    该层原子的坐标列表
    """
    atoms = []

    # 六边形网格参数
    a = bond_length  # 六边形边长
    dx = a * math.sqrt(3)  # x方向间距
    dy = a * 1.5  # y方向间距

    # 计算网格尺寸
    n = int(radius / dx) + 1

    # 生成六边形网格
    for i, j in product(range(-n, n + 1), range(-n, n + 1)):
        x = i * dx
        y = j * dy + (i % 2) * dy / 2  # 交错排列

        # 检查是否在圆形区域内
        if math.sqrt(x ** 2 + y ** 2) <= radius:
            atoms.append((x + center[0], y + center[1], center[2]))

    return atoms


def project_to_sphere(atoms, sphere_radius=800.0):
    """
    将平面上的原子投影到球面上

    参数:
    atoms: 平面原子坐标列表
    sphere_radius: 球面半径

    返回:
    投影到球面上的原子坐标列表
    """
    spherical_atoms = []

    for x, y, z in atoms:
        # 计算点到原点的水平距离
        r_xy = math.sqrt(x ** 2 + y ** 2)

        # 如果点在球内，计算对应的球面坐标
        if r_xy <= sphere_radius:
            # 计算球面上的z坐标
            z_sphere = math.sqrt(sphere_radius ** 2 - r_xy ** 2)

            # 保持原始x,y方向，但调整z坐标
            spherical_atoms.append((x, y, z_sphere))

    return spherical_atoms


def generate_spherical_fragment(num_atoms, radius, bond_length=1.42):
    """
    生成球面片段原子坐标

    参数:
    num_atoms: 目标原子数量
    radius: 球面半径
    bond_length: 原子间键长

    返回:
    原子坐标列表
    """
    # 计算所需的层数和每层半径
    # 使用六边形密堆积的原子密度公式
    area_per_atom = (bond_length ** 2) * math.sqrt(3) / 2
    target_area = num_atoms * area_per_atom
    disk_radius = math.sqrt(target_area / math.pi)

    print(f"目标原子数: {num_atoms}")
    print(f"所需圆盘半径: {disk_radius:.2f} Å")

    # 生成平面六边形网格
    plane_atoms = generate_hexagonal_layer((0, 0, 0), disk_radius, bond_length)

    # 如果原子数过多，随机抽样到目标数量
    if len(plane_atoms) > num_atoms:
        indices = np.random.choice(len(plane_atoms), num_atoms, replace=False)
        plane_atoms = [plane_atoms[i] for i in indices]

    # 投影到球面上
    spherical_atoms = project_to_sphere(plane_atoms, radius)

    # 如果原子数不足，添加更多层
    while len(spherical_atoms) < num_atoms:
        # 增加圆盘半径
        disk_radius += bond_length
        additional_atoms = generate_hexagonal_layer((0, 0, 0), disk_radius, bond_length)

        # 过滤掉已经存在的原子
        existing_positions = set((round(x, 2), round(y, 2)) for x, y, z in spherical_atoms)
        new_atoms = []

        for atom in additional_atoms:
            pos_key = (round(atom[0], 2), round(atom[1], 2))
            if pos_key not in existing_positions:
                new_atoms.append(atom)
                existing_positions.add(pos_key)

        # 投影新原子到球面
        new_spherical = project_to_sphere(new_atoms, radius)
        spherical_atoms.extend(new_spherical)

        # 如果原子数超过目标，随机抽样
        if len(spherical_atoms) > num_atoms:
            indices = np.random.choice(len(spherical_atoms), num_atoms, replace=False)
            spherical_atoms = [spherical_atoms[i] for i in indices]
            break

    return spherical_atoms[:num_atoms]


def save_xyz_file(atoms, element, filename=None):
    """
    保存为XYZ格式文件

    参数:
    atoms: 原子坐标列表
    element: 原子元素符号
    filename: 输出文件名（如果为None则自动生成）
    """
    if filename is None:
        filename = f"F:\dft_C\lyy\_800A\{element}{len(atoms)}_spherical_r{int(radius)}.xyz"

    with open(filename, 'w') as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"{element}{len(atoms)} spherical fragment with curvature radius {radius} A\n")

        for i, (x, y, z) in enumerate(atoms):
            f.write(f"{element} {x:.6f} {y:.6f} {z:.6f}\n")

    print(f"XYZ文件已保存: {filename}")


def save_cif_file(atoms, element, radius, filename=None):
    """
    保存为CIF格式文件

    参数:
    atoms: 原子坐标列表
    element: 原子元素符号
    radius: 球面半径
    filename: 输出文件名（如果为None则自动生成）
    """
    if filename is None:
        filename = f"F:\dft_C\lyy\_800A\{element}{len(atoms)}_spherical_r{int(radius)}.cif"

    # 计算晶胞大小（略大于结构尺寸）
    max_coord = max(max(abs(x), abs(y), abs(z)) for x, y, z in atoms)
    cell_size = 2.2 * max_coord  # 留出足够空间

    with open(filename, 'w') as f:
        f.write("# Generated CIF file for spherical fragment\n")
        f.write(f"# Element: {element}, Atoms: {len(atoms)}, Radius: {radius} Angstrom\n")
        f.write(f"data_{element}{len(atoms)}_spherical_r{int(radius)}\n")
        f.write("\n")

        # 晶体学信息
        f.write("_cell_length_a       {:.4f}\n".format(cell_size))
        f.write("_cell_length_b       {:.4f}\n".format(cell_size))
        f.write("_cell_length_c       {:.4f}\n".format(cell_size))
        f.write("_cell_angle_alpha    90.0000\n")
        f.write("_cell_angle_beta     90.0000\n")
        f.write("_cell_angle_gamma    90.0000\n")
        f.write("_cell_volume         {:.4f}\n".format(cell_size ** 3))
        f.write("\n")

        f.write("_symmetry_space_group_name_H-M    'P 1'\n")
        f.write("_symmetry_Int_Tables_number       1\n")
        f.write("\n")

        # 原子坐标
        f.write("loop_\n")
        f.write("_atom_site_label\n")
        f.write("_atom_site_type_symbol\n")
        f.write("_atom_site_fract_x\n")
        f.write("_atom_site_fract_y\n")
        f.write("_atom_site_fract_z\n")
        f.write("_atom_site_occupancy\n")

        for i, (x, y, z) in enumerate(atoms):
            # 转换为分数坐标
            frac_x = (x + cell_size / 2) / cell_size
            frac_y = (y + cell_size / 2) / cell_size
            frac_z = (z + cell_size / 2) / cell_size

            f.write("{}{:4d} {} {:.6f} {:.6f} {:.6f} 1.000000\n".format(
                element, i + 1, element, frac_x, frac_y, frac_z))

    print(f"CIF文件已保存: {filename}")


def calculate_statistics(atoms, target_radius, element, bond_length):
    """
    计算结构统计信息

    参数:
    atoms: 原子坐标列表
    target_radius: 目标球面半径
    element: 原子元素符号
    bond_length: 原子间键长
    """
    # 计算所有原子到球心的距离
    distances = [math.sqrt(x ** 2 + y ** 2 + z ** 2) for x, y, z in atoms]
    avg_distance = np.mean(distances)
    std_distance = np.std(distances)

    # 计算曲率
    actual_curvature = 1.0 / avg_distance
    target_curvature = 1.0 / target_radius

    # 计算键长统计
    bond_lengths = []
    for i, (x1, y1, z1) in enumerate(atoms):
        for j, (x2, y2, z2) in enumerate(atoms):
            if i < j:  # 避免重复计算
                dist = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
                if bond_length * 0.9 < dist < bond_length * 1.1:  # 合理的键长范围
                    bond_lengths.append(dist)
                # 只检查前几个原子以提高效率
                if j > min(50, len(atoms)):
                    break

    avg_bond_length = np.mean(bond_lengths) if bond_lengths else 0
    std_bond_length = np.std(bond_lengths) if bond_lengths else 0

    print(f"\n结构统计信息:")
    print(f"目标曲率半径: {target_radius:.2f} Å")
    print(f"实际平均曲率半径: {avg_distance:.2f} Å (标准差: {std_distance:.4f} Å)")
    print(f"目标曲率: {target_curvature:.6f} Å^-1")
    print(f"实际平均曲率: {actual_curvature:.6f} Å^-1")
    print(f"曲率误差: {abs(actual_curvature - target_curvature) / target_curvature * 100:.2f}%")
    print(f"平均{element}-{element}键长: {avg_bond_length:.4f} Å (标准差: {std_bond_length:.4f} Å)")
    print(f"总原子数: {len(atoms)}")


def visualize_structure(atoms, element):
    """
    生成简单的ASCII可视化（用于快速检查结构）

    参数:
    atoms: 原子坐标列表
    element: 原子元素符号
    """
    # 创建简化的2D投影
    print(f"\n{element}原子结构俯视图 (简化投影):")
    print("=" * 50)

    # 将3D坐标投影到2D平面
    projected = [(x, y) for x, y, z in atoms]

    # 计算边界
    x_coords = [p[0] for p in projected]
    y_coords = [p[1] for p in projected]
    x_min, x_max = min(x_coords), max(x_coords)
    y_min, y_max = min(y_coords), max(y_coords)

    # 创建简化的ASCII网格
    grid_size = min(60, len(atoms) // 5)  # 根据原子数调整网格大小
    grid = [[' ' for _ in range(grid_size)] for _ in range(grid_size // 2)]

    # 将原子映射到网格
    for x, y in projected:
        i = int((y - y_min) / (y_max - y_min) * (len(grid) - 1))
        j = int((x - x_min) / (x_max - x_min) * (len(grid[0]) - 1))
        if 0 <= i < len(grid) and 0 <= j < len(grid[0]):
            grid[i][j] = '•'

    # 打印网格
    for row in grid:
        print(''.join(row))


def main():
    """主函数"""
    print("可定制球面片段原子坐标生成器")
    print("=" * 50)

    # 获取用户输入
    args = get_user_input()
    global radius
    radius = args.radius
    num_atoms = args.num_atoms
    element = args.element
    bond_length = args.bond_length

    print(f"生成参数:")
    print(f"  原子种类: {element}")
    print(f"  原子数量: {num_atoms}")
    print(f"  曲率半径: {radius} Å")
    print(f"  键长: {bond_length} Å")
    print()

    # 验证输入
    if num_atoms <= 0:
        print("错误: 原子数量必须大于0")
        return

    if radius <= 0:
        print("错误: 曲率半径必须大于0")
        return

    if bond_length <= 0:
        print("错误: 键长必须大于0")
        return

    # 生成原子坐标
    atoms = generate_spherical_fragment(num_atoms, radius, bond_length)

    # 计算统计信息
    calculate_statistics(atoms, radius, element, bond_length)

    # 保存文件
    save_xyz_file(atoms, element)
    save_cif_file(atoms, element, radius)

    # 显示结构信息
    print(f"\n结构尺寸:")
    x_coords = [atom[0] for atom in atoms]
    y_coords = [atom[1] for atom in atoms]
    z_coords = [atom[2] for atom in atoms]

    x_range = max(x_coords) - min(x_coords)
    y_range = max(y_coords) - min(y_coords)
    z_range = max(z_coords) - min(z_coords)

    print(f"X方向范围: {x_range:.2f} Å")
    print(f"Y方向范围: {y_range:.2f} Å")
    print(f"Z方向范围: {z_range:.2f} Å")
    print(f"近似直径: {max(x_range, y_range):.2f} Å")

    # 简单可视化（仅对较小结构）
    if num_atoms <= 1000:
        visualize_structure(atoms, element)

    print("\n生成完成！")
    print("\n使用建议:")
    print("1. 用VMD、PyMOL等软件查看XYZ文件")
    print("2. 用VESTA、Mercury等软件查看CIF文件")
    print("3. 可用于分子动力学模拟或量子化学计算")
    print("\n命令行使用示例:")
    print("  python script.py --radius 500 --num_atoms 300 --element Si --bond_length 2.35")


if __name__ == "__main__":
    main()