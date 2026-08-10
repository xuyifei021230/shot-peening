# -*- coding: utf-8 -*-
"""
Abaqus 2025：从静态分析 ODB 自动测量单颗喷丸弹坑直径/半径

V3 测量定义：
1. 读取靶材最终帧顶部表面节点；
2. 使用远离弹坑的节点拟合变形后的参考平面；
3. 计算顶部节点相对参考平面的凹陷深度；
4. 沿通过冲击中心的 X、Y 两条剖面提取原始 U3 曲线；
5. 将中心左、右两侧各自的原始 U3 最高点定义为弹坑的两个边缘点；
6. 按变形后表面路径累计距离计算两峰间距，其一半作为弹坑半径。

运行：
Abaqus/CAE -> 文件 -> 运行脚本
或：abaqus python measure_crater_from_odb_abaqus2025.py
"""

from __future__ import print_function

import os
import csv
import math
from odbAccess import openOdb


# =============================================================================
# 1. 用户参数区
# =============================================================================

# True：使用 Abaqus 当前工作目录；False：使用 WORK_DIRECTORY。
USE_CURRENT_WORK_DIRECTORY = True
WORK_DIRECTORY = r"D:\temp\one-shot"

ODB_FILE_NAME = "Job-1-Static.odb"
INSTANCE_NAME = "TARGET-1"      # 忽略大小写
STEP_NAME = None                 # None：最后一个分析步；也可写 "Step-1"
FRAME_INDEX = -1                 # -1：最后一帧

# 冲击中心原始坐标；设为 None 时按最深凹陷位置自动识别。
IMPACT_CENTER_X = 0.0
IMPACT_CENTER_Y = 0.0

# None：自动取靶材较短边尺寸的 35%；2 mm 靶材约为 0.7 mm。
FAR_FIELD_RADIUS = None

# 中心剖面带宽。中心细网格尺寸约 0.01 mm 时，建议 0.005～0.008 mm。
PROFILE_HALF_WIDTH = 0.006

# 弹坑边缘固定采用中心两侧的原始 U3 最高点。
EDGE_METHOD = "U3_PEAKS"

# 移动平均仅作为剖面 CSV 中的辅助结果，不参与弹坑边缘和半径判断。
SMOOTH_WINDOW = 5

# True：生成 X/Y 剖面和顶部表面节点 CSV；False：不生成任何 CSV。
# 无论此开关为何值，测量摘要 TXT 都会始终生成。
WRITE_CSV_FILES = False

SUMMARY_FILE_NAME = "crater-job-1-05MPa_measurement_summary.txt"
X_PROFILE_FILE_NAME = "crater-job-1-05MPa-1_profile_X.csv"
Y_PROFILE_FILE_NAME = "crater-job-1-05MPa-1_profile_Y.csv"
TOP_SURFACE_FILE_NAME = "crate-job-1-05MPar-1_top_surface_nodes.csv"


# =============================================================================
# 2. 通用函数
# =============================================================================

def get_work_directory():
    directory = os.getcwd() if USE_CURRENT_WORK_DIRECTORY else os.path.abspath(WORK_DIRECTORY)
    if not os.path.isdir(directory):
        raise IOError("工作目录不存在：{}".format(directory))
    os.chdir(directory)
    return directory


def find_key_case_insensitive(mapping, requested_name):
    requested_upper = requested_name.upper()
    for key in mapping.keys():
        if key.upper() == requested_upper:
            return key
    raise KeyError("未找到 '{}'; 可用名称：{}".format(requested_name, list(mapping.keys())))


def solve_3x3(matrix, vector):
    a = [[float(matrix[i][j]) for j in range(3)] + [float(vector[i])] for i in range(3)]
    for column in range(3):
        pivot = max(range(column, 3), key=lambda row: abs(a[row][column]))
        if abs(a[pivot][column]) < 1.0e-20:
            raise ValueError("参考平面拟合失败：矩阵接近奇异。")
        if pivot != column:
            a[column], a[pivot] = a[pivot], a[column]
        pivot_value = a[column][column]
        for j in range(column, 4):
            a[column][j] /= pivot_value
        for row in range(3):
            if row == column:
                continue
            factor = a[row][column]
            for j in range(column, 4):
                a[row][j] -= factor * a[column][j]
    return a[0][3], a[1][3], a[2][3]


def fit_plane(points):
    """
    最小二乘拟合 z=a*x+b*y+c。

    注意：
    在 Abaqus/CAE 中通过“运行脚本”执行时，名称 sum 有时会被
    Abaqus 自身的函数覆盖。若直接使用 sum(generator)，会出现：
        TypeError: 已找到 'generator', 应当提供 a recognized type

    因此这里使用显式循环累计，避免调用可能被覆盖的 sum。
    """
    if len(points) < 3:
        raise ValueError("参考平面拟合至少需要 3 个点。")

    s_xx = 0.0
    s_yy = 0.0
    s_xy = 0.0
    s_x = 0.0
    s_y = 0.0
    s_z = 0.0
    s_xz = 0.0
    s_yz = 0.0

    for point in points:
        x = float(point[0])
        y = float(point[1])
        z = float(point[2])

        s_xx += x * x
        s_yy += y * y
        s_xy += x * y
        s_x += x
        s_y += y
        s_z += z
        s_xz += x * z
        s_yz += y * z

    n = float(len(points))

    matrix = (
        (s_xx, s_xy, s_x),
        (s_xy, s_yy, s_y),
        (s_x, s_y, n),
    )
    vector = (s_xz, s_yz, s_z)

    return solve_3x3(matrix, vector)


def moving_average(values, window):
    """
    移动平均。

    同样不使用 sum()，避免 Abaqus/CAE 全局命名空间中的 sum
    覆盖 Python 内置函数。
    """
    if window <= 1:
        return list(values)

    window = int(window)
    if window % 2 == 0:
        window += 1

    half = window // 2
    result = []

    for i in range(len(values)):
        start = max(0, i - half)
        end = min(len(values), i + half + 1)

        subtotal = 0.0
        count = 0

        for j in range(start, end):
            subtotal += float(values[j])
            count += 1

        result.append(subtotal / float(count))

    return result


def linear_intersection(x1, y1, x2, y2, target):
    denominator = y2 - y1
    if abs(denominator) < 1.0e-30:
        return 0.5 * (x1 + x2)
    ratio = (target - y1) / denominator
    return x1 + ratio * (x2 - x1)


def group_profile(records, axis, center_x, center_y, half_width):
    axis = axis.upper()
    if axis == "X":
        selected = [r for r in records if abs(r["y0"] - center_y) <= half_width]
        coordinate_name, group_name = "x_def", "x0"
    elif axis == "Y":
        selected = [r for r in records if abs(r["x0"] - center_x) <= half_width]
        coordinate_name, group_name = "y_def", "y0"
    else:
        raise ValueError("axis 必须为 X 或 Y。")

    if len(selected) < 5:
        raise ValueError("{} 向剖面节点过少（{} 个），请增大 PROFILE_HALF_WIDTH。".format(axis, len(selected)))

    groups = {}
    for rec in selected:
        groups.setdefault(round(rec[group_name], 12), []).append(rec)

    profile = []
    for key in sorted(groups.keys()):
        group = groups[key]
        coordinate_total = 0.0
        original_coordinate_total = 0.0
        x_deformed_total = 0.0
        y_deformed_total = 0.0
        z_deformed_total = 0.0
        u3_total = 0.0
        depth_total = 0.0
        relative_height_total = 0.0

        for record in group:
            coordinate_total += float(record[coordinate_name])
            original_coordinate_total += float(record[group_name])
            x_deformed_total += float(record["x_def"])
            y_deformed_total += float(record["y_def"])
            z_deformed_total += float(record["z_def"])
            u3_total += float(record["u3"])
            depth_total += float(record["depth"])
            relative_height_total += float(record["relative_height"])

        group_count = float(len(group))

        profile.append({
            "coordinate": coordinate_total / group_count,
            "original_coordinate": original_coordinate_total / group_count,
            "x_deformed": x_deformed_total / group_count,
            "y_deformed": y_deformed_total / group_count,
            "z_deformed": z_deformed_total / group_count,
            "u3_raw": u3_total / group_count,
            "depth_raw": depth_total / group_count,
            "relative_height_raw": relative_height_total / group_count,
            "node_count": len(group),
            "edge_role": "",
        })

    # 与 Abaqus 路径 XY 数据的 TRUE_DISTANCE 定义一致：沿变形后的表面剖面
    # 累加相邻采样点的三维欧氏距离。弹坑半径使用两个 U3 峰值点之间的
    # 这段路径距离之半，而不是只使用 X/Y 坐标差。
    cumulative_distance = 0.0
    for index, item in enumerate(profile):
        if index > 0:
            previous = profile[index - 1]
            dx = item["x_deformed"] - previous["x_deformed"]
            dy = item["y_deformed"] - previous["y_deformed"]
            dz = item["z_deformed"] - previous["z_deformed"]
            cumulative_distance += math.sqrt(dx * dx + dy * dy + dz * dz)
        item["path_distance"] = cumulative_distance

    smooth = moving_average([item["depth_raw"] for item in profile], SMOOTH_WINDOW)
    for item, value in zip(profile, smooth):
        item["depth_smooth"] = value
        item["relative_height_smooth"] = -value
    return profile


def find_u3_peak_edges(profile, center_coordinate):
    """以中心两侧各自的原始 U3 最高点作为坑口边缘。"""
    coordinates = [item["coordinate"] for item in profile]
    center_index = min(range(len(coordinates)), key=lambda i: abs(coordinates[i] - center_coordinate))

    if center_index < 1 or center_index >= len(profile) - 1:
        raise ValueError("U3 最低点位于剖面端部，无法在中心两侧寻找坑口峰值。")

    # 若出现相同峰值，优先选择更靠近弹坑中心的点。
    left_index = max(
        range(0, center_index),
        key=lambda i: (profile[i]["u3_raw"], i),
    )
    right_index = max(
        range(center_index + 1, len(profile)),
        key=lambda i: (profile[i]["u3_raw"], -i),
    )

    left = profile[left_index]
    right = profile[right_index]
    left["edge_role"] = "LEFT_U3_PEAK"
    right["edge_role"] = "RIGHT_U3_PEAK"

    diameter = right["path_distance"] - left["path_distance"]
    if diameter <= 0.0:
        raise ValueError("两个 U3 峰值点的路径距离无效，请检查中心剖面和节点顺序。")

    return {
        "left_coordinate": left["coordinate"],
        "right_coordinate": right["coordinate"],
        "left_original_coordinate": left["original_coordinate"],
        "right_original_coordinate": right["original_coordinate"],
        "left_path_distance": left["path_distance"],
        "right_path_distance": right["path_distance"],
        "left_u3": left["u3_raw"],
        "right_u3": right["u3_raw"],
        "diameter": diameter,
        "radius": 0.5 * diameter,
    }


def write_profile_csv(path, profile, axis_name):
    with open(path, "w", newline="") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow([axis_name + "_original", axis_name + "_deformed",
                         "deformed_path_distance", "u3_raw", "depth_raw", "depth_smoothed",
                         "relative_height_raw", "relative_height_smoothed",
                         "averaged_node_count", "edge_role"])
        for item in profile:
            writer.writerow(["{:.12e}".format(item["original_coordinate"]),
                             "{:.12e}".format(item["coordinate"]),
                             "{:.12e}".format(item["path_distance"]),
                             "{:.12e}".format(item["u3_raw"]),
                             "{:.12e}".format(item["depth_raw"]),
                             "{:.12e}".format(item["depth_smooth"]),
                             "{:.12e}".format(item["relative_height_raw"]),
                             "{:.12e}".format(item["relative_height_smooth"]),
                             item["node_count"], item["edge_role"]])


# =============================================================================
# 3. ODB 数据读取与测量
# =============================================================================

def main():
    work_directory = get_work_directory()
    odb_path = os.path.join(work_directory, ODB_FILE_NAME)
    if not os.path.isfile(odb_path):
        raise IOError("未找到 ODB 文件：{}".format(odb_path))

    print("=" * 78)
    print("Abaqus ODB 弹坑尺寸自动测量")
    print("=" * 78)
    print("工作目录：{}".format(work_directory))
    print("ODB 文件：{}".format(odb_path))

    odb = openOdb(path=odb_path, readOnly=True)
    try:
        instance_key = find_key_case_insensitive(odb.rootAssembly.instances, INSTANCE_NAME)
        instance = odb.rootAssembly.instances[instance_key]
        step_key = list(odb.steps.keys())[-1] if STEP_NAME is None else find_key_case_insensitive(odb.steps, STEP_NAME)
        step = odb.steps[step_key]
        frame = step.frames[FRAME_INDEX]

        if "U" not in frame.fieldOutputs:
            raise KeyError("当前帧没有位移场 U，请确认场输出请求包含 U。")

        displacement_by_label = {}
        for value in frame.fieldOutputs["U"].values:
            if value.instance is not None and value.instance.name.upper() == instance_key.upper():
                displacement_by_label[value.nodeLabel] = tuple(value.data)
        if not displacement_by_label:
            raise ValueError("没有读取到实例 {} 的节点位移。".format(instance_key))

        nodes = list(instance.nodes)
        z_values = [node.coordinates[2] for node in nodes]
        z_min, z_max = min(z_values), max(z_values)
        thickness = z_max - z_min
        top_tolerance = max(1.0e-9, 1.0e-6 * max(1.0, abs(thickness)))
        top_nodes = [node for node in nodes
                     if abs(node.coordinates[2] - z_max) <= top_tolerance
                     and node.label in displacement_by_label]
        if len(top_nodes) < 10:
            raise ValueError("顶部节点过少（{}），请检查厚度方向是否为全局 Z。".format(len(top_nodes)))

        x_values = [node.coordinates[0] for node in top_nodes]
        y_values = [node.coordinates[1] for node in top_nodes]
        x_mid = 0.5 * (min(x_values) + max(x_values))
        y_mid = 0.5 * (min(y_values) + max(y_values))
        provisional_x = x_mid if IMPACT_CENTER_X is None else IMPACT_CENTER_X
        provisional_y = y_mid if IMPACT_CENTER_Y is None else IMPACT_CENTER_Y
        far_radius = (0.35 * min(max(x_values) - min(x_values), max(y_values) - min(y_values))
                      if FAR_FIELD_RADIUS is None else float(FAR_FIELD_RADIUS))

        records = []
        for node in top_nodes:
            u = displacement_by_label[node.label]
            x0, y0, z0 = node.coordinates
            records.append({
                "label": node.label, "x0": x0, "y0": y0, "z0": z0,
                "u1": u[0], "u2": u[1], "u3": u[2],
                "x_def": x0 + u[0], "y_def": y0 + u[1], "z_def": z0 + u[2],
            })

        def fit_reference_plane(center_x, center_y):
            far_points = [(r["x_def"], r["y_def"], r["z_def"]) for r in records
                          if math.hypot(r["x0"] - center_x, r["y0"] - center_y) >= far_radius]
            if len(far_points) < 20:
                raise ValueError("远场节点过少（{}），请减小 FAR_FIELD_RADIUS。".format(len(far_points)))
            return fit_plane(far_points), len(far_points)

        def update_depths(plane_coefficients):
            a, b, c = plane_coefficients
            for rec in records:
                reference_z = a * rec["x_def"] + b * rec["y_def"] + c
                rec["reference_z"] = reference_z
                rec["relative_height"] = rec["z_def"] - reference_z
                rec["depth"] = reference_z - rec["z_def"]

        plane, far_count = fit_reference_plane(provisional_x, provisional_y)
        update_depths(plane)
        deepest = max(records, key=lambda r: r["depth"])
        center_x = deepest["x0"] if IMPACT_CENTER_X is None else float(IMPACT_CENTER_X)
        center_y = deepest["y0"] if IMPACT_CENTER_Y is None else float(IMPACT_CENTER_Y)

        plane, far_count = fit_reference_plane(center_x, center_y)
        update_depths(plane)
        deepest = max(records, key=lambda r: r["depth"])
        maximum_depth = deepest["depth"]
        if maximum_depth <= 0.0:
            raise ValueError("最大凹陷深度不为正，请检查实例、表面方向和帧。")

        x_profile = group_profile(records, "X", center_x, center_y, PROFILE_HALF_WIDTH)
        y_profile = group_profile(records, "Y", center_x, center_y, PROFILE_HALF_WIDTH)

        used_method = EDGE_METHOD.strip().upper()
        if used_method != "U3_PEAKS":
            raise ValueError("V3 的 EDGE_METHOD 必须为 U3_PEAKS。")

        center_x_def = center_x + deepest["u1"]
        center_y_def = center_y + deepest["u2"]
        x_edges = find_u3_peak_edges(x_profile, center_x_def)
        y_edges = find_u3_peak_edges(y_profile, center_y_def)

        x_left = x_edges["left_coordinate"]
        x_right = x_edges["right_coordinate"]
        y_lower = y_edges["left_coordinate"]
        y_upper = y_edges["right_coordinate"]
        diameter_x = x_edges["diameter"]
        diameter_y = y_edges["diameter"]
        radius_x = x_edges["radius"]
        radius_y = y_edges["radius"]
        mean_diameter = 0.5 * (diameter_x + diameter_y)
        mean_radius = 0.5 * mean_diameter
        equivalent_ellipse_diameter = math.sqrt(diameter_x * diameter_y)
        ellipticity = max(diameter_x, diameter_y) / min(diameter_x, diameter_y)
        a, b, c = plane

        summary_lines = [
            "=" * 78, "弹坑尺寸测量结果", "=" * 78,
            "ODB 文件：{}".format(ODB_FILE_NAME),
            "实例：{}".format(instance_key),
            "分析步：{}".format(step_key),
            "帧时间：{:.12e}".format(frame.frameValue), "",
            "顶部节点数：{}".format(len(top_nodes)),
            "参考平面远场节点数：{}".format(far_count),
            "远场半径：{:.12e}".format(far_radius),
            "参考平面：z=({:.12e})*x+({:.12e})*y+({:.12e})".format(a, b, c), "",
            "冲击中心原始坐标：X={:.12e}, Y={:.12e}".format(center_x, center_y),
            "最深节点：{}".format(deepest["label"]),
            "最大弹坑深度：{:.12e}".format(maximum_depth), "",
            "边缘方法：{}".format(used_method),
            "边缘定义：中心左右两侧的原始 U3 最高点",
            "半径定义：两个 U3 峰值点在变形后表面路径上的距离之半",
            "注意：移动平均结果不参与边缘判断", "",
            "X 向左峰值变形坐标：{:.12e}".format(x_left),
            "X 向右峰值变形坐标：{:.12e}".format(x_right),
            "X 向左峰值路径距离：{:.12e}".format(x_edges["left_path_distance"]),
            "X 向右峰值路径距离：{:.12e}".format(x_edges["right_path_distance"]),
            "X 向左峰值 U3：{:.12e}".format(x_edges["left_u3"]),
            "X 向右峰值 U3：{:.12e}".format(x_edges["right_u3"]),
            "X 向直径：{:.12e}".format(diameter_x),
            "X 向半径：{:.12e}".format(radius_x), "",
            "Y 向下峰值变形坐标：{:.12e}".format(y_lower),
            "Y 向上峰值变形坐标：{:.12e}".format(y_upper),
            "Y 向下峰值路径距离：{:.12e}".format(y_edges["left_path_distance"]),
            "Y 向上峰值路径距离：{:.12e}".format(y_edges["right_path_distance"]),
            "Y 向下峰值 U3：{:.12e}".format(y_edges["left_u3"]),
            "Y 向上峰值 U3：{:.12e}".format(y_edges["right_u3"]),
            "Y 向直径：{:.12e}".format(diameter_y),
            "Y 向半径：{:.12e}".format(radius_y), "",
            "X/Y 平均直径：{:.12e}".format(mean_diameter),
            "X/Y 平均半径：{:.12e}".format(mean_radius),
            "椭圆等效直径 sqrt(Dx*Dy)：{:.12e}".format(equivalent_ellipse_diameter),
            "椭圆率 max(Dx,Dy)/min(Dx,Dy)：{:.12e}".format(ellipticity),
            "=" * 78,
        ]

        summary_path = os.path.join(work_directory, SUMMARY_FILE_NAME)
        with open(summary_path, "w", encoding="utf-8") as summary_file:
            summary_file.write("\n".join(summary_lines))

        if WRITE_CSV_FILES:
            write_profile_csv(os.path.join(work_directory, X_PROFILE_FILE_NAME), x_profile, "X")
            write_profile_csv(os.path.join(work_directory, Y_PROFILE_FILE_NAME), y_profile, "Y")

            with open(os.path.join(work_directory, TOP_SURFACE_FILE_NAME), "w", newline="") as csv_file:
                writer = csv.writer(csv_file)
                writer.writerow(["node_label", "x0", "y0", "z0", "u1", "u2", "u3",
                                 "x_deformed", "y_deformed", "z_deformed", "reference_z",
                                 "relative_height", "depth"])
                for rec in records:
                    writer.writerow([rec["label"],
                                     "{:.12e}".format(rec["x0"]), "{:.12e}".format(rec["y0"]),
                                     "{:.12e}".format(rec["z0"]), "{:.12e}".format(rec["u1"]),
                                     "{:.12e}".format(rec["u2"]), "{:.12e}".format(rec["u3"]),
                                     "{:.12e}".format(rec["x_def"]), "{:.12e}".format(rec["y_def"]),
                                     "{:.12e}".format(rec["z_def"]), "{:.12e}".format(rec["reference_z"]),
                                     "{:.12e}".format(rec["relative_height"]), "{:.12e}".format(rec["depth"])])

        print("\n".join(summary_lines))
        print("\n结果文件已写入当前工作目录：")
        print("  {}".format(SUMMARY_FILE_NAME))
        if WRITE_CSV_FILES:
            print("  {}".format(X_PROFILE_FILE_NAME))
            print("  {}".format(Y_PROFILE_FILE_NAME))
            print("  {}".format(TOP_SURFACE_FILE_NAME))
        else:
            print("  WRITE_CSV_FILES=False：已跳过全部 CSV 文件。")

    finally:
        odb.close()


if __name__ == "__main__":
    main()
