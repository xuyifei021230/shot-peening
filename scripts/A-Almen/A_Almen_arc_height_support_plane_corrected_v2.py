#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A-Almen 试片弧高自动测量脚本——四支撑点基准平面法（修正版）

修正原理
--------
原脚本把“中心点到 X/Y 最外边缘的 U3 差”分别作为两个方向弧高，再相加。
这种算法不等同于 Almen 弧高仪的实际测量原理。

Almen 弧高仪的正确几何关系为：

    1. 试片由四个精密支撑球支承；
    2. 四个支撑接触点 A、B、C、D 确定测量基准平面；
    3. 指示表测头位于支撑矩形中心；
    4. 弧高是中心测量点相对于该基准平面的高度。

本脚本采用以下数值实现：

    1. 在试片未喷丸侧（默认下表面）提取位移；
    2. 在标准四支撑点的名义位置插值得到变形后坐标；
    3. 对四个支撑点进行最小二乘平面拟合：
           z = a*x + b*y + c
    4. 在中心测量点的 X、Y 坐标处求基准平面的 Z 高度；
    5. 计算中心点沿全局 Z 方向到基准平面的高度：
           h_Z = z_center - z_plane
    6. 同时输出中心点到平面的法向距离，作为倾斜检查值；
    7. 可选地用初始帧进行预弯补偿：
           h_corrected = h_final - h_initial

标准名义尺寸
------------
四个支撑球中心距默认采用：

    长度方向：31.750 mm
    宽度方向：15.875 mm

因此，相对于完整试片中心，四个名义支撑位置为：

    (+15.8750, +7.9375)
    (+15.8750, -7.9375)
    (-15.8750, +7.9375)
    (-15.8750, -7.9375)

当前用户模型是完整试片的 1/4：

    X = 0 ~ 38.0 mm
    Y = 0 ~ 9.5 mm
    x = 0、y = 0 为对称面

脚本会根据四分之一对称关系重建完整试片的四个支撑点。

重要说明
--------
1. 默认测量未喷丸侧，即 BOTTOM；当前建模脚本的喷丸影响层位于上表面。
2. 弧高使用 U3/变形后 Z 坐标，不使用 U, Magnitude。
3. 支撑点不一定正好落在网格节点上，因此采用局部二次曲面拟合插值。
4. 四个数值支撑点若存在微小不共面误差，使用最小二乘平面拟合。
5. 这是合理的后处理等效方法。最严格的物理复现仍需在有限元模型中
   建立四个刚性支撑球、压紧作用和中心测头接触。
6. Abaqus 使用一致单位制；本脚本假定几何和位移单位均为 mm。
7. 已显式调用 Python 内置的 sum/max/min，兼容 Abaqus/CAE Run Script
   环境中同名函数被覆盖的情况。

运行方法
--------
在 Abaqus/CAE 中：

    File -> Run Script -> 选择本脚本

或在命令行中：

    abaqus python A_Almen_arc_height_support_plane_corrected.py
"""

from __future__ import print_function

import os
import io
import math

# Abaqus/CAE 的脚本全局环境可能包含与 Python 内置函数同名的对象，
# 例如 sum、max、min。显式从 Python 内置模块取得这些函数，
# 避免出现 “found generator, expected a recognized type”。
try:
    import builtins as _python_builtins
except ImportError:
    import __builtin__ as _python_builtins

_py_sum = _python_builtins.sum
_py_max = _python_builtins.max
_py_min = _python_builtins.min

from odbAccess import openOdb


# =============================================================================
# 0. 用户参数区
# =============================================================================

# ---------- ODB ----------
ODB_FILE = 'Job-1.odb'
INSTANCE_NAME = 'A-Almen-1'
STEP_NAME = 'Step-1'

# 最终结果帧，-1 表示最后一帧。
FINAL_FRAME_INDEX = -1

# 是否用分析步第 0 帧进行预弯/初始形状补偿。
APPLY_INITIAL_FRAME_CORRECTION = True
INITIAL_FRAME_INDEX = 0

# ---------- 模型类型 ----------
# 当前模型为四分之一双对称模型。
# 可选：
#   'QUARTER_SYMMETRY'：X>=0、Y>=0，x=0/y=0 为对称面；
#   'FULL_MODEL'      ：完整试片，坐标中心在 (FULL_CENTER_X,FULL_CENTER_Y)。
MODEL_TYPE = 'QUARTER_SYMMETRY'

# 完整模型使用时的几何中心坐标。
FULL_CENTER_X = 0.0
FULL_CENTER_Y = 0.0

# ---------- 测量表面 ----------
# 当前模型的喷丸影响层在上表面，Almen 弧高仪应从未喷丸侧确定支撑平面，
# 因此默认选 BOTTOM。
MEASUREMENT_SURFACE = 'BOTTOM'

# ---------- Almen 弧高仪名义几何 ----------
# 四个支撑球的中心距/名义接触点矩形尺寸。
SUPPORT_SPACING_X = 31.750
SUPPORT_SPACING_Y = 15.875

# 支撑球直径。当前后处理算法不直接使用，仅写入报告便于追溯。
SUPPORT_BALL_DIAMETER = 4.760

# ---------- 插值参数 ----------
# 支撑点通常不与网格节点重合，使用邻近表面节点进行局部二次拟合。
LOCAL_FIT_NODE_NUMBER = 24

# 局部拟合最大搜索半径，单位 mm。
# 设为 None 时只按最近节点数选择，不限制半径。
LOCAL_FIT_MAX_RADIUS = 4.0

# 当目标点与某节点的平面距离小于此值时，直接使用该节点值。
EXACT_NODE_TOL = 1.0e-8

# 表面节点 Z 坐标筛选容差。
SURFACE_Z_TOL = 1.0e-5

# 二次拟合失效时，反距离加权的幂指数。
IDW_POWER = 2.0

# ---------- 输出 ----------
SUMMARY_FILE_NAME = 'A_Almen-02MPa.txt'
SUPPORT_POINTS_FILE_NAME = 'A_Almen-02MPa.csv'

# 是否把报告同时打印到 Abaqus 信息区。
PRINT_FULL_REPORT = True


# =============================================================================
# 1. 基础函数
# =============================================================================

def info(message):
    """在 Abaqus 信息区输出进度。"""
    print('[信息] ' + str(message))


def warn(message):
    """在 Abaqus 信息区输出警告。"""
    print('[警告] ' + str(message))


def resolve_odb_path(file_name):
    """
    查找 ODB，顺序为：
        1. 绝对路径；
        2. 脚本所在目录；
        3. 当前工作目录；
        4. 当前 CAE 所在目录。
    """
    candidates = []

    if os.path.isabs(file_name):
        candidates.append(file_name)
    else:
        if '__file__' in globals():
            script_dir = os.path.dirname(os.path.abspath(__file__))
            candidates.append(os.path.join(script_dir, file_name))

        candidates.append(os.path.join(os.getcwd(), file_name))

        try:
            from abaqus import mdb
            if mdb.pathName:
                cae_dir = os.path.dirname(os.path.abspath(mdb.pathName))
                candidates.append(os.path.join(cae_dir, file_name))
        except Exception:
            pass

    unique_paths = []
    for path in candidates:
        path = os.path.abspath(path)
        if path not in unique_paths:
            unique_paths.append(path)

    for path in unique_paths:
        if os.path.isfile(path):
            return path

    message = ['找不到 ODB 文件：%s' % file_name, '已搜索以下位置：']
    message.extend(['  - ' + path for path in unique_paths])
    raise RuntimeError('\n'.join(message))


def find_key_ignore_case(repository, requested_name, description):
    """忽略大小写查找 ODB repository 键。"""
    keys = list(repository.keys())

    if requested_name in keys:
        return requested_name

    requested_upper = requested_name.upper()
    for key in keys:
        if key.upper() == requested_upper:
            return key

    raise RuntimeError(
        '找不到%s：%s\n可用名称：%s' %
        (description, requested_name, ', '.join(keys))
    )


def normalize_frame_index(frame_count, requested_index, description):
    """把负帧编号转换为实际编号，并检查范围。"""
    index = requested_index

    if index < 0:
        index = frame_count + index

    if index < 0 or index >= frame_count:
        raise RuntimeError(
            '%s帧编号无效：%d；该分析步共有 %d 帧。' %
            (description, requested_index, frame_count)
        )

    return index


def get_step_and_frames(odb):
    """获取分析步、初始帧和最终帧。"""
    step_keys = list(odb.steps.keys())
    if not step_keys:
        raise RuntimeError('ODB 中没有分析步。')

    if STEP_NAME:
        step_key = find_key_ignore_case(odb.steps, STEP_NAME, '分析步')
    else:
        step_key = step_keys[-1]

    step = odb.steps[step_key]
    frame_count = len(step.frames)

    if frame_count == 0:
        raise RuntimeError('分析步 %s 中没有结果帧。' % step_key)

    final_index = normalize_frame_index(
        frame_count, FINAL_FRAME_INDEX, '最终'
    )
    final_frame = step.frames[final_index]

    initial_index = None
    initial_frame = None

    if APPLY_INITIAL_FRAME_CORRECTION:
        initial_index = normalize_frame_index(
            frame_count, INITIAL_FRAME_INDEX, '初始'
        )
        initial_frame = step.frames[initial_index]

    return (
        step_key,
        final_index,
        final_frame,
        initial_index,
        initial_frame
    )


def build_displacement_dictionary(frame, instance):
    """
    建立 nodeLabel -> (U1,U2,U3) 字典。

    若第 0 帧没有 U，则按零位移处理。
    """
    if 'U' not in frame.fieldOutputs.keys():
        warn('当前帧没有 U 输出，按零位移处理。')
        return dict(
            (node.label, (0.0, 0.0, 0.0))
            for node in instance.nodes
        )

    subset = frame.fieldOutputs['U'].getSubset(region=instance)
    displacement = {}

    for value in subset.values:
        data = value.data

        if len(data) < 3:
            raise RuntimeError('本脚本要求三维位移 U1、U2、U3。')

        displacement[value.nodeLabel] = (
            float(data[0]),
            float(data[1]),
            float(data[2])
        )

    if not displacement:
        raise RuntimeError('没有读取到实例的节点位移。')

    return displacement


def get_surface_nodes(instance):
    """
    按初始坐标选择上表面或下表面节点。

    对当前平板模型：
        TOP    -> 最大初始 Z；
        BOTTOM -> 最小初始 Z。
    """
    nodes = list(instance.nodes)

    if not nodes:
        raise RuntimeError('所选实例没有节点。')

    z_values = [float(node.coordinates[2]) for node in nodes]
    z_min = _py_min(z_values)
    z_max = _py_max(z_values)

    surface_mode = MEASUREMENT_SURFACE.upper()

    if surface_mode == 'TOP':
        target_z = z_max
        surface_cn = '上表面'
    elif surface_mode == 'BOTTOM':
        target_z = z_min
        surface_cn = '下表面（未喷丸侧）'
    else:
        raise RuntimeError(
            "MEASUREMENT_SURFACE 只能为 'TOP' 或 'BOTTOM'。"
        )

    surface_nodes = [
        node for node in nodes
        if abs(float(node.coordinates[2]) - target_z) <= SURFACE_Z_TOL
    ]

    if not surface_nodes:
        raise RuntimeError(
            '没有选到%s节点；目标 Z=%.9f，容差=%.3e。' %
            (surface_cn, target_z, SURFACE_Z_TOL)
        )

    info(
        '测量表面：%s；初始 Z=%.9f mm；节点数=%d。' %
        (surface_cn, target_z, len(surface_nodes))
    )

    return surface_nodes, target_z, surface_cn


# =============================================================================
# 2. 四分之一模型对称重建
# =============================================================================

def append_unique_sample(sample_list, sample_keys, sample):
    """按初始 XY 坐标去重添加采样点。"""
    key = (
        round(sample['x0'], 10),
        round(sample['y0'], 10)
    )

    if key not in sample_keys:
        sample_keys[key] = True
        sample_list.append(sample)


def build_full_surface_samples(surface_nodes, displacement, surface_z):
    """
    建立完整试片表面的采样点。

    QUARTER_SYMMETRY：
        将第一象限节点按 x=0、y=0 对称重建到四个象限。
        对称位移关系：
            U1(-x,y) = -U1(x,y)
            U2(x,-y) = -U2(x,y)
            U3 在镜像点保持相同。

    FULL_MODEL：
        直接使用完整模型节点，并把坐标平移到完整试片中心坐标系。
    """
    model_type = MODEL_TYPE.upper()
    samples = []
    sample_keys = {}

    for node in surface_nodes:
        x_global = float(node.coordinates[0])
        y_global = float(node.coordinates[1])
        z0 = float(node.coordinates[2])

        if node.label not in displacement:
            raise RuntimeError(
                '节点 %d 没有位移 U。' % node.label
            )

        u1, u2, u3 = displacement[node.label]

        if model_type == 'QUARTER_SYMMETRY':
            x_local = x_global
            y_local = y_global

            sign_x_values = (1.0,) if abs(x_local) <= EXACT_NODE_TOL else (1.0, -1.0)
            sign_y_values = (1.0,) if abs(y_local) <= EXACT_NODE_TOL else (1.0, -1.0)

            for sign_x in sign_x_values:
                for sign_y in sign_y_values:
                    sample = {
                        'source_label': node.label,
                        'x0': sign_x * x_local,
                        'y0': sign_y * y_local,
                        'z0': z0,
                        'u1': sign_x * u1,
                        'u2': sign_y * u2,
                        'u3': u3
                    }
                    append_unique_sample(samples, sample_keys, sample)

        elif model_type == 'FULL_MODEL':
            sample = {
                'source_label': node.label,
                'x0': x_global - FULL_CENTER_X,
                'y0': y_global - FULL_CENTER_Y,
                'z0': z0,
                'u1': u1,
                'u2': u2,
                'u3': u3
            }
            append_unique_sample(samples, sample_keys, sample)

        else:
            raise RuntimeError(
                "MODEL_TYPE 只能为 'QUARTER_SYMMETRY' 或 'FULL_MODEL'。"
            )

    if len(samples) < 6:
        raise RuntimeError('完整表面采样点少于 6 个，无法进行二次拟合。')

    return samples


# =============================================================================
# 3. 线性方程与局部插值
# =============================================================================

def solve_linear_system(matrix, vector, pivot_tolerance=1.0e-14):
    """
    使用带部分选主元的高斯消元求解线性方程组。

    返回：
        解向量；
    若矩阵奇异：
        返回 None。
    """
    n = len(vector)

    augmented = []
    for i in range(n):
        augmented.append(
            [float(value) for value in matrix[i]] + [float(vector[i])]
        )

    for column in range(n):
        pivot_row = _py_max(
            list(range(column, n)),
            key=lambda row: abs(augmented[row][column])
        )

        pivot_value = augmented[pivot_row][column]

        if abs(pivot_value) <= pivot_tolerance:
            return None

        if pivot_row != column:
            augmented[column], augmented[pivot_row] = (
                augmented[pivot_row],
                augmented[column]
            )

        pivot_value = augmented[column][column]

        for j in range(column, n + 1):
            augmented[column][j] /= pivot_value

        for row in range(n):
            if row == column:
                continue

            factor = augmented[row][column]

            if factor == 0.0:
                continue

            for j in range(column, n + 1):
                augmented[row][j] -= factor * augmented[column][j]

    return [augmented[i][n] for i in range(n)]


def idw_interpolate(selected_samples, target_x, target_y, field_name):
    """
    反距离加权插值，作为二次拟合失败时的备用算法。
    """
    numerator = 0.0
    denominator = 0.0

    for sample in selected_samples:
        dx = sample['x0'] - target_x
        dy = sample['y0'] - target_y
        distance = math.sqrt(dx * dx + dy * dy)

        if distance <= EXACT_NODE_TOL:
            return sample[field_name]

        weight = 1.0 / (distance ** IDW_POWER)
        numerator += weight * sample[field_name]
        denominator += weight

    if denominator == 0.0:
        raise RuntimeError('反距离加权插值失败：权重和为零。')

    return numerator / denominator


def local_quadratic_interpolate(samples, target_x, target_y, field_name):
    """
    在目标点附近拟合局部二次曲面：

        f = c0 + c1*X + c2*Y + c3*X^2 + c4*X*Y + c5*Y^2

    其中 X、Y 是以目标点为原点并经过尺度归一化的局部坐标。
    目标点的插值值就是 c0。

    返回：
        value             插值值
        method            使用的方法
        nearest_distance  最近采样点距离
        used_count        拟合使用节点数
        farthest_distance 最远使用节点距离
    """
    distances = []

    for sample in samples:
        dx = sample['x0'] - target_x
        dy = sample['y0'] - target_y
        d2 = dx * dx + dy * dy
        distances.append((d2, sample))

    distances.sort(key=lambda item: item[0])

    if not distances:
        raise RuntimeError('没有可用于插值的表面采样点。')

    nearest_distance = math.sqrt(distances[0][0])

    if nearest_distance <= EXACT_NODE_TOL:
        return (
            distances[0][1][field_name],
            '精确节点',
            nearest_distance,
            1,
            nearest_distance
        )

    selected = []

    for d2, sample in distances:
        distance = math.sqrt(d2)

        if LOCAL_FIT_MAX_RADIUS is not None and distance > LOCAL_FIT_MAX_RADIUS:
            continue

        selected.append((distance, sample))

        if len(selected) >= LOCAL_FIT_NODE_NUMBER:
            break

    # 半径限制导致点数不足时，退回最近的若干点。
    if len(selected) < 6:
        selected = [
            (math.sqrt(d2), sample)
            for d2, sample in distances[:_py_max(6, LOCAL_FIT_NODE_NUMBER)]
        ]

    used_count = len(selected)
    farthest_distance = _py_max([item[0] for item in selected])

    # 使用局部距离尺度改善正规方程的数值条件。
    scale = _py_max(farthest_distance, 1.0e-6)

    size = 6
    normal_matrix = [[0.0 for _ in range(size)] for _ in range(size)]
    right_vector = [0.0 for _ in range(size)]

    for distance, sample in selected:
        x_local = (sample['x0'] - target_x) / scale
        y_local = (sample['y0'] - target_y) / scale

        basis = [
            1.0,
            x_local,
            y_local,
            x_local * x_local,
            x_local * y_local,
            y_local * y_local
        ]

        # 距离越近权重越大；加入小量避免极端权重。
        weight = 1.0 / (
            (distance / scale) * (distance / scale) + 1.0e-6
        )

        value = sample[field_name]

        for i in range(size):
            right_vector[i] += weight * basis[i] * value

            for j in range(size):
                normal_matrix[i][j] += (
                    weight * basis[i] * basis[j]
                )

    # 极小正则化，仅用于抑制病态矩阵。
    trace = _py_sum([normal_matrix[i][i] for i in range(size)])
    regularization = _py_max(trace, 1.0) * 1.0e-12

    for i in range(size):
        normal_matrix[i][i] += regularization

    coefficients = solve_linear_system(normal_matrix, right_vector)

    if coefficients is None:
        value = idw_interpolate(
            [item[1] for item in selected],
            target_x,
            target_y,
            field_name
        )
        method = '反距离加权备用法'
    else:
        value = coefficients[0]
        method = '局部二次曲面拟合'

    return (
        value,
        method,
        nearest_distance,
        used_count,
        farthest_distance
    )


def interpolate_material_point(samples, target_x, target_y, surface_z, point_name):
    """
    在名义材料坐标 (target_x,target_y) 处插值 U1、U2、U3，
    并计算变形后坐标。
    """
    result = {
        'name': point_name,
        'x0': target_x,
        'y0': target_y,
        'z0': surface_z
    }

    metadata = None

    for field_name in ('u1', 'u2', 'u3'):
        (
            value,
            method,
            nearest_distance,
            used_count,
            farthest_distance
        ) = local_quadratic_interpolate(
            samples,
            target_x,
            target_y,
            field_name
        )

        result[field_name] = value

        if metadata is None:
            metadata = {
                'method': method,
                'nearest_distance': nearest_distance,
                'used_count': used_count,
                'farthest_distance': farthest_distance
            }

    result.update(metadata)

    result['xd'] = result['x0'] + result['u1']
    result['yd'] = result['y0'] + result['u2']
    result['zd'] = result['z0'] + result['u3']

    return result


# =============================================================================
# 4. 四点平面拟合与弧高
# =============================================================================

def fit_plane_z_xy(points):
    """
    用最小二乘拟合四个支撑点的平面：

        z = a*x + b*y + c

    返回：
        a, b, c, 各支撑点到拟合平面的 Z 残差、最大绝对残差。
    """
    n = float(len(points))

    sum_x = _py_sum([point['xd'] for point in points])
    sum_y = _py_sum([point['yd'] for point in points])
    sum_z = _py_sum([point['zd'] for point in points])

    sum_xx = _py_sum([point['xd'] * point['xd'] for point in points])
    sum_xy = _py_sum([point['xd'] * point['yd'] for point in points])
    sum_yy = _py_sum([point['yd'] * point['yd'] for point in points])

    sum_xz = _py_sum([point['xd'] * point['zd'] for point in points])
    sum_yz = _py_sum([point['yd'] * point['zd'] for point in points])

    matrix = [
        [sum_xx, sum_xy, sum_x],
        [sum_xy, sum_yy, sum_y],
        [sum_x,  sum_y,  n]
    ]

    vector = [sum_xz, sum_yz, sum_z]

    solution = solve_linear_system(matrix, vector)

    if solution is None:
        raise RuntimeError('四支撑点基准平面拟合失败。')

    a, b, c = solution

    residuals = []
    for point in points:
        z_plane = a * point['xd'] + b * point['yd'] + c
        residuals.append(point['zd'] - z_plane)

    max_abs_residual = _py_max([abs(value) for value in residuals])

    return a, b, c, residuals, max_abs_residual


def calculate_frame_arc_height(
    frame,
    instance,
    surface_nodes,
    surface_z,
    frame_description
):
    """
    对一个结果帧计算四支撑平面法弧高。
    """
    displacement = build_displacement_dictionary(frame, instance)

    samples = build_full_surface_samples(
        surface_nodes,
        displacement,
        surface_z
    )

    half_x = SUPPORT_SPACING_X / 2.0
    half_y = SUPPORT_SPACING_Y / 2.0

    nominal_supports = [
        ('支撑点A（+X,+Y）', +half_x, +half_y),
        ('支撑点B（-X,+Y）', -half_x, +half_y),
        ('支撑点C（-X,-Y）', -half_x, -half_y),
        ('支撑点D（+X,-Y）', +half_x, -half_y)
    ]

    support_points = []

    for name, x0, y0 in nominal_supports:
        point = interpolate_material_point(
            samples,
            x0,
            y0,
            surface_z,
            name
        )
        support_points.append(point)

    center = interpolate_material_point(
        samples,
        0.0,
        0.0,
        surface_z,
        '中心测量点'
    )

    a, b, c, residuals, max_residual = fit_plane_z_xy(
        support_points
    )

    # 沿全局 Z 方向测量：在中心点的变形后 X/Y 位置上计算基准面 Z。
    z_plane_at_center = (
        a * center['xd'] +
        b * center['yd'] +
        c
    )

    height_z_signed = center['zd'] - z_plane_at_center

    # 中心点到平面的法向距离。
    normal_denominator = math.sqrt(a * a + b * b + 1.0)
    height_normal_signed = height_z_signed / normal_denominator

    for index, point in enumerate(support_points):
        point['plane_residual_z'] = residuals[index]

    return {
        'frame_description': frame_description,
        'frame_value': float(frame.frameValue),
        'support_points': support_points,
        'center': center,
        'plane_a': a,
        'plane_b': b,
        'plane_c': c,
        'plane_max_residual_z': max_residual,
        'z_plane_at_center': z_plane_at_center,
        'height_z_signed': height_z_signed,
        'height_z_abs': abs(height_z_signed),
        'height_normal_signed': height_normal_signed,
        'height_normal_abs': abs(height_normal_signed)
    }


# =============================================================================
# 5. 输出
# =============================================================================

def format_point_line(point):
    """把关键点结果格式化为中文文本。"""
    return (
        '%-20s 初始坐标=(% .6f,% .6f,% .6f) mm；'
        '位移=(% .9e,% .9e,% .9e) mm；'
        '变形后坐标=(% .9f,% .9f,% .9f) mm'
    ) % (
        point['name'],
        point['x0'], point['y0'], point['z0'],
        point['u1'], point['u2'], point['u3'],
        point['xd'], point['yd'], point['zd']
    )


def write_points_csv(path, initial_result, final_result):
    """输出初始帧和最终帧的中心点、四支撑点及平面残差。"""
    with io.open(path, 'w', encoding='utf-8-sig', newline='') as file_object:
        header = [
            '结果帧',
            '点名称',
            '初始X_mm',
            '初始Y_mm',
            '初始Z_mm',
            'U1_mm',
            'U2_mm',
            'U3_mm',
            '变形后X_mm',
            '变形后Y_mm',
            '变形后Z_mm',
            '支撑点平面Z残差_mm',
            '插值方法',
            '最近节点距离_mm',
            '拟合节点数',
            '拟合最远节点距离_mm'
        ]
        file_object.write(','.join(header) + '\n')

        result_list = []
        if initial_result is not None:
            result_list.append(initial_result)
        result_list.append(final_result)

        for result in result_list:
            points = list(result['support_points']) + [result['center']]

            for point in points:
                residual = point.get('plane_residual_z', '')

                row = [
                    result['frame_description'],
                    point['name'],
                    '%.12g' % point['x0'],
                    '%.12g' % point['y0'],
                    '%.12g' % point['z0'],
                    '%.12g' % point['u1'],
                    '%.12g' % point['u2'],
                    '%.12g' % point['u3'],
                    '%.12g' % point['xd'],
                    '%.12g' % point['yd'],
                    '%.12g' % point['zd'],
                    '' if residual == '' else '%.12g' % residual,
                    point['method'],
                    '%.12g' % point['nearest_distance'],
                    str(point['used_count']),
                    '%.12g' % point['farthest_distance']
                ]

                file_object.write(','.join(row) + '\n')


def build_report(
    odb_path,
    instance_key,
    step_key,
    initial_index,
    final_index,
    surface_cn,
    initial_result,
    final_result,
    corrected_height_z_signed,
    corrected_height_normal_signed
):
    """生成中文 TXT 报告。"""
    lines = []

    lines.append('=' * 86)
    lines.append('A-Almen 试片弧高测量报告——四支撑点基准平面法')
    lines.append('=' * 86)
    lines.append('ODB结果文件              ：%s' % odb_path)
    lines.append('分析实例                  ：%s' % instance_key)
    lines.append('分析步                    ：%s' % step_key)
    lines.append('模型类型                  ：%s' % MODEL_TYPE)
    lines.append('测量表面                  ：%s' % surface_cn)
    lines.append('最终帧编号                ：%d' % final_index)

    if initial_result is not None:
        lines.append('初始补偿帧编号            ：%d' % initial_index)
    else:
        lines.append('初始帧补偿                ：未启用')

    lines.append('')
    lines.append('一、弧高仪名义几何')
    lines.append('-' * 86)
    lines.append(
        '四支撑点长度方向中心距    ：%.6f mm' % SUPPORT_SPACING_X
    )
    lines.append(
        '四支撑点宽度方向中心距    ：%.6f mm' % SUPPORT_SPACING_Y
    )
    lines.append(
        '支撑球名义直径            ：%.6f mm' % SUPPORT_BALL_DIAMETER
    )
    lines.append(
        '中心测量点名义坐标        ：(0.000000, 0.000000)'
    )
    lines.append(
        '四支撑点名义半间距        ：X=%.6f mm，Y=%.6f mm' %
        (SUPPORT_SPACING_X / 2.0, SUPPORT_SPACING_Y / 2.0)
    )

    if initial_result is not None:
        lines.append('')
        lines.append('二、初始帧基准结果')
        lines.append('-' * 86)
        lines.append(
            '初始帧值/时间             ：%.12g' %
            initial_result['frame_value']
        )
        lines.append(
            '初始基准平面              ：z = (% .9e)x + (% .9e)y + (% .9e)' %
            (
                initial_result['plane_a'],
                initial_result['plane_b'],
                initial_result['plane_c']
            )
        )
        lines.append(
            '初始四点最大不共面残差    ：%.9e mm' %
            initial_result['plane_max_residual_z']
        )
        lines.append(
            '初始沿Z方向有符号弧高     ：% .9e mm' %
            initial_result['height_z_signed']
        )
        lines.append(
            '初始沿平面法向有符号距离  ：% .9e mm' %
            initial_result['height_normal_signed']
        )

    lines.append('')
    lines.append('三、最终帧四支撑点与中心点')
    lines.append('-' * 86)

    for point in final_result['support_points']:
        lines.append(format_point_line(point))
        lines.append(
            '    相对拟合平面的Z残差：% .9e mm；插值：%s；最近节点距离：%.6f mm' %
            (
                point['plane_residual_z'],
                point['method'],
                point['nearest_distance']
            )
        )

    lines.append(format_point_line(final_result['center']))
    lines.append(
        '    插值：%s；最近节点距离：%.6f mm' %
        (
            final_result['center']['method'],
            final_result['center']['nearest_distance']
        )
    )

    lines.append('')
    lines.append('四、最终帧基准平面')
    lines.append('-' * 86)
    lines.append(
        '基准平面方程              ：z = (% .9e)x + (% .9e)y + (% .9e)' %
        (
            final_result['plane_a'],
            final_result['plane_b'],
            final_result['plane_c']
        )
    )
    lines.append(
        '四支撑点最大不共面残差    ：%.9e mm' %
        final_result['plane_max_residual_z']
    )
    lines.append(
        '中心点处基准平面Z高度     ：%.9f mm' %
        final_result['z_plane_at_center']
    )
    lines.append(
        '中心点变形后Z坐标         ：%.9f mm' %
        final_result['center']['zd']
    )

    lines.append('')
    lines.append('五、弧高结果')
    lines.append('-' * 86)
    lines.append(
        '最终沿全局Z方向有符号弧高 ：% .9f mm' %
        final_result['height_z_signed']
    )
    lines.append(
        '最终沿全局Z方向弧高绝对值 ：%.9f mm' %
        final_result['height_z_abs']
    )
    lines.append(
        '最终沿平面法向距离绝对值  ：%.9f mm' %
        final_result['height_normal_abs']
    )

    if initial_result is not None:
        lines.append(
            '预弯补偿公式              ：h修正 = h最终 - h初始'
        )
        lines.append(
            '补偿后沿Z方向有符号弧高   ：% .9f mm' %
            corrected_height_z_signed
        )
        lines.append(
            '补偿后沿Z方向弧高绝对值   ：%.9f mm' %
            abs(corrected_height_z_signed)
        )
        lines.append(
            '补偿后沿平面法向距离绝对值：%.9f mm' %
            abs(corrected_height_normal_signed)
        )
        primary_height = abs(corrected_height_z_signed)
    else:
        primary_height = final_result['height_z_abs']

    lines.append('')
    lines.append('最终采用的 Almen 弧高     ：%.9f mm' % primary_height)

    lines.append('')
    lines.append('六、结果判读')
    lines.append('-' * 86)
    lines.append(
        '1. 本结果不是中心点到试片最外边缘的高度差，也不再把 X/Y 两方向结果相加。'
    )
    lines.append(
        '2. 本结果先由四个标准支撑位置拟合基准平面，再测量中心点相对该平面的高度。'
    )
    lines.append(
        '3. “沿全局Z方向弧高”对应当前 Abaqus 坐标系中的厚度方向测量。'
    )
    lines.append(
        '4. 若基准平面存在明显倾斜，应同时检查“平面法向距离”和模型刚体转动。'
    )
    lines.append(
        '5. 四点最大不共面残差应远小于目标弧高；残差过大说明非对称、插值或网格需要检查。'
    )
    lines.append(
        '6. 若要与真实 Almen 弧高仪完全一致，应进一步建立支撑球、压紧力和测头接触模型。'
    )
    lines.append('=' * 86)

    return '\n'.join(lines), primary_height


# =============================================================================
# 6. 主程序
# =============================================================================

def main():
    """执行四支撑点基准平面法弧高测量。"""
    info('=' * 86)
    info('开始执行 A-Almen 四支撑点基准平面法弧高测量。')

    if SUPPORT_SPACING_X <= 0.0 or SUPPORT_SPACING_Y <= 0.0:
        raise RuntimeError('支撑点中心距必须大于零。')

    odb_path = resolve_odb_path(ODB_FILE)
    output_directory = os.path.dirname(odb_path)

    info('ODB：' + odb_path)

    odb = openOdb(path=odb_path, readOnly=True)

    try:
        instance_key = find_key_ignore_case(
            odb.rootAssembly.instances,
            INSTANCE_NAME,
            '实例'
        )
        instance = odb.rootAssembly.instances[instance_key]

        (
            step_key,
            final_index,
            final_frame,
            initial_index,
            initial_frame
        ) = get_step_and_frames(odb)

        surface_nodes, surface_z, surface_cn = get_surface_nodes(
            instance
        )

        initial_result = None

        if initial_frame is not None:
            info(
                '计算初始补偿帧：编号=%d，帧值/时间=%.12g。' %
                (initial_index, float(initial_frame.frameValue))
            )
            initial_result = calculate_frame_arc_height(
                initial_frame,
                instance,
                surface_nodes,
                surface_z,
                '初始帧'
            )

        info(
            '计算最终帧：编号=%d，帧值/时间=%.12g。' %
            (final_index, float(final_frame.frameValue))
        )
        final_result = calculate_frame_arc_height(
            final_frame,
            instance,
            surface_nodes,
            surface_z,
            '最终帧'
        )

        if initial_result is not None:
            corrected_height_z_signed = (
                final_result['height_z_signed'] -
                initial_result['height_z_signed']
            )
            corrected_height_normal_signed = (
                final_result['height_normal_signed'] -
                initial_result['height_normal_signed']
            )
        else:
            corrected_height_z_signed = final_result['height_z_signed']
            corrected_height_normal_signed = (
                final_result['height_normal_signed']
            )

        report, primary_height = build_report(
            odb_path=odb_path,
            instance_key=instance_key,
            step_key=step_key,
            initial_index=initial_index,
            final_index=final_index,
            surface_cn=surface_cn,
            initial_result=initial_result,
            final_result=final_result,
            corrected_height_z_signed=corrected_height_z_signed,
            corrected_height_normal_signed=corrected_height_normal_signed
        )

        summary_path = os.path.join(
            output_directory,
            SUMMARY_FILE_NAME
        )
        points_path = os.path.join(
            output_directory,
            SUPPORT_POINTS_FILE_NAME
        )

        with io.open(
            summary_path,
            'w',
            encoding='utf-8-sig',
            newline=''
        ) as file_object:
            file_object.write(report + '\n')

        write_points_csv(
            points_path,
            initial_result,
            final_result
        )

        if PRINT_FULL_REPORT:
            print(report)

        info('最终采用弧高：%.9f mm' % primary_height)
        info('中文汇总报告：' + summary_path)
        info('支撑点与中心点数据：' + points_path)
        info('弧高测量完成。')
        info('=' * 86)

    finally:
        odb.close()


if __name__ == '__main__':
    main()
