# -*- coding: utf-8 -*-
"""
Abaqus/CAE 随机多弹丸喷丸全流程建模脚本
============================================================
适用目标：复刻“随机喷丸-Python二次开发”PPT/视频中的建模流程，并将
原教程中需要手工完成的材料、网格、装配、分析步、接触、载荷、结果分层集、
作业创建等步骤尽量自动化。

建议运行方式：
    Abaqus/CAE -> 文件(File) -> 运行脚本(Run Script) -> 选择本文件

单位制（本脚本默认）：
    长度 mm；时间 s；质量 tonne；应力 MPa；密度 tonne/mm^3；速度 mm/s。

兼容性设计：
    1. 不使用 f-string、类型注解等新语法，兼容 Abaqus 2022 的 Python 2.7
       以及 Abaqus 2024/2025 的 Python 3。
    2. 脚本未在当前环境中调用 Abaqus 求解器实机运行；若不同版本的 API
       对某一网格控制参数存在差异，脚本会尽量给出明确的报错位置。

修复记录：
    2026-07-16：修复装配阶段创建 SURF-SHOT-ALL 时的“特征创建失败”。
    修复方式为：只筛选球体外表面，在 Part 层级创建表面，并在装配层级
    使用 SurfaceByBoolean 合并；若当前 Abaqus 版本不支持，则自动切换为
    General Contact 的 All with self。

    2026-08-03（V2）：
    1. 默认改用 General Contact 的 All with self，避免人工表面法向不一致造成
       improperly defined surface(s) 和 Explicit Packager 中止。
    2. 靶材固定区改为按网格节点坐标筛选，并用 SetFromNodeLabels 创建装配节点集，
       避免 assembly.Set(faces=...) 出现“特征创建失败”。
    3. 节点集创建失败时自动退回临时 Region，保证边界条件继续建立。

    2026-08-03（V3）：
    1. 将弹丸球心的 X、Y 随机范围与靶材中心加密区域直接关联。
    2. 在 X、Y 边界预留一个弹丸半径，确保弹丸沿 Z 方向下落后，
       整个弹丸投影而不仅是球心均落在中心区域内。
    3. 增加中心区域尺寸检查和生成后的逐颗弹丸越界复核。

    2026-08-03（V4）：
    1. 修复 MODEL_NAME='Model-1' 且它是数据库中唯一模型时，直接删除
       mdb.models['Model-1'] 触发 KeyError 的问题。
    2. 删除同名模型时先创建临时替代模型，再删除旧模型并重命名临时模型，
       保证模型数据库中始终至少存在一个模型。
    3. 模型存在性检查改为直接访问 Repository，避免 keys() 状态异常。

    2026-08-04（读取记录修正）：
    1. 默认创建全部弹丸外表面的装配级组合表面，并仅指定弹丸自接触及
       弹丸/靶材顶面接触，不再默认使用 All with self。
    2. 场输出改为 200 个时间间隔，并按节点位置输出。
    3. 结果分层单元集的 X、Y 范围直接跟随中心加密区域。

    2026-08-04（相互作用记录修正）：
    1. 按 Abaqus/CAE 2025 成功录制的 FaceArray 拼接方式创建 Surf-shot，
       避免把跨实例面转换成普通 tuple 后出现“特征创建失败”。
    2. 使用录制成功的拓扑掩码创建 Surf-shot 和 Surf-target。
    3. General Contact 仅包含 Surf-shot/Surf-target 接触对，并增加记录中的
       wearSurfacePropertyAssignments 设置。

    2026-08-05（V8）：
    1. 结果分层集合改用参考脚本末尾的独立包围盒范围：X、Y 均为
       [-0.4, 0.4] mm，并在各方向加入 1.0e-4 mm 容差。
    2. 集合范围不再跟随中心加密区的 ±0.5 mm 参数；层数、层厚和 Z 向
       逐层选取公式保持不变。

重要说明：
    原教程随机点判据为“球心距离 >= 0.5R”，这会允许弹丸严重相互穿透，
    且 50 个半径 0.15 mm 的弹丸无法无重叠地放入教程给定的小区域。
    为忠实复刻，默认 MIN_CENTER_DISTANCE_FACTOR = 0.5。
    若要得到物理上无初始穿透的弹丸云，请改为 2.02，并同步扩大喷射区域
    或提高弹丸云高度。
"""

from __future__ import print_function

from abaqus import mdb, session
from abaqusConstants import *
from caeModules import *
import mesh
import regionToolset

import os
import io
import math
import random
import traceback


# ============================================================================
#                         一、用户可修改参数区
# ============================================================================
# 绝大多数日常修改只需要在本区域完成。除非需要改变算法，否则不建议修改后文。

# ---------------------------------------------------------------------------
# 1. 模型、文件和作业控制
# ---------------------------------------------------------------------------
MODEL_NAME = 'Model-1'
JOB_NAME = 'Job-1'

# True：若同名模型已存在，删除旧模型后重新生成。
# False：同名模型存在时直接报错，防止误覆盖。
DELETE_EXISTING_MODEL = False

# 是否创建 Job；是否在脚本结束时直接提交计算。
CREATE_JOB = True
SUBMIT_JOB = False

# 是否保存 CAE。相对路径会相对于 Abaqus 当前工作目录转换为绝对路径。
SAVE_CAE = True
OUTPUT_CAE_PATH = 'Random-Shot.cae'

# 是否把本次生成的弹丸球心坐标输出为 CSV，便于复现实验与核查随机位置。
EXPORT_SHOT_COORDINATES = True
SHOT_COORDINATE_CSV = 'Random-Shot-Centers.csv'

# Job 计算资源。numDomains 通常设置为与 numCpus 相同。
NUM_CPUS = 18
MEMORY_PERCENT = 90
EXPLICIT_PRECISION = SINGLE       # 可改为 DOUBLE
NODAL_OUTPUT_PRECISION = SINGLE   # 可改为 FULL


# ---------------------------------------------------------------------------
# 2. 靶材几何参数（PPT 原值：2 mm × 2 mm × 1 mm）
# ---------------------------------------------------------------------------
TARGET_LENGTH = 2.0       # X 方向长度
TARGET_WIDTH = 2.0        # Y 方向宽度
TARGET_HEIGHT = 1.0       # Z 方向厚度；靶材顶面坐标为 Z = TARGET_HEIGHT

# 中心加密区域：X ∈ [-0.4, 0.4]，Y ∈ [-0.4, 0.4]
FINE_ZONE_HALF_X = 0.5
FINE_ZONE_HALF_Y = 0.5

# 靶材顶面向下的加密深度。PPT 的 30 层 × 0.01 mm 对应 0.30 mm。
FINE_ZONE_DEPTH = 0.30

# 靶材网格尺寸。0.01 mm 的模型规模较大，首次调试可先改成 0.02～0.03 mm。
TARGET_FINE_MESH_SIZE = 0.01
TARGET_COARSE_MESH_SIZE = 0.10

# STRUCTURED_HEX：复刻教程的分区六面体网格。
# FREE_TET：自由四面体网格，建模更稳健，但与教程网格拓扑不同。
TARGET_MESH_STRATEGY = 'STRUCTURED_HEX'
TARGET_TET_MESH_SIZE = 0.03


# ---------------------------------------------------------------------------
# 3. 弹丸几何、数量、随机区域和网格
# ---------------------------------------------------------------------------
SHOT_RADIUS = 0.15
SHOT_NUMBER = 50
SHOT_MESH_SIZE = 0.01

# 为确保所有弹丸沿 Z 方向下落后完整落在靶材中心加密区域内，
# 球心随机边界必须从中心区域边界向内缩进一个弹丸半径。
# 默认中心区域为 X、Y ∈ [-0.4, 0.4]，R=0.15 mm，因此球心范围为
# X、Y ∈ [-0.25, 0.25]。这样弹丸最外缘不会越过中心区域边界。
SHOT_X_MIN = -FINE_ZONE_HALF_X + SHOT_RADIUS
SHOT_X_MAX =  FINE_ZONE_HALF_X - SHOT_RADIUS
SHOT_Y_MIN = -FINE_ZONE_HALF_Y + SHOT_RADIUS
SHOT_Y_MAX =  FINE_ZONE_HALF_Y - SHOT_RADIUS

# 教程中 Z ∈ [1.01+R, 2.0]。这里写成与靶材顶面关联的形式。
SHOT_Z_MIN = TARGET_HEIGHT + 0.01 + SHOT_RADIUS
SHOT_Z_MAX = 2.0

# 固定随机种子可保证每次弹丸位置完全相同；设为 None 则每次不同。
RANDOM_SEED = 2026

# 最小球心距离 = SHOT_RADIUS × 此系数。
# 0.5：严格复刻教程代码，允许大范围初始穿透。
# 2.0：球面刚好相切；建议使用 2.02 留出少量几何间隙。
MIN_CENTER_DISTANCE_FACTOR = 0.5

# 随机放置的最大尝试次数。若无重叠要求过严，脚本会在达到此次数后停止。
MAX_PLACEMENT_ATTEMPTS = 1000000


# ---------------------------------------------------------------------------
# 4. 靶材 TC4 参数（按 PPT 截图录入）
# ---------------------------------------------------------------------------
TARGET_MATERIAL_NAME = 'Material-target-TC4'
TARGET_SECTION_NAME = 'Section-target-TC4'
TARGET_DENSITY = 4.5e-9       # tonne/mm^3
TARGET_E = 1.2e5              # MPa
TARGET_NU = 0.34

# Johnson-Cook：sigma = (A + B*eps_p^n)*(1 + C*ln(epsdot/epsdot0))*温度项
TARGET_JC_A = 961          # MPa
TARGET_JC_B = 902          # MPa
TARGET_JC_N = 0.87
TARGET_JC_M = 0.0
TARGET_MELTING_TEMPERATURE = 0.0
TARGET_TRANSITION_TEMPERATURE = 0.0
TARGET_JC_C = 0.01
TARGET_REFERENCE_STRAIN_RATE = 1.0

# 若没有可靠的温度参数，保持 False。设为 True 才启用绝热升温选项。
INCLUDE_ADIABATIC_HEATING = False


# ---------------------------------------------------------------------------
# 5. 弹丸 S110 参数（按 PPT 截图录入）
# ---------------------------------------------------------------------------
SHOT_MATERIAL_NAME = 'Material-shot-S110'
SHOT_SECTION_NAME = 'Section-shot-S110'
SHOT_DENSITY = 7.85e-9        # tonne/mm^3
SHOT_E = 210000.0             # MPa
SHOT_NU = 0.30

# 弹丸通过 RigidBody 约束作为刚体运动；密度用于计算刚体质量和转动惯量。


# ---------------------------------------------------------------------------
# 6. 喷丸速度、分析步和接触
# ---------------------------------------------------------------------------
STEP_NAME = 'Step-1'

# PPT/视频中为 V3 = -80000 mm/s，即 80 m/s，垂直冲击靶材顶面。
SHOT_VELOCITY_X = 0.0
SHOT_VELOCITY_Y = 0.0
SHOT_VELOCITY_Z = -59090.0

# AUTO：按“最高弹丸到靶面距离 / 速度 × 安全系数”自动计算。
# MANUAL：直接采用 MANUAL_STEP_TIME。
STEP_TIME_MODE = 'AUTO'
STEP_TIME_SAFETY_FACTOR = 1.2   # 视频中 1.4268e-5 最终取约 1.8e-5
MANUAL_STEP_TIME = 1.8e-5

# 接触属性。视频中摩擦系数为 0.2，法向为 Hard Contact，可分离。
FRICTION_COEFFICIENT = 0.20
CONTACT_PROPERTY_NAME = 'IntProp-1'
CONTACT_NAME = 'Int-1'

# TUTORIAL_SELECTED：仅包含操作记录指定的“弹丸-靶材顶面”接触对。
# ALL_EXTERIOR：使用 General Contact 的 All with self。
# 根据 PythonReader 记录，默认显式创建全部弹丸外表面和靶材顶面。
CONTACT_SCOPE = 'TUTORIAL_SELECTED'

# 以下名称和掩码直接来自“修改.txt”中已经在 Abaqus/CAE 2025 成功执行的命令。
# 本脚本生成的弹丸和靶材拓扑固定，因此各同源实例可安全复用相同掩码。
SHOT_ASSEMBLY_SURFACE_NAME = 'Surf-shot'
TARGET_TOP_SURFACE_NAME = 'Surf-target'
SHOT_EXTERIOR_FACE_MASK = ('[#e7030 ]',)
TARGET_TOP_FACE_MASK = ('[#0 #200000 ]',)


# ---------------------------------------------------------------------------
# 7. 约束、输出和分层网格集
# ---------------------------------------------------------------------------
# 教程中靶材底面和四个侧面完全固定，顶面保持自由。
FIX_TARGET_BOTTOM = True
FIX_TARGET_SIDES = True

# ODB 场输出时间间隔数量。PythonReader 记录中改为 200，并按节点输出。
FIELD_OUTPUT_INTERVALS = 200
FIELD_OUTPUT_POSITION = NODES
FIELD_OUTPUT_VARIABLES = ('S', 'PEEQ', 'LE', 'U', 'V', 'A', 'RF', 'CSTRESS')

# 结果提取分层集：复刻参考 random-shot.py 末尾的集合范围。
# X、Y 均为 [-0.4, 0.4] mm；顶面向下 30 层，每层 0.01 mm。
CREATE_LAYER_ELEMENT_SETS = True
LAYER_SET_COUNT = 30
LAYER_SET_THICKNESS = 0.01
LAYER_SET_HALF_X = 0.4
LAYER_SET_HALF_Y = 0.4
LAYER_SET_PREFIX = 'SET-LAYER-'
LAYER_SELECTION_TOLERANCE = 1.0e-4


# ---------------------------------------------------------------------------
# 8. 调试和显示
# ---------------------------------------------------------------------------
VERBOSE = True
DISPLAY_FINAL_ASSEMBLY = True


# ============================================================================
#                         二、内部辅助函数
# ============================================================================

def log(message):
    """统一输出运行信息。"""
    if VERBOSE:
        print('[RandomShot] {0}'.format(message))


def absolute_path(path_text):
    """将用户输入的相对路径转换成绝对路径。"""
    return os.path.abspath(os.path.expanduser(path_text))


def validate_parameters():
    """在真正建模前检查参数，避免运行到中途才发现明显错误。"""
    if SHOT_RADIUS <= 0.0:
        raise ValueError('SHOT_RADIUS 必须大于 0。')
    if SHOT_NUMBER <= 0:
        raise ValueError('SHOT_NUMBER 必须为正整数。')
    if CONTACT_SCOPE.upper() not in ('TUTORIAL_SELECTED', 'ALL_EXTERIOR'):
        raise ValueError("CONTACT_SCOPE 只能为 'TUTORIAL_SELECTED' 或 'ALL_EXTERIOR'。")
    if TARGET_LENGTH <= 0.0 or TARGET_WIDTH <= 0.0 or TARGET_HEIGHT <= 0.0:
        raise ValueError('靶材长、宽、高必须大于 0。')
    if not (0.0 < FINE_ZONE_HALF_X < TARGET_LENGTH / 2.0):
        raise ValueError('FINE_ZONE_HALF_X 必须位于 (0, TARGET_LENGTH/2) 内。')
    if not (0.0 < FINE_ZONE_HALF_Y < TARGET_WIDTH / 2.0):
        raise ValueError('FINE_ZONE_HALF_Y 必须位于 (0, TARGET_WIDTH/2) 内。')
    if not (0.0 < FINE_ZONE_DEPTH < TARGET_HEIGHT):
        raise ValueError('FINE_ZONE_DEPTH 必须位于 (0, TARGET_HEIGHT) 内。')
    if SHOT_RADIUS >= FINE_ZONE_HALF_X or SHOT_RADIUS >= FINE_ZONE_HALF_Y:
        raise ValueError(
            '弹丸半径过大：为保证整个弹丸落在中心区域内，必须满足 '
            'SHOT_RADIUS < FINE_ZONE_HALF_X 且 SHOT_RADIUS < FINE_ZONE_HALF_Y。'
        )
    if SHOT_X_MIN >= SHOT_X_MAX or SHOT_Y_MIN >= SHOT_Y_MAX or SHOT_Z_MIN >= SHOT_Z_MAX:
        raise ValueError('弹丸随机区域上下限设置错误。')
    if SHOT_VELOCITY_Z >= 0.0:
        raise ValueError('当前模型靶面法向为 +Z，向下冲击应设置 SHOT_VELOCITY_Z < 0。')
    if MIN_CENTER_DISTANCE_FACTOR < 2.0:
        log('警告：MIN_CENTER_DISTANCE_FACTOR < 2.0，初始弹丸会发生几何穿透。')
    if TARGET_FINE_MESH_SIZE <= 0.0 or TARGET_COARSE_MESH_SIZE <= 0.0:
        raise ValueError('网格尺寸必须大于 0。')
    if TARGET_MESH_STRATEGY not in ('STRUCTURED_HEX', 'FREE_TET'):
        raise ValueError("TARGET_MESH_STRATEGY 只能为 'STRUCTURED_HEX' 或 'FREE_TET'。")
    if CREATE_LAYER_ELEMENT_SETS:
        total_layer_depth = LAYER_SET_COUNT * LAYER_SET_THICKNESS
        if not (0.0 < LAYER_SET_HALF_X <= TARGET_LENGTH / 2.0):
            raise ValueError('LAYER_SET_HALF_X 必须位于 (0, TARGET_LENGTH/2] 内。')
        if not (0.0 < LAYER_SET_HALF_Y <= TARGET_WIDTH / 2.0):
            raise ValueError('LAYER_SET_HALF_Y 必须位于 (0, TARGET_WIDTH/2] 内。')
        if total_layer_depth > TARGET_HEIGHT + 1.0e-12:
            raise ValueError('分层网格集总深度超过靶材厚度。')


def _repository_contains(repository, key):
    """
    稳健检查 Abaqus Repository 中是否存在指定键。

    不直接依赖 ``key in repository.keys()``，因为部分 Abaqus 版本中
    Repository 的 keys() 视图在模型新建、删除或重命名后可能出现短暂不同步。
    """
    try:
        repository[key]
        return True
    except KeyError:
        return False


def _unique_temporary_model_name():
    """生成不会与现有模型重名的临时模型名称。"""
    base_name = '__RandomShot_Temporary_Model__'
    candidate = base_name
    index = 1

    while _repository_contains(mdb.models, candidate):
        candidate = '{0}_{1}'.format(base_name, index)
        index += 1

    return candidate


def prepare_model():
    """
    创建全新的目标模型，并稳健处理同名旧模型。

    Abaqus/CAE 新建模型数据库时通常只有默认模型 ``Model-1``。
    某些版本不允许先删除数据库中的唯一模型，因此下面采用“临时替代模型”流程：

        1. 创建一个临时空模型；
        2. 删除同名旧模型；
        3. 将临时模型重命名为 MODEL_NAME。

    该流程保证模型数据库中始终至少存在一个模型，可避免：
        KeyError: 'Model-1'
    """
    model_exists = _repository_contains(mdb.models, MODEL_NAME)

    if not model_exists:
        model = mdb.Model(name=MODEL_NAME)
        log('已创建模型：{0}'.format(MODEL_NAME))
        return model

    if not DELETE_EXISTING_MODEL:
        raise RuntimeError(
            '模型 {0} 已存在。请修改 MODEL_NAME 或启用 DELETE_EXISTING_MODEL。'.format(
                MODEL_NAME
            )
        )

    temporary_name = _unique_temporary_model_name()
    log('检测到同名旧模型：{0}'.format(MODEL_NAME))
    log('先创建临时替代模型：{0}'.format(temporary_name))

    # 先创建临时模型，确保即使 MODEL_NAME 是数据库中的唯一模型，
    # 删除它时模型数据库中仍保留至少一个模型。
    mdb.Model(name=temporary_name)

    try:
        del mdb.models[MODEL_NAME]
        log('已删除同名旧模型：{0}'.format(MODEL_NAME))
    except Exception:
        # 删除失败时清理刚创建的临时模型，避免污染当前数据库。
        try:
            if _repository_contains(mdb.models, temporary_name):
                del mdb.models[temporary_name]
        except Exception:
            pass
        raise RuntimeError(
            '无法删除同名旧模型 {0}。请关闭正在引用该模型的作业、视口或对话框后重试。'
            .format(MODEL_NAME)
        )

    try:
        mdb.models.changeKey(fromName=temporary_name, toName=MODEL_NAME)
    except Exception:
        # 极少数版本若不支持 changeKey，则退回到“新建目标模型后删除临时模型”。
        model = mdb.Model(name=MODEL_NAME)
        if _repository_contains(mdb.models, temporary_name):
            del mdb.models[temporary_name]
        log('已创建模型：{0}'.format(MODEL_NAME))
        return model

    model = mdb.models[MODEL_NAME]
    log('已创建全新模型：{0}'.format(MODEL_NAME))
    return model


def create_materials_and_sections(model):
    """创建 TC4 靶材和 S110 弹丸材料，并创建实体截面。"""
    log('创建靶材 TC4 的密度、弹性和 Johnson-Cook 塑性参数。')
    target_mat = model.Material(name=TARGET_MATERIAL_NAME)
    target_mat.Density(table=((TARGET_DENSITY,),))
    target_mat.Elastic(table=((TARGET_E, TARGET_NU),))
    target_mat.Plastic(
        hardening=JOHNSON_COOK,
        table=((TARGET_JC_A,
                TARGET_JC_B,
                TARGET_JC_N,
                TARGET_JC_M,
                TARGET_MELTING_TEMPERATURE,
                TARGET_TRANSITION_TEMPERATURE),)
    )
    target_mat.plastic.RateDependent(
        type=JOHNSON_COOK,
        table=((TARGET_JC_C, TARGET_REFERENCE_STRAIN_RATE),)
    )

    model.HomogeneousSolidSection(
        name=TARGET_SECTION_NAME,
        material=TARGET_MATERIAL_NAME,
        thickness=None
    )

    log('创建弹丸 S110 的密度和线弹性参数。')
    shot_mat = model.Material(name=SHOT_MATERIAL_NAME)
    shot_mat.Density(table=((SHOT_DENSITY,),))
    shot_mat.Elastic(table=((SHOT_E, SHOT_NU),))

    model.HomogeneousSolidSection(
        name=SHOT_SECTION_NAME,
        material=SHOT_MATERIAL_NAME,
        thickness=None
    )


def collect_exterior_geometry_faces(face_array):
    """
    从几何 FaceArray 中筛选实体外表面。

    球体被三个基准面分区后，part.faces / instance.faces 同时包含：
        1. 球体真正的外表面；
        2. 各 Cell 之间的内部公共分区面。

    外表面只邻接 1 个 Cell，内部公共面通常邻接 2 个 Cell。
    若把内部面也传给 Surface(side1Faces=...)，部分 Abaqus 版本会直接报
    "Feature creation failed / 特征创建失败"。
    """
    exterior_faces = ()
    unresolved_faces = 0

    for face in face_array:
        try:
            adjacent_cells = face.getCells()
            if len(adjacent_cells) == 1:
                exterior_faces = exterior_faces + (face,)
        except Exception:
            unresolved_faces += 1

    if unresolved_faces > 0:
        log('警告：有 {0} 个几何面无法读取相邻 Cell。'.format(unresolved_faces))

    return exterior_faces


def create_shot_part(model):
    """
    创建球形弹丸 Part。

    流程与 PPT 第 2～3 页一致：
        半圆截面绕中心轴旋转 360° -> 球体；
        球心创建参考点 -> 后续刚体控制点；
        通过 X=0、Y=0、Z=0 三个平面将球体划分为 8 个八分体；
        赋予 S110 截面并划分网格。
    """
    log('创建球形弹丸 Part。')

    sketch_name = '__shot_profile__'
    sketch = model.ConstrainedSketch(
        name=sketch_name,
        sheetSize=max(2.0, SHOT_RADIUS * 10.0)
    )
    sketch.ConstructionLine(
        point1=(0.0, SHOT_RADIUS * 2.0),
        point2=(0.0, -SHOT_RADIUS * 2.0)
    )
    sketch.ArcByCenterEnds(
        center=(0.0, 0.0),
        point1=(0.0, SHOT_RADIUS),
        point2=(0.0, -SHOT_RADIUS),
        direction=CLOCKWISE
    )
    sketch.Line(
        point1=(0.0, SHOT_RADIUS),
        point2=(0.0, -SHOT_RADIUS)
    )

    part = model.Part(
        name='shot',
        dimensionality=THREE_D,
        type=DEFORMABLE_BODY
    )
    part.BaseSolidRevolve(
        sketch=sketch,
        angle=360.0,
        flipRevolveDirection=OFF
    )
    del model.sketches[sketch_name]

    # 返回的 Feature.id 是稳定的参考点键值，避免硬编码原教程中的 rp[2]。
    rp_feature = part.ReferencePoint(point=(0.0, 0.0, 0.0))
    shot_reference_point_id = rp_feature.id

    # 三个主平面将球分成 8 个拓扑近似六面体区域，便于生成较规整网格。
    for principal_plane in (YZPLANE, XZPLANE, XYPLANE):
        datum_feature = part.DatumPlaneByPrincipalPlane(
            principalPlane=principal_plane,
            offset=0.0
        )
        part.PartitionCellByDatumPlane(
            datumPlane=part.datums[datum_feature.id],
            cells=part.cells[:]
        )

    shot_region = regionToolset.Region(cells=part.cells[:])
    part.SectionAssignment(
        region=shot_region,
        sectionName=SHOT_SECTION_NAME
    )

    # 为弹丸设置 Explicit 实体单元。球体作为刚体后不会发生材料变形，
    # 但其离散表面用于接触，密度用于计算质量/转动惯量。
    elem_type_hex = mesh.ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT)
    elem_type_wedge = mesh.ElemType(elemCode=C3D6, elemLibrary=EXPLICIT)
    elem_type_tet = mesh.ElemType(elemCode=C3D4, elemLibrary=EXPLICIT)
    part.setElementType(
        regions=(part.cells[:],),
        elemTypes=(elem_type_hex, elem_type_wedge, elem_type_tet)
    )

    # 优先使用结构化六面体；若某版本不接受，则退回扫掠六面体主导网格。
    try:
        part.setMeshControls(
            regions=part.cells[:],
            elemShape=HEX,
            technique=STRUCTURED
        )
    except Exception:
        log('弹丸 STRUCTURED HEX 控制未被当前版本接受，改用 HEX_DOMINATED + SWEEP。')
        part.setMeshControls(
            regions=part.cells[:],
            elemShape=HEX_DOMINATED,
            technique=SWEEP
        )

    part.seedPart(
        size=SHOT_MESH_SIZE,
        deviationFactor=0.1,
        minSizeFactor=0.1
    )
    part.generateMesh()

    log('弹丸网格完成：节点 {0}，单元 {1}。'.format(len(part.nodes), len(part.elements)))
    return part, shot_reference_point_id


def create_target_part(model):
    """
    创建靶材 Part，并按教程形成中心十字加密区域。

    分区平面：
        X = ±FINE_ZONE_HALF_X；
        Y = ±FINE_ZONE_HALF_Y；
        Z = TARGET_HEIGHT - FINE_ZONE_DEPTH。

    形成 3 × 3 × 2 = 18 个块状 Cell。结构化六面体模式下：
        中心 X 段和中心 Y 段采用细网格；
        顶部厚度段 Z 方向采用细网格；
        外围采用粗网格。
    """
    log('创建并分区靶材 Part。')

    sketch_name = '__target_profile__'
    sketch = model.ConstrainedSketch(
        name=sketch_name,
        sheetSize=max(TARGET_LENGTH, TARGET_WIDTH) * 2.0
    )
    sketch.rectangle(
        point1=(-TARGET_LENGTH / 2.0, -TARGET_WIDTH / 2.0),
        point2=( TARGET_LENGTH / 2.0,  TARGET_WIDTH / 2.0)
    )

    part = model.Part(
        name='target',
        dimensionality=THREE_D,
        type=DEFORMABLE_BODY
    )
    part.BaseSolidExtrude(sketch=sketch, depth=TARGET_HEIGHT)
    del model.sketches[sketch_name]

    # X 常量平面：YZPLANE；Y 常量平面：XZPLANE；Z 常量平面：XYPLANE。
    partition_specs = (
        (YZPLANE, -FINE_ZONE_HALF_X),
        (YZPLANE,  FINE_ZONE_HALF_X),
        (XZPLANE, -FINE_ZONE_HALF_Y),
        (XZPLANE,  FINE_ZONE_HALF_Y),
        (XYPLANE, TARGET_HEIGHT - FINE_ZONE_DEPTH),
    )
    for principal_plane, offset_value in partition_specs:
        datum_feature = part.DatumPlaneByPrincipalPlane(
            principalPlane=principal_plane,
            offset=offset_value
        )
        part.PartitionCellByDatumPlane(
            datumPlane=part.datums[datum_feature.id],
            cells=part.cells[:]
        )

    target_region = regionToolset.Region(cells=part.cells[:])
    part.SectionAssignment(
        region=target_region,
        sectionName=TARGET_SECTION_NAME
    )

    if TARGET_MESH_STRATEGY == 'STRUCTURED_HEX':
        mesh_target_structured_hex(part)
    else:
        mesh_target_free_tet(part)

    log('靶材网格完成：节点 {0}，单元 {1}。'.format(len(part.nodes), len(part.elements)))
    return part


def edge_end_coordinates(part, edge):
    """返回直线 Edge 两端点坐标；若无法获得两个端点则返回 None。"""
    vertex_ids = edge.getVertices()
    if len(vertex_ids) != 2:
        return None
    point_0 = part.vertices[vertex_ids[0]].pointOn[0]
    point_1 = part.vertices[vertex_ids[1]].pointOn[0]
    return point_0, point_1


def mesh_target_structured_hex(part):
    """按教程风格划分靶材结构化六面体网格。"""
    log('设置靶材结构化六面体网格。')

    elem_type_hex = mesh.ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT)
    elem_type_wedge = mesh.ElemType(elemCode=C3D6, elemLibrary=EXPLICIT)
    elem_type_tet = mesh.ElemType(elemCode=C3D4, elemLibrary=EXPLICIT)
    part.setElementType(
        regions=(part.cells[:],),
        elemTypes=(elem_type_hex, elem_type_wedge, elem_type_tet)
    )
    part.setMeshControls(
        regions=part.cells[:],
        elemShape=HEX,
        technique=STRUCTURED
    )

    # 先全局粗划分，再对“中心 X 段、中心 Y 段和顶部 Z 段”的边覆盖细网格。
    part.seedPart(
        size=TARGET_COARSE_MESH_SIZE,
        deviationFactor=0.1,
        minSizeFactor=0.1
    )

    fine_edges = []
    tolerance = 1.0e-8
    z_partition = TARGET_HEIGHT - FINE_ZONE_DEPTH

    for edge in part.edges:
        coordinates = edge_end_coordinates(part, edge)
        if coordinates is None:
            continue

        point_0, point_1 = coordinates
        dx = abs(point_1[0] - point_0[0])
        dy = abs(point_1[1] - point_0[1])
        dz = abs(point_1[2] - point_0[2])
        midpoint = (
            0.5 * (point_0[0] + point_1[0]),
            0.5 * (point_0[1] + point_1[1]),
            0.5 * (point_0[2] + point_1[2])
        )

        # X 向边：仅中心 X 分段（-hx 到 +hx）细化。
        if dx >= dy and dx >= dz:
            if abs(midpoint[0]) < FINE_ZONE_HALF_X - tolerance:
                fine_edges.append(edge)

        # Y 向边：仅中心 Y 分段（-hy 到 +hy）细化。
        elif dy >= dx and dy >= dz:
            if abs(midpoint[1]) < FINE_ZONE_HALF_Y - tolerance:
                fine_edges.append(edge)

        # Z 向边：仅顶面向下 FINE_ZONE_DEPTH 的分段细化。
        else:
            if midpoint[2] > z_partition + tolerance:
                fine_edges.append(edge)

    if fine_edges:
        part.seedEdgeBySize(
            edges=tuple(fine_edges),
            size=TARGET_FINE_MESH_SIZE,
            deviationFactor=0.1,
            minSizeFactor=0.1,
            constraint=FINER
        )
    else:
        raise RuntimeError('未识别到靶材细化边，请检查分区尺寸和 Abaqus 版本。')

    part.generateMesh()


def mesh_target_free_tet(part):
    """备用自由四面体网格模式，建模稳定但不复刻教程的六面体拓扑。"""
    log('设置靶材自由四面体网格。')

    elem_type_tet = mesh.ElemType(elemCode=C3D4, elemLibrary=EXPLICIT)
    part.setElementType(
        regions=(part.cells[:],),
        elemTypes=(elem_type_tet,)
    )
    part.setMeshControls(
        regions=part.cells[:],
        elemShape=TET,
        technique=FREE
    )
    part.seedPart(
        size=TARGET_TET_MESH_SIZE,
        deviationFactor=0.1,
        minSizeFactor=0.1
    )
    part.generateMesh()


def squared_distance(point_a, point_b):
    """计算两点距离平方，随机筛选时避免重复开平方。"""
    dx = point_a[0] - point_b[0]
    dy = point_a[1] - point_b[1]
    dz = point_a[2] - point_b[2]
    return dx * dx + dy * dy + dz * dz


def generate_random_shot_centers():
    """
    随机生成弹丸球心。

    与原教程不同的是：
        - 采用固定随机种子选项，保证可重复；
        - 设置最大尝试次数，避免 while 无限循环；
        - 使用距离平方提高效率；
        - 对无法满足的无重叠要求给出明确报错。
    """
    if RANDOM_SEED is not None:
        random.seed(RANDOM_SEED)
        log('随机种子：{0}'.format(RANDOM_SEED))

    minimum_distance = SHOT_RADIUS * MIN_CENTER_DISTANCE_FACTOR
    minimum_distance_squared = minimum_distance * minimum_distance

    points = []
    attempts = 0

    while len(points) < SHOT_NUMBER and attempts < MAX_PLACEMENT_ATTEMPTS:
        attempts += 1
        candidate = (
            random.uniform(SHOT_X_MIN, SHOT_X_MAX),
            random.uniform(SHOT_Y_MIN, SHOT_Y_MAX),
            random.uniform(SHOT_Z_MIN, SHOT_Z_MAX)
        )

        accepted = True
        for existing_point in points:
            if squared_distance(candidate, existing_point) < minimum_distance_squared:
                accepted = False
                break

        if accepted:
            points.append(candidate)
            if VERBOSE and (len(points) == 1 or len(points) % 10 == 0 or len(points) == SHOT_NUMBER):
                log('已生成弹丸球心：{0}/{1}'.format(len(points), SHOT_NUMBER))

    if len(points) != SHOT_NUMBER:
        raise RuntimeError(
            '随机放置失败：尝试 {0} 次后只生成 {1}/{2} 个弹丸。\n'
            '当前 X、Y 范围被严格限制在中心区域内，不建议向外扩大。请优先增大 '
            'SHOT_Z_MAX、减少 SHOT_NUMBER，或减小 MIN_CENTER_DISTANCE_FACTOR。'.format(
                attempts, len(points), SHOT_NUMBER
            )
        )

    # 二次复核：不仅球心，而且弹丸在 X-Y 平面上的完整投影都必须位于中心区域。
    containment_tolerance = 1.0e-12
    for index, point in enumerate(points, 1):
        x_coord = point[0]
        y_coord = point[1]
        if (x_coord - SHOT_RADIUS < -FINE_ZONE_HALF_X - containment_tolerance or
                x_coord + SHOT_RADIUS > FINE_ZONE_HALF_X + containment_tolerance or
                y_coord - SHOT_RADIUS < -FINE_ZONE_HALF_Y - containment_tolerance or
                y_coord + SHOT_RADIUS > FINE_ZONE_HALF_Y + containment_tolerance):
            raise RuntimeError(
                '第 {0} 个弹丸超出中心区域：球心=({1:.8g}, {2:.8g}, {3:.8g})。'.format(
                    index, point[0], point[1], point[2]
                )
            )

    log(
        '随机球心生成完成，总尝试次数：{0}；全部弹丸完整投影均位于中心区域内。'.format(
            attempts
        )
    )
    log(
        '中心区域：X=[-{0}, {0}]，Y=[-{1}, {1}] mm；球心范围：'
        'X=[{2}, {3}]，Y=[{4}, {5}] mm。'.format(
            FINE_ZONE_HALF_X, FINE_ZONE_HALF_Y,
            SHOT_X_MIN, SHOT_X_MAX, SHOT_Y_MIN, SHOT_Y_MAX
        )
    )
    return points


def export_shot_centers(points):
    """将球心坐标保存为 UTF-8 CSV。"""
    if not EXPORT_SHOT_COORDINATES:
        return

    csv_path = absolute_path(SHOT_COORDINATE_CSV)
    with io.open(csv_path, mode='w', encoding='utf-8') as csv_file:
        csv_file.write(u'序号,球心X_mm,球心Y_mm,球心Z_mm\n')
        for index, point in enumerate(points, 1):
            csv_file.write(u'{0},{1:.12g},{2:.12g},{3:.12g}\n'.format(
                index, point[0], point[1], point[2]
            ))
    log('弹丸球心坐标已输出：{0}'.format(csv_path))


def concatenate_objects(current_tuple, object_array):
    """把 Abaqus 几何数组逐项加入普通 tuple，便于跨多个实例创建集合/表面。"""
    result = current_tuple
    for item in object_array:
        result = result + (item,)
    return result


def create_assembly_and_rigid_bodies(model, shot_part, target_part,
                                      shot_reference_point_id, shot_centers):
    """装配靶材和所有弹丸，并为每个弹丸建立独立刚体约束。"""
    log('创建装配、弹丸实例和刚体约束。')

    assembly = model.rootAssembly
    assembly.DatumCsysByDefault(CARTESIAN)

    target_instance_name = 'target-1'
    assembly.Instance(
        name=target_instance_name,
        part=target_part,
        dependent=ON
    )

    shot_reference_points = ()
    # 必须保留为 Abaqus FaceArray（或多个 FaceArray 的相加结果）。
    # assembly.Surface 不接受这里原先构造的跨实例普通 Python tuple。
    shot_exterior_faces = None
    shot_surface_instance_count = 0
    shot_instance_names = []

    for index, center in enumerate(shot_centers):
        instance_name = 'shot-{0:03d}'.format(index + 1)
        shot_instance_names.append(instance_name)

        assembly.Instance(
            name=instance_name,
            part=shot_part,
            dependent=ON
        )
        assembly.translate(
            instanceList=(instance_name,),
            vector=center
        )

        instance = assembly.instances[instance_name]
        rp = instance.referencePoints[shot_reference_point_id]
        shot_reference_points = shot_reference_points + (rp,)

        # 操作记录是在装配层级为每个弹丸用同一掩码选择 8 个球面片，
        # 再将各 FaceArray 直接相加，最终一次性创建 50 × 8 = 400 个面的表面。
        if CONTACT_SCOPE.upper() == 'TUTORIAL_SELECTED':
            instance_exterior_faces = instance.faces.getSequenceFromMask(
                mask=SHOT_EXTERIOR_FACE_MASK
            )
            if len(instance_exterior_faces) != 8:
                raise RuntimeError(
                    '{0} 的外表面掩码应选中 8 个面，实际选中 {1} 个；'
                    '请重新录制 SHOT_EXTERIOR_FACE_MASK。'.format(
                        instance_name, len(instance_exterior_faces)
                    )
                )
            if shot_exterior_faces is None:
                shot_exterior_faces = instance_exterior_faces
            else:
                shot_exterior_faces = shot_exterior_faces + instance_exterior_faces
            shot_surface_instance_count += 1


        model.RigidBody(
            name='Rigid-{0:03d}'.format(index + 1),
            refPointRegion=regionToolset.Region(referencePoints=(rp,)),
            bodyRegion=regionToolset.Region(cells=instance.cells[:])
        )

    assembly.Set(
        name='SET-SHOT-RP',
        referencePoints=shot_reference_points
    )

    target_instance = assembly.instances[target_instance_name]

    # 指定接触模式下，直接把所有实例外表面创建成一个装配级表面。
    # 不使用 SurfaceByBoolean，避免旧版本中布尔合并表面失败。
    if CONTACT_SCOPE.upper() == 'TUTORIAL_SELECTED':
        if shot_surface_instance_count != len(shot_centers) or shot_exterior_faces is None:
            raise RuntimeError('未能收集全部弹丸实例的外表面。')

        assembly.Surface(
            name=SHOT_ASSEMBLY_SURFACE_NAME,
            side1Faces=shot_exterior_faces
        )
        log('已创建弹丸组合表面 {0}：{1} 个实例，共 {2} 个外表面。'.format(
            SHOT_ASSEMBLY_SURFACE_NAME,
            shot_surface_instance_count,
            len(shot_exterior_faces)
        ))

        target_top_faces = target_instance.faces.getSequenceFromMask(
            mask=TARGET_TOP_FACE_MASK
        )
        if len(target_top_faces) != 1:
            raise RuntimeError(
                '靶材顶面掩码应选中 1 个面，实际选中 {0} 个；'
                '请重新录制 TARGET_TOP_FACE_MASK。'.format(len(target_top_faces))
            )
        assembly.Surface(
            name=TARGET_TOP_SURFACE_NAME,
            side1Faces=target_top_faces
        )
        log('已创建靶材顶面接触表面：{0}。'.format(TARGET_TOP_SURFACE_NAME))
    else:
        log('ALL_EXTERIOR 模式：跳过装配级接触表面创建。')

    # 创建靶材全部单元集合，便于后处理和输出控制。
    assembly.Set(
        name='SET-TARGET-ALL',
        elements=target_instance.elements[:]
    )

    assembly.regenerate()
    log('已创建 {0} 个弹丸实例和 {0} 个刚体约束。'.format(len(shot_centers)))
    return assembly, target_instance_name, shot_instance_names


def calculate_step_time(shot_centers):
    """按最高弹丸的最上端到靶材顶面的距离自动计算分析时间。"""
    if STEP_TIME_MODE.upper() == 'MANUAL':
        if MANUAL_STEP_TIME <= 0.0:
            raise ValueError('MANUAL_STEP_TIME 必须大于 0。')
        return MANUAL_STEP_TIME

    vertical_speed = abs(SHOT_VELOCITY_Z)
    if vertical_speed <= 0.0:
        raise ValueError('自动计算分析时间时，SHOT_VELOCITY_Z 不能为 0。')

    # 视频教程的计算方式是：最高弹丸的最上端 Z 坐标减去靶面 Z 坐标，
    # 再除以喷射速度。也就是使用 center_z + R - target_top，而不是
    # 首次接触所需的 center_z - R - target_top。这样会留出更长的碰撞响应时间。
    maximum_surface_gap = 0.0
    for point in shot_centers:
        surface_gap = point[2] + SHOT_RADIUS - TARGET_HEIGHT
        maximum_surface_gap = max(maximum_surface_gap, surface_gap)

    if maximum_surface_gap < 0.0:
        log('警告：自动计算得到的弹丸云高度小于 0，请检查 Z 范围。')
        maximum_surface_gap = 0.0

    impact_time = maximum_surface_gap / vertical_speed
    calculated_time = impact_time * STEP_TIME_SAFETY_FACTOR

    # 防止所有弹丸几乎贴面时分析时间变成 0。
    calculated_time = max(calculated_time, 1.0e-9)

    log('最高弹丸最上端到靶面距离：{0:.8g} mm。'.format(maximum_surface_gap))
    log('自动分析时间：{0:.8g} s。'.format(calculated_time))
    return calculated_time


def create_explicit_step(model, step_time):
    """创建 Dynamic, Explicit 分析步。"""
    log('创建显式动力学分析步。')
    model.ExplicitDynamicsStep(
        name=STEP_NAME,
        previous='Initial',
        timePeriod=step_time,
        nlgeom=ON
    )


def create_contact(model, assembly):
    """
    创建 General Contact (Explicit)。

    “修改.txt”最终录制的接触对为：
        Surf-shot 与 Surf-target：弹丸-靶材。
    """
    log('创建接触属性和显式通用接触。')

    model.ContactProperty(CONTACT_PROPERTY_NAME)
    contact_property = model.interactionProperties[CONTACT_PROPERTY_NAME]
    contact_property.TangentialBehavior(
        formulation=PENALTY,
        directionality=ISOTROPIC,
        slipRateDependency=OFF,
        pressureDependency=OFF,
        temperatureDependency=OFF,
        dependencies=0,
        table=((FRICTION_COEFFICIENT,),),
        shearStressLimit=None,
        maximumElasticSlip=FRACTION,
        fraction=0.005,
        elasticSlipStiffness=None
    )
    contact_property.NormalBehavior(
        pressureOverclosure=HARD,
        allowSeparation=ON,
        constraintEnforcementMethod=DEFAULT
    )

    model.ContactExp(
        name=CONTACT_NAME,
        createStepName='Initial'
    )
    contact = model.interactions[CONTACT_NAME]

    if CONTACT_SCOPE.upper() == 'TUTORIAL_SELECTED':
        required_surfaces_exist = (
            _repository_contains(assembly.surfaces, SHOT_ASSEMBLY_SURFACE_NAME) and
            _repository_contains(assembly.surfaces, TARGET_TOP_SURFACE_NAME)
        )

        if required_surfaces_exist:
            # 与操作记录完全一致：在一次 setValuesInStep 调用中关闭 Allstar
            # 并通过 addPairs 加入唯一的弹丸/靶面接触对。
            contact.includedPairs.setValuesInStep(
                stepName='Initial',
                useAllstar=OFF,
                addPairs=((
                    assembly.surfaces[SHOT_ASSEMBLY_SURFACE_NAME],
                    assembly.surfaces[TARGET_TOP_SURFACE_NAME]
                ),)
            )
            log('General Contact 接触域：Surf-shot / Surf-target。')
        else:
            raise RuntimeError(
                '未找到 {0} 或 {1}，无法按操作记录设置接触对。'.format(
                    SHOT_ASSEMBLY_SURFACE_NAME, TARGET_TOP_SURFACE_NAME
                )
            )
    else:
        contact.includedPairs.setValuesInStep(
            stepName='Initial',
            useAllstar=ON
        )
        log('General Contact 接触域：All with self。')

    contact.contactPropertyAssignments.appendInStep(
        stepName='Initial',
        assignments=((GLOBAL, SELF, CONTACT_PROPERTY_NAME),)
    )

    # “修改.txt”最后一条相互作用命令。空字符串表示不指定磨损属性，
    # 但保留 Abaqus/CAE 为该通用接触生成的全局磨损表面属性赋值项。
    contact.wearSurfacePropertyAssignments.appendInStep(
        stepName='Initial',
        assignments=((GLOBAL, ''),)
    )


def create_boundary_conditions_and_velocity(model, assembly, target_instance_name):
    """
    固定靶材底面/侧面，并给所有弹丸参考点施加初速度。

    V2 修复：不再用几何 Face 创建固定 Set，而是按网格节点坐标筛选，
    再使用 SetFromNodeLabels 创建装配级节点集。
    """
    log('创建靶材固定约束和弹丸初速度。')

    target_instance = assembly.instances[target_instance_name]
    all_nodes = target_instance.nodes

    if len(all_nodes) == 0:
        raise RuntimeError('靶材实例没有网格节点，无法创建固定边界条件。')

    x_values = [node.coordinates[0] for node in all_nodes]
    y_values = [node.coordinates[1] for node in all_nodes]
    z_values = [node.coordinates[2] for node in all_nodes]

    x_min = min(x_values)
    x_max = max(x_values)
    y_min = min(y_values)
    y_max = max(y_values)
    z_min = min(z_values)
    z_max = max(z_values)

    coordinate_tolerance = max(
        1.0e-8,
        max(TARGET_LENGTH, TARGET_WIDTH, TARGET_HEIGHT) * 1.0e-6
    )

    fixed_node_labels = set()
    bottom_count = 0
    side_count = 0

    for node in all_nodes:
        x_coord = node.coordinates[0]
        y_coord = node.coordinates[1]
        z_coord = node.coordinates[2]

        is_bottom = abs(z_coord - z_min) <= coordinate_tolerance
        is_x_side = (
            abs(x_coord - x_min) <= coordinate_tolerance or
            abs(x_coord - x_max) <= coordinate_tolerance
        )
        is_y_side = (
            abs(y_coord - y_min) <= coordinate_tolerance or
            abs(y_coord - y_max) <= coordinate_tolerance
        )

        selected = False
        if FIX_TARGET_BOTTOM and is_bottom:
            selected = True
            bottom_count += 1

        if FIX_TARGET_SIDES and (is_x_side or is_y_side):
            selected = True
            side_count += 1

        if selected:
            fixed_node_labels.add(node.label)

    if len(fixed_node_labels) == 0:
        raise RuntimeError(
            '未识别到靶材固定节点。请检查 FIX_TARGET_BOTTOM、FIX_TARGET_SIDES和靶材坐标。'
        )

    fixed_node_labels = tuple(sorted(fixed_node_labels))
    fixed_set_name = 'SET-TARGET-FIXED-NODES'
    fixed_region = None

    try:
        if fixed_set_name in assembly.sets.keys():
            del assembly.sets[fixed_set_name]

        assembly.SetFromNodeLabels(
            name=fixed_set_name,
            nodeLabels=((target_instance_name, fixed_node_labels),),
            unsorted=True
        )
        fixed_region = assembly.sets[fixed_set_name]
        log('靶材固定节点集创建完成：{0}，唯一节点 {1} 个。'.format(
            fixed_set_name, len(fixed_node_labels)
        ))
    except Exception:
        log('警告：SetFromNodeLabels 创建持久节点集失败，改用临时节点 Region。')
        log(traceback.format_exc())

        selected_nodes = target_instance.nodes.sequenceFromLabels(
            labels=fixed_node_labels
        )
        if len(selected_nodes) == 0:
            raise RuntimeError('固定节点标签已识别，但 sequenceFromLabels 返回空节点数组。')
        fixed_region = regionToolset.Region(nodes=selected_nodes)

    model.EncastreBC(
        name='BC-TARGET-ENCASTRE',
        createStepName='Initial',
        region=fixed_region
    )

    if 'SET-SHOT-RP' not in assembly.sets.keys():
        raise RuntimeError('未找到弹丸参考点集合 SET-SHOT-RP，无法施加初速度。')

    model.Velocity(
        name='PREDEFINED-SHOT-VELOCITY',
        region=assembly.sets['SET-SHOT-RP'],
        velocity1=SHOT_VELOCITY_X,
        velocity2=SHOT_VELOCITY_Y,
        velocity3=SHOT_VELOCITY_Z,
        omega=0.0
    )

    log('固定节点筛选统计：底面命中 {0} 次，侧面命中 {1} 次，去重后 {2} 个节点。'.format(
        bottom_count, side_count, len(fixed_node_labels)
    ))
    log('靶材节点范围：X[{0:.8g}, {1:.8g}]，Y[{2:.8g}, {3:.8g}]，Z[{4:.8g}, {5:.8g}]。'.format(
        x_min, x_max, y_min, y_max, z_min, z_max
    ))

def create_layer_element_sets(assembly, target_instance_name):
    """
    按 PPT 第 11～12 页逻辑，在中心区域创建逐层单元集。

    第 1 层位于顶面下方 0～LAYER_SET_THICKNESS；
    第 2 层位于顶面下方 1～2 个层厚；以此类推。
    """
    if not CREATE_LAYER_ELEMENT_SETS:
        return

    log('创建结果提取分层单元集：X=[-{0}, {0}]，Y=[-{1}, {1}] mm。'.format(
        LAYER_SET_HALF_X, LAYER_SET_HALF_Y
    ))
    target_elements = assembly.instances[target_instance_name].elements
    tolerance = LAYER_SELECTION_TOLERANCE

    created_count = 0
    for layer_index in range(LAYER_SET_COUNT):
        z_max = TARGET_HEIGHT - layer_index * LAYER_SET_THICKNESS + tolerance
        z_min = TARGET_HEIGHT - (layer_index + 1) * LAYER_SET_THICKNESS - tolerance

        selected_elements = target_elements.getByBoundingBox(
            xMin=-LAYER_SET_HALF_X - tolerance,
            xMax= LAYER_SET_HALF_X + tolerance,
            yMin=-LAYER_SET_HALF_Y - tolerance,
            yMax= LAYER_SET_HALF_Y + tolerance,
            zMin=z_min,
            zMax=z_max
        )

        set_name = '{0}{1:02d}'.format(LAYER_SET_PREFIX, layer_index + 1)
        if len(selected_elements) > 0:
            assembly.Set(
                name=set_name,
                elements=selected_elements
            )
            created_count += 1
        else:
            log('警告：{0} 未选到单元；请检查网格尺寸和分层厚度。'.format(set_name))

    log('已创建分层单元集：{0}/{1}。'.format(created_count, LAYER_SET_COUNT))


def configure_output_requests(model):
    """设置场输出和整模型能量历史输出。"""
    log('配置 ODB 场输出和能量历史输出。')

    if 'F-Output-1' in model.fieldOutputRequests.keys():
        model.fieldOutputRequests['F-Output-1'].setValues(
            variables=FIELD_OUTPUT_VARIABLES,
            numIntervals=FIELD_OUTPUT_INTERVALS,
            position=FIELD_OUTPUT_POSITION
        )
    else:
        model.FieldOutputRequest(
            name='F-Output-1',
            createStepName=STEP_NAME,
            variables=FIELD_OUTPUT_VARIABLES,
            numIntervals=FIELD_OUTPUT_INTERVALS,
            position=FIELD_OUTPUT_POSITION
        )

    # 保留默认 H-Output-1 时不重复创建；另外创建能量监控，检查人工能量是否过大。
    if 'H-ENERGY' not in model.historyOutputRequests.keys():
        model.HistoryOutputRequest(
            name='H-ENERGY',
            createStepName=STEP_NAME,
            variables=('ALLIE', 'ALLKE', 'ALLAE', 'ALLVD', 'ALLWK', 'ETOTAL'),
            region=MODEL,
            numIntervals=FIELD_OUTPUT_INTERVALS
        )


def create_job_and_save(model):
    """创建分析作业、保存 CAE，并按用户选项提交。"""
    if CREATE_JOB:
        if JOB_NAME in mdb.jobs.keys():
            del mdb.jobs[JOB_NAME]

        log('创建 Job：{0}'.format(JOB_NAME))
        mdb.Job(
            name=JOB_NAME,
            model=MODEL_NAME,
            type=ANALYSIS,
            explicitPrecision=EXPLICIT_PRECISION,
            nodalOutputPrecision=NODAL_OUTPUT_PRECISION,
            memory=MEMORY_PERCENT,
            memoryUnits=PERCENTAGE,
            multiprocessingMode=DEFAULT,
            numCpus=NUM_CPUS,
            # Abaqus 2025 已移除 parallelizationMethodExplicit；
            # 显式并行由 numCpus、numDomains 和 multiprocessingMode 控制。
            numDomains=NUM_CPUS
        )

    if SAVE_CAE:
        cae_path = absolute_path(OUTPUT_CAE_PATH)
        log('保存 CAE：{0}'.format(cae_path))
        mdb.saveAs(pathName=cae_path)

    if CREATE_JOB and SUBMIT_JOB:
        log('提交 Job：{0}'.format(JOB_NAME))
        mdb.jobs[JOB_NAME].submit(consistencyChecking=OFF)
        log('Job 已提交。脚本不等待求解结束，可在 Job Monitor 中查看进度。')


def display_assembly(assembly):
    """在 CAE 主视口显示最终装配。"""
    if not DISPLAY_FINAL_ASSEMBLY:
        return

    try:
        viewport_name = session.currentViewportName
        viewport = session.viewports[viewport_name]
        viewport.setValues(displayedObject=assembly)
        viewport.view.fitView()
    except Exception:
        # 非图形模式运行时没有 viewport，忽略即可。
        pass


def print_summary(step_time, shot_centers, target_part, shot_part):
    """输出本次建模摘要。"""
    print('')
    print('=' * 72)
    print('随机喷丸模型创建完成')
    print('=' * 72)
    print('模型名称          : {0}'.format(MODEL_NAME))
    print('弹丸数量          : {0}'.format(len(shot_centers)))
    print('弹丸半径          : {0} mm'.format(SHOT_RADIUS))
    print('中心区域          : X=[-{0}, {0}], Y=[-{1}, {1}] mm'.format(
        FINE_ZONE_HALF_X, FINE_ZONE_HALF_Y))
    print('弹丸球心范围      : X=[{0}, {1}], Y=[{2}, {3}] mm'.format(
        SHOT_X_MIN, SHOT_X_MAX, SHOT_Y_MIN, SHOT_Y_MAX))
    print('初速度            : ({0}, {1}, {2}) mm/s'.format(
        SHOT_VELOCITY_X, SHOT_VELOCITY_Y, SHOT_VELOCITY_Z))
    print('分析步时间        : {0:.8g} s'.format(step_time))
    print('靶材节点/单元     : {0} / {1}'.format(len(target_part.nodes), len(target_part.elements)))
    print('单个弹丸节点/单元 : {0} / {1}'.format(len(shot_part.nodes), len(shot_part.elements)))
    print('摩擦系数          : {0}'.format(FRICTION_COEFFICIENT))
    print('接触域            : {0}'.format(CONTACT_SCOPE))
    if CREATE_JOB:
        print('Job 名称          : {0}'.format(JOB_NAME))
    if SAVE_CAE:
        print('CAE 路径          : {0}'.format(absolute_path(OUTPUT_CAE_PATH)))
    if EXPORT_SHOT_COORDINATES:
        print('球心坐标 CSV      : {0}'.format(absolute_path(SHOT_COORDINATE_CSV)))
    print('')
    print('计算前建议检查：')
    print('1. 在 Mesh 模块检查靶材结构化网格是否完整、无未划分区域。')
    print('2. 在 Interaction 模块检查 General Contact 和摩擦系数。')
    print('3. 若 MIN_CENTER_DISTANCE_FACTOR < 2，检查初始穿透及接触能量。')
    print('4. 求解后核对 ALLAE/ALLIE、ALLKE 和总能量平衡。')
    print('=' * 72)


# ============================================================================
#                         三、主程序
# ============================================================================

def main():
    """按既定顺序完成整个随机喷丸模型。"""
    validate_parameters()

    model = prepare_model()
    create_materials_and_sections(model)

    shot_part, shot_reference_point_id = create_shot_part(model)
    target_part = create_target_part(model)

    shot_centers = generate_random_shot_centers()
    export_shot_centers(shot_centers)

    assembly, target_instance_name, shot_instance_names = create_assembly_and_rigid_bodies(
        model=model,
        shot_part=shot_part,
        target_part=target_part,
        shot_reference_point_id=shot_reference_point_id,
        shot_centers=shot_centers
    )

    step_time = calculate_step_time(shot_centers)
    create_explicit_step(model, step_time)
    create_contact(model, assembly)
    create_boundary_conditions_and_velocity(model, assembly, target_instance_name)
    create_layer_element_sets(assembly, target_instance_name)
    configure_output_requests(model)

    assembly.regenerate()
    display_assembly(assembly)
    create_job_and_save(model)
    print_summary(step_time, shot_centers, target_part, shot_part)


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        print('')
        print('=' * 72)
        print('脚本执行失败：{0}'.format(error))
        print('以下为错误追踪信息：')
        traceback.print_exc()
        print('=' * 72)
        raise
