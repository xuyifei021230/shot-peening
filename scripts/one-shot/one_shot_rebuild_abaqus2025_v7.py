# -*- coding: utf-8 -*-
# V7 修订：当前 CAE 中存在同名模型或作业时立即停止，避免覆盖已有对象。
"""
Abaqus/CAE 2025 单颗喷丸冲击参数化建模脚本
================================================

用途
----
1. 创建直径可调的球形弹丸，并将其划分为 8 个可结构化网格的区域；
2. 创建 2 mm x 2 mm x 1 mm 的 TC4 靶材；
3. 将靶材划分为 3 x 3 x 2 共 18 个区域，中央区域细化、外侧区域粗化；
4. 定义 S110 弹丸和 TC4 Johnson-Cook 材料；
5. 创建装配、显式动力学分析步、接触、刚体约束、底面固定和初速度；
6. 创建双精度 Abaqus/Explicit 作业，可选择只写入 INP 或直接提交计算。

单位制
------
长度：mm
质量：tonne
时间：s
力：N
应力：MPa (= N/mm^2)
速度：mm/s
密度：tonne/mm^3

运行方法
--------
Abaqus/CAE -> 文件 -> 运行脚本 -> 选择本文件。

重要说明
--------
- 本脚本会在 RESET_DATABASE=True 时新建模型数据库，当前未保存模型会被清除。
- 建模开始前会检查当前 CAE；若 MODEL_NAME 或 JOB_NAME 已存在，脚本立即停止。
- 默认 SUBMIT_JOB=False，只创建 CAE 和 INP；检查模型后可改为 True 自动提交。
- 分析步真实时长必须写入 ExplicitDynamicsStep 的 timePeriod 参数。
- V5 默认根据弹丸速度、初始间距、弹丸尺寸和材料参数自动估算分析步时间，
  使分析过程覆盖弹丸飞行、撞击接触以及与靶材分离后的安全余量。
- 该时间由碰撞前飞行时间与 Hertz 接触时间估算共同确定；仍建议在 ODB 中检查
  接触力/接触压力是否已降为零以及弹丸速度方向是否已经反向。
"""

from abaqus import *
from abaqusConstants import *
from caeModules import *
import regionToolset
import mesh
import os
import math
import traceback
import abaqus as abaqus_module


# =============================================================================
# 1. 用户参数区：通常只需要修改这里
# =============================================================================

# ---------- 文件、模型和作业 ----------
# True：使用你在启动 Abaqus/CAE 时设置的当前工作目录。
#       脚本启动时通过 os.getcwd() 读取该目录。
# False：忽略 Abaqus/CAE 当前工作目录，改用下面 WORK_DIRECTORY 指定的目录。
USE_CURRENT_WORK_DIRECTORY = True

# 仅当 USE_CURRENT_WORK_DIRECTORY=False 时生效。
# 保留该参数，便于需要时将结果强制保存到指定目录。
WORK_DIRECTORY = r"D:\temp\one-shot"

CAE_FILE_NAME = "one_shot.cae"
MODEL_NAME = "Model-1"
JOB_NAME = "Job-1"
STEP_NAME = "Step-1"

# True：清空当前模型数据库并从头建模；False：在当前数据库中新增模型和作业
# 无论此开关为何值，只要当前 CAE 已有同名模型或作业，脚本都会立即停止。
RESET_DATABASE = False

# 是否写出 Job-1.inp
WRITE_INPUT = True

# 是否在脚本末尾自动提交分析。首次运行建议保持 False，先检查网格和接触。
SUBMIT_JOB = False

# SUBMIT_JOB=True 时，是否等待作业结束后脚本才退出
WAIT_FOR_COMPLETION = True

# ---------- 弹丸几何和运动 ----------
SHOT_RADIUS = 0.15                  # 弹丸半径，mm；直径为 0.30 mm
SHOT_TARGET_GAP = 0.02             # 弹丸最低点与靶材上表面的初始间距，mm
SHOT_VELOCITY_Z = -39030.0         # 沿全局 Z 负方向速度，mm/s；-80000 mm/s = -80 m/s
SHOT_MESH_SIZE = 0.01              # 弹丸网格尺寸，mm

# ---------- 靶材几何 ----------
TARGET_LENGTH = 2.0                 # X 方向尺寸，mm
TARGET_WIDTH = 2.0                  # Y 方向尺寸，mm
TARGET_THICKNESS = 1.0              # Z 方向厚度，mm

# 中央细化区尺寸。默认 x,y 均为 -0.5~0.5 mm，与原操作中的 3x3 分区一致。
TARGET_CENTER_LENGTH = 1.0
TARGET_CENTER_WIDTH = 1.0

# 上部细化层厚度。默认分界面 z=0.7 mm，即上部 0.3 mm 细化层。
TARGET_TOP_FINE_THICKNESS = 0.30

# ---------- 靶材网格 ----------
TARGET_FINE_SIZE = 0.01             # 中央 X/Y 段和上部厚度方向网格尺寸，mm
TARGET_COARSE_SIZE = 0.10           # 外围 X/Y 段粗网格尺寸，mm
TARGET_BOTTOM_MIN_SIZE = 0.01       # 下部厚度方向偏置网格靠近上表面的最小尺寸，mm
TARGET_BOTTOM_MAX_SIZE = 0.10       # 下部厚度方向偏置网格靠近底面的最大尺寸，mm
USE_BIASED_THICKNESS_MESH = True    # True：下部厚度方向由细到粗；失败时自动退回均匀网格
TARGET_BOTTOM_FALLBACK_SIZE = 0.05  # 偏置布种失败后的均匀尺寸，mm

# ---------- 分析步时间 ----------
# True：根据用户填写的 SHOT_VELOCITY_Z 自动计算分析步时间，推荐保持 True。
# False：使用下面的 MANUAL_STEP_TIME。
AUTO_CALCULATE_STEP_TIME = True

# 仅当 AUTO_CALCULATE_STEP_TIME=False 时使用，单位 s。
MANUAL_STEP_TIME = 7.0e-6

# 自动时间计算采用：
#   总时间 = 弹丸飞行时间 + 接触/分离时间窗口
#
# 弹丸飞行时间：
#   t_flight = SHOT_TARGET_GAP / abs(SHOT_VELOCITY_Z)
#
# 接触时间使用球体撞击弹性半空间的 Hertz 近似进行估算。
# Hertz 接触时间本身近似覆盖“首次接触到首次分离”的过程。
# 为考虑 TC4 塑性变形、有限元离散误差和分离后的安全余量，
# 再乘以下面的安全系数。
CONTACT_AND_SEPARATION_TIME_FACTOR = 3.0

# 自动计算时，接触后的时间窗口不得小于该值，单位 s。
# 对当前 0.30 mm 弹丸可保留 1.0e-6 s。
MIN_POST_CONTACT_TIME = 1.0e-6

# 场输出在整个自动计算后的分析步内均匀划分的区间数。
OUTPUT_INTERVALS = 200

# ---------- 接触 ----------
FRICTION_COEFFICIENT = 0.2

# ----------  弹丸材料 ----------
SHOT_DENSITY = 7.85e-9              # tonne/mm^3
SHOT_YOUNG_MODULUS = 210000.0       # MPa
SHOT_POISSON_RATIO = 0.30

# ---------- TC4 靶材材料 ----------
TARGET_DENSITY = 4.50e-9            # tonne/mm^3
TARGET_YOUNG_MODULUS = 120000.0     # MPa
TARGET_POISSON_RATIO = 0.34

# Johnson-Cook 塑性参数：sigma = (A+B*eps_p^n)*(1+C*ln(epsdot/epsdot0))*热软化项
JC_A = 961.0                       # MPa
JC_B = 902.0                       # MPa
JC_N = 0.87
JC_M = 0.0                          # 原 PDF/RPY 设置为 0，即不考虑热软化
JC_MELTING_TEMPERATURE = 0.0
JC_TRANSITION_TEMPERATURE = 0.0
JC_C = 0.01
JC_REFERENCE_STRAIN_RATE = 1.0      # 1/s

# ---------- 作业计算资源 ----------
NUM_CPUS = 16                        # 按电脑实际核心数修改
NUM_DOMAINS = 16                     # DOMAIN 并行时应为 NUM_CPUS 的整数倍
MEMORY_PERCENT = 90

# 使用“打包器+分析器双精度”，避免显式分析因增量很多而被单精度舍入误差终止。
EXPLICIT_PRECISION = DOUBLE_PLUS_PACK
NODAL_OUTPUT_PRECISION = FULL


# =============================================================================
# 2. 通用辅助函数：一般不需要修改
# =============================================================================

def _print_title(text):
    print("\n" + "=" * 78)
    print(text)
    print("=" * 78)


def _estimate_analysis_time():
    """
    根据弹丸速度自动估算显式分析步时间。

    组成
    ----
    1. 飞行时间：
           t_flight = gap / velocity

    2. Hertz 接触时间：
           t_contact = 2.94 * [m^2 / (R * E_eff^2 * velocity)]^(1/5)

       其中：
           m      为球形弹丸质量；
           R      为弹丸半径；
           E_eff  为弹丸和靶材的等效接触模量。

       Hertz 接触时间近似表示从首次接触到首次分离的持续时间。
       当前模型存在 TC4 塑性变形，因此脚本再乘安全系数，并设置
       最小接触后时间窗口。

    返回
    ----
    step_time, flight_time, hertz_contact_time, post_contact_time,
    shot_mass, effective_modulus
    """
    speed = abs(float(SHOT_VELOCITY_Z))
    flight_time = float(SHOT_TARGET_GAP) / speed

    # 球形弹丸质量。当前一致单位制为 mm、N、s、tonne。
    shot_mass = (
        4.0
        / 3.0
        * math.pi
        * float(SHOT_RADIUS) ** 3
        * float(SHOT_DENSITY)
    )

    # 两种材料的等效接触模量：
    # 1/E_eff = (1-v1^2)/E1 + (1-v2^2)/E2
    compliance = (
        (1.0 - float(SHOT_POISSON_RATIO) ** 2)
        / float(SHOT_YOUNG_MODULUS)
        + (1.0 - float(TARGET_POISSON_RATIO) ** 2)
        / float(TARGET_YOUNG_MODULUS)
    )
    effective_modulus = 1.0 / compliance

    hertz_contact_time = 2.94 * (
        shot_mass ** 2
        / (
            float(SHOT_RADIUS)
            * effective_modulus ** 2
            * speed
        )
    ) ** 0.2

    estimated_post_contact_time = (
        float(CONTACT_AND_SEPARATION_TIME_FACTOR)
        * hertz_contact_time
    )
    post_contact_time = max(
        estimated_post_contact_time,
        float(MIN_POST_CONTACT_TIME),
    )

    automatic_step_time = flight_time + post_contact_time

    return (
        automatic_step_time,
        flight_time,
        hertz_contact_time,
        post_contact_time,
        shot_mass,
        effective_modulus,
    )


def _validate_parameters():
    """
    在正式建模前检查参数，并返回本次实际采用的分析步时间。

    当 AUTO_CALCULATE_STEP_TIME=True 时，修改 SHOT_VELOCITY_Z 后，
    分析步时间会自动更新，不再需要手动同步修改固定 STEP_TIME。
    """
    positive_values = {
        "SHOT_RADIUS": SHOT_RADIUS,
        "SHOT_TARGET_GAP": SHOT_TARGET_GAP,
        "SHOT_MESH_SIZE": SHOT_MESH_SIZE,
        "SHOT_DENSITY": SHOT_DENSITY,
        "SHOT_YOUNG_MODULUS": SHOT_YOUNG_MODULUS,
        "TARGET_LENGTH": TARGET_LENGTH,
        "TARGET_WIDTH": TARGET_WIDTH,
        "TARGET_THICKNESS": TARGET_THICKNESS,
        "TARGET_CENTER_LENGTH": TARGET_CENTER_LENGTH,
        "TARGET_CENTER_WIDTH": TARGET_CENTER_WIDTH,
        "TARGET_TOP_FINE_THICKNESS": TARGET_TOP_FINE_THICKNESS,
        "TARGET_FINE_SIZE": TARGET_FINE_SIZE,
        "TARGET_COARSE_SIZE": TARGET_COARSE_SIZE,
        "TARGET_DENSITY": TARGET_DENSITY,
        "TARGET_YOUNG_MODULUS": TARGET_YOUNG_MODULUS,
        "CONTACT_AND_SEPARATION_TIME_FACTOR":
            CONTACT_AND_SEPARATION_TIME_FACTOR,
        "OUTPUT_INTERVALS": OUTPUT_INTERVALS,
    }

    for name, value in positive_values.items():
        if value <= 0.0:
            raise ValueError(
                "参数 %s 必须大于 0，当前值为 %s"
                % (name, value)
            )

    if MIN_POST_CONTACT_TIME < 0.0:
        raise ValueError("MIN_POST_CONTACT_TIME 不能小于 0。")

    if SHOT_VELOCITY_Z >= 0.0:
        raise ValueError(
            "SHOT_VELOCITY_Z 应为负值，使弹丸沿全局 Z 负方向撞击靶材。"
        )

    if not (-1.0 < SHOT_POISSON_RATIO < 0.5):
        raise ValueError("SHOT_POISSON_RATIO 必须位于 -1 和 0.5 之间。")
    if not (-1.0 < TARGET_POISSON_RATIO < 0.5):
        raise ValueError("TARGET_POISSON_RATIO 必须位于 -1 和 0.5 之间。")

    if TARGET_CENTER_LENGTH >= TARGET_LENGTH:
        raise ValueError("TARGET_CENTER_LENGTH 必须小于 TARGET_LENGTH。")
    if TARGET_CENTER_WIDTH >= TARGET_WIDTH:
        raise ValueError("TARGET_CENTER_WIDTH 必须小于 TARGET_WIDTH。")
    if TARGET_TOP_FINE_THICKNESS >= TARGET_THICKNESS:
        raise ValueError(
            "TARGET_TOP_FINE_THICKNESS 必须小于 TARGET_THICKNESS。"
        )

    if NUM_CPUS < 1 or NUM_DOMAINS < 1:
        raise ValueError("NUM_CPUS 和 NUM_DOMAINS 必须为正整数。")
    if NUM_DOMAINS % NUM_CPUS != 0:
        raise ValueError(
            "采用 DOMAIN 并行时，NUM_DOMAINS 必须是 NUM_CPUS 的整数倍。"
        )

    (
        automatic_step_time,
        flight_time,
        hertz_contact_time,
        post_contact_time,
        shot_mass,
        effective_modulus,
    ) = _estimate_analysis_time()

    if AUTO_CALCULATE_STEP_TIME:
        actual_step_time = automatic_step_time
        time_source = "根据弹丸速度自动计算"
    else:
        if MANUAL_STEP_TIME <= 0.0:
            raise ValueError("MANUAL_STEP_TIME 必须大于 0。")
        actual_step_time = float(MANUAL_STEP_TIME)
        time_source = "使用 MANUAL_STEP_TIME"

        if actual_step_time <= flight_time:
            raise ValueError(
                "MANUAL_STEP_TIME=%.6e s 不足以使弹丸到达靶材；"
                "飞行时间约为 %.6e s。"
                % (actual_step_time, flight_time)
            )

        if actual_step_time < automatic_step_time:
            print(
                "警告：手动分析步时间 %.6e s 小于自动建议值 %.6e s，"
                "可能无法完整覆盖弹丸接触和分离过程。"
                % (actual_step_time, automatic_step_time)
            )

    print("分析步时间来源：%s" % time_source)
    print("弹丸速度绝对值：%.6e mm/s" % abs(SHOT_VELOCITY_Z))
    print("弹丸质量：%.6e tonne" % shot_mass)
    print("等效接触模量：%.6e MPa" % effective_modulus)
    print("弹丸到达靶材前的飞行时间：%.6e s" % flight_time)
    print("Hertz 估算首次接触至首次分离时间：%.6e s" % hertz_contact_time)
    print("采用的接触/分离时间窗口：%.6e s" % post_contact_time)
    print("最终分析步时间：%.6e s" % actual_step_time)

    return actual_step_time


def _resolve_work_directory():
    """
    确定本次脚本实际使用的工作目录。

    USE_CURRENT_WORK_DIRECTORY=True：
        使用 Abaqus/CAE 当前工作目录，即脚本启动时的 os.getcwd()。
    USE_CURRENT_WORK_DIRECTORY=False：
        使用用户参数区中的 WORK_DIRECTORY。
    """
    if USE_CURRENT_WORK_DIRECTORY:
        work_directory = os.getcwd()
        source_text = "Abaqus/CAE 当前工作目录"
    else:
        work_directory = WORK_DIRECTORY
        source_text = "脚本参数 WORK_DIRECTORY"

    if work_directory is None or str(work_directory).strip() == "":
        raise ValueError("工作目录为空，请检查 USE_CURRENT_WORK_DIRECTORY 和 WORK_DIRECTORY。")

    work_directory = os.path.abspath(str(work_directory))
    if not os.path.isdir(work_directory):
        os.makedirs(work_directory)

    print("本次使用的工作目录：%s" % work_directory)
    print("工作目录来源：%s" % source_text)
    return work_directory


def _check_existing_model_and_job_names():
    """若当前 CAE 中已存在同名模型或作业，则在任何建模操作前停止脚本。"""
    conflicts = []

    if MODEL_NAME in list(mdb.models.keys()):
        conflicts.append("模型 %s" % MODEL_NAME)
    if JOB_NAME in list(mdb.jobs.keys()):
        conflicts.append("作业 %s" % JOB_NAME)

    if conflicts:
        raise RuntimeError(
            "当前 CAE 中已存在同名%s。脚本已停止，未删除或覆盖已有模型和作业。"
            % "、".join(conflicts)
        )

    print("名称检查通过：当前 CAE 中不存在模型 %s 或作业 %s。" % (MODEL_NAME, JOB_NAME))


def _prepare_database():
    """
    创建干净的模型数据库，避免旧特征、旧作业和同名对象干扰复现。

    关键修正：
    Mdb() 会创建新的当前模型数据库。若脚本仍继续使用旧的 mdb Python
    引用，可能出现模型树中已经显示 Model-1，但访问 mdb.models['Model-1']
    或删除同名模型时仍触发 KeyError 的情况。因此这里显式重新绑定全局 mdb。
    """
    global mdb

    if RESET_DATABASE:
        new_database = Mdb()

        # Abaqus 版本/执行环境不同，Mdb() 通常返回新 ModelDatabase；
        # 若返回 None，则从 abaqus 模块重新取得当前 mdb。
        if new_database is not None:
            mdb = new_database
        else:
            mdb = abaqus_module.mdb

        model_names = list(mdb.models.keys())

        # 不假定新数据库中一定存在 Model-1，避免再次触发 KeyError。
        if MODEL_NAME in model_names:
            model = mdb.models[MODEL_NAME]
        elif "Model-1" in model_names:
            if MODEL_NAME != "Model-1":
                mdb.models.changeKey(fromName="Model-1", toName=MODEL_NAME)
            model = mdb.models[MODEL_NAME]
        elif len(model_names) > 0:
            # 极少数环境可能使用不同的默认模型名，直接重命名第一个模型。
            first_model_name = model_names[0]
            if first_model_name != MODEL_NAME:
                mdb.models.changeKey(fromName=first_model_name, toName=MODEL_NAME)
            model = mdb.models[MODEL_NAME]
        else:
            model = mdb.Model(name=MODEL_NAME)

    else:
        # 同名检查已在 main() 开始处完成，此处只新增模型，不删除已有对象。
        model = mdb.Model(name=MODEL_NAME)

    print("当前模型数据库中的模型：%s" % list(mdb.models.keys()))
    return model


def _partition_by_principal_plane(part, principal_plane, offset):
    """使用基准面切分当前部件的全部实体单元，避免 RPY 中不稳定的 mask 索引。"""
    datum_feature = part.DatumPlaneByPrincipalPlane(
        principalPlane=principal_plane,
        offset=offset,
    )
    part.PartitionCellByDatumPlane(
        datumPlane=part.datums[datum_feature.id],
        cells=part.cells[:],
    )


def _create_shot(model):
    """创建球形弹丸、中心参考点以及三个正交切分面。"""
    sketch_name = "__shot_profile__"
    sketch = model.ConstrainedSketch(
        name=sketch_name,
        sheetSize=max(10.0 * SHOT_RADIUS, 1.0),
    )

    # 绕 Y 轴旋转半圆形成球体。闭合直线位于旋转轴上。
    sketch.ConstructionLine(
        point1=(0.0, -2.0 * SHOT_RADIUS),
        point2=(0.0, 2.0 * SHOT_RADIUS),
    )
    sketch.ArcByCenterEnds(
        center=(0.0, 0.0),
        point1=(0.0, SHOT_RADIUS),
        point2=(0.0, -SHOT_RADIUS),
        direction=CLOCKWISE,
    )
    sketch.Line(
        point1=(0.0, SHOT_RADIUS),
        point2=(0.0, -SHOT_RADIUS),
    )

    shot = model.Part(
        name="shot",
        dimensionality=THREE_D,
        type=DEFORMABLE_BODY,
    )
    shot.BaseSolidRevolve(
        sketch=sketch,
        angle=360.0,
        flipRevolveDirection=OFF,
    )
    del model.sketches[sketch_name]

    # 参考点位于球心，后续用于刚体约束和初速度。
    rp_feature = shot.ReferencePoint(point=(0.0, 0.0, 0.0))
    rp_id = rp_feature.id
    shot.Set(
        name="SET_SHOT_RP",
        referencePoints=(shot.referencePoints[rp_id],),
    )

    # 三个相互正交的中心平面将球体切成 8 个区域，便于六面体结构网格。
    _partition_by_principal_plane(shot, XYPLANE, 0.0)
    _partition_by_principal_plane(shot, XZPLANE, 0.0)
    _partition_by_principal_plane(shot, YZPLANE, 0.0)

    return shot, rp_id


def _create_target(model):
    """创建靶材，并用 x/y 两对分区面和一个厚度分区面切成 18 个区域。"""
    sketch_name = "__target_profile__"
    half_l = TARGET_LENGTH / 2.0
    half_w = TARGET_WIDTH / 2.0

    sketch = model.ConstrainedSketch(
        name=sketch_name,
        sheetSize=2.5 * max(TARGET_LENGTH, TARGET_WIDTH),
    )
    sketch.rectangle(
        point1=(-half_l, -half_w),
        point2=(half_l, half_w),
    )

    target = model.Part(
        name="target",
        dimensionality=THREE_D,
        type=DEFORMABLE_BODY,
    )
    target.BaseSolidExtrude(sketch=sketch, depth=TARGET_THICKNESS)
    del model.sketches[sketch_name]

    center_half_x = TARGET_CENTER_LENGTH / 2.0
    center_half_y = TARGET_CENTER_WIDTH / 2.0
    z_interface = TARGET_THICKNESS - TARGET_TOP_FINE_THICKNESS

    # YZPLANE 的偏置控制 x 坐标；XZPLANE 的偏置控制 y 坐标。
    _partition_by_principal_plane(target, YZPLANE, -center_half_x)
    _partition_by_principal_plane(target, YZPLANE, center_half_x)
    _partition_by_principal_plane(target, XZPLANE, -center_half_y)
    _partition_by_principal_plane(target, XZPLANE, center_half_y)

    # XYPLANE 偏置形成上部细化层。
    _partition_by_principal_plane(target, XYPLANE, z_interface)

    return target


def _assign_materials_and_sections(model, shot, target):
    """定义材料、实体截面，并对两个部件的全部实体区域赋予截面。"""
    shot_mat_name = "Material-shot-S110"
    shot_sec_name = "Section-shot-S110"
    target_mat_name = "Material-target-TC4"
    target_sec_name = "Section-target-TC4"

    shot_material = model.Material(name=shot_mat_name)
    shot_material.Density(table=((SHOT_DENSITY,),))
    shot_material.Elastic(
        table=((SHOT_YOUNG_MODULUS, SHOT_POISSON_RATIO),)
    )

    target_material = model.Material(name=target_mat_name)
    target_material.Density(table=((TARGET_DENSITY,),))
    target_material.Elastic(
        table=((TARGET_YOUNG_MODULUS, TARGET_POISSON_RATIO),)
    )
    target_material.Plastic(
        hardening=JOHNSON_COOK,
        scaleStress=None,
        table=((
            JC_A,
            JC_B,
            JC_N,
            JC_M,
            JC_MELTING_TEMPERATURE,
            JC_TRANSITION_TEMPERATURE,
        ),),
    )
    target_material.plastic.RateDependent(
        type=JOHNSON_COOK,
        table=((JC_C, JC_REFERENCE_STRAIN_RATE),),
    )

    model.HomogeneousSolidSection(
        name=shot_sec_name,
        material=shot_mat_name,
        thickness=None,
    )
    model.HomogeneousSolidSection(
        name=target_sec_name,
        material=target_mat_name,
        thickness=None,
    )

    shot.SectionAssignment(
        region=regionToolset.Region(cells=shot.cells[:]),
        sectionName=shot_sec_name,
        offset=0.0,
        offsetType=MIDDLE_SURFACE,
        offsetField="",
        thicknessAssignment=FROM_SECTION,
    )
    target.SectionAssignment(
        region=regionToolset.Region(cells=target.cells[:]),
        sectionName=target_sec_name,
        offset=0.0,
        offsetType=MIDDLE_SURFACE,
        offsetField="",
        thicknessAssignment=FROM_SECTION,
    )


def _edge_coordinates(part, edge):
    """返回直线边两个端点坐标；若无法取得两个端点则返回 None。"""
    vertex_ids = edge.getVertices()
    if len(vertex_ids) != 2:
        return None
    p1 = part.vertices[vertex_ids[0]].pointOn[0]
    p2 = part.vertices[vertex_ids[1]].pointOn[0]
    return p1, p2


def _seed_target_edges(target):
    """
    按几何坐标自动识别 X/Y/Z 方向边并布种。

    这种方法替代 RPY 中 getSequenceFromMask(...) 的索引选择，因此修改尺寸后仍较稳定。
    X/Y 中央段使用细网格，外围段使用粗网格；厚度上层细化，下层偏置。
    """
    center_half_x = TARGET_CENTER_LENGTH / 2.0
    center_half_y = TARGET_CENTER_WIDTH / 2.0
    z_interface = TARGET_THICKNESS - TARGET_TOP_FINE_THICKNESS
    tol = 1.0e-7 * max(TARGET_LENGTH, TARGET_WIDTH, TARGET_THICKNESS, 1.0)

    x_fine, x_coarse = [], []
    y_fine, y_coarse = [], []
    z_top = []
    z_bottom_end1 = []
    z_bottom_end2 = []
    z_bottom_all = []

    for edge in target.edges:
        coords = _edge_coordinates(target, edge)
        if coords is None:
            continue

        p1, p2 = coords
        dx = abs(p2[0] - p1[0])
        dy = abs(p2[1] - p1[1])
        dz = abs(p2[2] - p1[2])

        # X 方向边：根据该边是否处于中央 x 段决定细/粗网格。
        if dx > tol and dy <= tol and dz <= tol:
            xmin = min(p1[0], p2[0])
            xmax = max(p1[0], p2[0])
            if xmin >= -center_half_x - tol and xmax <= center_half_x + tol:
                x_fine.append(edge)
            else:
                x_coarse.append(edge)

        # Y 方向边：根据该边是否处于中央 y 段决定细/粗网格。
        elif dy > tol and dx <= tol and dz <= tol:
            ymin = min(p1[1], p2[1])
            ymax = max(p1[1], p2[1])
            if ymin >= -center_half_y - tol and ymax <= center_half_y + tol:
                y_fine.append(edge)
            else:
                y_coarse.append(edge)

        # Z 方向边：上部细化层均匀细化；下部从上向下逐渐变粗。
        elif dz > tol and dx <= tol and dy <= tol:
            zmin = min(p1[2], p2[2])
            zmax = max(p1[2], p2[2])

            if zmin >= z_interface - tol:
                z_top.append(edge)
            elif zmax <= z_interface + tol:
                z_bottom_all.append(edge)

                # seedEdgeByBias 需要区分边的内部方向，使最小尺寸位于 z_interface 一侧。
                # getVertices() 返回顺序对应 edge 的 end1/end2。
                if p1[2] > p2[2]:
                    z_bottom_end1.append(edge)
                else:
                    z_bottom_end2.append(edge)

    # 先给整个部件粗网格，再对中央段和顶部层覆盖为细网格。
    target.seedPart(
        size=TARGET_COARSE_SIZE,
        deviationFactor=0.1,
        minSizeFactor=0.1,
    )

    if x_coarse:
        target.seedEdgeBySize(
            edges=tuple(x_coarse),
            size=TARGET_COARSE_SIZE,
            deviationFactor=0.1,
            minSizeFactor=0.1,
            constraint=FINER,
        )
    if y_coarse:
        target.seedEdgeBySize(
            edges=tuple(y_coarse),
            size=TARGET_COARSE_SIZE,
            deviationFactor=0.1,
            minSizeFactor=0.1,
            constraint=FINER,
        )
    if x_fine:
        target.seedEdgeBySize(
            edges=tuple(x_fine),
            size=TARGET_FINE_SIZE,
            deviationFactor=0.1,
            minSizeFactor=0.1,
            constraint=FINER,
        )
    if y_fine:
        target.seedEdgeBySize(
            edges=tuple(y_fine),
            size=TARGET_FINE_SIZE,
            deviationFactor=0.1,
            minSizeFactor=0.1,
            constraint=FINER,
        )
    if z_top:
        target.seedEdgeBySize(
            edges=tuple(z_top),
            size=TARGET_FINE_SIZE,
            deviationFactor=0.1,
            minSizeFactor=0.1,
            constraint=FINER,
        )

    if z_bottom_all:
        if USE_BIASED_THICKNESS_MESH:
            try:
                kwargs = {
                    "biasMethod": SINGLE,
                    "minSize": TARGET_BOTTOM_MIN_SIZE,
                    "maxSize": TARGET_BOTTOM_MAX_SIZE,
                    "constraint": FINER,
                }
                if z_bottom_end1:
                    kwargs["end1Edges"] = tuple(z_bottom_end1)
                if z_bottom_end2:
                    kwargs["end2Edges"] = tuple(z_bottom_end2)
                target.seedEdgeByBias(**kwargs)
                print("靶材下部厚度方向：已采用从上部细到下部粗的偏置网格。")
            except Exception:
                print("警告：厚度方向偏置布种失败，自动改用均匀网格。")
                target.seedEdgeBySize(
                    edges=tuple(z_bottom_all),
                    size=TARGET_BOTTOM_FALLBACK_SIZE,
                    deviationFactor=0.1,
                    minSizeFactor=0.1,
                    constraint=FINER,
                )
        else:
            target.seedEdgeBySize(
                edges=tuple(z_bottom_all),
                size=TARGET_BOTTOM_FALLBACK_SIZE,
                deviationFactor=0.1,
                minSizeFactor=0.1,
                constraint=FINER,
            )

    print(
        "靶材边分类：X细=%d, X粗=%d, Y细=%d, Y粗=%d, Z上层=%d, Z下层=%d"
        % (
            len(x_fine), len(x_coarse),
            len(y_fine), len(y_coarse),
            len(z_top), len(z_bottom_all),
        )
    )


def _mesh_parts(shot, target):
    """指定 C3D8R 显式单元、结构化六面体网格并生成网格。"""
    explicit_hex = mesh.ElemType(
        elemCode=C3D8R,
        elemLibrary=EXPLICIT,
        kinematicSplit=AVERAGE_STRAIN,
        hourglassControl=DEFAULT,
        distortionControl=DEFAULT,
    )

    # 弹丸：8 分区结构化六面体网格。
    shot.setMeshControls(
        regions=shot.cells[:],
        elemShape=HEX,
        technique=STRUCTURED,
    )
    shot.setElementType(
        regions=(shot.cells[:],),
        elemTypes=(explicit_hex,),
    )
    shot.seedPart(
        size=SHOT_MESH_SIZE,
        deviationFactor=0.1,
        minSizeFactor=0.1,
    )
    shot.generateMesh()

    # 靶材：18 个块体均使用结构化六面体网格。
    target.setMeshControls(
        regions=target.cells[:],
        elemShape=HEX,
        technique=STRUCTURED,
    )
    target.setElementType(
        regions=(target.cells[:],),
        elemTypes=(explicit_hex,),
    )
    _seed_target_edges(target)
    target.generateMesh()

    print("弹丸网格：节点 %d，单元 %d" % (len(shot.nodes), len(shot.elements)))
    print("靶材网格：节点 %d，单元 %d" % (len(target.nodes), len(target.elements)))


def _get_geom_sequence_from_indices(entity_array, indices, description):
    """将普通 Python 索引列表转换为 Abaqus GeomSequence。

    Abaqus 的 ``Region(side1Faces=...)``、``Set(faces=...)`` 等接口不接受
    Python 的 list/tuple，而要求 FaceArray（错误信息中称为 GeomSequence）。
    ``getSequenceFromMask`` 返回的正是 Abaqus 原生几何序列，因此这里根据
    实体索引自动生成压缩掩码，避免继续把 tuple 传给 Abaqus API。

    参数
    ----
    entity_array : FaceArray / EdgeArray / CellArray
        Abaqus 几何数组，例如 ``instance.faces``。
    indices : sequence of int
        需要选取的实体索引。
    description : str
        用于异常信息的区域说明。
    """
    unique_indices = sorted(set(int(index) for index in indices))
    if not unique_indices:
        raise RuntimeError("%s：没有可转换的几何实体索引。" % description)

    # Abaqus 压缩掩码每个十六进制字段表示 32 个连续实体。
    word_count = unique_indices[-1] // 32 + 1
    words = [0] * word_count
    for index in unique_indices:
        words[index // 32] |= 1 << (index % 32)

    mask_text = "[" + " ".join("#%x" % word for word in words) + " ]"
    sequence = entity_array.getSequenceFromMask(mask=(mask_text,))

    if len(sequence) != len(unique_indices):
        raise RuntimeError(
            "%s：GeomSequence 转换数量不一致，期望 %d，实际 %d；mask=%s"
            % (description, len(unique_indices), len(sequence), mask_text)
        )
    return sequence


def _select_shot_outer_faces(shot_instance, center):
    """按球心距离选出球形外表面，并返回 Abaqus FaceArray。

    球形外表面上任意内部点到球心的距离都应接近 SHOT_RADIUS；
    三个内部切分平面上的 pointOn 点到球心距离明显更小。
    这里采用 0.90R 的阈值，可排除球体内部的三个切分平面。
    """
    selected_indices = []
    threshold = 0.90 * SHOT_RADIUS
    for face in shot_instance.faces:
        point = face.pointOn[0]
        radius = math.sqrt(
            (point[0] - center[0]) ** 2
            + (point[1] - center[1]) ** 2
            + (point[2] - center[2]) ** 2
        )
        if radius >= threshold:
            selected_indices.append(face.index)

    if not selected_indices:
        raise RuntimeError("未识别到弹丸外表面，请检查 SHOT_RADIUS 或几何分区。")

    return _get_geom_sequence_from_indices(
        shot_instance.faces,
        selected_indices,
        "弹丸外表面",
    )


def _select_faces_at_z(instance, z_value, tolerance):
    """按 Z 坐标选取水平表面，并返回 Abaqus FaceArray。"""
    selected_indices = []
    for face in instance.faces:
        point = face.pointOn[0]
        if abs(point[2] - z_value) <= tolerance:
            selected_indices.append(face.index)

    if not selected_indices:
        raise RuntimeError("在 z=%.6g 附近未找到目标表面。" % z_value)

    return _get_geom_sequence_from_indices(
        instance.faces,
        selected_indices,
        "z=%.6g 的水平表面" % z_value,
    )


def _create_assembly_step_contact_and_loads(
    model,
    shot,
    target,
    shot_rp_id,
    actual_step_time,
):
    """创建装配、分析步、接触、刚体约束、固定边界和弹丸初速度。"""
    assembly = model.rootAssembly
    assembly.DatumCsysByDefault(CARTESIAN)

    shot_instance = assembly.Instance(
        name="shot-1",
        part=shot,
        dependent=ON,
    )
    target_instance = assembly.Instance(
        name="target-1",
        part=target,
        dependent=ON,
    )

    # 靶材上表面位于 z=TARGET_THICKNESS。
    # 弹丸中心高度 = 靶材厚度 + 弹丸半径 + 初始间距。
    shot_center_z = TARGET_THICKNESS + SHOT_RADIUS + SHOT_TARGET_GAP
    assembly.translate(
        instanceList=("shot-1",),
        vector=(0.0, 0.0, shot_center_z),
    )
    assembly.regenerate()

    # 显式动力学分析步。
    # actual_step_time 已根据弹丸速度自动计算，覆盖飞行、接触和分离过程。
    model.ExplicitDynamicsStep(
        name=STEP_NAME,
        previous="Initial",
        timePeriod=actual_step_time,
        improvedDtMethod=ON,
    )

    # 默认场输出只调整输出间隔，保留 Abaqus/Explicit 的预选变量。
    model.fieldOutputRequests["F-Output-1"].setValues(
        numIntervals=OUTPUT_INTERVALS,
    )

    # ---------------------- 接触属性 ----------------------
    contact_property_name = "IntProp-1"
    model.ContactProperty(contact_property_name)
    model.interactionProperties[contact_property_name].TangentialBehavior(
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
        elasticSlipStiffness=None,
    )
    model.interactionProperties[contact_property_name].NormalBehavior(
        pressureOverclosure=HARD,
        allowSeparation=ON,
        constraintEnforcementMethod=DEFAULT,
    )

    tol = 1.0e-6 * max(TARGET_LENGTH, TARGET_WIDTH, TARGET_THICKNESS, 1.0)
    shot_center = (0.0, 0.0, shot_center_z)
    shot_outer_faces = _select_shot_outer_faces(shot_instance, shot_center)
    target_top_faces = _select_faces_at_z(
        target_instance,
        TARGET_THICKNESS,
        tol,
    )

    # 直接创建 Region，不额外创建装配级 Surface 特征。
    # 注意：side1Faces 必须是 Abaqus FaceArray/GeomSequence，不能是 Python tuple。
    # 上面的选择函数已经通过 getSequenceFromMask 完成类型转换。
    print("识别到弹丸外表面数量：%d" % len(shot_outer_faces))
    print("识别到靶材上表面数量：%d" % len(target_top_faces))

    shot_contact_region = regionToolset.Region(
        side1Faces=shot_outer_faces,
    )
    target_contact_region = regionToolset.Region(
        side1Faces=target_top_faces,
    )

    model.SurfaceToSurfaceContactExp(
        name="Int-1",
        createStepName="Initial",
        main=shot_contact_region,
        secondary=target_contact_region,
        mechanicalConstraint=KINEMATIC,
        sliding=FINITE,
        interactionProperty=contact_property_name,
        initialClearance=OMIT,
        datumAxis=None,
        clearanceRegion=None,
    )

    # ---------------------- 弹丸刚体约束 ----------------------
    shot_rp = shot_instance.referencePoints[shot_rp_id]
    assembly.Set(
        name="SET_SHOT_RP",
        referencePoints=(shot_rp,),
    )
    assembly.Set(
        name="SET_SHOT_BODY",
        cells=shot_instance.cells[:],
    )

    model.RigidBody(
        name="Constraint-shot-rigid",
        refPointRegion=assembly.sets["SET_SHOT_RP"],
        bodyRegion=assembly.sets["SET_SHOT_BODY"],
    )

    # ---------------------- 靶材底面固定 ----------------------
    target_bottom_faces = _select_faces_at_z(target_instance, 0.0, tol)
    assembly.Set(
        name="SET_TARGET_BOTTOM",
        faces=target_bottom_faces,
    )
    model.EncastreBC(
        name="BC-target-bottom",
        createStepName="Initial",
        region=assembly.sets["SET_TARGET_BOTTOM"],
        localCsys=None,
    )

    # ---------------------- 弹丸初速度 ----------------------
    model.Velocity(
        name="Predefined-Field-shot-velocity",
        region=assembly.sets["SET_SHOT_RP"],
        field="",
        distributionType=MAGNITUDE,
        velocity1=0.0,
        velocity2=0.0,
        velocity3=SHOT_VELOCITY_Z,
        omega=0.0,
    )

    return assembly


def _create_job_and_save(model, work_directory):
    """创建双精度显式作业，保存 CAE，按开关写入 INP 或提交。"""
    # 切换到已经解析并验证过的工作目录，使 INP、ODB、STA 等文件均写入此处。
    os.chdir(work_directory)

    job = mdb.Job(
        name=JOB_NAME,
        model=MODEL_NAME,
        description="Single-shot peening impact model",
        type=ANALYSIS,
        atTime=None,
        waitMinutes=0,
        waitHours=0,
        queue=None,
        memory=MEMORY_PERCENT,
        memoryUnits=PERCENTAGE,
        getMemoryFromAnalysis=True,
        explicitPrecision=EXPLICIT_PRECISION,
        nodalOutputPrecision=NODAL_OUTPUT_PRECISION,
        echoPrint=OFF,
        modelPrint=OFF,
        contactPrint=OFF,
        historyPrint=OFF,
        userSubroutine="",
        scratch="",
        resultsFormat=ODB,
        multiprocessingMode=DEFAULT,
        numCpus=NUM_CPUS,
        numDomains=NUM_DOMAINS,
        activateLoadBalancing=False,
    )

    cae_path = os.path.join(work_directory, CAE_FILE_NAME)
    mdb.saveAs(pathName=cae_path)
    print("CAE 文件已保存：%s" % cae_path)

    if WRITE_INPUT:
        job.writeInput(consistencyChecking=OFF)
        print("INP 文件已写出：%s" % os.path.join(work_directory, JOB_NAME + ".inp"))

    if SUBMIT_JOB:
        print("正在提交作业 %s ..." % JOB_NAME)
        job.submit(consistencyChecking=OFF)
        if WAIT_FOR_COMPLETION:
            job.waitForCompletion()
            print("作业最终状态：%s" % job.status)
            mdb.save()
    else:
        print("SUBMIT_JOB=False：模型和 INP 已创建，但未自动提交计算。")


def main():
    _print_title("Abaqus/Explicit 单颗喷丸冲击参数化建模")

    # 必须首先检查当前 CAE，确保脚本不会删除或覆盖同名模型和作业。
    _check_existing_model_and_job_names()

    # 必须在 Mdb() 新建数据库之前读取当前目录，确保采用用户在 CAE 中设置的工作目录。
    work_directory = _resolve_work_directory()
    actual_step_time = _validate_parameters()

    _print_title("1/6 初始化模型数据库")
    model = _prepare_database()

    _print_title("2/6 创建弹丸和靶材几何及分区")
    shot, shot_rp_id = _create_shot(model)
    target = _create_target(model)
    print("弹丸实体分区数：%d（目标为 8）" % len(shot.cells))
    print("靶材实体分区数：%d（目标为 18）" % len(target.cells))

    _print_title("3/6 定义材料和截面")
    _assign_materials_and_sections(model, shot, target)

    _print_title("4/6 划分结构化六面体网格")
    _mesh_parts(shot, target)

    _print_title("5/6 创建装配、分析步、接触和载荷")
    _create_assembly_step_contact_and_loads(
        model,
        shot,
        target,
        shot_rp_id,
        actual_step_time,
    )

    _print_title("6/6 创建作业、保存并按设置写入/提交")
    _create_job_and_save(model, work_directory)

    _print_title("脚本执行完成")
    print("模型：%s" % MODEL_NAME)
    print("作业：%s" % JOB_NAME)
    print("分析步时间：%.6e s" % actual_step_time)
    print("显式精度：DOUBLE_PLUS_PACK")
    print("结果目录：%s" % work_directory)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print("\n脚本执行失败：%s" % str(error))
        traceback.print_exc()
        raise
