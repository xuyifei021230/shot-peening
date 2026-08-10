# -*- coding: utf-8 -*-
"""
Abaqus/CAE 2025 单颗喷丸冲击后的静态分析参数化脚本
======================================================

脚本用途
--------
本脚本以已经完成计算的 one-shot 显式冲击模型为基础，自动建立冲击后的
Abaqus/Standard 静态模型，用于计算靶材在保留冲击残余应力、塑性应变等
状态变量后的静力平衡或回弹状态。

主要流程
--------
1. 打开已有 one-shot.cae；
2. 将源模型复制为静态模型，不修改原始显式模型；
3. 抑制弹丸实例；
4. 抑制弹丸—靶材接触、弹丸刚体约束和弹丸初速度；
5. 从 Job-1.odb 的最后一步、最后增量导入 target-1 的初始状态；
6. 删除复制模型中的显式动力学分析步；
7. 创建 Static, General 静力分析步，并启用几何非线性；
8. 创建 Abaqus/Standard 作业；
9. 保存新的 CAE，并按设置写出 INP 或直接提交计算。

单位制
------
本脚本不重新定义材料或几何，因此单位制必须与源 one-shot 模型完全一致。
原 one-shot 模型采用：mm、N、s、MPa、tonne。

运行方法
--------
Abaqus/CAE 2025 -> 文件 -> 运行脚本 -> 选择本文件。

重要说明
--------
- 默认使用 Abaqus/CAE 当前工作目录；也可以在参数区指定固定目录。
- 默认根据源 CAE 文件名自动生成静态 CAE 名称，例如 one-shot.cae -> one-shot-static.cae。
- 默认只创建模型并写出 Job-1-Static.inp，不自动提交。
- 默认自动令域数等于 CPU 核心数，避免并行作业参数不匹配。
- Job-1.odb 必须来自与 one-shot.cae 中源模型网格完全一致的显式分析。
"""

from abaqus import *
from abaqusConstants import *
from caeModules import *
import os
import traceback
import abaqus as abaqus_module


# =============================================================================
# 1. 用户参数区：通常只需要修改这里
# =============================================================================

# -----------------------------------------------------------------------------
# 1.1 工作目录和输入文件
# -----------------------------------------------------------------------------

# True：使用启动 Abaqus/CAE 后设置的当前工作目录，即 os.getcwd()。
# False：使用下面 WORK_DIRECTORY 指定的目录。
USE_CURRENT_WORK_DIRECTORY = True

# 仅当 USE_CURRENT_WORK_DIRECTORY=False 时生效。
WORK_DIRECTORY = r"D:\temp\one-shot"

# True：脚本主动打开下面的源 CAE 文件。
# False：直接使用当前已经在 Abaqus/CAE 中打开的模型数据库。
# 为保证复现稳定，通常建议保持 True。
OPEN_SOURCE_CAE = True

# 源显式分析模型数据库文件。使用文件名时，默认位于工作目录中。
SOURCE_CAE_FILE_NAME = "one_shot.cae"

# 显式计算结果文件。初始状态将从该 ODB 的最后一步、最后增量导入。
SOURCE_ODB_FILE_NAME = "Job-1.odb"

# True：找不到 ODB 时立即终止，避免生成无法计算的静态模型。
# False：仅给出警告并继续建模；写入 INP 时 Abaqus 可能仍会报错。
REQUIRE_SOURCE_ODB = True


# -----------------------------------------------------------------------------
# 1.2 模型、实例、分析步和作业名称
# -----------------------------------------------------------------------------

# one-shot.cae 中已有的显式源模型名称。
SOURCE_MODEL_NAME = "Model-1"

# 新建的静态模型名称。必须与 SOURCE_MODEL_NAME 不同。
STATIC_MODEL_NAME = "Model-1-Static"

# 源装配中的实例名称。
SHOT_INSTANCE_NAME = "shot-1"
TARGET_INSTANCE_NAME = "target-1"

# 新静态分析步和作业名称。
STATIC_STEP_NAME = "Step-1"
STATIC_JOB_NAME = "Job-1-Static"

# 导入初始状态的预定义场名称。
INITIAL_STATE_NAME = "Initial-State-from-Job-1"

# 若静态模型或静态作业已经存在，是否自动删除并重新创建。
REPLACE_EXISTING_STATIC_MODEL = True
REPLACE_EXISTING_STATIC_JOB = True


# -----------------------------------------------------------------------------
# 1.3 输出 CAE 文件
# -----------------------------------------------------------------------------

# True：将包含显式模型和静态模型的数据库另存为新 CAE，不覆盖源文件。
# False：直接保存到已经打开的源 CAE 文件中。
SAVE_AS_NEW_CAE = False

# 仅当 SAVE_AS_NEW_CAE=True 时使用。
# 默认根据 SOURCE_CAE_FILE_NAME 自动生成：
#   one-shot.cae -> one-shot-static.cae
#   impact_model.cae -> impact_model-static.cae
# 因此修改源 CAE 文件名后，无需再同步修改输出 CAE 文件名。
# 需要自定义名称时，也可以直接改成字符串，例如：
# OUTPUT_CAE_FILE_NAME = "my-static-result.cae"
_source_cae_name_without_ext, _source_cae_ext = os.path.splitext(
    SOURCE_CAE_FILE_NAME
)
OUTPUT_CAE_FILE_NAME = (
    _source_cae_name_without_ext
    + "-static"
    + (_source_cae_ext if _source_cae_ext else ".cae")
)


# -----------------------------------------------------------------------------
# 1.4 从显式 ODB 导入初始状态
# -----------------------------------------------------------------------------

# False：保留原始参考构型，只导入应力、应变等状态；对应 RPY 中的设置。
# True：将 ODB 最终变形构型更新为新的参考构型。
# 对于冲击后的残余应力平衡复现，通常保持 False。
UPDATE_REFERENCE_CONFIGURATION = False

# 当前脚本使用 ODB 的最后分析步和该步最后增量：
# endStep=LAST_STEP, endIncrement=STEP_END。
# 这两个设置与原 RPY 操作一致，通常无需修改。


# -----------------------------------------------------------------------------
# 1.5 需要从复制模型中抑制的对象
# -----------------------------------------------------------------------------

# 原始 RPY 的名称与参数化 one-shot 脚本的名称不同，因此同时列出两套候选名。
# 找不到某个名称时脚本会跳过，不会因此终止。
INTERACTIONS_TO_SUPPRESS = (
    "Int-1",
)

CONSTRAINTS_TO_SUPPRESS = (
    "Constraint-shot-rigid",  # 参数化 one-shot 脚本中的弹丸刚体约束
    "Constraint-1",           # 原始 RPY/GUI 操作中可能采用的名称
)

PREDEFINED_FIELDS_TO_SUPPRESS = (
    "Predefined-Field-shot-velocity",  # 参数化脚本名称
    "Predefined Field-1",              # 原始 RPY 名称
)

# True：除上述明确名称外，还自动抑制名称中包含 shot、rigid、velocity 的
# 相关对象，以兼容不同版本脚本产生的命名差异。
AUTO_SUPPRESS_SHOT_RELATED_OBJECTS = True

# 靶材底面边界条件必须保留。脚本用这些候选名称检查它是否存在。
TARGET_BOTTOM_BC_CANDIDATES = (
    "BC-target-bottom",
    "BC-1",
)

# True：若找不到可能的靶材固定边界条件则终止；False：只警告。
REQUIRE_TARGET_BOTTOM_BC = True


# -----------------------------------------------------------------------------
# 1.6 静力分析步参数
# -----------------------------------------------------------------------------

# 是否考虑大变形，即 Static, General 中的 Nlgeom。
STATIC_NLGEOM = True

# True：完全采用 Abaqus 创建 StaticStep 时的默认增量参数，最接近原 RPY：
#       StaticStep(name='Step-1', previous='Initial', nlgeom=ON)
# False：使用下面显式给出的时间和增量设置。
USE_ABAQUS_DEFAULT_INCREMENT_SETTINGS = True

# 仅当 USE_ABAQUS_DEFAULT_INCREMENT_SETTINGS=False 时生效。
STATIC_TIME_PERIOD = 1.0
STATIC_INITIAL_INCREMENT = 0.01
STATIC_MIN_INCREMENT = 1.0e-8
STATIC_MAX_INCREMENT = 0.10
STATIC_MAX_INCREMENTS = 1000


# -----------------------------------------------------------------------------
# 1.7 作业设置
# -----------------------------------------------------------------------------

# 用于计算的 CPU 核心数。
NUM_CPUS = 16

# True：自动令域数等于 CPU 核心数，这是本脚本的推荐设置。
# False：使用下面 NUM_DOMAINS 手动指定域数。
#
# Abaqus 要求域数必须是用于域计算的处理器数的整数倍。若默认域数为 1、
# 而 NUM_CPUS 设置为 12，就会出现“域的个数必须是用于域计算的处理器
# 个数的倍数”的错误。因此一般保持 True 即可。
AUTO_SET_NUM_DOMAINS = True

# 仅当 AUTO_SET_NUM_DOMAINS=False 时生效。
# 通常设置为与 NUM_CPUS 相同；也可设置为 NUM_CPUS 的整数倍。
NUM_DOMAINS = 16

MEMORY_PERCENT = 90
NUM_GPUS = 0
NODAL_OUTPUT_PRECISION = SINGLE

# True：在脚本末尾写出 Job-1-Static.inp。
WRITE_INPUT = True

# True：自动提交 Abaqus/Standard 分析。
# 首次运行建议保持 False，先检查初始状态、边界条件和静态分析步。
SUBMIT_JOB = False

# SUBMIT_JOB=True 时，是否等待计算结束后脚本才退出。
WAIT_FOR_COMPLETION = True

# 写 INP 和提交时是否执行一致性检查。
CONSISTENCY_CHECKING = OFF


# =============================================================================
# 2. 通用辅助函数：一般不需要修改
# =============================================================================


def _print_title(text):
    print("\n" + "=" * 78)
    print(text)
    print("=" * 78)


def _repository_names(repository):
    """将 Abaqus repository 的键转换成普通 Python 列表，便于打印和判断。"""
    return list(repository.keys())


def _validate_parameters():
    """运行前检查明显错误，尽早给出可读提示。"""
    if SOURCE_MODEL_NAME == STATIC_MODEL_NAME:
        raise ValueError("STATIC_MODEL_NAME 必须与 SOURCE_MODEL_NAME 不同。")

    if str(SOURCE_CAE_FILE_NAME).strip() == "":
        raise ValueError("SOURCE_CAE_FILE_NAME 不能为空。")

    if str(SOURCE_ODB_FILE_NAME).strip() == "":
        raise ValueError("SOURCE_ODB_FILE_NAME 不能为空。")

    if SAVE_AS_NEW_CAE and str(OUTPUT_CAE_FILE_NAME).strip() == "":
        raise ValueError(
            "SAVE_AS_NEW_CAE=True 时，OUTPUT_CAE_FILE_NAME 不能为空。"
        )

    if NUM_CPUS < 1:
        raise ValueError("NUM_CPUS 必须为正整数。")

    if not AUTO_SET_NUM_DOMAINS:
        if NUM_DOMAINS < 1:
            raise ValueError("NUM_DOMAINS 必须为正整数。")
        if NUM_DOMAINS % NUM_CPUS != 0:
            raise ValueError(
                "NUM_DOMAINS 必须是 NUM_CPUS 的整数倍。当前 NUM_CPUS=%d, "
                "NUM_DOMAINS=%d。通常将二者设置为相同值。"
                % (NUM_CPUS, NUM_DOMAINS)
            )

    if MEMORY_PERCENT <= 0 or MEMORY_PERCENT > 100:
        raise ValueError("MEMORY_PERCENT 必须位于 1~100。")

    if not USE_ABAQUS_DEFAULT_INCREMENT_SETTINGS:
        if STATIC_TIME_PERIOD <= 0.0:
            raise ValueError("STATIC_TIME_PERIOD 必须大于 0。")
        if STATIC_INITIAL_INCREMENT <= 0.0:
            raise ValueError("STATIC_INITIAL_INCREMENT 必须大于 0。")
        if STATIC_MIN_INCREMENT <= 0.0:
            raise ValueError("STATIC_MIN_INCREMENT 必须大于 0。")
        if STATIC_MAX_INCREMENT <= 0.0:
            raise ValueError("STATIC_MAX_INCREMENT 必须大于 0。")
        if STATIC_MAX_INCREMENTS < 1:
            raise ValueError("STATIC_MAX_INCREMENTS 必须为正整数。")
        if STATIC_MIN_INCREMENT > STATIC_INITIAL_INCREMENT:
            raise ValueError(
                "STATIC_MIN_INCREMENT 不应大于 STATIC_INITIAL_INCREMENT。"
            )
        if STATIC_INITIAL_INCREMENT > STATIC_MAX_INCREMENT:
            raise ValueError(
                "STATIC_INITIAL_INCREMENT 不应大于 STATIC_MAX_INCREMENT。"
            )


def _resolve_work_directory():
    """
    确定实际工作目录。

    USE_CURRENT_WORK_DIRECTORY=True：
        使用 Abaqus/CAE 当前工作目录。
    USE_CURRENT_WORK_DIRECTORY=False：
        使用 WORK_DIRECTORY。
    """
    if USE_CURRENT_WORK_DIRECTORY:
        work_directory = os.getcwd()
        source_text = "Abaqus/CAE 当前工作目录"
    else:
        work_directory = WORK_DIRECTORY
        source_text = "脚本参数 WORK_DIRECTORY"

    if work_directory is None or str(work_directory).strip() == "":
        raise ValueError("工作目录为空。")

    work_directory = os.path.abspath(str(work_directory))
    if not os.path.isdir(work_directory):
        raise IOError("工作目录不存在：%s" % work_directory)

    print("本次使用的工作目录：%s" % work_directory)
    print("工作目录来源：%s" % source_text)
    return work_directory


def _absolute_path(work_directory, file_name):
    """允许用户填写绝对路径，也允许只填写工作目录中的文件名。"""
    file_name = str(file_name)
    if os.path.isabs(file_name):
        return os.path.abspath(file_name)
    return os.path.abspath(os.path.join(work_directory, file_name))


def _open_or_get_database(source_cae_path):
    """
    打开源模型数据库，并重新绑定全局 mdb。

    与新建 Mdb() 类似，openMdb() 会改变 Abaqus 当前模型数据库。
    显式重新绑定 mdb 可避免继续访问旧 Python 引用。
    """
    global mdb

    if OPEN_SOURCE_CAE:
        if not os.path.isfile(source_cae_path):
            raise IOError("找不到源 CAE 文件：%s" % source_cae_path)

        print("正在打开源模型数据库：%s" % source_cae_path)
        opened_database = openMdb(pathName=source_cae_path)

        if opened_database is not None:
            mdb = opened_database
        else:
            mdb = abaqus_module.mdb
    else:
        mdb = abaqus_module.mdb
        print("使用当前已经打开的模型数据库。")

    print("当前模型数据库中的模型：%s" % _repository_names(mdb.models))
    return mdb


def _delete_existing_static_objects():
    """按参数删除同名旧静态作业和旧静态模型。"""
    # 先删除作业，避免作业仍引用即将删除的模型。
    if STATIC_JOB_NAME in mdb.jobs.keys():
        if REPLACE_EXISTING_STATIC_JOB:
            del mdb.jobs[STATIC_JOB_NAME]
            print("已删除旧作业：%s" % STATIC_JOB_NAME)
        else:
            raise ValueError(
                "作业 %s 已存在；请修改名称或启用 REPLACE_EXISTING_STATIC_JOB。"
                % STATIC_JOB_NAME
            )

    if STATIC_MODEL_NAME in mdb.models.keys():
        if REPLACE_EXISTING_STATIC_MODEL:
            del mdb.models[STATIC_MODEL_NAME]
            print("已删除旧静态模型：%s" % STATIC_MODEL_NAME)
        else:
            raise ValueError(
                "模型 %s 已存在；请修改名称或启用 REPLACE_EXISTING_STATIC_MODEL。"
                % STATIC_MODEL_NAME
            )


def _copy_source_model():
    """复制显式源模型，原模型保持不变。"""
    if SOURCE_MODEL_NAME not in mdb.models.keys():
        raise KeyError(
            "源模型 %s 不存在。当前模型：%s"
            % (SOURCE_MODEL_NAME, _repository_names(mdb.models))
        )

    model = mdb.Model(
        name=STATIC_MODEL_NAME,
        objectToCopy=mdb.models[SOURCE_MODEL_NAME],
    )
    print("已复制模型：%s -> %s" % (SOURCE_MODEL_NAME, STATIC_MODEL_NAME))
    return model


def _suppress_named_objects(repository, candidate_names, object_label):
    """抑制 repository 中存在的候选对象，并返回已抑制名称列表。"""
    suppressed_names = []

    for name in candidate_names:
        if name in repository.keys() and name not in suppressed_names:
            try:
                repository[name].suppress()
                suppressed_names.append(name)
                print("已抑制%s：%s" % (object_label, name))
            except Exception as error:
                print("警告：抑制%s %s 失败：%s" % (object_label, name, error))

    return suppressed_names


def _auto_suppress_by_keywords(repository, keywords, object_label, excluded_names):
    """按名称关键词补充抑制对象，兼容不同脚本的命名方式。"""
    suppressed_names = []
    excluded_lower = set([str(name).lower() for name in excluded_names])

    for name in list(repository.keys()):
        lower_name = str(name).lower()
        if lower_name in excluded_lower:
            continue

        matched = False
        for keyword in keywords:
            if str(keyword).lower() in lower_name:
                matched = True
                break

        if not matched:
            continue

        try:
            repository[name].suppress()
            suppressed_names.append(name)
            print("已自动抑制%s：%s" % (object_label, name))
        except Exception as error:
            print("警告：自动抑制%s %s 失败：%s" % (object_label, name, error))

    return suppressed_names


def _suppress_shot_instance(model):
    """在静态模型中抑制弹丸实例，只保留靶材参与静力平衡。"""
    assembly = model.rootAssembly

    if SHOT_INSTANCE_NAME not in assembly.features.keys():
        raise KeyError(
            "装配特征中找不到弹丸实例 %s。现有特征：%s"
            % (SHOT_INSTANCE_NAME, _repository_names(assembly.features))
        )

    try:
        assembly.features[SHOT_INSTANCE_NAME].suppress()
        print("已抑制弹丸实例：%s" % SHOT_INSTANCE_NAME)
    except Exception as error:
        # 如果实例已经处于抑制状态，Abaqus 可能返回异常；此时继续检查目标实例。
        print("警告：抑制弹丸实例时返回：%s" % error)

    assembly.regenerate()

    if TARGET_INSTANCE_NAME not in assembly.instances.keys():
        raise KeyError(
            "装配中找不到靶材实例 %s。当前活动实例：%s"
            % (TARGET_INSTANCE_NAME, _repository_names(assembly.instances))
        )

    return assembly


def _suppress_explicit_only_objects(model):
    """
    抑制只属于显式弹丸冲击阶段的对象。

    靶材底面固定边界条件不在此处抑制，因为静态平衡仍需要该约束。
    """
    suppressed_interactions = _suppress_named_objects(
        model.interactions,
        INTERACTIONS_TO_SUPPRESS,
        "相互作用",
    )

    suppressed_constraints = _suppress_named_objects(
        model.constraints,
        CONSTRAINTS_TO_SUPPRESS,
        "约束",
    )

    suppressed_fields = _suppress_named_objects(
        model.predefinedFields,
        PREDEFINED_FIELDS_TO_SUPPRESS,
        "预定义场",
    )

    if AUTO_SUPPRESS_SHOT_RELATED_OBJECTS:
        suppressed_interactions += _auto_suppress_by_keywords(
            model.interactions,
            ("shot", "contact", "int-"),
            "相互作用",
            suppressed_interactions,
        )

        suppressed_constraints += _auto_suppress_by_keywords(
            model.constraints,
            ("shot", "rigid"),
            "约束",
            suppressed_constraints,
        )

        suppressed_fields += _auto_suppress_by_keywords(
            model.predefinedFields,
            ("shot", "velocity"),
            "预定义场",
            suppressed_fields,
        )

    print("已抑制相互作用数量：%d" % len(set(suppressed_interactions)))
    print("已抑制约束数量：%d" % len(set(suppressed_constraints)))
    print("已抑制原初速度场数量：%d" % len(set(suppressed_fields)))


def _check_target_boundary_condition(model):
    """确认静态模型仍保留靶材底面固定边界条件。"""
    found = []
    for name in TARGET_BOTTOM_BC_CANDIDATES:
        if name in model.boundaryConditions.keys():
            found.append(name)

    if found:
        print("检测到靶材约束边界条件：%s" % found)
        return

    message = (
        "未找到候选靶材底面边界条件 %s。当前边界条件：%s"
        % (TARGET_BOTTOM_BC_CANDIDATES, _repository_names(model.boundaryConditions))
    )

    if REQUIRE_TARGET_BOTTOM_BC:
        raise KeyError(message)
    print("警告：%s" % message)


def _create_initial_state(model, assembly, source_odb_path):
    """从显式 ODB 最后一步、最后增量导入靶材初始状态。"""
    if not os.path.isfile(source_odb_path):
        message = "找不到显式结果文件：%s" % source_odb_path
        if REQUIRE_SOURCE_ODB:
            raise IOError(message)
        print("警告：%s" % message)

    # InitialState 的 fileName 通常使用不带 .odb 的作业名或路径。
    # 这里先切换到 ODB 所在目录，再传入不带扩展名的文件名，兼容性最好。
    odb_directory = os.path.dirname(source_odb_path)
    odb_job_name = os.path.splitext(os.path.basename(source_odb_path))[0]
    os.chdir(odb_directory)

    if INITIAL_STATE_NAME in model.predefinedFields.keys():
        del model.predefinedFields[INITIAL_STATE_NAME]
        print("已删除旧初始状态预定义场：%s" % INITIAL_STATE_NAME)

    target_instance = assembly.instances[TARGET_INSTANCE_NAME]
    instances = (target_instance,)

    model.InitialState(
        name=INITIAL_STATE_NAME,
        createStepName="Initial",
        fileName=odb_job_name,
        endStep=LAST_STEP,
        endIncrement=STEP_END,
        updateReferenceConfiguration=(
            ON if UPDATE_REFERENCE_CONFIGURATION else OFF
        ),
        instances=instances,
    )

    print("已创建初始状态：%s" % INITIAL_STATE_NAME)
    print("初始状态来源：%s，最后一步、最后增量" % source_odb_path)
    print(
        "更新参考构型：%s"
        % ("是" if UPDATE_REFERENCE_CONFIGURATION else "否")
    )


def _delete_source_analysis_steps(model):
    """
    删除复制模型中的所有非 Initial 分析步。

    原 one-shot 模型采用 ExplicitDynamicsStep，Abaqus/Standard 静态作业中不能
    保留该显式分析步，因此必须先删除，再创建 StaticStep。
    """
    step_names = [name for name in list(model.steps.keys()) if name != "Initial"]

    for step_name in reversed(step_names):
        del model.steps[step_name]
        print("已删除源分析步：%s" % step_name)

    remaining_steps = _repository_names(model.steps)
    print("删除后模型中的分析步：%s" % remaining_steps)


def _create_static_step(model):
    """创建 Static, General 分析步。"""
    nlgeom_value = ON if STATIC_NLGEOM else OFF

    if USE_ABAQUS_DEFAULT_INCREMENT_SETTINGS:
        model.StaticStep(
            name=STATIC_STEP_NAME,
            previous="Initial",
            nlgeom=nlgeom_value,
        )
        print("静态步采用 Abaqus 默认时间与增量设置。")
    else:
        model.StaticStep(
            name=STATIC_STEP_NAME,
            previous="Initial",
            nlgeom=nlgeom_value,
            timePeriod=STATIC_TIME_PERIOD,
            maxNumInc=STATIC_MAX_INCREMENTS,
            initialInc=STATIC_INITIAL_INCREMENT,
            minInc=STATIC_MIN_INCREMENT,
            maxInc=STATIC_MAX_INCREMENT,
        )
        print(
            "静态步增量：timePeriod=%g, initial=%g, min=%g, max=%g, maxNumInc=%d"
            % (
                STATIC_TIME_PERIOD,
                STATIC_INITIAL_INCREMENT,
                STATIC_MIN_INCREMENT,
                STATIC_MAX_INCREMENT,
                STATIC_MAX_INCREMENTS,
            )
        )

    print("已创建静态分析步：%s" % STATIC_STEP_NAME)
    print("几何非线性 NLGEOM：%s" % ("开启" if STATIC_NLGEOM else "关闭"))


def _get_num_domains():
    """返回实际使用的域数，并保证满足 Abaqus 的并行设置要求。"""
    if AUTO_SET_NUM_DOMAINS:
        return int(NUM_CPUS)

    domains = int(NUM_DOMAINS)
    if domains < 1:
        raise ValueError("NUM_DOMAINS 必须为正整数。")
    if domains % int(NUM_CPUS) != 0:
        raise ValueError(
            "NUM_DOMAINS 必须是 NUM_CPUS 的整数倍。当前 NUM_CPUS=%d, "
            "NUM_DOMAINS=%d。" % (NUM_CPUS, domains)
        )
    return domains


def _create_static_job(model):
    """创建 Abaqus/Standard 静态作业。"""
    if STATIC_JOB_NAME in mdb.jobs.keys():
        if REPLACE_EXISTING_STATIC_JOB:
            del mdb.jobs[STATIC_JOB_NAME]
        else:
            raise ValueError("作业已存在：%s" % STATIC_JOB_NAME)

    actual_num_domains = _get_num_domains()

    print("并行计算设置：CPU 核心数=%d，域数=%d" % (NUM_CPUS, actual_num_domains))

    job = mdb.Job(
        name=STATIC_JOB_NAME,
        model=STATIC_MODEL_NAME,
        description="Static equilibrium after single-shot impact",
        type=ANALYSIS,
        atTime=None,
        waitMinutes=0,
        waitHours=0,
        queue=None,
        memory=MEMORY_PERCENT,
        memoryUnits=PERCENTAGE,
        getMemoryFromAnalysis=True,
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
        numDomains=actual_num_domains,
        numGPUs=NUM_GPUS,
    )

    print("已创建 Abaqus/Standard 作业：%s" % STATIC_JOB_NAME)
    return job


def _save_database(source_cae_path, output_cae_path):
    """根据参数另存或保存模型数据库。"""
    if SAVE_AS_NEW_CAE:
        # 防止用户将输出文件误设为源文件，却以为保留了原始数据库。
        if os.path.normcase(output_cae_path) == os.path.normcase(source_cae_path):
            raise ValueError(
                "SAVE_AS_NEW_CAE=True 时，OUTPUT_CAE_FILE_NAME 不能与源 CAE 相同。"
            )
        mdb.saveAs(pathName=output_cae_path)
        print("模型数据库已另存为：%s" % output_cae_path)
    else:
        mdb.save()
        print("模型数据库已保存到当前源 CAE。")


def _write_or_submit_job(job, work_directory):
    """按开关写入 INP 或提交静态作业。"""
    os.chdir(work_directory)

    if WRITE_INPUT:
        job.writeInput(consistencyChecking=CONSISTENCY_CHECKING)
        print("作业输入文件已写入：%s.inp" % STATIC_JOB_NAME)

    if SUBMIT_JOB:
        print("正在提交静态作业：%s" % STATIC_JOB_NAME)
        job.submit(consistencyChecking=CONSISTENCY_CHECKING)

        if WAIT_FOR_COMPLETION:
            job.waitForCompletion()
            print("作业计算结束：%s" % STATIC_JOB_NAME)


def main():
    """主流程。"""
    _print_title("Abaqus/Standard 单颗喷丸冲击后静态分析参数化建模")

    _validate_parameters()
    work_directory = _resolve_work_directory()
    source_cae_path = _absolute_path(work_directory, SOURCE_CAE_FILE_NAME)
    source_odb_path = _absolute_path(work_directory, SOURCE_ODB_FILE_NAME)
    output_cae_path = _absolute_path(work_directory, OUTPUT_CAE_FILE_NAME)

    print("源 CAE：%s" % source_cae_path)
    print("源 ODB：%s" % source_odb_path)
    if SAVE_AS_NEW_CAE:
        print("输出 CAE：%s" % output_cae_path)

    _print_title("1/7 打开显式 one-shot 模型数据库")
    _open_or_get_database(source_cae_path)

    _print_title("2/7 复制源模型并创建静态模型")
    _delete_existing_static_objects()
    static_model = _copy_source_model()

    _print_title("3/7 抑制弹丸及显式冲击专用对象")
    assembly = _suppress_shot_instance(static_model)
    _suppress_explicit_only_objects(static_model)
    _check_target_boundary_condition(static_model)

    _print_title("4/7 从显式 ODB 导入靶材初始状态")
    _create_initial_state(static_model, assembly, source_odb_path)

    _print_title("5/7 删除显式分析步并创建静态分析步")
    _delete_source_analysis_steps(static_model)
    _create_static_step(static_model)

    _print_title("6/7 创建静态作业并保存 CAE")
    static_job = _create_static_job(static_model)
    _save_database(source_cae_path, output_cae_path)

    _print_title("7/7 写入 INP 或提交分析")
    _write_or_submit_job(static_job, work_directory)

    # 最后再保存一次，将作业和输入文件写入状态记录进 CAE。
    mdb.save()

    _print_title("脚本执行完成")
    print("静态模型：%s" % STATIC_MODEL_NAME)
    print("静态分析步：%s" % STATIC_STEP_NAME)
    print("静态作业：%s" % STATIC_JOB_NAME)
    print("SUBMIT_JOB=%s" % SUBMIT_JOB)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print("\n脚本执行失败：%s" % error)
        traceback.print_exc()
        raise
