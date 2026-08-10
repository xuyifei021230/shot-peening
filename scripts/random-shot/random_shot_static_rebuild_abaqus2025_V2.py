# -*- coding: utf-8 -*-
"""
Abaqus/CAE 2025 随机多弹丸喷丸后的静态弹性/静力平衡分析脚本 V2
===================================================================

用途
----
本脚本以已经完成 Abaqus/Explicit 随机多弹丸冲击计算的 CAE 与 ODB 为基础，
自动复制显式模型，并建立 Abaqus/Standard 的 Static, General 静态分析模型。

V2 说明
-------
V2 以“random_shot_static_elastic_rebuild_abaqus2025改.py”中已经通过
Job-5-Static 实算测试的参数为功能基线，只补充更细的中文解释和阅读导航，
不改变原修改版的模型转换逻辑。默认参数之间的对应关系为：

    显式源模型      Model-5
    显式结果        Job-5.odb
             ↓  Initial State（最后一步、最后增量）
    静态复制模型    Model-5-Static
    静态作业        Job-5-Static
    静态结果        Job-5-Static.odb
             ↓  extract_random_shot_PE_for_Almen.py
    分层 PE 文本    Job-5-Static-PE.txt

术语说明
--------
“源模型”指包含靶材、全部弹丸、刚体约束、接触和显式分析步的模型。
“静态模型”指源模型的副本；脚本只修改副本，不直接删改源模型中的对象。
“抑制”不是删除：对象仍保留在模型树中，但不参与静态作业，便于追溯。
“Initial State”把显式 ODB 中的应力、塑性应变等材料状态映射到静态模型。
Static, General 步随后在保留靶材约束的条件下求解静力平衡/卸载状态。

脚本复现的视频流程
------------------
1. 复制随机喷丸显式模型，原模型保持不变；
2. 抑制全部随机弹丸实例；
3. 抑制全部弹丸刚体约束；
4. 抑制显式通用接触；
5. 抑制弹丸初速度预定义场；
6. 保留靶材的底面/侧面固定边界条件；
7. 从显式 ODB 的最后一步、最后增量向靶材实例导入 Initial State；
8. 不更新参考构型，即 updateReferenceConfiguration=OFF；
9. 删除复制模型中的 Dynamic, Explicit 分析步；
10. 创建 Static, General，默认 NLGEOM=ON、时间周期为 1；
11. 创建 Abaqus/Standard 静态作业；
12. 另存新的 CAE，并按开关写出 INP 或提交计算。

与单喷丸静态脚本相比
--------------------
单喷丸只需要抑制一个 shot 实例和一个刚体约束；本脚本会自动识别并批量抑制：
    shot-001、shot-002、...、shot-050
    Rigid-001、Rigid-002、...、Rigid-050
同时兼容视频中的 shot-0～shot-49、rigid-0～rigid-49 命名。

单位制
------
本脚本不重新创建几何、材料和网格，单位制必须与源随机喷丸模型一致。
前述随机喷丸脚本采用：mm、s、tonne、MPa。

运行方式
--------
Abaqus/CAE -> File（文件）-> Run Script（运行脚本）-> 选择本文件。

重要说明
--------
1. SOURCE_ODB_FILE_NAME 必须来自与 SOURCE_CAE_FILE_NAME 中源模型完全一致的网格。
2. 初始状态默认读取 ODB 的最后一步和最后增量。
3. 默认不更新参考构型，只导入最终变形网格和材料状态用于静力再平衡；
   该设置与视频和先前单喷丸静态脚本一致。
4. “静态弹性分析”在这里指 Abaqus/Standard 的静态平衡步骤。脚本不会把
   TC4 的 Johnson-Cook 塑性材料强制改成纯弹性材料，以免丢失显式冲击后的
   塑性应变和残余状态。
5. 默认只写出 INP，不自动提交。首次运行后应先检查静态模型、Initial State、
   边界条件和分析步，再决定是否提交。
6. 脚本必须在 Abaqus/CAE 中运行，不能用系统 Python 直接执行；因为 mdb、
   session、InitialState、StaticStep 等对象只存在于 Abaqus Python 环境。
7. 若 OPEN_SOURCE_CAE=False，运行前必须先在 CAE 中打开正确的数据库；
   SOURCE_CAE_FILE_NAME 此时主要用于日志和推导另存文件名，不负责打开文件。
8. 若 SAVE_AS_NEW_CAE=False，mdb.save() 会保存当前已打开的 CAE。该设置符合
   你的测试版，但正式批量计算时建议改为 True，以免静态副本写回源 CAE。
"""

from __future__ import print_function

# Abaqus 主对象：mdb 管理模型/作业，session 管理当前 CAE 会话和视口。
from abaqus import *
# ON、OFF、LAST_STEP、STEP_END、ANALYSIS 等符号常量均来自这里。
from abaqusConstants import *
# 导入 CAE 各模块所需对象；脚本虽然主要操作模型数据库，但仍需 CAE 环境。
from caeModules import *

import os
import traceback
# 单独保留 abaqus 模块引用。openMdb 后用 abaqus_module.mdb 重新取得当前数据库，
# 可避免某些 Abaqus 版本中通配符导入的 mdb 引用没有及时更新。
import abaqus as abaqus_module


# =============================================================================
# 一、用户可修改参数区
# =============================================================================
# 日常使用通常只需要修改本区域。建议按以下顺序核对：
#   ① 工作目录和 ODB 文件名；② 源/静态模型名；③ 靶材实例名；
#   ④ 是否另存 CAE；⑤ CPU 数；⑥ WRITE_INPUT / SUBMIT_JOB 开关。

# -----------------------------------------------------------------------------
# 1. 工作目录与源文件
# -----------------------------------------------------------------------------

# True：使用 Abaqus/CAE 当前工作目录，即 os.getcwd()。
#       适合你当前“先设置工作目录，再运行脚本”的操作方式。
# False：忽略当前目录，固定使用下面 WORK_DIRECTORY 指定的目录。
USE_CURRENT_WORK_DIRECTORY = True

# 仅当 USE_CURRENT_WORK_DIRECTORY=False 时生效。
WORK_DIRECTORY = r"D:\temp\random-shot"

# True：脚本主动打开 SOURCE_CAE_FILE_NAME 指定的源 CAE。
#       适合从空白 CAE 会话开始运行。
# False：直接使用 Abaqus/CAE 当前已经打开的模型数据库。
#        你的测试记录采用此模式；务必确认 Model-5 位于当前数据库中。
OPEN_SOURCE_CAE = False

# SOURCE_CAE_FILE_NAME 与 SOURCE_ODB_FILE_NAME 可写文件名或绝对路径。
# 当只写文件名时，它们会与前面解析出的工作目录拼接。
# ODB 必须来自同一套靶材网格，否则 Initial State 导入会因标签/网格不匹配失败。
SOURCE_CAE_FILE_NAME = "Random-Shot.cae"
SOURCE_ODB_FILE_NAME = "Job-1.odb"

# 找不到源 ODB 时是否立即终止。
REQUIRE_SOURCE_ODB = True


# -----------------------------------------------------------------------------
# 2. 源模型、静态模型、实例与作业名称
# -----------------------------------------------------------------------------

# 显式源模型名。你的测试数据库包含 Model-1～Model-5，此处明确选择 Model-5。
SOURCE_MODEL_NAME = "Model-1"

# 若 SOURCE_MODEL_NAME 不存在，是否根据 target/shot 实例自动识别源模型。
AUTO_DETECT_SOURCE_MODEL = True

# 下面三个名称共同决定静态模型、分析步以及最终 ODB 的名称：
#   STATIC_JOB_NAME="Job-5-Static" -> 求解后生成 Job-5-Static.odb。
STATIC_MODEL_NAME = "Model-1-Static"
STATIC_STEP_NAME = "Step-1"
STATIC_JOB_NAME = "Job-1-Static"
# Initial State 会出现在静态模型的 Predefined Fields 下。
INITIAL_STATE_NAME = "Initial-State-from-Random-Shot"

# 源随机喷丸装配中的靶材实例。V4 脚本使用 target-1。
TARGET_INSTANCE_NAME = "target-1"

# 找不到上述名称时，是否自动从实例名中搜索唯一的 target 实例。
AUTO_DETECT_TARGET_INSTANCE = True

# 重新运行脚本时，是否删除同名旧静态模型和旧静态作业。
REPLACE_EXISTING_STATIC_MODEL = True
REPLACE_EXISTING_STATIC_JOB = True


# -----------------------------------------------------------------------------
# 3. 输出 CAE
# -----------------------------------------------------------------------------

# True：通过 mdb.saveAs() 另存新 CAE，不覆盖显式源 CAE。
# False：通过 mdb.save() 保存当前数据库。你的测试版采用 False。
SAVE_AS_NEW_CAE = False

# None：根据 SOURCE_CAE_FILE_NAME 自动生成，例如：
#   Random-Shot-Model.cae -> Random-Shot-Model-static.cae
# 也可以直接填写字符串，例如 "my-random-shot-static.cae"。
OUTPUT_CAE_FILE_NAME = None


# -----------------------------------------------------------------------------
# 4. 弹丸实例、刚体约束、接触和速度场识别规则
# -----------------------------------------------------------------------------

# 弹丸实例前缀。比较时不区分大小写；既兼容当前脚本生成的 shot-001，
# 也兼容视频旧模型中的 shot-0。只要实例名以其中任一字符串开头就会被抑制。
SHOT_INSTANCE_PREFIXES = (
    "shot-",
    "shot_",
)

# 刚体约束前缀/关键词。除名称匹配外，程序还会检查对象类名中是否含 rigid，
# 因此 V4 的 Rigid-001 和视频中的 rigid-0 均可识别。
RIGID_CONSTRAINT_PREFIXES = (
    "rigid-",
    "rigid_",
    "constraint-shot",
    "constraint-rigid",
)

# True：抑制复制模型中的全部相互作用。
# 随机喷丸模型中通常只有显式 General Contact：Int-1。
SUPPRESS_ALL_INTERACTIONS = True

# SUPPRESS_ALL_INTERACTIONS=False 时，按这些名称和关键词抑制。
INTERACTION_NAMES_OR_KEYWORDS = (
    "Int-1",
    "contact",
    "shot",
)

# 需要抑制的弹丸初速度候选名称。
# 即使视频中名称为 Predefined Field-1，脚本也会根据对象类型 Velocity 识别。
VELOCITY_FIELD_NAMES_OR_KEYWORDS = (
    "PREDEFINED-SHOT-VELOCITY",
    "Predefined-Field-shot-velocity",
    "Predefined Field-1",
    "velocity",
    "shot",
)

# True：找不到任何弹丸实例或刚体约束时终止，防止选错源模型。
REQUIRE_SHOT_INSTANCES = True
REQUIRE_RIGID_CONSTRAINTS = True


# -----------------------------------------------------------------------------
# 5. 靶材边界条件
# -----------------------------------------------------------------------------

# 此前随机喷丸 V4 脚本使用 BC-TARGET-ENCASTRE。
TARGET_FIXED_BC_CANDIDATES = (
    "BC-TARGET-ENCASTRE",
    "BC-target-bottom",
    "BC-1",
)

# True：找不到靶材固定边界条件时终止。
REQUIRE_TARGET_FIXED_BC = True


# -----------------------------------------------------------------------------
# 6. Initial State 导入设置
# -----------------------------------------------------------------------------

# False：不更新参考构型；与视频和单喷丸静态脚本一致。
#        原始参考坐标保持不变，但显式末态中的应力/塑性应变会被导入。
# True：将显式计算最终变形构型作为新的参考构型。
#       两者物理含义不同，除非已确定固有应变流程需要，否则不要随意切换。
UPDATE_REFERENCE_CONFIGURATION = False

# 当前固定采用最后一步、最后增量：
#   endStep=LAST_STEP
#   endIncrement=STEP_END


# -----------------------------------------------------------------------------
# 7. Static, General 分析步设置
# -----------------------------------------------------------------------------

# ON 允许静态再平衡考虑几何非线性；喷丸后靶材可能已有变形，因此默认开启。
STATIC_NLGEOM = True

# True：采用 Abaqus 创建 StaticStep 时的默认增量参数，最接近视频操作。
# False：使用下面的手动参数。
USE_ABAQUS_DEFAULT_INCREMENT_SETTINGS = True

STATIC_TIME_PERIOD = 1.0
STATIC_INITIAL_INCREMENT = 0.01
STATIC_MIN_INCREMENT = 1.0e-8
STATIC_MAX_INCREMENT = 0.10
STATIC_MAX_INCREMENTS = 1000


# -----------------------------------------------------------------------------
# 8. 静态作业设置
# -----------------------------------------------------------------------------

# 你的修改版设置为 12 核。该参数影响静态求解，不影响模型转换结果。
NUM_CPUS = 16

# True：实际域数自动等于 NUM_CPUS，此时下方 NUM_DOMAINS=4 不参与计算。
# False：使用手工 NUM_DOMAINS，并要求它是 NUM_CPUS 的整数倍。
AUTO_SET_NUM_DOMAINS = True
NUM_DOMAINS = 4

MEMORY_PERCENT = 90
NUM_GPUS = 0
NODAL_OUTPUT_PRECISION = SINGLE

# WRITE_INPUT=True：生成 Job-5-Static.inp，便于提交前检查 *Import 和集合。
# SUBMIT_JOB=True：脚本自动提交 Standard 求解；False 时只建模/写 INP。
# WAIT_FOR_COMPLETION 只在 SUBMIT_JOB=True 时有效。
# 首次运行建议 WRITE_INPUT=True、SUBMIT_JOB=False。
WRITE_INPUT = True
SUBMIT_JOB = False
WAIT_FOR_COMPLETION = True
CONSISTENCY_CHECKING = OFF


# -----------------------------------------------------------------------------
# 9. 日志与显示
# -----------------------------------------------------------------------------

WRITE_SUMMARY_LOG = False
SUMMARY_LOG_FILE_NAME = "Random-Shot-Static-Summary.txt"
DISPLAY_STATIC_MODEL = True


# =============================================================================
# 二、通用辅助函数
# =============================================================================

_ACTION_LOG = []


def log(message):
    """
    向 Abaqus 命令区输出带统一前缀的信息，并缓存到 _ACTION_LOG。

    缓存内容会在流程结束时写入 SUMMARY_LOG_FILE_NAME，方便关闭 CAE 后复核
    脚本实际选择了哪个模型、ODB、实例和作业参数。
    """
    text = str(message)
    print("[RandomShotStatic] " + text)
    _ACTION_LOG.append(text)


def print_title(text):
    """输出阶段标题；主流程用 1/8～8/8 标识执行位置。"""
    line = "=" * 78
    print("\n" + line)
    print(text)
    print(line)
    _ACTION_LOG.append("")
    _ACTION_LOG.append(text)


def repository_names(repository):
    """
    将 Abaqus Repository 的键转换成普通 Python 列表。

    mdb.models、model.steps、assembly.instances 等都不是普通 dict；不同
    Abaqus 版本的 keys() 返回类型略有差异，转换后更便于记录和遍历。
    """
    try:
        return list(repository.keys())
    except Exception:
        return []


def repository_contains(repository, name):
    """通过直接索引稳健检查 Repository 中是否存在指定名称。"""
    try:
        repository[name]
        return True
    except Exception:
        return False


def lower_text(value):
    """安全转为小写字符串，用于不区分大小写的名称匹配。"""
    try:
        return str(value).lower()
    except Exception:
        return ""


def starts_with_any(name, prefixes):
    """判断名称是否以任一候选前缀开头。"""
    lower_name = lower_text(name)
    for prefix in prefixes:
        if lower_name.startswith(lower_text(prefix)):
            return True
    return False


def contains_any(name, keywords):
    """判断名称中是否包含任一候选关键词。"""
    lower_name = lower_text(name)
    for keyword in keywords:
        if lower_text(keyword) in lower_name:
            return True
    return False


def validate_parameters():
    """
    在修改模型数据库前检查明显的参数错误。

    此处只检查名称冲突、CPU/域关系、内存百分比和手动增量参数。尽早终止
    可以避免脚本复制大模型以后才因简单参数错误失败。
    """
    if SOURCE_MODEL_NAME == STATIC_MODEL_NAME:
        raise ValueError("SOURCE_MODEL_NAME 与 STATIC_MODEL_NAME 不能相同。")

    if str(SOURCE_CAE_FILE_NAME).strip() == "":
        raise ValueError("SOURCE_CAE_FILE_NAME 不能为空。")

    if str(SOURCE_ODB_FILE_NAME).strip() == "":
        raise ValueError("SOURCE_ODB_FILE_NAME 不能为空。")

    if NUM_CPUS < 1:
        raise ValueError("NUM_CPUS 必须为正整数。")

    if not AUTO_SET_NUM_DOMAINS:
        if NUM_DOMAINS < 1:
            raise ValueError("NUM_DOMAINS 必须为正整数。")
        if NUM_DOMAINS % NUM_CPUS != 0:
            raise ValueError(
                "NUM_DOMAINS 必须是 NUM_CPUS 的整数倍。当前 NUM_CPUS=%d，"
                "NUM_DOMAINS=%d。" % (NUM_CPUS, NUM_DOMAINS)
            )

    if MEMORY_PERCENT <= 0 or MEMORY_PERCENT > 100:
        raise ValueError("MEMORY_PERCENT 必须位于 1～100。")

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
            raise ValueError("STATIC_MIN_INCREMENT 不应大于 STATIC_INITIAL_INCREMENT。")
        if STATIC_INITIAL_INCREMENT > STATIC_MAX_INCREMENT:
            raise ValueError("STATIC_INITIAL_INCREMENT 不应大于 STATIC_MAX_INCREMENT。")


def resolve_work_directory():
    """
    确定源 CAE、源 ODB、输出 INP 和摘要日志依赖的基准目录。

    返回值始终转换成绝对路径，并确认目录真实存在；函数不会自动创建目录。
    """
    if USE_CURRENT_WORK_DIRECTORY:
        directory = os.getcwd()
        source = "Abaqus/CAE 当前工作目录"
    else:
        directory = WORK_DIRECTORY
        source = "参数 WORK_DIRECTORY"

    if directory is None or str(directory).strip() == "":
        raise ValueError("工作目录为空。")

    directory = os.path.abspath(str(directory))
    if not os.path.isdir(directory):
        raise IOError("工作目录不存在：%s" % directory)

    log("工作目录：%s" % directory)
    log("工作目录来源：%s" % source)
    return directory


def absolute_path(base_directory, file_name):
    """把“仅文件名”与基准目录拼接；绝对路径则保持其目录含义。"""
    text = str(file_name)
    if os.path.isabs(text):
        return os.path.abspath(text)
    return os.path.abspath(os.path.join(base_directory, text))


def derive_output_cae_name():
    """未显式指定 OUTPUT_CAE_FILE_NAME 时，在源文件名后追加 -static。"""
    if OUTPUT_CAE_FILE_NAME is not None:
        if str(OUTPUT_CAE_FILE_NAME).strip() == "":
            raise ValueError("OUTPUT_CAE_FILE_NAME 不能是空字符串。")
        return str(OUTPUT_CAE_FILE_NAME)

    source_base = os.path.basename(str(SOURCE_CAE_FILE_NAME))
    root, extension = os.path.splitext(source_base)
    if extension == "":
        extension = ".cae"
    return root + "-static" + extension


def open_or_get_database(source_cae_path):
    """
    根据 OPEN_SOURCE_CAE 打开源 CAE，或取得当前已打开数据库。

    OPEN_SOURCE_CAE=False 时不会核对当前 CAE 文件名；后续依靠模型名称及
    target/shot 装配结构识别来降低选错数据库的风险。
    """
    global mdb

    if OPEN_SOURCE_CAE:
        if not os.path.isfile(source_cae_path):
            raise IOError("找不到源 CAE 文件：%s" % source_cae_path)

        log("正在打开源 CAE：%s" % source_cae_path)
        opened_database = openMdb(pathName=source_cae_path)
        if opened_database is not None:
            mdb = opened_database
        else:
            mdb = abaqus_module.mdb
    else:
        mdb = abaqus_module.mdb
        log("使用 Abaqus/CAE 当前已经打开的模型数据库。")

    log("当前数据库中的模型：%s" % repository_names(mdb.models))
    return mdb


def model_has_random_shot_structure(model):
    """
    判断候选模型是否同时具有 target 实例和至少一个 shot 实例。

    这只用于源模型自动识别；ODB 与模型的网格一致性仍由 Abaqus 检查。
    """
    assembly = model.rootAssembly
    instance_names = repository_names(assembly.instances)

    has_target = False
    has_shot = False
    for name in instance_names:
        lower_name = lower_text(name)
        if "target" in lower_name:
            has_target = True
        if starts_with_any(name, SHOT_INSTANCE_PREFIXES):
            has_shot = True

    return has_target and has_shot


def resolve_source_model_name():
    """
    确定源随机喷丸模型名称。

    优先使用 SOURCE_MODEL_NAME；若不存在且允许自动识别，则排除名称含
    static 的模型，再寻找唯一具有 target+shot 结构的候选模型。
    """
    if repository_contains(mdb.models, SOURCE_MODEL_NAME):
        log("采用指定源模型：%s" % SOURCE_MODEL_NAME)
        return SOURCE_MODEL_NAME

    if not AUTO_DETECT_SOURCE_MODEL:
        raise KeyError(
            "找不到源模型 %s。当前模型：%s"
            % (SOURCE_MODEL_NAME, repository_names(mdb.models))
        )

    candidates = []
    for name in repository_names(mdb.models):
        # 避免把之前创建的静态模型再次当成显式源模型。
        if "static" in lower_text(name):
            continue
        try:
            if model_has_random_shot_structure(mdb.models[name]):
                candidates.append(name)
        except Exception:
            pass

    if len(candidates) == 1:
        log("指定源模型不存在，自动识别为：%s" % candidates[0])
        return candidates[0]

    raise KeyError(
        "无法唯一识别随机喷丸源模型。指定名称=%s，候选模型=%s，全部模型=%s。"
        % (SOURCE_MODEL_NAME, candidates, repository_names(mdb.models))
    )


def delete_existing_static_objects():
    """
    清理上一次运行遗留的同名静态作业和静态模型。

    先删作业再删模型，避免作业仍引用旧模型；显式源模型不会在此删除。
    """
    if repository_contains(mdb.jobs, STATIC_JOB_NAME):
        if REPLACE_EXISTING_STATIC_JOB:
            del mdb.jobs[STATIC_JOB_NAME]
            log("已删除旧静态作业：%s" % STATIC_JOB_NAME)
        else:
            raise ValueError("静态作业已存在：%s" % STATIC_JOB_NAME)

    if repository_contains(mdb.models, STATIC_MODEL_NAME):
        if REPLACE_EXISTING_STATIC_MODEL:
            del mdb.models[STATIC_MODEL_NAME]
            log("已删除旧静态模型：%s" % STATIC_MODEL_NAME)
        else:
            raise ValueError("静态模型已存在：%s" % STATIC_MODEL_NAME)


def copy_source_model(source_model_name):
    """
    用 Model(objectToCopy=...) 深复制显式源模型，原模型保持不变。

    零件、装配、材料、网格、集合、边界条件和分析对象都会进入静态副本；
    后续所有 suppress/delete 操作只针对 STATIC_MODEL_NAME。
    """
    source_model = mdb.models[source_model_name]
    static_model = mdb.Model(
        name=STATIC_MODEL_NAME,
        objectToCopy=source_model,
    )
    log("已复制模型：%s -> %s" % (source_model_name, STATIC_MODEL_NAME))
    return static_model


def resolve_target_instance_name(model):
    """优先采用 target-1；否则按 target 关键词寻找唯一靶材实例。"""
    assembly = model.rootAssembly

    if repository_contains(assembly.instances, TARGET_INSTANCE_NAME):
        log("采用指定靶材实例：%s" % TARGET_INSTANCE_NAME)
        return TARGET_INSTANCE_NAME

    if not AUTO_DETECT_TARGET_INSTANCE:
        raise KeyError(
            "找不到靶材实例 %s。当前实例：%s"
            % (TARGET_INSTANCE_NAME, repository_names(assembly.instances))
        )

    candidates = []
    for name in repository_names(assembly.instances):
        if "target" in lower_text(name):
            candidates.append(name)

    if len(candidates) == 1:
        log("指定靶材实例不存在，自动识别为：%s" % candidates[0])
        return candidates[0]

    raise KeyError(
        "无法唯一识别靶材实例。指定名称=%s，候选=%s，全部实例=%s。"
        % (TARGET_INSTANCE_NAME, candidates, repository_names(assembly.instances))
    )


def find_shot_feature_names(model, target_instance_name):
    """
    从 assembly.features 与 assembly.instances 双重识别全部弹丸实例。

    suppress() 实际作用于装配特征；补充读取 instances 是为了兼容部分版本中
    两个 Repository 的名称视图暂时不一致的情况。
    """
    assembly = model.rootAssembly
    candidates = []

    for name in repository_names(assembly.features):
        if lower_text(name) == lower_text(target_instance_name):
            continue
        if starts_with_any(name, SHOT_INSTANCE_PREFIXES):
            candidates.append(name)

    # 部分模型中的 feature/instance 名称视图可能不同，补充从 instances 搜索。
    for name in repository_names(assembly.instances):
        if lower_text(name) == lower_text(target_instance_name):
            continue
        if starts_with_any(name, SHOT_INSTANCE_PREFIXES) and name not in candidates:
            candidates.append(name)

    return sorted(candidates)


def suppress_shot_instances(model, target_instance_name):
    """
    批量抑制全部随机弹丸装配特征，只保留靶材参与静态求解。

    弹丸已完成冲击作用，静态卸载阶段继续保留它们会增加自由度和接触问题。
    抑制后 regenerate() 用于刷新装配状态。
    """
    assembly = model.rootAssembly
    shot_names = find_shot_feature_names(model, target_instance_name)

    if len(shot_names) == 0:
        message = "未识别到任何弹丸实例。"
        if REQUIRE_SHOT_INSTANCES:
            raise RuntimeError(message + " 当前装配特征：%s" % repository_names(assembly.features))
        log("警告：" + message)
        return []

    suppressed = []
    for name in shot_names:
        if not repository_contains(assembly.features, name):
            log("警告：实例 %s 不在 assembly.features 中，无法执行 suppress。" % name)
            continue
        try:
            assembly.features[name].suppress()
            suppressed.append(name)
        except Exception as error:
            # 重复运行或对象已被抑制时，某些版本会抛出异常。
            log("警告：抑制弹丸实例 %s 时返回：%s" % (name, error))

    try:
        assembly.regenerate()
    except Exception:
        pass

    log("识别到弹丸实例 %d 个；成功执行抑制 %d 个。" % (len(shot_names), len(suppressed)))
    if len(suppressed) > 0:
        log("弹丸实例范围：%s ... %s" % (suppressed[0], suppressed[-1]))
    return suppressed


def find_rigid_constraint_names(model):
    """综合约束名称、候选前缀和对象类名识别全部弹丸刚体约束。"""
    candidates = []
    for name in repository_names(model.constraints):
        class_name = ""
        try:
            class_name = model.constraints[name].__class__.__name__
        except Exception:
            pass

        if starts_with_any(name, RIGID_CONSTRAINT_PREFIXES):
            candidates.append(name)
        elif "rigid" in lower_text(class_name):
            candidates.append(name)
        elif "rigid" in lower_text(name) and "shot" in lower_text(name):
            candidates.append(name)

    return sorted(set(candidates))


def suppress_rigid_constraints(model):
    """
    批量抑制所有弹丸刚体约束。

    即使弹丸实例已抑制，也应抑制引用其参考点/单元的 RigidBody 约束，
    防止 Standard 输入文件中出现悬空引用。
    """
    names = find_rigid_constraint_names(model)

    if len(names) == 0:
        message = "未识别到任何弹丸刚体约束。"
        if REQUIRE_RIGID_CONSTRAINTS:
            raise RuntimeError(message + " 当前约束：%s" % repository_names(model.constraints))
        log("警告：" + message)
        return []

    suppressed = []
    for name in names:
        try:
            model.constraints[name].suppress()
            suppressed.append(name)
        except Exception as error:
            log("警告：抑制刚体约束 %s 时返回：%s" % (name, error))

    log("识别到刚体约束 %d 个；成功执行抑制 %d 个。" % (len(names), len(suppressed)))
    return suppressed


def suppress_interactions(model):
    """
    抑制显式冲击阶段的 General Contact 等相互作用。

    静态模型只剩靶材且没有外部接触对象，继续保留 Int-1 没有物理意义。
    默认抑制全部相互作用，避免漏掉名称不同的接触。
    """
    names = []
    for name in repository_names(model.interactions):
        if SUPPRESS_ALL_INTERACTIONS or contains_any(name, INTERACTION_NAMES_OR_KEYWORDS):
            names.append(name)

    suppressed = []
    for name in names:
        try:
            model.interactions[name].suppress()
            suppressed.append(name)
        except Exception as error:
            log("警告：抑制相互作用 %s 时返回：%s" % (name, error))

    log("已抑制相互作用：%s" % suppressed)
    return suppressed


def is_velocity_predefined_field(field_object, field_name):
    """根据对象类名和候选关键词判断预定义场是否为弹丸初速度。"""
    class_name = ""
    try:
        class_name = field_object.__class__.__name__
    except Exception:
        pass

    if "velocity" in lower_text(class_name):
        return True
    if contains_any(field_name, VELOCITY_FIELD_NAMES_OR_KEYWORDS):
        return True
    return False


def suppress_velocity_fields(model):
    """
    抑制显式阶段给弹丸施加的初速度，避免它进入 Standard 输入文件。

    此函数在创建 Initial State 之前执行，不会误抑制随后新建的初始状态场。
    """
    names = []
    for name in repository_names(model.predefinedFields):
        try:
            field_object = model.predefinedFields[name]
            if is_velocity_predefined_field(field_object, name):
                names.append(name)
        except Exception:
            pass

    suppressed = []
    for name in names:
        try:
            model.predefinedFields[name].suppress()
            suppressed.append(name)
        except Exception as error:
            log("警告：抑制速度预定义场 %s 时返回：%s" % (name, error))

    log("已抑制速度预定义场：%s" % suppressed)
    return suppressed


def check_target_fixed_boundary_condition(model):
    """
    确认靶材固定边界条件存在，并保留它用于消除刚体自由度。

    若找不到固定条件仍继续计算，Standard 通常会出现数值奇异，因此默认
    REQUIRE_TARGET_FIXED_BC=True 并立即终止。
    """
    found = []
    for name in TARGET_FIXED_BC_CANDIDATES:
        if repository_contains(model.boundaryConditions, name):
            found.append(name)

    if len(found) == 0:
        # 名称不同的情况下，进一步按关键词识别。
        for name in repository_names(model.boundaryConditions):
            lower_name = lower_text(name)
            if "target" in lower_name or "encastre" in lower_name or "fixed" in lower_name:
                found.append(name)

    if len(found) == 0:
        message = (
            "未找到靶材固定边界条件。候选名称=%s，当前边界条件=%s。"
            % (TARGET_FIXED_BC_CANDIDATES, repository_names(model.boundaryConditions))
        )
        if REQUIRE_TARGET_FIXED_BC:
            raise RuntimeError(message)
        log("警告：" + message)
        return []

    log("检测到并保留靶材固定边界条件：%s" % found)
    return found


def create_initial_state(model, target_instance_name, source_odb_path):
    """
    从显式 ODB 的最后一步、最后增量导入靶材初始状态。

    InitialState 的 fileName 传入不带 .odb 的作业基名，且 Abaqus 会从当前
    工作目录查找该 ODB，所以函数先切换到 ODB 所在目录。instances 只包含
    target-1，避免把已经抑制的弹丸状态导入静态模型。

    导入内容由源 ODB 和材料模型决定，通常包括应力、塑性应变等状态变量；
    这不是简单复制位移结果，而是 Standard 后续再平衡的初始材料状态。
    """
    if not os.path.isfile(source_odb_path):
        message = "找不到显式 ODB：%s" % source_odb_path
        if REQUIRE_SOURCE_ODB:
            raise IOError(message)
        log("警告：" + message)

    # 例如 D:\work\Job-5.odb -> 目录 D:\work、作业基名 Job-5。
    odb_directory = os.path.dirname(source_odb_path)
    odb_job_name = os.path.splitext(os.path.basename(source_odb_path))[0]

    # InitialState 的 fileName 使用不带 .odb 的作业名时兼容性最好。
    # 因此切换到 ODB 所在目录，再传入作业基名。
    os.chdir(odb_directory)

    if repository_contains(model.predefinedFields, INITIAL_STATE_NAME):
        del model.predefinedFields[INITIAL_STATE_NAME]
        log("已删除旧 Initial State：%s" % INITIAL_STATE_NAME)

    assembly = model.rootAssembly
    if not repository_contains(assembly.instances, target_instance_name):
        raise KeyError("创建 Initial State 时找不到靶材实例：%s" % target_instance_name)

    target_instance = assembly.instances[target_instance_name]

    # LAST_STEP + STEP_END 精确对应“最后一步、最后增量”。
    model.InitialState(
        name=INITIAL_STATE_NAME,
        createStepName="Initial",
        fileName=odb_job_name,
        endStep=LAST_STEP,
        endIncrement=STEP_END,
        updateReferenceConfiguration=(
            ON if UPDATE_REFERENCE_CONFIGURATION else OFF
        ),
        instances=(target_instance,),
    )

    log("已创建 Initial State：%s" % INITIAL_STATE_NAME)
    log("初始状态来源：%s" % source_odb_path)
    log("导入位置：最后一步、最后增量。")
    log(
        "更新参考构型：%s"
        % ("是" if UPDATE_REFERENCE_CONFIGURATION else "否")
    )


def delete_source_analysis_steps(model):
    """
    删除静态副本中除 Initial 外的全部显式/旧分析步。

    Dynamic, Explicit 步不能直接作为 Standard 静态步的前置步；逆序删除可
    避免多步模型中后一步依赖前一步所造成的删除顺序问题。
    """
    names = [name for name in repository_names(model.steps) if name != "Initial"]

    # 逆序删除，避免多分析步之间的前后依赖造成问题。
    for name in reversed(names):
        del model.steps[name]
        log("已删除源分析步：%s" % name)

    log("删除后保留的分析步：%s" % repository_names(model.steps))


def create_static_step(model):
    """
    在 Initial 之后创建 Static, General 静态平衡分析步。

    默认增量设置最接近视频的 GUI 操作；若收敛困难，可关闭默认设置并在参数区
    减小 STATIC_INITIAL_INCREMENT / STATIC_MIN_INCREMENT。
    """
    nlgeom_value = ON if STATIC_NLGEOM else OFF

    if USE_ABAQUS_DEFAULT_INCREMENT_SETTINGS:
        model.StaticStep(
            name=STATIC_STEP_NAME,
            previous="Initial",
            nlgeom=nlgeom_value,
        )
        log("静态步采用 Abaqus 默认时间与增量设置。")
    else:
        model.StaticStep(
            name=STATIC_STEP_NAME,
            previous="Initial",
            nlgeom=nlgeom_value,
            timePeriod=STATIC_TIME_PERIOD,
            initialInc=STATIC_INITIAL_INCREMENT,
            minInc=STATIC_MIN_INCREMENT,
            maxInc=STATIC_MAX_INCREMENT,
            maxNumInc=STATIC_MAX_INCREMENTS,
        )
        log(
            "静态步参数：timePeriod=%g，initialInc=%g，minInc=%g，"
            "maxInc=%g，maxNumInc=%d。"
            % (
                STATIC_TIME_PERIOD,
                STATIC_INITIAL_INCREMENT,
                STATIC_MIN_INCREMENT,
                STATIC_MAX_INCREMENT,
                STATIC_MAX_INCREMENTS,
            )
        )

    log("已创建 Static, General：%s" % STATIC_STEP_NAME)
    log("NLGEOM：%s" % ("ON" if STATIC_NLGEOM else "OFF"))


def actual_num_domains():
    """计算实际域数；自动模式下直接返回 CPU 数，手工模式再检查整倍数。"""
    if AUTO_SET_NUM_DOMAINS:
        return int(NUM_CPUS)

    domains = int(NUM_DOMAINS)
    if domains < 1:
        raise ValueError("NUM_DOMAINS 必须为正整数。")
    if domains % int(NUM_CPUS) != 0:
        raise ValueError(
            "NUM_DOMAINS 必须是 NUM_CPUS 的整数倍。当前 NUM_CPUS=%d，"
            "NUM_DOMAINS=%d。" % (NUM_CPUS, domains)
        )
    return domains


def create_static_job():
    """
    为 STATIC_MODEL_NAME 创建 Abaqus/Standard 分析作业。

    作业名决定 INP/ODB 文件基名；因此 Job-5-Static 会输出
    Job-5-Static.inp 和 Job-5-Static.odb，并与后续 PE 提取脚本默认值衔接。
    """
    if repository_contains(mdb.jobs, STATIC_JOB_NAME):
        if REPLACE_EXISTING_STATIC_JOB:
            del mdb.jobs[STATIC_JOB_NAME]
        else:
            raise ValueError("静态作业已存在：%s" % STATIC_JOB_NAME)

    domains = actual_num_domains()

    job = mdb.Job(
        name=STATIC_JOB_NAME,
        model=STATIC_MODEL_NAME,
        description="Static equilibrium after random multi-shot peening",
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
        numDomains=domains,
        numGPUs=NUM_GPUS,
    )

    log("已创建 Abaqus/Standard 作业：%s" % STATIC_JOB_NAME)
    log("并行设置：CPU=%d，Domains=%d，GPU=%d。" % (NUM_CPUS, domains, NUM_GPUS))
    return job


def save_database(source_cae_path, output_cae_path):
    """
    根据 SAVE_AS_NEW_CAE 另存或原位保存数据库。

    另存模式会检查输出路径不能等于源路径；原位保存模式则直接保存当前打开的
    数据库，符合你的测试操作，但会把 Model-5-Static 写入当前 CAE。
    """
    if SAVE_AS_NEW_CAE:
        if os.path.normcase(os.path.abspath(output_cae_path)) == os.path.normcase(os.path.abspath(source_cae_path)):
            raise ValueError("输出 CAE 不能与源 CAE 相同。")
        mdb.saveAs(pathName=output_cae_path)
        log("模型数据库已另存为：%s" % output_cae_path)
    else:
        mdb.save()
        log("已保存当前模型数据库。")


def write_or_submit_job(job, analysis_directory):
    """
    在分析目录写入 INP，并按开关选择是否提交/等待静态作业。

    WRITE_INPUT 与 SUBMIT_JOB 相互独立：可以只写 INP，也可以写 INP 后自动
    提交。WAIT_FOR_COMPLETION=True 时脚本会阻塞到求解结束，再保存 CAE。
    """
    os.chdir(analysis_directory)

    if WRITE_INPUT:
        job.writeInput(consistencyChecking=CONSISTENCY_CHECKING)
        log("已写出输入文件：%s" % os.path.join(analysis_directory, STATIC_JOB_NAME + ".inp"))

    if SUBMIT_JOB:
        log("正在提交静态作业：%s" % STATIC_JOB_NAME)
        job.submit(consistencyChecking=CONSISTENCY_CHECKING)
        if WAIT_FOR_COMPLETION:
            job.waitForCompletion()
            log("静态作业结束，状态：%s" % job.status)
            mdb.save()
    else:
        log("SUBMIT_JOB=False：未自动提交计算。")


def display_static_model(model):
    """图形模式下把当前视口切换到静态模型；无图形模式失败时仅记录警告。"""
    if not DISPLAY_STATIC_MODEL:
        return

    try:
        if len(session.viewports.keys()) > 0:
            viewport_name = session.currentViewportName
            session.viewports[viewport_name].setValues(
                displayedObject=model.rootAssembly
            )
            log("当前视口已切换到静态模型装配。")
    except Exception as error:
        log("警告：切换视口显示失败：%s" % error)


def write_summary_log(output_directory):
    """
    把 _ACTION_LOG 写为带 UTF-8 BOM 的中文 TXT。

    BOM 使 Windows 记事本可直接识别中文编码；写日志失败不会撤销已完成的模型。
    """
    if not WRITE_SUMMARY_LOG:
        return

    path = absolute_path(output_directory, SUMMARY_LOG_FILE_NAME)
    try:
        # UTF-8 BOM，便于 Windows 记事本直接正确显示中文。
        with open(path, "wb") as handle:
            text = "\r\n".join(_ACTION_LOG) + "\r\n"
            handle.write(text.encode("utf-8-sig"))
        print("[RandomShotStatic] 摘要日志：%s" % path)
    except Exception as error:
        print("[RandomShotStatic] 警告：写摘要日志失败：%s" % error)


# =============================================================================
# 三、主流程
# =============================================================================


def main():
    """
    执行随机喷丸显式模型到静态模型的八阶段转换。

    各阶段顺序具有依赖关系：必须先复制模型，再抑制弹丸/显式对象；必须在删除
    显式步之前创建 Initial State；静态步和作业只能在清理完成后创建。
    """
    print_title("Abaqus/Standard 随机多弹丸喷丸后静态分析")

    # 阶段 0：只做参数与路径准备，尚未修改模型数据库。
    validate_parameters()
    work_directory = resolve_work_directory()

    # 统一转换为绝对路径，后续切换工作目录时不会丢失文件定位。
    source_cae_path = absolute_path(work_directory, SOURCE_CAE_FILE_NAME)
    source_odb_path = absolute_path(work_directory, SOURCE_ODB_FILE_NAME)
    output_cae_name = derive_output_cae_name()
    output_cae_path = absolute_path(work_directory, output_cae_name)
    analysis_directory = os.path.dirname(source_odb_path)

    log("源 CAE：%s" % source_cae_path)
    log("源 ODB：%s" % source_odb_path)
    log("输出 CAE：%s" % output_cae_path)
    log("静态作业目录：%s" % analysis_directory)

    if REQUIRE_SOURCE_ODB and not os.path.isfile(source_odb_path):
        raise IOError("找不到源 ODB：%s" % source_odb_path)

    print_title("1/8 打开源 CAE 并复制显式模型")
    # 复制是安全边界：从这一行以后只修改 static_model，不修改 source model。
    open_or_get_database(source_cae_path)
    source_model_name = resolve_source_model_name()
    delete_existing_static_objects()
    static_model = copy_source_model(source_model_name)
    target_instance_name = resolve_target_instance_name(static_model)

    print_title("2/8 抑制全部随机弹丸实例")
    # 从静态装配中移除全部 shot，只让 target 进入 Standard 作业。
    suppress_shot_instances(static_model, target_instance_name)

    print_title("3/8 抑制显式阶段对象")
    # 刚体、接触和初速度属于冲击阶段；靶材固定边界条件则必须保留。
    suppress_rigid_constraints(static_model)
    suppress_interactions(static_model)
    suppress_velocity_fields(static_model)
    check_target_fixed_boundary_condition(static_model)

    print_title("4/8 从显式 ODB 导入靶材 Initial State")
    # 静态再平衡的核心：把 Job-5.odb 末态映射到 target-1。
    create_initial_state(static_model, target_instance_name, source_odb_path)

    print_title("5/8 删除显式分析步并创建 Static, General")
    # Initial State 已建立后，才能安全移除 Dynamic, Explicit 步并创建静态步。
    delete_source_analysis_steps(static_model)
    create_static_step(static_model)

    print_title("6/8 创建 Abaqus/Standard 静态作业")
    # 作业只引用 Model-5-Static，最终结果名为 Job-5-Static.odb。
    static_job = create_static_job()

    print_title("7/8 保存静态 CAE")
    # 当前默认 SAVE_AS_NEW_CAE=False，因此保存到当前已经打开的 CAE。
    save_database(source_cae_path, output_cae_path)
    display_static_model(static_model)

    print_title("8/8 写出 INP 或提交计算")
    # 默认只写 Job-5-Static.inp；确认无误后可把 SUBMIT_JOB 改为 True。
    write_or_submit_job(static_job, analysis_directory)

    log("随机喷丸静态分析模型创建完成。")
    log("静态模型：%s" % STATIC_MODEL_NAME)
    log("静态分析步：%s" % STATIC_STEP_NAME)
    log("静态作业：%s" % STATIC_JOB_NAME)

    write_summary_log(work_directory)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print("\n" + "=" * 78)
        print("脚本执行失败：%s" % error)
        print("以下为错误追踪信息：")
        print(traceback.format_exc())
        print("=" * 78)

        _ACTION_LOG.append("脚本执行失败：%s" % error)
        _ACTION_LOG.append(traceback.format_exc())
        try:
            directory = os.getcwd()
            write_summary_log(directory)
        except Exception:
            pass
        raise
