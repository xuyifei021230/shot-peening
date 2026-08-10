# -*- coding: utf-8 -*-
"""
Abaqus 随机喷丸静态结果分层 PE 提取脚本（连接 Almen 固有应变分析）
=====================================================================

脚本在三段流程中的位置
----------------------
第一段：Random_Shot_Peening_Full_Auto_V5.1.py
    创建随机多弹丸显式模型和 SET-LAYER-01～SET-LAYER-30 分层单元集，
    求解后得到 Job-5.odb（实际作业名可在参数区修改）。

第二段：random_shot_static_elastic_rebuild_abaqus2025_V2.py
    把 Job-5.odb 的最后状态导入 Model-5-Static，抑制弹丸、刚体和接触，
    进行 Abaqus/Standard 静态再平衡，得到 Job-5-Static.odb。

第三段：本脚本
    从 Job-5-Static.odb 的 Step-1 最后一帧读取每个 SET-LAYER-xx，计算
    PE11、PE22、PE33 的算术平均值，并输出：

        Job-5-Static-PE.txt
            必然生成。
            无表头、三列，列顺序固定为 PE11 PE22 PE33；可直接供旧工程
            “赋值.py”中的 numpy.loadtxt() 读取，并作为 Almen 分层材料的
            ORTHOTROPIC Expansion 三个方向系数。

        Job-5-Static-PE-details.csv
            可通过用户参数选择是否生成。
            带层号、深度、集合名、数据点数量和 PE11+PE22+PE33 校核列，
            用于人工核查，不作为 Almen 赋值脚本的直接输入。

        Job-5-Static-PE-summary.txt
            可通过用户参数选择是否生成，记录本次提取的运行信息和逐层摘要。

为什么默认读取静态 ODB
----------------------
旧工程的“提取PE.py”读取 Job-4.odb；从 Job-4.inp 可确认 Job-4 是通过
*Import, state=yes, update=no 建立的 Static 分析，而不是原始 Explicit 作业。
因此当前连贯流程应读取第二段产生的 Job-5-Static.odb 最后一帧。

平均方法
--------
脚本保持旧“提取PE.py”的算法：对集合内 PE 场的全部数值做算术平均。
当前靶材主要使用 C3D8R 单元，每个单元一个积分点，因此等价于对层内单元值
等权平均。若以后改成多积分点单元或需要严格体积加权，应另行输出 IVOL 并改用
积分体积加权，不能把本结果误称为体积加权平均。

运行方式
--------
推荐在 Abaqus/CAE 中：File -> Run Script -> 选择本文件。
也可以在命令行使用：abaqus python extract_random_shot_PE_for_Almen.py

脚本兼容 Abaqus 2022 Python 2.7 与 Abaqus 2025 Python 3，不使用 f-string、
类型注解或依赖 NumPy 的文件写出功能。
"""

from __future__ import print_function

from odbAccess import openOdb

import io
import os
import traceback


# =============================================================================
# 一、用户参数区
# =============================================================================

# True：使用运行脚本时 Abaqus 的当前工作目录。
# False：使用 WORK_DIRECTORY。建议与前两个随机喷丸脚本使用同一工作目录。
USE_CURRENT_WORK_DIRECTORY = True
WORK_DIRECTORY = r"D:\temp\random-shot"

# 第二段静态脚本默认作业名为 Job-5-Static，因此静态结果文件为下列名称。
# 也可以填写绝对路径；绝对路径不再与工作目录拼接。
SOURCE_ODB_FILE_NAME = "Job-1-Static.odb"

# 优先读取指定分析步。若名称不存在且 AUTO_USE_LAST_STEP=True，自动采用
# ODB 中最后一个分析步，并在日志中明确输出实际名称。
SOURCE_STEP_NAME = "Step-1"
AUTO_USE_LAST_STEP = True

# LAST：最后一帧，代表静态再平衡完成后的结果（推荐）。
# FIRST：该分析步的第 0 帧；INDEX：使用 FRAME_INDEX，可填负数。
FRAME_SELECTION = "LAST"
FRAME_INDEX = -1

# Abaqus 塑性应变张量场及需要输出的三个正应变分量。
FIELD_OUTPUT_NAME = "PE"
PE_COMPONENTS = ("PE11", "PE22", "PE33")

# 与第一段随机喷丸脚本保持一致：SET-LAYER-01～SET-LAYER-30。
LAYER_SET_PREFIX = "SET-LAYER-"
LAYER_SET_COUNT = 30
LAYER_INDEX_WIDTH = 2
LAYER_THICKNESS_MM = 0.01

# 第一段脚本在装配层创建分层集合，通常可直接从 rootAssembly.elementSets
# 读取。若找不到，脚本还会尝试 target-1 实例内部的同名集合。
TARGET_INSTANCE_NAME = "target-1"
AUTO_DETECT_TARGET_INSTANCE = True
SEARCH_INSTANCE_ELEMENT_SETS = True

# None：根据 ODB 基名自动生成 Job-5-Static-PE.txt / -PE-details.csv / -PE-summary.txt。
# 也可填写文件名或绝对路径。OUTPUT_TEXT_FILE_NAME 对应的 TXT 必然生成，
# 且必须保持无表头三列，才能直接衔接赋值.py。
OUTPUT_TEXT_FILE_NAME = None
OUTPUT_CSV_FILE_NAME = None
OUTPUT_SUMMARY_FILE_NAME = None

# True：生成对应辅助文件；False：不生成。主 TXT 不设置开关，始终生成。
# 默认保持旧脚本行为，即两个辅助文件均生成。
GENERATE_CSV_FILE = False
GENERATE_SUMMARY_FILE = False

# 旧脚本使用 %.10f；保持一致可以直接比较新旧结果。
OUTPUT_FLOAT_FORMAT = "%.10f"

# 提取结束后关闭脚本打开的 ODB，避免同名 ODB 文件持续被 CAE 占用。
CLOSE_ODB_AFTER_EXTRACTION = False


# =============================================================================
# 二、通用工具函数
# =============================================================================

_LOG_LINES = []


def log(message):
    """输出统一日志并缓存，最后写入摘要文件。"""
    text = str(message)
    print("[ExtractPE] " + text)
    _LOG_LINES.append(text)


def repository_names(repository):
    """把 Abaqus Repository 的键转换为普通列表。"""
    try:
        return list(repository.keys())
    except Exception:
        return []


def repository_contains(repository, name):
    """通过直接索引检查 Repository 是否含指定键。"""
    try:
        repository[name]
        return True
    except Exception:
        return False


def resolve_key_case_insensitive(repository, requested_name):
    """优先精确匹配，再进行不区分大小写匹配；找不到时返回 None。"""
    if repository_contains(repository, requested_name):
        return requested_name

    requested_lower = str(requested_name).lower()
    for name in repository_names(repository):
        if str(name).lower() == requested_lower:
            return name
    return None


def absolute_path(base_directory, file_name):
    """允许参数使用文件名或绝对路径。"""
    text = str(file_name)
    if os.path.isabs(text):
        return os.path.abspath(text)
    return os.path.abspath(os.path.join(base_directory, text))


def resolve_work_directory():
    """确定并校验工作目录。"""
    if USE_CURRENT_WORK_DIRECTORY:
        directory = os.getcwd()
        source = "Abaqus 当前工作目录"
    else:
        directory = WORK_DIRECTORY
        source = "WORK_DIRECTORY"

    directory = os.path.abspath(str(directory))
    if not os.path.isdir(directory):
        raise IOError("工作目录不存在：%s" % directory)

    log("工作目录：%s" % directory)
    log("工作目录来源：%s" % source)
    return directory


def validate_parameters():
    """在打开大型 ODB 以前检查参数，避免无效读取。"""
    if LAYER_SET_COUNT < 1:
        raise ValueError("LAYER_SET_COUNT 必须大于 0。")
    if LAYER_INDEX_WIDTH < 1:
        raise ValueError("LAYER_INDEX_WIDTH 必须大于 0。")
    if LAYER_THICKNESS_MM <= 0.0:
        raise ValueError("LAYER_THICKNESS_MM 必须大于 0。")
    if len(PE_COMPONENTS) != 3:
        raise ValueError("PE_COMPONENTS 必须恰好包含三个分量。")
    if len(set(PE_COMPONENTS)) != len(PE_COMPONENTS):
        raise ValueError("PE_COMPONENTS 中存在重复分量。")
    if not isinstance(GENERATE_CSV_FILE, bool):
        raise ValueError("GENERATE_CSV_FILE 必须是 True 或 False。")
    if not isinstance(GENERATE_SUMMARY_FILE, bool):
        raise ValueError("GENERATE_SUMMARY_FILE 必须是 True 或 False。")

    mode = str(FRAME_SELECTION).upper()
    if mode not in ("FIRST", "LAST", "INDEX"):
        raise ValueError("FRAME_SELECTION 只能是 FIRST、LAST 或 INDEX。")

    try:
        OUTPUT_FLOAT_FORMAT % 0.0
    except Exception:
        raise ValueError("OUTPUT_FLOAT_FORMAT 不是有效的单浮点数格式。")


def default_output_paths(work_directory, source_odb_path):
    """根据 ODB 基名生成三个默认输出路径。"""
    odb_base = os.path.splitext(os.path.basename(source_odb_path))[0]

    text_name = OUTPUT_TEXT_FILE_NAME
    if text_name is None:
        text_name = odb_base + "-PE.txt"

    csv_name = OUTPUT_CSV_FILE_NAME
    if csv_name is None:
        csv_name = odb_base + "-PE-details.csv"

    summary_name = OUTPUT_SUMMARY_FILE_NAME
    if summary_name is None:
        summary_name = odb_base + "-PE-summary.txt"

    return (
        absolute_path(work_directory, text_name),
        absolute_path(work_directory, csv_name),
        absolute_path(work_directory, summary_name),
    )


def layer_set_name(layer_number):
    """把第 1 层转换为 SET-LAYER-01 形式。"""
    return LAYER_SET_PREFIX + str(int(layer_number)).zfill(LAYER_INDEX_WIDTH)


def resolve_step(odb):
    """取得指定分析步；必要时自动退回 ODB 中最后一个分析步。"""
    key = resolve_key_case_insensitive(odb.steps, SOURCE_STEP_NAME)
    if key is not None:
        log("采用指定分析步：%s" % key)
        return odb.steps[key], key

    step_names = repository_names(odb.steps)
    if AUTO_USE_LAST_STEP and len(step_names) > 0:
        key = step_names[-1]
        log("指定分析步 %s 不存在，自动采用最后一步：%s" % (SOURCE_STEP_NAME, key))
        return odb.steps[key], key

    raise KeyError(
        "找不到分析步 %s。ODB 中的分析步：%s" % (SOURCE_STEP_NAME, step_names)
    )


def resolve_frame(step):
    """根据 FRAME_SELECTION 取得实际帧和帧索引。"""
    frames = step.frames
    if len(frames) == 0:
        raise RuntimeError("所选分析步不包含任何结果帧。")

    mode = str(FRAME_SELECTION).upper()
    if mode == "FIRST":
        index = 0
    elif mode == "LAST":
        index = len(frames) - 1
    else:
        index = int(FRAME_INDEX)
        if index < 0:
            index = len(frames) + index

    if index < 0 or index >= len(frames):
        raise IndexError(
            "FRAME_INDEX 超出范围：实际索引=%d，帧数=%d。" % (index, len(frames))
        )

    frame = frames[index]
    log("采用结果帧：%d/%d，frameValue=%g" % (index, len(frames) - 1, frame.frameValue))
    return frame, index


def resolve_target_instance(root_assembly):
    """为实例级集合回退搜索确定 target-1 实例。"""
    key = resolve_key_case_insensitive(root_assembly.instances, TARGET_INSTANCE_NAME)
    if key is not None:
        return root_assembly.instances[key], key

    if AUTO_DETECT_TARGET_INSTANCE:
        candidates = []
        for name in repository_names(root_assembly.instances):
            if "target" in str(name).lower():
                candidates.append(name)
        if len(candidates) == 1:
            key = candidates[0]
            log("指定靶材实例不存在，自动识别为：%s" % key)
            return root_assembly.instances[key], key

    return None, None


def resolve_layer_element_set(root_assembly, target_instance, requested_name):
    """
    先查装配级 elementSets，再查靶材实例 elementSets。

    第一段脚本使用 assembly.Set()，所以正常情况下会在装配级命中。
    """
    key = resolve_key_case_insensitive(root_assembly.elementSets, requested_name)
    if key is not None:
        return root_assembly.elementSets[key], "ASSEMBLY/%s" % key

    if SEARCH_INSTANCE_ELEMENT_SETS and target_instance is not None:
        key = resolve_key_case_insensitive(target_instance.elementSets, requested_name)
        if key is not None:
            return target_instance.elementSets[key], "INSTANCE/%s" % key

    return None, None


def available_layer_like_sets(root_assembly, target_instance):
    """列出与分层前缀相似的集合，帮助定位命名不一致。"""
    prefix_lower = str(LAYER_SET_PREFIX).lower()
    found = []

    for name in repository_names(root_assembly.elementSets):
        if str(name).lower().startswith(prefix_lower):
            found.append("ASSEMBLY/%s" % name)

    if target_instance is not None:
        for name in repository_names(target_instance.elementSets):
            if str(name).lower().startswith(prefix_lower):
                found.append("INSTANCE/%s" % name)

    return found


def get_field_output(frame):
    """取得 PE 张量场，并在缺失时列出当前帧可用变量。"""
    key = resolve_key_case_insensitive(frame.fieldOutputs, FIELD_OUTPUT_NAME)
    if key is None:
        raise KeyError(
            "当前帧没有 %s 场输出。可用变量：%s。请在前序模型中输出 PE。"
            % (FIELD_OUTPUT_NAME, repository_names(frame.fieldOutputs))
        )
    log("采用场输出：%s" % key)
    return frame.fieldOutputs[key]


def arithmetic_component_average(field_subset, component_label):
    """提取一个张量分量并对全部 FieldValue.data 做算术平均。"""
    try:
        scalar_field = field_subset.getScalarField(componentLabel=component_label)
    except Exception as error:
        raise RuntimeError(
            "无法从 %s 提取分量 %s：%s" % (FIELD_OUTPUT_NAME, component_label, error)
        )

    values = scalar_field.values
    if len(values) == 0:
        raise RuntimeError("分量 %s 在当前分层集合中没有数值。" % component_label)

    total = 0.0
    for value in values:
        total += float(value.data)
    return total / float(len(values)), len(values)


# =============================================================================
# 三、分层提取与文件写出
# =============================================================================

def extract_layer_rows(odb, frame):
    """依次提取 30 个分层集合，返回可写入 TXT/CSV 的记录列表。"""
    root_assembly = odb.rootAssembly
    target_instance, target_instance_key = resolve_target_instance(root_assembly)
    if target_instance_key is not None:
        log("靶材实例：%s" % target_instance_key)
    else:
        log("警告：未识别靶材实例；仅搜索装配级分层集合。")

    pe_field = get_field_output(frame)
    rows = []

    for layer_number in range(1, LAYER_SET_COUNT + 1):
        requested_set_name = layer_set_name(layer_number)
        element_set, resolved_location = resolve_layer_element_set(
            root_assembly, target_instance, requested_set_name
        )

        if element_set is None:
            raise KeyError(
                "找不到第 %d 层集合 %s。相似集合：%s"
                % (
                    layer_number,
                    requested_set_name,
                    available_layer_like_sets(root_assembly, target_instance),
                )
            )

        # 一次取得该层 PE 子集，再从同一子集中分解三个正应变分量。
        field_subset = pe_field.getSubset(region=element_set)
        component_averages = []
        component_counts = []
        for component_label in PE_COMPONENTS:
            average, count = arithmetic_component_average(
                field_subset, component_label
            )
            component_averages.append(average)
            component_counts.append(count)

        if len(set(component_counts)) != 1:
            raise RuntimeError(
                "%s 三个分量的数据数量不一致：%s"
                % (requested_set_name, component_counts)
            )

        depth_top = (layer_number - 1) * LAYER_THICKNESS_MM
        depth_bottom = layer_number * LAYER_THICKNESS_MM
        depth_mid = 0.5 * (depth_top + depth_bottom)
        trace_value = sum(component_averages)

        row = {
            "layer": layer_number,
            "set_name": requested_set_name,
            "resolved_location": resolved_location,
            "depth_top_mm": depth_top,
            "depth_mid_mm": depth_mid,
            "depth_bottom_mm": depth_bottom,
            "value_count": component_counts[0],
            "PE11": component_averages[0],
            "PE22": component_averages[1],
            "PE33": component_averages[2],
            "trace": trace_value,
        }
        rows.append(row)

        log(
            "%s：N=%d，PE11=% .6e，PE22=% .6e，PE33=% .6e，trace=% .3e"
            % (
                requested_set_name,
                row["value_count"],
                row["PE11"],
                row["PE22"],
                row["PE33"],
                row["trace"],
            )
        )

    return rows


def write_almen_text(rows, output_path):
    """写出赋值.py 可直接 numpy.loadtxt() 的无表头三列 PE 文本。"""
    with io.open(output_path, mode="w", encoding="ascii", newline="") as handle:
        for row in rows:
            line = " ".join(
                (
                    OUTPUT_FLOAT_FORMAT % row["PE11"],
                    OUTPUT_FLOAT_FORMAT % row["PE22"],
                    OUTPUT_FLOAT_FORMAT % row["PE33"],
                )
            )
            handle.write(line + "\n")
    log("已写出 Almen 三列 PE 文本：%s" % output_path)


def write_detail_csv(rows, output_path):
    """写出带深度、数量和 trace 校核信息的 UTF-8 CSV。"""
    header = (
        "layer,set_name,resolved_location,depth_top_mm,depth_mid_mm,"
        "depth_bottom_mm,value_count,PE11,PE22,PE33,PE_trace"
    )
    with io.open(output_path, mode="w", encoding="utf-8-sig", newline="") as handle:
        handle.write(header + "\n")
        for row in rows:
            handle.write(
                "%d,%s,%s,%.10f,%.10f,%.10f,%d,%.10e,%.10e,%.10e,%.10e\n"
                % (
                    row["layer"],
                    row["set_name"],
                    row["resolved_location"],
                    row["depth_top_mm"],
                    row["depth_mid_mm"],
                    row["depth_bottom_mm"],
                    row["value_count"],
                    row["PE11"],
                    row["PE22"],
                    row["PE33"],
                    row["trace"],
                )
            )
    log("已写出 PE 详细 CSV：%s" % output_path)


def write_summary(output_path):
    """把运行选择和每层摘要写成便于记事本查看的中文日志。"""
    with io.open(output_path, mode="w", encoding="utf-8-sig", newline="") as handle:
        handle.write(u"\r\n".join(_LOG_LINES) + u"\r\n")
    print("[ExtractPE] 已写出摘要日志：%s" % output_path)


# =============================================================================
# 四、主流程
# =============================================================================

def main():
    """打开静态 ODB、提取 30 层 PE，并按开关写出结果文件。"""
    print("=" * 78)
    print("随机喷丸静态结果 -> Almen 固有应变：分层 PE 提取")
    print("=" * 78)

    validate_parameters()
    work_directory = resolve_work_directory()
    source_odb_path = absolute_path(work_directory, SOURCE_ODB_FILE_NAME)
    text_path, csv_path, summary_path = default_output_paths(
        work_directory, source_odb_path
    )

    log("源静态 ODB：%s" % source_odb_path)
    log("目标集合：%s01 ～ %s%s" % (
        LAYER_SET_PREFIX,
        LAYER_SET_PREFIX,
        str(LAYER_SET_COUNT).zfill(LAYER_INDEX_WIDTH),
    ))
    log("分量顺序：%s" % (PE_COMPONENTS,))
    log("平均方法：集合内 FieldValue 算术平均")
    log("详细 CSV：%s" % ("生成" if GENERATE_CSV_FILE else "不生成"))
    log("摘要日志：%s" % ("生成" if GENERATE_SUMMARY_FILE else "不生成"))

    if not os.path.isfile(source_odb_path):
        raise IOError("找不到静态 ODB：%s" % source_odb_path)

    odb = None
    try:
        odb = openOdb(path=source_odb_path, readOnly=True)
        log("ODB 已打开。")

        step, step_name = resolve_step(odb)
        frame, frame_index = resolve_frame(step)
        log("实际读取：step=%s，frame=%d" % (step_name, frame_index))

        rows = extract_layer_rows(odb, frame)
        if len(rows) != LAYER_SET_COUNT:
            raise RuntimeError(
                "提取层数不完整：%d/%d。" % (len(rows), LAYER_SET_COUNT)
            )

        write_almen_text(rows, text_path)
        if GENERATE_CSV_FILE:
            write_detail_csv(rows, csv_path)
        else:
            log("GENERATE_CSV_FILE=False：未生成 PE 详细 CSV。")

        log("分层 PE 提取完成：%d 层。" % len(rows))
        if GENERATE_SUMMARY_FILE:
            write_summary(summary_path)
        else:
            print("[ExtractPE] GENERATE_SUMMARY_FILE=False：未生成摘要日志。")

    finally:
        if odb is not None and CLOSE_ODB_AFTER_EXTRACTION:
            try:
                odb.close()
                print("[ExtractPE] ODB 已关闭。")
            except Exception as error:
                print("[ExtractPE] 警告：关闭 ODB 失败：%s" % error)

    print("=" * 78)
    print("输出给 Almen 赋值脚本的文件：%s" % text_path)
    if GENERATE_CSV_FILE:
        print("PE 详细 CSV：%s" % csv_path)
    if GENERATE_SUMMARY_FILE:
        print("摘要日志：%s" % summary_path)
    print("=" * 78)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print("")
        print("=" * 78)
        print("PE 提取失败：%s" % error)
        print("错误追踪：")
        traceback.print_exc()
        print("=" * 78)
        raise
