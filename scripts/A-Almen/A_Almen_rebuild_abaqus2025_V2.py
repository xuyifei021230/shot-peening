#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A-Almen 阿尔门试片喷丸等效热膨胀法建模全流程脚本

功能：
    1. 新建/重建 Abaqus 模型 Model-1；
    2. 创建 A-Almen 试片：38 mm × 9.5 mm × 1.29 mm；
    3. 在距离上表面 0.30 mm 处进行分区；
    4. 网格划分：
        - X 方向约 1 mm；
        - Y 方向约 1 mm；
        - 表层 0.30 mm 区域沿厚度方向每层 0.01 mm；
        - 下部基体区域沿厚度方向约 0.1 mm；
    5. 创建 set-1 ~ set-30，每个集合对应表层 0.01 mm；
    6. 创建 Set-0，对应未喷丸影响的下部基体区域；
    7. 为 Set-0 赋普通弹性材料；
    8. 为 set-1 ~ set-30 读取/使用 80ms-PE 数据并赋予正交热膨胀系数；
    9. 创建装配、分析步、边界条件和温度预定义场；
    10. 创建 Job-1。

运行方式：
    在 Abaqus/CAE 中：
        File -> Run Script -> 选择本脚本

重要说明：
    - 本脚本是根据 Abaqus/CAE 的 .rpy 操作记录整理的“稳定复现版”。
    - 已删除大量视角调整、显示开关、手动点击 mask 等不稳定命令。
    - 关键几何/网格/材料参数集中在“用户参数区”，后续换试片或换喷丸数据时只改那里即可。
    - 脚本默认不自动提交计算，只创建 Job；需要计算时可在脚本末尾取消 submit_job() 的注释。
"""

from __future__ import print_function

from abaqus import *
from abaqusConstants import *
from caeModules import *

import os
import mesh
import regionToolset


# =============================================================================
# 0. 用户参数区：以后最常改的地方都在这里
# =============================================================================

# ---------- 模型、零件、实例、作业名称 ----------
MODEL_NAME = 'Model-1'
PART_NAME = 'A-Almen'
INSTANCE_NAME = 'A-Almen-1'
JOB_NAME = 'Job-1'

# 是否重建模型。
# True ：每次运行都会删除已有 Model-1，然后重新建模，最适合“完整复现视频操作”。
# False：若 Model-1 已存在，则继续在已有模型上操作，适合二次修改。
RECREATE_MODEL = False

# ---------- 几何尺寸，单位要与你 Abaqus 模型单位一致 ----------
# 视频中 A 型阿尔门试片尺寸为 38 mm × 9.5 mm × 1.29 mm。
LENGTH_X = 38.0      # X 方向长度
WIDTH_Y = 9.5        # Y 方向宽度
THICKNESS_Z = 1.29   # 总厚度

# 喷丸影响层深度。
# 视频中取上表面向下 0.30 mm。
AFFECTED_DEPTH = 0.30

# 分区平面的 Z 坐标。
# 默认底面 z=0，上表面 z=1.29，所以分区位置 z=1.29-0.30=0.99。
PARTITION_Z = THICKNESS_Z - AFFECTED_DEPTH

# ---------- 网格参数 ----------
GLOBAL_MESH_SIZE = 1.0        # 整体种子尺寸，主要控制 X/Y 方向
TOP_LAYER_SIZE_Z = 0.01       # 喷丸影响层厚度方向网格尺寸
BASE_LAYER_SIZE_Z = 0.10      # 下部基体厚度方向网格尺寸

# 表层分层数量。
# 0.30 mm / 0.01 mm = 30 层。
LAYER_NUMBER = 30

# 框选单元、查找几何时使用的容差。
TOL = 1.0e-4

# ---------- 材料参数 ----------
# 如果你的单位是 N-mm-MPa，205000 表示 205000 MPa。
ELASTIC_MODULUS = 205000.0
POISSON_RATIO = 0.29

BASE_MATERIAL_NAME = 'Material-SAE'
BASE_SECTION_NAME = 'Section-SAE-0'

LAYER_MATERIAL_PREFIX = 'SAE-'
LAYER_SECTION_PREFIX = 'Section-SAE-'
LAYER_SET_PREFIX = 'set-'
BASE_SET_NAME = 'Set-0'

# ---------- PE 数据文件 ----------
# 本脚本要求必须读取外部 PE 数据文件。
# 如果找不到 PE_DATA_FILE，脚本会在建模前直接终止，并打印已经搜索过的路径。
#
# 推荐做法：
#   1. 将 80ms-PE.txt 与本脚本放在同一个文件夹；
#   2. 或者把 PE_DATA_FILE 改成绝对路径，例如：
#      PE_DATA_FILE = r'D:/AbaqusWorks/A_shipian/80ms-PE.txt'
#
# 文件格式要求：
#   每行 3 个数字，分别对应 ORTHOTROPIC 热膨胀系数的 1、2、3 方向；
#   共至少 LAYER_NUMBER 行。
PE_DATA_FILE = '80ms-PE.txt'


# ---------- 边界条件参数 ----------
# 视频中使用了 XsymmBC、YsymmBC 和一个 EncastreBC 顶点约束。
# 默认脚本将：
#   XsymmBC 施加在 x=0 的面；
#   YsymmBC 施加在 y=0 的面；
#   EncastreBC 施加在 (0, 0, 0) 顶点。
#
# 如果你的模型选面和视频相反，可改成：
#   XSYMM_FACE = 'X_MAX'
#   YSYMM_FACE = 'Y_MAX'
XSYMM_FACE = 'X_MIN'        # 可选：'X_MIN' 或 'X_MAX'
YSYMM_FACE = 'Y_MIN'        # 可选：'Y_MIN' 或 'Y_MAX'
FIX_VERTEX_COORD = (0.0, 0.0, 0.0)

# 温度场：Initial 中温度为 0，Step-1 中温度变为 1。
# 这样热膨胀系数数值就直接等效为应变量。
INITIAL_TEMPERATURE = 0.0
STEP_TEMPERATURE = 1.0

# ---------- 分析步和 Job 参数 ----------
STEP_NAME = 'Step-1'
NUM_CPUS = 16

# ---------- CAE 保存参数 ----------
# 是否自动保存 CAE。
SAVE_CAE = True

# 保存 CAE 文件名。
# 以后只想改保存出来的 .cae 名称，主要改这里即可。
# 例如：CAE_SAVE_FILE_NAME = 'A_shipian_80ms.cae'
CAE_SAVE_FILE_NAME = 'A-Almen.cae'

# 保存 CAE 的文件夹。
# 为空字符串 '' 时，默认保存到 Abaqus 当前工作目录。
# 也可以改成绝对路径，例如：
# CAE_SAVE_DIRECTORY = r'D:/AbaqusWorks/A_shipian'
CAE_SAVE_DIRECTORY = ''

# 如果你更习惯直接写完整路径，也可以在这里填写完整 .cae 路径。
# 注意：CAE_SAVE_PATH 非空时，会优先使用它，忽略上面的 CAE_SAVE_DIRECTORY 和 CAE_SAVE_FILE_NAME。
# 例如：CAE_SAVE_PATH = r'D:/AbaqusWorks/A_shipian/A_shipian.cae'
CAE_SAVE_PATH = ''


# =============================================================================
# 1. 通用工具函数
# =============================================================================

def info(msg):
    """统一打印信息，便于在 Abaqus Message Area 中查看进度。"""
    print('[INFO] ' + str(msg))


def warn(msg):
    """统一打印警告。"""
    print('[WARNING] ' + str(msg))


def delete_if_exists(repo, name):
    """如果 Abaqus 仓库 repo 中存在 name，则删除，避免重复运行时报重名错误。"""
    if name in repo.keys():
        del repo[name]


def safe_make_dir(path):
    """如果路径不存在则创建。"""
    if path and (not os.path.isdir(path)):
        os.makedirs(path)


def reset_model(model_name):
    """
    新建或重建模型。
    Abaqus 至少需要保留一个模型，因此删除唯一模型时先创建一个临时模型。
    """
    if RECREATE_MODEL and model_name in mdb.models.keys():
        temp_name = '__temp_model_for_delete__'
        if len(mdb.models.keys()) == 1:
            mdb.Model(name=temp_name)
        del mdb.models[model_name]
        mdb.Model(name=model_name)
        if temp_name in mdb.models.keys() and len(mdb.models.keys()) > 1:
            del mdb.models[temp_name]
        info('Recreated model: ' + model_name)
    elif model_name not in mdb.models.keys():
        mdb.Model(name=model_name)
        info('Created model: ' + model_name)
    else:
        info('Use existing model: ' + model_name)

    return mdb.models[model_name]


def resolve_file(file_name):
    """
    查找外部 PE 数据文件。

    返回：
        found_path, searched_paths

    found_path:
        找到文件时为完整路径；找不到时为 None。

    searched_paths:
        实际搜索过的路径列表，用于报错提示。

    查找顺序：
        1. 绝对路径；
        2. 脚本所在目录；
        3. Abaqus 当前工作目录；
        4. 当前 CAE 文件所在目录。
    """
    candidate_paths = []

    if os.path.isabs(file_name):
        candidate_paths.append(file_name)
    else:
        if '__file__' in globals():
            script_dir = os.path.dirname(os.path.abspath(__file__))
            candidate_paths.append(os.path.join(script_dir, file_name))

        candidate_paths.append(os.path.join(os.getcwd(), file_name))

        try:
            if mdb.pathName:
                cae_dir = os.path.dirname(os.path.abspath(mdb.pathName))
                candidate_paths.append(os.path.join(cae_dir, file_name))
        except Exception:
            pass

    # 去掉重复路径，避免提示信息太乱
    unique_paths = []
    for path in candidate_paths:
        if path and path not in unique_paths:
            unique_paths.append(path)

    for path in unique_paths:
        if path and os.path.isfile(path):
            return path, unique_paths

    return None, unique_paths


def read_pe_data_from_txt(file_path):
    """
    读取 PE 数据 txt 文件。
    文件格式要求：
        每行 3 个数字，可以用空格、Tab 或英文逗号分隔。
    """
    data = []

    f = open(file_path, 'r')
    try:
        line_number = 0
        for raw_line in f:
            line_number += 1
            line = raw_line.strip()

            if not line:
                continue
            if line.startswith('#'):
                continue

            line = line.replace(',', ' ')
            items = line.split()

            if len(items) < 3:
                raise RuntimeError(
                    'PE file line %d has fewer than 3 numbers: %s' %
                    (line_number, line)
                )

            try:
                row = (float(items[0]), float(items[1]), float(items[2]))
            except ValueError:
                raise RuntimeError(
                    'PE file line %d contains non-numeric value: %s' %
                    (line_number, line)
                )

            data.append(row)
    finally:
        f.close()

    return data


def get_pe_data():
    """
    获取 PE 数据。

    本版本要求必须找到外部 PE_DATA_FILE。
    如果找不到文件，脚本会立即终止，不再使用任何内置备用数据，
    这样可以避免误用旧的 PE 数据。
    """
    file_path, searched_paths = resolve_file(PE_DATA_FILE)

    if not file_path:
        message = []
        message.append('PE data file was not found: %s' % PE_DATA_FILE)
        message.append('Script stopped before creating/modifying the model.')
        message.append('Please put the txt file in one of the following paths,')
        message.append('or change PE_DATA_FILE to the correct absolute path:')
        for path in searched_paths:
            message.append('  - %s' % path)
        raise RuntimeError('\n'.join(message))

    info('Use external PE data file: ' + file_path)
    data = read_pe_data_from_txt(file_path)

    if len(data) < LAYER_NUMBER:
        raise RuntimeError(
            'PE data rows are not enough. Required %d rows, but found %d rows.\n'
            'Please check the file: %s' %
            (LAYER_NUMBER, len(data), file_path)
        )

    if len(data) > LAYER_NUMBER:
        warn('PE data has %d rows, but only first %d rows will be used.' %
             (len(data), LAYER_NUMBER))

    return data[:LAYER_NUMBER]


def unique_geometry_items(items):
    """
    去除 Face/Edge/Vertex 等几何对象中的重复项。
    Abaqus 对象通常有 index 属性。
    """
    result = []
    used = {}

    for item in items:
        try:
            key = item.index
        except Exception:
            key = str(item)

        if key not in used:
            used[key] = True
            result.append(item)

    return tuple(result)


# =============================================================================
# 2. 建几何与分区
# =============================================================================

def create_part_geometry(model):
    """
    创建 38 × 9.5 × 1.29 的三维实体阿尔门试片，
    并在 z = PARTITION_Z 处建立水平分区面。
    """
    info('Creating part geometry...')

    delete_if_exists(model.parts, PART_NAME)

    sketch = model.ConstrainedSketch(name='__profile__', sheetSize=200.0)
    sketch.rectangle(point1=(0.0, 0.0), point2=(LENGTH_X, WIDTH_Y))

    part = model.Part(
        name=PART_NAME,
        dimensionality=THREE_D,
        type=DEFORMABLE_BODY
    )
    part.BaseSolidExtrude(sketch=sketch, depth=THICKNESS_Z)

    delete_if_exists(model.sketches, '__profile__')

    # 用 DatumPlane 进行分区，比 rpy 中的 Edge Parameter + PointNormal 更稳定。
    datum = part.DatumPlaneByPrincipalPlane(
        principalPlane=XYPLANE,
        offset=PARTITION_Z
    )
    part.PartitionCellByDatumPlane(
        datumPlane=part.datums[datum.id],
        cells=part.cells[:]
    )

    info('Geometry created. Partition plane z = %.6f.' % PARTITION_Z)
    return part


# =============================================================================
# 3. 网格划分
# =============================================================================

def append_edges_from_box(part, target_list, x, y, z_min, z_max, label):
    """
    用 BoundingBox 选择指定角点处的一段竖直边，并追加到 target_list。

    为什么不用 findAt：
        不同 Abaqus 版本中 findAt 在只给一个点时，可能返回 Edge，
        也可能返回 GeomSequence。如果把 GeomSequence 再放进 tuple 传给
        seedEdgeBySize，就会出现类似：
            TypeError: edges[0]; found GeomSequence, expected Edge/Vertex
        的错误。

    本函数直接用 getByBoundingBox 返回 EdgeArray，再逐个取出 Edge 对象，
    可以避免“嵌套 GeomSequence”导致的报错。
    """
    edges = part.edges.getByBoundingBox(
        xMin=x - TOL,
        xMax=x + TOL,
        yMin=y - TOL,
        yMax=y + TOL,
        zMin=z_min - TOL,
        zMax=z_max + TOL
    )

    if len(edges) == 0:
        raise RuntimeError(
            'No vertical edge was found for %s at x=%.6f, y=%.6f, z=%.6f~%.6f.\n'
            'Please check LENGTH_X, WIDTH_Y, THICKNESS_Z, PARTITION_Z and TOL.' %
            (label, x, y, z_min, z_max)
        )

    for edge in edges:
        target_list.append(edge)


def seed_vertical_edges(part):
    """
    对厚度方向的竖直边进行局部种子设置：
        - 上部喷丸影响层：0.01 mm；
        - 下部基体区域：0.10 mm。

    改进点：
        原版本通过 findAt 找边，某些 Abaqus 版本会把 EdgeArray/GeomSequence
        嵌套传入 seedEdgeBySize，从而报：
            TypeError: edges[0]; found GeomSequence, expected Vertex/Edge
        这里改为 BoundingBox 选边，并确保传入的是普通 Edge 对象 tuple。
    """
    info('Seeding vertical edges...')

    x_values = (0.0, LENGTH_X)
    y_values = (0.0, WIDTH_Y)

    top_edges = []
    base_edges = []

    for x in x_values:
        for y in y_values:
            # 上部喷丸影响层竖直边：z = PARTITION_Z ~ THICKNESS_Z
            append_edges_from_box(
                part=part,
                target_list=top_edges,
                x=x,
                y=y,
                z_min=PARTITION_Z,
                z_max=THICKNESS_Z,
                label='top affected layer edge'
            )

            # 下部基体竖直边：z = 0 ~ PARTITION_Z
            append_edges_from_box(
                part=part,
                target_list=base_edges,
                x=x,
                y=y,
                z_min=0.0,
                z_max=PARTITION_Z,
                label='base layer edge'
            )

    top_edges = unique_geometry_items(top_edges)
    base_edges = unique_geometry_items(base_edges)

    if len(top_edges) == 0:
        raise RuntimeError('No top vertical edges were collected. Mesh seeding stopped.')
    if len(base_edges) == 0:
        raise RuntimeError('No base vertical edges were collected. Mesh seeding stopped.')

    part.seedEdgeBySize(
        edges=top_edges,
        size=TOP_LAYER_SIZE_Z,
        deviationFactor=0.1,
        minSizeFactor=0.1,
        constraint=FINER
    )

    part.seedEdgeBySize(
        edges=base_edges,
        size=BASE_LAYER_SIZE_Z,
        deviationFactor=0.1,
        minSizeFactor=0.1,
        constraint=FINER
    )

    info('Top vertical edges seeded: %d, size = %.6f.' %
         (len(top_edges), TOP_LAYER_SIZE_Z))
    info('Base vertical edges seeded: %d, size = %.6f.' %
         (len(base_edges), BASE_LAYER_SIZE_Z))

def mesh_part(part):
    """
    对零件划分结构化六面体网格，并设置 C3D8R 单元类型。
    """
    info('Meshing part...')

    part.seedPart(
        size=GLOBAL_MESH_SIZE,
        deviationFactor=0.1,
        minSizeFactor=0.1
    )

    seed_vertical_edges(part)

    cells = part.cells[:]

    # 对简单长方体分区，STRUCTURED + HEX 通常可以稳定生成六面体网格。
    try:
        part.setMeshControls(regions=cells, elemShape=HEX, technique=STRUCTURED)
    except Exception:
        warn('Structured mesh control failed. Abaqus will use default mesh control.')

    elem_type_1 = mesh.ElemType(
        elemCode=C3D8R,
        elemLibrary=STANDARD,
        kinematicSplit=AVERAGE_STRAIN,
        secondOrderAccuracy=OFF,
        hourglassControl=DEFAULT,
        distortionControl=DEFAULT
    )
    elem_type_2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)
    elem_type_3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)

    part.setElementType(
        regions=(cells, ),
        elemTypes=(elem_type_1, elem_type_2, elem_type_3)
    )

    part.generateMesh()

    info('Mesh generated. Element number = %d.' % len(part.elements))


# =============================================================================
# 4. 创建集合
# =============================================================================

def create_element_set_by_box(part, set_name, z_min, z_max):
    """
    根据 BoundingBox 创建单元集合。
    这里默认选取整个试片 X/Y 范围内、指定 Z 厚度范围内的单元。
    """
    elements = part.elements.getByBoundingBox(
        xMin=0.0 - TOL,
        xMax=LENGTH_X + TOL,
        yMin=0.0 - TOL,
        yMax=WIDTH_Y + TOL,
        zMin=z_min - TOL,
        zMax=z_max + TOL
    )

    if len(elements) == 0:
        raise RuntimeError(
            'No elements selected for %s. z range = %.6f ~ %.6f. '
            'Please check mesh size and geometry coordinates.' %
            (set_name, z_min, z_max)
        )

    delete_if_exists(part.sets, set_name)
    part.Set(elements=elements, name=set_name)

    info('Created %-10s | elements = %-6d | z = %.6f ~ %.6f' %
         (set_name, len(elements), z_min, z_max))


def create_layer_sets(part):
    """
    创建：
        set-1 ~ set-30：上表面向下 0.30 mm 的喷丸影响层；
        Set-0          ：下部未喷丸影响区域。
    """
    info('Creating element sets...')

    # 下部基体区域 Set-0
    create_element_set_by_box(
        part=part,
        set_name=BASE_SET_NAME,
        z_min=0.0,
        z_max=PARTITION_Z
    )

    # 表层 30 个分层集合
    for i in range(LAYER_NUMBER):
        layer_id = i + 1
        z_max = THICKNESS_Z - i * TOP_LAYER_SIZE_Z
        z_min = THICKNESS_Z - (i + 1) * TOP_LAYER_SIZE_Z

        set_name = LAYER_SET_PREFIX + str(layer_id)

        create_element_set_by_box(
            part=part,
            set_name=set_name,
            z_min=z_min,
            z_max=z_max
        )


# =============================================================================
# 5. 材料、截面、材料方向
# =============================================================================

def create_base_material_and_section(model, part):
    """
    为 Set-0 创建普通弹性材料和截面。
    """
    info('Creating base material and section...')

    delete_if_exists(model.sections, BASE_SECTION_NAME)
    delete_if_exists(model.materials, BASE_MATERIAL_NAME)

    model.Material(name=BASE_MATERIAL_NAME)
    model.materials[BASE_MATERIAL_NAME].Elastic(
        table=((ELASTIC_MODULUS, POISSON_RATIO), )
    )

    model.HomogeneousSolidSection(
        name=BASE_SECTION_NAME,
        material=BASE_MATERIAL_NAME,
        thickness=None
    )

    region = part.sets[BASE_SET_NAME]
    part.SectionAssignment(
        region=region,
        sectionName=BASE_SECTION_NAME,
        offset=0.0,
        offsetType=MIDDLE_SURFACE,
        offsetField='',
        thicknessAssignment=FROM_SECTION
    )

    info('Assigned %s to %s.' % (BASE_SECTION_NAME, BASE_SET_NAME))


def create_layer_material_and_section(model, part, layer_id, expansion_row):
    """
    为某一表层集合创建材料、截面，并赋值给对应集合。
    """
    set_name = LAYER_SET_PREFIX + str(layer_id)
    material_name = LAYER_MATERIAL_PREFIX + str(layer_id)
    section_name = LAYER_SECTION_PREFIX + str(layer_id)

    delete_if_exists(model.sections, section_name)
    delete_if_exists(model.materials, material_name)

    model.Material(name=material_name)
    material = model.materials[material_name]

    material.Elastic(table=((ELASTIC_MODULUS, POISSON_RATIO), ))
    material.Expansion(type=ORTHOTROPIC, table=(expansion_row, ))

    model.HomogeneousSolidSection(
        name=section_name,
        material=material_name,
        thickness=None
    )

    region = part.sets[set_name]
    part.SectionAssignment(
        region=region,
        sectionName=section_name,
        offset=0.0,
        offsetType=MIDDLE_SURFACE,
        offsetField='',
        thicknessAssignment=FROM_SECTION
    )

    info('Assigned %-10s -> %-15s | alpha = %s' %
         (set_name, section_name, str(expansion_row)))


def create_layer_materials_and_sections(model, part, pe_data):
    """
    为 set-1 ~ set-30 分别创建材料和截面。
    """
    info('Creating layer materials and sections...')

    for i in range(LAYER_NUMBER):
        create_layer_material_and_section(
            model=model,
            part=part,
            layer_id=i + 1,
            expansion_row=pe_data[i]
        )


def assign_material_orientation(part):
    """
    为所有 cell 指定全局材料方向。
    对 ORTHOTROPIC 热膨胀材料，材料方向会影响 1、2、3 三个方向对应到全局坐标的方式。
    """
    info('Assigning material orientation...')

    region = regionToolset.Region(cells=part.cells[:])
    part.MaterialOrientation(
        region=region,
        orientationType=GLOBAL,
        axis=AXIS_3,
        additionalRotationType=ROTATION_NONE,
        localCsys=None,
        fieldName='',
        stackDirection=STACK_3
    )

    info('Material orientation assigned: GLOBAL, AXIS_3.')


# =============================================================================
# 6. 装配、分析步、边界条件、温度场
# =============================================================================

def create_assembly(model, part):
    """
    创建装配和实例。
    """
    info('Creating assembly instance...')

    assembly = model.rootAssembly
    assembly.DatumCsysByDefault(CARTESIAN)

    delete_if_exists(assembly.instances, INSTANCE_NAME)
    instance = assembly.Instance(
        name=INSTANCE_NAME,
        part=part,
        dependent=ON
    )

    assembly.regenerate()
    info('Assembly instance created: ' + INSTANCE_NAME)

    return assembly, instance


def create_static_step(model):
    """
    创建静力分析步 Step-1。
    """
    info('Creating static step...')

    # 删除已有同名分析步。
    if STEP_NAME in model.steps.keys():
        del model.steps[STEP_NAME]

    model.StaticStep(name=STEP_NAME, previous='Initial')
    info('Static step created: ' + STEP_NAME)


def get_side_faces_by_box(instance, side_name):
    """
    用 BoundingBox 直接选取装配实例上的对称面。

    为什么不用 findAt + tuple：
        某些 Abaqus 版本中，Region(faces=...) 不接受普通 Python tuple，
        而要求 Abaqus 自己返回的 GeomSequence/FaceArray。
        因此这里直接使用 instance.faces.getByBoundingBox(...)，
        其返回值可以直接传给 regionToolset.Region(faces=...)。

    side_name 可选：
        X_MIN：x = 0 面
        X_MAX：x = LENGTH_X 面
        Y_MIN：y = 0 面
        Y_MAX：y = WIDTH_Y 面
    """
    if side_name == 'X_MIN':
        faces = instance.faces.getByBoundingBox(
            xMin=0.0 - TOL,
            xMax=0.0 + TOL,
            yMin=0.0 - TOL,
            yMax=WIDTH_Y + TOL,
            zMin=0.0 - TOL,
            zMax=THICKNESS_Z + TOL
        )
    elif side_name == 'X_MAX':
        faces = instance.faces.getByBoundingBox(
            xMin=LENGTH_X - TOL,
            xMax=LENGTH_X + TOL,
            yMin=0.0 - TOL,
            yMax=WIDTH_Y + TOL,
            zMin=0.0 - TOL,
            zMax=THICKNESS_Z + TOL
        )
    elif side_name == 'Y_MIN':
        faces = instance.faces.getByBoundingBox(
            xMin=0.0 - TOL,
            xMax=LENGTH_X + TOL,
            yMin=0.0 - TOL,
            yMax=0.0 + TOL,
            zMin=0.0 - TOL,
            zMax=THICKNESS_Z + TOL
        )
    elif side_name == 'Y_MAX':
        faces = instance.faces.getByBoundingBox(
            xMin=0.0 - TOL,
            xMax=LENGTH_X + TOL,
            yMin=WIDTH_Y - TOL,
            yMax=WIDTH_Y + TOL,
            zMin=0.0 - TOL,
            zMax=THICKNESS_Z + TOL
        )
    else:
        raise RuntimeError(
            'Unknown side name: %s. Use X_MIN, X_MAX, Y_MIN or Y_MAX.' %
            side_name
        )

    if len(faces) == 0:
        raise RuntimeError(
            'No faces were found for %s. Please check LENGTH_X, WIDTH_Y, '
            'THICKNESS_Z, instance position and TOL.' % side_name
        )

    info('Selected %d face(s) for %s.' % (len(faces), side_name))
    return faces


def get_vertex_by_box(instance, coord):
    """
    用 BoundingBox 选取固定顶点，返回 VertexArray/GeomSequence。
    这样可避免 Region(vertices=(vertex,)) 在某些 Abaqus 版本中因 tuple 报错。
    """
    x, y, z = coord
    vertices = instance.vertices.getByBoundingBox(
        xMin=x - TOL,
        xMax=x + TOL,
        yMin=y - TOL,
        yMax=y + TOL,
        zMin=z - TOL,
        zMax=z + TOL
    )

    if len(vertices) == 0:
        raise RuntimeError(
            'No fixed vertex was found near %s. Please check FIX_VERTEX_COORD and TOL.' %
            (str(coord),)
        )
    if len(vertices) > 1:
        warn('More than one vertex was found near %s. All of them will be fixed.' %
             (str(coord),))

    return vertices


def create_boundary_conditions(model, instance):
    """
    创建边界条件：
        BC-1：YsymmBC，默认 y=0 面；
        BC-2：XsymmBC，默认 x=0 面；
        BC-3：EncastreBC，默认固定 (0,0,0) 顶点。

    注意：
        如果你的视频操作中选的是另一侧面，请修改用户参数区的 XSYMM_FACE、YSYMM_FACE。
    """
    info('Creating boundary conditions...')

    # 先删除已有同名边界条件，避免重复运行脚本冲突。
    for bc_name in ('BC-1', 'BC-2', 'BC-3'):
        delete_if_exists(model.boundaryConditions, bc_name)

    ysymm_faces = get_side_faces_by_box(instance, YSYMM_FACE)
    ysymm_region = regionToolset.Region(faces=ysymm_faces)
    model.YsymmBC(
        name='BC-1',
        createStepName='Initial',
        region=ysymm_region,
        localCsys=None
    )
    info('BC-1 YsymmBC applied on ' + YSYMM_FACE)

    xsymm_faces = get_side_faces_by_box(instance, XSYMM_FACE)
    xsymm_region = regionToolset.Region(faces=xsymm_faces)
    model.XsymmBC(
        name='BC-2',
        createStepName='Initial',
        region=xsymm_region,
        localCsys=None
    )
    info('BC-2 XsymmBC applied on ' + XSYMM_FACE)

    fixed_vertices = get_vertex_by_box(instance, FIX_VERTEX_COORD)
    fixed_region = regionToolset.Region(vertices=fixed_vertices)
    model.EncastreBC(
        name='BC-3',
        createStepName='Initial',
        region=fixed_region,
        localCsys=None
    )
    info('BC-3 EncastreBC applied on vertex %s' % (str(FIX_VERTEX_COORD),))


def create_temperature_field(model, instance):
    """
    创建温度预定义场：
        Initial：温度为 0；
        Step-1 ：温度为 1。

    由于各层材料的热膨胀系数就是 PE 数据，温度增量为 1 时，
    热应变增量 = alpha × ΔT = alpha。
    """
    info('Creating predefined temperature field...')

    field_name = 'Predefined Field-1'
    delete_if_exists(model.predefinedFields, field_name)

    region = regionToolset.Region(cells=instance.cells[:])

    model.Temperature(
        name=field_name,
        createStepName='Initial',
        region=region,
        distributionType=UNIFORM,
        crossSectionDistribution=CONSTANT_THROUGH_THICKNESS,
        magnitudes=(INITIAL_TEMPERATURE, )
    )

    model.predefinedFields[field_name].setValuesInStep(
        stepName=STEP_NAME,
        magnitudes=(STEP_TEMPERATURE, )
    )

    info('Temperature field created: Initial = %.6f, %s = %.6f.' %
         (INITIAL_TEMPERATURE, STEP_NAME, STEP_TEMPERATURE))


# =============================================================================
# 7. Job 与保存
# =============================================================================

def create_job():
    """
    创建 Job-1。
    默认只创建，不提交。
    """
    info('Creating job...')

    delete_if_exists(mdb.jobs, JOB_NAME)

    mdb.Job(
        name=JOB_NAME,
        model=MODEL_NAME,
        description='A-Almen equivalent peening simulation by thermal expansion method',
        type=ANALYSIS,
        atTime=None,
        waitMinutes=0,
        waitHours=0,
        queue=None,
        memory=90,
        memoryUnits=PERCENTAGE,
        getMemoryFromAnalysis=True,
        explicitPrecision=SINGLE,
        nodalOutputPrecision=SINGLE,
        echoPrint=OFF,
        modelPrint=OFF,
        contactPrint=OFF,
        historyPrint=OFF,
        userSubroutine='',
        scratch='',
        resultsFormat=ODB,
        numThreadsPerMpiProcess=0,
        numCpus=NUM_CPUS,
        numDomains=NUM_CPUS,
        numGPUs=0
    )

    info('Job created: %s, CPUs = %d.' % (JOB_NAME, NUM_CPUS))


def save_cae_file():
    """
    保存 CAE 文件。

    路径优先级：
        1. 如果 CAE_SAVE_PATH 非空，直接保存到这个完整路径；
        2. 否则使用 CAE_SAVE_DIRECTORY + CAE_SAVE_FILE_NAME；
        3. 如果 CAE_SAVE_DIRECTORY 为空，则保存到 Abaqus 当前工作目录。
    """
    if not SAVE_CAE:
        return

    if CAE_SAVE_PATH:
        save_path = CAE_SAVE_PATH
    else:
        if not CAE_SAVE_FILE_NAME.lower().endswith('.cae'):
            raise RuntimeError(
                'CAE_SAVE_FILE_NAME must end with .cae, current value: %s' %
                CAE_SAVE_FILE_NAME
            )

        if CAE_SAVE_DIRECTORY:
            save_dir = CAE_SAVE_DIRECTORY
        else:
            save_dir = os.getcwd()

        save_path = os.path.join(save_dir, CAE_SAVE_FILE_NAME)

    save_path = os.path.abspath(save_path)
    save_dir = os.path.dirname(save_path)
    safe_make_dir(save_dir)

    mdb.saveAs(pathName=save_path)
    info('CAE saved to: ' + save_path)


def submit_job():
    """
    提交计算。
    默认主程序中不调用，避免脚本运行后直接开始计算。
    如需自动提交，可在 main() 末尾取消注释。
    """
    mdb.jobs[JOB_NAME].submit(consistencyChecking=OFF)
    mdb.jobs[JOB_NAME].waitForCompletion()
    info('Job finished: ' + JOB_NAME)


# =============================================================================
# 8. 主程序
# =============================================================================

def main():
    """
    主流程：
        建模 -> 分区 -> 网格 -> 集合 -> 材料截面 -> 装配 -> 分析步 -> 约束 -> 温度场 -> Job -> 保存
    """
    # 记录几何时采用坐标方式，通常比 INDEX 更稳定。
    session.journalOptions.setValues(
        replayGeometry=COORDINATE,
        recoverGeometry=COORDINATE
    )

    info('=' * 72)
    info('Start A-Almen full workflow script.')

    # 基本参数检查
    if abs(AFFECTED_DEPTH - LAYER_NUMBER * TOP_LAYER_SIZE_Z) > 1.0e-8:
        warn(
            'AFFECTED_DEPTH is not equal to LAYER_NUMBER * TOP_LAYER_SIZE_Z. '
            'AFFECTED_DEPTH = %.6f, layers total = %.6f.' %
            (AFFECTED_DEPTH, LAYER_NUMBER * TOP_LAYER_SIZE_Z)
        )

    if PARTITION_Z <= 0.0 or PARTITION_Z >= THICKNESS_Z:
        raise RuntimeError(
            'Invalid PARTITION_Z = %.6f. It must be between 0 and THICKNESS_Z.' %
            PARTITION_Z
        )

    pe_data = get_pe_data()

    model = reset_model(MODEL_NAME)

    part = create_part_geometry(model)
    mesh_part(part)
    create_layer_sets(part)

    create_base_material_and_section(model, part)
    create_layer_materials_and_sections(model, part, pe_data)
    assign_material_orientation(part)

    assembly, instance = create_assembly(model, part)

    create_static_step(model)
    create_boundary_conditions(model, instance)
    create_temperature_field(model, instance)

    create_job()
    save_cae_file()

    # 如需脚本运行后自动提交计算，请取消下面两行注释。
    # submit_job()

    info('Full workflow finished successfully.')
    info('=' * 72)


if __name__ == '__main__':
    main()
