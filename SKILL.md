---
name: shot-peening
description: Abaqus 2025 喷丸强化仿真全流程：单丸粒标定→随机多丸粒→固有应变提取→Almen 弧高，含参数化、命名规范与覆盖率公式
---

# 喷丸强化仿真全流程（Abaqus 2025）：单丸粒 → 随机多丸粒 → Almen 试片

在 Abaqus/CAE 2025 中用 Python 脚本完成喷丸强化仿真：单丸粒冲击标定 → 随机多丸粒喷丸 → 固有应变提取 → Almen 试片等效喷丸 → 弧高测量。对应研究报告《空心叶片喷丸仿真研究报告》第三章。

## 0. 单位制、目录与执行总则

- 单位制：长度 mm，时间 s，质量 tonne，力 N，应力 MPa，密度 tonne/mm³，速度 mm/s。
- **脚本副本已随 skill 存放（自包含）**：`C:\Users\win11\AppData\Roaming\reasonix\skills\shot-peening/scripts/` 下的 `one-shot/`、`random-shot/`、`A-Almen/`（含 8 个主流程 .py 与配套计算器 exe）。执行时用**绝对路径**调用这些副本，原脚本不动。
- **工作目录（结果输出目录）**：用户指定一个根文件夹，在其中新建 `one-shot/`、`random-shot/`、`A-Almen/` 三个子文件夹；**每部分工作目录 = 各自子文件夹**。执行前 `Set-Location` 到该文件夹（脚本 `USE_CURRENT_WORK_DIRECTORY=True` → 产物输出到工作目录，与脚本所在位置分离）。
- 每个速度/气压跑**一条完整链路**；多速度的 model/job 共存于**同一 CAE 文件**（每部分一个 CAE），靠命名区分。
- 所有脚本**默认只创建 Job 不提交**。执行者**不改原脚本**，复制出调整版副本（`*_submit.py`）修改后运行：
  - 建模/静态脚本：`SUBMIT_JOB = False` → `True`；`OPEN_SOURCE_CAE` 为 False 的改 `True`（noGUI 空会话需自动打开源 CAE）；`RESET_DATABASE`/`DELETE_EXISTING_MODEL`/`RECREATE_MODEL` 保持 `False`（多 model 共存）。
  - `A_Almen_rebuild_abaqus2025_V2.py`：取消 `# submit_job()` 注释启用提交。
  - 弧高脚本：**注释 `write_points_csv(...)` 调用**（默认不生成 CSV）。
- 运行：建模/静态脚本 `abaqus cae noGUI=<脚本绝对路径>`；纯后处理（弹坑、PE、弧高）`abaqus python <脚本绝对路径>`。提交后 `waitForCompletion` 或轮询 `.lck` 消失。
  - 调用示例（以单丸粒建模为例）：`Set-Location <工作目录>\one-shot; abaqus cae noGUI="C:\Users\win11\AppData\Roaming\reasonix\skills\shot-peening\scripts\one-shot\one_shot_rebuild_abaqus2025_v7_submit.py"`——脚本副本在 skill 内，产物（ODB/CAE/报告）输出到工作目录。
- **每个任务/计算完成后保存 CAE**（`mdb.save()`/`saveAs` 保存回本部分 CAE）。
- Abaqus 2025 要点：打开已有 CAE 用 `openMdb(pathName=...)`（不是 `mdb.open()`）；stdout 中文在 PowerShell 可能乱码，不影响结果文件。

## 1. 命名规范（贯穿全流程，用户已确认）

| 变量类型 | 标识格式 | 示例 |
|---|---|---|
| 弹丸速度 | `39_1ms`（速度去小数点下划线 + ms） | 39.1 m/s → `39_1ms` |
| 喷丸气压 | `05MPa` | 0.5 MPa → `05MPa` |

- 三部分 model/job 名**统一不加前缀**（分属不同 CAE，不冲突）：`Model-<标识>` / `Job-<标识>`。
- 静态：`Model-<标识>-Static` / `Job-<标识>-Static`。
- CAE 文件每部分固定：`one_shot.cae`、`Random-Shot.cae`、`A-Almen.cae`。
- 弹坑摘要：`crater-<标识>_measurement_summary.txt`；PE 文件：`Job-<标识>-Static-PE.txt`；弧高报告：`A_Almen-<标识>.txt`。
- ODB 由作业名决定。**不用逗号、点号**（下划线允许）。

## 2. 弹丸数计算（覆盖率公式，与“弹丸数目计算器”一致）

```
C = 100 × (1 − exp(−N·π·a²/S))   →   N = ceil(−S·ln(1−C/100) / (π·a²))
```
- `a`：单丸粒弹坑 X/Y 平均半径（mm），取自弹坑测量摘要。
- `S`：中心加密区面积 = (2×FINE_ZONE_HALF_X)×(2×FINE_ZONE_HALF_Y)（默认 1×1 = 1 mm²）。
- `C`：目标覆盖率（用户给定，如 50、98）。结果**向上取整**。
- 验证案例：N=85、a=0.05098、S=1 → C=50%（与 50fugai 案例吻合）。
- 材料速查（报告表 3.1~3.5）：S110 弹丸 7.85e-9 / 2.1e5 / 0.3；TC4 4.5e-9 / 1.2e5 / 0.34，JC A=961 B=902 C=0.01 n=0.87；SAE1070 7.85e-9 / 2.1e5 / 0.3，JC A=961.06 B=901.62 C=0.02 n=0.87。
- 材料参数**只在一开始（第一个速度）设置一次，同一 CAE 内后续速度继承**；用户提出更改才改。

## 3. 逐脚本参数规则

### 阶段 A · 单丸粒（`one-shot/`）

**A1 `one_shot_rebuild_abaqus2025_v7.py`（显式建模）**
- 命名：`CAE_FILE_NAME='one_shot.cae'`（固定）、`MODEL_NAME='Model-<标识>'`、`JOB_NAME='Job-<标识>'`、`RESET_DATABASE=False`。
- 弹丸：`SHOT_VELOCITY_Z`（mm/s 负值，每次速度改）；`SHOT_RADIUS=0.15`、`SHOT_MESH_SIZE=0.01` 固定；`SHOT_TARGET_GAP=0.02`（**副本已按报告改好**）；材料 `SHOT_DENSITY=7.85e-9`、`SHOT_YOUNG_MODULUS=210000`、`SHOT_POISSON_RATIO=0.30`（首次设置，后续继承）。
- 靶材：2×2×1mm、网格 0.01/0.10mm、细化层 0.3mm 固定；`TARGET_DENSITY=4.5e-9`、`TARGET_YOUNG_MODULUS=120000`、`TARGET_POISSON_RATIO=0.34`；JC **副本已按报告改好** `JC_A=961 JC_B=902 JC_C=0.01 JC_N=0.87 JC_M=0`（首次设置，后续继承）。
- 接触：`FRICTION_COEFFICIENT=0.2`（**副本已按报告改好**）。
- 作业：`NUM_CPUS`/`NUM_DOMAINS` 按机器最大核心留余量（如 20 核→16/18）；`EXPLICIT_PRECISION=DOUBLE_PLUS_PACK`、`NODAL_OUTPUT_PRECISION=FULL`；分析步时间 AUTO。
- 副本：`SUBMIT_JOB=True`。

**A2 `one_shot_static_rebuild_abaqus2025_v3.py`（静态再平衡）**
- 命名衔接：`SOURCE_ODB_FILE_NAME='Job-<标识>.odb'`、`SOURCE_MODEL_NAME='Model-<标识>'`、`STATIC_MODEL_NAME='Model-<标识>-Static'`、`STATIC_JOB_NAME='Job-<标识>-Static'`、`INITIAL_STATE_NAME='Initial-State-from-Job-<标识>'`。
- 实例名 `shot-1`/`target-1` 固定（与 A1 对应）；`SAVE_AS_NEW_CAE=False`（保存回 `one_shot.cae`）。
- 固定：`OPEN_SOURCE_CAE=True`、`UPDATE_REFERENCE_CONFIGURATION=False`、`STATIC_NLGEOM=True`、`REPLACE_EXISTING_STATIC_*=True`、`NUM_CPUS` 按机器。
- 副本：`SUBMIT_JOB=True`。

**A3 `measure_crater_from_odb_abaqus2025_v3.py`（弹坑测量，abaqus python）**
- `ODB_FILE_NAME='Job-<标识>-Static.odb'`、`SUMMARY_FILE_NAME='crater-<标识>_measurement_summary.txt'`。
- `WRITE_CSV_FILES=False`（**默认不生成 CSV**）；剖面/表面 CSV 文件名含标识。
- 固定：`IMPACT_CENTER_X/Y=0`、`FAR_FIELD_RADIUS=None`、`EDGE_METHOD='U3_PEAKS'`。
- 产物：摘要含 **X/Y 平均半径** → 供 §2 算弹丸数。

### 阶段 B · 随机多丸粒（`random-shot/`）

**B1 `random_shot_rebuild_abaqus2025_V9.py`（显式建模）**
- 命名：`MODEL_NAME='Model-<标识>'`、`JOB_NAME='Job-<标识>'`（不加前缀）、`OUTPUT_CAE_PATH='Random-Shot.cae'`、`DELETE_EXISTING_MODEL=False`。
- `SHOT_NUMBER` ← §2 覆盖率公式计算值。
- `SHOT_VELOCITY_Z` 与单丸粒同速度；`RANDOM_SEED=2026`（可复现）。
- **弹丸材料副本已搬取单丸粒值**（原默认 2.3e-9/300000 已改为）：`SHOT_DENSITY=7.85e-9`、`SHOT_E=210000`（**变量名是 `SHOT_E`/`SHOT_NU`，不是 `SHOT_YOUNG_MODULUS`**）、`SHOT_NU=0.30`。
- 靶材 TC4 默认已与单丸粒一致（4.5e-9/1.2e5/0.34 + JC 961/902/0.01/0.87），无需改；几何/网格/摩擦 0.2 一致。
- 固定：`MIN_CENTER_DISTANCE_FACTOR=0.5`、`CREATE_LAYER_ELEMENT_SETS`（30 层 SET-LAYER-01~30）、`EXPORT_SHOT_COORDINATES=True`、`CONTACT_SCOPE`、`FIELD_OUTPUT_INTERVALS=200`、`NUM_CPUS` 按机器。
- 副本：`SUBMIT_JOB=True`。

**B2 `random_shot_static_rebuild_abaqus2025_V2.py`（静态再平衡）**
- 命名衔接：`SOURCE_ODB_FILE_NAME='Job-<标识>.odb'`、`SOURCE_MODEL_NAME='Model-<标识>'`、`STATIC_MODEL_NAME='Model-<标识>-Static'`、`STATIC_JOB_NAME='Job-<标识>-Static'`、`INITIAL_STATE_NAME='Initial-State-from-Random-Shot'`（固定名）。
- 实例识别自动：`SHOT_INSTANCE_PREFIXES=('shot-',...)` 匹配 `shot-001` 等、`TARGET_INSTANCE_NAME='target-1'`、`AUTO_DETECT_*=True`。
- `SAVE_AS_NEW_CAE=False`（保存回 `Random-Shot.cae`）；`OPEN_SOURCE_CAE` 副本改 `True`。
- 固定：`UPDATE_REFERENCE_CONFIGURATION=False`、`STATIC_NLGEOM=True`、`REPLACE_EXISTING_STATIC_*=True`、`NUM_CPUS` 按机器。
- 副本：`SUBMIT_JOB=True`。

**B3 `extract_random_shot_PE_for_Almen.py`（PE 提取，abaqus python）**
- `SOURCE_ODB_FILE_NAME='Job-<标识>-Static.odb'`（静态后多喷丸模型）、`OUTPUT_TEXT_FILE_NAME=None` → 自动 `Job-<标识>-Static-PE.txt`（30 行无表头三列）。
- `GENERATE_CSV_FILE=False`、`GENERATE_SUMMARY_FILE=False`（**默认不生成 CSV/摘要**）。
- 固定：`FRAME_SELECTION='LAST'`、`PE_COMPONENTS=('PE11','PE22','PE33')`、`LAYER_SET_*`（SET-LAYER-01~30、0.01mm/层）、`OUTPUT_FLOAT_FORMAT='%.10f'`。

### 阶段 C · Almen 试片（`A-Almen/`）

**C1 `A_Almen_rebuild_abaqus2025_V2.py`（等效建模）**
- 命名：`MODEL_NAME='Model-<标识>'`、`JOB_NAME='Job-<标识>'`（不加前缀）、`PE_DATA_FILE='Job-<标识>-Static-PE.txt'`（**从 random-shot 复制到本文件夹**）、`CAE_SAVE_FILE_NAME='A-Almen.cae'`、`RECREATE_MODEL=False`。
- 固定：`PART_NAME='A-Almen'`、`INSTANCE_NAME='A-Almen-1'`、材料 `Material-SAE`、截面 `Section-SAE-*`、集合 `set-1~30`/`Set-0`；试片 38×9.5×1.29mm、表层 0.30mm/30 层（0.01mm/层）、SAE1070（205000/0.29）、`NUM_CPUS` 按机器。
- 副本：启用 `submit_job()`（默认被注释）。PE 行数必须 ≥ `LAYER_NUMBER=30`，否则脚本终止。

**C2 `A_Almen_arc_height_support_plane_corrected_v2.py`（弧高，abaqus python）**
- `ODB_FILE='Job-<标识>.odb'`（**与 C1 作业名对应**）、`SUMMARY_FILE_NAME='A_Almen-<标识>.txt'`。
- **默认不生成 CSV**：副本中注释 `write_points_csv(...)` 调用（脚本无开关）。
- 固定：`INSTANCE_NAME='A-Almen-1'`、`STEP_NAME='Step-1'`、`MODEL_TYPE='QUARTER_SYMMETRY'`、`MEASUREMENT_SURFACE='BOTTOM'`（下表面未喷丸侧）、支撑距 31.75×15.875、`APPLY_INITIAL_FRAME_CORRECTION=True`。

## 4. 完成检查清单

- 每步：`*.sta` 显示 `THE ANALYSIS HAS COMPLETED SUCCESSFULLY`；`.lck` 不存在；CAE 已保存。
- 弹坑摘要含 X/Y 平均半径；PE 文件恰好 30 行三列且随深度衰减到 0；弧高报告生成、四支撑点残差≈0。
- 多速度场景：同一 CAE 内各速度 model/job 命名可区分，无覆盖；材料参数仅首次设置、后续继承一致。
