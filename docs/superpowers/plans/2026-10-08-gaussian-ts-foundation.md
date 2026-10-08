# Gaussian TS Foundation Implementation Plan

> **For agentic workers:** Use superpowers:executing-plans in this session; write failing behavior tests before implementation.

**Goal:** 搭建供 Codex/GPT 定义 TS 研究问题、调用 Gaussian/几何工具与验证证据的辅助工具包。

**Architecture:** Codex/GPT 自身负责研究决策，类型化动作经确定性工具执行，Gaussian 原始日志转换为证据。问题定义/context/schema 辅助模型，持久化状态、预算和独立验证器约束计算。

**Tech Stack:** Python >=3.9，标准库，unittest，setuptools。

**Spec:** `docs/superpowers/specs/2026-10-08-gaussian-ts-foundation-design.md`

## Global Constraints

- 原子使用稳定、正整数映射；Gaussian 输入编号按当前坐标顺序转换。
- 不执行 LLM 生成的 shell 或原始 route。
- 默认保存所有结构、日志和动作；无真实证据不得声称成功。
- 无 Gaussian 环境通过离线用例验收工程部分。

## Review Focus

- 不连续原子映射和错误自旋：模型与输入生成测试拒绝不一致数据。
- Gaussian Link1、多频率块、截断日志：解析测试保留完整性信息，验证器拒绝混用证据。
- symlink 和目录穿越：存储与作业执行测试阻止越界。
- 进程超时或中断：执行与状态测试终止自己的进程、保守保留预算和动作。
- 相同频率但错误机理：验证测试核对模式方向与 IRC/端点链路。

### Task 1: 结构和几何

Files: `gaussian_ts/models.py`, `gaussian_ts/geometry.py`, `tests/test_geometry.py`。
Interfaces: `Structure.from_xyz`, `Reaction`, `edit_geometry`, `inspect_clashes`。

- [x] 写 XYZ、映射、电荷自旋、键长/角度/二面角、刚体片段与碰撞行为测试。
- [x] 运行测试确认缺失功能，再实现并通过测试。

### Task 2: Gaussian 计算边界

Files: `gaussian_ts/gaussian.py`, `gaussian_ts/parser.py`, `gaussian_ts/runner.py`, `tests/test_gaussian.py`, `tests/fixtures/*`。
Interfaces: `JobSpec`, `prepare_input`, `parse_log`, `LocalRunner.start/wait/cancel`。

- [x] 写合法输入、参数拒绝、频率模式/收敛/日志截断、进程超时测试并确认失败。
- [x] 实现并验证；每条日志保留来源摘要，不混合 Link1 作业。

### Task 3: 独立验证器

Files: `gaussian_ts/validation.py`, `tests/test_validation.py`。
Interfaces: `validate_ts(reaction, evidence)` 返回五种手册标签和分项依据。

- [x] 写正常结束不等于成功、数值小虚频、错误模式、缺 IRC、端点错误、完整证据测试并确认失败。
- [x] 实现结构/计算水平和证据父子链检查；合成结果不能通过。

### Task 4: 持久化工具引擎

Files: `gaussian_ts/store.py`, `gaussian_ts/engine.py`, `tests/test_engine.py`。
Interfaces: `Session.create/open/execute/observe`，`Action` 字典。

- [x] 写审计、撤销、恢复、去重、预算、路径边界和错误动作测试并确认失败。
- [x] 实现事务式 checkpoint、单写者锁、记录失败与中断；运行全套测试。

### Task 5: 入口、LLM 和交付

Files: `gaussian_ts/agent.py`, `gaussian_ts/cli.py`, `README.md`, `docs/API.md`, `examples/*`, `pyproject.toml`。
Interfaces: `DecisionPolicy.decide(observation)`，`python -m gaussian_ts`。

- [x] 写反馈依赖策略、多轮/暂停恢复、CLI 和 HTTP 错误边界测试并确认失败。
- [x] 实现、运行全部测试及命令行离线演示；记录真实计算未验收的限制。

### Task 6: Codex 工具包入口（用户澄清定位）

Files: `gaussian_ts/toolkit.py`, `AGENTS.md`, `docs/CODEX_GUIDE.md`, `examples/problem.json`, `tests/test_toolkit.py`。
Interfaces: `define_problem`, `tool_catalog`, `function_tools`, `export_context`；CLI `define/context/tools`。

- [x] 写问题缺项、固定条件/假设分离、上下文反馈及工具参数 schema 测试并确认缺失功能。
- [x] 实现 Codex 主入口与仓库说明；独立 Agent/HTTP 保留为可选宿主示例。
- [x] 修复独立审查发现的频率自由度/坐标系、复合截断日志、审计恢复窗口与 QST2 真实输入去重；补入回归。
