# 蛋白／酶设计上下文入口

这里整合了此前为后续酶设计收集的多工作流知识，供 Agent_Gaussian_TS 的 Codex 决策层读取。来源快照截至 2026-10-17，整合日期为 2026-10-19；共 17 篇唯一研究、13 类综合路线，重点为 RFD。本目录是文献与设计判断储备，不是完整文献审计或独立计算复现。后续开发的可执行能力和边界另见 [酶设计工具接口](../ENZYME_DESIGN.md)。

## 阅读顺序

1. [项目衔接与判断规则](AGENT_INTEGRATION.md)：先确定 Gaussian 在当前任务中的作用、输入和交接条件。
2. [酶设计路线与 RFD 实用指南](enzyme_design_guide_zh.md)：路线选择、实验判据与失败诊断。
3. [研究来源目录](source_catalog_zh.md)与[详细登记](source_register.json)：核对研究、版本、数值分母及未决问题。
4. [路线比较](route_selection.json)：已有亲本、天然家族、从头生成、匹配及其他路线的入口条件。
5. 按所选路线读取 [RFD 丝氨酸水解工作流](rfd_serine_workflow.json)、[RFD3 工作流](rfd3_workflow.json)、[RFD 代码来源](rfd_code_provenance.json)与[补充材料核查](rfd3_supplement_audit.json)。
6. 需要规划从量化到蛋白设计的交接时，读取[历史 V1 工程方案](enzyme_design_workflow_v1.md)。

[原知识库索引](knowledge_index.json)、[后续核查队列](rfd_focus.json)、[历史交付审计](delivery_audit.json)与[研究记录](research_history.jsonl)完整保留。它们描述原收集任务的进度，不是当前计算任务的完成标签；历史待办与后来交付不一致时以具体文件及证据核查范围为准。

## 使用边界

- 用户当前研究条件、仓库现有工具 API 和独立验证器约束优先。本文献集不改变 Gaussian CLI、验证标签或作业预算。
- 历史 V1 中“完整游离氨基酸、不显式加水、不截断/封端”等约束属于当时项目范围；不得自动覆盖新任务明确采用的侧链簇、金属、水或截断模型。也不能因文献中出现某组分就擅自加入当前计算。
- V1 的服务器配置、软件可用性、样本数和并发量是历史记录或示例；执行时重新核实并采用本次预算。论文阈值不是通用验收标准。
- [导入清单](import_manifest.json)保存所有导入文本的 SHA-256；文本按原样保存，仅 V1 规范化了文末空行，其原始与导入哈希均已登记，本次没有重新核查外部论文。`papers/` 原文、补充材料及上游代码缓存不随 Git 提交，清单中保存其名称、大小和 SHA-256；实际使用前按登记来源取得并核查。已有“部分核查”不能扩大解释成全文或原始数据均已审计。

## 读取方式

在仓库根目录直接读取本目录的文档及按需 JSON。`python3 -m gaussian_ts context --run-dir ...` 仍只输出该运行的状态和工具反馈，不会自动注入本文献库；可选 HTTPPolicy 也不自动加载这些文件。主入口 AGENTS.md 与 CODEX_GUIDE.md 提醒 Codex 在相关研究任务开始时主动读取，并将采用的来源及假设记录到 `problem.json` 和每次动作的 `reason` 中。
