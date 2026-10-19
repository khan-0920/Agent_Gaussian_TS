# Theozyme 到 de novo enzyme：实现与验收

目标是让当前 Codex/GPT 通过可审计工具，从明确的 theozyme 生成、设计并计算评估候选酶，交付序列、结构、分阶段指标、淘汰原因与实验候选清单。计算筛选不等于证明实验活性。保留既有 Gaussian TS 接口和独立验证标准。

## 完成条件

- 输入：PDB/mmCIF theozyme、催化残基和固定原子、配体化学表示与稳定原子映射、反应状态、几何约束、TS/假设证据范围、资源预算；拒绝歧义身份和缺失必要映射。
- 编排：持久化 campaign、阶段依赖、不可覆盖的尝试目录、超时/进程组清理、并发互斥、失败重试、断点恢复、原始日志、输入/输出哈希与版本记录；不得把旧输出当新成功。
- 生成：RFD3 1.0.0 显式入口与模式、完整输出/催化残基重编号映射、生成后几何和拓扑检查。
- 序列：LigandMPNN 固定催化残基，配体上下文和侧链打包，序列/结构身份核对；每个骨架与每条序列分别登记。
- 优化：Rosetta 受约束结构优化、配体参数及原子名契约、金属/催化位点几何复核；缺失化学参数不能静默丢弃配体。
- 回预测：RF3 从序列与明确配体状态重新预测，记录置信度、全局/催化位点几何恢复与配体身份；不能把输入设计结构改名充当预测结构。
- 活性位点评估：PLACER 构象集合、各状态几何通过比例、模型分数和未覆盖化学状态；缺失数据不算通过。
- 筛选交付：显式阈值、结构/序列多样性、候选 FASTA/CSV/结构/映射/来源及失败清单；区分计算候选与实验验证。
- 智能体接口：CLI、Python 和 function-tool schema，按需 context 与带 reason 的单步动作，支持从完整 campaign 自动执行到报告。
- 验证：有意义的失败/映射/恢复测试和既有 Gaussian 回归；在 guest 上用真实工具完成小规模完整流程。无命中可以是有效筛选结果，但必须有真实全流程评估证据。

## 版本基线（2026-10-19 核查）

旧 `/home/guest/miniconda3/envs/rfd3` 为 rc-foundry 0.2.0 / torch 2.7.1+cu128。新 `/data2/guest/rfd3_foundry1_20261018/env_validated` 为 rc-foundry 1.0.0 / torch 2.11.0 / atomworks 3.0.1 / biotite 1.6.0，安装源为服务器 `runtime_source` 快照。官方 production 的 pyproject 同样为 1.0.0，但未据此宣称源码逐字相同。

新环境的既有 `benchmark_pinned_fast` 和 `benchmark_pinned_legacy` 记录各 5 个成功作业、33 个输出且身份检查通过。这只支持既有运行，不替代新适配器验收。新流程使用新版绝对路径，显式记录 inference_mode 与 checkpoint。原环境、全局 PATH、其他在跑任务不作升级或替换。

官方来源：https://github.com/RosettaCommons/foundry/blob/production/pyproject.toml

## 实施顺序

先实现输入契约、结构/几何工具和可恢复编排；随后逐个接入真实生成、设计、优化、回预测和集合评估适配器；最后验证完整小规模运行、文档和交付。此计划的勾选或测试通过不替代完成条件中的原始证据。

## 本轮交付证据

七阶段工具、CLI/Python/action 接口、示例和实际支持边界见 [ENZYME_DESIGN.md](ENZYME_DESIGN.md)。guest `run02` 已完成 2 骨架 → 2 序列 → 2 优化 → 2 回预测 → 4 PLACER 样本 → 筛选；因催化几何未恢复，最终 0 个合格候选。逐阶段版本、原始结果哈希和失败重试见 [验证摘要](validation/enzyme_smoke_20261019.json)。不把运行成功解释为新酶活性，也不声称已覆盖所有历史路线或化学状态。
