# 开发记录

2026-10-08：目录只有需求手册，无 Git 仓库、Gaussian 或已有实现。直接在用户提供的目录开发，不创建工作树、不自动提交。使用标准库支持当前 Python 3.9.6，避免以依赖安装阻塞基础验证。

Ruling: 用户明确要求开始开发，需求手册作为设计依据；采用当前会话直接执行并落盘设计/计划，不重复请求启动授权。

Pre-flight: models → geometry/input/parser → validation → engine → agent/CLI；所有接口遵循设计文件，真实计算与离线演示证据分离。

用户随后要求当前目录初始化 Git，已建立 main 分支；Git 忽略运行结果、二进制 Gaussian 文件、本地配置和密钥。

用户澄清产品是 Codex/GPT 辅助工具包：主入口调整为 define/context/tools/apply/report，增加 AGENTS.md 和 CODEX_GUIDE；独立 Agent/HTTP 仅保留作可选宿主示例。需求手册原文保留。

Task 1–6: 已实现并通过对应回归。独立代码审查发现的 5 个问题均复现后修复；追加恢复缺 journal 时原动作审计、丢失日志保守判定的回归。

验收分层：标准库工具和离线闭环验收；真实 Gaussian、真实 LLM 化学决策及已知反应的科学成功率尚未验收。当前主机 doctor 未检测到 Gaussian。
