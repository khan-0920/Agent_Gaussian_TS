# Codex 使用本工具包

本仓库是给 Codex/GPT 的 Gaussian 过渡态搜索辅助工具包。Codex 负责定义研究问题、解释计算反馈、选择动作；Python 模块提供可审计的执行工具与独立验证。无需启动独立 LLM 服务或 `Agent.run` 才能使用。

开始新的计算化学任务时，读 `docs/CODEX_GUIDE.md`；工具参数见 `docs/API.md`。用户的研究问题、指定条件与预算优先，假设/未知信息要与已知输入分开记录。

常用入口（从仓库根目录运行）：

```bash
python3 -m gaussian_ts define --problem examples/problem.json --config examples/hydrogen_exchange/reaction.json
python3 -m gaussian_ts tools
python3 -m gaussian_ts init --config REACTION_JSON --run-dir runs/NEW_RUN
python3 -m gaussian_ts context --run-dir runs/NEW_RUN
python3 -m gaussian_ts apply --run-dir runs/NEW_RUN --action ACTION_JSON --actor codex
python3 -m gaussian_ts report --run-dir runs/NEW_RUN
```

读取 context，形成一个带依据的动作，写入 JSON，再调用 apply；执行后重新读取 context。用户条件不全时使用 define 列出缺项，不根据示例猜电荷、自旋、映射或计算水平。无需把动作顺序固定为某一流水线。

工具参数中的原子是稳定 ID；Gaussian 行号只在生成输入时转换。几何动作显式指定移动片段。计算输出的 `ok` 表示工具执行成功，实际结果见 job.status、反馈和独立验证标签。

只有原始真实证据经验证器得到 `VALIDATED_TS` 时才能称为已验证目标过渡态。模拟日志、正常结束、单个虚频都不足以单独证明成功；不替 LLM 的判断制造频率或 IRC 证据。

工程开发验证：`python3 -m unittest discover -v`。核心无第三方运行依赖。真实 Gaussian/许可/主机环境由用户提供；没有 g16/g09 时可运行 `python3 -m examples.offline_demo` 检查工程链路，科学验收保持未完成。
