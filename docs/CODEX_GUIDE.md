# 用 Codex 定义问题和调用工具

本工具包把 Codex/GPT 放在研究决策层。它提供三类辅助：把自然语言问题整理成显式研究条件；把历史计算压缩成结构化上下文；执行经过验证的几何/Gaussian 动作并保存证据。

## 定义问题

先写一个 `problem.json`，包含 question、fixed_conditions、hypotheses、uncertainties 和可选 notes，示例见 `examples/problem.json`。问题应明确目标基元步骤、哪些条件固定、哪些机理是待检验假设及什么证据才算成功。

```bash
python3 -m gaussian_ts define --problem examples/problem.json
python3 -m gaussian_ts define --problem examples/problem.json \
  --config examples/hydrogen_exchange/reaction.json
```

define 返回 specified_inputs、missing_inputs、ready_for_calculation 和不可绕过的 success_criteria。ready 只表示问题定义字段齐全；不会验证文件是否存在、真实化学正确性、Gaussian 许可或运行主机。init 负责原子/电荷/自旋验证，doctor 与真实主机测试负责环境。

可把同样的 problem 对象放入 reaction.json 的 `problem` 字段，初始化时将问题定义持久化到 session。若没有该字段，context 会明确返回 problem_definition=null，Codex 应补齐研究问题。程序的 HF/STO-3G 等默认参数是接口默认，不能自动视为用户认可的研究设置。

## 准备输入与环境

根据用户给定条件准备 Reactant/Product XYZ、相同的有序原子映射、charge、multiplicity、bond_changes、job_defaults 和 limits。对缺失或矛盾条件，优先询问影响物理问题的内容；结构图、日志、构象生成等独立工作可继续。此阶段不自动猜映射。

执行 doctor 检查 Gaussian 可执行文件；在装有 Gaussian 的计算主机加载环境，再 init 新的运行目录。无需配置第二个 LLM：当前 Codex 就是决策模型。

## 调用工具

```bash
python3 -m gaussian_ts context --run-dir runs/reaction-001
python3 -m gaussian_ts tools
```

context 包括问题、映射、反应坐标、当前结构、登记结构、作业状态/诊断/模式、历史依据、假设、预算、验证结果和完整工具参数。大型轨迹保留在 parsed.json，日常上下文只带相关结构化反馈。

选择动作后写文件：

```json
{
  "tool": "set_distance",
  "params": {"atoms": [1, 2], "value": 1.5, "moving": [2]},
  "reason": "根据上一轮诊断，把目标转移原子朝反应中心移动；保留其余锚点"
}
```

```bash
python3 -m gaussian_ts apply --run-dir runs/reaction-001 \
  --action ACTION_JSON --actor codex
python3 -m gaussian_ts context --run-dir runs/reaction-001
```

Codex 自行选择下一动作，既可通过 CLI，也可直接使用 `Session.execute`。每动作有依据；记录计算反馈后再决定调整几何、约束、数值策略、验证或切换候选。`use_structure` 和保存的 hypothesis 支持串行探索分支。

`tools --format functions` 输出常见 function-tools JSON 描述，可给自己的 GPT 宿主注册。该输出是工具 schema；宿主仍须把选中的工具名/参数/依据分发到 `Session.execute`。当前未注册全局 MCP 服务；Codex 直接用 shell/Python 调用即可。

## 判定与汇报

```bash
python3 -m gaussian_ts report --run-dir runs/reaction-001
```

引用 report 的真实标签、分项依据和关联作业，说明未完成的验证及预算消耗。原始日志路径在 context/jobs 中；可用它们进一步理解失败，但不能以语言推断覆盖独立验证结果。

当前端点身份/模式检查支持人工给定的目标键集合，并非通用复杂机理识别。证据不适配当前规则时明确为不确定，保留原始结果供人工/扩展验证器分析。

独立 HTTPPolicy/Agent 是可选的宿主集成和测试入口；主工作流不需要它们。研究级验收仍需要真实 Gaussian、真实反应数据及完整频率/IRC/端点证据。
