# 工具 API（v0.1）

坐标单位 Å，角度 °，Gaussian 频率 cm⁻¹、能量 Hartree。工具原子参数都是稳定的正整数 `atom_ids`；生成 Gaussian 时转换为 1 起始行号。

## 公共接口

```python
from gaussian_ts.engine import Session
from gaussian_ts.models import Structure, Reaction, Limits
from gaussian_ts.agent import Agent, DecisionPolicy

# 初次 Session.create(root, reaction, limits, command=["g16"], job_defaults={...})
session = Session.open("runs/reaction-001")
observation = session.observe()
result = session.execute({
    "tool": "set_distance",
    "params": {"atoms": [1, 2], "value": 1.5, "moving": [2]},
    "reason": "依据当前反应坐标调整初猜",
    "actor": "chemist"
})
report = session.evaluate()
```

`execute` 返回 `ok`、`action_id` 和 `result` 或 `error`。科学/数值失败仍可产生 `ok=true` 的已执行动作，必须继续检查 `job.status`、诊断和验证标签。无效动作也消耗动作预算并保存失败记录；暂停状态下 agent 动作和预算耗尽后的动作在入口拒绝。

## 几何与决策动作

| tool | params | 行为 |
|---|---|---|
| `inspect_geometry` | `{}` | 当前坐标及粗碰撞检查 |
| `set_distance` | `atoms:[a,b], value, moving:[IDs]` | a 固定，移动含 b 的片段；缺省只移动 b |
| `set_angle` | `atoms:[a,b,c], value, moving:[IDs]` | 绕 b 旋转含 c 的片段，a/b 固定 |
| `set_dihedral` | `atoms:[a,b,c,d], value, moving:[IDs]` | 绕 b→c 轴旋转含 d 的片段，a/b/c 固定 |
| `translate_fragment` | `moving:[IDs], vector:[dx,dy,dz]` | 片段平移 |
| `rotate_fragment` | `moving:[IDs], origin:[x,y,z], axis:[x,y,z], degrees` | Rodrigues 刚体旋转 |
| `use_structure` | `structure_id` | 从登记结构分支/切换；不接受任意文件路径 |
| `undo_geometry_edit` | `{}` | 恢复上一个当前结构（也可撤销切换/优化输出） |
| `save_hypothesis` | `hypothesis:string` | 保存假设，关联结构及动作 |
| `request_human_review` | `question:string` | 记录审阅问题并暂停 |
| `pause` / `resume` | `{}` | 检查点暂停/继续 |
| `verify_ts` | `ts_job_id`（可选） | 根据关联原始证据验证选定或最新 TS |

移动片段必须排除几何参考锚点。不存在自动推断移动原子；LLM/用户显式给出 fragment。当前粗碰撞阈值为 0.5 Å，不能替代成键、价态或局部化学合理性判断。

## Gaussian 动作

```json
{
  "tool": "run_gaussian",
  "params": {
    "kind": "freq",
    "structure_id": "job_0001_output",
    "parent_job": "job_0001",
    "method": "UHF",
    "basis": "STO-3G",
    "timeout_seconds": 300
  },
  "reason": "在已收敛的无约束 TS 上检查完整频率谱"
}
```

`structure_id` 缺省当前结构，`parent_job` 缺省该结构的来源作业。Freq/IRC 的父作业必须为 TS/QST；endpoint_forward/reverse 的父作业必须为相应方向 IRC。即使父作业不是已验证 TS，也允许收集后续证据；最终独立验证器统一裁决。

| 字段 | 默认 | 含义 |
|---|---|---|
| `kind` | 必填 | opt/frozen_opt/scan/ts/qst2/qst3/freq/irc_forward/irc_reverse/endpoint_forward/endpoint_reverse |
| `method` / `basis` | HF / STO-3G | 由 session.job_defaults 覆盖；只允许安全 token，不支持 Gen/GenECP |
| `cpus` / `memory_mb` | 2 / 1024 | Gaussian 请求，受 Limits 限制 |
| `max_cycles` | 100 | Opt MaxCycles |
| `irc_max_points` / `irc_step_size` | 50 / 10 | Gaussian MaxPoints/StepSize；StepSize 保持 Gaussian 原生整数单位 |
| `scf_xqc` | false | 显式受控 SCF 恢复选项 |
| `solvent` | 空 | 无溶剂，或 Water/Acetonitrile/Methanol/Ethanol/Toluene/Dichloromethane 的 SMD |
| `constraints` | 空 | 只允许 frozen_opt / scan；TS/QST/Freq/IRC 不接受约束 |
| `timeout_seconds` | Limits.job_timeout_seconds | 必须在每作业时限内且完整预留额度适配总预算 |

所有输入使用 `NoSymm` 和 `Integral=UltraFine`；TS 默认 `CalcFC` 和 `Tight`，保持本征值检查。

冻结键示例：

```json
{"kind":"frozen_opt","constraints":[{"type":"B","atoms":[1,2],"operation":"F"}]}
```

松弛扫描示例：

```json
{"kind":"scan","constraints":[{"type":"B","atoms":[1,2],"operation":"S","steps":10,"step_size":0.1}]}
```

冻结当前值；若要冻结另一数值先编辑结构。A/D 坐标分别要求 3/4 个原子。扫描可同时带冻结约束，但至少一个 S 坐标；Freeze 不接受 steps/step_size。

输出保存即使计算失败时得到的结构版本；只有成功且优化判据满足的优化作业会把输出设为当前结构。Freq/IRC 输出保留但不会自动改变当前结构。结构文件名来自内部编号，作业路径固定在 session 根目录。

## 解析与执行接口

```python
from gaussian_ts.gaussian import JobSpec, prepare_input
from gaussian_ts.parser import parse_log
from gaussian_ts.runner import LocalRunner

input_text = prepare_input(reaction, structure, JobSpec("ts"))
runner = LocalRunner(["g16"])
handle = runner.start(job_directory, input_text)  # 新目录，禁止覆盖
status = runner.get_job_status(handle)
runner.cancel(handle)                          # 只能操作本 runner 持有的 handle
result = runner.wait(handle, timeout=300)
parsed = parse_log(raw_text, atom_ids=structure.atom_ids)
```

`run(directory, input_text, timeout)` 是 start/wait 的同步便利入口。状态有 completed/failed/timed_out/cancelled，进程 returncode=0 仅表示进程层面完成；Session 继续核查正常结束。

`parse_log` 返回完整轨迹、最终结构、能量序列、收敛字典、频率、位移模式、IRC 点/完成状态、diagnostics、parse_warnings、原始文本 SHA256。不支持的日志布局保守地缺失证据；不要拿压缩 observation 替代原始 parsed.json。

## 策略适配

实现 `decide(observation)` 返回动作字典或 None。`Agent(session, policy).run(max_steps=20)` 每轮重新观察 checkpoint，再执行一个动作；None 表示策略结束，不能当作成功。Agent 强制操作者为 `agent`，模型不能冒充人工跳过暂停。

历史动作的依据、失败与假设在 observation 中可见；原始日志保留在作业目录，结构化数据用于日常决策。供应商适配器应使用已校验的工具参数，不扩展成任意代码执行。

## 验证标签

- `VALIDATED_TS`：全部真实证据满足当前验证协议。
- `CANDIDATE_TS`：驻点/频率/模式满足，但完整 IRC 或端点证据不足。
- `WRONG_CHANNEL`：模式或两端与目标反应不匹配。
- `FAILED_SEARCH`：目标候选 TS 作业失败/超时。
- `INCONCLUSIVE`：无候选、模拟数据、结构/水平/模式证据不完整或不确定。

原始日志被改写后不能复用先前通过的证据。完整报告包含逐项检查、关联作业、预算、假设、错误和明确限制；不提供没有 thermochemistry 证据的能垒数值。
