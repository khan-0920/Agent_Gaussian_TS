# Gaussian TS Toolkit for Codex / GPT

按照 [项目手册](LLM_Gaussian_TS_Agent_Project_Goals.md) 搭建的 **v0.1 计算化学辅助工具包**。用于帮助当前 Codex/GPT 定义研究问题、读取实验上下文，并调用结构编辑、Gaussian 和证据验证工具。Codex/GPT 负责研究判断，工具包负责确定性执行和记录。

当前交付是人工引导 MVP 的工具基础，尚未通过真实 Gaussian 反应案例的科学验收。自动原子映射、机理发现和通用构象生成属于后续阶段。

## Codex 使用入口

仓库的 [AGENTS.md](AGENTS.md) 提供 Codex 操作指引；具体用法见 [Codex 使用说明](docs/CODEX_GUIDE.md)。主工作流不需要另一个 API key 或独立 Agent 服务。

```bash
# 先定义问题，分开固定条件、机理假设与未知信息
python3 -m gaussian_ts define --problem examples/problem.json \
  --config examples/hydrogen_exchange/reaction.json

# 获取工具清单；也可用 --format functions 导出 function-tools schema
python3 -m gaussian_ts tools

# 初始化后，读取给 Codex 的结构化上下文，并执行它选择的动作
python3 -m gaussian_ts context --run-dir runs/reaction-001
python3 -m gaussian_ts apply --run-dir runs/reaction-001 \
  --action ACTION_JSON --actor codex
```

`define` 检查问题定义是否齐全；`context` 提供映射、结构、反馈、历史依据、预算和验证结果；`apply` 执行一个带依据的动作。Codex 根据新的 context 再选择下一步，无需固定流程。

## 已有能力

| 部分 | 实现 |
|---|---|
| 输入 | XYZ、稳定原子映射、原子守恒、电荷/自旋一致性检查 |
| 几何 | 键长、键角、二面角、片段刚体平移/旋转、碰撞检查、结构版本与撤销 |
| Gaussian | Opt、Frozen Opt、松弛扫描、TS/CalcFC、QST2/QST3、Freq、双向 IRC、端点优化 |
| 执行 | 本机进程、独立作业/scratch 目录、状态/取消 Python API、超时、预算预留、重复计算阻止 |
| 反馈 | SCF 能量、四项优化判据、坐标轨迹、频率/位移模式、IRC 点、失败诊断 |
| Codex 接口 | 问题定义、上下文导出、工具 JSON schema、单动作 CLI/Python 调用；可选独立宿主/回放示例 |
| 记录 | 每动作 `action.json`、原始输入/日志、解析 JSON、结构 XYZ、原子写入检查点、暂停/继续 |
| 验证 | 五类手册标签；同一计算水平和结构、目标振动方向、双向 IRC 与优化端点的关联检查 |

框架不会让 LLM 写 shell 或原始 Gaussian route，不提供 `NoEigenTest` 开关。

## 环境

- Python **3.9+**。核心代码和测试只使用标准库；可直接运行，无须先安装依赖。
- 真实计算需要用户已有的合法 Gaussian 09/16 安装、许可、环境和计算资源。框架不下载或安装 Gaussian。
- 在计算主机上先按安装说明加载 Gaussian 环境，再启动框架。`gaussian_command` 是可信配置中的 **argv 数组**；如果站点需要初始化脚本，请由用户创建/检查包装程序，并将其绝对路径放在该数组中。
- 真实 LLM 需要支持 JSON 输出的 chat-completions 兼容接口；配置完整 endpoint 和模型名。

```bash
python3 -m gaussian_ts doctor
python3 -m unittest discover -v
```

测试包含一个绑定 `127.0.0.1` 临时端口的 HTTP 接口集成测试；受限执行环境需允许本机端口绑定。测试不调用外部 LLM 或真实 Gaussian。

可选安装（需要 setuptools 构建环境）：

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
gaussian-ts --help
```

## 先运行离线演示

```bash
python3 -m examples.offline_demo --run-dir runs/offline-demo
python3 -m gaussian_ts status --run-dir runs/offline-demo
python3 -m gaussian_ts report --run-dir runs/offline-demo
```

演示通过独立进程产生明确标记的**模拟日志**：第一次 TS 尝试失败，策略读到 Hessian 反馈后修改几何，再收集频率、正反 IRC 和端点优化。它证明动作、进程、反馈和记录接口可以联通；演示策略是确定性示例策略，并非真实 LLM。

所有演示作业标记 `synthetic: true`。最终应为 `INCONCLUSIVE`，永远不能成为 `VALIDATED_TS`。`examples/hydrogen_exchange` 是映射和调用格式示例，结构/能量不是优化后的基准数据。

运行目录必须不存在；重复演示请使用新的目录，避免覆盖已有实验。

## 人工引导的真实计算

先修改 `examples/hydrogen_exchange/reaction.json`，换成你的反应结构、映射、成断键、方法/基组、Gaussian 命令及预算。路径相对于配置文件。`change=-1` 表示断键，`change=1` 表示成键；原子 ID 按 `atom_ids` 映射，不是随意改变的坐标行号。

```bash
python3 -m gaussian_ts init \
  --config examples/hydrogen_exchange/reaction.json \
  --run-dir runs/reaction-001

python3 -m gaussian_ts apply \
  --run-dir runs/reaction-001 --action examples/edit_guess.json

python3 -m gaussian_ts apply \
  --run-dir runs/reaction-001 --action examples/run_ts.json

python3 -m gaussian_ts report --run-dir runs/reaction-001
```

`apply` 的 `ok` 表示工具动作执行成功；Gaussian 是否计算成功需查看 `result.job.status` 和 `parsed`。只跑 TS 优化不能完成验证。

Freq 和双向 IRC 必须选择 TS 作业的 `output_structure_id` 和 `parent_job`，端点优化选择相应 IRC 的输出。输入结构与父作业输出及计算水平不一致时拒绝提交；从状态输出获取实际作业 ID，详见 [工具 API](docs/API.md)。

## 可选：独立 LLM 宿主示例

本节供另行构建 LLM 宿主时使用；当前 Codex 直接通过上述工具入口操作即可。

复制 `examples/llm.example.json` 到忽略的 `config.local.json`，填写你的 endpoint 和模型名。在环境变量中设置 `TS_LLM_API_KEY`；框架不会自动读取 `.env`，不要把密钥写入 JSON 或提交到 Git。

```bash
python3 -m gaussian_ts run \
  --run-dir runs/reaction-001 --llm-config config.local.json --steps 15
```

LLM 每轮接收当前结构、反应键、已有作业反馈、历史动作、假设和预算，返回一个带依据的 JSON 动作。模型可以选择编辑、扫描、TS、频率、IRC、回退或人工审阅。provider 出错会暂停并留下审阅请求。

API 返回格式为 `choices[0].message.content`，其中 content 是一个 JSON 对象字符串。请求包含 `response_format: {"type":"json_object"}`；不支持该格式的服务需自行实现 `DecisionPolicy` 适配器。LLM 每次 HTTP 请求默认上限 45 秒/1500 输出 tokens，轮数由 `--steps` 和动作预算限制；当前不统计供应商费用或输入 token 总额。

显式人工动作计划可用 `run --actions path/to/actions.json` 回放（JSON 数组）；回放是复现/基线接口，不能当作自主 LLM 策略。

## 人工接管与恢复

```bash
python3 -m gaussian_ts pause --run-dir runs/reaction-001
python3 -m gaussian_ts apply --run-dir runs/reaction-001 \
  --action examples/edit_guess.json --actor chemist
python3 -m gaussian_ts resume --run-dir runs/reaction-001
```

重新运行 `run` 即从 checkpoint 继续。`use_structure` 可选择任意已登记结构建立候选分支，`undo_geometry_edit` 恢复前一个结构。每次覆盖保留操作者及前后结构。

当前执行器是**串行、前台**本机执行；人机接管发生在动作检查点。Python `LocalRunner` 提供启动、查询、等待和取消自己作业的接口，CLI 不提供跨进程后台取消。中断时终止当前执行器拥有的进程并保守保留预算预留；若整个主进程被强制杀死，重开 session 会标记该动作/作业 interrupted，不盲目重新提交或杀其他进程。恢复前需在计算主机检查遗留进程。

保存 `job.chk` 以便追溯和后续扩展。当前 checkpoint 恢复是 **Agent 状态恢复**；二进制 Gaussian `Opt=Restart`、SLURM 和并行执行尚未实现。

## 证据规则与范围

`VALIDATED_TS` 要求所有证据都来自真实、成功结束的同一计算水平作业：无约束 TS 四项优化判据通过、完整频率谱中仅一个频率 ≤ −30 cm⁻¹ 的显著负曲率、其位移与所有目标成断键方向一致、完成双向多点 IRC、优化两端并匹配反应物/产物。额外小负频不会被悄悄忽略。

模式分析采用目标键长导数投影；端点身份采用固定映射的全距离矩阵（RMS ≤ 0.25 Å）和关键键长（偏差 ≤ 0.20 Å），容许两端交换。它们是当前人工指定反应中心的保守检查器，尚不是通用机理/化学身份识别器。不同构象、复合物解离和复杂重排可能需要人工审阅或更好的身份适配器。

验证时重新读取原始日志并核对摘要；日志被修改则拒绝原证据。日志解析目前只支持单作业、常规 Cartesian orientation、标准频率/位移表；Link1、缺失数据或不认识的 IRC 完成标记会留下不确定性。方法、基组、SMD 溶剂等来自生成输入的元数据，首次接入具体 Gaussian 版本时仍需检查实际输出是否符合目标设置。

预算包括动作次数、Gaussian 作业数、累计作业墙钟及**分配核秒（墙钟×请求核数）**。提交前预留完整超时额度；实际结束后结算，不把失去跟踪的作业当作零消耗。`%nprocshared`/`%mem` 限制 Gaussian 的每作业请求，并不是对任意包装程序的 OS 级内存/CPU 隔离。没有 ZPE/热修正/参考态证据时不会报告能垒或自由能。

## 运行产物

```text
runs/reaction-001/
  state.json                 # 可恢复状态、预算、结构与作业关系
  report.json                # 独立验证器结果、证据、失败与限制
  structures/*.xyz           # 映射见 state.json
  actions/000001/action.json # 依据、操作者、输入/输出、错误和资源
  jobs/job_0001/
    input.gjf
    output.log
    parsed.json
    job.json
    job.chk                  # 真实 Gaussian 生成时存在
    scratch/
```

默认 Git 忽略运行目录、密钥、本地配置和 Gaussian 二进制文件；需要归档真实基准时，应专门筛选并登记数据来源。

后续优先工作：在实际 Gaussian 主机用精选已知反应完成全链验收，补入脱敏真实日志回归集，再扩展自动映射/构象、调度器和标准基准评估。详见 [设计](docs/superpowers/specs/2026-10-08-gaussian-ts-foundation-design.md) 与 [实现计划](docs/superpowers/plans/2026-10-08-gaussian-ts-foundation.md)。
