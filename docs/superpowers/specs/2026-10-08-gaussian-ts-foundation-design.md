# Gaussian TS 研究框架：首版设计

依据：`LLM_Gaussian_TS_Agent_Project_Goals.md`。用户已要求开始开发和搭建；本次直接落实需求手册中的人工引导 MVP 基础，不额外等待设计授权。

## 目标与验收

提供给 **Codex/GPT 使用的 Python 辅助工具包**：辅助它定义研究问题、区分固定条件/机理假设/未知信息，导出实验上下文，并提供受控动作接口编辑结构、运行 Gaussian 和验证证据。用户提供映射一致的反应物、产物、电荷、自旋和成断键；当前 Codex/GPT 自身是研究决策者。允许回退、分支、暂停和恢复。离线演示只证明工程闭环，不作为真实化学成功证据。

## 架构

- `models.py`：原子映射、XYZ、反应定义和计算/预算配置，使用 Python 3.9+ 标准库。
- `geometry.py`：确定性键长、角度、二面角、片段平移与旋转，保持稳定映射；碰撞检查。
- `gaussian.py`：限定计算类型与白名单参数，生成 Opt、ModRedundant、Scan、TS、QST2/QST3、Freq、IRC 输入。禁止任意 route 和 shell 拼接。
- `parser.py`：能量、优化判据、坐标轨迹、频率和模式、IRC 进度、数值失败；对截断与复合日志保守处理。
- `runner.py`：本机 Gaussian 进程启动、状态、取消、超时，独立工作目录和 scratch。`%mem`/`%nprocshared` 限制每作业请求；不声称对任意外部程序进行操作系统级内存隔离。
- `store.py` / `engine.py`：原子写入状态、文件沙盒、单写者锁、每动作 `action.json`、结构版本、预算预留、重复检测和恢复。
- `validation.py`：独立验证器；同一 TS/计算水平的无约束收敛、唯一显著虚频、目标键变化投影、双向 IRC 与端点优化/身份匹配。
- `toolkit.py` / `cli.py`：问题定义、工具 schema、结构化 context、单动作 CLI/Python 调用；`AGENTS.md` 为 Codex 仓库使用指引。
- `agent.py`：可选宿主示例，包括可注入决策策略、JSON 动作回放、兼容 chat-completions JSON 的 HTTP LLM 接口，不是主产品入口。

## 验证规则

无正常结束、四项收敛判据不全、约束 TS、缺模式数据、微小负频、计算水平/结构不一致、IRC 缺少路径完成证据均不得升级成功。映射保持不变；端点以映射的全距离矩阵和关键反应键共同比较，可容许正反向交换。所有证据关联引擎执行记录和原始日志摘要。测试/演示后端的结果永远不能成为已验证 TS。

## 边界

本次不实现自动原子映射/机理发现、通用构象采样、SLURM、跨进程后台作业恢复、并行候选执行、溶剂与自由能报告。串行候选可用结构版本分支；人工可从 checkpoint 继续并记录操作人。MVP 完整科学验收需真实 Gaussian 日志及已知反应，当前无 Gaussian，不宣称通过该验收。

## 参考

Gaussian 官方关键词文档：https://gaussian.com/opt/ 、https://gaussian.com/irc/ 、https://gaussian.com/freq/ 。关键词页面抓取失败；生成策略交叉核对 ASE 官方 Gaussian 接口 https://wiki.fysik.dtu.dk/ase/ase/calculators/gaussian.html 与 Gaussian 官方振动分析 https://gaussian.com/wp-content/uploads/dl/vib.pdf 。未知日志形式按证据不足处理。
