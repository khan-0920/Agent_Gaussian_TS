# Theozyme 到候选酶的工具接口

`enzyme_design` 为当前 Codex/GPT 提供可恢复的计算工具。它从明确的 PDB/mmCIF 催化模型生成骨架、设计序列、优化、独立回预测和采样活性位点，最后交付筛选表。模型推理使用各自已有的 Python 环境；编排、PDB 几何、测试和 CLI 使用 Python 3.9+ 标准库。不需要另一个 LLM 服务。

该模块实现的是计算筛选。`COMPUTATION_COMPLETE` 只表示流程执行完毕；`NO_CANDIDATES_PASSED` 是有效筛选结果。输出绝不自动成为已验证过渡态或有催化活性的酶。Gaussian 的 `VALIDATED_TS` 标准保持独立。

## 输入与版本

完整示例见 [campaign.json](../examples/enzyme_smoke/campaign.json) 和 [服务器配置](../examples/enzyme_smoke/profile.json)。这是 Zn–His3 假设几何的工程实例，不是可直接用于新反应的科学默认设置。配体 SDF 使用稳定反应物键图，不能用它表示部分成断键的 TS 电子结构；Rosetta 参数未经物理标定。

新任务应明确以下字段，缺失身份或映射时停止，而不是按最近距离猜测：

| 字段 | 内容 |
|---|---|
| `theozyme` | 单模型 PDB/mmCIF；元素、链、残基编号、原子名明确 |
| `evidence` | `hypothesis`、`external_evidence` 或 `gaussian_session`，说明和可选原始报告路径；标签本身不证明 TS |
| `motif` | 源残基 ID → 固定原子名，如 `A12:NE2` |
| `ligands` | 每个组分的源残基、名称、CCD 或 SDF、PLACER 是否预测；SDF 必须给出全部重原子的零基索引映射 |
| `rosetta_params` / `rosetta_atom_map` | 自定义配体的显式 Rosetta 类型文件与原子名映射；由研究者核查化学类型、电荷和质子化 |
| `constraints` | 距离 Å 或角度 ° 的上下界，明确适用的 `design` / 预测状态 |
| `states` | 化学解释和配体子集；不同连接状态使用不同显式配体。apo 状态需设 `assess: false` 且仍提供内部几何约束 |
| `generation` | 长度区间、骨架数、种子、`fast`/`legacy`、步数及生成时保留的配体 |
| `sequence` / `relax` / `refold` / `assess` | 每骨架序列数和温度、优化次数及坐标约束、回预测循环与步数、集合样本数 |
| `selection` | pLDDT、主链 RMSD、集合几何通过率、序列身份阈值、交付数量；可选每骨架上限和最小结构 RMSD |
| `limits` | 尝试总数、每阶段超时、总墙钟预算；失败和重试也计入 |

2026-10-19 guest 环境核查：旧 `miniconda3/envs/rfd3` 为 rc-foundry 0.2.0；新版位于 `/data2/guest/rfd3_foundry1_20261018/env_validated`，为 rc-foundry 1.0.0 / torch 2.11.0 / atomworks 3.0.1。适配器使用新版 Python 和 `rfd3` 绝对路径，执行前核验版本；保留 venv Python 的符号链接路径，不能解析到旧基础解释器。RF3 使用另一个独立环境，不随 RFD3 自动升级。

`profile.json` 是可信站点配置，含各阶段 Python、模型入口、权重、仓库、GPU、CPU 和允许的缓存环境变量。它不应包含账号、令牌或密码。示例中的 guest 绝对路径只适用于该服务器，迁移时替换。权重 SHA-256 和实际导入的软件版本写入每次 `provenance.json`；这不等于已下载最新发布的权重。

## CLI 和智能体调用

在仓库根目录运行；安装包后也可使用 `enzyme-design`。已有入口 `python3 -m gaussian_ts enzyme ...` 转发到相同接口。

```bash
python3 -m enzyme_design validate \
  --config examples/enzyme_smoke/campaign.json \
  --profile examples/enzyme_smoke/profile.json
python3 -m enzyme_design init \
  --config examples/enzyme_smoke/campaign.json \
  --profile examples/enzyme_smoke/profile.json \
  --run-dir campaigns/example
python3 -m enzyme_design doctor --run-dir campaigns/example
python3 -m enzyme_design context --run-dir campaigns/example
python3 -m enzyme_design stage prepare --run-dir campaigns/example \
  --reason '核对明确的化学身份、映射和假设几何'
python3 -m enzyme_design run --run-dir campaigns/example \
  --reason '按已指定的样本数和预算执行计算筛选'
python3 -m enzyme_design report --run-dir campaigns/example
```

`validate` 只检查配置契约；`doctor` 只检查路径和可执行权限。模型导入、显存、权重兼容性和真实输出必须由实际运行验证。示例的完整推理需要 guest 上对应 GPU/模型环境，本地标准库测试不会运行模型。

按阶段决策可使用以下动作文件：

```json
{
  "tool": "run_stage",
  "params": {"stage": "generate", "retry": false},
  "reason": "输入身份和设计几何已检查，按指定预算生成骨架"
}
```

```bash
python3 -m enzyme_design apply --run-dir campaigns/example --action action.json
python3 -m enzyme_design tools --format functions
```

宿主将 function 名称/参数与研究依据包装成上述动作，调用 `enzyme_design.cli.execute(Campaign(path), action)`。Python 也可直接调用 `Campaign.create(...)`、`context()`、`run_stage(stage, reason, retry=False)`、`run(reason)` 和 `report()`。工具不接受任意 shell 字符串。

## 阶段证据与筛选

1. `prepare`：保存输入快照，检查重原子身份与 SDF 映射，规范 PDB 名称，生成 RFD3 条件和初始几何报告。
2. `generate`：RFD3 1.0.0 输出原始 CIF、JSON 和重编号映射；检查主链、催化几何及 RFD3 配体碰撞。缺失检查数据不能算通过。
3. `sequence`：LigandMPNN 固定催化序列、使用配体上下文并打包侧链。打包器固定扭转角仍会重建坐标，因此显式保留输入固定残基的笛卡尔坐标，保存原始打包结构及逐原子位移审计。
4. `relax`：PyRosetta 使用明确配体参数和可选催化残基类型，约束 CA、催化原子、配体和指定几何。重查序列/身份/几何。REU 仅作结构优化指标，不解释为自由能或速率。
5. `refold`：RF3 输入序列和化学组分，不提供设计蛋白模板；自定义配体独立生成构象，SDF 原子标签跨移氢过程保留。输出各状态的序列匹配、CA RMSD、模型置信度和催化几何。
6. `assess`：PLACER 每状态采样。显式化学图和原子名保持一致；去除 OpenBabel 在声明为独立组分之间猜出的共价键并记录审计，金属配位保留为几何关系。每个样本单独核对身份和约束。
7. `select`：按用户阈值筛选，重算集合通过比例，记录失败原因与序列/结构多样性比较。交付 FASTA、候选结构、CSV、映射、输入哈希和实验记录模板；实验状态初始为 `not_tested`。

`steps/select/NNNN/` 下的 `candidate_scores.csv`、`rejected.csv`、`all_evaluations.json`、`delivery_manifest.json` 保留无命中情况；空 FASTA 不会被报告为成功设计。RF3、PLACER 原始结构和日志保留在相应阶段目录。

## 恢复、预算与边界

初始化后配置和输入不可变；改变研究条件应新建 campaign。每次执行保留独立尝试目录、源码快照、命令、日志、版本和权重哈希。读取完成结果时复核输出哈希；重复调用完成阶段返回原结果。失败后先读日志，再显式 `--retry`，不会覆盖旧失败证据。控制器异常退出时，对失联尝试保留预留预算；仍活跃的 worker 禁止重复提交。

运行中的控制器持有文件锁。Ctrl-C 和阶段超时会终止所属进程组；`pause`/`resume` 在阶段边界生效，不能用 `pause` 中断正在持锁的模型。执行中不要从另一进程修改运行文件。跨主机恢复应先核实原 worker 状态。

当前输出蛋白限一个链、标准氨基酸；未实现共价蛋白–配体连接、非天然氨基酸、多链复合物及任意插入码/替代构象。输入催化 motif 可以来自多个源链，但生成输出须为单链。所有状态目前使用同一蛋白序列；apo 可回预测，PLACER apo 采样未实现。全局 RMSD 为序列对应 CA 的刚体拟合；结构多样性是序列对齐后的 CA 比较，不是通用折叠分类。未实现独立通道分析、表达预测、实验反馈学习或催化速率预测。

测试：`python3 -m unittest discover -v`。单元/集成测试验证身份、映射、几何、反射排除、固定坐标、输入/输出篡改拒绝、超时清理、锁、预算入口和显式重试；它们不替代 guest 上真实模型运行或实验验证。

## 真实运行记录：2026-10-19

guest 的 `/data2/guest/agent_ts_enzyme_smoke_20261019/run02` 已完成七阶段真实流程。2 个骨架、2 条序列均完成受约束优化、RF3 回预测和每状态 2 次 PLACER 采样。最终 0 个候选通过，阈值未放宽。各阶段原始结果哈希、版本、权重哈希、失败重试和逐候选指标见 [验证摘要](validation/enzyme_smoke_20261019.json)。第一次开发运行 `run01` 的结果也保留在服务器。

| 候选 | RF3 pLDDT（0–1） | CA RMSD（Å） | PLACER 几何通过 | 最终结果 |
|---|---:|---:|---:|---|
| backbone_00001_seq0001 | 0.8048 | 0.799 | 0/2 | 催化几何未恢复，淘汰 |
| backbone_00002_seq0001 | 0.8349 | 2.696 | 0/2 | 主链恢复及催化几何未通过，淘汰 |

验证中发现并修复了 LigandMPNN 固定扭转角后坐标重建、RF3 移氢后的索引、`LIG` 别名触发 CCD 化学替换、PLACER 环境没有 RDKit 四类接口问题。错误输出被检查拒绝，原始失败记录没有被删除或改为通过。最终流程分别使用既有环境；PLACER 自定义配体的 RDKit 准备由独立子进程在 prepare 环境执行。这个实例证明指定单链/非共价组分接口可运行，不代表其他反应、legacy 模式、多状态化学或实验性能已经独立验证。
