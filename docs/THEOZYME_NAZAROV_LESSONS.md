# Nazarov Theozyme 构建过程与经验

记录日期：2026-10-12。用途：为未来反应提供决策依据，不把此案例当所有体系的参数模板。复用入口：[指导Prompt](THEOZYME_REUSABLE_PROMPT.md)。

## 研究演进

最初目标是从中性二乙烯基酮出发，探索酸活化及关环，而不是只计算预先质子化底物。先研究Glu侧链：保留Cβ–Cγ–Cδ–COOH、去主链、断口补氢，无水，中性单重态。经历摆放、受限扫描、无约束TS、频率、双向IRC和端点验证。成环端点出现限制性到非限制性电子不稳定，说明找到几何鞍点并不等于电子态验证完成。用户随后明确停止Glu路线。

并行AlCl3对照采用单体AlCl3与羰基O配位，无水气相、中性单重态。同水平结果保存在 `/data2/guest/enzyme_design/nazarov_alcl3/validated_ring_closure`；历史记录报告相对配位反应物的ΔE‡=26.93355、ΔH‡=25.86718、ΔG‡=27.55267 kcal/mol（298.15K、1atm、RRHO）。复用数值时重新查包内原始证据。该模型不包含Al2Cl6、溶剂物种、产物释放或完整周转。

随后建立ZnII–His3模型：49原子，二乙烯基酮+Zn2++三个中性4-methylimidazole侧链；保留His CB及咪唑环、去主链补氢，总电荷+2、单重态，无显式水。Zn13、O6、供体N19/N31/N43，关环C1/C5。统一wB97XD/6-31+G(d)、UltraFine，8CPU/8GB，单作业7200秒。并发从一个按用户授权增加到最多两个；这些都是本案例设置，不是其他金属体系通用推荐。

## Zn路线与关键判断

1. 手动配位摆放后独立优化反应物，得到141个正频率；C1–C5约3.170Å，Zn–O约2.000Å。反应物电子稳定性后续通过。
2. 以C1–C5进行松弛扫描，8个收敛点约从3.170到1.770Å，最大扫描电子能差约25.83kcal/mol，峰附近约1.970Å。超时保留原作业与已收敛点，再接续，未把扫描峰当作最终TS。
3. 初次TS出现“Desired1 Actual2”负特征值错误；重算Hessian后曾正常结束但位移未达严格阈值，因此进一步精修，没有放宽门槛。用户允许必要时NoEigenTest；这只是中途搜索策略，不能替代最终验证。
4. ts_polish_01_session/job_0001严格四项收敛；job_0002完整141模，唯一虚频−567.9266cm−1，目标关环模式。C1–C5约2.05639Å，Zn–O约1.95047Å。job_0005电子稳定性通过。
5. 正反IRC均因时限分段checkpoint续算：正向最终119点检测极小值，反向134点。原session超时job仍保留超时状态；完整续算证据与验证器整合尚需审计，不能改写原job状态。
6. 两端无约束优化及限制性频率完成。成环端点出现微小R→U不稳定。GuessMix优化曾坍缩为限制性解，因此改读Stable=Opt得到的U checkpoint。U核优化后141频率均正，最低23.2906cm−1，S²约0.0270，固定几何电子稳定性通过。
7. 独立反应物来源是初始模型；独立成环参考来源是预IRC扫描点07，未复制IRC端点。独立成环限制性优化四项通过。与对应限制性IRC端点的距离矩阵差RMS：反应物0.001568Å，成环0.002297Å。仅说明几何匹配，不证明电子态或极小值全套验证完成。
8. 独立成环参考也有R→U不稳定；能量降低约0.003715kcal/mol，S²约0.0209。启动读取其自身U checkpoint的独立核优化，随后用户结束任务，该作业被停止；没有借用IRC端点轨道冒充独立来源。
9. 初步端基二面角与正向路径支持顺转，但完整双向续算路径的顺转审计、R/U电子分支衔接及最终独立参考验证未完成。不能将初步判断写为全套通过。

## 停止时状态与可追溯文件

用户于2026-10-12要求结束任务。自动监测已暂停，运行中的独立U成环参考优化已停止，已有证据保留。记录为STOPPED_BY_USER，overall_theozyme_valid=false；“任务结束”与“科学验证通过”分开。

主目录：`/data2/guest/enzyme_design/nazarov_zn_his3`。

- `status.json`、`history.jsonl`：状态和决策链。
- `ts_polish_01_session`：TS、频率、原始IRC和稳定性。
- `irc_forward_restart_01/02`、`irc_reverse_restart_01/02`及manifest：续算。
- `forward_u_read_opt_01`、`forward_u_read_freq_01`、`forward_u_read_stable_01`：U端点。
- `independent_product_session`、`independent_product_u_opt_01`：独立成环参考。
- `independent_geometry_comparison.json`、`independent_product_electronic_audit.json`：参考匹配和电子态证据。
- `conrotation_audit.json`、`rotation_raw_audit.json`：初步顺转检查，不是完整最终认证。

## 最值得复用的经验

正常终止不保证四项收敛；单虚频不保证正确反应；几何相近不保证同一电子分支；Stable=Opt不等于核几何优化；GuessMix不保证保留U解；独立参考不能来自IRC端点；气相催化簇单步势垒不能直接解释酶催化效率。以上问题要在状态中显式保留，不因用户希望尽快得到结果而隐藏。
