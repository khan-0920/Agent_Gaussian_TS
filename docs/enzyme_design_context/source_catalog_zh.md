# 酶设计研究来源目录

按唯一研究计数，镜像、作者稿、期刊网页和代码仓库不另算论文。目录用于找到证据，不代表所有研究已完整审计。

目前共17篇研究：8篇部分正文核查，8篇摘要或元数据，1篇摘要及扩展图注。路线表覆盖13类综合路径，类别可以重叠。

## 研究与核查范围

| 编号 | 研究 | 年份 | 核查范围 | 使用重点 |
| --- | --- | --- | --- | --- |
| E001 | [Kemp elimination catalysts by computational enzyme design](https://www.nature.com/articles/nature06879) | 2008 | 摘要或元数据 | 经典催化模型匹配；优化前后分列 |
| E002 | [De novo computational design of retro-aldol enzymes](https://pubmed.ncbi.nlm.nih.gov/18323453/) | 2008 | 摘要或元数据 | 逆醛缩酶；详细实验待核查 |
| E003 | [Computational design of an enzyme catalyst for a stereoselective bimolecular Diels-Alder reaction](https://pmc.ncbi.nlm.nih.gov/articles/PMC3241958/) | 2010 | 部分正文已核查 | 双底物与立体选择性；分母区分订购和可溶表达 |
| E004 | [Precision is essential for efficient catalysis in an evolved Kemp eliminase](https://pubmed.ncbi.nlm.nih.gov/24132235/) | 2013 | 摘要或元数据 | 进化后精确几何；不归给初始设计 |
| E005 | [Computational design of serine hydrolases](https://pmc.ncbi.nlm.nih.gov/articles/PMC12288761/) | 2025 | 部分正文已核查 | RFD丝氨酸水解；区分单次反应与周转 |
| E006 | [Computational design of metallohydrolases](https://www.nature.com/articles/s41586-025-09746-w) | 2025 | 部分正文已核查 | Zn水解；配位变化和动力学条件 |
| E007 | [Large language models generate functional protein sequences across diverse families](https://doi.org/10.1038/s41587-022-01618-2) | 2023 | 部分正文已核查 | 家族生成；新实验与既有数据评估分开 |
| E008 | [Machine learning-assisted directed protein evolution with combinatorial libraries](https://pmc.ncbi.nlm.nih.gov/articles/PMC6500146/) | 2019 | 部分正文已核查 | ML进化；结合基准与催化实验分开 |
| E009 | [Biocatalytic asymmetric synthesis of chiral amines from ketones applied to sitagliptin manufacture](https://pubmed.ncbi.nlm.nih.gov/20558668/) | 2010 | 摘要或元数据 | 底物逐步适配后进化；工艺指标单列 |
| E010 | [Design and evolution of an enzyme with a non-canonical organocatalytic mechanism](https://www.nature.com/articles/s41586-019-1262-8) | 2019 | 摘要或元数据 | 非天然氨基酸催化；细节待全文 |
| E011 | [Directed evolution of artificial metalloenzymes for in vivo metathesis](https://www.nature.com/articles/nature19114) | 2016 | 摘要及扩展图注 | 人工金属辅因子；游离辅因子对照 |
| E012 | [Engineering highly functional thermostable proteins using ancestral sequence reconstruction](https://www.nature.com/articles/s41929-018-0159-5) | 2018 | 摘要或元数据 | 祖先重建；热稳定与活性分开 |
| E013 | [Automated Structure- and Sequence-Based Design of Proteins for High Bacterial Expression and Stability](https://pubmed.ncbi.nlm.nih.gov/27425410/) | 2016 | 摘要或元数据 | 表达和稳定性工程；脚本更正已核查 |
| E014 | [De novo Design of All-atom Biomolecular Interactions with RFdiffusion3](https://pmc.ncbi.nlm.nih.gov/articles/PMC12458353/) | 2025 | 部分正文已核查 | RFD3预印本；部分补充方法已读 |
| E015 | [Atom-level enzyme active site scaffolding using RFdiffusion2](https://www.nature.com/articles/s41592-025-02975-x) | 2026 | 部分正文已核查 | RFD2基准；与E006共享金属酶实验 |
| E016 | [De novo design of luciferases using deep learning](https://www.nature.com/articles/s41586-023-05696-3) | 2023 | 部分正文已核查 | 生成骨架库与RIF；文库覆盖未知 |
| E017 | [Computational design of metalloproteases](https://www.nature.com/articles/s41586-026-11167-2) | 2026 | 摘要或元数据 | 金属蛋白酶；模型分支和正文待核查 |

## 计数与复用规则

E006与E015的共享金属水解酶实验不能重复累计。E007对已有实验数据的预测也不能计为新生成序列的湿实验。其余跨论文实验重叠尚未全面审计，独立实验总数留空。

正文已读并不等于补充材料、原始数据和代码均已复核。RFD3预印本版本、论文更正、初始设计与优化变体均保留在source_register.json中。

后续使用时先读enzyme_design_guide_zh.md，再按本表进入原始研究。精确数字、分母、待核查项和代码提交分别见source_register.json及rfd_code_provenance.json。
