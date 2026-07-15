"""Single-language message table for the interactive lipid agent.

交互式脂质 agent 的单语言消息表。

After the user picks a language at step 1, every prompt and progress message is shown ONLY in that
language (not bilingual). This module holds every string in both English and Chinese and returns the
chosen one.

用户在第 1 步选择语言后，所有提示与过程信息都只用该语言显示（不再中英并列）。本模块保存所有中英文
字符串并返回所选语言的那一份。
"""


MESSAGES = {
    "ask_language": {
        "en": "Which language should we use for this session? Type 'en' for English or 'zh' for Chinese: ",
        "zh": "本次对话使用哪种语言？输入 'en' 用英文，输入 'zh' 用中文：",
    },
    "lang_set": {
        "en": "Language set to English. All following messages will be in English.",
        "zh": "语言已设为中文，后续所有信息都用中文显示。",
    },
    "banner_params": {
        "en": "=== Step 2 | Parameters and data-folder detection ===",
        "zh": "=== 第二步 | 参数与数据文件夹检测 ===",
    },
    "ask_data_root": {
        "en": "Enter the data root folder (it should contain pos/ and/or neg/ subfolders of mzML): ",
        "zh": "请输入数据根目录（其中应包含 pos/ 和/或 neg/ 存放 mzML 的子文件夹）：",
    },
    "ask_library_dir": {
        "en": "Enter the spectral library folder (with pos_common.pkl / neg_common.pkl): ",
        "zh": "请输入谱图库文件夹（内含 pos_common.pkl / neg_common.pkl）：",
    },
    "detect_found_both": {
        "en": "Detected BOTH polarities: pos ({n_pos} mzML) and neg ({n_neg} mzML).",
        "zh": "检测到两种极性：pos（{n_pos} 个 mzML）和 neg（{n_neg} 个 mzML）。",
    },
    "detect_found_one": {
        "en": "Detected ONLY {polarity} ({n} mzML files).",
        "zh": "只检测到 {polarity}（{n} 个 mzML 文件）。",
    },
    "detect_found_none": {
        "en": "No mzML files were found under pos/ or neg/ of the given root.",
        "zh": "在该根目录的 pos/ 或 neg/ 下都没有找到 mzML 文件。",
    },
    "ask_run_both": {
        "en": "Run both pos and neg? Type 'both', 'pos', or 'neg': ",
        "zh": "要都跑 pos 和 neg 吗？输入 'both'、'pos' 或 'neg'：",
    },
    "ask_run_one": {
        "en": "Only {polarity} is available. Run it? Type 'y' to run or 'n' to abort: ",
        "zh": "只有 {polarity} 可用，是否运行？输入 'y' 运行，'n' 取消：",
    },
    "params_summary": {
        "en": "Parameters: MS2_filter={MS2_filter}, MS1 ppm={ppm_err1}, MS1 Da={da_err1}, "
              "MS2 ppm={ppm_err2}, MS2 Da={da_err2}, EQ3 workers={eq3_workers}.",
        "zh": "参数：MS2_filter={MS2_filter}，MS1 ppm={ppm_err1}，MS1 Da={da_err1}，"
              "MS2 ppm={ppm_err2}，MS2 Da={da_err2}，EQ3 线程={eq3_workers}。",
    },
    "ask_use_default_params": {
        "en": "Use these default parameters? Type 'y' to accept or 'n' to edit: ",
        "zh": "使用以上默认参数吗？输入 'y' 接受，'n' 修改：",
    },
    "ask_param_value": {
        "en": "  {name} (current {current}): press Enter to keep, or type a new value: ",
        "zh": "  {name}（当前 {current}）：直接回车保留，或输入新值：",
    },
    "banner_pipeline": {
        "en": "=== Step 3 | Running the qualitative + quantitative pipeline ===",
        "zh": "=== 第三步 | 运行定性 + 定量全流程 ===",
    },
    "pipeline_polarity_start": {
        "en": "-- Running polarity: {polarity} --",
        "zh": "-- 正在运行极性：{polarity} --",
    },
    "pipeline_done": {
        "en": "Pipeline finished for {polarities}.",
        "zh": "{polarities} 的流程已完成。",
    },
    "banner_sampleinfo": {
        "en": "=== Step 4 | Sample information and quality control ===",
        "zh": "=== 第四步 | 样本信息与质量控制 ===",
    },
    "detected_samples": {
        "en": "Detected samples ({n}): {samples}",
        "zh": "检测到的样本（{n} 个）：{samples}",
    },
    "auto_groups": {
        "en": "Auto-inferred groups (by name): {groups}",
        "zh": "按名称自动推断的分组：{groups}",
    },
    "ask_edit_groups": {
        "en": "Accept these groups? Type 'y' to accept, or 'n' to type them manually: ",
        "zh": "接受这些分组吗？输入 'y' 接受，'n' 手动填写：",
    },
    "ask_group_for": {
        "en": "  Group label for sample {sample}: ",
        "zh": "  样本 {sample} 的分组标签：",
    },
    "ask_case_control": {
        "en": "Which two groups are the comparison? Enter as CASE,CONTROL (e.g. B,A): ",
        "zh": "用于比较的两组是哪两个？按 病例组,对照组 输入（例如 B,A）：",
    },
    "qc_running": {
        "en": "Running QC (total-area correction, CV, PCA, PLS-DA) ...",
        "zh": "正在进行质控（总峰面积校正、CV、PCA、PLS-DA）……",
    },
    "ask_sample_context": {
        "en": "Describe the samples / study background (disease, tissue, groups, etc.). "
              "This guides the biological conclusions: ",
        "zh": "请描述样本/研究背景（疾病、组织、分组等），用于指导生物学结论生成：",
    },
    "banner_filter": {
        "en": "=== Annotation & quantification filtering (applied before differential analysis) ===",
        "zh": "=== 注释与定量筛选（在差异分析之前应用）===",
    },
    "ask_include_ms1only": {
        "en": "Include MS1-only inferred annotations in the analysis? (y/n): ",
        "zh": "是否让仅有 MS1 的推断注释参与分析？(y/n)：",
    },
    "ask_include_unannotated": {
        "en": "Include MS1 unannotated-search annotations? (y/n): ",
        "zh": "是否让 MS1 未注释搜库结果参与分析？(y/n)：",
    },
    "ask_custom_grades": {
        "en": "Restrict to confidence grades? Enter a comma list of 1/2/3 (e.g. 1,2), "
              "or press Enter to keep all: ",
        "zh": "是否限定置信度等级？输入 1/2/3 的逗号列表（如 1,2），或直接回车表示全部：",
    },
    "ask_cv_filter": {
        "en": "Drop features whose QC CV is too high? (y/n): ",
        "zh": "是否剔除 QC CV 太差（过高）的特征？(y/n)：",
    },
    "ask_cv_threshold": {
        "en": "  Maximum acceptable QC CV percent (default 30): ",
        "zh": "  可接受的 QC CV 上限（百分数，默认 30）：",
    },
    "filter_applied": {
        "en": "Filtering ({stage}): {n0} -> {n1} features ({dropped} removed).",
        "zh": "筛选（{stage}）：{n0} -> {n1} 个特征（移除 {dropped} 个）。",
    },
    "confirm_conclusions_header": {
        "en": "The following biological conclusions were generated:",
        "zh": "已生成以下生物学结论：",
    },
    "ask_confirm_conclusions": {
        "en": "Accept these for verification? 'y'=accept all, 'r'=regenerate, "
              "'q'=skip verification, or numbers to drop (e.g. 2,4): ",
        "zh": "确认进入验证吗？'y'=全部接受，'r'=重新生成，'q'=跳过验证，"
              "或输入要删除的编号（如 2,4）：",
    },
    "conclusions_updated": {
        "en": "Kept {n} conclusions.",
        "zh": "保留 {n} 条结论。",
    },
    "no_conclusions_kept": {
        "en": "No conclusions kept; skipping verification.",
        "zh": "未保留任何结论，跳过验证。",
    },
    "banner_analysis": {
        "en": "=== Step 5 | Lipid differential analysis ===",
        "zh": "=== 第五步 | 脂质差异分析 ===",
    },
    "analysis_running": {
        "en": "Building volcano plot, heatmap, FA-chain change, class bubble, and random-forest classifier ...",
        "zh": "正在生成火山图、热图、脂肪酸链变化、脂质分类气泡图和随机森林分类器……",
    },
    "banner_conclusions": {
        "en": "=== Step 6 | Biological conclusions and multi-agent verification ===",
        "zh": "=== 第六步 | 生物学结论与多智能体验证 ===",
    },
    "conclusions_generated": {
        "en": "Generated {n} biological conclusions; launching {n} verification sub-agents ...",
        "zh": "已生成 {n} 条生物学结论；正在启动 {n} 个验证子智能体……",
    },
    "banner_report": {
        "en": "=== Step 7 | Writing the full Markdown report ===",
        "zh": "=== 第七步 | 生成完整 Markdown 报告 ===",
    },
    "report_saved": {
        "en": "Full report saved to: {path}",
        "zh": "完整报告已保存到：{path}",
    },
    "ask_llm_provider": {
        "en": "LLM provider for conclusions/verification (anthropic / openai / deepseek / custom / none): ",
        "zh": "用于结论/验证的大模型提供方（anthropic / openai / deepseek / custom / none）：",
    },
    "ask_llm_key": {
        "en": "API key (leave empty to read from environment): ",
        "zh": "API 密钥（留空则从环境变量读取）：",
    },
    "ask_llm_model": {
        "en": "Model name (Enter for the provider default): ",
        "zh": "模型名称（回车使用该提供方默认值）：",
    },
    "ask_llm_baseurl": {
        "en": "Base URL (only for 'custom'; Enter to use the provider default): ",
        "zh": "Base URL（仅 custom 需要；回车使用默认）：",
    },
    "llm_none": {
        "en": "No LLM configured. Conclusions and verification will use the offline rule-based engine.",
        "zh": "未配置大模型，结论与验证将使用离线规则引擎。",
    },
    "aborted": {
        "en": "Aborted by user.",
        "zh": "已被用户取消。",
    },
    "summary_header": {
        "en": "===== Short summary =====",
        "zh": "===== 简要总结 =====",
    },
}


class Messages:
    def __init__(self, language="en"):
        self.language = language if language in ("en", "zh") else "en"

    def t(self, key, **fmt):
        entry = MESSAGES.get(key)
        if entry is None:
            return key
        text = entry.get(self.language, entry.get("en", key))
        if fmt:
            try:
                return text.format(**fmt)
            except (KeyError, IndexError):
                return text
        return text
