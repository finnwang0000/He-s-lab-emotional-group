% =========================================================================
% Mvpa_1_0_searchlight_state_singletrial_Aug18.m  —— 状态版 / 单试次(beta)输入
% 变更记录(Aug18, 基于 Aug5 提速, 不改变计算结果):
%   ① cosmo_searchlight 加 nproc 并行开关(n_proc, 默认4, 结果与串行完全一致);
%   ② helper_load_singletrial_sub 把逐 trial 调 cosmo_fmri_dataset 改为按被试一次性
%      批量读取(cellstr 传入), 顺序与原逐个读取严格一致, 只是少了重复的 mask/仿射开销。
%   未改动: SVM 核函数、mask 范围、精度(仍 double) —— 结果应与 Aug5 版逐位可比。
% 变更记录(Aug5): 由 Aug4 searchlight 骨架 + MVPA 参考 Mvpa_1_1_basic_singletrial_Aug04.m
%   的单试次读取逻辑合并而成。输入不再是 con 平均图, 而是每个 trial 的 beta_XXXX.nii,
%   经 beta_info.mat 按【条件名】捞取; 每 trial = 1 个样本(样本量远大于 con 版)。
%   4 状态(全 -clear 口径, 依据 1material\condition名称.md):
%       1 全Focus   : focus-clear
%       2 自然Nature: nature-clear
%       3 屏蔽Shield: {distancing-clear, distraction-clear, inhibition-clear} 三合一(同 target=3)
%       4 闭Close   : close-clear
%   骨架保留 Aug4 的: modal 开关 / dry_run / 时间戳 / dbclear / radius 批量 / LOSO(chunk=被试)。
%   分类器用 cosmo_classify_libsvm(球内, 快); 单试次 beta 须已在与 mask 相同的 MNI 网格。
% ⚠️ 待你确认的路径: single_trial 根目录、全脑 GM mask(必须和 beta 同网格)、每模态 sub_list。
% 小规模复现修改：本文件已改为情绪调节项目的小规模测试版, 改动处均以 [小规模复现修改：] 标注, 汇总如下:
%   (1) 工具箱路径改为本机 CoSMoMVPA; 本机 libsvm 为 3.17(CoSMo 要求>=3.18), 暂不加载;
%   (2) 时间窗 0-9 不存在 → 0-6; radius 只跑 3; n_proc=1(不开并行池); dry_run=false;
%   (3) 4 状态 → 3 条件: natural_viewing / perspective_substitution / cognitive_reappraisal;
%   (4) 数据根目录、mask、被试列表(S1-S3)改为本项目; 被试文件夹名加 'S' 前缀;
%   (5) 结果写到 Codex outputs 目录, 不写入数据盘; 输出文件名改为 3conds_smalltest;
%   (6) 分类器改为 cosmo_classify_lda; 只计算每隔 2000 个体素取 1 个的中心点(center_ids)。
%   (7) helper 中批量读取(cellstr)改回逐 trial 读取, 因本机 CoSMoMVPA 版本不支持 cellstr 输入。
%   正式运行时: 恢复全部被试、删除 center_ids 行、按需改回 libsvm(需装 >=3.18)、换组水平 mask。
% =========================================================================
clear; clc;
dbclear all;   % 清残留断点, 避免长循环停在 K>>

% ------- 工具箱路径(请确保正确) -------
addpath(genpath('D:\tools\Matlab\CoSMoMVPA-master'));   % 小规模复现修改：原为相对路径 'CoSMoMVPA-v.1.1.0/', 改为本机 CoSMoMVPA 绝对路径
% addpath(genpath('libsvm-3.25/'));   % 小规模复现修改：注释掉。本机只有 libsvm 3.17, CoSMo 要求 >=3.18, 本测试改用 LDA, 不需要 libsvm

% =========================================================================
% 0. 总开关 —— 只改这里
% =========================================================================
modal       = 'view';        % 'view' 或 'voice'
dry_run     = false;         % true=只核对文件/维度并汇总, 不跑 searchlight   % 小规模复现修改：true → false, 直接跑小规模 searchlight(维度核对信息仍会先打印)
time_window = '0-6';         % 小规模复现修改：'0-9' → '0-6'。本项目单试次只有 0-1.5/0-3/0-6/1.5-3/3-4.5/3-6/4.5-6 七个窗, 没有 0-9
radius_list = 3;             % cosmo 半径单位=体素; 数据2mm → radius=5 即半径10mm 的球   % 小规模复现修改：[3, 5] → 3, 只测一个半径
n_proc      = 1;             % searchlight 并行 worker 数(=1 则不并行); 只提速, 不改变结果   % 小规模复现修改：4 → 1, 小规模测试不启动 parpool

% ---- 4 状态 → 条件名(beta_info.mat 里的 condition 名; 屏蔽三合一) ----
% 顺序 = target 1..4: 1全 2自然 3屏蔽 4闭
% 小规模复现修改：原 4 状态(focus/nature/shield三合一/close)改为本项目 3 条件, 顺序与一级模型 NV, PS, CR 一致
state_labels = {'NV(natural_viewing)', 'PS(perspective_substitution)', 'CR(cognitive_reappraisal)'};
cond_names   = cell(3,1);
cond_names{1} = 'natural_viewing';
cond_names{2} = 'perspective_substitution';
cond_names{3} = 'cognitive_reappraisal';
nStates       = numel(state_labels);

ts = @() datestr(now, 'yyyy-mm-dd HH:MM:SS');
fprintf('[%s] ====== searchlight LOSO [状态版/单试次] 开始: modal=%s, dry_run=%d ======\n', ts(), modal, dry_run);
fprintf('[%s] 4 状态: %s\n', ts(), strjoin(state_labels, ' | '));

% =========================================================================
% 1. 按模态设置路径 / 被试 / mask
% =========================================================================
cfg = struct();
switch lower(modal)
    case 'view'
        cfg.base_path = 'F:\HeLAB\EEEEEEEmotional\617_single trial model';   % 单试次 beta 根目录   % 小规模复现修改：改为本项目单试次模型根目录
        cfg.mask_fn   = 'F:\HeLAB\EEEEEEEmotional\617_single trial model\S1\0-6\mask.nii';  % ⚠️须与 beta 同 MNI 网格, 不符请换   % 小规模复现修改：暂用 S1 一级模型 mask(91x109x91, 2mm, 与 beta 同网格); 正式分析建议换组水平/灰质 mask
        cfg.sub_list  = [1, 2, 3];    % 已做 single_trial 的视觉被试(有新增自行加)   % 小规模复现修改：只用 S1-S3 做 LOSO 测试(避开 covariate 数异常的 S12)
    case 'voice'
        cfg.base_path = 'H:\deinfo\deinfo_voice\single_trial';
        cfg.mask_fn   = 'H:\deinfo_voice\mask\GM\mean_gm_mask_2mm.nii';  % ⚠️同上, 须与 beta 同网格
        cfg.sub_list  = [];   % ⚠️填入做了 single_trial 的语音被试编号
    otherwise
        error('modal 必须是 ''view'' 或 ''voice'', 现在是: %s', modal);
end

if isempty(cfg.sub_list)
    error('modal=%s 的 sub_list 为空, 请先填入做了 single_trial 的被试编号。', modal);
end

save_path = 'C:\Users\ASUS\Documents\Codex\2026-11-18\new-chat-2\outputs\results_searchlight_LOSO_3conds_smalltest';   % 小规模复现修改：测试结果不写入 F 盘数据目录, 改写到 Codex outputs
if ~dry_run && ~exist(save_path, 'dir'); mkdir(save_path); end

% =========================================================================
% 2. 组装跨被试单试次 Dataset (chunk = 被试编号 → LOSO)
% =========================================================================
fprintf('[%s] --- 加载所有被试单试次 beta 并拼装 ---\n', ts());
ds_all       = [];
loaded_subs  = [];
skipped_subs = {};
trial_counts = [];   % 每被试样本数, 供核对
for s = cfg.sub_list
    subj_dir = fullfile(cfg.base_path, ['S' num2str(s)], time_window);   % 小规模复现修改：本项目被试文件夹名为 S1, S2..., 加 'S' 前缀
    ds_sub   = helper_load_singletrial_sub(subj_dir, cfg.mask_fn, cond_names, s);
    if isempty(ds_sub)
        fprintf('  [跳过] 被试 %d —— 无 beta_info.mat 或无匹配 trial (%s)\n', s, subj_dir);
        skipped_subs{end+1} = num2str(s); %#ok<AGROW>
        continue;
    end
    trial_counts(end+1) = numel(ds_sub.sa.targets); %#ok<AGROW>
    fprintf('  [加载] 被试 %d —— %d 个 trial (各状态: %s)\n', s, trial_counts(end), ...
            mat2str(arrayfun(@(k) sum(ds_sub.sa.targets==k), 1:nStates)));
    if isempty(ds_all); ds_all = ds_sub; else; ds_all = cosmo_stack({ds_all, ds_sub}); end
    loaded_subs(end+1) = s; %#ok<AGROW>
end

if isempty(ds_all)
    fprintf('[%s] !! 没有任何被试被成功加载, 请检查 single_trial 路径/beta_info/mask 网格。终止。\n', ts());
    return;
end
ds_all = cosmo_remove_useless_data(ds_all);

% ---- 维度核对 ----
n_subj = numel(unique(ds_all.sa.chunks));
fprintf('\n[%s] === 维度核对 ===\n', ts());
fprintf('  模态          : %s\n', modal);
fprintf('  时间窗        : %s\n', time_window);
fprintf('  加载被试数    : %d  (%s)\n', n_subj, num2str(loaded_subs));
fprintf('  跳过被试数    : %d  (%s)\n', numel(skipped_subs), strjoin(skipped_subs, ', '));
fprintf('  状态数        : %d  (%s)\n', nStates, strjoin(state_labels, ', '));
fprintf('  总样本数(行)  : %d  (每被试 trial 数不等, 屏蔽约为其它状态3倍)\n', size(ds_all.samples,1));
fprintf('  体素数(列)    : %d\n', size(ds_all.samples,2));
fprintf('  各状态样本数  : %s\n', mat2str(arrayfun(@(k) sum(ds_all.sa.targets==k), 1:nStates)));
fprintf('  targets 取值  : %s\n', num2str(unique(ds_all.sa.targets)'));

% =========================================================================
% 3. 批量跑不同 radius 的 Searchlight —— dry_run 时到此为止
% =========================================================================
measure                    = @cosmo_crossvalidation_measure;
measure_args               = struct();
measure_args.classifier    = @cosmo_classify_lda;               % libsvm 远快于 fitcecoc   % 小规模复现修改：libsvm → LDA(本机 libsvm 3.17 版本过低, CoSMo 会报 out of date)
measure_args.partitions    = cosmo_nfold_partitioner(ds_all);   % LOSO 交叉验证(chunk=被试)
measure_args.normalization = 'zscore';                          % 每折内 z 标准化(无泄漏)
measure_args.center_ids    = 1:2000:size(ds_all.samples, 2);   % 小规模复现修改：新增。只计算每 2000 个体素取 1 个的中心点, 正式运行时删除此行

if dry_run
    fprintf('\n[%s] (dry_run) 空跑结束: 未跑 searchlight, 未写结果。核对无误后把 dry_run 改为 false 再跑。\n', ts());
    return;
end

if n_proc > 1
    pool = gcp('nocreate');
    if isempty(pool)
        fprintf('[%s] 启动 parpool(%d) 用于 searchlight 并行(只提速, 不改变结果)...\n', ts(), n_proc);
        parpool(n_proc);
    end
    measure_args.nproc = n_proc;
end

fprintf('\n[%s] === 开始 %d 个 radius: %s ===\n', ts(), numel(radius_list), num2str(radius_list));
for r_idx = 1:numel(radius_list)
    radius = radius_list(r_idx);
    fprintf('[%s] --- radius = %g (第 %d/%d) ---\n', ts(), radius, r_idx, numel(radius_list));

    nbrhood = cosmo_spherical_neighborhood(ds_all, 'radius', radius);
    fprintf('[%s] Searchlight 运行中(单试次样本多, 可能很久)...\n', ts());
    acc_map = cosmo_searchlight(ds_all, nbrhood, measure, measure_args);

    save_filename = fullfile(save_path, sprintf('Group_LOSO_Searchlight_Acc_3conds_smalltest_%s_%s_r%g.nii', modal, time_window, radius));   % 小规模复现修改：文件名 4states_singletrial → 3conds_smalltest, 并加入时间窗
    cosmo_map2fmri(acc_map, save_filename);
    fprintf('[%s] radius=%g 完成 → %s\n', ts(), radius, save_filename);
end

fprintf('[%s] ====== 全部 radius 跑完 [状态版/单试次] (modal=%s) ======\n', ts(), modal);


% =========================================================================
% 单试次读取辅助函数 (解析 beta_info.mat) —— 搬自 Mvpa_1_1_basic_singletrial_Aug04.m
% Aug18: 逐 trial 单独调 cosmo_fmri_dataset 改为按被试一次性批量读取(cellstr), 顺序不变, 只提速。
% =========================================================================
function ds_sub = helper_load_singletrial_sub(subj_dir, mask_fn, cond_names, sub_id)
    ds_sub = [];
    beta_info_file = fullfile(subj_dir, 'beta_info.mat');
    if ~exist(beta_info_file, 'file'); return; end

    tmp = load(beta_info_file); f = fieldnames(tmp); raw = tmp.(f{1});
    if iscell(raw), tab = cell2table(raw(2:end,:), 'VariableNames', strtrim(raw(1,:)));
    elseif istable(raw), tab = raw; elseif isstruct(raw), tab = struct2table(raw); end

    tab.Properties.VariableNames = strtrim(tab.Properties.VariableNames);
    if isstring(tab.condition) || iscategorical(tab.condition), tab.condition = cellstr(tab.condition); end
    tab.condition = strtrim(strrep(tab.condition, '''', ''));

    if ismember('beta_number', tab.Properties.VariableNames)
        if iscell(tab.beta_number), tab.beta_number = cell2mat(tab.beta_number); end
        tab.beta_name = arrayfun(@(x) sprintf('beta_%04d.nii', x), tab.beta_number, 'UniformOutput', false);
    end

    % ---- 先收集该被试所有匹配 trial 的路径(顺序=原逐条件/逐 trial 顺序), 再一次性批量读取 ----
    all_paths   = {};
    targets_sub = [];
    for c = 1:length(cond_names)
        c_names = cond_names{c};
        if ischar(c_names) || isstring(c_names), c_names = {char(c_names)}; end

        idx = find(ismember(tab.condition, c_names));
        for tr = 1:length(idx)
            b_name = tab.beta_name{idx(tr)};
            b_path = fullfile(subj_dir, strrep(b_name, '''', ''));
            if exist(b_path, 'file')
                all_paths{end+1}     = b_path; %#ok<AGROW>
                targets_sub(end+1,1) = c;      %#ok<AGROW>
            end
        end
    end
    if isempty(all_paths), return; end

    % 小规模复现修改：本机 CoSMoMVPA-master 的 cosmo_fmri_dataset 不接受 cellstr(报 Could not find image format for <cell> input),
    %   改回逐 trial 读取再 cosmo_stack, 顺序与 all_paths 相同, 结果不变
    ds_sub      = cosmo_stack(cellfun(@(p) cosmo_fmri_dataset(p, 'mask', mask_fn), all_paths, 'UniformOutput', false));
    samples_sub = ds_sub.samples;

    valid_mask = ~any(isnan(samples_sub), 1) & ~any(samples_sub == 0, 1);
    if any(valid_mask); samples_sub(:, valid_mask) = zscore(samples_sub(:, valid_mask), 0, 2); end

    ds_sub.samples    = samples_sub;
    ds_sub.sa.targets = targets_sub; ds_sub.sa.chunks = repmat(sub_id, length(targets_sub), 1);
end
