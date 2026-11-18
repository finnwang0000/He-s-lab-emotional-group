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
% =========================================================================
clear; clc;
dbclear all;   % 清残留断点, 避免长循环停在 K>>

% ------- 工具箱路径(请确保正确) -------
addpath(genpath('CoSMoMVPA-v.1.1.0/'));
addpath(genpath('libsvm-3.25/'));

% =========================================================================
% 0. 总开关 —— 只改这里
% =========================================================================
modal       = 'view';        % 'view' 或 'voice'
dry_run     = true;          % true=只核对文件/维度并汇总, 不跑 searchlight
time_window = '0-6';
radius_list = [3, 5];        % cosmo 半径单位=体素; 数据2mm → radius=5 即半径10mm 的球
n_proc      = 4;             % searchlight 并行 worker 数(=1 则不并行); 只提速, 不改变结果

% ---- 4 状态 → 条件名(beta_info.mat 里的 condition 名; 屏蔽三合一) ----
% 顺序 = target 1..4: 1全 2自然 3屏蔽 4闭
state_labels = {'Focus(全)', 'Nature(自然)', 'Shielding(屏蔽)', 'Close(闭)'};
cond_names   = cell(4,1);
cond_names{1} = 'focus-clear';
cond_names{2} = 'nature-clear';
cond_names{3} = {'distancing-clear', 'distraction-clear', 'inhibition-clear'};  % 三合一 → target=3
cond_names{4} = 'close-clear';
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
        cfg.base_path = 'H:\deinfo\deinfo_view\single_trial';   % 单试次 beta 根目录
        cfg.mask_fn   = 'H:\deinfo_view\mask\GM\rmean_gm_mask.nii';  % ⚠️须与 beta 同 MNI 网格, 不符请换
        cfg.sub_list  = [4, 5, 6, 7, 9, 10, 11, 13, 14, 15];    % 已做 single_trial 的视觉被试(有新增自行加)
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

save_path = fullfile(cfg.base_path, 'results_searchlight_LOSO_4states_singletrial');
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
    subj_dir = fullfile(cfg.base_path, num2str(s), time_window);
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
measure_args.classifier    = @cosmo_classify_libsvm;            % libsvm 远快于 fitcecoc
measure_args.partitions    = cosmo_nfold_partitioner(ds_all);   % LOSO 交叉验证(chunk=被试)
measure_args.normalization = 'zscore';                          % 每折内 z 标准化(无泄漏)

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

    save_filename = fullfile(save_path, sprintf('Group_LOSO_Searchlight_Acc_4states_singletrial_%s_r%g.nii', modal, radius));
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

    ds_sub      = cosmo_fmri_dataset(all_paths, 'mask', mask_fn);
    samples_sub = ds_sub.samples;

    valid_mask = ~any(isnan(samples_sub), 1) & ~any(samples_sub == 0, 1);
    if any(valid_mask); samples_sub(:, valid_mask) = zscore(samples_sub(:, valid_mask), 0, 2); end

    ds_sub.samples    = samples_sub;
    ds_sub.sa.targets = targets_sub; ds_sub.sa.chunks = repmat(sub_id, length(targets_sub), 1);
end
