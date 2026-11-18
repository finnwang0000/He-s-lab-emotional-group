clear;
clc;

%% 0. 初始化
cd(fileparts(mfilename('fullpath')));
restoredefaultpath;
rehash toolboxcache;
addpath('D:\tools\Matlab\spm12');
spm('defaults','FMRI');
spm_jobman('initcfg');

%% 1. 参数配置

% 填写需要运行的被试
subjects_to_run = {'S17'};

parts_to_run = {'run1','run2','run3'};

preproc_base_dir = 'F:\HeLAB\EEEEEEEmotional\preprocessing';
timing_base_dir  = 'F:\HeLAB\EEEEEEEmotional\onset_in_different_format\onsets_new';
output_base_dir  = 'F:\HeLAB\EEEEEEEmotional\temporal';

conditions = {
    'natural_viewing'
    'perspective_substitution'
    'cognitive_reappraisal'
};

% 原始 onset 已对应删除前 10 个 scan 后的时间轴，不再减 15 秒。
% 仅在建模时添加下列偏移，不修改原始 onset MAT。
% 六个 1.5 秒连续时间窗，加五个指定时间窗；共 11 个时间窗，分别独立建模。
% 每行依次为：时间窗标签、相对原始 onset 的偏移（秒）、duration（秒）。
durations_cfg = {
    '0-1.5', 0, 1.5;
    '1.5-3', 1.5, 1.5;
    '3-4.5', 3, 1.5;
    '4.5-6', 4.5, 1.5;
    '6-7.5', 6, 1.5;
    '7.5-9', 7.5, 1.5;
    '0-6', 0, 6;
    '0-9', 0, 9;
    '3-6', 3, 3;
    '6-9', 6, 3;
    '0-3', 0, 3
};

TR = 1.5;
FMRI_T = 72;
FMRI_T0 = 17;

% 与其他普通被试保持一致
HRF_DERIVS = [0 0];

SKIP_COMPLETED = true;

%% 2. 主程序

for i = 1:numel(subjects_to_run)

    subject_number = subjects_to_run{i};
    subject_numeric = str2double(regexprep(subject_number,'^[Ss]',''));

    if isnan(subject_numeric)
        error('无法识别被试编号：%s',subject_number);
    end

    timing_subject = sprintf('S%d',subject_numeric);

    fprintf('\n====================================================\n');
    fprintf('开始处理被试：%s（onset文件夹：%s）\n', ...
        subject_number,timing_subject);
    fprintf('====================================================\n');

    try
        for d = 1:size(durations_cfg,1)

            dur_label   = durations_cfg{d,1};
            onset_shift = durations_cfg{d,2};
            event_dur   = durations_cfg{d,3};

            fprintf('\n--- 时间窗：%s｜onset + %.1fs｜duration %.1fs ---\n', ...
                dur_label,onset_shift,event_dur);

            output_dir = fullfile(output_base_dir, ...
                subject_number,dur_label);

            if ~exist(output_dir,'dir')
                mkdir(output_dir);
            end

            expected_last_con = fullfile(output_dir, ...
                sprintf('con_%04d.nii',7));

            existing_spm = fullfile(output_dir,'SPM.mat');

            if SKIP_COMPLETED && exist(expected_last_con,'file')
                fprintf('  [跳过] 已完成：%s\n',output_dir);
                continue;
            end

            if exist(existing_spm,'file') && ...
                    ~exist(expected_last_con,'file')
                error(['发现未完成或旧版本结果：%s。' ...
                    '请先删除该时间窗文件夹后重新运行。'], ...
                    output_dir);
            end

            clear matlabbatch

            %% 3. Model specification

            matlabbatch{1}.spm.stats.fmri_spec.dir = {output_dir};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = FMRI_T;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = FMRI_T0;

            nuisance_dir = fullfile(output_dir,'nuisance');

            if ~exist(nuisance_dir,'dir')
                mkdir(nuisance_dir);
            end

            %% 4. 遍历run

            for sess = 1:numel(parts_to_run)

                part_name = parts_to_run{sess};

                subject_dir = fullfile(preproc_base_dir, ...
                    subject_number);

                % 兼容硬盘中的四种run命名：
                % run1、run1_new、run1_art、run1_art_new。
                run_variant_dirs = {
                    fullfile(subject_dir,part_name)
                    fullfile(subject_dir,[part_name '_new'])
                    fullfile(subject_dir,[part_name '_art'])
                    fullfile(subject_dir,[part_name '_art_new'])
                };

                fprintf('\n  正在载入%s...\n',part_name);

                %% 4.1 自动查找功能像

                candidate_nii_dirs = {};

                for rd = 1:numel(run_variant_dirs)
                    current_run_root = run_variant_dirs{rd};
                    candidate_nii_dirs = [candidate_nii_dirs; {
                        fullfile(current_run_root,'sw')
                        fullfile(current_run_root,'swuaf')
                        current_run_root
                    }]; %#ok<AGROW>
                end

                % 兼容S1/S3/S5/S11的sw/run1布局。
                candidate_nii_dirs{end+1,1} = ...
                    fullfile(subject_dir,'sw',part_name);
                candidate_nii_dirs = unique( ...
                    candidate_nii_dirs,'stable');

                nii_struct = [];
                nii_dir = '';

                % 必须且只能找到一套swuaf*.nii，避免误选旧数据。
                for nd = 1:numel(candidate_nii_dirs)

                    current_dir = candidate_nii_dirs{nd};

                    if ~exist(current_dir,'dir')
                        continue;
                    end

                    temp_nii = dir(fullfile(current_dir,'swuaf*.nii'));

                    if ~isempty(temp_nii)
                        hidden_mask = startsWith({temp_nii.name},'._');
                        temp_nii = temp_nii(~hidden_mask);
                    end

                    if ~isempty(temp_nii)
                        if ~isempty(nii_struct)
                            error(['%s找到多套swuaf NIfTI：%s 和 %s。' ...
                                '请先人工确认应使用哪套。'], ...
                                part_name,nii_dir,current_dir);
                        end
                        nii_struct = temp_nii;
                        nii_dir = current_dir;
                    end
                end

                if isempty(nii_struct)
                    error('%s未找到swuaf NIfTI。已搜索：%s', ...
                        part_name,strjoin(candidate_nii_dirs,' | '));
                end

                [~,order] = sort({nii_struct.name});
                nii_struct = nii_struct(order);

                scans = cell(numel(nii_struct),1);

                for n = 1:numel(nii_struct)

                    nii_file = fullfile( ...
                        nii_struct(n).folder,nii_struct(n).name);

                    scans{n,1} = sprintf('%s,1',nii_file);
                end

                num_scans = numel(scans);
                run_duration = num_scans * TR;

                matlabbatch{1}.spm.stats.fmri_spec ...
                    .sess(sess).scans = scans;

                fprintf('  功能像目录：%s\n',nii_dir);
                fprintf('  %s：找到%d个swuaf文件。\n', ...
                    part_name,num_scans);

                %% 4.2 读取onset

                % 兼容 run1.mat 和 run1/onset.mat；多个候选时禁止猜选。
                timing_candidates = {
                    fullfile(timing_base_dir,timing_subject,[part_name '.mat'])
                    fullfile(timing_base_dir,timing_subject,part_name,'onset.mat')
                };
                timing_exists = cellfun(@(p) exist(p,'file') == 2,timing_candidates);
                if sum(timing_exists) ~= 1
                    error('onset文件缺失或存在多个候选：%s / %s', ...
                        timing_candidates{1},timing_candidates{2});
                end
                timing_file = timing_candidates{find(timing_exists,1)};

                timing_data = load(timing_file);

                if ~isfield(timing_data,'names') || ...
                        ~isfield(timing_data,'onsets')
                    error('onset文件缺少names或onsets变量：%s', ...
                        timing_file);
                end

                if ~iscell(timing_data.onsets)
                    error('onset文件中的onsets不是cell数组：%s', ...
                        timing_file);
                end

                if iscell(timing_data.names)
                    loaded_names = cellfun( ...
                        @(x) strtrim(char(x)), ...
                        timing_data.names, ...
                        'UniformOutput',false);
                elseif isstring(timing_data.names)
                    loaded_names = cellstr(timing_data.names);
                elseif ischar(timing_data.names)
                    loaded_names = cellstr(timing_data.names);
                else
                    error('无法识别names变量格式：%s',timing_file);
                end

                loaded_names = loaded_names(:)';

                fprintf('  %s onset文件：%s\n', ...
                    part_name,timing_file);

                fprintf('  %s onset条件：%s\n', ...
                    part_name,strjoin(loaded_names,', '));

                %% 4.3 条件建模

                for c = 1:numel(conditions)

                    cond_name = conditions{c};
                    match_idx = find(strcmp(loaded_names,cond_name),1);

                    if isempty(match_idx)
                        error('%s缺少condition：%s', ...
                            part_name,cond_name);
                    end

                    base_onset = double( ...
                        timing_data.onsets{match_idx});

                    base_onset = base_onset(isfinite(base_onset));
                    base_onset = base_onset(:)';

                    if isempty(base_onset)
                        error('%s的%s onset为空。', ...
                            part_name,cond_name);
                    end

                    onset_data = base_onset + onset_shift;

                    if any(onset_data < 0)
                        error('%s的%s存在负数onset。', ...
                            part_name,cond_name);
                    end

                    if any(onset_data + event_dur > ...
                            run_duration + 1e-6)
                        error(['%s的%s超出扫描时长。' ...
                            '扫描时长=%.3fs，最大事件结束时间=%.3fs。'], ...
                            part_name,cond_name,run_duration, ...
                            max(onset_data + event_dur));
                    end

                    duration_data = repmat( ...
                        event_dur,size(onset_data));

                    matlabbatch{1}.spm.stats.fmri_spec ...
                        .sess(sess).cond(c).name = cond_name;

                    matlabbatch{1}.spm.stats.fmri_spec ...
                        .sess(sess).cond(c).onset = onset_data;

                    matlabbatch{1}.spm.stats.fmri_spec ...
                        .sess(sess).cond(c).duration = duration_data;

                    matlabbatch{1}.spm.stats.fmri_spec ...
                        .sess(sess).cond(c).tmod = 0;

                    matlabbatch{1}.spm.stats.fmri_spec ...
                        .sess(sess).cond(c).orth = 1;

                    matlabbatch{1}.spm.stats.fmri_spec ...
                        .sess(sess).cond(c).pmod = ...
                        struct('name',{},'param',{},'poly',{});

                    fprintf(['    %s：%d个事件｜' ...
                        '首个onset=%.3fs｜末个onset=%.3fs\n'], ...
                        cond_name,numel(onset_data), ...
                        min(onset_data),max(onset_data));
                end

                %% 4.4 读取RP文件

                rp_search_dirs = {nii_dir};

                for rd = 1:numel(run_variant_dirs)
                    current_run_root = run_variant_dirs{rd};
                    rp_search_dirs = [rp_search_dirs; {
                        fullfile(current_run_root,'rp')
                        current_run_root
                    }]; %#ok<AGROW>
                end

                % 兼容S1/S3/S5/S11的rp/run1布局。
                rp_search_dirs{end+1,1} = ...
                    fullfile(subject_dir,'rp',part_name);
                rp_search_dirs = unique(rp_search_dirs,'stable');

                rp_paths = {};

                for sd = 1:numel(rp_search_dirs)

                    current_dir = rp_search_dirs{sd};

                    if ~exist(current_dir,'dir')
                        continue;
                    end

                    temp_rp = dir(fullfile(current_dir,'rp_*.txt'));

                    for k = 1:numel(temp_rp)

                        if ~startsWith(temp_rp(k).name,'._')
                            rp_paths{end+1} = fullfile( ...
                                temp_rp(k).folder, ...
                                temp_rp(k).name); %#ok<SAGROW>
                        end
                    end
                end

                rp_paths = unique(rp_paths,'stable');

                if isempty(rp_paths)
                    error('%s未找到rp头动文件。已搜索：%s', ...
                        part_name,strjoin(rp_search_dirs,' | '));
                end

                if numel(rp_paths) > 1
                    error('%s找到多个rp文件，禁止自动选择：%s', ...
                        part_name,strjoin(rp_paths,' | '));
                end

                rp_file = rp_paths{1};
                rp_data = double(load(rp_file));

                if isempty(rp_data) || ~isnumeric(rp_data)
                    error('rp文件不是有效数字矩阵：%s',rp_file);
                end

                if size(rp_data,1) ~= num_scans
                    error(['%s功能像数量与rp行数不一致：' ...
                        '功能像=%d，rp行数=%d'], ...
                        part_name,num_scans,size(rp_data,1));
                end

                if any(~isfinite(rp_data(:)))
                    error('rp文件中存在NaN或Inf：%s',rp_file);
                end

                %% 4.5 自动查找ART回归量

                art_search_dirs = unique( ...
                    [rp_search_dirs;candidate_nii_dirs], ...
                    'stable');

                art_paths = {};

                for sd = 1:numel(art_search_dirs)

                    current_dir = art_search_dirs{sd};

                    if ~exist(current_dir,'dir')
                        continue;
                    end

                    temp_art_mat = dir(fullfile(current_dir,'art*.mat'));
                    temp_art_txt = dir(fullfile(current_dir,'art*.txt'));
                    temp_art = [temp_art_mat;temp_art_txt];

                    for k = 1:numel(temp_art)

                        if ~startsWith(temp_art(k).name,'._')
                            art_paths{end+1} = fullfile( ...
                                temp_art(k).folder, ...
                                temp_art(k).name); %#ok<SAGROW>
                        end
                    end
                end

                art_paths = unique(art_paths,'stable');
                art_data = [];
                art_file = '';

                if ~isempty(art_paths)

                    % 给ART文件打分，优先选择完整回归量文件
                    art_scores = zeros(1,numel(art_paths));

                    for a = 1:numel(art_paths)

                        art_name_lower = lower(art_paths{a});

                        if contains(art_name_lower,'and_movement')
                            art_scores(a) = art_scores(a) + 100;
                        end

                        if contains(art_name_lower,'regression_outliers')
                            art_scores(a) = art_scores(a) + 50;
                        end

                        if endsWith(art_name_lower,'.mat')
                            art_scores(a) = art_scores(a) + 10;
                        end
                    end

                    [~,best_art_idx] = max(art_scores);
                    art_file = art_paths{best_art_idx};

                    fprintf('  找到ART文件：%s\n',art_file);

                    [~,~,art_ext] = fileparts(art_file);

                    if strcmpi(art_ext,'.mat')

                        art_loaded = load(art_file);

                        if isfield(art_loaded,'R') && ...
                                isnumeric(art_loaded.R)

                            art_data = double(art_loaded.R);

                        else
                            field_names = fieldnames(art_loaded);
                            candidate_names = {};

                            for f = 1:numel(field_names)

                                current_value = ...
                                    art_loaded.(field_names{f});

                                if isnumeric(current_value) && ...
                                        ismatrix(current_value) && ...
                                        size(current_value,1) == num_scans

                                    candidate_names{end+1} = ...
                                        field_names{f}; %#ok<SAGROW>
                                end
                            end

                            if isempty(candidate_names)
                                error(['ART文件中没有找到行数等于' ...
                                    'scan数量的数字矩阵：%s'], ...
                                    art_file);
                            end

                            if numel(candidate_names) > 1
                                error(['ART文件中存在多个可能的' ...
                                    '回归量矩阵，请检查：%s'], ...
                                    art_file);
                            end

                            art_data = double( ...
                                art_loaded.(candidate_names{1}));
                        end

                    else
                        art_data = double(load(art_file));
                    end

                    if size(art_data,1) ~= num_scans
                        error(['ART行数与scan数量不一致：' ...
                            'scan=%d，ART行数=%d\n%s'], ...
                            num_scans,size(art_data,1),art_file);
                    end

                    if any(~isfinite(art_data(:)))
                        error('ART文件中存在NaN或Inf：%s',art_file);
                    end
                end

                %% 4.6 生成最终协变量

                if isempty(art_data)

                    nuisance_R = rp_data;
                    nuisance_source = '普通RP';

                else
                    % 判断ART是否已经包含6列运动参数
                    art_has_motion = ...
                        contains(lower(art_file),'and_movement') || ...
                        contains(lower(art_file),'motion');

                    if ~art_has_motion && ...
                            size(art_data,2) >= 6 && ...
                            size(rp_data,2) >= 6

                        motion_difference = ...
                            art_data(:,1:6) - rp_data(:,1:6);

                        if max(abs(motion_difference(:))) < 1e-6
                            art_has_motion = true;
                        end
                    end

                    if art_has_motion
                        nuisance_R = art_data;
                        nuisance_source = ...
                            'ART完整回归量（已含运动参数）';
                    else
                        nuisance_R = [rp_data,art_data];
                        nuisance_source = ...
                            'RP + ART异常点回归量';
                    end
                end

                nuisance_file = fullfile(nuisance_dir, ...
                    sprintf('%s_%s_nuisance.txt', ...
                    subject_number,part_name));

                save(nuisance_file,'nuisance_R','-ascii','-tabs');

                matlabbatch{1}.spm.stats.fmri_spec ...
                    .sess(sess).multi_reg = {nuisance_file};

                matlabbatch{1}.spm.stats.fmri_spec ...
                    .sess(sess).multi = {''};

                matlabbatch{1}.spm.stats.fmri_spec ...
                    .sess(sess).regress = ...
                    struct('name',{},'val',{});

                matlabbatch{1}.spm.stats.fmri_spec ...
                    .sess(sess).hpf = 128;

                fprintf(['  [完成] %s：%d个scan，%d个条件，' ...
                    '%d列协变量（%s）。\n'], ...
                    part_name,num_scans,numel(conditions), ...
                    size(nuisance_R,2),nuisance_source);
            end

            %% 5. 模型设置

            matlabbatch{1}.spm.stats.fmri_spec.fact = ...
                struct('name',{},'levels',{});

            matlabbatch{1}.spm.stats.fmri_spec ...
                .bases.hrf.derivs = HRF_DERIVS;

            matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
            matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

            %% 6. 模型估计

            matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = ...
                cfg_dep( ...
                'fMRI model specification: SPM.mat File', ...
                substruct('.', 'val','{}',{1}, ...
                          '.', 'val','{}',{1}, ...
                          '.', 'val','{}',{1}), ...
                substruct('.','spmmat'));

            matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

            %% 7. Contrast：跨Session复制并缩放

            matlabbatch{3}.spm.stats.con.spmmat(1) = ...
                cfg_dep( ...
                'Model estimation: SPM.mat File', ...
                substruct('.', 'val','{}',{2}, ...
                          '.', 'val','{}',{1}, ...
                          '.', 'val','{}',{1}), ...
                substruct('.','spmmat'));

            for c = 1:numel(conditions)

                condition_weights = zeros(1,numel(conditions));
                condition_weights(c) = 1;

                matlabbatch{3}.spm.stats.con ...
                    .consess{c}.tcon.name = conditions{c};

                matlabbatch{3}.spm.stats.con ...
                    .consess{c}.tcon.weights = ...
                    condition_weights;

                matlabbatch{3}.spm.stats.con ...
                    .consess{c}.tcon.sessrep = 'replsc';
            end

            matlabbatch{3}.spm.stats.con.delete = 0;

            %% 8. 保存并运行

            batch_dir = fullfile(output_dir,'batch');

            if ~exist(batch_dir,'dir')
                mkdir(batch_dir);
            end

            batch_filename = fullfile(batch_dir, ...
                sprintf('1st_batch_sub-%s_dur-%s.mat', ...
                subject_number,dur_label));

            save(batch_filename,'matlabbatch');

            fprintf('\n>> 开始运行SPM Job...\n');
            spm_jobman('run',matlabbatch);

            %% 8.1 合并 add.m 的四个条件差异 contrast（原权重：跨 run 求和）
            result_spm = fullfile(output_dir,'SPM.mat');
            loaded_model = load(result_spm,'SPM');
            SPM = loaded_model.SPM;
            assert(numel(SPM.xCon)==3,'预期原有三个 activation contrast');
            condition_columns = zeros(size(SPM.xX.X,2),3);
            for cc=1:3
                for rr=1:numel(SPM.Sess)
                    reg_name = sprintf('Sn(%d) %s*bf(1)',rr,conditions{cc});
                    col = find(strcmp(SPM.xX.name,reg_name));
                    assert(numel(col)==1,'条件设计列缺失或重复：%s',reg_name);
                    condition_columns(col,cc)=1;
                end
            end
            extra_names = {'Substitution > Reappraisal', ...
                'Reappraisal > Substitution','Substitution > Passive', ...
                'Reappraisal > Passive'};
            extra_weights = [0 1 -1;0 -1 1;-1 1 0;-1 0 1];
            for cc=1:4
                SPM.xCon(end+1)=spm_FcUtil('Set',extra_names{cc},'T','c', ...
                    condition_columns*extra_weights(cc,:)',SPM.xX.xKXs);
            end
            save(result_spm,'SPM');
            SPM=spm_contrasts(SPM,4:7);
            assert(numel(SPM.xCon)==7);
            fprintf('>> 已生成全部七个 contrast：%s %s\n',subject_number,dur_label);


            if exist(expected_last_con,'file')
                fprintf('>> 完成被试%s，时间窗%s。\n', ...
                    subject_number,dur_label);
            else
                warning(['SPM运行结束，但没有检测到%s。' ...
                    '请检查输出目录。'],expected_last_con);
            end
        end

    catch ME

        fprintf(2,'\n!!!! 被试%s处理失败：%s\n', ...
            subject_number,ME.message);

        if ~isempty(ME.stack)
            fprintf(2,'出错脚本：%s\n',ME.stack(1).name);
            fprintf(2,'出错行号：%d\n',ME.stack(1).line);
        end

        fprintf(2,'该被试后续时间窗已停止，继续下一个被试。\n');
        continue;
    end
end

fprintf('\n====================================================\n');
fprintf('所有任务处理完毕。\n');
fprintf('====================================================\n');
