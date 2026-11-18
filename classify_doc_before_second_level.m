win = '0-9';      % 时间窗口
con = 'con_0005'; % 文件编号

src = 'F:\HeLAB\EEEEEEEmotional\first_level_new';
dst = fullfile('F:\HeLAB\EEEEEEEmotional\temporal', [win '_' con]);
mkdir(dst);
subs = dir(fullfile(src, 'S*'));
for i = 1:numel(subs)
    copyfile(fullfile(src, subs(i).name, win, [con '.nii']), ...
             fullfile(dst, [subs(i).name '_' con '.nii']));
end
