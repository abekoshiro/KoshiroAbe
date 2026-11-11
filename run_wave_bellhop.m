%% ============================================================
%  run_wave_bellhop.m
%  make_wave_ati.m が生成した 8枚のスナップショット .ati を順に
%  Bellhop で実行し、波位相ごとの ARR（_waveXX.arr）を生成する。
%
%  前提:
%   ・JamTank_Arr_wave01.ati … _wave08.ati が存在すること
%   ・JamTank_Arr.env のトップオプション5桁目に '*'（例 'CVW *'）
%   ・SSP最上点が波の山をカバー（z <= -H/2）していること
%   ・bellhop が MATLAB パス上にあること（Acoustics Toolbox）
%
%  出力: JamTank_Arr_wave01.arr … _wave08.arr
%        → Wave_TimeVarying_simulator_ver1.m が読み込む
% ============================================================
base   = 'JamTank_Arr';
n_snap = 8;

for k = 1:n_snap
    ati_k = sprintf('%s_wave%02d.ati', base, k);
    if ~isfile(ati_k)
        error('見つかりません: %s （先に make_wave_ati.m を実行）', ati_k);
    end
    % スナップショット k を本体名にコピー（Bellhopは base.ati しか読まない）
    copyfile(ati_k, sprintf('%s.ati', base));

    % Bellhop 実行
    bellhop(base);

    % 生成された base.arr を波位相連番に退避
    arr_out = sprintf('%s_wave%02d.arr', base, k);
    copyfile(sprintf('%s.arr', base), arr_out);

    fprintf('スナップショット %d/%d 完了 → %s\n', k, n_snap, arr_out);
end

fprintf('\n全%d枚の ARR 生成完了。\n', n_snap);
fprintf('次: Wave_TimeVarying_simulator_ver1.m を実行（f_wave を実験条件に設定してから）\n');
