%% ============================================================
%  make_wave_ati.m
%  海面波（正弦波）を表す Bellhop Altimetry (.ati) ファイルを生成
%  ・単一スナップショット、または時間位相をずらした複数枚を出力
%  ・生成した .ati を使うには .env のトップオプション4桁目に '*' を追加
%    例）'SVWT'  →  'SVWT*'
% ============================================================

%% ==== パラメータ（実環境に合わせて変更）====
base_name = 'JamTank';   % .env / .arr のベース名（'JamTank.env'なら'JamTank'）
H         = 0.5;         % 波高 peak-to-peak [m]（振幅は H/2）
lambda    = 5.0;         % 波長 [m]
R_max_m   = 120;         % 海面を定義する最大距離 [m]（受信機の最遠より少し広く）
pts_per_wl = 20;         % 1波長あたりの折れ線点数（形状の滑らかさ）
n_snap    = 8;           % 時間スナップショット数（波の位相を等分割）
interp    = 'C';         % 補間: 'L'=線形折れ線, 'C'=曲線(curvilinear/推奨)

% 深度の符号規約：Bellhop の .ati は「海面の深度[m]」
%   平水面 z=0 を基準に、山は上（負）、谷は下（正）
%   z(r) = -(H/2) * sin(2*pi*r/lambda + phase)
%   （符号は環境の座標系に合わせて確認。通常 下向き正）

%% ==== 距離グリッド ====
R_max_km = R_max_m / 1000;                 % .ati の距離単位は km
dr_m     = lambda / pts_per_wl;            % 点間隔 [m]
r_m      = (0 : dr_m : R_max_m).';         % 距離ベクトル [m]
r_km     = r_m / 1000;
N        = numel(r_m);

%% ==== スナップショットごとに .ati を書き出し ====
for k = 1:n_snap
    phase = 2*pi*(k-1)/n_snap;             % 位相を等分割（波の進行を表現）
    z_m   = -(H/2) * sin(2*pi*r_m/lambda + phase);   % 海面深度 [m]

    if n_snap == 1
        fname = sprintf('%s.ati', base_name);
    else
        fname = sprintf('%s_wave%02d.ati', base_name, k);  % 時系列は連番
    end

    fid = fopen(fname, 'w');
    fprintf(fid, '''%s''\n', interp);      % 補間方式
    fprintf(fid, '%d\n', N);               % 海面点数
    for i = 1:N
        fprintf(fid, '%.6f %.6f\n', r_km(i), z_m(i));  % 距離(km)  深度(m)
    end
    fclose(fid);

    fprintf('生成: %s  (波高%.2fm, 波長%.2fm, 位相%.0f°, %d点)\n', ...
            fname, H, lambda, rad2deg(phase), N);
end

fprintf('\n次のステップ:\n');
fprintf(' 1) %s.env のトップオプション文字列4桁目に ''*'' を追加\n', base_name);
fprintf('    例)  ''SVWT''  ->  ''SVWT*''\n');
if n_snap == 1
    fprintf(' 2) bellhop(''%s'') を実行 → %s.arr が波込みで生成される\n', base_name, base_name);
else
    fprintf(' 2) 各スナップショットで .ati をリネームして Bellhop 実行:\n');
    fprintf('    for k=1:%d\n', n_snap);
    fprintf('        copyfile(sprintf(''%s_wave%%02d.ati'',k), ''%s.ati'');\n', base_name, base_name);
    fprintf('        bellhop(''%s'');\n', base_name);
    fprintf('        copyfile(''%s.arr'', sprintf(''%s_wave%%02d.arr'',k));\n', base_name, base_name);
    fprintf('    end\n');
    fprintf('    → 各 %s_waveXX.arr を SCM_PTR_simulator で順に解析し、\n', base_name);
    fprintf('      波の位相による OSNR/BER 変動を追える\n');
end
