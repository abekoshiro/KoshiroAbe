%% ============================================================
%  arrival_angle_delay_map.m
%  到来波を「到来角 × 遅延 × パワー」で可視化する
%    横軸: 遅延 [s]（最速到来を0とした相対遅延）
%    縦軸: 受信到来角 [deg]
%    色  : パワー [dB]（最大到来を0 dBとした相対値）
%  → 「何度の波が、何秒遅れて、どれだけのエネルギーで届くか」を一目で把握
%
%   ※ARRファイルに受信到来角が含まれている必要があります。
%    read_arrivals_asc のバージョン差に対応するため角度フィールド名を自動判定。
% ============================================================
[Arr, Pos] = read_arrivals_asc('JamTank_Arr.arr');

Nrr = size(Arr,1); Nrz = size(Arr,2); Nsd = size(Arr,3);

%% ==== 表示対象（1対1リンク）を選択 ====
idx_sd = 1;    % 音源番号
idx_rr = 1;    % 受信距離インデックス（1..Nrr）
idx_rz = 1;    % 受信深度インデックス（1..Nrz）

%% ==== 表示オプション ====
dyn_range = 40;      % 表示ダイナミックレンジ [dB]（-dyn_range dB 未満は非表示）
use_abs_delay = false;  % true=絶対到来時刻[s], false=最速到来を0とした相対遅延

%% ==== 到来データの取り出し ====
el = Arr(idx_rr, idx_rz, idx_sd);
A     = el.A(:);                     % 複素振幅
delay = real(el.delay(:));           % 遅延 [s]
ang   = get_rcvr_angle(el);          % 受信到来角 [deg]（自動判定）

% 有効な到来のみ
v = (delay>0) & isfinite(delay) & (abs(A)>0);
A = A(v); delay = delay(v);
if ~isempty(ang), ang = ang(v); end

if isempty(ang)
    error(['受信到来角フィールドが見つかりません。read_arrivals_asc が角度を', ...
           '読み込むバージョンか確認してください（RcvrDeclAngle 等）。']);
end

%% ==== パワー[dB]・遅延[s] の整形 ====
powerdB = 20*log10(abs(A)/max(abs(A)));   % 最大到来を0 dB
if ~use_abs_delay
    delay = delay - min(delay);           % 最速到来を0に
end

% ダイナミックレンジで足切り
keep = powerdB > -dyn_range;
A=A(keep); delay=delay(keep); ang=ang(keep); powerdB=powerdB(keep);

%% ==== 情報表示 ====
fprintf('=== 到来角-遅延マップ: 音源%d, 受信(距離idx%d, 深度idx%d) ===\n', idx_sd, idx_rr, idx_rz);
if isfield(Pos,'r')&&isfield(Pos.r,'r'), fprintf('  受信距離: %.1f m\n', Pos.r.r(idx_rr)); end
if isfield(Pos,'r')&&isfield(Pos.r,'z'), fprintf('  受信深度: %.2f m\n', Pos.r.z(idx_rz)); end
fprintf('  到来波数(表示): %d 本（%.0f dBダイナミックレンジ）\n', numel(A), dyn_range);
fprintf('  到来角範囲: %.1f 〜 %.1f deg\n', min(ang), max(ang));
fprintf('  遅延範囲  : %.3f 〜 %.3f ms\n\n', min(delay)*1e3, max(delay)*1e3);

%% ==== 描画：到来角 × 遅延 × パワー ====
figure('Name','Arrival angle-delay-power map','Position',[100 100 760 560]);
scatter(delay*1e3, ang, 60, powerdB, 'filled', 'MarkerEdgeColor',[0.2 0.2 0.2]);
grid on; box on;
xlabel('遅延 [ms]'); ylabel('受信到来角 [deg]');
cb = colorbar; ylabel(cb, 'パワー [dB]（最大到来=0 dB）');
colormap(jet); caxis([-dyn_range 0]);
if use_abs_delay
    ttl_x = '絶対到来時刻';
else
    ttl_x = '相対遅延（最速到来=0）';
end
title(sprintf('到来角-遅延-パワー  音源%d 受信(%d,%d)   横軸:%s', ...
      idx_sd, idx_rr, idx_rz, ttl_x));

%% ============================================================
function ang = get_rcvr_angle(el)
% ARR要素から受信到来角[deg]を取り出す（バージョン差に対応）
    cands = {'RcvrDeclAngle','RcvrAngle','rcvrAngle','RcvrDeclAngleDeg', ...
             'RecvrAngle','AngleRcvr','rd_angle','beta'};
    ang = [];
    for k=1:numel(cands)
        if isfield(el, cands{k})
            ang = el.(cands{k})(:);
            return;
        end
    end
    % 見つからない場合は空を返す（呼び出し側でエラー表示）
end
