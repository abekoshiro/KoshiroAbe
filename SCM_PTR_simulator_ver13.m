%% ============================================================
%  SCM-PTR (Subarray Combining Passive Time Reversal)  ver13
%  ★ver11ベース：距離群（サブアレイ）ごとの星座図描画を追加
%    例）75m に8深度 + 100m に8深度 = 計16ch
%  通常PTR（全ch直接合成）と SCM-PTR（距離ごとサブアレイ→MRC合成）を比較
%  ※DFEなし。PTR/アレイ合成のみでOSNRを評価
% ============================================================
[Arr, Pos] = read_arrivals_asc('JamTank_Arr.arr');

%% ARRの次元を自動取得（ハードコードしない）
Nrr = size(Arr, 1);     % 受信距離の数（例：2 → 75m, 100m）
Nrz = size(Arr, 2);     % 受信深度の数（例：8 → 1.25〜4.75m）
Nsd = size(Arr, 3);     % 音源の数
NUM_RX = Nrr * Nrz;     % 全受信機チャネル数（例：2×8 = 16）

fprintf('ARR構成: 受信距離 %d点 × 受信深度 %d点 = %dch, 音源 %d個\n', ...
        Nrr, Nrz, NUM_RX, Nsd);
if isfield(Pos,'r') && isfield(Pos.r,'r')
    fprintf('  受信距離: '); fprintf('%.1f ', Pos.r.r); fprintf('m\n');
end
if isfield(Pos,'r') && isfield(Pos.r,'z')
    fprintf('  受信深度: '); fprintf('%.2f ', Pos.r.z); fprintf('m\n\n');
end

%% ★SCM設定：サブアレイの分割
% 受信距離ごとに1サブアレイとするのが物理的に自然（75m群 / 100m群）
% → NUM_SUB = Nrr, 各サブアレイ = その距離の全深度(Nrz ch)
NUM_SUB  = Nrr;             % サブアレイ数（= 受信距離の数）
SUB_SIZE = Nrz;            % 1サブアレイあたりのch数（= 深度数）
% ※距離が1つしかない旧ARRの場合は深度を等分する設定に切り替える
if NUM_SUB < 2
    NUM_SUB  = 2;
    SUB_SIZE = NUM_RX / NUM_SUB;
end

OSNR_matrix      = zeros(NUM_RX, Nsd);
BER_matrix       = zeros(NUM_RX, Nsd);
PTR_OSNR_vector  = zeros(Nsd,1);
PTR_BER_vector   = zeros(Nsd,1);
SCM_OSNR_vector  = zeros(Nsd,1);
SCM_BER_vector   = zeros(Nsd,1);

PTR_syms_all = cell(Nsd,1);
SCM_syms_all = cell(Nsd,1);
SUB_syms_all = cell(Nsd, NUM_SUB);   % ★各音源・各サブアレイの補正後シンボル
SUB_OSNR     = nan(Nsd, NUM_SUB);    % ★各サブアレイ単体のOSNR
SUB_BER      = nan(Nsd, NUM_SUB);    % ★各サブアレイ単体のBER

%% チャネル番号(1..NUM_RX) ⇔ (距離idx, 深度idx) の対応表
% 距離メジャー順（同一距離の深度を連続配置）→ サブアレイ=距離 が自然に成立
ch_rr = zeros(NUM_RX,1);   % 各chの距離インデックス
ch_rz = zeros(NUM_RX,1);   % 各chの深度インデックス
cc = 0;
for irr = 1:Nrr
    for irz = 1:Nrz
        cc = cc + 1;
        ch_rr(cc) = irr;
        ch_rz(cc) = irz;
    end
end

for idx_sd = 1:Nsd

    %% エリア①：送信信号
    Rs = 500; Sps = 16; Fs = Rs * Sps; numSymbols = 2000;

    bits       = randi([0 1], 2 * numSymbols, 1);
    inPhase    = 2*bits(1:2:end) - 1;
    quadrature = 2*bits(2:2:end) - 1;
    symbols    = inPhase + 1i*quadrature;
    sigPower   = mean(abs(symbols).^2);

    txSignal = zeros(numSymbols * Sps, 1);
    for i = 1:numSymbols
        txSignal((i-1)*Sps + 1 : i*Sps) = symbols(i);
    end

    %% 全受信機共通の最小遅延（全距離×全深度を走査）
    global_min_delay = inf;
    for ch = 1:NUM_RX
        d = real(Arr(ch_rr(ch), ch_rz(ch), idx_sd).delay); d = d(d>0);
        if ~isempty(d), global_min_delay = min(global_min_delay, min(d)); end
    end

    ptr_buf_len = length(txSignal) * 4;

    total_ptr_signal = zeros(ptr_buf_len,1);
    total_autocorr   = zeros(ptr_buf_len,1);
    sub_ptr_signal   = zeros(ptr_buf_len, NUM_SUB);
    sub_autocorr     = zeros(ptr_buf_len, NUM_SUB);

    for ch = 1:NUM_RX

        %% エリア②：チャネル適用・ノイズ付加
        amplitudes = Arr(ch_rr(ch), ch_rz(ch), idx_sd).A;
        delays     = real(Arr(ch_rr(ch), ch_rz(ch), idx_sd).delay);
        v = (delays>0); amplitudes=amplitudes(v); delays=delays(v);
        delays = delays - global_min_delay;
        v = (delays>=0); amplitudes=amplitudes(v); delays=delays(v);
        if isempty(delays), continue; end

        maxDelaySamples = round(max(delays)*Fs);
        rxBase = zeros(length(txSignal)+maxDelaySamples,1);
        for p=1:length(delays)
            ds = round(delays(p)*Fs);
            rxBase(ds+(1:length(txSignal))) = rxBase(ds+(1:length(txSignal))) + amplitudes(p)*txSignal;
        end
        SNR_dB=20; sp=mean(abs(rxBase).^2); np=sp/10^(SNR_dB/10);
        rxSignal = rxBase + sqrt(np/2)*(randn(size(rxBase))+1i*randn(size(rxBase)));

        %% --- 単一chの参考OSNR/BER ---
        so = round(Sps/2);
        rxSym = rxSignal(so : Sps : so + Sps*(numSymbols-1));
        a1 = (symbols'*rxSym)/(symbols'*symbols); rxc = rxSym/a1;
        OSNR_matrix(ch,idx_sd) = 10*log10(sigPower/mean(abs(symbols-rxc).^2));
        db = zeros(2*numSymbols,1); db(1:2:end)=real(rxc)>0; db(2:2:end)=imag(rxc)>0;
        BER_matrix(ch,idx_sd) = sum(bits~=db)/length(bits);

        %% --- PTR個別処理 ---
        ds2 = round(delays*Fs)+1;
        h = zeros(max(ds2),1); h(ds2)=amplitudes;
        g = conj(flipud(h));
        z = conv(rxSignal,g);
        hac = conv(h,g);

        Lz = min(length(z),ptr_buf_len);
        La = min(length(hac),ptr_buf_len);

        % 通常PTR：全ch合算
        total_ptr_signal(1:Lz) = total_ptr_signal(1:Lz) + z(1:Lz);
        total_autocorr(1:La)   = total_autocorr(1:La)   + hac(1:La);

        % SCM：このchが属するサブアレイ（=距離）に振り分け
        sub_idx = ch_rr(ch);          % 距離インデックスがそのままサブアレイ番号
        if sub_idx > NUM_SUB, sub_idx = ceil(ch/SUB_SIZE); end  % 旧ARR用フォールバック
        sub_ptr_signal(1:Lz, sub_idx) = sub_ptr_signal(1:Lz, sub_idx) + z(1:Lz);
        sub_autocorr(1:La, sub_idx)   = sub_autocorr(1:La, sub_idx)   + hac(1:La);

    end % ◀ 受信機ループ

    %% エリア③-A：通常PTR（全ch直接合成）
    [~, peak_idx] = max(abs(total_autocorr));
    ptr_offset    = peak_idx + round(Sps/2) - 1;
    ptrSymbols    = total_ptr_signal(ptr_offset : Sps : ptr_offset + Sps*(numSymbols-1));
    a_ptr = (symbols'*ptrSymbols)/(symbols'*symbols); ptrc = ptrSymbols/a_ptr;
    PTR_OSNR_vector(idx_sd) = 10*log10(sigPower/mean(abs(symbols-ptrc).^2));
    db = zeros(2*numSymbols,1); db(1:2:end)=real(ptrc)>0; db(2:2:end)=imag(ptrc)>0;
    PTR_BER_vector(idx_sd) = sum(bits~=db)/length(bits);
    PTR_syms_all{idx_sd} = ptrc;

    %% エリア③-B：SCM-PTR（サブアレイごとにPTR→MRC合成）
    scm_combined = zeros(numSymbols,1);
    weight_sum   = 0;
    for s = 1:NUM_SUB
        [~, pk_s] = max(abs(sub_autocorr(:,s)));
        off_s     = pk_s + round(Sps/2) - 1;
        if off_s + Sps*(numSymbols-1) > ptr_buf_len, continue; end
        subSym    = sub_ptr_signal(off_s : Sps : off_s + Sps*(numSymbols-1), s);
        a_s   = (symbols'*subSym)/(symbols'*symbols);
        subC  = subSym / a_s;
        err_s = mean(abs(symbols - subC).^2);
        w_s   = 1 / err_s;                 % MRC重み（残留誤差の逆数）
        scm_combined = scm_combined + w_s * subC;
        weight_sum   = weight_sum + w_s;

        % ★各サブアレイ単体の星座図用にシンボル・OSNR・BERを保存
        SUB_syms_all{idx_sd, s} = subC;
        SUB_OSNR(idx_sd, s) = 10*log10(sigPower/err_s);
        dbs = zeros(2*numSymbols,1); dbs(1:2:end)=real(subC)>0; dbs(2:2:end)=imag(subC)>0;
        SUB_BER(idx_sd, s)  = sum(bits~=dbs)/length(bits);
    end
    scm_combined = scm_combined / weight_sum;

    a_scm = (symbols'*scm_combined)/(symbols'*symbols);
    scmC  = scm_combined / a_scm;
    SCM_OSNR_vector(idx_sd) = 10*log10(sigPower/mean(abs(symbols-scmC).^2));
    db = zeros(2*numSymbols,1); db(1:2:end)=real(scmC)>0; db(2:2:end)=imag(scmC)>0;
    SCM_BER_vector(idx_sd) = sum(bits~=db)/length(bits);
    SCM_syms_all{idx_sd} = scmC;

end % ◀ 音源ループ

%% 選択ダイバーシティ（参考）
[selected_BER,  best_ber_idx]  = min(BER_matrix, [], 1);
[selected_OSNR, best_osnr_idx] = max(OSNR_matrix, [], 1);

%% サブアレイ（距離群）のラベル作成
sub_label = cell(NUM_SUB,1);
for s = 1:NUM_SUB
    if NUM_SUB == Nrr && isfield(Pos,'r') && isfield(Pos.r,'r') && numel(Pos.r.r) >= s
        sub_label{s} = sprintf('%.0fm群', Pos.r.r(s));
    else
        sub_label{s} = sprintf('サブアレイ%d', s);
    end
end

%% 結果出力
disp('==========================================================================');
disp(' 比較：【選択ダイバーシティ】 vs 【通常PTR】 vs 【SCM-PTR】');
fprintf(' （サブアレイ数 NUM_SUB=%d, 各%dch, シンボルレート%dbaud）\n', NUM_SUB, SUB_SIZE, Rs);
disp('==========================================================================');
for idx_sd = 1:Nsd
    fprintf('【音源 %d】\n', idx_sd);
    fprintf('  ➔ 選択ダイバーシティ (Best 1ch):\n');
    fprintf('     - 最良BER : %f (ch %d)\n',      selected_BER(idx_sd),  best_ber_idx(idx_sd));
    fprintf('     - 最高OSNR: %.2f dB (ch %d)\n', selected_OSNR(idx_sd), best_osnr_idx(idx_sd));
    fprintf('  ➔ 通常PTR (全ch直接合成):\n');
    fprintf('     - OSNR: %.2f dB   BER: %f\n', PTR_OSNR_vector(idx_sd), PTR_BER_vector(idx_sd));
    fprintf('  ➔ サブアレイ単体PTR:\n');
    for s = 1:NUM_SUB
        fprintf('     - %-10s OSNR: %6.2f dB   BER: %f\n', sub_label{s}, SUB_OSNR(idx_sd,s), SUB_BER(idx_sd,s));
    end
    fprintf('  ➔ SCM-PTR (距離ごとサブアレイMRC合成):\n');
    fprintf('     - OSNR: %.2f dB   BER: %f\n\n', SCM_OSNR_vector(idx_sd), SCM_BER_vector(idx_sd));
end

%% 星座図の描画①：通常PTR vs SCM-PTR
qpsk_ideal = [1+1i, 1-1i, -1+1i, -1-1i];
figure('Name', 'Constellation: PTR vs SCM', 'Position', [100 100 300*Nsd 560]);
for idx_sd = 1:Nsd
    ptrc = PTR_syms_all{idx_sd};
    scmC = SCM_syms_all{idx_sd};

    subplot(2, Nsd, idx_sd);
    plot(real(ptrc), imag(ptrc), '.', 'Color', [0.3 0.6 1.0], 'MarkerSize', 3); hold on;
    plot(real(qpsk_ideal), imag(qpsk_ideal), 'r+', 'MarkerSize', 12, 'LineWidth', 2); hold off;
    axis equal; grid on;
    lim = max(3, ceil(max(abs([real(ptrc); imag(ptrc)]))));
    xlim([-lim lim]); ylim([-lim lim]); xlabel('I'); ylabel('Q');
    title(sprintf('PTR 音源%d\nOSNR=%.1f dB, BER=%.3f', idx_sd, PTR_OSNR_vector(idx_sd), PTR_BER_vector(idx_sd)));

    subplot(2, Nsd, Nsd + idx_sd);
    plot(real(scmC), imag(scmC), '.', 'Color', [0.1 0.8 0.3], 'MarkerSize', 3); hold on;
    plot(real(qpsk_ideal), imag(qpsk_ideal), 'r+', 'MarkerSize', 12, 'LineWidth', 2); hold off;
    axis equal; grid on;
    lim = max(3, ceil(max(abs([real(scmC); imag(scmC)]))));
    xlim([-lim lim]); ylim([-lim lim]); xlabel('I'); ylabel('Q');
    title(sprintf('SCM-PTR 音源%d\nOSNR=%.1f dB, BER=%.3f', idx_sd, SCM_OSNR_vector(idx_sd), SCM_BER_vector(idx_sd)));
end
sgtitle(sprintf('QPSK 星座図比較  (NUM\\_SUB=%d, Rs=%d baud)', NUM_SUB, Rs), 'FontSize', 13, 'FontWeight', 'bold');

%% ★星座図の描画②：距離群（サブアレイ）ごとに1枚、各図に全音源をまとめる
%   1枚 = 1距離群（例：75m群／100m群 = 2枚）。図内に音源1〜Nsdをサブプロット配置
sub_colors = lines(NUM_SUB);
% 音源サブプロットのレイアウト（なるべく正方形に近い行列）
n_col = ceil(sqrt(Nsd));
n_row = ceil(Nsd / n_col);
for s = 1:NUM_SUB
    figure('Name', sprintf('%s 各音源の星座図', sub_label{s}), ...
           'Position', [100 + 40*s, 100, 320*n_col, 300*n_row]);
    for idx_sd = 1:Nsd
        subC = SUB_syms_all{idx_sd, s};
        subplot(n_row, n_col, idx_sd);
        if isempty(subC)
            axis off;
            title(sprintf('音源%d（データなし）', idx_sd));
            continue;
        end
        plot(real(subC), imag(subC), '.', 'Color', sub_colors(s,:), 'MarkerSize', 4); hold on;
        plot(real(qpsk_ideal), imag(qpsk_ideal), 'r+', 'MarkerSize', 14, 'LineWidth', 2.5); hold off;
        axis equal; grid on;
        lim = max(3, ceil(max(abs([real(subC); imag(subC)]))));
        xlim([-lim lim]); ylim([-lim lim]); xlabel('I'); ylabel('Q');
        title(sprintf('音源%d\nOSNR=%.1f dB, BER=%.3f', ...
              idx_sd, SUB_OSNR(idx_sd,s), SUB_BER(idx_sd,s)), 'FontSize', 11);
    end
    sgtitle(sprintf('%s の QPSK 星座図（全音源）  (Rs=%d baud)', sub_label{s}, Rs), ...
            'FontSize', 13, 'FontWeight', 'bold');
end
