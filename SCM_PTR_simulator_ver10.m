%% ============================================================
%  SCM-PTR (Subarray Combining Passive Time Reversal)
%  通常PTR（16ch直接合成）と SCM-PTR（サブアレイ分割→MRC合成）を比較
%  ※DFEなし。PTR/アレイ合成のみでOSNRを評価する版
% ============================================================
[Arr, Pos] = read_arrivals_asc('JamTank_Arr.arr');

idx_r = 1;
NUM_RX = 16;            % 全受信機数

%% ★SCM設定：サブアレイの分割
% 16chを NUM_SUB 個のサブアレイに均等分割（16はNUM_SUBで割り切れること）
NUM_SUB  = 4;                       % サブアレイ数（例：4 → 各4ch）
SUB_SIZE = NUM_RX / NUM_SUB;        % 1サブアレイあたりの受信機数

OSNR_matrix      = zeros(NUM_RX,4); % 各単一chのOSNR（参考）
BER_matrix       = zeros(NUM_RX,4);
PTR_OSNR_vector  = zeros(4,1);      % 通常PTR（16ch直接合成）
PTR_BER_vector   = zeros(4,1);
SCM_OSNR_vector  = zeros(4,1);      % SCM-PTR（サブアレイMRC合成）
SCM_BER_vector   = zeros(4,1);

% 星座図描画用：各音源の復調後シンボルを保存
PTR_syms_all = cell(4,1);
SCM_syms_all = cell(4,1);

for idx_sd = 1:4

    %% エリア①：送信信号
    Rs = 500; Sps = 16; Fs = Rs * Sps; numSymbols = 2000;   % ★低レートで評価

    bits       = randi([0 1], 2 * numSymbols, 1);
    inPhase    = 2*bits(1:2:end) - 1;
    quadrature = 2*bits(2:2:end) - 1;
    symbols    = inPhase + 1i*quadrature;
    sigPower   = mean(abs(symbols).^2);

    txSignal = zeros(numSymbols * Sps, 1);
    for i = 1:numSymbols
        txSignal((i-1)*Sps + 1 : i*Sps) = symbols(i);
    end

    %% 全受信機共通の最小遅延
    global_min_delay = inf;
    for k = 1:NUM_RX
        d = real(Arr(idx_r,k,idx_sd).delay); d = d(d>0);
        if ~isempty(d), global_min_delay = min(global_min_delay, min(d)); end
    end

    ptr_buf_len = length(txSignal) * 4;

    % 通常PTR用バケツ（全ch合算）
    total_ptr_signal = zeros(ptr_buf_len,1);
    total_autocorr   = zeros(ptr_buf_len,1);

    % SCM用バケツ（サブアレイごとに分けて保持）
    sub_ptr_signal = zeros(ptr_buf_len, NUM_SUB);
    sub_autocorr   = zeros(ptr_buf_len, NUM_SUB);

    for idx_rd = 1:NUM_RX

        %% エリア②：チャネル適用・ノイズ付加
        amplitudes = Arr(idx_r,idx_rd,idx_sd).A;
        delays     = real(Arr(idx_r,idx_rd,idx_sd).delay);
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
        OSNR_matrix(idx_rd,idx_sd) = 10*log10(sigPower/mean(abs(symbols-rxc).^2));
        db = zeros(2*numSymbols,1); db(1:2:end)=real(rxc)>0; db(2:2:end)=imag(rxc)>0;
        BER_matrix(idx_rd,idx_sd) = sum(bits~=db)/length(bits);

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

        % SCM：このchが属するサブアレイに振り分け
        sub_idx = ceil(idx_rd / SUB_SIZE);
        sub_ptr_signal(1:Lz, sub_idx) = sub_ptr_signal(1:Lz, sub_idx) + z(1:Lz);
        sub_autocorr(1:La, sub_idx)   = sub_autocorr(1:La, sub_idx)   + hac(1:La);

    end % ◀ 受信機ループ

    %% エリア③-A：通常PTR（16ch直接合成）の復調
    [~, peak_idx] = max(abs(total_autocorr));
    ptr_offset    = peak_idx + round(Sps/2) - 1;
    ptrSymbols    = total_ptr_signal(ptr_offset : Sps : ptr_offset + Sps*(numSymbols-1));
    a_ptr = (symbols'*ptrSymbols)/(symbols'*symbols); ptrc = ptrSymbols/a_ptr;
    PTR_OSNR_vector(idx_sd) = 10*log10(sigPower/mean(abs(symbols-ptrc).^2));
    db = zeros(2*numSymbols,1); db(1:2:end)=real(ptrc)>0; db(2:2:end)=imag(ptrc)>0;
    PTR_BER_vector(idx_sd) = sum(bits~=db)/length(bits);
    PTR_syms_all{idx_sd} = ptrc;   % 星座図用に保存

    %% エリア③-B：SCM-PTR（サブアレイごとにPTR→MRC合成）
    scm_combined = zeros(numSymbols,1);   % MRC合成後シンボル
    weight_sum   = 0;

    for s = 1:NUM_SUB
        % サブアレイsのq関数ピークでシンボル抽出
        [~, pk_s]   = max(abs(sub_autocorr(:,s)));
        off_s       = pk_s + round(Sps/2) - 1;
        last_s      = off_s + Sps*(numSymbols-1);
        if last_s > ptr_buf_len, continue; end
        subSym      = sub_ptr_signal(off_s : Sps : off_s + Sps*(numSymbols-1), s);

        % サブアレイの複素ゲインalpha_sと残留誤差を推定
        a_s   = (symbols'*subSym)/(symbols'*symbols);
        subC  = subSym / a_s;                     % 送信シンボルにスケール合わせ
        err_s = mean(abs(symbols - subC).^2);     % 残留誤差電力（≒雑音+ISI）

        % ★MRC重み：信号対雑音比に比例（誤差が小さいサブアレイを重視）
        w_s = 1 / err_s;

        scm_combined = scm_combined + w_s * subC;
        weight_sum   = weight_sum + w_s;
    end

    scm_combined = scm_combined / weight_sum;   % 正規化

    % SCM合成後のOSNR/BER（最終ゲイン補正）
    a_scm = (symbols'*scm_combined)/(symbols'*symbols);
    scmC  = scm_combined / a_scm;
    SCM_OSNR_vector(idx_sd) = 10*log10(sigPower/mean(abs(symbols-scmC).^2));
    db = zeros(2*numSymbols,1); db(1:2:end)=real(scmC)>0; db(2:2:end)=imag(scmC)>0;
    SCM_BER_vector(idx_sd) = sum(bits~=db)/length(bits);
    SCM_syms_all{idx_sd} = scmC;   % 星座図用に保存

end % ◀ 音源ループ

%% 選択ダイバーシティ（参考）
[selected_BER,  best_ber_idx]  = min(BER_matrix, [], 1);
[selected_OSNR, best_osnr_idx] = max(OSNR_matrix, [], 1);

%% 結果出力
clc;
disp('==========================================================================');
disp(' 比較：【選択ダイバーシティ】 vs 【通常PTR(16ch)】 vs 【SCM-PTR】');
fprintf(' （サブアレイ数 NUM_SUB=%d, 各%dch, シンボルレート%dbaud）\n', NUM_SUB, SUB_SIZE, Rs);
disp('==========================================================================');
for idx_sd = 1:4
    fprintf('【音源 %d】\n', idx_sd);
    fprintf('  ➔ 選択ダイバーシティ (Best 1ch):\n');
    fprintf('     - 最良BER : %f (受信機 %d)\n',      selected_BER(idx_sd),  best_ber_idx(idx_sd));
    fprintf('     - 最高OSNR: %.2f dB (受信機 %d)\n', selected_OSNR(idx_sd), best_osnr_idx(idx_sd));
    fprintf('  ➔ 通常PTR (16ch直接合成):\n');
    fprintf('     - OSNR: %.2f dB   BER: %f\n', PTR_OSNR_vector(idx_sd), PTR_BER_vector(idx_sd));
    fprintf('  ➔ SCM-PTR (サブアレイMRC合成):\n');
    fprintf('     - OSNR: %.2f dB   BER: %f\n\n', SCM_OSNR_vector(idx_sd), SCM_BER_vector(idx_sd));
end

%% 星座図の描画（4音源 × 2手法 = 8プロット、2×4レイアウト）
% 理想QPSK点（参照）
qpsk_ideal = [1+1i, 1-1i, -1+1i, -1-1i];

figure('Name', 'Constellation Diagrams', 'Position', [100 100 1200 560]);

for idx_sd = 1:4
    ptrc = PTR_syms_all{idx_sd};
    scmC = SCM_syms_all{idx_sd};

    % --- 通常PTR ---
    subplot(2, 4, idx_sd);
    plot(real(ptrc), imag(ptrc), '.', 'Color', [0.3 0.6 1.0], 'MarkerSize', 3);
    hold on;
    plot(real(qpsk_ideal), imag(qpsk_ideal), 'r+', 'MarkerSize', 12, 'LineWidth', 2);
    hold off;
    axis equal; grid on;
    lim = max(3, ceil(max(abs([real(ptrc); imag(ptrc)]))));
    xlim([-lim lim]); ylim([-lim lim]);
    xlabel('I'); ylabel('Q');
    title(sprintf('PTR 音源%d\nOSNR=%.1f dB, BER=%.3f', ...
        idx_sd, PTR_OSNR_vector(idx_sd), PTR_BER_vector(idx_sd)));

    % --- SCM-PTR ---
    subplot(2, 4, 4 + idx_sd);
    plot(real(scmC), imag(scmC), '.', 'Color', [0.1 0.8 0.3], 'MarkerSize', 3);
    hold on;
    plot(real(qpsk_ideal), imag(qpsk_ideal), 'r+', 'MarkerSize', 12, 'LineWidth', 2);
    hold off;
    axis equal; grid on;
    lim = max(3, ceil(max(abs([real(scmC); imag(scmC)]))));
    xlim([-lim lim]); ylim([-lim lim]);
    xlabel('I'); ylabel('Q');
    title(sprintf('SCM-PTR 音源%d\nOSNR=%.1f dB, BER=%.3f', ...
        idx_sd, SCM_OSNR_vector(idx_sd), SCM_BER_vector(idx_sd)));
end

% 上段・下段にラベル
annotation('textbox', [0.01 0.93 0.1 0.05], 'String', '【通常PTR】', ...
    'FontSize', 11, 'FontWeight', 'bold', 'EdgeColor', 'none');
annotation('textbox', [0.01 0.45 0.1 0.05], 'String', '【SCM-PTR】', ...
    'FontSize', 11, 'FontWeight', 'bold', 'EdgeColor', 'none', 'Color', [0 0.5 0]);

sgtitle(sprintf('QPSK 星座図比較  (NUM\\_SUB=%d, Rs=%d baud)', NUM_SUB, Rs), ...
    'FontSize', 13, 'FontWeight', 'bold');
