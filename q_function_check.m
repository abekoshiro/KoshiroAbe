%% ============================================================
%  q_function_check.m
%  PTRの集束品質を q-function で検証する（技術者提案の手法）
%    ① 各受波器で観測CIR h_m に、自身の時間反転(共役) を畳み込む
%         = 自己相関 conv(h_m, conj(flip(h_m)))
%    ② 全受波器で足し合わせる  → q(t) = Σ_m h_m ⊛ h_m*(-t)
%    ③ デルタ関数状になるか検証（中央だけ鋭く、他が沈むほど良い）
%  ※各自己相関の中心(ゼロ遅延)を揃えるため、全CIRを共通長に整えてから加算
%  ※時間変化なし(合成データ)なので、理想はデルタ。ならなければ集束不良
% ============================================================
arr_file = 'JamTank_Arr.arr';     % ★使用中のARRに合わせる
[Arr, Pos] = read_arrivals_asc(arr_file);
Nrr=size(Arr,1); Nrz=size(Arr,2); Nsd=size(Arr,3); NUM_RX=Nrr*Nrz;
fprintf('ARR: 距離%d × 深度%d = %dch, 音源%d\n', Nrr,Nrz,NUM_RX,Nsd);

%% ==== 設定 ====
Fs = 80000;          % CIR構築の標本化周波数[Hz]（遅延分解能）
Rs = 5000;           % 参考: シンボル間隔PSRを測るレート
zoom_ms = 5;         % q関数の表示範囲（焦点±zoom_ms）

ch_rr=zeros(NUM_RX,1); ch_rz=zeros(NUM_RX,1); cc=0;
for irr=1:Nrr, for irz=1:Nrz, cc=cc+1; ch_rr(cc)=irr; ch_rz(cc)=irz; end, end

figure('Name','q-function check','Position',[80 80 360*Nsd 640]);

for idx_sd = 1:Nsd
    %% 全受波器共通の最小遅延・共通長を決める（★中心を揃えるため）
    gmin=inf; maxlen=0;
    for ch=1:NUM_RX
        dl=real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay); dl=dl(dl>0);
        if ~isempty(dl), gmin=min(gmin,min(dl)); end
    end
    for ch=1:NUM_RX
        dl=real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay); dl=dl(dl>0)-gmin; dl=dl(dl>=0);
        if ~isempty(dl), maxlen=max(maxlen, round(max(dl)*Fs)+1); end
    end
    Lmax=maxlen;                       % 全受波器共通のCIR長

    %% ①② 各受波器：CIR ⊛ 時間反転共役 → 全受波器で加算
    q = zeros(2*Lmax-1, 1);
    for ch=1:NUM_RX
        amp=Arr(ch_rr(ch),ch_rz(ch),idx_sd).A;
        dl =real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay);
        v=(dl>0); amp=amp(v); dl=dl(v); dl=dl-gmin; v=(dl>=0); amp=amp(v); dl=dl(v);
        if isempty(dl), continue; end
        h = zeros(Lmax,1);             % ★共通長Lmaxに揃える（ゼロ埋め）
        ds = round(dl*Fs)+1; ds=ds(ds>=1 & ds<=Lmax);
        for p=1:numel(ds), h(ds(p)) = h(ds(p)) + amp(p); end
        q = q + conv(h, conj(flipud(h)));   % 自己相関（中心が常に index=Lmax）
    end

    %% ③ 検証：デルタらしさの指標
    [qpk, pk] = max(abs(q));           % ピーク（理想は中心=Lmax）
    qn = q/qpk;                        % ピーク正規化
    % 全サンプルでの ピーク対サイドローブ（積分PSR）
    side_all = sum(abs(q).^2) - qpk^2;
    PSR_samp = 10*log10(qpk^2 / max(side_all,eps));
    % ピーク対「最大サイドローブ」比（PSLR）: 中央1点を除く最大
    tmp=abs(qn); tmp(pk)=0; PSLR = -20*log10(max(tmp)+eps);
    % シンボル間隔PSR（通信で効くISIの目安）
    Sps=round(Fs/Rs); lags=-round(Lmax/Sps):round(Lmax/Sps);
    qs=zeros(numel(lags),1);
    for j=1:numel(lags), idx=pk+lags(j)*Sps; if idx>=1&&idx<=length(q), qs(j)=q(idx); end, end
    PSR_sym = 10*log10(abs(qs(lags==0))^2 / max(sum(abs(qs(lags~=0)).^2),eps));

    fprintf('\n=== 音源%d ===\n', idx_sd);
    fprintf('  ピーク位置: index %d（中心=%d なら整列OK）\n', pk, Lmax);
    fprintf('  最大サイドローブ比 PSLR   = %6.2f dB（大きいほどデルタに近い）\n', PSLR);
    fprintf('  積分PSR(全サンプル)       = %6.2f dB\n', PSR_samp);
    fprintf('  シンボル間隔PSR(Rs=%d)  = %6.2f dB（通信ISIの目安=OSNR上限）\n', Rs, PSR_sym);

    %% 描画（時間, dB）
    t = ((0:length(q)-1)-(pk-1))/Fs*1e3;   % 焦点を0msに
    subplot(2,Nsd,idx_sd);
    plot(t, 20*log10(abs(qn)+eps)); grid on; xlim([-zoom_ms zoom_ms]); ylim([-60 2]);
    xlabel('焦点からの時間 [ms]'); ylabel('|q|/|q_{peak}| [dB]');
    title(sprintf('音源%d q関数（dB）\nPSLR=%.1fdB, sym-PSR=%.1fdB', idx_sd, PSLR, PSR_sym));

    subplot(2,Nsd,Nsd+idx_sd);
    plot(t, abs(qn)); grid on; xlim([-zoom_ms zoom_ms]); ylim([0 1.05]);
    xlabel('焦点からの時間 [ms]'); ylabel('|q|/|q_{peak}|（線形）');
    title(sprintf('音源%d q関数（線形）\n中央だけ1・他0が理想（デルタ）', idx_sd));
end
sgtitle(sprintf('q-function 検証（%dch合成, 雑音なし）', NUM_RX), 'FontWeight','bold');
