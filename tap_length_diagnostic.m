%% ============================================================
%  tap_length_diagnostic.m
%  DFE のタップ数(Nf/Nb)の根拠を数値化する診断
%  ・SCM-PTR 合成後の「等価シンボル応答」ρ(k) を相関で実測
%    ρ(k) = E[ y[n]·conj(s[n-k]) ] / E[|s|^2]
%    （y=合成シンボル列, s=送信シンボル。sは白色なのでρがそのまま等価応答）
%  ・|ρ(k)| が閾値を超える lag の範囲から
%      Nf ≳ 前方(lag<0)の有意範囲 + 余裕
%      Nb ≳ 後方(lag>0)の有意範囲
%    を音源ごとに推奨値として出力
% ============================================================
[Arr, Pos] = read_arrivals_asc('JamTank_Arr_ver3.arr');

Nrr = size(Arr,1); Nrz = size(Arr,2); Nsd = size(Arr,3);
NUM_RX = Nrr*Nrz;

ch_rr=zeros(NUM_RX,1); ch_rz=zeros(NUM_RX,1); cc=0;
for irr=1:Nrr, for irz=1:Nrz, cc=cc+1; ch_rr(cc)=irr; ch_rz(cc)=irz; end, end

%% サブアレイ分割とch→サブアレイ対応（距離が1つでも動くようにする）
if Nrr >= 2
    NUM_SUB = Nrr; SUB_SIZE = Nrz;   % 距離ごとに1サブアレイ
    ch_sub  = ch_rr;                 % chの距離idx = サブアレイ番号
else
    NUM_SUB = 2; SUB_SIZE = ceil(NUM_RX/NUM_SUB);  % 単一距離：深度で2分割
    ch_sub  = min(NUM_SUB, ceil((1:NUM_RX)'/SUB_SIZE));
end
fprintf('サブアレイ数 NUM_SUB=%d（距離%d, 深度%d）\n', NUM_SUB, Nrr, Nrz);

Rs = 500; Sps = 16; Fs = Rs*Sps; numSymbols = 2000;

%% ==== 診断パラメータ ====
Lmax   = 30;         % 調べる最大 lag（シンボル）
thr_dB = -20;        % 有意サイドローブの閾値（ピーク比 dB）。-20dB=10%
margin = 3;          % FF に加える余裕タップ

fprintf('=== DFE タップ長診断（閾値 %d dB, 最大lag %d）===\n\n', thr_dB, Lmax);

Nf_rec = zeros(Nsd,1);
Nb_rec = zeros(Nsd,1);
rho_all = cell(Nsd,1);

for idx_sd = 1:Nsd

    %% 送信信号
    bits=randi([0 1],2*numSymbols,1);
    symbols=(2*bits(1:2:end)-1)+1i*(2*bits(2:2:end)-1);
    txSignal=zeros(numSymbols*Sps,1);
    for i=1:numSymbols, txSignal((i-1)*Sps+1:i*Sps)=symbols(i); end

    %% 共通最小遅延
    gmin=inf;
    for ch=1:NUM_RX
        d=real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay); d=d(d>0);
        if ~isempty(d), gmin=min(gmin,min(d)); end
    end

    ptr_buf_len=length(txSignal)*4;
    sub_ptr=zeros(ptr_buf_len,NUM_SUB);
    sub_ac =zeros(ptr_buf_len,NUM_SUB);

    for ch=1:NUM_RX
        amp=Arr(ch_rr(ch),ch_rz(ch),idx_sd).A;
        dl =real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay);
        v=(dl>0); amp=amp(v); dl=dl(v); dl=dl-gmin;
        v=(dl>=0); amp=amp(v); dl=dl(v);
        if isempty(dl), continue; end
        mds=round(max(dl)*Fs);
        rxBase=zeros(length(txSignal)+mds,1);
        for p=1:length(dl)
            ds=round(dl(p)*Fs);
            rxBase(ds+(1:length(txSignal)))=rxBase(ds+(1:length(txSignal)))+amp(p)*txSignal;
        end
        sp=mean(abs(rxBase).^2); np=sp/10^(20/10);
        rx=rxBase+sqrt(np/2)*(randn(size(rxBase))+1i*randn(size(rxBase)));
        ds2=round(dl*Fs)+1; h=zeros(max(ds2),1); h(ds2)=amp;
        g=conj(flipud(h)); z=conv(rx,g); hac=conv(h,g);
        Lz=min(length(z),ptr_buf_len); La=min(length(hac),ptr_buf_len);
        si=ch_sub(ch);
        sub_ptr(1:Lz,si)=sub_ptr(1:Lz,si)+z(1:Lz);
        sub_ac(1:La,si) =sub_ac(1:La,si)+hac(1:La);
    end

    %% SCM-MRC 合成（ver13と同じ）
    scm=zeros(numSymbols,1); wsum=0;
    for s=1:NUM_SUB
        [~,pk]=max(abs(sub_ac(:,s)));
        off=pk+round(Sps/2)-1;
        if off+Sps*(numSymbols-1)>ptr_buf_len, continue; end
        ss=sub_ptr(off:Sps:off+Sps*(numSymbols-1),s);
        a=(symbols'*ss)/(symbols'*symbols); sc=ss/a;
        e=mean(abs(symbols-sc).^2); w=1/e;
        scm=scm+w*sc; wsum=wsum+w;
    end
    scm=scm/wsum;
    a=(symbols'*scm)/(symbols'*symbols); y=scm/a;   % 合成後シンボル（利得正規化）

    %% 等価シンボル応答 ρ(k) を相関で実測
    lags=-Lmax:Lmax;
    rho=zeros(numel(lags),1);
    for ii=1:numel(lags)
        k=lags(ii);
        n1=max(1,1+k); n2=min(numSymbols,numSymbols+k);
        idxN=n1:n2;
        rho(ii)=mean( y(idxN).*conj(symbols(idxN-k)) );
    end
    rho=rho/mean(abs(symbols).^2);
    rho_all{idx_sd}=struct('lags',lags,'rho',rho);

    %% 有意 lag の判定 → Nf/Nb 推奨
    rho_dB = 20*log10(abs(rho)/max(abs(rho))+eps);
    sig = (rho_dB > thr_dB);               % 閾値超えの lag
    sig = sig(:).';                        % lags(行ベクトル)と向きを揃える
    fwd_lags = lags(sig & lags<0);         % 前方（負lag）
    bwd_lags = lags(sig & lags>0);         % 後方（正lag）
    if isempty(fwd_lags), max_fwd = 0; else, max_fwd = max(abs(fwd_lags)); end
    if isempty(bwd_lags), max_bwd = 0; else, max_bwd = max(bwd_lags);      end

    Nf_rec(idx_sd) = max_fwd + margin + 1;  % +1は現在タップ, marginは余裕
    Nb_rec(idx_sd) = max_bwd;

    fprintf('音源%d:\n', idx_sd);
    fprintf('  前方(lag<0)有意範囲: lag -%d まで  → FF寄与\n', max_fwd);
    fprintf('  後方(lag>0)有意範囲: lag +%d まで  → FB寄与\n', max_bwd);
    fprintf('  推奨 Nf = %d,  Nb = %d\n\n', Nf_rec(idx_sd), Nb_rec(idx_sd));
end

%% 全音源で安全側（最大）を採用
fprintf('----------------------------------------------------\n');
fprintf('全音源をカバーする推奨タップ数: Nf = %d,  Nb = %d\n', max(Nf_rec), max(Nb_rec));
fprintf('（PTR_DFE_simulator_ver1.m の Nf/Nb にこの値を設定）\n');

%% ==== 可視化：等価シンボル応答 ρ(k) ====
figure('Name','等価シンボル応答 ρ(k)（DFEタップ根拠）','Position',[80 80 340*Nsd 420]);
for idx_sd=1:Nsd
    L=rho_all{idx_sd}.lags; r=rho_all{idx_sd}.rho;
    rdb=20*log10(abs(r)/max(abs(r))+eps);
    subplot(1,Nsd,idx_sd);
    stem(L, rdb, 'filled','MarkerSize',3); hold on;
    yline(thr_dB,'r--','閾値');
    xline(0,'k:');
    grid on; xlabel('lag [シンボル]'); ylabel('|ρ(k)| [dB]');
    ylim([-40 2]); xlim([-Lmax Lmax]);
    title(sprintf('音源%d  Nf=%d,Nb=%d', idx_sd, Nf_rec(idx_sd), Nb_rec(idx_sd)));
    hold off;
end
sgtitle('SCM-PTR 合成後の残留ISI（等価シンボル応答）＝DFEタップ数の根拠','FontWeight','bold');
