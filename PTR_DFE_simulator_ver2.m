%% ============================================================
%  PTR_DFE_simulator_ver2.m
%  SCM-PTR + RLS-DFE。ver1に「トレーニング区間を除いたOSNR/BER」を追加。
%  ★過学習チェック：
%    DFEはトレーニング区間(先頭Ntrain)で係数を学習するため、全区間でOSNRを
%    測ると学習済みシンボルに合わせ込んだ「見かけの改善」が混ざりうる。
%    そこで判定指向区間(Ntrain+1 以降＝未知データ相当)だけでも評価し、
%    そこでも改善が残るかを確認する。→ 残れば「本物の改善」。
%  比較: 通常PTR / SCM-PTR / SCM-PTR + DFE を [全区間] と [判定指向区間] で表示
% ============================================================
[Arr, Pos] = read_arrivals_asc('JamTank_Arr_ver3.arr');

%% ARR次元の自動取得
Nrr = size(Arr,1); Nrz = size(Arr,2); Nsd = size(Arr,3);
NUM_RX = Nrr*Nrz;

fprintf('ARR構成: 受信距離 %d点 × 受信深度 %d点 = %dch, 音源 %d個\n', Nrr, Nrz, NUM_RX, Nsd);
if isfield(Pos,'r')&&isfield(Pos.r,'r'), fprintf('  受信距離: '); fprintf('%.1f ',Pos.r.r); fprintf('m\n'); end
if isfield(Pos,'r')&&isfield(Pos.r,'z'), fprintf('  受信深度: '); fprintf('%.2f ',Pos.r.z); fprintf('m\n\n'); end

%% ch ⇔ (距離,深度) 対応
ch_rr=zeros(NUM_RX,1); ch_rz=zeros(NUM_RX,1); cc=0;
for irr=1:Nrr, for irz=1:Nrz, cc=cc+1; ch_rr(cc)=irr; ch_rz(cc)=irz; end, end

%% サブアレイ分割とch→サブアレイ対応（距離が1つでも動くようにする）
if Nrr >= 2
    NUM_SUB = Nrr; SUB_SIZE = Nrz;
    ch_sub  = ch_rr;
else
    NUM_SUB = 2; SUB_SIZE = ceil(NUM_RX/NUM_SUB);
    ch_sub  = min(NUM_SUB, ceil((1:NUM_RX)'/SUB_SIZE));
end
fprintf('サブアレイ数 NUM_SUB=%d（距離%d, 深度%d）\n', NUM_SUB, Nrr, Nrz);

%% ★DFEパラメータ（tap_length_diagnostic.m の推奨値に合わせる）
Nf     = 8;       % フィードフォワードタップ数
Nb     = 6;       % フィードバックタップ数
lambda = 0.999;   % RLS忘却係数
delta  = 0.01;    % RLS正則化
Ntrain = 400;     % トレーニングシンボル数
% ※目安: Ntrain ≳ 数×(Nf+Nb)。タップを増やしたらNtrainも増やすこと

%% 結果格納  [全区間] と [判定指向区間(test)] の2種
PTR_OSNR=zeros(Nsd,1);      PTR_BER=zeros(Nsd,1);
SCM_OSNR=zeros(Nsd,1);      SCM_BER=zeros(Nsd,1);
DFE_OSNR=zeros(Nsd,1);      DFE_BER=zeros(Nsd,1);
SCM_OSNR_t=zeros(Nsd,1);    SCM_BER_t=zeros(Nsd,1);   % test区間
DFE_OSNR_t=zeros(Nsd,1);    DFE_BER_t=zeros(Nsd,1);   % test区間
SCM_syms=cell(Nsd,1); DFE_syms=cell(Nsd,1);

for idx_sd = 1:Nsd

    %% 送信信号
    Rs=500; Sps=16; Fs=Rs*Sps; numSymbols=2000;
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
    total_ptr=zeros(ptr_buf_len,1); total_ac=zeros(ptr_buf_len,1);
    sub_ptr=zeros(ptr_buf_len,NUM_SUB); sub_ac=zeros(ptr_buf_len,NUM_SUB);

    for ch=1:NUM_RX
        %% チャネル適用・雑音
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
        SNR_dB=20; sp=mean(abs(rxBase).^2); np=sp/10^(SNR_dB/10);
        rx=rxBase+sqrt(np/2)*(randn(size(rxBase))+1i*randn(size(rxBase)));

        %% PTR個別
        ds2=round(dl*Fs)+1; h=zeros(max(ds2),1); h(ds2)=amp;
        g=conj(flipud(h)); z=conv(rx,g); hac=conv(h,g);
        Lz=min(length(z),ptr_buf_len); La=min(length(hac),ptr_buf_len);
        total_ptr(1:Lz)=total_ptr(1:Lz)+z(1:Lz);
        total_ac(1:La) =total_ac(1:La)+hac(1:La);
        si=ch_sub(ch);
        sub_ptr(1:Lz,si)=sub_ptr(1:Lz,si)+z(1:Lz);
        sub_ac(1:La,si) =sub_ac(1:La,si)+hac(1:La);
    end

    all_idx  = 1:numSymbols;               % 全区間
    test_idx = (Ntrain+1):numSymbols;      % 判定指向区間（未知データ相当）

    %% --- 通常PTR ---
    [~,pk]=max(abs(total_ac)); off=pk+round(Sps/2)-1;
    ptrSym=total_ptr(off:Sps:off+Sps*(numSymbols-1));
    [PTR_OSNR(idx_sd),PTR_BER(idx_sd)] = eval_region(symbols, ptrSym, all_idx);

    %% --- SCM-PTR（MRC合成） ---
    scm=zeros(numSymbols,1); wsum=0;
    for s=1:NUM_SUB
        [~,pks]=max(abs(sub_ac(:,s))); offs=pks+round(Sps/2)-1;
        if offs+Sps*(numSymbols-1)>ptr_buf_len, continue; end
        ss=sub_ptr(offs:Sps:offs+Sps*(numSymbols-1),s);
        as=(symbols'*ss)/(symbols'*symbols); sc=ss/as;
        es=mean(abs(symbols-sc).^2); ws=1/es;
        scm=scm+ws*sc; wsum=wsum+ws;
    end
    scm=scm/wsum;                                  % ← DFE入力（生の合成シンボル）
    [SCM_OSNR(idx_sd),  SCM_BER(idx_sd)]   = eval_region(symbols, scm, all_idx);
    [SCM_OSNR_t(idx_sd),SCM_BER_t(idx_sd)] = eval_region(symbols, scm, test_idx);
    [~,~,scmC] = eval_region(symbols, scm, all_idx);
    SCM_syms{idx_sd}=scmC;

    %% --- SCM-PTR + RLS-DFE ---
    dfe_out = rls_dfe(scm, symbols, Nf, Nb, lambda, delta, Ntrain);
    [DFE_OSNR(idx_sd),  DFE_BER(idx_sd)]   = eval_region(symbols, dfe_out, all_idx);
    [DFE_OSNR_t(idx_sd),DFE_BER_t(idx_sd)] = eval_region(symbols, dfe_out, test_idx);
    [~,~,dfeC] = eval_region(symbols, dfe_out, all_idx);
    DFE_syms{idx_sd}=dfeC;

end % ◀ 音源ループ

%% 結果出力
disp('==========================================================================');
disp(' 比較：通常PTR / SCM-PTR / SCM-PTR+DFE');
fprintf(' NUM_SUB=%d, Rs=%dbaud, DFE: Nf=%d Nb=%d λ=%.3f Ntrain=%d\n', NUM_SUB, Rs, Nf, Nb, lambda, Ntrain);
fprintf(' ※[全]=全区間, [test]=判定指向区間(先頭%dシンボルのトレーニングを除外)\n', Ntrain);
fprintf(' ※[test]は未知データ相当。ここで改善が残れば過学習でない「本物」\n');
disp('==========================================================================');
for idx_sd=1:Nsd
    fprintf('【音源 %d】\n', idx_sd);
    fprintf('  通常PTR             : OSNR %6.2f dB   BER %f\n', PTR_OSNR(idx_sd), PTR_BER(idx_sd));
    fprintf('  SCM-PTR      [全]   : OSNR %6.2f dB   BER %f\n', SCM_OSNR(idx_sd), SCM_BER(idx_sd));
    fprintf('  SCM-PTR      [test] : OSNR %6.2f dB   BER %f\n', SCM_OSNR_t(idx_sd), SCM_BER_t(idx_sd));
    fprintf('  SCM-PTR+DFE  [全]   : OSNR %6.2f dB   BER %f\n', DFE_OSNR(idx_sd), DFE_BER(idx_sd));
    fprintf('  SCM-PTR+DFE  [test] : OSNR %6.2f dB   BER %f\n', DFE_OSNR_t(idx_sd), DFE_BER_t(idx_sd));
    fprintf('  ── DFE改善量(test基準, SCM[test]比): ΔOSNR = %+.2f dB\n', ...
            DFE_OSNR_t(idx_sd)-SCM_OSNR_t(idx_sd));
    fprintf('  ── 過学習チェック(DFE 全 − test): %+.2f dB  ', DFE_OSNR(idx_sd)-DFE_OSNR_t(idx_sd));
    if DFE_OSNR(idx_sd)-DFE_OSNR_t(idx_sd) > 1.0
        fprintf('⚠ 全区間が良すぎ→過学習の疑い(Ntrain増やす)\n\n');
    else
        fprintf('✓ 差が小さく健全\n\n');
    end
end

%% 星座図：SCM-PTR vs SCM-PTR+DFE（タイトルにtest区間OSNRを併記）
qpsk_ideal=[1+1i,1-1i,-1+1i,-1-1i];
figure('Name','星座図: SCM vs SCM+DFE','Position',[100 100 300*Nsd 560]);
for idx_sd=1:Nsd
    scmC=SCM_syms{idx_sd}; dfeC=DFE_syms{idx_sd};
    subplot(2,Nsd,idx_sd);
    plot(real(scmC),imag(scmC),'.','Color',[0.1 0.8 0.3],'MarkerSize',3); hold on;
    plot(real(qpsk_ideal),imag(qpsk_ideal),'r+','MarkerSize',12,'LineWidth',2); hold off;
    axis equal; grid on; lim=max(3,ceil(max(abs([real(scmC);imag(scmC)]))));
    xlim([-lim lim]); ylim([-lim lim]); xlabel('I'); ylabel('Q');
    title(sprintf('SCM-PTR 音源%d\ntest OSNR=%.1fdB BER=%.3f',idx_sd,SCM_OSNR_t(idx_sd),SCM_BER_t(idx_sd)));

    subplot(2,Nsd,Nsd+idx_sd);
    plot(real(dfeC),imag(dfeC),'.','Color',[0.85 0.3 0.1],'MarkerSize',3); hold on;
    plot(real(qpsk_ideal),imag(qpsk_ideal),'r+','MarkerSize',12,'LineWidth',2); hold off;
    axis equal; grid on; lim=max(3,ceil(max(abs([real(dfeC);imag(dfeC)]))));
    xlim([-lim lim]); ylim([-lim lim]); xlabel('I'); ylabel('Q');
    title(sprintf('SCM+DFE 音源%d\ntest OSNR=%.1fdB BER=%.3f',idx_sd,DFE_OSNR_t(idx_sd),DFE_BER_t(idx_sd)));
end
sgtitle(sprintf('DFE等化の効果（test区間OSNR表示） Nf=%d Nb=%d Ntrain=%d',Nf,Nb,Ntrain),'FontSize',13,'FontWeight','bold');

%% ================= ローカル関数 =================
function [osnr, ber, comp] = eval_region(symbols, raw, idx)
% 指定区間idxで OSNR/BER を評価（αは区間内で最小二乗推定）
%  symbols: 送信シンボル, raw: 評価対象(未補正), idx: 評価するシンボル番号
    s = symbols(idx); r = raw(idx);
    a = (s'*r)/(s'*s);          % 最小二乗複素利得
    c = r / a;
    osnr = 10*log10(mean(abs(s).^2)/mean(abs(s-c).^2));
    % QPSK象限判定でBER（I/Qビット = 実部/虚部の符号）
    dec = [real(c)>0, imag(c)>0];
    ref = [real(s)>0, imag(s)>0];
    ber = mean(dec(:) ~= ref(:));
    comp = raw / a;             % 全体を同じαで補正して返す（星座図用）
end

function out = rls_dfe(x, symbols, Nf, Nb, lambda, delta, Ntrain)
% RLS判定帰還等化器（QPSK, T間隔）
    N = length(x);
    Ntap = Nf + Nb;
    w = zeros(Ntap,1); w(1)=1;
    P = (1/delta)*eye(Ntap);
    ff = zeros(Nf,1); fb = zeros(Nb,1);
    out = zeros(N,1);
    for n=1:N
        ff = [x(n); ff(1:end-1)];
        u  = [ff; -fb];
        y  = w' * u;
        out(n) = y;
        d  = sign(real(y)) + 1i*sign(imag(y));
        if n <= Ntrain
            s_ref = symbols(n); fb_new = symbols(n);
        else
            s_ref = d;          fb_new = d;
        end
        Pu    = P*u;
        kappa = Pu/(lambda + u'*Pu);
        e     = s_ref - y;
        w     = w + kappa*conj(e);
        P     = (P - kappa*(Pu'))/lambda;
        fb    = [fb_new; fb(1:end-1)];
    end
end
