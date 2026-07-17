%% ============================================================
%  PTR_DFE_simulator_ver1.m
%  SCM-PTR（距離ごとサブアレイ→MRC合成）に RLS-DFE 等化を追加
%  ・比較: 通常PTR / SCM-PTR / SCM-PTR + DFE
%  ・DFE は合成後シンボル列（シンボルレート, T間隔）に適用
%  ・RLS採用理由: LMSは収束に多数の反復が必要で本規模では発散するため
%  ※タップ数 Nf/Nb は tap_length_diagnostic.m の推奨値を設定すること
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
    NUM_SUB = Nrr; SUB_SIZE = Nrz;   % 距離ごとに1サブアレイ
    ch_sub  = ch_rr;
else
    NUM_SUB = 2; SUB_SIZE = ceil(NUM_RX/NUM_SUB);  % 単一距離：深度で2分割
    ch_sub  = min(NUM_SUB, ceil((1:NUM_RX)'/SUB_SIZE));
end
fprintf('サブアレイ数 NUM_SUB=%d（距離%d, 深度%d）\n', NUM_SUB, Nrr, Nrz);

%% ★DFEパラメータ（tap_length_diagnostic.m の推奨値に合わせる）
Nf     = 8;       % フィードフォワードタップ数
Nb     = 6;       % フィードバックタップ数
lambda = 0.999;   % RLS忘却係数（1に近いほど過去重視／定常チャネル向き）
delta  = 0.01;    % RLS正則化（P初期値スケール）
Ntrain = 400;     % トレーニングシンボル数（RLSは少数で収束）

%% 結果格納
PTR_OSNR=zeros(Nsd,1);      PTR_BER=zeros(Nsd,1);
SCM_OSNR=zeros(Nsd,1);      SCM_BER=zeros(Nsd,1);
SCMDFE_OSNR=zeros(Nsd,1);   SCMDFE_BER=zeros(Nsd,1);
PTR_syms=cell(Nsd,1); SCM_syms=cell(Nsd,1); DFE_syms=cell(Nsd,1);
RX_syms=cell(Nsd,NUM_RX);              % ★合成前：各受信機単体の生受信シンボル
RX_OSNR=nan(Nsd,NUM_RX); RX_BER=nan(Nsd,NUM_RX);

for idx_sd = 1:Nsd

    %% 送信信号
    Rs=500; Sps=16; Fs=Rs*Sps; numSymbols=2000;
    bits=randi([0 1],2*numSymbols,1);
    symbols=(2*bits(1:2:end)-1)+1i*(2*bits(2:2:end)-1);
    sigPower=mean(abs(symbols).^2);
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

        %% ★合成前：この受信機単体の生受信星座図（PTRなし・中央サンプリング）
        so=round(Sps/2);
        rxSym=rx(so:Sps:so+Sps*(numSymbols-1));
        a1=(symbols'*rxSym)/(symbols'*symbols); rxc=rxSym/a1;
        RX_syms{idx_sd,ch}=rxc;
        RX_OSNR(idx_sd,ch)=10*log10(sigPower/mean(abs(symbols-rxc).^2));
        RX_BER(idx_sd,ch)=ber_qpsk(symbols,rxc,bits);

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

    %% --- 通常PTR ---
    [~,pk]=max(abs(total_ac)); off=pk+round(Sps/2)-1;
    ptrSym=total_ptr(off:Sps:off+Sps*(numSymbols-1));
    a=(symbols'*ptrSym)/(symbols'*symbols); ptrc=ptrSym/a;
    PTR_OSNR(idx_sd)=10*log10(sigPower/mean(abs(symbols-ptrc).^2));
    PTR_BER(idx_sd)=ber_qpsk(symbols,ptrc,bits);
    PTR_syms{idx_sd}=ptrc;

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
    a=(symbols'*scm)/(symbols'*symbols); scmC=scm/a;
    SCM_OSNR(idx_sd)=10*log10(sigPower/mean(abs(symbols-scmC).^2));
    SCM_BER(idx_sd)=ber_qpsk(symbols,scmC,bits);
    SCM_syms{idx_sd}=scmC;

    %% --- SCM-PTR + RLS-DFE ---
    %  入力は α補正前の scm（DFE自身が利得を学習。α二重補正を避ける）
    dfe_out = rls_dfe(scm, symbols, Nf, Nb, lambda, delta, Ntrain);
    ad=(symbols'*dfe_out)/(symbols'*symbols); dfeC=dfe_out/ad;
    SCMDFE_OSNR(idx_sd)=10*log10(sigPower/mean(abs(symbols-dfeC).^2));
    SCMDFE_BER(idx_sd)=ber_qpsk(symbols,dfeC,bits);
    DFE_syms{idx_sd}=dfeC;

end % ◀ 音源ループ

%% サブアレイ距離ラベル
sub_label=cell(NUM_SUB,1);
for s=1:NUM_SUB
    if NUM_SUB==Nrr && isfield(Pos,'r')&&isfield(Pos.r,'r')&&numel(Pos.r.r)>=s
        sub_label{s}=sprintf('%.0fm群',Pos.r.r(s));
    else, sub_label{s}=sprintf('サブアレイ%d',s); end
end

%% 結果出力
disp('==========================================================================');
disp(' 比較：【通常PTR】 vs 【SCM-PTR】 vs 【SCM-PTR + RLS-DFE】');
fprintf(' （NUM_SUB=%d, Rs=%dbaud, DFE: Nf=%d Nb=%d λ=%.3f Ntrain=%d）\n', ...
        NUM_SUB, Rs, Nf, Nb, lambda, Ntrain);
disp('==========================================================================');
for idx_sd=1:Nsd
    fprintf('【音源 %d】\n', idx_sd);
    fprintf('  通常PTR        : OSNR %6.2f dB   BER %f\n', PTR_OSNR(idx_sd), PTR_BER(idx_sd));
    fprintf('  SCM-PTR        : OSNR %6.2f dB   BER %f\n', SCM_OSNR(idx_sd), SCM_BER(idx_sd));
    fprintf('  SCM-PTR + DFE  : OSNR %6.2f dB   BER %f', SCMDFE_OSNR(idx_sd), SCMDFE_BER(idx_sd));
    fprintf('   （ΔOSNR = %+.2f dB）\n\n', SCMDFE_OSNR(idx_sd)-SCM_OSNR(idx_sd));
end

%% 星座図：SCM-PTR vs SCM-PTR+DFE
qpsk_ideal=[1+1i,1-1i,-1+1i,-1-1i];
figure('Name','星座図: SCM vs SCM+DFE','Position',[100 100 300*Nsd 560]);
for idx_sd=1:Nsd
    scmC=SCM_syms{idx_sd}; dfeC=DFE_syms{idx_sd};
    subplot(2,Nsd,idx_sd);
    plot(real(scmC),imag(scmC),'.','Color',[0.1 0.8 0.3],'MarkerSize',3); hold on;
    plot(real(qpsk_ideal),imag(qpsk_ideal),'r+','MarkerSize',12,'LineWidth',2); hold off;
    axis equal; grid on; lim=max(3,ceil(max(abs([real(scmC);imag(scmC)]))));
    xlim([-lim lim]); ylim([-lim lim]); xlabel('I'); ylabel('Q');
    title(sprintf('SCM-PTR 音源%d\nOSNR=%.1fdB BER=%.3f',idx_sd,SCM_OSNR(idx_sd),SCM_BER(idx_sd)));

    subplot(2,Nsd,Nsd+idx_sd);
    plot(real(dfeC),imag(dfeC),'.','Color',[0.85 0.3 0.1],'MarkerSize',3); hold on;
    plot(real(qpsk_ideal),imag(qpsk_ideal),'r+','MarkerSize',12,'LineWidth',2); hold off;
    axis equal; grid on; lim=max(3,ceil(max(abs([real(dfeC);imag(dfeC)]))));
    xlim([-lim lim]); ylim([-lim lim]); xlabel('I'); ylabel('Q');
    title(sprintf('SCM+DFE 音源%d\nOSNR=%.1fdB BER=%.3f',idx_sd,SCMDFE_OSNR(idx_sd),SCMDFE_BER(idx_sd)));
end
sgtitle(sprintf('DFE等化の効果  (Nf=%d, Nb=%d, RLS λ=%.3f)',Nf,Nb,lambda),'FontSize',13,'FontWeight','bold');

%% ★合成前：各受信機単体の生受信星座図（音源ごとに1枚）
%   PTRも合成も通さない「そのまま復調」の星座図。合成でどれだけ改善したか比較用
nc=ceil(sqrt(NUM_RX)); nr=ceil(NUM_RX/nc);
for idx_sd=1:Nsd
    figure('Name',sprintf('合成前 各受信機 星座図 音源%d',idx_sd), ...
           'Position',[80+30*idx_sd, 80, 200*nc, 190*nr]);
    for ch=1:NUM_RX
        rxc=RX_syms{idx_sd,ch};
        subplot(nr,nc,ch);
        if isempty(rxc), axis off; title(sprintf('ch%d なし',ch)); continue; end
        plot(real(rxc),imag(rxc),'.','Color',[0.4 0.4 0.9],'MarkerSize',2); hold on;
        plot(real(qpsk_ideal),imag(qpsk_ideal),'r+','MarkerSize',8,'LineWidth',1.5); hold off;
        axis equal; grid on;
        lim=max(3,ceil(max(abs([real(rxc);imag(rxc)]))));
        xlim([-lim lim]); ylim([-lim lim]);
        % ラベル：距離・深度が取れれば表示
        if isfield(Pos,'r')&&isfield(Pos.r,'r')&&isfield(Pos.r,'z')
            lab=sprintf('%.0fm/%.2fm',Pos.r.r(ch_rr(ch)),Pos.r.z(ch_rz(ch)));
        else
            lab=sprintf('ch%d',ch);
        end
        title(sprintf('%s\nOSNR=%.1f BER=%.2f',lab,RX_OSNR(idx_sd,ch),RX_BER(idx_sd,ch)),'FontSize',8);
    end
    sgtitle(sprintf('合成前：各受信機単体の生受信星座図（音源%d, PTR/合成なし）',idx_sd), ...
            'FontSize',12,'FontWeight','bold');
end

%% ================= ローカル関数 =================
function out = rls_dfe(x, symbols, Nf, Nb, lambda, delta, Ntrain)
% RLS判定帰還等化器（QPSK, T間隔）
%  x       : 入力シンボル列（合成後, α補正前）
%  symbols : 既知送信シンボル（トレーニング用）
%  返り値 out: 等化後シンボル列
    N = length(x);
    Ntap = Nf + Nb;
    w = zeros(Ntap,1); w(1)=1;          % FF先頭タップ=1で初期化（恒等から開始）
    P = (1/delta)*eye(Ntap);            % 逆相関行列初期化
    ff = zeros(Nf,1);                   % FFバッファ（受信シンボル）
    fb = zeros(Nb,1);                   % FBバッファ（判定シンボル）
    out = zeros(N,1);
    for n=1:N
        ff = [x(n); ff(1:end-1)];       % 新しい受信シンボルを先頭へ
        u  = [ff; -fb];                 % 入力ベクトル（FBは減算のため符号反転）
        y  = w' * u;                    % 等化出力  y = w^H u
        out(n) = y;
        d  = sign(real(y)) + 1i*sign(imag(y));   % QPSK硬判定
        if n <= Ntrain
            s_ref = symbols(n); fb_new = symbols(n);   % トレーニング
        else
            s_ref = d;          fb_new = d;            % 判定指向
        end
        % RLS更新
        Pu    = P*u;
        kappa = Pu/(lambda + u'*Pu);    % カルマンゲイン
        e     = s_ref - y;              % 事前誤差
        w     = w + kappa*conj(e);      % 係数更新
        P     = (P - kappa*(Pu'))/lambda;  % 逆相関行列更新
        fb    = [fb_new; fb(1:end-1)];  % 判定シンボルをFBバッファへ
    end
end

function ber = ber_qpsk(symbols, comp, bits)
% QPSK象限判定でBERを計算
    db = zeros(2*length(symbols),1);
    db(1:2:end)=real(comp)>0; db(2:2:end)=imag(comp)>0;
    ber = sum(bits~=db)/length(bits);
end
