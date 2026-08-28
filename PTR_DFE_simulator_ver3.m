%% ============================================================
%  PTR_DFE_simulator_ver3.m
%  ver2に Monte Carlo 平均を追加。
%  ・単発だと送信ビット・雑音の乱数でOSNR/BERが数dB揺れるため、
%    同一チャネル(Arr)で Ntrials 回繰り返し、平均±標準偏差で評価する。
%  ・各方式を [全区間] と [判定指向区間(test)] で評価（過学習チェック維持）
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

%% サブアレイ分割
if Nrr >= 2
    NUM_SUB = Nrr; SUB_SIZE = Nrz; ch_sub = ch_rr;
else
    NUM_SUB = 2; SUB_SIZE = ceil(NUM_RX/NUM_SUB);
    ch_sub  = min(NUM_SUB, ceil((1:NUM_RX)'/SUB_SIZE));
end
fprintf('サブアレイ数 NUM_SUB=%d（距離%d, 深度%d）\n', NUM_SUB, Nrr, Nrz);

%% ★DFEパラメータ
Nf     = 8;
Nb     = 6;
lambda = 0.999;
delta  = 0.01;
Ntrain = 400;

%% ★Monte Carlo 設定
Ntrials = 30;     % 試行回数（多いほど平均が安定。目安20〜50）
fprintf('Monte Carlo: %d 回試行して平均±標準偏差を算出\n\n', Ntrials);

%% 結果格納 [Nsd × Ntrials]
PTR_O=zeros(Nsd,Ntrials);  PTR_B=zeros(Nsd,Ntrials);
SCM_O=zeros(Nsd,Ntrials);  SCM_B=zeros(Nsd,Ntrials);   % 全区間
SCM_Ot=zeros(Nsd,Ntrials); SCM_Bt=zeros(Nsd,Ntrials);  % test区間
DFE_O=zeros(Nsd,Ntrials);  DFE_B=zeros(Nsd,Ntrials);   % 全区間
DFE_Ot=zeros(Nsd,Ntrials); DFE_Bt=zeros(Nsd,Ntrials);  % test区間
SCM_syms=cell(Nsd,1); DFE_syms=cell(Nsd,1);            % 星座図（最終試行）

Rs=500; Sps=16; Fs=Rs*Sps; numSymbols=2000;

for idx_sd = 1:Nsd

    %% チャネル依存量は試行間で不変 → 事前に1回だけ算出
    gmin=inf;
    for ch=1:NUM_RX
        d=real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay); d=d(d>0);
        if ~isempty(d), gmin=min(gmin,min(d)); end
    end

    for it = 1:Ntrials

        %% 送信信号（試行ごとに乱数生成）
        bits=randi([0 1],2*numSymbols,1);
        symbols=(2*bits(1:2:end)-1)+1i*(2*bits(2:2:end)-1);
        txSignal=zeros(numSymbols*Sps,1);
        for i=1:numSymbols, txSignal((i-1)*Sps+1:i*Sps)=symbols(i); end

        ptr_buf_len=length(txSignal)*4;
        total_ptr=zeros(ptr_buf_len,1); total_ac=zeros(ptr_buf_len,1);
        sub_ptr=zeros(ptr_buf_len,NUM_SUB); sub_ac=zeros(ptr_buf_len,NUM_SUB);

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
            SNR_dB=20; sp=mean(abs(rxBase).^2); np=sp/10^(SNR_dB/10);
            rx=rxBase+sqrt(np/2)*(randn(size(rxBase))+1i*randn(size(rxBase)));

            ds2=round(dl*Fs)+1; h=zeros(max(ds2),1); h(ds2)=amp;
            g=conj(flipud(h)); z=conv(rx,g); hac=conv(h,g);
            Lz=min(length(z),ptr_buf_len); La=min(length(hac),ptr_buf_len);
            total_ptr(1:Lz)=total_ptr(1:Lz)+z(1:Lz);
            total_ac(1:La) =total_ac(1:La)+hac(1:La);
            si=ch_sub(ch);
            sub_ptr(1:Lz,si)=sub_ptr(1:Lz,si)+z(1:Lz);
            sub_ac(1:La,si) =sub_ac(1:La,si)+hac(1:La);
        end

        all_idx  = 1:numSymbols;
        test_idx = (Ntrain+1):numSymbols;

        %% 通常PTR
        [~,pk]=max(abs(total_ac)); off=pk+round(Sps/2)-1;
        ptrSym=total_ptr(off:Sps:off+Sps*(numSymbols-1));
        [PTR_O(idx_sd,it),PTR_B(idx_sd,it)] = eval_region(symbols, ptrSym, all_idx);

        %% SCM-PTR
        scm=zeros(numSymbols,1); wsum=0;
        for s=1:NUM_SUB
            [~,pks]=max(abs(sub_ac(:,s))); offs=pks+round(Sps/2)-1;
            if offs+Sps*(numSymbols-1)>ptr_buf_len, continue; end
            ss=sub_ptr(offs:Sps:offs+Sps*(numSymbols-1),s);
            as=(symbols'*ss)/(symbols'*symbols); sc=ss/as;
            es=mean(abs(symbols-sc).^2); ws=1/es;
            scm=scm+ws*sc; wsum=wsum+ws;
        end
        scm=scm/wsum;
        [SCM_O(idx_sd,it), SCM_B(idx_sd,it)]  = eval_region(symbols, scm, all_idx);
        [SCM_Ot(idx_sd,it),SCM_Bt(idx_sd,it)] = eval_region(symbols, scm, test_idx);

        %% SCM-PTR + RLS-DFE
        dfe_out = rls_dfe(scm, symbols, Nf, Nb, lambda, delta, Ntrain);
        [DFE_O(idx_sd,it), DFE_B(idx_sd,it)]  = eval_region(symbols, dfe_out, all_idx);
        [DFE_Ot(idx_sd,it),DFE_Bt(idx_sd,it)] = eval_region(symbols, dfe_out, test_idx);

        %% 星座図用（最終試行のみ保存）
        if it==Ntrials
            [~,~,SCM_syms{idx_sd}] = eval_region(symbols, scm, all_idx);
            [~,~,DFE_syms{idx_sd}] = eval_region(symbols, dfe_out, all_idx);
        end
    end % ◀ 試行ループ
end % ◀ 音源ループ

%% 平均±標準偏差の算出（試行方向 dim=2）
m=@(X) mean(X,2); s=@(X) std(X,0,2);

%% 結果出力
disp('==========================================================================');
fprintf(' Monte Carlo 平均（%d試行）  NUM_SUB=%d, Rs=%dbaud, DFE:Nf=%d Nb=%d Ntrain=%d\n', ...
        Ntrials, NUM_SUB, Rs, Nf, Nb, Ntrain);
fprintf(' 表記: 平均 ± 標準偏差   [全]=全区間  [test]=判定指向区間\n');
disp('==========================================================================');
for idx_sd=1:Nsd
    fprintf('【音源 %d】\n', idx_sd);
    fprintf('  通常PTR             : OSNR %5.2f ± %4.2f dB   BER %.4f ± %.4f\n', ...
        m(PTR_O(idx_sd,:)), s(PTR_O(idx_sd,:)), m(PTR_B(idx_sd,:)), s(PTR_B(idx_sd,:)));
    fprintf('  SCM-PTR      [全]   : OSNR %5.2f ± %4.2f dB   BER %.4f ± %.4f\n', ...
        m(SCM_O(idx_sd,:)), s(SCM_O(idx_sd,:)), m(SCM_B(idx_sd,:)), s(SCM_B(idx_sd,:)));
    fprintf('  SCM-PTR      [test] : OSNR %5.2f ± %4.2f dB   BER %.4f ± %.4f\n', ...
        m(SCM_Ot(idx_sd,:)), s(SCM_Ot(idx_sd,:)), m(SCM_Bt(idx_sd,:)), s(SCM_Bt(idx_sd,:)));
    fprintf('  SCM-PTR+DFE  [全]   : OSNR %5.2f ± %4.2f dB   BER %.4f ± %.4f\n', ...
        m(DFE_O(idx_sd,:)), s(DFE_O(idx_sd,:)), m(DFE_B(idx_sd,:)), s(DFE_B(idx_sd,:)));
    fprintf('  SCM-PTR+DFE  [test] : OSNR %5.2f ± %4.2f dB   BER %.4f ± %.4f\n', ...
        m(DFE_Ot(idx_sd,:)), s(DFE_Ot(idx_sd,:)), m(DFE_Bt(idx_sd,:)), s(DFE_Bt(idx_sd,:)));
    fprintf('  ── DFE改善量(test基準): ΔOSNR = %+.2f dB\n', m(DFE_Ot(idx_sd,:))-m(SCM_Ot(idx_sd,:)));
    fprintf('  ── 過学習チェック(DFE 全−test): %+.2f dB\n\n', m(DFE_O(idx_sd,:))-m(DFE_Ot(idx_sd,:)));
end

%% 収束の目安：標準偏差 / √Ntrials ≒ 平均値の標準誤差
fprintf('（平均値の標準誤差 ≒ 標準偏差/√%d。誤差が十分小さくなければ Ntrials を増やす）\n', Ntrials);

%% 星座図（最終試行, test区間の平均OSNRをタイトルに併記）
qpsk_ideal=[1+1i,1-1i,-1+1i,-1-1i];
figure('Name','星座図: SCM vs SCM+DFE (最終試行)','Position',[100 100 300*Nsd 560]);
for idx_sd=1:Nsd
    scmC=SCM_syms{idx_sd}; dfeC=DFE_syms{idx_sd};
    subplot(2,Nsd,idx_sd);
    plot(real(scmC),imag(scmC),'.','Color',[0.1 0.8 0.3],'MarkerSize',3); hold on;
    plot(real(qpsk_ideal),imag(qpsk_ideal),'r+','MarkerSize',12,'LineWidth',2); hold off;
    axis equal; grid on; lim=max(3,ceil(max(abs([real(scmC);imag(scmC)]))));
    xlim([-lim lim]); ylim([-lim lim]); xlabel('I'); ylabel('Q');
    title(sprintf('SCM-PTR 音源%d\ntest OSNR=%.1f±%.1fdB',idx_sd,m(SCM_Ot(idx_sd,:)),s(SCM_Ot(idx_sd,:))));

    subplot(2,Nsd,Nsd+idx_sd);
    plot(real(dfeC),imag(dfeC),'.','Color',[0.85 0.3 0.1],'MarkerSize',3); hold on;
    plot(real(qpsk_ideal),imag(qpsk_ideal),'r+','MarkerSize',12,'LineWidth',2); hold off;
    axis equal; grid on; lim=max(3,ceil(max(abs([real(dfeC);imag(dfeC)]))));
    xlim([-lim lim]); ylim([-lim lim]); xlabel('I'); ylabel('Q');
    title(sprintf('SCM+DFE 音源%d\ntest OSNR=%.1f±%.1fdB',idx_sd,m(DFE_Ot(idx_sd,:)),s(DFE_Ot(idx_sd,:))));
end
sgtitle(sprintf('DFE等化（%d試行平均, test区間OSNR） Nf=%d Nb=%d',Ntrials,Nf,Nb),'FontSize',13,'FontWeight','bold');

%% ================= ローカル関数 =================
function [osnr, ber, comp] = eval_region(symbols, raw, idx)
    s = symbols(idx); r = raw(idx);
    a = (s'*r)/(s'*s);
    c = r / a;
    osnr = 10*log10(mean(abs(s).^2)/mean(abs(s-c).^2));
    dec = [real(c)>0, imag(c)>0];
    ref = [real(s)>0, imag(s)>0];
    ber = mean(dec(:) ~= ref(:));
    comp = raw / a;
end

function out = rls_dfe(x, symbols, Nf, Nb, lambda, delta, Ntrain)
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
