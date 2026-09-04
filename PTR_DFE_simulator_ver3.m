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
% ★線形平均OSNR用：各試行の信号電力Ps・誤差電力Pe（線形）を保存
PTR_Ps =zeros(Nsd,Ntrials); PTR_Pe =zeros(Nsd,Ntrials);
SCM_Ps =zeros(Nsd,Ntrials); SCM_Pe =zeros(Nsd,Ntrials);   % 全区間
SCM_Pst=zeros(Nsd,Ntrials); SCM_Pet=zeros(Nsd,Ntrials);   % test区間
DFE_Ps =zeros(Nsd,Ntrials); DFE_Pe =zeros(Nsd,Ntrials);   % 全区間
DFE_Pst=zeros(Nsd,Ntrials); DFE_Pet=zeros(Nsd,Ntrials);   % test区間
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
        [PTR_O(idx_sd,it),PTR_B(idx_sd,it),~,PTR_Ps(idx_sd,it),PTR_Pe(idx_sd,it)] = eval_region(symbols, ptrSym, all_idx);

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
        [SCM_O(idx_sd,it), SCM_B(idx_sd,it), ~,SCM_Ps(idx_sd,it), SCM_Pe(idx_sd,it)]  = eval_region(symbols, scm, all_idx);
        [SCM_Ot(idx_sd,it),SCM_Bt(idx_sd,it),~,SCM_Pst(idx_sd,it),SCM_Pet(idx_sd,it)] = eval_region(symbols, scm, test_idx);

        %% SCM-PTR + RLS-DFE
        dfe_out = rls_dfe(scm, symbols, Nf, Nb, lambda, delta, Ntrain);
        [DFE_O(idx_sd,it), DFE_B(idx_sd,it), ~,DFE_Ps(idx_sd,it), DFE_Pe(idx_sd,it)]  = eval_region(symbols, dfe_out, all_idx);
        [DFE_Ot(idx_sd,it),DFE_Bt(idx_sd,it),~,DFE_Pst(idx_sd,it),DFE_Pet(idx_sd,it)] = eval_region(symbols, dfe_out, test_idx);

        %% 星座図用（最終試行のみ保存）
        if it==Ntrials
            [~,~,SCM_syms{idx_sd}] = eval_region(symbols, scm, all_idx);
            [~,~,DFE_syms{idx_sd}] = eval_region(symbols, dfe_out, all_idx);
        end
    end % ◀ 試行ループ
end % ◀ 音源ループ

%% 集計関数（試行方向 dim=2）
m =@(X) mean(X,2); s=@(X) std(X,0,2); md=@(X) median(X,2);
% ★線形平均OSNR：各試行の誤差電力を平均してからdB化（Jensen上振れを排除）
%   OSNR_lin = 10*log10( mean_trials(Psig) / mean_trials(Perr) )
linO=@(Ps,Pe) 10*log10(mean(Ps,2)./mean(Pe,2));

%% 結果出力
disp('==========================================================================');
fprintf(' Monte Carlo（%d試行）  NUM_SUB=%d, Rs=%dbaud, DFE:Nf=%d Nb=%d Ntrain=%d\n', ...
        Ntrials, NUM_SUB, Rs, Nf, Nb, Ntrain);
fprintf(' ★OSNRは「誤差電力を平均してからdB化」した線形平均値を主表示（保守的・物理的）\n');
fprintf(' 参考: dB平均(Jensen上振れ) と 中央値 も併記。 [test]=判定指向区間\n');
disp('==========================================================================');
for idx_sd=1:Nsd
    fprintf('【音源 %d】\n', idx_sd);
    fprintf('  通常PTR             : OSNR(線形平均) %6.2f dB   [dB平均 %5.2f / 中央値 %5.2f]  BER %.4f\n', ...
        linO(PTR_Ps(idx_sd,:),PTR_Pe(idx_sd,:)), m(PTR_O(idx_sd,:)), md(PTR_O(idx_sd,:)), m(PTR_B(idx_sd,:)));
    fprintf('  SCM-PTR      [全]   : OSNR(線形平均) %6.2f dB   [dB平均 %5.2f / 中央値 %5.2f]  BER %.4f\n', ...
        linO(SCM_Ps(idx_sd,:),SCM_Pe(idx_sd,:)), m(SCM_O(idx_sd,:)), md(SCM_O(idx_sd,:)), m(SCM_B(idx_sd,:)));
    fprintf('  SCM-PTR      [test] : OSNR(線形平均) %6.2f dB   [dB平均 %5.2f / 中央値 %5.2f]  BER %.4f\n', ...
        linO(SCM_Pst(idx_sd,:),SCM_Pet(idx_sd,:)), m(SCM_Ot(idx_sd,:)), md(SCM_Ot(idx_sd,:)), m(SCM_Bt(idx_sd,:)));
    fprintf('  SCM-PTR+DFE  [全]   : OSNR(線形平均) %6.2f dB   [dB平均 %5.2f / 中央値 %5.2f]  BER %.4f\n', ...
        linO(DFE_Ps(idx_sd,:),DFE_Pe(idx_sd,:)), m(DFE_O(idx_sd,:)), md(DFE_O(idx_sd,:)), m(DFE_B(idx_sd,:)));
    fprintf('  SCM-PTR+DFE  [test] : OSNR(線形平均) %6.2f dB   [dB平均 %5.2f / 中央値 %5.2f]  BER %.4f\n', ...
        linO(DFE_Pst(idx_sd,:),DFE_Pet(idx_sd,:)), m(DFE_Ot(idx_sd,:)), md(DFE_Ot(idx_sd,:)), m(DFE_Bt(idx_sd,:)));
    dfe_gain = linO(DFE_Pst(idx_sd,:),DFE_Pet(idx_sd,:)) - linO(SCM_Pst(idx_sd,:),SCM_Pet(idx_sd,:));
    fprintf('  ── DFE改善量(線形平均, test基準): ΔOSNR = %+.2f dB\n\n', dfe_gain);
end

fprintf('（線形平均 vs dB平均の差 = Jensen上振れ。差が大きいほど試行ごとのばらつき大）\n');

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
    title(sprintf('SCM-PTR 音源%d\ntest OSNR(線形)=%.1fdB (±%.1f)',idx_sd,linO(SCM_Pst(idx_sd,:),SCM_Pet(idx_sd,:)),s(SCM_Ot(idx_sd,:))));

    subplot(2,Nsd,Nsd+idx_sd);
    plot(real(dfeC),imag(dfeC),'.','Color',[0.85 0.3 0.1],'MarkerSize',3); hold on;
    plot(real(qpsk_ideal),imag(qpsk_ideal),'r+','MarkerSize',12,'LineWidth',2); hold off;
    axis equal; grid on; lim=max(3,ceil(max(abs([real(dfeC);imag(dfeC)]))));
    xlim([-lim lim]); ylim([-lim lim]); xlabel('I'); ylabel('Q');
    title(sprintf('SCM+DFE 音源%d\ntest OSNR(線形)=%.1fdB (±%.1f)',idx_sd,linO(DFE_Pst(idx_sd,:),DFE_Pet(idx_sd,:)),s(DFE_Ot(idx_sd,:))));
end
sgtitle(sprintf('DFE等化（%d試行平均, test区間OSNR） Nf=%d Nb=%d',Ntrials,Nf,Nb),'FontSize',13,'FontWeight','bold');

%% ================= ローカル関数 =================
function [osnr, ber, comp, psig, perr] = eval_region(symbols, raw, idx)
    s = symbols(idx); r = raw(idx);
    a = (s'*r)/(s'*s);
    c = r / a;
    psig = mean(abs(s).^2);              % 信号電力（線形）
    perr = mean(abs(s-c).^2);            % 残留誤差電力（線形）
    osnr = 10*log10(psig/perr);          % この試行のOSNR(dB)
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
