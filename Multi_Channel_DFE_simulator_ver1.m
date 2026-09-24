%% ============================================================
%  Multi_Channel_DFE_simulator_ver1.m
%  PTR_DFE_simulator_ver3 をベースに「マルチチャネルDFE」を実装。
%
%  従来(ver3):  サブアレイPTR → MRCで1本に合成 → 単一チャネルDFE
%  本コード  :  サブアレイPTR → 少数ビーム(K本) → マルチチャネルDFE
%               （各ビームに個別FFフィルタ＋共有FB＋PLLをRLSで同時最適化）
%  → 合成と等化をMMSE基準で同時に行い、空間ダイバーシティをISIヌルに活用
%  参考: Stojanovic et al. 1994 (IEEE JOE), Sensors 2017 (TR-MC-DFE)
%
%  比較: 【SCM-MRC + 単一DFE】 vs 【マルチチャネルDFE】
%  評価: Monte Carlo平均。OSNR は誤差電力を平均してからdB化（線形平均）
% ============================================================
arr_file = 'JamTank_Arr_ver3.arr';    % ★使用中のARRに合わせる
[Arr, Pos] = read_arrivals_asc(arr_file);
Nrr=size(Arr,1); Nrz=size(Arr,2); Nsd=size(Arr,3); NUM_RX=Nrr*Nrz;
fprintf('ARR: 距離%d × 深度%d = %dch, 音源%d\n', Nrr,Nrz,NUM_RX,Nsd);

%% ==== 信号パラメータ ====
Rs=5000; Sps=16; Fs=Rs*Sps; numSymbols=2000;

%% ==== サブアレイ（ビーム）設定 ====
K_SUB = min(8, NUM_RX);          % ビーム数（=MC-DFEの入力チャネル数）
% 全受信機を K_SUB 個の連続グループに分割 → 各群をPTR合成して1ビーム
ch_sub = min(K_SUB, ceil((1:NUM_RX)'/ceil(NUM_RX/K_SUB)))';

%% ==== DFEパラメータ ====
% ★低次で良い：PTRが粗く集束済みなので、MC-DFEは残差の掃除に専念
Nf     = 8;      % 各ビームのフィードフォワードタップ数
Nb     = 15;     % 共有フィードバックタップ数
lambda = 0.999;  % RLS忘却係数
delta  = 0.01;   % RLS正則化
mu_pll = 0.01;   % PLL利得
Ntrain = 600;    % トレーニングシンボル数
SNR_dB = 20;

%% ==== Monte Carlo ====
Ntrials = 20;
fprintf('K_SUB(ビーム数)=%d, MC-DFE: Nf=%d/ch Nb=%d, Rs=%d, %d試行\n\n', K_SUB,Nf,Nb,Rs,Ntrials);

ch_rr=zeros(NUM_RX,1); ch_rz=zeros(NUM_RX,1); cc=0;
for irr=1:Nrr, for irz=1:Nrz, cc=cc+1; ch_rr(cc)=irr; ch_rz(cc)=irz; end, end

% 結果（線形平均用に信号/誤差電力を保存, test区間）
SCM_Ps=zeros(Nsd,Ntrials); SCM_Pe=zeros(Nsd,Ntrials); SCM_B=zeros(Nsd,Ntrials);
MC_Ps =zeros(Nsd,Ntrials); MC_Pe =zeros(Nsd,Ntrials); MC_B =zeros(Nsd,Ntrials);
SCM_syms=cell(Nsd,1); MC_syms=cell(Nsd,1);

for idx_sd=1:Nsd
    gmin=inf;
    for ch=1:NUM_RX
        d=real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay); d=d(d>0);
        if ~isempty(d), gmin=min(gmin,min(d)); end
    end

    for it=1:Ntrials
        %% 送信 QPSK
        bits=randi([0 1],2*numSymbols,1);
        symbols=(2*bits(1:2:end)-1)+1i*(2*bits(2:2:end)-1);
        txSignal=zeros(numSymbols*Sps,1);
        for i=1:numSymbols, txSignal((i-1)*Sps+1:i*Sps)=symbols(i); end

        ptr_buf_len=length(txSignal)*4;
        sub_ptr=zeros(ptr_buf_len,K_SUB); sub_ac=zeros(ptr_buf_len,K_SUB);
        for ch=1:NUM_RX
            amp=Arr(ch_rr(ch),ch_rz(ch),idx_sd).A; dl=real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay);
            v=(dl>0); amp=amp(v); dl=dl(v); dl=dl-gmin; v=(dl>=0); amp=amp(v); dl=dl(v);
            if isempty(dl), continue; end
            mds=round(max(dl)*Fs); rxBase=zeros(length(txSignal)+mds,1);
            for p=1:length(dl), ds=round(dl(p)*Fs); rxBase(ds+(1:length(txSignal)))=rxBase(ds+(1:length(txSignal)))+amp(p)*txSignal; end
            sp=mean(abs(rxBase).^2); np=sp/10^(SNR_dB/10);
            rx=rxBase+sqrt(np/2)*(randn(size(rxBase))+1i*randn(size(rxBase)));
            ds2=round(dl*Fs)+1; h=zeros(max(ds2),1); h(ds2)=amp;
            g=conj(flipud(h)); z=conv(rx,g); hac=conv(h,g);
            Lz=min(length(z),ptr_buf_len); La=min(length(hac),ptr_buf_len);
            si=ch_sub(ch);
            sub_ptr(1:Lz,si)=sub_ptr(1:Lz,si)+z(1:Lz); sub_ac(1:La,si)=sub_ac(1:La,si)+hac(1:La);
        end

        %% 各ビームを焦点でサンプリングし、複素利得で正規化（RLSの数値安定化に必須）
        X=zeros(numSymbols,K_SUB); valid=false(K_SUB,1); err_beam=zeros(K_SUB,1);
        for s=1:K_SUB
            if all(sub_ac(:,s)==0), continue; end
            [~,pks]=max(abs(sub_ac(:,s))); off=pks+round(Sps/2)-1;
            if off+Sps*(numSymbols-1)>ptr_buf_len, continue; end
            beam=sub_ptr(off:Sps:off+Sps*(numSymbols-1),s);
            a=(symbols'*beam)/(symbols'*symbols);          % 複素利得
            if abs(a)<eps || any(~isfinite(beam)), continue; end
            X(:,s)=beam/a;                                 % ★正規化（〜単位振幅シンボル）
            valid(s)=true; err_beam(s)=mean(abs(symbols-X(:,s)).^2);
        end
        X=X(:,valid); err_beam=err_beam(valid); Kv=size(X,2);
        if Kv==0, continue; end

        %% --- 参考: SCM-MRC合成 → 単一DFE ---（ビームは正規化済み）
        scm=zeros(numSymbols,1); wsum=0;
        for k=1:Kv
            w=1/err_beam(k); scm=scm+w*X(:,k); wsum=wsum+w;
        end
        scm=scm/wsum;
        scm_dfe=rls_dfe(scm,symbols,Nf*Kv,Nb,lambda,delta,Ntrain);   % 単一DFE(タップ総数を公平化)
        ti=(Ntrain+1):numSymbols;
        [~,SCM_B(idx_sd,it),cS,SCM_Ps(idx_sd,it),SCM_Pe(idx_sd,it)]=eval_region(symbols,scm_dfe,ti);

        %% --- 本命: マルチチャネルDFE ---
        mc_out=mc_dfe(X,symbols,Nf,Nb,lambda,delta,Ntrain,mu_pll);
        [~,MC_B(idx_sd,it),cM,MC_Ps(idx_sd,it),MC_Pe(idx_sd,it)]=eval_region(symbols,mc_out,ti);

        if it==Ntrials, SCM_syms{idx_sd}=cS; MC_syms{idx_sd}=cM; end
    end
end

%% ==== 集計（線形平均OSNR）====
linO=@(Ps,Pe) 10*log10(mean(Ps,2)./mean(Pe,2));
mB=@(B) mean(B,2);

disp('==========================================================================');
fprintf(' 比較: SCM-MRC+単一DFE  vs  マルチチャネルDFE（%d試行, Rs=%d, K=%d）\n', Ntrials,Rs,K_SUB);
fprintf(' OSNRは線形平均（誤差電力平均→dB化）, test区間(トレーニング除外)\n');
disp('==========================================================================');
for idx_sd=1:Nsd
    oS=linO(SCM_Ps(idx_sd,:),SCM_Pe(idx_sd,:)); oM=linO(MC_Ps(idx_sd,:),MC_Pe(idx_sd,:));
    fprintf('【音源%d】\n',idx_sd);
    fprintf('  SCM-MRC + 単一DFE : OSNR %6.2f dB   BER %.4f\n', oS, mB(SCM_B(idx_sd,:)));
    fprintf('  マルチチャネルDFE : OSNR %6.2f dB   BER %.4f   （改善 %+.2f dB）\n\n', oM, mB(MC_B(idx_sd,:)), oM-oS);
end

%% ==== 星座図（最終試行）====
qpsk=[1+1i,1-1i,-1+1i,-1-1i];
figure('Name','SCM+単一DFE vs マルチチャネルDFE','Position',[100 100 300*Nsd 560]);
for idx_sd=1:Nsd
    oS=linO(SCM_Ps(idx_sd,:),SCM_Pe(idx_sd,:)); oM=linO(MC_Ps(idx_sd,:),MC_Pe(idx_sd,:));
    subplot(2,Nsd,idx_sd);
    plot(real(SCM_syms{idx_sd}),imag(SCM_syms{idx_sd}),'.','Color',[0.1 0.6 0.9],'MarkerSize',3); hold on;
    plot(real(qpsk),imag(qpsk),'r+','MarkerSize',12,'LineWidth',2); axis equal; grid on;
    lim=max(3,ceil(max(abs([real(SCM_syms{idx_sd});imag(SCM_syms{idx_sd})])))); xlim([-lim lim]); ylim([-lim lim]);
    title(sprintf('SCM+単一DFE 音源%d\nOSNR=%.1fdB',idx_sd,oS));
    subplot(2,Nsd,Nsd+idx_sd);
    plot(real(MC_syms{idx_sd}),imag(MC_syms{idx_sd}),'.','Color',[0.85 0.3 0.1],'MarkerSize',3); hold on;
    plot(real(qpsk),imag(qpsk),'r+','MarkerSize',12,'LineWidth',2); axis equal; grid on;
    lim=max(3,ceil(max(abs([real(MC_syms{idx_sd});imag(MC_syms{idx_sd})])))); xlim([-lim lim]); ylim([-lim lim]);
    title(sprintf('マルチチャネルDFE 音源%d\nOSNR=%.1fdB',idx_sd,oM));
end
sgtitle(sprintf('マルチチャネルDFEの効果（K=%d, Nf=%d/ch, Nb=%d, Rs=%d）',K_SUB,Nf,Nb,Rs),'FontWeight','bold');

%% ================= ローカル関数 =================
function out = mc_dfe(X, symbols, Nf, Nb, lambda, delta, Ntrain, mu_pll)
% マルチチャネルDFE（各チャネル個別FF＋共有FB＋共有PLL, RLS）
%  X: numSymbols × K（K本の入力ビーム, α補正前）
    [N,K]=size(X); Ntap=K*Nf+Nb;
    w=zeros(Ntap,1);
    for k=1:K, w((k-1)*Nf+1)=1/K; end     % 初期値：各chの先頭FFを1/K（初期出力≈平均）
    P=(1/delta)*eye(Ntap);
    ff=zeros(Nf,K); fb=zeros(Nb,1); theta=0; out=zeros(N,1);
    for n=1:N
        ff=[X(n,:); ff(1:end-1,:)];        % 各チャネルのFFバッファを更新
        ffvec=reshape(ff,[],1)*exp(-1i*theta);   % FF部（PLLで位相回転）K*Nf×1
        u=[ffvec; -fb];
        y=w'*u; out(n)=y;
        d=sign(real(y))+1i*sign(imag(y));
        if n<=Ntrain, s_ref=symbols(n); fb_new=symbols(n); else, s_ref=d; fb_new=d; end
        e=s_ref-y;
        Pu=P*u; kap=Pu/(lambda+u'*Pu);
        w=w+kap*conj(e); P=(P-kap*(Pu'))/lambda;
        phi=imag(y*conj(s_ref)); theta=theta+mu_pll*phi;
        fb=[fb_new; fb(1:end-1)];
    end
end

function out = rls_dfe(x, symbols, Nf, Nb, lambda, delta, Ntrain)
% 単一チャネルRLS-DFE（比較基準用）
    N=length(x); Ntap=Nf+Nb; w=zeros(Ntap,1); w(1)=1; P=(1/delta)*eye(Ntap);
    ff=zeros(Nf,1); fb=zeros(Nb,1); out=zeros(N,1);
    for n=1:N
        ff=[x(n); ff(1:end-1)]; u=[ff; -fb]; y=w'*u; out(n)=y;
        d=sign(real(y))+1i*sign(imag(y));
        if n<=Ntrain, s_ref=symbols(n); fb_new=symbols(n); else, s_ref=d; fb_new=d; end
        e=s_ref-y; Pu=P*u; kap=Pu/(lambda+u'*Pu);
        w=w+kap*conj(e); P=(P-kap*(Pu'))/lambda; fb=[fb_new; fb(1:end-1)];
    end
end

function [osnr, ber, comp, psig, perr] = eval_region(symbols, raw, idx)
    s=symbols(idx); r=raw(idx); a=(s'*r)/(s'*s); c=r/a;
    psig=mean(abs(s).^2); perr=mean(abs(s-c).^2); osnr=10*log10(psig/perr);
    ber=mean([real(c)>0,imag(c)>0]~=[real(s)>0,imag(s)>0],'all');
    comp=raw/a;
end
