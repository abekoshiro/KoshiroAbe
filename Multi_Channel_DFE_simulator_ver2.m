%% ============================================================
%  Multi_Channel_DFE_simulator_ver2.m
%  ver1 に【RRCパルス整形（ロールオフフィルタ）】を導入した版。
%
%  専門家指摘への対応:
%   ⓪ 送信を矩形→RRCベースバンド整形（帯域制限, ナイキストでISIなし）
%   ② 受信にも整合RRCを適用（帯域外雑音除去 → SNRは帯域内で定義）
%   ① PTR合成の整列を共通長CIRに統一（遅延依存の焦点ズレを除去）
%     ※遅延の丸め誤差は Sps を上げると低減（本版は時間領域, 将来は周波数領域/sinc補間）
%
%  構成: サブアレイPTR → 少数ビーム → マルチチャネルDFE
%  比較: SCM-MRC+単一DFE vs マルチチャネルDFE（Monte Carlo平均, 線形平均OSNR）
%  依存: rrc_time_domain.m（同ディレクトリ）
% ============================================================
arr_file = 'JamTank_Arr_ver3.arr';    % ★使用中のARRに合わせる
[Arr, Pos] = read_arrivals_asc(arr_file);
Nrr=size(Arr,1); Nrz=size(Arr,2); Nsd=size(Arr,3); NUM_RX=Nrr*Nrz;
fprintf('ARR: 距離%d × 深度%d = %dch, 音源%d\n', Nrr,Nrz,NUM_RX,Nsd);

%% ==== 信号パラメータ ====
Rs=5000; Sps=16; Fs=Rs*Sps; numSymbols=2000;   % BW=Rs, Fs=80kHz
alpha  = 0.25;    % ★ロールオフ率 (0<α<1)。帯域=Rs*(1+α)
Ncycle = 20;      % RRCのシンボルスパン（±Ncycle。20以上推奨）
SNR_dB = 20;      % 帯域内SNR（受信整合RRC後に有効）

%% ==== サブアレイ（ビーム）設定 ====
K_SUB = min(8, NUM_RX);
ch_sub = min(K_SUB, ceil((1:NUM_RX)'/ceil(NUM_RX/K_SUB)))';

%% ==== DFEパラメータ ====
Nf=8; Nb=15; lambda=0.999; delta=0.01; mu_pll=0.01; Ntrain=600;

%% ==== Monte Carlo ====
Ntrials = 10;    % RRC+長い畳込みで重いので控えめ（必要なら増やす）
fprintf('RRC: α=%.2f, Ncycle=%d | K_SUB=%d, MC-DFE Nf=%d/ch Nb=%d | Rs=%d, %d試行\n\n', ...
        alpha,Ncycle,K_SUB,Nf,Nb,Rs,Ntrials);

ch_rr=zeros(NUM_RX,1); ch_rz=zeros(NUM_RX,1); cc=0;
for irr=1:Nrr, for irz=1:Nrz, cc=cc+1; ch_rr(cc)=irr; ch_rz(cc)=irz; end, end

%% ==== RRCフィルタ（送信・受信整合で共用）====
rrc = rrc_time_domain(Rs, Ncycle, alpha, Fs);   % 提供関数
% ★特異点 t=±Ts/(4α) で浮動小数点判定が外れInf/NaNになる場合を閉形式値で補修
bad = ~isfinite(rrc);
if any(bad)
    val = alpha/sqrt(2)*((1+2/pi)*sin(pi/(4*alpha))+(1-2/pi)*cos(pi/(4*alpha)));
    rrc(bad) = val;
    fprintf('  （RRC特異点 %d 点を閉形式値 %.4f で補修）\n', nnz(bad), val);
end
rrc = rrc(:) / sqrt(sum(abs(rrc).^2));          % エネルギー正規化
gd  = Ncycle*Sps;                               % RRC群遅延[サンプル]

% 結果（test区間, 線形平均用）
SCM_Ps=zeros(Nsd,Ntrials); SCM_Pe=zeros(Nsd,Ntrials); SCM_B=zeros(Nsd,Ntrials);
MC_Ps =zeros(Nsd,Ntrials); MC_Pe =zeros(Nsd,Ntrials); MC_B =zeros(Nsd,Ntrials);
SCM_syms=cell(Nsd,1); MC_syms=cell(Nsd,1);

for idx_sd=1:Nsd
    %% 共通最小遅延・共通CIR長（★PTR焦点の整列に必須）
    gmin=inf; Lmax=0;
    for ch=1:NUM_RX
        d=real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay); d=d(d>0);
        if ~isempty(d), gmin=min(gmin,min(d)); end
    end
    for ch=1:NUM_RX
        d=real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay); d=d(d>0)-gmin; d=d(d>=0);
        if ~isempty(d), Lmax=max(Lmax, round(max(d)*Fs)+1); end
    end

    for it=1:Ntrials
        %% 送信 QPSK → RRCベースバンド整形
        bits=randi([0 1],2*numSymbols,1);
        symbols=(2*bits(1:2:end)-1)+1i*(2*bits(2:2:end)-1);
        up=zeros(numSymbols*Sps,1); up(1:Sps:end)=symbols;    % シンボルをSps間隔に配置
        tx=conv(up,rrc);                                       % RRC整形（ベースバンド）
        Ltx=length(tx);

        ptr_buf_len=Ltx+2*Lmax+4*length(rrc);
        sub_ptr=zeros(ptr_buf_len,K_SUB); sub_ac=zeros(ptr_buf_len,K_SUB);

        for ch=1:NUM_RX
            amp=Arr(ch_rr(ch),ch_rz(ch),idx_sd).A; dl=real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay);
            v=(dl>0); amp=amp(v); dl=dl(v); dl=dl-gmin; v=(dl>=0); amp=amp(v); dl=dl(v);
            if isempty(dl), continue; end

            %% チャネル適用：各パスを遅延・複素振幅で重畳
            mds=round(max(dl)*Fs); rxBase=zeros(Ltx+mds,1);
            for p=1:length(dl)
                ds=round(dl(p)*Fs);
                rxBase(ds+(1:Ltx))=rxBase(ds+(1:Ltx))+amp(p)*tx;
            end
            %% 帯域内SNRで雑音付加 → 受信整合RRC（帯域外雑音除去）
            sp=mean(abs(rxBase).^2); np=sp/10^(SNR_dB/10);
            rxn=rxBase+sqrt(np/2)*(randn(size(rxBase))+1i*randn(size(rxBase)));
            rxmf=conv(rxn,rrc);                               % 受信整合RRC

            %% PTR：共通長Lmaxの h で整合フィルタ（★焦点整列）
            h=zeros(Lmax,1); ds2=round(dl*Fs)+1; ds2=ds2(ds2>=1&ds2<=Lmax);
            for p=1:numel(ds2), h(ds2(p))=h(ds2(p))+amp(p); end
            g=conj(flipud(h));
            z=conv(rxmf,g); hac=conv(h,g);
            Lz=min(length(z),ptr_buf_len); La=min(length(hac),ptr_buf_len);
            si=ch_sub(ch);
            sub_ptr(1:Lz,si)=sub_ptr(1:Lz,si)+z(1:Lz);
            sub_ac(1:La,si) =sub_ac(1:La,si)+hac(1:La);
        end

        %% 各ビームを焦点でサンプリング（整合後の最良タイミングを探索）し正規化
        X=zeros(numSymbols,K_SUB); valid=false(K_SUB,1); err_beam=zeros(K_SUB,1);
        for s=1:K_SUB
            if all(sub_ac(:,s)==0), continue; end
            [~,pk]=max(abs(sub_ac(:,s)));
            % 焦点はq関数ピークからRRC群遅延(送受×2)ぶん後方 → 広域を粗→微で探索
            beam=find_sample(sub_ptr(:,s),pk,Sps,numSymbols,symbols,3*gd+4*Sps);
            if isempty(beam), continue; end
            a=(symbols'*beam)/(symbols'*symbols);
            if abs(a)<eps||any(~isfinite(beam)), continue; end
            X(:,s)=beam/a; valid(s)=true; err_beam(s)=mean(abs(symbols-X(:,s)).^2);
        end
        X=X(:,valid); err_beam=err_beam(valid); Kv=size(X,2);
        if Kv==0, continue; end
        ti=(Ntrain+1):numSymbols;

        %% 参考: SCM-MRC → 単一DFE（従来型）
        scm=zeros(numSymbols,1); wsum=0;
        for k=1:Kv, w=1/err_beam(k); scm=scm+w*X(:,k); wsum=wsum+w; end
        scm=scm/wsum;
        %  ※単一ストリームにNf*Kvタップは過大で発散するため、等倍Nf/Nbを使用
        scm_dfe=rls_dfe(scm,symbols,Nf,Nb,lambda,delta,Ntrain);
        if all(isfinite(scm_dfe))
            [~,SCM_B(idx_sd,it),cS,SCM_Ps(idx_sd,it),SCM_Pe(idx_sd,it)]=eval_region(symbols,scm_dfe,ti);
        else
            SCM_Ps(idx_sd,it)=1; SCM_Pe(idx_sd,it)=1; SCM_B(idx_sd,it)=0.5; cS=zeros(numSymbols,1); % 発散→OSNR0dB相当
        end

        %% 本命: マルチチャネルDFE
        mc_out=mc_dfe(X,symbols,Nf,Nb,lambda,delta,Ntrain,mu_pll);
        [~,MC_B(idx_sd,it),cM,MC_Ps(idx_sd,it),MC_Pe(idx_sd,it)]=eval_region(symbols,mc_out,ti);

        if it==Ntrials, SCM_syms{idx_sd}=cS; MC_syms{idx_sd}=cM; end
    end
end

%% ==== 集計（線形平均OSNR）====
linO=@(Ps,Pe) 10*log10(mean(Ps,2)./mean(Pe,2)); mB=@(B) mean(B,2);
disp('==========================================================================');
fprintf(' RRC導入版: SCM-MRC+単一DFE vs マルチチャネルDFE（%d試行, Rs=%d, α=%.2f, K=%d）\n', Ntrials,Rs,alpha,K_SUB);
fprintf(' OSNRは線形平均（誤差電力平均→dB化）, test区間\n');
disp('==========================================================================');
for idx_sd=1:Nsd
    oS=linO(SCM_Ps(idx_sd,:),SCM_Pe(idx_sd,:)); oM=linO(MC_Ps(idx_sd,:),MC_Pe(idx_sd,:));
    fprintf('【音源%d】\n',idx_sd);
    fprintf('  SCM-MRC + 単一DFE : OSNR %6.2f dB   BER %.4f\n', oS, mB(SCM_B(idx_sd,:)));
    fprintf('  マルチチャネルDFE : OSNR %6.2f dB   BER %.4f   （改善 %+.2f dB）\n\n', oM, mB(MC_B(idx_sd,:)), oM-oS);
end

%% ==== 星座図 ====
qpsk=[1+1i,1-1i,-1+1i,-1-1i];
figure('Name','RRC導入: SCM+単一DFE vs マルチチャネルDFE','Position',[100 100 300*Nsd 560]);
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
sgtitle(sprintf('RRC導入版（α=%.2f, K=%d, Nf=%d/ch, Nb=%d, Rs=%d）',alpha,K_SUB,Nf,Nb,Rs),'FontWeight','bold');

%% ================= ローカル関数 =================
function sym = find_sample(sig, start, Sps, numSymbols, symbols, maxoff)
% 焦点を粗→微の2段で探索し、送信シンボルと最も相関する点でサンプリング
%  群遅延（RRC送受）により焦点が start から広く後方にずれるため広域探索
    % --- 粗探索（Sps刻みで symbol整合の当たりを付ける）---
    bestd=0; best=-inf;
    for d=0:Sps:maxoff
        o=start+d; idx=o+(0:numSymbols-1)*Sps;
        if idx(1)<1||idx(end)>length(sig), continue; end
        v=osnr_of(sig(idx),symbols); if v>best, best=v; bestd=d; end
    end
    % --- 微探索（±Sps を1刻み）---
    sym=[]; best=-inf;
    for d=bestd-Sps:bestd+Sps
        o=start+d; idx=o+(0:numSymbols-1)*Sps;
        if idx(1)<1||idx(end)>length(sig), continue; end
        cand=sig(idx); v=osnr_of(cand,symbols); if v>best, best=v; sym=cand; end
    end
end

function v = osnr_of(cand, symbols)
    g0=(symbols'*cand)/(symbols'*symbols);
    if abs(g0)<eps, v=-inf; return; end
    cc=cand/g0; v=sum(abs(symbols).^2)/max(sum(abs(symbols-cc).^2),eps);
end

function out = mc_dfe(X, symbols, Nf, Nb, lambda, delta, Ntrain, mu_pll)
% マルチチャネルDFE（各ch個別FF＋共有FB＋共有PLL, RLS）
    [N,K]=size(X); Ntap=K*Nf+Nb; w=zeros(Ntap,1);
    for k=1:K, w((k-1)*Nf+1)=1/K; end
    P=(1/delta)*eye(Ntap); ff=zeros(Nf,K); fb=zeros(Nb,1); theta=0; out=zeros(N,1);
    for n=1:N
        ff=[X(n,:); ff(1:end-1,:)];
        ffvec=reshape(ff,[],1)*exp(-1i*theta);
        u=[ffvec; -fb]; y=w'*u; out(n)=y;
        d=sign(real(y))+1i*sign(imag(y));
        if n<=Ntrain, s_ref=symbols(n); fb_new=symbols(n); else, s_ref=d; fb_new=d; end
        e=s_ref-y; Pu=P*u; kap=Pu/(lambda+u'*Pu);
        w=w+kap*conj(e); P=(P-kap*(Pu'))/lambda;
        phi=imag(y*conj(s_ref)); theta=theta+mu_pll*phi;
        fb=[fb_new; fb(1:end-1)];
    end
end

function out = rls_dfe(x, symbols, Nf, Nb, lambda, delta, Ntrain)
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
    ber=mean([real(c)>0,imag(c)>0]~=[real(s)>0,imag(s)>0],'all'); comp=raw/a;
end
