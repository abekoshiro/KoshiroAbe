%% ============================================================
%  Kida2017_PTR_DFE_replication.m   (検証反映版 rev2)
%  Kida et al., Jpn. J. Appl. Phys. 56, 07JG04 (2017) の計算再現
%  "Performance analysis of passive time reversal ... shallow sea"
%
%  論文の流れ:
%    ① Pekeris導波路のCIR生成（像法で自己完結。論文はnormal mode法）
%    ② 受信波動場行列 → ③ T-SVDで直達波 XD と 多重波 XM を分離（式6,7）
%    ④ X = XD + α·XM + N でSIR/SNRをパラメトリックに制御（式8）
%    ⑤ BPSK + RRC送信信号を生成、CIRで畳み込み
%    ⑥ PTR（全受信機で q関数集束, 式1）→ RLS-DFE(FF+FB+PLL) で残留ISI除去
%    ⑦ OSNR（式5）/ BER を評価
%
%  ★意図的な近似（論文との差, 開示事項）:
%    - チャネル: 論文=normal mode法, 本コード=像法Pekeris（物理的に同等な多重波
%      CIRだが固有音線の振幅位相は厳密には一致しない）。ARR差し替えも可。
%    - 信号: 論文=20kHz通過帯域+DDC, 本コード=等価ベースバンド（DDCで最終的に
%      ベースバンド化されるため処理結果は保存。搬送波依存効果は非対象）。
%    - PTRプローブ: 本コードは既知CIRを整合フィルタに使用（プローブ推定は省略）。
% ============================================================
clear; close all;

%% ===== 環境パラメータ（論文 Sect.3.2）=====
c1=1500; c2=2000; rho1=1000; rho2=1800;   % 音速[m/s]・密度[kg/m^3]（rho2は代表値）
D=20; zs=10; R=200;                        % 水深/音源深度/水平距離 [m]

%% ===== 信号パラメータ（論文 Sect.3.2）=====
fc=20e3; BW=8e3;              % 中心周波数20kHz, 帯域8kHz
beta=0.5;                    % RRCロールオフ
Rs=round(BW/(1+beta));       % シンボルレート（占有帯域≈BW）
sps=8; Fs=Rs*sps;            % 標本化（等価ベースバンド）
span=8;                      % RRCスパン[symbol]
isBPSK=true;                 % 論文：BPSK

% 論文の時間長：送信信号150ms（データ）, 全記録500ms
Tdata=0.150;                 % [s]
Ndata =round(Tdata*Rs);      % データシンボル数（≈150ms分）
Ntrain=round(0.4*Ndata);     % トレーニング前置（DFE収束用）
numSymbols=Ntrain+Ndata;

%% ===== 受信アレイ（論文 Fig.6：長さ一定で間引き）=====
z_top=0.1; z_bot=19.9;       % アレイ両端[m]（長さ19.8m固定）
NchList=[199 100 41 11 6];

%% ===== パラメトリック掃引（論文の3変数: No.Ch, SIR(α), SNR）=====
alphaList=[0 0.2 0.5 1 2 4];     % 多重波スケール（式8, 0→SIR=inf）
SNRList  =[20 6 0];              % 入力SNR [dB]（論文 Fig9/11のパネル）

%% ===== DFEパラメータ（論文：RLS, FF+FB+PLL）=====
Nf=20; Nb=10; lambda=0.999; delta=0.01; mu_pll=0.01;

%% ===== 送受信RRC =====
rrcF=rrc_filter(beta,span,sps);
gd_rrc=span*sps;             % 送信+受信整合の合計群遅延[サンプル]

%% ===== 結果格納 (Nch × α × SNR) =====
nN=numel(NchList); nA=numel(alphaList); nS=numel(SNRList);
OSNR=nan(nN,nA,nS); BER=nan(nN,nA,nS); SIRv=nan(nN,nA,nS); ISNR=nan(nN,nA,nS);

fprintf('Kida2017 PTR-DFE 再現 rev2 (BPSK, fc=%.0fkHz BW=%.0fkHz Rs=%d baud sps=%d)\n', fc/1e3,BW/1e3,Rs,sps);
fprintf('像法Pekeris: D=%.0fm zs=%.0fm R=%.0fm c1=%.0f c2=%.0f\n', D,zs,R,c1,c2);
fprintf('掃引: No.Ch=%s, α=%s, SNR=%s dB\n\n', mat2str(NchList), mat2str(alphaList), mat2str(SNRList));

for iN=1:nN
    Nch=NchList(iN);
    zr=linspace(z_top,z_bot,Nch).';

    %% ① 各受信機のCIR（像法Pekeris）
    Hraw=pekeris_cir_array(R,zs,zr,D,c1,c2,rho1,rho2,Fs,BW);

    %% ③ T-SVDで直達/多重分離（音源深度zs基準で整列）
    [HD,HM]=tsvd_separate_direct(Hraw,Fs,c1,R,zr,zs);

    for ia=1:nA
        a=alphaList(ia);
        Heff=HD+a*HM;                                  % ④ SIR制御

        % SIR（式4）: 直達RMS / 多重RMS
        AD=sqrt(mean(sum(abs(HD).^2,2)));
        AM=sqrt(mean(sum(abs(a*HM).^2,2)));
        sir = (AM>0)*20*log10(AD/max(AM,eps)) + (AM==0)*inf;

        for iS=1:nS
            SNR_dB=SNRList(iS);

            %% ⑤ 送信 BPSK + RRC
            bits=randi([0 1],numSymbols,1);
            symbols=2*bits-1;                          % BPSK ±1
            tx=upfirdn(symbols,rrcF,sps);

            %% ⑥-a 各受信機で畳込み＋雑音 → PTR相関蓄積
            rxset=cell(Nch,1); sigpow=0;
            for m=1:Nch
                r=conv(tx,Heff(m,:).'); rxset{m}=r; sigpow=sigpow+mean(abs(r).^2);
            end
            sigpow=sigpow/Nch; npow=sigpow/10^(SNR_dB/10);

            bufLen=length(tx)+2*size(Heff,2);
            sback=zeros(bufLen,1); qsum=zeros(2*size(Heff,2)-1,1);
            for m=1:Nch
                r=rxset{m}+sqrt(npow/2)*(randn(size(rxset{m}))+1i*randn(size(rxset{m})));
                g=conj(flipud(Heff(m,:).'));           % PTR整合フィルタ
                z=conv(r,g); L=min(length(z),bufLen); sback(1:L)=sback(1:L)+z(1:L);
                qac=conv(Heff(m,:).',g); Lq=min(length(qac),length(qsum)); qsum(1:Lq)=qsum(1:Lq)+qac(1:Lq);
            end

            % ⑥-b 受信整合RRC（RRC×RRC=レイズドコサインでNyquist化）
            sback=conv(sback,rrcF);

            % ISNR（式10）: (AD+AM)/An  ← 算術和
            ISNR(iN,ia,iS)=20*log10((AD+AM)/sqrt(npow));
            SIRv(iN,ia,iS)=sir;

            %% ⑥-c 焦点でシンボル抽出（q関数ピーク+群遅延、近傍探索で最良点）
            [~,qpk]=max(abs(qsum));
            off0=qpk+gd_rrc;
            ptrSym=sample_at_focus(sback,off0,sps,numSymbols,symbols,2*sps);

            %% ⑥-d RLS-DFE + PLL
            dfe=rls_dfe_pll(ptrSym,symbols,Nf,Nb,lambda,delta,Ntrain,mu_pll,isBPSK);

            %% ⑦ OSNR（式5）/ BER：判定指向区間
            idx=(Ntrain+1):numSymbols; s=symbols(idx); y=dfe(idx);
            g0=(s'*y)/(s'*s); yc=y/g0;
            OSNR(iN,ia,iS)=10*log10(sum(abs(s).^2)/sum(abs(s-yc).^2));
            BER(iN,ia,iS)=mean((real(yc)>0)~=(s>0));
        end
    end
    fprintf('Nch=%3d 完了\n', Nch);
end

%% ===== コンソール要約 =====
fprintf('\n===== OSNR[dB] 要約 =====\n');
for iS=1:nS
    fprintf('--- SNR=%d dB ---\n', SNRList(iS));
    fprintf('%6s', 'Nch\\α'); fprintf('%8.1f',alphaList); fprintf('\n');
    for iN=1:nN
        fprintf('%6d', NchList(iN)); fprintf('%8.2f',squeeze(OSNR(iN,:,iS))); fprintf('\n');
    end
end

%% ===== 図9相当：SIR–OSNR（SNRパネル別, No.Ch重ね）=====
figure('Name','SIR-OSNR (Fig.9相当)','Position',[60 60 380*nS 420]); mk='osd^v><';
for iS=1:nS
    subplot(1,nS,iS);
    for iN=1:nN
        plot(squeeze(SIRv(iN,:,iS)),squeeze(OSNR(iN,:,iS)),['-' mk(mod(iN-1,7)+1)],...
            'LineWidth',1.3,'DisplayName',sprintf('%dCh',NchList(iN))); hold on;
    end
    grid on; xlabel('SIR [dB]'); ylabel('OSNR [dB]');
    title(sprintf('SNR=%d dB',SNRList(iS))); if iS==1, legend('Location','best'); end
end
sgtitle('SIR–OSNR 関係（Fig.9相当）','FontWeight','bold');

%% ===== 図12相当：ISNR–OSNR と 各No.Ch の理論限界 [ISNR+Ga] =====
iS0=1;  % SNR=20dB代表
figure('Name','ISNR-OSNR (Fig.12相当)','Position',[80 80 640 480]);
for iN=1:nN
    plot(squeeze(ISNR(iN,:,iS0)),squeeze(OSNR(iN,:,iS0)),['-' mk(mod(iN-1,7)+1)],...
        'LineWidth',1.3,'DisplayName',sprintf('%dCh',NchList(iN))); hold on;
end
xg=linspace(min(ISNR(:),[],'omitnan')-2,max(ISNR(:),[],'omitnan')+2,50);
for iN=1:nN
    plot(xg,xg+10*log10(NchList(iN)),'--','HandleVisibility','off');
end
grid on; xlabel('ISNR [dB]'); ylabel('OSNR [dB]'); legend('Location','best');
title(sprintf('ISNR–OSNR と理論限界 OSNR=ISNR+Ga（SNR=%d dB, 破線=各Ch）',SNRList(iS0)));

%% ================= ローカル関数 =================
function h = rrc_filter(beta, span, sps)
    N=span*sps; t=(-N/2:N/2)/sps; h=zeros(size(t));
    for i=1:numel(t)
        tt=t(i);
        if tt==0
            h(i)=1-beta+4*beta/pi;
        elseif abs(abs(tt)-1/(4*beta))<1e-8
            h(i)=(beta/sqrt(2))*((1+2/pi)*sin(pi/(4*beta))+(1-2/pi)*cos(pi/(4*beta)));
        else
            num=sin(pi*tt*(1-beta))+4*beta*tt.*cos(pi*tt*(1+beta));
            den=pi*tt.*(1-(4*beta*tt).^2); h(i)=num/den;
        end
    end
    h=h/sqrt(sum(h.^2));
end

function H = pekeris_cir_array(R,zs,zr,D,c1,c2,rho1,rho2,Fs,BW)
% 像法によるPekeris導波路CIR（各受信機）
    Nch=numel(zr); Nimg=60;
    arr=cell(Nch,1); tmin=inf; tmax=0;
    for m=1:Nch
        [amp,tau]=pekeris_eigenrays(R,zs,zr(m),D,c1,c2,rho1,rho2,Nimg);
        arr{m}=struct('amp',amp,'tau',tau);
        tmin=min(tmin,min(tau)); tmax=max(tmax,max(tau));
    end
    bl=blimpulse(Fs,BW); klen=(numel(bl)-1)/2;
    Lh=round((tmax-tmin)*Fs)+2*klen+8; H=zeros(Nch,Lh);
    for m=1:Nch
        a=arr{m}.amp; tau=arr{m}.tau-tmin; row=zeros(1,Lh+2*klen);
        for p=1:numel(tau)
            n0=round(tau(p)*Fs)+klen+1; idx=n0+(-klen:klen);
            ok=idx>=1&idx<=numel(row); row(idx(ok))=row(idx(ok))+a(p)*bl(ok);
        end
        H(m,:)=row(klen+1:klen+Lh);
    end
end

function [amp,tau]=pekeris_eigenrays(R,zs,zr,D,c1,c2,rho1,rho2,Nimg)
% Pekeris導波路の固有音線（像法, 2系統）
%  表面反射=-1（圧力解放）, 海底反射=Rb（レイリー, grazing角依存）
    amp=[]; tau=[];
    for n=-Nimg:Nimg
        % 系統A: 像深度 2nD+zs, 反射回数 s=|n|, b=|n|
        addray(2*n*D+zs, abs(n), abs(n));
        % 系統B: 像深度 2nD-zs
        if n>=1, s=n-1; b=n; else, s=abs(n)+1; b=abs(n); end
        addray(2*n*D-zs, s, b);
    end
    [tau,si]=sort(tau); amp=amp(si);
    function addray(imgz, ns, nb)
        dz=imgz-zr; Ls=sqrt(R^2+dz^2); t=Ls/c1;
        grazing=atan2(abs(dz),R);                     % 水平からの角（=grazing）
        Rb=rayleigh_refl(grazing,c1,c2,rho1,rho2);
        A=((-1)^ns)*(Rb^nb)/(4*pi*Ls);
        if abs(A)>1e-6, amp(end+1,1)=A; tau(end+1,1)=t; end
    end
end

function Rb=rayleigh_refl(grazing,c1,c2,rho1,rho2)
% 平面波レイリー反射係数（水→堆積層, 入力=grazing角）
    st1=sin(grazing);
    s2=(c2/c1)*cos(grazing);           % スネル則（臨界: grazing<acos(c1/c2)で全反射）
    if s2>=1
        Rb=1;
    else
        st2=sqrt(1-s2^2);
        Z1=rho1*c1/max(st1,eps); Z2=rho2*c2/st2;
        Rb=(Z2-Z1)/(Z2+Z1);
    end
end

function bl=blimpulse(Fs,BW)
    K=round(4*Fs/BW); n=-K:K;
    bl=(BW/Fs)*sinc((BW/Fs)*n).*hamming(2*K+1)'; bl=bl/max(abs(bl));
end

function [HD,HM]=tsvd_separate_direct(H,Fs,c1,R,zr,zs)
% T-SVDで直達(低ランク)/多重を分離（式6,7）。音源深度zs基準で直達フロント整列
    [Nch,L]=size(H);
    tau_d=sqrt(R^2+(zr-zs).^2)/c1; shift=round((tau_d-min(tau_d))*Fs);
    Hal=zeros(Nch,L); for m=1:Nch, Hal(m,:)=circshift(H(m,:),-shift(m)); end
    [U,S,V]=svd(Hal,'econ'); sv=diag(S);
    % 直達ランク：特異値の明確な落ち込み（elbow）で選択。既定は最大1〜2
    if numel(sv)>=2 && sv(2)<0.3*sv(1), m_rank=1; else, m_rank=min(2,numel(sv)); end
    HDal=U(:,1:m_rank)*S(1:m_rank,1:m_rank)*V(:,1:m_rank)'; HMal=Hal-HDal;
    HD=zeros(Nch,L); HM=zeros(Nch,L);
    for m=1:Nch, HD(m,:)=circshift(HDal(m,:),shift(m)); HM(m,:)=circshift(HMal(m,:),shift(m)); end
end

function sym=sample_at_focus(sback,off0,sps,numSymbols,symbols,srch)
    best=-inf; sym=zeros(numSymbols,1);
    for d=-srch:srch
        off=off0+d; idx=off+(0:numSymbols-1)*sps;
        if idx(1)<1||idx(end)>length(sback), continue; end
        cand=sback(idx); g0=(symbols'*cand)/(symbols'*symbols);
        if abs(g0)<eps, continue; end
        cc=cand/g0; osnr=sum(abs(symbols).^2)/sum(abs(symbols-cc).^2);
        if osnr>best, best=osnr; sym=cand; end
    end
end

function out=rls_dfe_pll(x,symbols,Nf,Nb,lambda,delta,Ntrain,mu_pll,isBPSK)
% RLS判定帰還等化器 + PLL（論文 Fig.2）。PLL回転はFFブランチのみに作用
    N=length(x); Ntap=Nf+Nb;
    w=zeros(Ntap,1); w(1)=1; P=(1/delta)*eye(Ntap);
    ff=zeros(Nf,1); fb=zeros(Nb,1); theta=0; out=zeros(N,1);
    for n=1:N
        ff=[x(n); ff(1:end-1)];
        u=[ff*exp(-1i*theta); -fb];      % FFのみ位相回転（FBは復調済みシンボル領域）
        y=w'*u; out(n)=y;
        if isBPSK, d=sign(real(y)); else, d=sign(real(y))+1i*sign(imag(y)); end
        if n<=Ntrain, s_ref=symbols(n); else, s_ref=d; end
        e=s_ref-y;
        Pu=P*u; kap=Pu/(lambda+u'*Pu);
        w=w+kap*conj(e); P=(P-kap*(Pu'))/lambda;   % RLS（regressorが回転済なので追加回転不要）
        phi_err=imag(y*conj(s_ref)); theta=theta+mu_pll*phi_err;  % PLL更新（w更新の後）
        fb=[s_ref; fb(1:end-1)];
    end
end
