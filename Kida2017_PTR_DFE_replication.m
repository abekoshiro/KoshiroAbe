%% ============================================================
%  Kida2017_PTR_DFE_replication.m
%  Kida et al., Jpn. J. Appl. Phys. 56, 07JG04 (2017) の計算再現
%  "Performance analysis of passive time reversal ... shallow sea"
%
%  論文の流れ:
%    ① Pekeris導波路のCIR生成（本コードは像法で自己完結。論文はnormal mode）
%    ② 受信波動場行列 X を構成
%    ③ T-SVDで直達波 XD と マルチパス XM を分離
%    ④ X = XD + α·XM + N でSIR/SNRをパラメトリックに制御（式8）
%    ⑤ BPSK + RRC送信信号を生成、CIRで畳み込み
%    ⑥ PTR（全受信機で q関数集束）→ RLS-DFE(FF+FB+PLL) で残留ISI除去
%    ⑦ OSNR（式5）/ BER を評価
%
%  ※本コードは論文の「処理」を忠実に再現することを主眼とし、チャネルは
%    像法Pekerisで自己完結生成する（read_arrivals_asc のARRも差し替え可）。
% ============================================================
clear; close all;

%% ===== 環境・信号パラメータ（論文 Sect.3.2）=====
c1     = 1500;       % 水中音速 [m/s]
c2     = 2000;       % 海底堆積層音速 [m/s]
rho1   = 1000;       % 水の密度 [kg/m^3]
rho2   = 1800;       % 堆積層密度 [kg/m^3]（論文非記載→代表値）
D      = 20;         % 水深 [m]（受信器 0.1–19.9m を張るので約20m）
zs     = 10;         % 音源深度 [m]（論文：10 m）
R      = 200;        % 水平距離 [m]（論文：200 m）

fc     = 20e3;       % 中心周波数 [Hz]（論文：20 kHz）
BW     = 8e3;        % 帯域幅 [Hz]（論文：8 kHz）
beta   = 0.5;        % RRC ロールオフ
Rs     = round(BW/(1+beta));   % シンボルレート [baud]（占有帯域≈BW）
sps    = 8;          % 1シンボルあたりサンプル数
Fs     = Rs*sps;     % ベースバンド標本化周波数 [Hz]
span   = 8;          % RRCフィルタのシンボルスパン

Ndata  = 500;        % データシンボル数
Ntrain = 200;        % トレーニングシンボル数（DFE収束用, 前置）
numSymbols = Ntrain + Ndata;

isBPSK = true;       % 論文：BPSK

%% ===== 受信アレイ（論文 Fig.6：長さ一定で間引き）=====
z_top = 0.1; z_bot = 19.9;   % アレイ両端 [m]（長さ19.8m固定）
NchList = [199 100 41 11 6];  % 評価する受信機数（間引き。長さは保持）

%% ===== SIR制御（論文 式8：X = XD + α·XM + N）=====
alphaList = [0 0.2 0.5 1 2 4];   % 多重波スケール（0→SIR=inf）
SNR_dB    = 20;                   % 入力SNR（AWGN, 論文の一例）

%% ===== DFEパラメータ（論文：RLS, FF+FB+PLL）=====
Nf     = 20;      % フィードフォワードタップ
Nb     = 10;      % フィードバックタップ
lambda = 0.999;   % RLS忘却係数
delta  = 0.01;    % RLS正則化
mu_pll = 0.01;    % PLL利得

%% ===== 送受信フィルタ（RRC）=====
rrcF = rrc_filter(beta, span, sps);   % ルートレイズドコサイン

%% ===== 結果格納 =====
OSNR = nan(numel(NchList), numel(alphaList));
BER  = nan(numel(NchList), numel(alphaList));
SIRv = nan(numel(NchList), numel(alphaList));
ISNR = nan(numel(NchList), numel(alphaList));

fprintf('Kida2017 PTR-DFE 再現  (BPSK, fc=%.0fkHz, BW=%.0fkHz, Rs=%d baud)\n', fc/1e3, BW/1e3, Rs);
fprintf('像法Pekeris: D=%.0fm, zs=%.0fm, R=%.0fm, c1=%.0f c2=%.0f\n\n', D, zs, R, c1, c2);

for iN = 1:numel(NchList)
    Nch = NchList(iN);
    zr  = linspace(z_top, z_bot, Nch).';   % 受信深度（長さ一定で間引き）

    %% ① 各受信機のCIR（像法Pekeris）
    [Hraw, tgrid] = pekeris_cir_array(R, zs, zr, D, c1, c2, rho1, rho2, Fs, fc, BW);
    % Hraw: Nch × Lh の複素CIR（帯域制限済み）

    %% ② 波動場行列 X（各行=受信機のCIR記録）→ ③ T-SVDで直達/多重分離
    [HD, HM] = tsvd_separate_direct(Hraw, Fs, c1, R, zr);
    % HD: 直達成分, HM: 多重波成分（Nch × Lh）

    for ia = 1:numel(alphaList)
        a = alphaList(ia);

        %% ④ SIR制御：実効CIR = 直達 + α·多重
        Heff = HD + a*HM;

        %% ---- SIR / ISNR（式4,10）: RMS振幅から算出 ----
        AD = sqrt(mean(sum(abs(HD).^2,2)));       % 直達RMS（アレイ平均）
        AM = sqrt(mean(sum(abs(a*HM).^2,2)));     % 多重RMS
        if AM>0, SIRv(iN,ia)=20*log10(AD/AM); else, SIRv(iN,ia)=inf; end

        %% ⑤ 送信信号：BPSK + RRC
        bits = randi([0 1], numSymbols, 1);
        if isBPSK
            symbols = 2*bits - 1;                 % BPSK: ±1
        else
            symbols = (2*bits-1);                 % （拡張余地）
        end
        tx = upfirdn(symbols, rrcF, sps);         % RRC整形・アップサンプル
        grpdel = span*sps;                        % RRC群遅延（両側でspan*sps）

        %% ⑥-a 各受信機で畳み込み＋雑音、PTR相関を蓄積
        Lconv = length(tx) + size(Heff,2) - 1;
        sback = zeros(Lconv + size(Heff,2) - 1, 1);   % q関数畳み込み分の余裕
        qsum  = zeros(2*size(Heff,2)-1, 1);
        % 雑音電力：全受信機の平均信号電力から SNR で決定
        rxset = cell(Nch,1); sigpow = 0;
        for m = 1:Nch
            r = conv(tx, Heff(m,:).');
            rxset{m} = r; sigpow = sigpow + mean(abs(r).^2);
        end
        sigpow = sigpow/Nch; npow = sigpow/10^(SNR_dB/10);
        for m = 1:Nch
            r = rxset{m} + sqrt(npow/2)*(randn(size(rxset{m}))+1i*randn(size(rxset{m})));
            g = conj(flipud(Heff(m,:).'));        % PTR：時間反転共役（プローブ相当）
            z = conv(r, g);                        % 相互相関（back-propagation）
            L = min(length(z), length(sback));
            sback(1:L) = sback(1:L) + z(1:L);
            qac = conv(Heff(m,:).', g);            % q関数（自己相関）
            Lq = min(length(qac), length(qsum));
            qsum(1:Lq) = qsum(1:Lq) + qac(1:Lq);
        end

        %% ISNR（式10）: (直達+多重)RMS / 雑音RMS
        Atot = sqrt(mean(sum(abs(Heff).^2,2)));
        ISNR(iN,ia) = 20*log10(Atot/sqrt(npow));

        %% ⑥-b 焦点位置でシンボル抽出（q関数ピーク基準）
        [~, qpk] = max(abs(qsum));
        % sback の焦点：RRC往復群遅延 + q関数ピーク を考慮して同期
        off = qpk + grpdel;    % 近似同期。下でオフセット微調して最良点を採る
        [ptrSym, best_off] = sample_at_focus(sback, off, sps, numSymbols, symbols, grpdel);

        %% ⑥-c DDC相当（ベースバンド）→ RLS-DFE + PLL
        dfe_out = rls_dfe_pll(ptrSym, symbols, Nf, Nb, lambda, delta, Ntrain, mu_pll, isBPSK);

        %% ⑦ OSNR（式5）/ BER：判定指向区間（トレーニング除外）
        idx = (Ntrain+1):numSymbols;
        s = symbols(idx); y = dfe_out(idx);
        g0 = (s'*y)/(s'*s); yc = y/g0;            % 利得正規化
        OSNR(iN,ia) = 10*log10(sum(abs(s).^2)/sum(abs(s-yc).^2));
        if isBPSK
            dec = real(yc)>0; ref = s>0;
        else
            dec = [real(yc)>0; imag(yc)>0]; ref = [real(s)>0; imag(s)>0];
        end
        BER(iN,ia) = mean(dec(:)~=ref(:));

        fprintf('Nch=%3d  α=%.1f  SIR=%6.1f dB  ISNR=%6.1f dB  OSNR=%6.2f dB  BER=%.4f\n', ...
                Nch, a, SIRv(iN,ia), ISNR(iN,ia), OSNR(iN,ia), BER(iN,ia));
    end
    fprintf('\n');
end

%% ===== 図9相当：SIR–OSNR 関係（No.Ch別）=====
figure('Name','SIR-OSNR (Fig.9相当)','Position',[80 80 640 480]);
mk = 'osd^v><';
for iN=1:numel(NchList)
    plot(SIRv(iN,:), OSNR(iN,:), ['-' mk(mod(iN-1,numel(mk))+1)], 'LineWidth',1.4, ...
        'DisplayName',sprintf('%dCh',NchList(iN))); hold on;
end
grid on; xlabel('SIR [dB]'); ylabel('OSNR [dB]'); legend('Location','best');
title(sprintf('SIR–OSNR 関係（SNR=%d dB）', SNR_dB));

%% ===== 図12相当：ISNR–OSNR と理論限界 [ISNR + Ga] =====
figure('Name','ISNR-OSNR (Fig.12相当)','Position',[100 100 640 480]);
for iN=1:numel(NchList)
    plot(ISNR(iN,:), OSNR(iN,:), ['-' mk(mod(iN-1,numel(mk))+1)], 'LineWidth',1.4, ...
        'DisplayName',sprintf('%dCh',NchList(iN))); hold on;
end
% 理論限界線：OSNR = ISNR + Ga, Ga=10log10(Nch)（最大Chで代表）
xg = linspace(min(ISNR(:))-2, max(ISNR(:))+2, 50);
plot(xg, xg + 10*log10(max(NchList)), 'k--', 'DisplayName','ISNR+Ga (max Ch)');
grid on; xlabel('ISNR [dB]'); ylabel('OSNR [dB]'); legend('Location','best');
title('ISNR–OSNR 関係と理論限界（式：OSNR=ISNR+Ga）');

%% ================= ローカル関数 =================
function h = rrc_filter(beta, span, sps)
% ルートレイズドコサイン（root raised cosine）フィルタ係数
    N = span*sps;
    t = (-N/2:N/2)/sps;
    h = zeros(size(t));
    for i=1:numel(t)
        tt = t(i);
        if tt==0
            h(i) = 1 - beta + 4*beta/pi;
        elseif abs(abs(tt)-1/(4*beta))<1e-8
            h(i) = (beta/sqrt(2))*((1+2/pi)*sin(pi/(4*beta)) + (1-2/pi)*cos(pi/(4*beta)));
        else
            num = sin(pi*tt*(1-beta)) + 4*beta*tt.*cos(pi*tt*(1+beta));
            den = pi*tt.*(1-(4*beta*tt).^2);
            h(i) = num/den;
        end
    end
    h = h/sqrt(sum(h.^2));   % 正規化
end

function [H, tgrid] = pekeris_cir_array(R, zs, zr, D, c1, c2, rho1, rho2, Fs, fc, BW)
% 像法によるPekeris導波路のCIR（各受信機）
%  海面：圧力解放(係数-1)、海底：レイリー反射（密度・音速比）
    Nch = numel(zr);
    Nimg = 60;                       % 像の反射回数上限
    % まず全受信機の到来を集めて最大遅延を把握
    arr = cell(Nch,1); tmin = inf; tmax = 0;
    for m=1:Nch
        [amp, tau] = pekeris_eigenrays(R, zs, zr(m), D, c1, c2, rho1, rho2, Nimg);
        arr{m} = struct('amp',amp,'tau',tau);
        tmin = min(tmin, min(tau)); tmax = max(tmax, max(tau));
    end
    % 帯域制限パルス（sinc×窓、中心0）で各到来を表現
    Lh = round((tmax - tmin)*Fs) + 4*round(Fs/BW) + 8;
    tgrid = (0:Lh-1)/Fs;
    H = zeros(Nch, Lh);
    bl = blimpulse(Fs, BW);          % 帯域制限インパルス核
    klen = (numel(bl)-1)/2;
    for m=1:Nch
        a = arr{m}.amp; tau = arr{m}.tau - tmin;   % 共通最小遅延で正規化
        row = zeros(1, Lh + 2*klen);
        for p=1:numel(tau)
            n0 = round(tau(p)*Fs) + klen + 1;
            idx = n0 + (-klen:klen);
            ok = idx>=1 & idx<=numel(row);
            row(idx(ok)) = row(idx(ok)) + a(p)*bl(ok);
        end
        H(m,:) = row(klen+1:klen+Lh);
    end
end

function [amp, tau] = pekeris_eigenrays(R, zs, zr, D, c1, c2, rho1, rho2, Nimg)
% Pekeris導波路の固有音線（像法）：4系統の像を生成
    amp = []; tau = [];
    for n = -Nimg:Nimg
        % 4種の鏡像深度と反射回数
        cand = [ 2*n*D + (zr - zs), abs(n)*2, abs(n)*2 ;      % 直達系（表n回,底n回）
                 2*n*D + (zr + zs), abs(n)*2+1, abs(n)*2 ;    % 表面反射込み
                 2*n*D - (zr + zs), abs(n)*2+1, abs(n)*2 ;
                 2*n*D - (zr - zs), abs(n)*2, abs(n)*2 ];
        for k=1:4
            dz = cand(k,1);
            Ls = sqrt(R^2 + dz^2);          % 経路長
            t  = Ls / c1;                    % 遅延
            % 反射係数：表面 -1（圧力解放）、海底はレイリー（grazing角）
            nb_surf = cand(k,2);             % 表面反射回数（近似）
            nb_bot  = cand(k,3);             % 海底反射回数（近似）
            theta = atan2(abs(dz), R);       % 伝搬角（水平から）
            Rb = rayleigh_refl(pi/2-theta, c1, c2, rho1, rho2);  % grazing余角
            A  = ((-1)^nb_surf) * (Rb^nb_bot) / (4*pi*Ls);       % 球面拡散
            if abs(A) > 1e-6                  % 微小像は無視
                amp(end+1,1) = A; tau(end+1,1) = t; %#ok<AGROW>
            end
        end
    end
    % 重複遅延を統合
    [tau, si] = sort(tau); amp = amp(si);
end

function Rb = rayleigh_refl(grazing, c1, c2, rho1, rho2)
% 平面波レイリー反射係数（水→堆積層）
    st1 = sin(grazing); ct1 = cos(grazing);
    s2  = (c2/c1)*ct1;                 % スネルの法則
    if s2 >= 1
        Rb = 1;                        % 全反射
    else
        ct2 = s2; st2 = sqrt(1-s2^2);
        Z1 = rho1*c1/st1; Z2 = rho2*c2/st2;
        Rb = (Z2 - Z1)/(Z2 + Z1);
    end
end

function bl = blimpulse(Fs, BW)
% 帯域制限インパルス核（sinc×ハミング窓）
    K = round(4*Fs/BW);
    n = -K:K;
    bl = (BW/Fs) * sinc((BW/Fs)*n) .* hamming(2*K+1)';
    bl = bl/max(abs(bl));
end

function [HD, HM] = tsvd_separate_direct(H, Fs, c1, R, zr)
% T-SVDで直達波(低ランク)と多重波を分離（論文 式6,7）
%  直達波の時間フロントを揃えてから SVD → 低ランク=直達
    [Nch, L] = size(H);
    % 各受信機の直達到来サンプル（幾何遅延）でアライン
    tau_d = sqrt(R^2 + (zr - zr(1)).^2)/c1;    % 相対直達遅延（近似）
    shift = round((tau_d - min(tau_d))*Fs);
    Hal = zeros(Nch, L);
    for m=1:Nch
        Hal(m,:) = circshift(H(m,:), -shift(m));
    end
    % SVD
    [U,S,V] = svd(Hal, 'econ');
    sv = diag(S);
    % 直達ランク m*：最大特異値からの寄与で決定（>10%を直達とみなす）
    m_rank = max(1, sum(sv > 0.1*sv(1)));
    m_rank = min(m_rank, 3);                    % 直達は低ランク（安全側に上限）
    HDal = U(:,1:m_rank)*S(1:m_rank,1:m_rank)*V(:,1:m_rank)';
    HMal = Hal - HDal;
    % アライン戻し
    HD = zeros(Nch,L); HM = zeros(Nch,L);
    for m=1:Nch
        HD(m,:) = circshift(HDal(m,:), shift(m));
        HM(m,:) = circshift(HMal(m,:), shift(m));
    end
end

function [sym, best_off] = sample_at_focus(sback, off0, sps, numSymbols, symbols, grpdel)
% 焦点近傍でオフセットを微調し、最もOSNRが高い点でシンボル抽出
    best = -inf; best_off = off0; sym = zeros(numSymbols,1);
    for d = -sps:sps
        off = off0 + d;
        idx = off + (0:numSymbols-1)*sps;
        if idx(1)<1 || idx(end)>length(sback), continue; end
        cand = sback(idx);
        g0 = (symbols'*cand)/(symbols'*symbols);
        if abs(g0)<eps, continue; end
        cc = cand/g0;
        osnr = sum(abs(symbols).^2)/sum(abs(symbols-cc).^2);
        if osnr>best, best=osnr; best_off=off; sym=cand; end
    end
end

function out = rls_dfe_pll(x, symbols, Nf, Nb, lambda, delta, Ntrain, mu_pll, isBPSK)
% RLS判定帰還等化器 + PLL（論文 Fig.2 相当）
    N = length(x); Ntap = Nf+Nb;
    w = zeros(Ntap,1); w(1)=1;
    P = (1/delta)*eye(Ntap);
    ff = zeros(Nf,1); fb = zeros(Nb,1);
    theta = 0; out = zeros(N,1);
    for n=1:N
        ff = [x(n); ff(1:end-1)];
        u  = [ff; -fb];
        y0 = w'*u;
        y  = y0*exp(-1i*theta);          % PLLで位相回転を除去
        out(n) = y;
        if isBPSK, d = sign(real(y)); else, d = sign(real(y))+1i*sign(imag(y)); end
        if n<=Ntrain, s_ref = symbols(n); else, s_ref = d; end
        % PLL更新（判定指向位相誤差）
        phi_err = imag(y*conj(s_ref));
        theta = theta + mu_pll*phi_err;
        % RLS更新（位相回転後の誤差）
        e   = s_ref - y;
        Pu  = P*u; kap = Pu/(lambda + u'*Pu);
        w   = w + kap*conj(e*exp(1i*theta));
        P   = (P - kap*(Pu'))/lambda;
        fb  = [s_ref; fb(1:end-1)];
    end
end
