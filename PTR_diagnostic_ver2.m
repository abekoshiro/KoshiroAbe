%% ============================================================
%  PTR 診断スクリプト ver2（新ARR：複数距離×複数深度 対応）
%  目的：全音源について q関数の鋭さ(PSR)とPTR単体OSNRを一覧比較し、
%        どの音源が物理的に厳しいかを切り分ける
% ============================================================
[Arr, Pos] = read_arrivals_asc('JamTank_Arr.arr');

Nrr = size(Arr,1);   % 受信距離数
Nrz = size(Arr,2);   % 受信深度数
Nsd = size(Arr,3);   % 音源数
NUM_RX = Nrr * Nrz;

fprintf('ARR構成: 距離%d × 深度%d = %dch, 音源%d\n', Nrr, Nrz, NUM_RX, Nsd);
if isfield(Pos,'r')&&isfield(Pos.r,'r'), fprintf('  距離: '); fprintf('%.1f ',Pos.r.r); fprintf('m\n'); end
if isfield(Pos,'r')&&isfield(Pos.r,'z'), fprintf('  深度: '); fprintf('%.2f ',Pos.r.z); fprintf('m\n\n'); end

% チャネル ⇔ (距離,深度) 対応
ch_rr=zeros(NUM_RX,1); ch_rz=zeros(NUM_RX,1); cc=0;
for irr=1:Nrr, for irz=1:Nrz, cc=cc+1; ch_rr(cc)=irr; ch_rz(cc)=irz; end, end

Rs = 500; Sps = 16; Fs = Rs*Sps; numSymbols = 2000;

PSR_vec  = zeros(Nsd,1);
OSNR_vec = zeros(Nsd,1);

for idx_sd = 1:Nsd

    %% 送信信号
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
    total_ptr=zeros(ptr_buf_len,1);
    total_ac =zeros(ptr_buf_len,1);

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
        g=conj(flipud(h));
        z=conv(rx,g); hac=conv(h,g);
        Lz=min(length(z),ptr_buf_len); La=min(length(hac),ptr_buf_len);
        total_ptr(1:Lz)=total_ptr(1:Lz)+z(1:Lz);
        total_ac(1:La) =total_ac(1:La)+hac(1:La);
    end

    %% q関数 PSR（シンボル間隔サンプリング）
    [qpk,pk]=max(abs(total_ac));
    lags=(-numSymbols:numSymbols); qsym=zeros(length(lags),1);
    for ii=1:length(lags)
        idx=pk+lags(ii)*Sps;
        if idx>=1&&idx<=length(total_ac), qsym(ii)=total_ac(idx); end
    end
    main=abs(qsym(lags==0))^2; side=sum(abs(qsym(lags~=0)).^2);
    PSR_vec(idx_sd)=10*log10(main/side);

    %% PTR単体OSNR
    off=pk+round(Sps/2)-1;
    ps=total_ptr(off:Sps:off+Sps*(numSymbols-1));
    a=(symbols'*ps)/(symbols'*symbols); pc=ps/a;
    OSNR_vec(idx_sd)=10*log10(sigPower/mean(abs(symbols-pc).^2));

    fprintf('音源%d:  q関数PSR = %6.2f dB,   PTR単体OSNR = %6.2f dB\n', ...
            idx_sd, PSR_vec(idx_sd), OSNR_vec(idx_sd));
end

%% PSR vs OSNR の関係を可視化
figure('Name','PSR vs OSNR per source');
bar([PSR_vec, OSNR_vec]);
set(gca,'XTickLabel',arrayfun(@(k)sprintf('音源%d',k),1:Nsd,'UniformOutput',false));
legend({'q関数 PSR (上限の目安)','PTR単体 OSNR (実測)'},'Location','best');
ylabel('dB'); grid on;
title('音源ごとの PSR(理論上限) vs 実測OSNR');
