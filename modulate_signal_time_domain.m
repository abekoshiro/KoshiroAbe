function [output,papr,rrc]=modulate_signal_time_domain(symbol,Fs,TLength,BW,Fc,shift,Alpha)
%- Input
% BW: Symbol Rate (Hz)
% Ncycle : RRC cycle 20< recommended
% Alpha : Roll-off rate 0 < Alpha < 1
% Fs : Sampling Rate (Hz)
% symbol : Tx. communiation symbols
% shift : pre zero-padding length (samples)
%- Output
% output: output passband signal (for baseband, Fc=0 )
% papr : peak to average power ratio of modulated signal
% rrc : time domain raised root cosine filter
%%
SPS=Fs/BW;
ns=length(symbol);
if TLength ==0
    % insert zeros for 10 symbols
    N=ns*SPS+10*SPS;
    temp=zeros(N,1);
    TLength=N/Fs;
else
    N=round(TLength*Fs);
    temp=zeros(N,1);
end
dt=1/Fs;
% df=1/TLength;
TAxis=0:dt:TLength-dt;
% FAxis=0:df:Fs-df;
Ncycle=30; %for this value, practically with Alpha=0.15, ISI ~ 54 dB
% Ncycle=20; %for this value, practically with Alpha=0.15, ISI ~ 36 dB
%% gen. root raised cosine filter
rrc = rrc_time_domain(BW, Ncycle, Alpha, Fs);
%% gen base-band signal
temp(shift + (1 : SPS : ns * SPS))=symbol; % insert symbol in time-series
tempr = conv(temp, rrc); % apply RRC
tempr=tempr(Ncycle * SPS + (1 : N)); % remove offset
%% gen passband signal
if Fc~=0
    tempe=tempr .* exp(2i * pi * Fc * TAxis.'); % frequency modulation
    output=real(tempe); % passband signal is real
end
output=output/rms(output(shift+(1:ns*SPS))); % amplitude normalization
%
papr=20*log10(max(abs(output))/rms(output(shift+(1:ns*SPS))));