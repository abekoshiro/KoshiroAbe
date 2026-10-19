function [rrc,tTs] = rrc_time_domain(BW, Ncycle, alpha, Fs)
Ts = 1/BW;
dt=1/Fs;
t = -Ncycle*Ts:dt:Ncycle*Ts;
tTs = t/Ts;
rrc = (sin(pi*tTs*(1-alpha)) + 4*alpha*tTs.*cos(pi*tTs*(1+alpha)))./(pi*tTs.*(1-(4*alpha*tTs).^2));
rrc(t==0)=(1+alpha*(4/pi-1));
rrc(abs(t)==Ts/(4*alpha))=alpha/sqrt(2)*((1+2/pi)*sin(pi/(4*alpha))+(1-2/pi)*cos(pi/(4*alpha)));
end