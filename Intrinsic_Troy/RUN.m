Tau= 0 ; %1/(6*D);
p=0;
% Set magnetic parameters in frequency units
g0=2.0067;
g2=-.00658;
gpar=g0 - 2*g2/3;
gprp=g0 + g2/3;


defaults_troy;
stv_troy;
mat_troy_diff;
spec_troy;
dderiv_troy;
unique_deriv=pD;

unique=absorption;

c=frequency;
figure
plot(c,pD)
title('First derivative with respect to diffusion spectra')
xlabel('Frequency rad/sec')
ylabel('Intensity')
% figure
% plot(c,unique)
% %flipped=absorption(end:-1:1);
% %BATMAN=flipped+absorption;
% %plot(frequency,BATMAN)
% title('Absorption with rotational diffusion parameter: 1e8')
% xlabel('Frequency rad/sec')
% ylabel('Intensity')

% Set magnetic parameters in frequency units
g0=2.0067;
g2=-.00658;
gpar=g0 + 2*g2/3;
gprp=g0 - g2/3;


defaults_troy;
stv_troy;
mat_troy_diff;
spec_troy;
absorption2=absorption;
dderiv_troy;

d=frequency;
% hold on
% plot(d, absorption2)
% %flipped=absorption(end:-1:1);
% %BATMAN=flipped+absorption;
% %plot(frequency,BATMAN)
% title('Absorption with rotational diffusion parameter: 1e7')
% xlabel('Frequency rad/sec')
% ylabel('Intensity')
hold on
plot(d,pD)

% figure
% plot(frequency, pD)
% title('pD')
% figure
% plot(frequency, pDel)
% title('pDel')
% figure
% plot(frequency, pLam)
% title('pLam')
% figure
% plot(frequency, pDD)
% title('pDD')


