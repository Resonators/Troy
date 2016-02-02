Tau= 0 ; %1/(6*D);
p=0;
defaults_troy;
stv_troy;
mat_troy_diff;
spec_troy;
% dderiv_troy;
figure
plot(frequency, absorption)
title('Absorption')
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

