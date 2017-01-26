%perform check on derivatives with the foundamental theorem of calculus.
% dSpectrum = [spectrum(x+e) - spectrum(x)]/e as e tends to zero.
% Code should calculate a spectrum with the variable x and the spectrum
% with variable x+e which we call deltaSpectrum
Tau=0;
p=0;
RUN;
absorption1=absorption;
figure;
plot(frequency, absorption1)
k=.00000001;
stv_troy;
mat_troy_diff_SmallChange;
spec_troy;
absorption2=absorption;
hold on
plot(frequency, absorption2,'red')
dabsorption= (absorption2 - absorption1)./k;
figure;
plot(frequency, dabsorption)
title('pD Fun-da-MENTAL Theorem of Calculus')

