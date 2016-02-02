% needs defaults.m stv.m mat.m

% Initialize data arrays

absorption = [];
frequency  = [];

tempom = om;

for n = 1:pts

  shift(1:length(alpha)) = 1i*(tempom-om0);

  C = diag(shift) + diag(alpha) + diag(beta,1) + diag(beta,-1) ...
    + diag(gamma,2) + diag(gamma,-2);

  temp = svec*(C\svec')/pi;
  absorption(n) = real(temp);

  frequency(n)  = tempom;
  tempom   = tempom + 2.5*F/pts;
end

if (tempom < om)
  frequency = frequency(end:-1:1,:);
  absorption = absorption(end:-1:1,:);
end

%Derivative spectra:

% dDeltempom = om;
% for n = 1:pts
% 
%   shift(1:length(alphaDel)) = 1i*(dDeltempom-om0);
% 
%   dDelC = diag(shift) + diag(alphaDel) + diag(betaDel,1) + diag(betaDel,-1) ...
%     + diag(gammaDel,2) + diag(gammaDel,-2);
% 
%   dDeltemp = -svec*C'*(dDelC)*(C'\svec')/pi;
%   dDelabsorption(n) = real(dDeltemp);
% 
%   dDelfrequency(n)  = dDeltempom;
%   dDeltempom   = dDeltempom + 2.5*F/pts;
% end
% 
% if (dDeltempom < om)
%   dDelfrequency = dDelfrequency(end:-1:1,:);
%   dDelabsorption = dDelabsorption(end:-1:1,:);
% end
% 
% 

