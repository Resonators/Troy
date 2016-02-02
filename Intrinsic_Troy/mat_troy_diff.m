% Set magnetic inhomogeneity diffusion rate potential coefficient and maxl in
% octave then call mat.  There is a script defaults.m that does this.
%

% Tau=0;
% p=0;
Del = Delta;
D   = Diff;
Lam = lambda;
Lmx = maxl;
T2i = T2inv;
plflag = 0;

% Initialize temporary arrays

tempd       = [];
tempd1      = [];
tempd2      = [];

% Initialize arrays for diagonal elements and their derivatives

alpha       = [];
alphaD      = [];
alphaDLam   = [];
alphaDel    = [];
alphaLam    = [];
alphaLamLam = [];
alphap      = [];
alphaTau    = [];

% Initialize arrays for superdiagonal elements and their derivatives

beta        = [];
betaD       = [];
betaDLam    = [];
betaDel     = [];
betaLam     = [];
betaLamLam  = [];

% Initialize arrays for superduperdiagonal elements and their derivatives

gamma       = [];
gammaD      = [];
gammaDLam   = [];
gammaLam    = [];
gammaLamLam = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute Diagonal matrix elements of the secular g tensor problem in an axial
% orienting potential and their derivatives.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%


for n = 0:ceil(Lmx/2.)
  tempd(n+1) = (2*n*(2*n+1)*D)/(1+2*n*(2*n+1)*D*Tau)^p;
  tempd(n+1) = tempd(n+1) + 1i*(Del/6)*n*(n+0.5)/((n+0.75)*(n-0.25));
  tempd(n+1) = tempd(n+1) + T2i;
%   tempd(n+1) = tempd(n+1) + 0.3*D*Lam^2;
%   tempd(n+1) = tempd(n+1) ...
%                -(4*n+1)*3*Lam*D*(1-Lam/14)*wigner3j(4*n,4,4*n,0,0,0)^2 ...
%                -(4*n+1)*(18/35)*Lam^2*D*wigner3j(4*n,8,4*n,0,0,0)^2;
end
alpha    = tempd;

for n = 0:ceil(Lmx/2.)
  tempd(n+1) = -2*n*(2*n+1)*D*log(1+2*n*(2*n+1)*D*Tau)/(1+2*n*(2*n+1)*D*Tau)^p;
end
alphap   = tempd;

for n = 0:ceil(Lmx/2.)
  tempd(n+1) = -p*(2*n*(2*n+1)*D)^2/(1+2*n*(2*n+1)*D*Tau)^(p+1);
end
alphaTau   = tempd;

for n = 0:ceil(Lmx/2.)
  tempd(n+1) = (1i/6)*n*(n+0.75)/((n+0.75)*(n-0.25));
end
alphaDel   = tempd;

for n = 0:ceil(Lmx/2.)
  tempd(n+1) = [2*n*(2*n+1)*(1+2*n*(2*n+1)*D*Tau)^(-p)+(2*n*(2*n+1)*D*(-p*(1+2*n*(2*n+1)*D*Tau)^(-p-1)*2*n*(2*n+1)*Tau))];
%   +0.3*Lam^2 - (4*n+1)* ...
%               ( 3*Lam*(1-Lam/14)*wigner3j(4*n,4,4*n,0,0,0)^2 ...
% 	       + (18/35)*Lam^2*wigner3j(4*n,8,4*n,0,0,0)^2);
end
alphaD     = tempd;
% 
% for n = 0:ceil(Lmx/2.)
%   tempd(n+1) = 0.6*Lam - (4*n+1)* ...
%               ( 3*(1-Lam/7)*wigner3j(4*n,4,4*n,0,0,0)^2 ...
% 	       + (36/35)*Lam*wigner3j(4*n,8,4*n,0,0,0)^2);
% end
 alphaDLam  = zeros(1,length(tempd));

% for n = 0:ceil(Lmx/2.)
%   tempd(n+1) = 0.6*Lam*D - (4*n+1)*D* ...
% 	     ( 3*(1-Lam/7)*wigner3j(4*n,4,4*n,0,0,0)^2 ...
% 	      +(36/35)*Lam*wigner3j(4*n,8,4*n,0,0,0)^2);
% end
alphaLam   = zeros(1,length(tempd));
% 
% for n = 0:ceil(Lmx/2.)
%   tempd(n+1) = 0.6*D - (4*n+1)*D* ...
% 	     (-(3/7)*wigner3j(4*n,4,4*n,0,0,0)^2 ...
% 	       + (36/35)*wigner3j(4*n,8,4*n,0,0,0)^2);
% end
 alphaLamLam   = zeros(1,length(tempd));

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute superdiagonal matrix elements and their derivatives
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

for n = 1:ceil(Lmx/2.)
  tempd1(n) = 1i*(Del/4)*n*(n-0.5)/(n-0.25);
  tempd1(n) = tempd1(n)/((n+0.25)*(n-0.75))^(0.5);
  tempd1(n) = tempd1(n) ...
             -sqrt((4*n-3)*(4*n+1))*D* ...
              (3*Lam*(1-Lam/14)*wigner3j(4*n-4,4,4*n,0,0,0)^2 ...
              +(18/35)*Lam^2*wigner3j(4*n-4,8,4*n,0,0,0)^2);
end
beta        = tempd1;

for n = 1:ceil(Lmx/2.)
  tempd1(n) = 1i*(1/4)*n*(n-0.5)/(n-0.25);
  tempd1(n) = tempd1(n)/((n+0.25)*(n-0.75))^(0.5);
end
betaDel     = tempd1;

for n = 1:ceil(Lmx/2.)
  tempd1(n) = -sqrt((4*n-3)*(4*n+1))* ...
               (3*Lam*(1-Lam/14)*wigner3j(4*n-4,4,4*n,0,0,0)^2 ...
              +(18/35)*Lam^2*wigner3j(4*n-4,8,4*n,0,0,0)^2);
end
betaD     = tempd1;

for n = 1:ceil(Lmx/2.)
  tempd1(n) = -sqrt((4*n-3)*(4*n+1))* ...
               (3*(1-Lam/7)*wigner3j(4*n-4,4,4*n,0,0,0)^2 ...
		+(36/35)*Lam*wigner3j(4*n-4,8,4*n,0,0,0)^2);
end
betaDLam    = tempd1;

for n = 1:ceil(Lmx/2.)
  tempd1(n) = -sqrt((4*n-3)*(4*n+1))*D* ...
              (3*(1-Lam/7)*wigner3j(4*n-4,4,4*n,0,0,0)^2 ...
              +(36/35)*Lam*wigner3j(4*n-4,8,4*n,0,0,0)^2);
end
betaLam     = tempd1;

for n = 1:ceil(Lmx/2.)
  tempd1(n) = -sqrt((4*n-3)*(4*n+1))*D* ...
	    (-(3/7)*wigner3j(4*n-4,4,4*n,0,0,0)^2 ...
              +(36/35)^wigner3j(4*n-4,8,4*n,0,0,0)^2);
end
betaLamLam  = tempd1;

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute superduper diagonal matrix elements and their derivatives.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

for n = 2:ceil(Lmx/2.)
  tempd2(n-1) = -sqrt((4*n-7)*(4*n+1))*D*(18/35)*Lam^2 ...
	     * wigner3j(4*n-8,8,4*n,0,0,0)^2;
end
gamma       = tempd2;
gammaD      = tempd2/D;
gammaDLam   = 2*tempd2/(D*Lam);
gammaLam    = 2*tempd2/Lam;
gammaLamLam = 2*gammaLam/(Lam^2);
gammaDel    = zeros(1,length(tempd2));
