% Compute derivatives with respect to the parameters of the spectral function
% call defaults.m stv.m mat.m spec.m before using deriv.m

% Initialize data arrays

pD         = [];
pDel       = [];
pLam       = [];
frequency  = [];

tempom = om;

CD   =   diag(alphaD) + diag(betaD,1) + diag(betaD,-1) ...
       + diag(gammaD,2) + diag(gammaD,-2);

% Magic is necessary to take account of the fact that the lineshape is the
% real part of the resolvent, but we're taking the derivative with respect
% to a pure imaginary quantity.

CDel =   diag(alphaDel) + diag(betaDel,1) + diag(betaDel,-1);
CDel = -1i*CDel;

% Magic is necessary to scale the derivative so that we don't crash into
% dynamic range issues.

CLam =   diag(alphaLam) + diag(betaLam,1) + diag(betaLam,-1) ...
       + diag(gammaLam,2) + diag(gammaLam,-2);
CLam = CLam/om0;

for n = 1:pts

  shift(1:length(alpha)) = 1i*(tempom-om0);

  C = diag(shift) + diag(alpha) + diag(beta,1) + diag(beta,-1) ...
    + diag(gamma,2) + diag(gamma,-2);

% Magic is necessary to scale lambda derivatives to similar values as
% dynamic derivatives.
  tempv   = C\svec';
  tempd   = (C\dvec')/om0;
  wLam    = svec + (dvec/om0);
  tempw   = C\wLam';

  pD(n)   = -real(tempv'*CD*tempv)/pi;
  pDel(n) = -real(tempv'*CDel*tempv)/pi;
  pLam(n) =  real(wLam*tempw)/pi - real(tempv'*CLam*tempv)/pi ...
           - real(svec*tempv)/pi - real(dvec*tempd)/pi;

  frequency(n)  = tempom;
  tempom   = tempom + 2.2*F/pts;
end

if (tempom < om)
  frequency = frequency(end:-1:1,:) ;
  pD        = pD(end:-1:1,:) ;
  pDel      = swapcols(pDel);
end
