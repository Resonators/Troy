# Compute derivatives with respect to the parameters of the spectral function
# call defaults.m stv.m mat.m spec.m before using deriv.m

# Initialize data arrays

pD         = [];
pDel       = [];
pLam       = [];
frequency  = [];

tempom = om;

CD   =   spdiag(alphaD) + spdiag(betaD,1) + spdiag(betaD,-1) ...
       + spdiag(gammaD,2) + spdiag(gammaD,-2);

# Magic is necessary to take account of the fact that the lineshape is the
# real part of the resolvent, but we're taking the derivative with respect
# to a pure imaginary quantity.

CDel =   spdiag(alphaDel) + spdiag(betaDel,1) + spdiag(betaDel,-1);
CDel = -1i*CDel;

# Magic is necessary to scale the derivative so that we don't crash into
# dynamic range issues.

CLam =   spdiag(alphaLam) + spdiag(betaLam,1) + spdiag(betaLam,-1) ...
       + spdiag(gammaLam,2) + spdiag(gammaLam,-2);
CLam = CLam/om0;

for n = 1:pts

  shift(1:size(alpha)(2)) = 1i*(tempom-om0);

  C = spdiag(shift) + spdiag(alpha) + spdiag(beta,1) + spdiag(beta,-1) ...
    + spdiag(gamma,2) + spdiag(gamma,-2);

# Magic is necessary to scale lambda derivatives to similar values as
# dynamic derivatives.
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
endfor

if (tempom < om)
  frequency = swapcols(frequency);
  pD        = swapcols(pD);
  pDel      = swapcols(pDel);
endif
