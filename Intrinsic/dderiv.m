# Compute derivatives with respect to the parameters of the spectral function
# call defaults.m stv.m mat.m spec.m before using deriv.m

# Initialize data arrays

pD         = [];
pDel       = [];
pLam       = [];
pDD        = [];
pDDel      = [];
pDLam      = [];
pDelLam    = [];
pLamLam    = [];
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

CLam  =   spdiag(alphaLam) + spdiag(betaLam,1) + spdiag(betaLam,-1) ...
        + spdiag(gammaLam,2) + spdiag(gammaLam,-2);
CLam  =   CLam/om0;

CDLam =   spdiag(alphaDLam) + spdiag(betaDLam,1) + spdiag(betaDLam,-1) ...
        + spdiag(gammaDLam,2) + spdiag(gammaDLam,-2);
CDLam =   CDLam/om0;

CLamLam =  spdiag(alphaLamLam)+spdiag(betaLamLam,1)+spdiag(betaLamLam,-1) ...
          +spdiag(gammaLamLam,2) + spdiag(gammaLamLam,-2);

CLamLam = CLamLam/(om0*om0);

for n = 1:pts

  shift(1:size(alpha)(2)) = 1i*(tempom-om0);

  C = spdiag(shift) + spdiag(alpha) + spdiag(beta,1) + spdiag(beta,-1) ...
    + spdiag(gamma,2) + spdiag(gamma,-2);

# Magic is necessary to scale lambda derivatives to similar values as
# dynamic derivatives.
  tempv      = C\svec';
  tempd      = (C\dvec')/om0;
  tempdd     = (C\ddvec')/(om0*om0);
  wLam       = svec + (dvec/om0);
  tempw      = C\wLam';
  wLamLam    = svec + ddvec/(om0*om0);
  tempww     = C\wLamLam';

  pD(n)      = -real(tempv'*CD*tempv)/pi;
  pDel(n)    = -real(tempv'*CDel*tempv)/pi;
  pLam(n)    =  real(wLam*tempw)/pi - real(tempv'*CLam*tempv)/pi ...
              - real(svec*tempv)/pi - real(dvec*tempd)/pi;

# Now we evaluate second derivatives

  tDv        = C\(CD*tempv);
  tDelv      = C\(CDel*tempv);
  tDDelv     = C\(CDel*tempv + CD*tempv);
  tDelLamv   = C\(CDel*tempv + CLam*tempv);
  tDLamv     = C\(CD*tempv + CLam*tempv);
  tLamv      = C\(CLam*tempv);

  pDD(n)     =  2*real((CD*tempv)'*tDv)/pi;
  pDelDel(n) =  2*real((CDel*tempv)'*tDelv)/pi;
  pDDel(n)   =    real((CDel*tempv + CD*tempv)'*tDDelv)/pi ...
              - 0.5*(pDD(n)+pDelDel(n));

# The second derivatives involving Lambda are more complicated

  pDelLam(n) = -real(wLam*CDel*wLam')/pi +real(tempd'*CDel*tempd)/pi ...
              + real(tempv'*CDel*tempv)/pi ...
              + real((CDel*tempv + CLam*tempv)'*tDelLamv)/pi ...
              - real((CDel*tempv)'*tDelv)/pi ...
              - real((CLam*tempv)'*tLamv)/pi;

  pDLam(n)   = -real(wLam*CD*wLam')/pi +real(tempd'*CD*tempd)/pi ...
              + real(tempv'*CD*tempv)/pi ...
              + real((CD*tempv + CLam*tempv)'*tDLamv)/pi ...
              - real((CD*tempv)'*tDv)/pi ...
              - real((CLam*tempv)'*tLamv)/pi ...
              - real(tempv'*CLamLam*tempv)/pi;

  pLamLam(n) =    real(wLamLam*tempww)/pi ...
              -   absorption(n) ...
              -   real((ddvec/(om0*om0))*tempdd)/pi ...
              + 2*real((dvec/om0)*tempd)/pi ...
              - 2*real(tempw'*CLam*tempw)/pi ...
              + 2*real(tempd'*CLam*tempd)/pi ...
              + 2*real((CLam*tempv)'*tLamv)/pi ...
              -   real(tempv'*CLamLam*tempv)/pi;                

  frequency(n)  = tempom;
  tempom   = tempom + 2.2*F/pts;
endfor

if (tempom < om)
  frequency = swapcols(frequency);
  pD        = swapcols(pD);
  pDel      = swapcols(pDel);
  pLam      = swapcols(pLam);
  pDD       = swapcols(pDD);
  pDelDel   = swapcols(pDelDel);
  pDDel     = swapcols(pDDel);
  pDelLam   = swapcols(pDelLam);
  pLamLam   = swapcols(pLamLam);
endif
