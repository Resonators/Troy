# needs defaults.m stv.m mat.m

# Initialize data arrays

absorption = [];
frequency  = [];

tempom = om;

for n = 1:pts

  shift(1:size(alpha)(2)) = 1i*(tempom-om0);

  C = spdiag(shift) + spdiag(alpha) + spdiag(beta,1) + spdiag(beta,-1) ...
    + spdiag(gamma,2) + spdiag(gamma,-2);

  temp = svec*(C\svec')/pi;
  absorption(n) = real(temp);

  frequency(n)  = tempom;
  tempom   = tempom + 2.5*F/pts;
endfor

if (tempom < om)
  frequency = swapcols(frequency);
  absorption = swapcols(absorption);
endif
