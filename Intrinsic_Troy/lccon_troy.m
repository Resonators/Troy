% Compute the Levi-Civita Connection
% run defaults, stv, mat, spec, dderiv, metric


l1 = [pD./absorption; pDel./absorption; pLam./absorption];
l2 = [pDD./absorption; pDDel./absorption; pDLam./absorption; ...
      pDDel./absorption; pDelDel./absorption; pDelLam./absorption; ...
      pDLam./absorption; pDelLam./absorption; pLamLam./absorption];

Gamma = [];

for l = 0:2
   for m = 0:2
      for n = 0:2
        temp = trapz(frequency, ...
               (l2(3*m+n+1)-0.5*l1(m+1).*l1(n+1)).*l1(l+1).*absorption);
        Gamma(9*l+3*m+n+1) = temp;
      end
   end
end
