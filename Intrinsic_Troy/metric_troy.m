% Compute the Fisher information from the spectral function and its derivatives
% call defaults.m stv.m mat.m spec.m deriv.m before calling metric.m

dlogpD   = pD./absorption;
dlogpDel = pDel./absorption;
dlogpLam = pLam./absorption;

gDD      = trapz(frequency,dlogpD.*dlogpD.*absorption);
gDDel    = trapz(frequency,dlogpD.*dlogpDel.*absorption);
gDLam    = trapz(frequency,dlogpD.*dlogpLam.*absorption);
gDelDel  = trapz(frequency,dlogpDel.*dlogpDel.*absorption);
gDelLam  = trapz(frequency,dlogpDel.*dlogpLam.*absorption);
gLamLam  = trapz(frequency,dlogpLam.*dlogpLam.*absorption);

gcov     = [gDD  , gDDel  , gDLam; ...
            gDDel, gDelDel, gDelLam; ...
            gDLam, gDelLam, gLamLam];

gcon     = inv(gcov);
