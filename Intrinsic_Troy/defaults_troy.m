lambda = 0.1;        % potential parameter
maxl   = 16;       % maximum L value

Diff   = 1e7;      % rotational diffusion parameter
T2inv  = 3e5;      % intrinsic linewidth

% Set magnetic parameters in frequency units

gpar = 2.002;      % parallel component of the g tensor
gprp = 2.007;      % perpendicular component of the g tensor
B0   = 1.2;        % magnetic field in Tesla

pts  = 512;        % number of spectral points to compute

% These variables determin spectral extent and center frequency
% F is the Magnetic anisotropy
% Delta is the frequency anisotropy defined by Moro and Segre
% om0 is the Larmor Frequency.
% om is the Starting frequency for the frequency sweep

F     = (2*(gpar-gprp)/3)*(9.274009e-24/1.054e-34)*B0;  
Delta = (3/2)*F;                                        
om0   = ((gpar+2*gprp)/3)*(9.274009e-24/1.054e-34)*B0;   
om    = om0 - 1.5*F;  
