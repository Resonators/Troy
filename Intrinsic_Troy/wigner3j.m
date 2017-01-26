%  Version 1.0 23 Apr 16
% ------------------------------------------------------------------------------
%  Copyright (C) 2016 Keith Earle
%  
%  This program is free software; you can redistribute it and/or modify it
%  under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 3 of the License, or
%  (at your option) any later version.
%  
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
%  Computes Wigner 3j symbols according to the algorithm developed by Wei.  The
%  arguments are specified by the tableau
% 
%  ( j1 j2 j3 )
%  ( m1 m2 m3 )
% 
%  subject to the constraints
% 
%  j1 + j2 + j3 is an integer, m1 + m2 + m3 = 0, j1 + j2 >= j3, j2 + j3 >= j1,
%  j3 + j1 >= j2, -j1 <= m1 <= j1, -j2 <= m2 <= j2, -j3 <= m3 <= j3.
% 
%  Usage:
% 
%  [VCC, status] = Wig3j (j1, j2, j3, m1, m2, m3)
% 
%  where VCC is the vector coupling coefficient and the error state is returned
%  in status: a successful computation returns status = true, an unsuccessful
%  computation returns status = 0.
% 
%  Literature:
% 
%  Liqiang Wei. "New formula for 9-j symbols and their direct calculation" 
%  Computers in Physics, 12(6):632 -- 634 (1998)
% 
%  Liqiang Wei. "Unified approach for exact calculation of angular momentum
%  coupling and recoupling coefficients" Computer Physics Communications,
%  120:222 -- 230 (1999).
% 
%  Calls:    Delta, recbin, Triangle
%  Author:   Keith Earle <kearle@albany.edu>
%  Created:  2016-04-23
%  Modified: 
%  Future:
% 

function [retval, nerror] = wigner3j (j1, j2, j3, m1, m2, m3)
  j1pm1 = ( floor ( j1 + m1 ) == j1 + m1);
  j2pm2 = ( floor ( j2 + m2 ) == j2 + m2);
  j3pm3 = ( floor ( j3 + m3 ) == j3 + m3);
  
  m1rng = ( abs ( m1 ) <= j1 );
  m2rng = ( abs ( m2 ) <= j2 );
  m3rng = ( abs ( m3 ) <= j3 );
  
  mflag = ( m1 + m2 + m3 == 0 );
  iflag = ( j1pm1 & j2pm2 & j3pm3 );
  rflag = ( m1rng & m2rng & m3rng );
  tflag = triangle ( j1, j2, j3 );
  
  l1 = max ([0, j1 - j3 + m2, j2 - j3 - m1]);
  l2 = min ([j1 + j2 - j3, j1 - m1, j2 + m2]);
  lflag = (l1 <= l2);
 
  nerror = iflag & lflag & mflag & rflag & tflag;
  temp = 0;
   
  if (nerror)
    for k = l1:1:l2
      temp = (temp + ((-1)^k)*recbin(j1+j2-j3,k)*recbin(j1-j2+j3,j1-m1-k)*recbin(-j1+j2+j3,j2+m2-k));
    end
    temp = (temp*((-1)^(j1-j2-m3))*delta(j1,j2,j3)*sqrt((recbin(2*j1, j1-j2+j3)/recbin(2*j1,j1+m1))*(recbin(2*j2, j1+j2-j3)/recbin(2*j2,j2+m2))*(recbin(2*j3,-j1+j2+j3)/recbin(2*j3,j3+m3))));
  end
  retval = temp;
end
