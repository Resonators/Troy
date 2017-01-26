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
%  Computes the sum over products of binomial coefficients appearing in Wei's
%  expression for Wigner 9j symbols.  The arguments are represented by the
%  tablesu
% 
%  [ m1 m2 m3 ]
%  [ m4 m5 m6 ]
% 
%  subject to the constraints
% 
%  m1 + m2 + m3, m1 + m5 + m6, m2 + m4 + m6, m3 + m4 + m5 are all integers
% 
%  Usage:
% 
%  [binprod, flag] = sixbin (m1, m2, m3, m4, m5, m6)
% 
%  where binprod is the required sum over binomial coefficients and flag
%  reports the error status: successful = true, unsuccessful = false.  In the
%  case of an unsuccessful computation, sixbin returns binprod = 0.  Note that
%  a value of binprod = 0 is not necessarily an error.
% 
%  Literature:
% 
%  Liqiang Wei. "New formula for 9-j symbols and their direct calculation" 
%  Computers in Physics, 12(6):632 -- 634 (1998)
% 
%  Calls:    recbin
%  Author:   Keith Earle <kearle@albany.edu>
%  Created:  2016-04-22
%  Modified:
%  Future:
% 

function [retval, nerror] = sixbin (m1, m2, m3, m4, m5, m6)

  m156 = (floor (m1 + m5 + m6) == m1 + m5 + m6 );
  m246 = (floor (m2 + m4 + m6) == m2 + m4 + m6 );
  m345 = (floor (m3 + m4 + m5) == m3 + m4 + m5 );
  m123 = (floor (m1 + m2 + m3) == m1 + m2 + m3 );
  
  mflag = m156 & m246 & m345 & m123;
  
  p     = max ([m1+m5+m6,m2+m4+m6,m3+m4+m5,m1+m2+m3]);
  q     = min ([m1+m2+m4+m5,m1+m3+m4+m6,m2+m3+m5+m6]);
  
  iflag = (floor ( p + q ) == ( p + q ));
  oflag = (p <= q);
  
  nerror = mflag & iflag & oflag;
  temp   = 0;
  
  if (nerror)
    for n = p:1:q
      temp = (temp + ((-1)^n)*recbin(n+1,n-m1-m5-m6)*recbin(m1+m5-m6,n-m2-m4-m6)*recbin(m1-m5+m6,n-m3-m4-m5)*recbin(-m1+m5+m6,n-m1-m2-m3));
    end
    end
  retval = temp;
    end
