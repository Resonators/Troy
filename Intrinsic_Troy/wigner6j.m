## Version 1.0 23 Apr 16
##------------------------------------------------------------------------------
## Copyright (C) 2016 Keith Earle
## 
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
## Computes the Wigner 6j symbol with arguments arranged in the following
## tableau
##
## { a b c }
## { d e f }
##
## subject to the following constraints on certain unordered triples
##
## (a,b,c), (c,f,j), (b,f,k), (a,j,k) satisfy the triangle conditions tested
## in Triangle (j1, j2, j3).
##
## Usage:
##
## [W6j, status] = Wig6j ( a, b, c, d, e, f )
##
## The value of the recoupling coefficient is returned in W6j. If the
## computation executes properly, status = true. If the computation executes
## improperly, W6j = 0 and status = false.  Note that W6j = 0 is not necessarily
## an error condition.
##
## Literature:
##
## Liqiang Wei. "New formula for 9-j symbols and their direct calculation" 
## Computers in Physics, 12(6):632 -- 634 (1998)
##
## Liqiang Wei. "Unified approach for exact calculation of angular momentum
## coupling and recoupling coefficients" Computer Physics Communications,
## 120:222 -- 230 (1999).
##
## Author:   Keith Earle <kearle@albany.edu>
## Created:  2016-04-23
## Calls:    Delta, recbin, Triangle
## Modified:
## Future:

function [retval, nerror] = Wig6j (a, b, c, d, e, f)
  abc   = Triangle (a, b, c);
  aef   = Triangle (a, e, f);
  bdf   = Triangle (b, d, f);
  cde   = Triangle (c, d, e);
  tflag = (abc & aef & bdf & cde);
  
  l1    = max ([a+b+c, c+d+e, b+d+f, a+e+f]);
  l2    = min ([a+b+d+e, b+c+e+f, a+c+d+f]);
  lflag = (l1 <= l2);
  
  nerror = (tflag & lflag);
  temp   = 0;
  
  if (nerror)
    for k = l1:1:l2
      temp = (temp + ((-1)**k)*recbin(k+1,k-a-b-c)*recbin(a+b-c,k-c-d-e)*
              recbin(a-b+c,k-b-d-f)*recbin(-a+b+c,k-a-e-f));
    endfor
    temp = temp * Delta(a,e,f)*Delta(b,d,f)*Delta(c,d,e)/Delta(a,b,c);
  endif
  retval = temp;
endfunction

