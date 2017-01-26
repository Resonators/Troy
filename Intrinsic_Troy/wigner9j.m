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
## Computes the Wigner 9j symbol with arguments arranged in the following
## tableau
##
## { a b c }
## { d e f }
## { g h j }
##
## subject to the following constraints on the following triples
##
## (a,b,c), (d,e,f), (g,h,j), (a,d,g), (b,e,h), (c,f,j) satisfy the triangle
## conditions tested in Triangle (j1, j2, j3).
##
## Usage:
##
## [W9j, status] = Wig9j ( a, b, c, d, e, f, g, h. j )
##
## The value of the recoupling coefficient is returned in W9j. If the
## computation executes properly, status = true. If the computation executes
## improperly, W9j = 0 and status = false.  Note that W9j = 0 is not necessarily
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
## Author:   Keith Earle <keith@keith-Dell-System-XPS-L322X>
## Created:  2016-04-23
## Calls:    Delta, sixbin, Triangle
## Modified:
## Future:

function [retval, nerror] = Wig9j ( a, b, c, d, e, f, g, h, j )

  abc = Triangle (a, b, c);
  def = Triangle (d, e, f);
  ghj = Triangle (g, h, j);
  
  adg = Triangle (a, d, g);
  beh = Triangle (b, e, h);
  cfj = Triangle (c, f, j);
  
  tflag = ( abc & def & ghj & adg & beh & cfj );
  
  l1 = max ([ abs(h-d), abs(b-f), abs(a-j)]);
  l2 = min ([     h+d ,     b+f ,     a+j ]);
  
  lflag = (l1 <= l2);
  
  nerror = (tflag & lflag);
  temp   = 0;
  
  if (nerror)
    for k = l1:1:l2
      temp = (temp+((-1)**(2*k))*(2*k+1)*sixbin(a,b,c,f,j,k)*sixbin(f,d,e,h,b,k)*
              sixbin(h,j,g,a,d,k));
    endfor
    temp = (temp * Delta(a,b,c)*Delta(d,e,f)*Delta(g,h,j)*
                   Delta(a,d,g)*Delta(b,e,h)*Delta(c,f,j));
  endif
  retval = temp;
endfunction
