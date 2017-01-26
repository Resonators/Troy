%  version 1.0 22 Apr 16
% -------------------------------------------------------------------------------
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
%  Function: [bincof, flag] = recbin ( n, k )
% 
%  Computes the binomial coefficient
% 
%  ( n )
%  (   )
%  { k }
% 
%  recursively and returns the result in bincof.  Recursion is used to provide 
%  stability against overflow. A successful calculation returns flag = true.
% 
%  The caclculation is based on the recurrence relation
% 
%  ( n )     (   n   ) (n + 1) - k
%  (   )  =  (       ) ------------
%  ( k )     ( k - 1 )      k
% 
%  which is readily derived from the definition of the binomial coefficient
% 
%  ( n )        n!
%  (   ) =  -----------
%  ( k )    k! (n - k)!
% 
%  Usage: n <= k are non-negative integers. 
%  Literature: Liqiang Wei. "New formula for 9-j symbols and their direct caclcu-
%  lation" Computers in Physics, 12(6) 632 -- 634 (1998).
%  Uses:
%  Calls:
%  Modified:
%  Author:   Keith Earle <kearle@albany.edu>
% 

function [bincof, flag] = recbin (n, k)
  
  bincof = 0.;
  flag   = false;
  nflag  = false;
  iflag  = false;
  rflag  = false;
  
  if ((k >= 0) & (n >= 0))
    nflag = true;
  end
  
  if (( floor (k) == k ) & ( floor (n) == n ))
    iflag = true;
  end
  
  if ( n >= k )
    rflag = true;
  end
  
  flag = nflag & iflag & rflag;
  temp = 1;
  m    = 1;
  
  if (flag)
    while ( m <= k )
      temp = temp * (n+1-m)/m;
      m=m+1;
    end
  end
  
  bincof = temp;
  end
