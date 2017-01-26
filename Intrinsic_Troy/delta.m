%  Version 1.0 22 Apr 16
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
%  function [D, check] = Delta (a,b,c)
% 
%  Returns
% 
%  ( a+b-c )! ( a-b+c )! ( -a+b+c )!
%  ---------------------------------
%          ( a + b + c + 1)!
% 
%  in D if (a,b,c) satisfy the triangle constraints:
% 
%  a+b+c must be an integer
%  the sum of any two arguments must not be less than the third.
% 
%  A successful evaluation returns D and check = true.
%  An unsuccessful evaluation returns D = 0 and check = false.
% 
%  Literature: Liqiang Wei. "New formula for 9-j symbols and their direct calcu-
%  lation" Computers in Physics, 12(6):632 -- 634 (1998)
%  Uses:     factorial
%  Calls:    Triangle
%  Author:   Keith Earle <kearle@albany.edu>
%  Modified:
%  Future:


function [retval, nerror] =  delta(a,b,c)
  nerror = triangle(a,b,c);
  D  = 0;
  
  if (nerror==1)
    D = sqrt(factorial(a+b-c)*factorial(a-b+c)*factorial(-a+b+c)/factorial(a+b+c+1));
  end
  retval =  D;
  
end

