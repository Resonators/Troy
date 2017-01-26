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
%  GNU General Public Licenset for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
%  function: retval = Triangle (a,b,c)
% 
%  Tests the triangle inequalities in Messiah, A. "Quantum Mechanics" v. 2,
%  pg. 1062 (North-Holland, Amsterdam) 1962.  Returns true if the inequalities
%  are satisfied, false if not.  Also tests if the triad sums to an integer
% 
%  Author:   Keith Earle <kearle@albany.edu>
%  Ported:   Octave 17 Apr 16
% 

function retval = triangle(x,y,z)

    retval=logical((abs(x-y) <= z) && (z <= x+y) && (floor(x+y+z) == x+y+z));

    end

