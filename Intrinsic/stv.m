# Set lambda (the axial potential strength) and maxl in octave then call stv
#

lam = lambda;
max = floor(maxl/2)+1;
vec = [];
dvec = [];
ddvec = [];

# Use Millers algorithm to evaluate the recurrence stably downwards.

vec(1:max+2) = 0.;
vec(max+1)=1.0;
norm = 0.;
for j = max+1:-1:2
  n = 2*(j-1);
  tempn2 = (2*n+1)*(1.-(3.*lam/2.)/((2*n+3)*(2*n-1)))*vec(j);
  tempn4 = (3.*lam/2.)*((n+2)/(2*n+3))*vec(j+1);
  tempn  = (tempn2+tempn4)*2*(2*n-1)/(3*lam*(n-1));
  vec(j-1) = tempn;
  norm   = norm + tempn^2;
endfor

norm = sqrt(norm);
vec = vec/norm;
svec = vec(1:max);

# compute the derivative of the starting vector in terms of the starting
# vector, using the contraction formula for Legendre Polynomials:
#
#   l     m     l + m           k   (l m k)2
#  P (z) P (z) = sum  (2k + 1) P (z)(     )
#               l - m               (0 0 0)
#

# The first derivative is a little easier to deal with.
# Take care of boundary case first

dvec(1:max) = 0;
dvec(1) = vec(2);
for j = 2:max
  n = 2*(j-1);
  dvec(j) = ((2*n+5)*coupling_3j(2*n,4,2*(n+2),0,0,0)^2*vec(j+1) ...
           +(2*n+1)*coupling_3j(2*n,4,2*n,0,0,0)^2*vec(j) ...
           +(2*n-3)*coupling_3j(2*n,4,2*(n-2),0,0,0)^2*vec(j-1))/2;
endfor

# compute the second derivative of the starting vector in terms of the
# starting vector, using the same technique as for the first derivative,
# there are just more terms.
# Take care of boundary cases first.

ddvec(1:max) = 0;
ddvec(1) = (81*coupling_3j(0,8,8,0,0,0)^2*coupling_3j(4,4,8,0,0,0)^2*vec(3) ...
          +25*coupling_3j(0,4,4,0,0,0)^2*coupling_3j(4,4,4,0,0,0)^2*vec(2) ...
	    + 1*coupling_3j(4,4,0,0,0,0)^2*vec(1))/4;

ddvec(2) = (( 13*coupling_3j(4,8,12,0,0,0)^2*vec(4) ...
            + 9*coupling_3j(4,8,8,0,0,0)^2*vec(3) ...
            + 5*coupling_3j(4,8,4,0,0,0)^2*vec(2)) ...
          *9*coupling_3j(4,4,8,0,0,0)^2 ...
         +(9*coupling_3j(4,4,8,0,0,0)^2*vec(3) ...
          +5*coupling_3j(4,4,4,0,0,0)^2*vec(2) ...
          +1*coupling_3j(4,4,0,0,0,0)^2*vec(1)) ...
          *5*coupling_3j(4,4,4,0,0,0)^2 ...
	    +1*coupling_3j(4,4,0,0,0,0)^2*vec(1))/4;

ddvec(3) = (( 17*coupling_3j(8,8,16,0,0,0)^2*vec(5) ...
            +13*coupling_3j(8,8,12,0,0,0)^2*vec(4) ...
            + 9*coupling_3j(8,8,8,0,0,0)^2*vec(3) ...
            + 5*coupling_3j(8,8,4,0,0,0)^2*vec(2)) ...
          *9*coupling_3j(4,4,8,0,0,0)^2 ...
         +(13*coupling_3j(8,4,12,0,0,0)^2*vec(4) ...
           +9*coupling_3j(8,4,8,0,0,0)^2*vec(3) ...
           +5*coupling_3j(8,4,4,0,0,0)^2*vec(2)) ...
           *5*coupling_3j(4,4,4,0,0,0) ...
	    +1*coupling_3j(4,4,0,0,0,0)^2)/4;

ddvec(4) = (( 21*coupling_3j(12,8,20,0,0,0)^2*vec(6) ...
            +17*coupling_3j(12,8,16,0,0,0)^2*vec(5) ...
            +13*coupling_3j(12,8,12,0,0,0)^2*vec(4) ...
            + 9*coupling_3j(12,8,8,0,0,0)^2*vec(3) ...
            + 5*coupling_3j(12,8,4,0,0,0)^2*vec(2)) ...
          *9*coupling_3j(4,4,8,0,0,0)^2 ...
         +(17*coupling_3j(12,4,16,0,0,0)^2*vec(5) ...
	  +13*coupling_3j(12,4,12,0,0,0)^2*vec(4) ...
	  + 9*coupling_3j(12,4,8,0,0,0)^2*vec(3)) ...
           *5*coupling_3j(4,4,4,0,0,0) ...
	    +1*coupling_3j(4,4,0,0,0,0)^2*vec(4))/4;

for j = 5:max
  n = 2*j;
  temp4 =  (2*n+9)*coupling_3j(2*n,8,2*(n+4),0,0,0)^2*vec(j+2) ...
       	  +(2*n+5)*coupling_3j(2*n,8,2*(n+2),0,0,0)^2*vec(j+1) ...
       	  +(2*n+1)*coupling_3j(2*n,8,2*(n+0),0,0,0)^2*vec(j+0) ...
       	  +(2*n-3)*coupling_3j(2*n,8,2*(n-2),0,0,0)^2*vec(j-1) ...
       	  +(2*n-7)*coupling_3j(2*n,8,2*(n-4),0,0,0)^2*vec(j-2);

  temp2 =  (2*n+5)*coupling_3j(2*n,4,2*(n+2),0,0,0)^2*vec(j+1) ...
          +(2*n+1)*coupling_3j(2*n,4,2*(n+0),0,0,0)^2*vec(j+0) ...
          +(2*n-3)*coupling_3j(2*n,4,2*(n-2),0,0,0)^2*vec(j-1);

  temp0 = vec(j);

  ddvec(j) =  (temp4*9*coupling_3j(4,4,8,0,0,0)^2 ...
            + temp2*5*coupling_3j(4,4,4,0,0,0)^2 ...
	     + temp0*1*coupling_3j(4,4,0,0,0,0)^2)/4;
endfor

