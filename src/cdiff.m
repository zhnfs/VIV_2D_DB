%=========================================
% cdiff.m:  2nd-order centered differences
%=========================================

%================================================
% Copyright Franz Hover and Michael Triantafyllou
% MIT/AMT November 2001
%================================================

function y = cdiff(x) ;

n = length(x) ;
[r,c] = size(x) ;
if r > c,
	trans = 1 ;
else,
	trans = 0 ;
end;

y(1) = -3*x(1) + 4*x(2) - x(3) ;
y(n) =  x(n-2) - 4*x(n-1) + 3*x(n) ;
y(2:n-1) = x(3:n) - x(1:n-2) ;

y = y/2 ;

if trans == 1,
	y = y' ;
end;
