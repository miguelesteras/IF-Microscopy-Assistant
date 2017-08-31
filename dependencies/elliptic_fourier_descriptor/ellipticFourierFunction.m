function xt = ellipticFourierFunction(a,b,T)
% Returns the elliptic fourier polynomial as function. a and b are the
% coefficients for cos and sin (in this very order), T is the period
% Bibliography:
% F.P. Kuhl and C.R. Giadina, Elliptic Fourier features of a closed
% contour, Comput. Vis. Graphics Image Process. 18, 236-258 (1982)
% W. Dahmen and A. Reusken, Numerik fuer Ingenieure und 
% Naturwissenschaftler, 2. korrigerte Auflage, Springer Verlag, Aachen 2007
% A. Angermann, M. Beuschel, M. Rau and U. Wohlfarth, MATLAB - Simulink - 
% Stateflow, 7. Auflage, Oldenbourg Verlag, M?nchen 2011
% Remark: 
% - We don't use a fast fourier transform. Hence, the calculation will take
%   O(N)
% - There are no parameter checks
%
% Copyright (c) 2013 by Tobias Pohlen (tobias.pohlen@gmail.com)
% www.geekstack.net
    N = size(a);
    xt = @(x) evalX(x,a,b,N,T);
end

function y = evalX(x,a,b,N,T)
% Returns the value of the polynomial at position x
    y = 0;
    for it = 1:N
        y = y + a(it)*cos(2*it*pi*x/T) + b(it) * sin(2*it*pi*x/T);
    end
end