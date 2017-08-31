function [a,b,c,d,T] = ellipticFourierDescriptor(filename, N, p, rotation)
% Returns the coefficients for the elliptic fourier descriptor
% The image referenced by filename is read and the fist contour is 
% extracted using MATLAB's internal contourc-algorithm. 
% The calculated fourier coefficients are phase-, scale and rotation-
% invariant. 
% Parameters: 
% filename - the filename of the image which is being processed. The image
% is supposed to be black/white and its only content must be the shape
% N - The number of fourier coefficients which are calculated. The higher
% this is, the better the approximation will be. 
% rotation - boolean value: If false, the rotation invariance will not be
% performed
% p - A percentage value in (0,1] which indicates how many of the contour
% points shall be used
% Return values:
% T - The period
% a,b,c,d - The coefficient vectors
%
% Bibliography:
% [1] F.P. Kuhl and C.R. Giadina, Elliptic Fourier features of a closed
% contour, Comput. Vis. Graphics Image Process. 18, 236-258 (1982)
%
% [2] W. Dahmen and A. Reusken, Numerik fuer Ingenieure und 
% Naturwissenschaftler, 2. korrigerte Auflage, Springer Verlag, Aachen 2007
%
% [3] A. Angermann, M. Beuschel, M. Rau and U. Wohlfarth, MATLAB - Simulink
% - Stateflow, 7. Auflage, Oldenbourg Verlag, Muenchen 2011
%
% Remark: 
% - We don't use a fast fourier transform. Hence, the calculation will take
%   O(N^2)
% - There are no parameter checks
%
% Copyright (c) 2013 by Tobias Pohlen (tobias.pohlen@gmail.com)
% www.geekstack.net
    
    % Get the contour
    [conX, conY] = getContour(filename);
    
    % Get the total number of points
    points = size(conX, 1) - 1;
    
    % Determine the steps size
    h = 1/p;
    
    x = conX(floor(1:h:points));
    y = conY(floor(1:h:points));
    
    m = size(x, 1);
    
    % Calculate the delta t's and the total period
    dt = calcDeltaT(x,y);
    T = sum(dt);
    
    % Calculate the coefficiencts
    % Please note: This is not a fast fourier transform!
    a = zeros(N,1);
    b = zeros(N,1);
    c = zeros(N,1);
    d = zeros(N,1);
    for it = 1:N
        a(it) = coeffA(it,m,T,dt,x,y);
        b(it) = coeffB(it,m,T,dt,x,y);
        c(it) = coeffC(it,m,T,dt,x,y);
        d(it) = coeffD(it,m,T,dt,x,y);
    end
    
    
    % Normalize the coefficients
    % --------------------------
    
    % 1. Phase invariance
    % This is needed to receive the same fourier descriptor no matter which
    % order of points you pick
    [a,b,c,d] = normalizePhase(a,b,c,d,N);
    
    % Extract the scale-invariance to use it later
    E = sqrt(a(1)^2+c(1)^2);
    
    % 2. Rotation invariance
    if rotation
        [a,b,c,d] = normalizeRotation(a,b,c,d,N);
    end
    
    % 3. Scale invariance
    a = a./E;
    b = b./E;
    c = c./E;
    d = d./E;
end

function [x,y] = getContour(image)
% Returns two vectors x,y which contain the outer contour of the object
    
    % Read the image and extract the contour using MATLAB's algorithms
    img = double(image);
    con = contourc(img);
    
    % Extract the first contour
    l = con(2,1);
    
    x = con(1, 2:l+1)';
    y = con(2, 2:l+1)';
end

function [a,b,c,d] = normalizePhase(a,b,c,d,N)
% Normalizes the phase according to [1]
    delta = 0.5*atan(2*(a(1)*b(1)+c(1)*d(1))/sqrt(a(1)^2+b(1)^2+c(1)^2+d(1)^2));
    for it = 1:N
        M = [a(it), b(it);c(it), d(it)];
        K = [cos(it*delta), -sin(it*delta);sin(it*delta), cos(it*delta)];
        B = M*K;
        a(it) = B(1,1);
        b(it) = B(1,2);
        c(it) = B(2,1);
        d(it) = B(2,2);
    end
end

function [a,b,c,d] = normalizeRotation(a,b,c,d,N)
% Normalizes the phase according to [1]
    psi = atan(c(1)/a(1));
    K = [cos(psi), sin(psi);-sin(psi), cos(psi)];
    for it = 1:N
        M = [a(it), b(it);c(it), d(it)];
        B = K*M;
        a(it) = B(1,1);
        b(it) = B(1,2);
        c(it) = B(2,1);
        d(it) = B(2,2);
    end
end

function dt = calcDeltaT(x,y)
% Calcualtes the delta t_i according to [1]
    n = size(x, 1);
    dt = zeros(n,1);
    for jt = 1:n
        jt1 = mod(jt-2, n) + 1;
        deltaX = x(jt) - x(jt1);
        deltaY = y(jt) - y(jt1);
        dt(jt) = sqrt(deltaX^2 + deltaY^2);
    end
end

function a = coeffA(n, m, T, dt, x, y)
% Calculates the coefficient a_n
    a = 0;
    for jt = 1:m
        deltaX = x(jt) - x(mod(jt-2, m)+1);
        a = a + deltaX/dt(jt) * (cos(2*n*pi*ti(jt, dt)/T) - cos(2*n*pi*ti(jt-1, dt)/T));
    end
    a = a*T/(2*n^2*pi^2);
end

function b = coeffB(n, m, T, dt, x, y)
% Calculates the coefficient b_n
    b = 0;
    for jt = 2:m
        deltaX = x(jt) - x(mod(jt-2, m)+1);
        b = b + deltaX/dt(jt) * (sin(2*n*pi*ti(jt, dt)/T) - sin(2*n*pi*ti(jt-1, dt)/T));
    end
    b = b*T/(2*n^2*pi^2);
end

function c = coeffC(n, m, T, dt, x, y)
% Calculates the coefficient c_n
    c = 0;
    for jt = 1:m
        deltaY = y(jt) - y(mod(jt-2, m)+1);
        c = c + deltaY/dt(jt) * (cos(2*n*pi*ti(jt, dt)/T) - cos(2*n*pi*ti(jt-1, dt)/T));
    end
    c = c*T/(2*n^2*pi^2);
end

function d = coeffD(n, m, T, dt, x, y)
% Calculates the coefficient d_n
    d = 0;
    for jt = 1:m
        deltaY = y(jt) - y(mod(jt-2, m)+1);
        d = d + deltaY/dt(jt) * (sin(2*n*pi*ti(jt, dt)/T) - sin(2*n*pi*ti(jt-1, dt)/T));
    end
    d = d*T/(2*n^2*pi^2);
end


function t = ti(m, dt)
    t = 0;
    for it = 1:m
        t = t + dt(it);
    end
end