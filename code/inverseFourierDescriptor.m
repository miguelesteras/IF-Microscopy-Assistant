function [image] = inverseFourierDescriptor(a,b,c,d,T,s)
% receives the fourier coefficients a,b,c,d; the estimated contour length
% T; and the step size s.
    
    % Receive the fourier functions
    xt = ellipticFourierFunction(a,b,T);
    yt = ellipticFourierFunction(c,d,T);
    
    k = 1;
    for it = 0:s:T
        pl(k) = xt(it) + i*yt(it);
        k = k + 1;
    end
    % plot reconstructed image with center [100,100] and rescale by scaleF.
    scaleF = 50;
    x2 = round(real(pl)*scaleF)';
    y2 = round(imag(pl)*scaleF)';
    ind2 = sub2ind([200,200], y2+100, x2+100);
    image = false(200,200);
    image(ind2) = true;
end