function [image] = inverseFourierDescriptor(a,b,c,d,T,s,ImgSize)
% receives the fourier coefficients a,b,c,d; the estimated contour length
% T; and the step size s.
    
    % Receive the fourier functions
    xt = ellipticFourierFunction(a,b,T);
    yt = ellipticFourierFunction(c,d,T);
    
    k = 1;
    for it = 0:s:T
        pl(k) = xt(it) + 1i*yt(it);
        k = k + 1;
    end
    addX = ceil(ImgSize(2)/2);
    addY = floor(ImgSize(1)/2);
    scaleF = 60;
    x2 = ceil(real(pl)*scaleF)';
    y2 = floor(imag(pl)*scaleF)';
    ind2 = sub2ind(ImgSize, y2+addY-3, x2+addX);
    image = false(ImgSize);
    image(ind2) = true;
end