function cm = centerobject(bw)
%CENTEROBJECT is a function which centers an object in a binary images to
%the image center.
% CM = CENTEROBJECT(BW) gets and MxN binary image BW as input and returns
% and RxR binary image CM, where the object in BW is centered to the image
% center. Due to numerical reasons, there can be an error of max +-0.5px.
% A requirement is, that the input image must contain only one object.
%
% by Frederik Kratzert 18. Aug 2015
% contact f.kratzert(at)gmail.com
%

if ~islogical(bw)
    error('Input BW must be of type ''logical''');
end

cc = bwconncomp(bw,8);

if cc.NumObjects > 1
    error('There is more than one object in the binary image');
end

%adds coloums/rows of zeros to create a square image
sz = size(bw);
if sz(1) < sz(2)
    temp = (sz(2) - sz(1))*0.5;
    if mod(temp,1) > 0
        bw = padarray(bw,[0 1], 'post');
    end
    bw = padarray(bw,[round(temp) 0]);
elseif sz(2) < sz(1)
    temp = (sz(1) - sz(2))*0.5;
    if mod(temp,1) > 0
        bw = padarray(bw, [1 0], 'post');
    end
    bw = padarray(bw,[0 round(temp)]);
end

%gets translation factors in x and y direction
state = regionprops(bw,'Centroid');
sz = size(bw);
delta_y = round(sz(1)/2-state.Centroid(2));
delta_x = round(sz(2)/2-state.Centroid(1));

%extend image, so that translation fit in any case in the boundarys
delta_max = max(abs(delta_y),abs(delta_x));
bw = padarray(bw,[delta_max+10, delta_max+10]);

%execute translation
cm = circshift(bw,[delta_y,delta_x]);

end

