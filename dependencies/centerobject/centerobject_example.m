%% Example to show the functionality of the centerobject-function
% Green cross marks the centroid of the object
% Red line shows the distance between objects centroid and image center
%
% Your need to have the carriage-17.gif and the centerobject.m file in the
% same folder as the example-script.
%
% The image "carriage-17.gif" is from the MPEG-7 dataset
%
% by Frederik Kratzert 18. Aug 2015
% contact f.kratzert(at)gmail.com

%% Create a shape and get image/object information

    bw = imread('carriage-17.gif');
    bw = double(bw);
    bw = im2bw(bw);
    
    sz_bw = size(bw);
    state_bw = regionprops(bw,'Centroid');
    
%% call centerobject function and get image/object information of result

    cm = centerobject(bw);
    sz_cm = size(cm);
    state_cm = regionprops(cm,'Centroid');
    
%% Plot Results

    figure();
    %input image
    subplot(1,2,1);
    imshow(bw);
    hold on;
    plot(state_bw.Centroid(1),state_bw.Centroid(2),'+','color','green');
    line([sz_bw(2)/2 state_bw.Centroid(1)],...
        [sz_bw(1)/2 state_bw.Centroid(2)],...
        'color','red');
    title(['Image size: ',num2str(sz_bw(2)),...
        'px x ',num2str(sz_bw(1)),'px']);
    hold off;
    
    %output image
    subplot(1,2,2);
    imshow(cm);
    hold on;
    plot(state_cm.Centroid(1),state_cm.Centroid(2),'+','color','green');
    line([sz_cm(2)/2 state_cm.Centroid(1)],...
        [sz_cm(1)/2 state_cm.Centroid(2)],...
        'color','red');
    title(['Image size: ',num2str(sz_cm(2)),...
        'px x ',num2str(sz_cm(1)),'px']);
    hold off;
    
