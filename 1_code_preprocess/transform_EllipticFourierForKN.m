%% Phase 4.1.b Shape descriptors for KNN
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   
%   ======================================================================

% Plot Elliptic Fourier for KNN

No = [200 201 20 21 650 651];
image = cell(1,numel(No));
for i = 1:numel(No)
    n = No(i);
    canvas = cat(3, dynamicImages{n,1}, dynamicImages{n,1}, dynamicImages{n,1});
    canvas(:,:,1) = canvas(:,:,1) + dynamicImages{n,2};
    canvas(:,:,2) = canvas(:,:,2) - dynamicImages{n,2};
    canvas(:,:,3) = canvas(:,:,3) - dynamicImages{n,2};
    image{i} = imresize(canvas,2);
end

figure;
subplot(2,3,1), imshow(image{1})
subplot(2,3,2), imshow(image{2})
subplot(2,3,3), imshow(image{3})
subplot(2,3,4), imshow(image{4})
subplot(2,3,5), imshow(image{5})
subplot(2,3,6), imshow(image{6})

cell1 = fourierDescriptorKNN{No(1),1:end-2};
cell2 = fourierDescriptorKNN{No(2),1:end-2};
cell3 = fourierDescriptorKNN{No(3),1:end-2};
cell4 = fourierDescriptorKNN{No(4),1:end-2};
cell5 = fourierDescriptorKNN{No(5),1:end-2};
cell6 = fourierDescriptorKNN{No(6),1:end-2};

% Compute Pairwise Euclidean distance between morphologies
dist = pdist([cell1,cell2,cell3,cell4,cell5,cell6]','euclidean');
distSqrF = squareform(dist);

%visualize results
cmap = [0.2622,    0.1804,    0.7603
        0.2803,    0.2782,    0.9221
        0.2440,    0.3358,    0.9988
        0.2120,    0.6314,    0.6901
        0.1709,    0.7154,    0.6342
        0.1938,    0.7758,    0.6251
        0.5044,    0.7993,    0.3480
        0.8634,    0.7406,    0.1596
        0.9892,    0.8136,    0.1885];   
figure;
imagesc(distSqrF);
colormap(cmap);
txt = num2str(distSqrF(:),'%0.2f');
txt = strtrim(cellstr(txt));
[x,y] = meshgrid(1:6);
hstr  = text(x(:),y(:),txt(:),'HorizontalAlignment','center','color','white');
names = {'cell 1','cell 2','cell 3','cell 4','cell 5','cell 6'};
set(gca,'FontSize', 16,'XTickLabel',names,'YTickLabel',names)

%% Boundary for KNN

No = [200 201 20 21 650 651];
image = cell(1,numel(No));
for i = 1:numel(No)
    n = No(i);
    canvas = cat(3, dynamicImages{n,1}, dynamicImages{n,1}, dynamicImages{n,1});
    canvas(:,:,1) = canvas(:,:,1) + dynamicImages{n,2};
    canvas(:,:,2) = canvas(:,:,2) - dynamicImages{n,2};
    canvas(:,:,3) = canvas(:,:,3) - dynamicImages{n,2};
    image{i} = imresize(canvas,2);
end

figure;
subplot(2,3,1), imshow(image{1})
subplot(2,3,2), imshow(image{2})
subplot(2,3,3), imshow(image{3})
subplot(2,3,4), imshow(image{4})
subplot(2,3,5), imshow(image{5})
subplot(2,3,6), imshow(image{6})

cell1 = BoundaryDescriptors{No(1),1:end-1};
cell2 = BoundaryDescriptors{No(2),1:end-1};
cell3 = BoundaryDescriptors{No(3),1:end-1};
cell4 = BoundaryDescriptors{No(4),1:end-1};
cell5 = BoundaryDescriptors{No(5),1:end-1};
cell6 = BoundaryDescriptors{No(6),1:end-1};

% Compute Pairwise Euclidean distance between morphologies
dist = pdist([cell1,cell2,cell3,cell4,cell5,cell6]','euclidean');
distSqrF = squareform(dist);

%visualize results
cmap = [0.2622,    0.1804,    0.7603
        0.2803,    0.2782,    0.9221
        0.2440,    0.3358,    0.9988
        0.2120,    0.6314,    0.6901
        0.1709,    0.7154,    0.6342
        0.1938,    0.7758,    0.6251
        0.5044,    0.7993,    0.3480
        0.8634,    0.7406,    0.1596
        0.9892,    0.8136,    0.1885];   
figure;
imagesc(distSqrF);
colormap(cmap);
txt = num2str(distSqrF(:),'%0.2f');
txt = strtrim(cellstr(txt));
[x,y] = meshgrid(1:6);
hstr  = text(x(:),y(:),txt(:),'HorizontalAlignment','center','color','white');
names = {'cell 1','cell 2','cell 3','cell 4','cell 5','cell 6'};
set(gca,'FontSize', 16,'XTickLabel',names,'YTickLabel',names)

