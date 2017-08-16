%% Phase 4.3. Transform cell masks into a conforming 2-D boundary
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   1. Extract sequences of frames (f1, f2, ... , fx), from which
%   input sequence = (f1, f2, ... , fx-1), and target sequence = fx
%   2. Transform cell masks into a conforming 2-D boundary of fix size
%   given by (x,y) coordenates. The boundary can shrink towards the 
%   interior of the hull. 
%   ======================================================================

%% Build data set of single cells

seqLength = 8;              % Number of frames in learning sequence (x input + 1 target)
vecSize   = 100;            % Boundary descriptor vector size

files = dir('*_metadata.mat');      
num_files = length(files);
for i = 1:num_files
    load(files(i).name,'metadata');                                  
    load(strcat(metadata.name,'_cellSequences.mat'),'cellSequences');   
    load(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');
    load(strcat(metadata.name,'_centerMass.mat'),'centerMass');     
    
    cellCoordenates = cell(size(cellSequences));
    index = find(~cellfun(@isempty,cellSequences)); % non empty cells in cell sequences array
    for h = 1:size(index)    
        selection = cellSequences{index(h)};
        center = centerMass{index(h)};
        image = false(metadata.imageSize);
        image(selection) = true;                % only show selected cell in binary image
        stats = regionprops(image,'PixelList');
        coordenates = stats.PixelList - center;
        [theta, rho] = cart2pol(coordenates(:,1),coordenates(:,2));   % polar coordenates
        cellCoordenates{index(h)} = [theta rho]; 
    end

    % detect sequence of length equal or greater than seqLength
    noFrames = sum(double(~cellfun(@isempty,cellCoordenates)),2);
    idx = find(noFrames >= seqLength);
    BoundaryDataSet = cell(100,seqLength); 
    count = 1;
    
    for j = 1:size(idx,1)
        idx2 = find(~cellfun(@isempty,cellCoordenates(idx(j),:)));
        for k = 1:size(idx2,2) - (seqLength-1)
            rotation = rotationUp{idx(j),idx2(k)} - pi/2;   % -pi/2 added so rotation points west
            input = int8(zeros((metadata.maxRadious*2)+1, (metadata.maxRadious*2)+1));
            for m = 1:seqLength                
                selection = cellCoordenates{idx(j),idx2(k+m-1)};
                selection(:,1) = selection(:,1) + rotation;          % apply rotation
                [x,y] = pol2cart(selection(:,1),selection(:,2));
                cell = [round(x) round(y)];
                
                % Start of Binary Search algorithm implementation for
                % values of the shrinking factor 's'
                from = 0; to = 1; s = 0.5;
                sumInd = boundary(cell,s);
                while numel(sumInd) ~= vecSize
                    if numel(sumInd) < vecSize
                        from = s;
                        s = (s + to)/2;                        
                    elseif numel(sumInd) > vecSize 
                        to = s;
                        s = (s + from)/2; 
                    end
                    sumInd = boundary(cell,s);
                end
                
                boundary = cell(sumInd(1:end-1),:);        % remove last value (because it is equal to first)
                BoundaryDataSet{count,m} = boundary;
            end           
            count = count+1;
        end               
    end
    save(strcat(metadata.name,'_BoundaryDataSet.mat'),'BoundaryDataSet');     
    save(strcat(metadata.name,'_cellCoordenates.mat'),'cellCoordenates');     
    
%     clearvars -except files num_files i
end

%%

figure

ind = sub2ind([200,200], contour(:,2)+100, contour(:,1)+100);
canvas = false(200,200);
canvas(ind) = true;

ind2 = sub2ind([200,200], cell(:,2)+100, cell(:,1)+100);
canvas2 = false(200,200);
canvas2(ind2) = true;

imshowpair(canvas,canvas2,'montage');
% 
% figure
% indSu = sub2ind([193,193], summary(1:50,2), summary(1:50,1));
% canvasSu = false(193,193);
% canvasSu(indSu) = true;
% imshowpair(I,canvasSu,'montage');

% figure
% plot(c(:,1),c(:,2),':','LineWidth',8,'Color',[.9 .9 .9]); hold on;
% scatter(summary(:,1),summary(:,2),75,'filled','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0])
