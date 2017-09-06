function [vector] = normBoundary(set,vecSize)
% receives a 2-D boundary matrix (m x 2) with coordinates x and y.
% returns the same matrix reduced or augmented to match the size (vecSize x 2).
    
    while size(set,1) ~= vecSize
    % measure euclidean distance between 2 adjacent points and add a new 
    % point between the points with max distance    
    if size(set,1) < vecSize
        v1 = set;
        v2 = [set(2:end,:);set(1,:)];
        distance = sqrt(sum((v1-v2).^2,2));
        [~,idx3] = max(distance);
        if idx3 == size(set,1)
            newPoint = (set(end,:)+set(1,:))./2;
            set = [set;newPoint];
        else
            newPoint = (set(idx3,:)+set(idx3+1,:))./2;
            set = [set(1:idx3,:);newPoint;set(idx3+1:end,:)];
        end
    % measure euclidean distance between 3 adjacent points sequence and delete 
    % the point in the middle of the shortest distance   
    elseif size(set,1) > vecSize
        v1 = set;
        v2 = [set(3:end,:);set(1:2,:)];
        distance = sqrt(sum((v1-v2).^2,2));
        [~,idx3] = min(distance);
        if idx3 == size(set,1)
            set(1,:) = [];
        else
            set(idx3+1,:) = [];
        end
    end
    vector = set;
end
