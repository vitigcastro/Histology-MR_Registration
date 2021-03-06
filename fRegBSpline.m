function imageRegBSpline = fRegBSpline(pointsMRI, pointsHistology, imagesMRI, sizeHistolImageGUI, imageHistology, hist2MRI)
% This function registers the MRI images to the histology image, or
% vice-versa, using the non-affine B-spline registration method. 
% The function will estimate the points in one image or in the other if 
% both had the same size, i.e. If the registration is histology to MRI, the
% the function would estimate the coordinates if the MRI cropped image had 
% a similar size as the histology one. If the registration is MRI to
% histology, the function would estimate the coordinates if the histology
% had a similar size as the MRI cropped ones.
%
% INPUT
%   - pointsMRI: Array Nx2 with the list of N 2D points (row,col) in the MRI 
%   image.
%   - pointsHistology: Array Nx2 with the list of N 2D points (row,col) in the 
%   MRI image. 
%   - imagesMRI: Cell array with the MRI images (cropped).
%   - sizeHistolImageGUI: Size of the histology image shown in GUI.
%   - imageHistology: Histology image to be registered. It is usually resized
%   to prevent memory issues.
%   - hist2MRI: Bool variable to indicate if the registration will be MRI
%   to histology (true) or histology to MRI (false)
%       
% OUTPUT
%   - imageRegBSpline: Registered image with B-Splines. If the registration
%   is MRI to histology, it is a cell array containing the different 
%   (registered) modalities.
%
% Created by:	Victor Gonzalez Castro
% Funded by:    Row Fogo Charitable Trust
%

    %sizeMRIImage = size(imagesMRI{1}); % There will be at least one MRI image
    sizeHistolImage = size(imageHistology);
    
    % GET THE RATIO BETWEEN IMAGES (i.e. How should we multiply the MRI
    % size to get (approximately) the original histology image)
%     ratX = sizeHistolImage(2)/sizeMRIImage(2);
%     ratY = sizeHistolImage(1)/sizeMRIImage(1);
%     ratio = mean([ratX, ratY]);
    % Ratio between the distances of the 2 more distant points in histology and MRI
    % It should be more robust than taking into account the rows and
    % columns of the image (the image may not have the same proportion as
    % the actual tissues)
    dists = squareform(pdist(pointsHistology));
    [points,~] = find(dists==max(max(dists)));
    distPointsHistol = dists(points(1), points(2));
    distPointsMRI = pdist2(pointsMRI(points(1),:), pointsMRI(points(2),:));
    ratio = distPointsHistol/distPointsMRI;
    
    % Ratio to resize GUI Histology to original histology
    %ratioHistol = sizeHistolImage(1)/sizeHistolImageGUI(1);
    %ratioHistol = mean([sizeHistolImage(1)/sizeHistolImageGUI(1), sizeHistolImage(2)/sizeHistolImageGUI(2)]);
    ratioHistolX = sizeHistolImage(2)/sizeHistolImageGUI(2);
    ratioHistolY = sizeHistolImage(1)/sizeHistolImageGUI(1);

    % IF THE REGISTRATION IS HISTOLOGY -> MRI, RESIZE MRI POINTS
    if(hist2MRI)
        pointsHistology(:,1) = pointsHistology(:,1)*ratioHistolY;
        pointsHistology(:,2) = pointsHistology(:,2)*ratioHistolX;
        pointsMRI = pointsMRI*ratio;
        %pointsMRI(:,1) = pointsMRI(:,1)*ratY; % rows == y
        %pointsMRI(:,2) = pointsMRI(:,2)*ratX; % cols == x
        
        %imageHistology = imread(pathHistolImage);
        options.Verbose = true;
        options.MaxRef = 10;
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%% FIRST: Affine registration %%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Needs (x,y) coordinates
%         tFormAff = fitgeotrans([pointsHistology(:,2), pointsHistology(:,1)], [pointsMRI(:,2), pointsMRI(:,1)], 'affine'); 
%         imageRegAffine = imwarp(imageHistology, tFormAff, 'nearest');
% 
%         % Apply affine registration to points
%         % PROVISIONAL: WOULD BE BETTER TO APPLY THE TRANSFORM TO THE POINTS
%         % ANALITICALLY
%         %[xHistAff, yHistAff] = tformfwd(tFormAff, pointsHistology(:,2),pointsHistology(:,1));
%         %[xHistAff, yHistAff,zHistAff] = transformPointsForward(tFormAff, pointsHistology(:,2), pointsHistology(:,1), zeros(size(pointsHistology,1),1));
% %         [xHistAff, yHistAff] = transformPointsForward(tFormAff, pointsHistology(:,2), pointsHistology(:,1));
% %         pointsHistology2 = [yHistAff, xHistAff];
%         imPointsHistol = im2uint8(zeros(size(imageHistology,1), size(imageHistology,2)));
%         for i=1:size(pointsHistology,1)
%             p = round(pointsHistology(i,:));
%             listPoints = [p(1)-50:p(1)+50; p(2)-50:p(2)+50]';
%             [rDiscard, ~] = find(listPoints(:,1)<=0 | listPoints(:,2)<=0 | listPoints(:,1)>sizeHistolImage(1) | listPoints(:,2)>sizeHistolImage(2));
%             rDiscard = unique(rDiscard);
%             listPoints(rDiscard,:) = [];
%             %imPointsHistol(p(1)-50:p(1)+50, p(2)-50:p(2)+50) = i;
%             %%%imPointsHistol(sub2ind(sizeHistolImage, listPoints(:,1), listPoints(:,2))) = i;
%             imPointsHistol(listPoints(:,1), listPoints(:,2)) = i;
%             %imPointsHistol(p(1), p(2)) = i;
%         end
% %         figure; imshow(imPointsHistol, []);
%         imPointsHistol_affine = imwarp(imPointsHistol, tFormAff, 'nearest');
% %         figure; imshow(imPointsHistol_affine,[]);
%         pointsHistologyAff = zeros(size(pointsHistology));
%         centroidsAff = regionprops(logical(imPointsHistol_affine), 'Centroid');
%         for i=1:size(pointsHistology,1)
%             cent = [centroidsAff(i).Centroid(2), centroidsAff(i).Centroid(1)];
%             index = imPointsHistol_affine(round(cent(1)), round(cent(2)));
%             pointsHistologyAff(index,:) = cent;
%         end
%         
% %         imHist_affine2 = im2uint8(imHist_affine);
% %         for i=1:size(pointsHistologyAff,1)
% %             p = round(pointsHistologyAff(i,:));
% %             imHist_affine2(p(1)-60:p(1)+60, p(2)-60:p(2)+60,1) = 0;
% %             imHist_affine2(p(1)-60:p(1)+60, p(2)-60:p(2)+60,2) = 255;
% %             imHist_affine2(p(1)-60:p(1)+60, p(2)-60:p(2)+60,3) = 0;
% %         end
% %         figure;imshow(imHist_affine2);
% %         for i=1:size(pointsHistologyAff,1)
% %             p = round(pointsHistologyAff(i,:));
% %             text(p(2),p(1), num2str(i), 'Color', 'k', 'FontSize', 16);
% %         end
% %         clear imageHistology

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% SECOND: B-Spline registration %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%         [O_trans, Spacing] = point_registration(size(imageRegAffine), pointsHistologyAff, pointsMRI, options); % TEST POINTS, IF IT IS (x,y) OR (row,col)
        [O_trans, Spacing] = point_registration(size(imageHistology), pointsHistology, pointsMRI, options); % TEST POINTS, IF IT IS (x,y) OR (row,col)
%         imageRegBSpline = bspline_transform(O_trans, imageRegAffine, Spacing, 3);
        imageRegBSpline = bspline_transform(O_trans, imageHistology, Spacing, 3);
        
    else % IF THE REGISTRATION IS MRI -> HISTOLOGY, RESIZE MRI POINTS
        imageRegBSpline = cell(size(imagesMRI));
        %imageRegAffine = cell(size(imagesMRI));
        %pointsHistology = (pointsHistology*ratioHistol)/ratio; 
        pointsHistology(:,1) = (pointsHistology(:,1)*ratioHistolY)/ratio; % rows == y 
        pointsHistology(:,2) = (pointsHistology(:,2)*ratioHistolX)/ratio; % cols == x
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%% FIRST: Affine registration %%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Needs (x,y) coordinates
%         tFormAff = fitgeotrans([pointsMRI(:,2), pointsMRI(:,1)], [pointsHistology(:,2), pointsHistology(:,1)], 'affine'); 
%         for i=1:length(imagesMRI)
%             imageRegAffine{i} = imwarp(imagesMRI{i}, tFormAff, 'nearest');
%         end
%         
%         % Apply affine registration to points
%         % PROVISIONAL: WOULD BE BETTER TO APPLY THE TRANSFORM TO THE POINTS
%         % ANALITICALLY
%         imPointsMRI = im2uint8(zeros(size(imagesMRI{i},1), size(imagesMRI{i},2)));
%         for i=1:size(pointsMRI,1)
%             p = round(pointsMRI(i,:));
%             listPoints = [p(1)-1:p(1)+1; p(2)-1:p(2)+1]';
%             [rDiscard, ~] = find(listPoints(:,1)<=0 | listPoints(:,2)<=0 | listPoints(:,1)>sizeMRIImage(1) | listPoints(:,2)>sizeMRIImage(2));
%             rDiscard = unique(rDiscard);
%             listPoints(rDiscard,:) = [];
%             %imPointsHistol(p(1)-50:p(1)+50, p(2)-50:p(2)+50) = i;
%             %%%imPointsHistol(sub2ind(sizeHistolImage, listPoints(:,1), listPoints(:,2))) = i;
%             imPointsMRI(listPoints(:,1), listPoints(:,2)) = i;
%             %imPointsHistol(p(1), p(2)) = i;
%         end
% 
%         imPointsMRI_affine = imwarp(imPointsMRI, tFormAff, 'nearest');
%         pointsMRIAff = zeros(size(pointsMRI));
%         centroidsAff = regionprops(logical(imPointsMRI_affine), 'Centroid');
%         for i=1:size(pointsMRI,1)
%             cent = [centroidsAff(i).Centroid(2), centroidsAff(i).Centroid(1)];
%             index = imPointsMRI_affine(round(cent(1)), round(cent(2)));
%             pointsMRIAff(index,:) = cent;
%         end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% SECOND: B-Spline registration %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        options.Verbose = true;
        options.MaxRef = 10;
        [O_trans, Spacing] = point_registration(sizeHistolImage, pointsMRI, pointsHistology, options); % TEST POINTS, IF IT IS (x,y) OR (row,col)
        for i=1:length(imagesMRI)
            imageRegBSpline{i} = bspline_transform(O_trans, imagesMRI{i}, Spacing, 3);
        end
    end

end