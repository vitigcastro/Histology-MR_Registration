function pointHistol = fInferPointInHistol(pMRI, refPointsMRI, refPointsHistol, sizeHistology)
% Given a point in the MRI image, this function infers where the equivalent
% point is in the registered histology image. This is done by estimating
% distances to a set of reference points in the original image (MRI) and in
% the target image (Histology). 
% TO TEST: TRY WITH THE TRANSFORMATION MATRIX
%
% INPUT
%   - pMRI: Point selected in the MRI image
%   - refPointsMRI: Reference points in the MRI image (i.e. the points
%   selected in the MR image to make the registration)
%   - refPointsHistol: Reference points in the Histology image (i.e. the 
%   points selected in the histology image to make the registration).
%   - sizeHistology: Size of the Histology image.
%
% OUTPUT
%   pointsHistol: Points in the Histology image that corresponds to the 
%   input MRI points.
%   
% Software Version: 1.0
% Created: 07/11/15
% Written by Dr Victor Gonzalez Castro
% Centre for Clinical Brain Sciences (CCBS) - Edinburgh
% email: victor.gonzalez@ed.ac.uk
%

    % Margin to estimate the distance from the reference points
    %delta = 50;
    delta = 100;
    
    % 1. GET THE RATIO OF THE DISTANCES BETWEEN THE MOST DISTANT REFERENCE
    % POINTS IN THAT IMAGE
    dists = squareform(pdist(refPointsHistol));
    [points,~] = find(dists==max(max(dists)));
    distPointsHistol = dists(points(1), points(2));
    distPointsMRI = pdist2(refPointsMRI(points(1),:), refPointsMRI(points(2),:));
    ratio = distPointsHistol/distPointsMRI;
    
    % 2. FIND THE DISTANCES FROM THE POINT TO EACH OF THE REFERENCE POINTS
    % IN THE MR IMAGE
    distsPoint2RefsMRI = pdist2(pMRI, refPointsMRI);
    distsPoint2RefsHistol_hat = distsPoint2RefsMRI * ratio;
    
    % 3. GET THE DISTANCE TRANSFORMS FROM THE REFERENCE POINTS OF THE 
    % HISTOLOGY IMAGES 
    numRefPoints = size(refPointsHistol,1);
    rndRefPointsHistol = round(refPointsHistol);
%     distTransfs = zeros([sizeHistology(1:2), numRefPoints]);
    distsEstimated = zeros([sizeHistology(1:2), numRefPoints]);
    for numIm=1:numRefPoints
        % Compute distance transform of a binary image where the numIm-th
        % histology reference point is the only white pixel
        binImageRefHistPoint = false(sizeHistology(1:2));
        %binImageRefHistPoint(rndRefPointsHistol(numIm,:)) = 1;
        binImageRefHistPoint(rndRefPointsHistol(numIm,1)-10:rndRefPointsHistol(numIm,1)+10,rndRefPointsHistol(numIm,2)-10:rndRefPointsHistol(numIm,2)+10) = 1;
%         distTransfs(:,:,numIm) = bwdist(binImageRefHistPoint);
        distTransf = bwdist(binImageRefHistPoint);
%         figure; imshow(binImageRefHistPoint);
%         figure; imshow(distTransfs(:,:,numIm),[]);
        
        % Image with only the pixels at a distance
        % distsPoint2RefsHistol_hat are white
%         distsEstimated(:,:,numIm) = (distTransfs(:,:,numIm)>=distsPoint2RefsHistol_hat(numIm)-delta & distTransfs(:,:,numIm)<=distsPoint2RefsHistol_hat(numIm)+delta);
        distsEstimated(:,:,numIm) = (distTransf>=distsPoint2RefsHistol_hat(numIm)-delta & distTransf<=distsPoint2RefsHistol_hat(numIm)+delta);
%         figure; imshow(distsEstimated(:,:,numIm));
    end
    
    % Find the region in the image where there are more coincidence of
    % circles
    combinedDists = sum(distsEstimated,3);
    
    % The output point will be taken as the centroid of this region
    [r,c] = find(combinedDists==max(combinedDists(:)));
    pointHistol = mean([r,c],1);
    
    %disp('done');
end