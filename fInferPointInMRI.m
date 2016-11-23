function pointMRI = fInferPointInMRI(pHistol, refPointsMRI, refPointsHistol, sizeMRI)
% Given a point in the histology image, this function infers where the 
% equivalent point is in the registered histology image. This is done by 
% estimating distances to a set of reference points in the original image 
% (MRI) and in the target image (Histology). 
% TO TEST: TRY WITH THE TRANSFORMATION MATRIX
%
% INPUT
%   - pHistol: Point selected in the Histology image
%   - refPointsMRI: Reference points in the MRI image (i.e. the points
%   selected in the MR image to make the registration)
%   - refPointsHistol: Reference points in the Histology image (i.e. the 
%   points selected in the histology image to make the registration).
%   - sizeMRI: Size of the MRI image
%
% OUTPUT
%   pointsMRI: Points in the Histology image that corresponds to the 
%   input MRI points.
%   
% Software Version: 1.0
% Created: 07/11/15
% Written by Dr Victor Gonzalez Castro
% Centre for Clinical Brain Sciences (CCBS) - Edinburgh
% email: victor.gonzalez@ed.ac.uk
%

    % Margin to estimate the distance from the reference points
    delta = 2;
    
    % 1. GET THE RATIO OF THE DISTANCES BETWEEN THE MOST DISTANT REFERENCE
    % POINTS IN THAT IMAGE
    dists = squareform(pdist(refPointsMRI));
    [points,~] = find(dists==max(max(dists)));
    distPointsMRI = dists(points(1), points(2));
    distPointsHistol = pdist2(refPointsHistol(points(1),:), refPointsHistol(points(2),:));
    ratio = distPointsMRI/distPointsHistol;
    
    % 2. FIND THE DISTANCES FROM THE POINT TO EACH OF THE REFERENCE POINTS
    % IN THE MR IMAGE
    distsPoint2RefsHistol = pdist2(pHistol, refPointsHistol);
    distsPoint2RefsMRI_hat = distsPoint2RefsHistol * ratio;
    
    % 3. GET THE DISTANCE TRANSFORMS FROM THE REFERENCE POINTS OF THE 
    % MR IMAGE 
    numRefPoints = size(refPointsMRI,1);
    rndRefPointsMRI = round(refPointsMRI);
%     distTransfs = zeros([sizeHistology(1:2), numRefPoints]);
    distsEstimated = zeros([sizeMRI(1:2), numRefPoints]);
    for numIm=1:numRefPoints
        % Compute distance transform of a binary image where the numIm-th
        % histology reference point is the only white pixel
        binImageRefMRIPoint = false(sizeMRI(1:2));
        %binImageRefMRIPoint(rndRefPointsMRI(numIm,1)-1:rndRefPointsMRI(numIm,1)+1,rndRefPointsMRI(numIm,2)-1:rndRefPointsMRI(numIm,2)+1) = 1;
        binImageRefMRIPoint(rndRefPointsMRI(numIm,1), rndRefPointsMRI(numIm,2)) = 1;
%         distTransfs(:,:,numIm) = bwdist(binImageRefHistPoint);
        distTransf = bwdist(binImageRefMRIPoint);
%         figure; imshow(binImageRefHistPoint);
%         figure; imshow(distTransfs(:,:,numIm),[]);
        
        % Image with only the pixels at a distance
        % distsPoint2RefsHistol_hat are white
%         distsEstimated(:,:,numIm) = (distTransfs(:,:,numIm)>=distsPoint2RefsHistol_hat(numIm)-delta & distTransfs(:,:,numIm)<=distsPoint2RefsHistol_hat(numIm)+delta);
        distsEstimated(:,:,numIm) = (distTransf>=distsPoint2RefsMRI_hat(numIm)-delta & distTransf<=distsPoint2RefsMRI_hat(numIm)+delta);
%         figure; imshow(distsEstimated(:,:,numIm));
    end
    
    % Find the region in the image where there are more coincidence of
    % circles
    combinedDists = sum(distsEstimated,3);
    
    % The output point will be taken as the centroid of this region
    [r,c] = find(combinedDists==max(combinedDists(:)));
    pointMRI = mean([r,c],1);
    
    %disp('done');
end