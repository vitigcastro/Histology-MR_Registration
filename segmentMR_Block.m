function [imgCropped, maskCrop] = segmentMR_Block(imagMRI, modality, varargin)
% This function crops the MR block img. The algorithm is suitable for the
% FLAIR and T1W image. The T2W and T2star are cropped using the mask
% obtained for the FLAIR image (we are assuming that all modalities are
% co-registered)
%
% INPUT
%   - imagMRI: MR block of the image to be segmented.
%   - modality: String continaing the modality of img. It can be either
%   'FLAIR', 'T1W', 'T2W' or 'T2star'.
%   - mask (OPTIONAL): Mask of the FLAIR image, if available.
%
% OUTPUT
%   - imgCropped: Cropped image.
%   - maskCrop: Mask of the cropped image.

    % CHECK INPUTS 
    
    % Get mask, if available
    if(nargin>2)
        mask = varargin{1};
    end
    
    % Check the modality
    switch lower(modality)
        case 'flair'
            % Nothing to check
        case 't1w'
            % Nothing to check
        case 't2w' % Error if mask does not exist
            if ~exist('mask','var')
                error('You must give a mask as an input image');
            end
        case 't2star' % Error if mask does not exist
            if ~exist('mask','var')
                error('You must give a mask as an input image');
            end
        otherwise
            error('Modality not known');
    end
    
    % SEGMENTATION OF THE FLAIR OR T1
    if(strcmpi(modality, 'FLAIR') || (strcmpi(modality, 'T1W') && ~exist('mask','var')))
            
        %Convert into uint8
        imagMRI8 = uint8((double(imagMRI)./double(intmax(class(imagMRI)))).*double(intmax('uint8')));
        %imagMRI8Adj = adapthisteq(imagMRI8);
        imagMRI8Adj = imadjust(imagMRI8);

        % Thresholding of the tissue sample
        %imBW = im2bw(imag8Adj, graythresh(imag8Adj)); % Otsu
        imBW = im2bw(imagMRI8Adj, 0.5);
        %imBW = im2bw(imagMRI8, 0.5);
        
        % Remove regions smaller than the bigger one
        %labs = bwlabel(imBW);
        %props = regionprops(labs, 'Area', 'PixelIdxList');
        props = regionprops(logical(imBW), 'Area', 'PixelIdxList');
        maxArea = max([props.Area]);
        smallRegions = [props(:).Area] < maxArea;
        imBWMask = imBW; 
        for k=1:length(smallRegions)
            if(smallRegions(k))
                imBWMask([props(k).PixelIdxList]) = 0;
            end
        end

        % Fill the holes
        imBWMaskNoHoles = imfill(imBWMask, 'holes');
        
        se = strel('square',3);
        imBWMaskDil = imdilate(imBWMaskNoHoles, se);

        % Extract bounding box
        %propsBB = regionprops(imBWMaskNoHoles, 'BoundingBox');
        propsBB = regionprops(imBWMaskDil, 'BoundingBox');

%         imBWCrop = imcrop(imBWMaskNoHoles, propsBB(1).BoundingBox);
        %imCrop = imcrop(imagMRI, propsBB(1).BoundingBox);
        imgCropped = imcrop(imagMRI8, propsBB(1).BoundingBox);

        % Create the mask
        [r,c] = find(imBWMaskDil);
        maskCrop = false(size(imagMRI));
        maskCrop(min(r):max(r), min(c):max(c)) = true;
        
    else
        
        imagMRI8 = uint8((double(imagMRI)./double(intmax(class(imagMRI)))).*double(intmax('uint8')));
        
        % Crop image with existing mask
        propsBB = regionprops(mask, 'BoundingBox');
        imgCropped = imcrop(imagMRI8, propsBB(1).BoundingBox);
        
        % The output maskCrop is the mask
        maskCrop = mask;
        
    end

end