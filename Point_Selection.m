function varargout = Point_Selection(varargin)
% POINT_SELECTION MATLAB code for Point_Selection.fig
%      POINT_SELECTION, by itself, creates a new POINT_SELECTION or raises the existing
%      singleton*.
%
%      H = POINT_SELECTION returns the handle to a new POINT_SELECTION or the handle to
%      the existing singleton*.
%
%      POINT_SELECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POINT_SELECTION.M with the given input arguments.
%
%      POINT_SELECTION('Property','Value',...) creates a new POINT_SELECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Point_Selection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Point_Selection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Point_Selection

% Last Modified by GUIDE v2.5 23-Nov-2016 15:13:36

% QUESTIONS FOR CAT:
%
% 1. When New reference points are loaded... should the reference points be
% deleted from the table of the regions of interest?
%   Yes, they shoul. The reason is that the coordinates may change in
%   different modalities, even if they are registered. TO BE TESTED
%
%
%
%

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @Point_Selection_OpeningFcn, ...
                       'gui_OutputFcn',  @Point_Selection_OutputFcn, ...
                       'gui_LayoutFcn',  [] , ...
                       'gui_Callback',   []);
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
% End initialization code - DO NOT EDIT

end
% *************************************************************************

% --- Executes just before Point_Selection is made visible.
function Point_Selection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Point_Selection (see VARARGIN)

    
    maxMRIImages = 4;
    
    % Choose default command line output for Point_Selection
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes Point_Selection wait for user response (see UIRESUME)
    % uiwait(handles.figure1);

     % Add additional data as a new field inside axesMRI. Specifically:
    % - A cell array to store the MRI images
    data = guidata(handles.axesMRI);
    data.data.MRIImage = cell(3,maxMRIImages); % Cell array for the MR images
    data.data.tableCoordinatesMRI = zeros(1,2);
    data.data.tableCoordinatesHistol = zeros(1,2);
    % Cell array for the selected labelled points. 
    %   · The first cell stores the array with the coordinates of the 
    %   points in the MR images
    %   · The second cell stores the array with the coordinates of the
    %   points in the histology image
    data.data.tableCoordsPointsLab = cell(1,2); 
    
    % Add additional data as a new field inside axesHist. Specifically:
    % - A cell array to store the histology images
    data.data.histologyImage = cell(1,3); % Cell array for the MR images
    
    guidata(handles.axesMRI, data);
    
    % Set uitables to empty
    dataTableMRI = get(handles.uitablePointsMRI, 'Data');
%     dataTableMRI(4,:) = []; 
%     dataTableMRI(3,:) = [];
%     dataTableMRI(2,:) = [];
%     dataTableMRI(1,:) = [];
    numPointsMRI = size(dataTableMRI,1);
    for nP=1:numPointsMRI
        dataTableMRI(1,:) = []; 
    end
    set(handles.uitablePointsMRI, 'Data', dataTableMRI);
    
    dataTableHistol = get(handles.uitablePointsHistol, 'Data');
%     dataTableHistol(4,:) = [];
%     dataTableHistol(3,:) = [];
%     dataTableHistol(2,:) = [];
%     dataTableHistol(1,:) = [];
    numPointsHistol = size(dataTableHistol,1);
    for nP=1:numPointsHistol
        dataTableHistol(1,:) = [];
    end
    set(handles.uitablePointsHistol, 'Data', dataTableHistol);
    
    
    dataTablePointsLab = get(handles.uitablePointsLab, 'Data');
    numPointsLab = size(dataTablePointsLab,1);
    for nP=1:numPointsLab
        dataTablePointsLab(1,:) = []; 
    end
    set(handles.uitablePointsLab, 'Data', dataTablePointsLab);
    
    % Put actual directory path in editWorkDir 
    set(handles.editWorkDir, 'String', pwd);
     
end
% *************************************************************************


% --- Outputs from this function are returned to the command line.
function varargout = Point_Selection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;

end
% *************************************************************************

function editWorkDir_Callback(hObject, eventdata, handles)
% hObject    handle to editWorkDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of editWorkDir as text
    %        str2double(get(hObject,'String')) returns contents of editWorkDir as a double

end
% *************************************************************************

% --- Executes during object creation, after setting all properties.
function editWorkDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editWorkDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

end
% *************************************************************************

% --- Executes on button press in pushbuttonSelDir.
function pushbuttonSelDir_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSelDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    % Get the current working directory
    currentWorkDir = get(handles.editWorkDir, 'String');
    
    % Choose a new directory
    folderName = uigetdir(currentWorkDir, 'Choose working directory');
    
    % Set the new working directory in edit3
    set(handles.editWorkDir, 'String', folderName);
end
% *************************************************************************

% --- Executes on button press in pushbuttonHist.
function pushbuttonHist_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonHist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    valueHistology = 1;
    
    % GET WORKING DIRECTORY
    workDir = get(handles.editWorkDir, 'String');
    
    % Get MRI file name
    [fileImage, pathImage] = uigetfile({'*.jpg;*.jpeg;*.png',...
    'Image Files (*.jpg,*.jpeg,*.png)';
    '*.jpg;*.jpeg','JPEG image (*.jpg;*.jpeg)'; ...
    '*.png',  'Portable Network Graphics (*.png)'; ...
    '*.tif;*.tiff',  'Tagged Image File Format (*.tig;*tiff)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Choose the H&E histology image', ...
    workDir);

    try
        % SHOW HELP DIALOG
        hDlg = helpdlg('Loading Histology image.', 'Loading image');
        
        % Load image
        imageHistology = imread([pathImage, fileImage]);
        
%         % Resize image (TEST: Set scale to 0.25 (i.e. 1/4 of the original
%         % size))
        %imageHistologyResize = imresize(imageHistology, 0.25);
        % Resize image (TEST: Set scale to 0.5 (i.e. 1/2 of the original
        % size))
%         imageHistologyResize = imresize(imageHistology, 0.5);
        sizeHistol = size(imageHistology);
%         clear imageHistology;
        
        % ADD TO THE AXIS. IF IT IS NOT THE FIRST ONE, CONTROL THE SLIDER
        % Store in the cell array inside axis1
        data = guidata(handles.axesMRI);
        data.data.histologyImage{valueHistology, 1} = imageHistology;
%         data.data.histologyImage{valueHistology,1} = imageHistologyResize;
        % Store also path of the original histology image
        data.data.histologyImage{valueHistology, 2} = [pathImage, fileImage];
        data.data.histologyImage{valueHistology, 3} = sizeHistol;
        guidata(handles.axesMRI, data);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% SHOW IN THE AXIS %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Next, set your image display axes as the current axes:
        axes(handles.axesHist);
        
        % Get axesHist tag
        tagHist = get(handles.axesHist, 'Tag');
        
        % Display image in axesMRI, using bicubic interpolation
        hHistology = imshow(data.data.histologyImage{valueHistology},[]);
%         resizePos = get(handles.axesHist, 'Position');
%         scaleRes = resizePos(3)/max(size(data.data.histologyImage{valueHistology})); % Axis is squared
%         %imageRes = imresize(data.data.MRIImage{2,valueFLAIR}, [resizePos(4) resizePos(3)]);
%         imageRes = imresize(data.data.histologyImage{valueHistology}, scaleRes, 'bicubic');
%         hHistology = imshow(imageRes);

        % CHECK IF THERE ARE POINTS ALREADY SELECTED, AND SHOW THEM IF
        % NECESSARY
        selPointsHistol = data.data.tableCoordinatesHistol;
        dataTableHistol = get(handles.uitablePointsHistol, 'Data'); % Gets current data of the table
        numPointsHistol = size(dataTableHistol,1);
        if(numPointsHistol>0)
            axes(handles.axesHist);
            for nP=1:size(selPointsHistol,1)
                coordinates = [selPointsHistol(nP,2), selPointsHistol(nP,1)]; %(x,y)
                % Put a marker on the point
                hold on;
                hPoint = plot(coordinates(1), coordinates(2), 'b*', 'MarkerSize', 9);
%                 hPoint = impoint(handles.axesMRI, coordinates);
                set(hPoint, 'UserData', nP);
% %                 % Construct boundary constraint function
% %                 fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
% %                 % Enforce boundary constraint function using setPositionConstraintFcn
% %                 setPositionConstraintFcn(hPoint, fcn);
%                 setColor(hPoint, 'r');
% %                 addNewPositionCallback(hPoint, @pointDragMRI);
%                 %setString(hPoint, num2str(k));
                %hTxt = 
                text(coordinates(1), coordinates(2), sprintf('%d',nP), ...
                    'Color', 'r', 'FontSize',8, ...
                    'HorizontalAlignment','left', 'VerticalAlignment','top');
            end
        end
        
        selPointsHistolLab = data.data.tableCoordsPointsLab{2};
        dataTableHistolLab = get(handles.uitablePointsLab, 'Data');
        numPointsHistolLab = size(dataTableHistolLab,1);
        if(numPointsHistolLab>0)
            for nP=1:size(selPointsHistolLab, 1)
                coordinates = [selPointsHistolLab(nP,2),selPointsHistolLab(nP,1)];
                hold on;
                % Put a marker on the point
                hPoint = impoint(handles.axesHist, coordinates);
                set(hPoint, 'UserData', nP);
                % Construct boundary constraint function
                fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
                % Enforce boundary constraint function using setPositionConstraintFcn
                setPositionConstraintFcn(hPoint, fcn);
                setColor(hPoint, 'g');
                % TO-DO: DECIDE WHETHER WE WANT TO LET THE POINTS TO BE DRAGGED OR
                % NOT
                %addNewPositionCallback(hPoint, @pointDragMRI);
                %hTxt = 
                text(coordinates(1), coordinates(2), sprintf('%d',k), ...
                    'Color', 'g', 'FontSize',8, ...
                    'HorizontalAlignment','left', 'VerticalAlignment','top');
            end
        end

        % Set again the properties lost in axesHist
        set(handles.axesHist, 'Tag', tagHist);
        %set(handles.axesMRI, 'ButtonDownFcn', @axesMRI_ButtonDownFcn);
        set(hHistology, 'ButtonDownFcn', {@getClicksHistology, handles});
        
        % ACTIVATE THE CORRESPONDING RADIOBUTTON
        set(handles.rbHE, 'Value', 1);
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%% DELETE ALL POINTS FROM TABLES AND AXES, IF ANY %%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         % 1 - DELETE FROM TABLES IN AXES
%         % Get the table with the points from handles.axis1
%         dataHandles = guidata(handles.axesMRI);
%     
%         % Copy an empty matrix into dataHandles.data.tableCoordinatesMRI
%         % and dataHandles.data.tableCoordinatesHistol
%         dataHandles.data.tableCoordinatesMRI = [];
%         dataHandles.data.tableCoordinatesHistol = [];
%     
%         % Save the points into handles.axisMRI
%         guidata(handles.axesMRI, dataHandles);
%         
%         % 2 - REMOVE FROM GUI TABLES
%         set(handles.uitablePointsMRI, 'Data', []);
%         set(handles.uitablePointsHistol, 'Data', []);
%     
%         % 3 - DELETE FROM IMAGES AXES
%         % If there are po THERE ARE POINTS IN THE IMAGE, REMOVE THEM AND PAINT THEM AGAIN
%         % axesMRI
%         axes(handles.axesMRI);
%         axesMRI_chil = get(handles.axesMRI, 'Children');
%         class_axesMRI = arrayfun(@class, axesMRI_chil, 'UniformOutput', false);
%         isImage = strcmpi('matlab.graphics.primitive.Image', class_axesMRI);
%         delete(axesMRI_chil(~isImage));
% 
%         % axesHist
%         axes(handles.axesHist);
%         axesHist_chil = get(handles.axesHist, 'Children');
%         class_axesHist = arrayfun(@class, axesHist_chil, 'UniformOutput', false);
%         isImage = strcmpi('matlab.graphics.primitive.Image', class_axesHist);
%         delete(axesHist_chil(~isImage));
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        close(hDlg);
    catch
        errordlg('Choose a JPEG image', 'Error reading the image');
    end

end
% *************************************************************************

% --- Executes on button press in pushbuttonFLAIR.
function pushbuttonFLAIR_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonFLAIR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    valueFLAIR = 1;

    % GET WORKING DIRECTORY
    workDir = get(handles.editWorkDir, 'String');
    
    % Get MRI file name
    [fileImage, pathImage] = uigetfile({'*.dcm;*.jpg;*.jpeg',...
    'Image Files (*.dcm,*.jpg,*.jpeg)';
    '*.dcm',  'DICOM (*.dcm)'; ...
    '*.jpg;*.jpeg','JPEG image (*.jpg;*.jpeg)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Choose the FLAIR image', ...
    workDir);
        
    % READ IMAGE FILE 
    try
        [~,~,ext] = fileparts(fileImage);
        
        if(strcmp(ext,'.dcm'))
            image = dicomread([pathImage, fileImage]);
        else
            image = imread([pathImage, fileImage]);
        end
        
        % SEGMENT THE FLAIR IMAGE (in a external function)
        [imgCropped, maskCrop] = segmentMR_Block(image, 'FLAIR');
                
        % ADD TO THE AXIS. IF IT IS NOT THE FIRST ONE, CONTROL THE SLIDER
        % Store in the cell array inside axis1
        data = guidata(handles.axesMRI);
        %data.MRIImage{valueFLAIR} = image;
        data.data.MRIImage{1,valueFLAIR} = image;
        data.data.MRIImage{2,valueFLAIR} = imgCropped;
        data.data.MRIImage{3,valueFLAIR} = maskCrop;
        guidata(handles.axesMRI, data);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% SHOW IN THE AXIS %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Next, set your image display axes as the current axes:
        axes(handles.axesMRI);
        
        % Get axesMRI tag
        tag = get(handles.axesMRI, 'Tag');
        
%         % Display image in axesMRI using bicubic interpolation
        hFLAIR = imshow(data.data.MRIImage{2,valueFLAIR},[]);
%         resizePos = get(handles.axesMRI, 'Position');
%         scaleRes = resizePos(3)/max(size(data.data.MRIImage{2,valueFLAIR})); % Axis is squared
%         %imageRes = imresize(data.data.MRIImage{2,valueFLAIR}, [resizePos(4) resizePos(3)]);
%         imageRes = imresize(data.data.MRIImage{2,valueFLAIR}, scaleRes, 'bicubic');
%         hFLAIR = imshow(imageRes);

        % CHECK IF THERE ARE POINTS ALREADY SELECTED, AND SHOW THEM IF
        % NECESSARY
        selPointsMRI = data.data.tableCoordinatesMRI;
        dataTableMRI = get(handles.uitablePointsMRI, 'Data'); % Gets current data of the table
        numPointsMRI = size(dataTableMRI,1);
        if(numPointsMRI>0)
            axes(handles.axesMRI);
            for nP=1:size(selPointsMRI,1)
                coordinates = [selPointsMRI(nP,2), selPointsMRI(nP,1)]; %(x,y)
                % Put a marker on the point
                hold on;
                hPoint = plot(coordinates(1), coordinates(2), 'r*', 'MarkerSize', 9);
%                 hPoint = impoint(handles.axesMRI, coordinates);
                set(hPoint, 'UserData', nP);
% %                 % Construct boundary constraint function
% %                 fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
% %                 % Enforce boundary constraint function using setPositionConstraintFcn
% %                 setPositionConstraintFcn(hPoint, fcn);
%                 setColor(hPoint, 'r');
% %                 addNewPositionCallback(hPoint, @pointDragMRI);
%                 %setString(hPoint, num2str(k));
                %hTxt = 
                text(coordinates(1), coordinates(2), sprintf('%d',nP), ...
                    'Color', 'r', 'FontSize',8, ...
                    'HorizontalAlignment','left', 'VerticalAlignment','top');
            end
        end
        
        selPointsMRILab = data.data.tableCoordsPointsLab{1};
        dataTableMRILab = get(handles.uitablePointsLab, 'Data');
        numPointsMRILab = size(dataTableMRILab,1);
        if(numPointsMRILab>0)
            for nP=1:size(selPointsMRILab, 1)
                coordinates = [selPointsMRILab(nP,2),selPointsMRILab(nP,1)];
                hold on;
                % Put a marker on the point
                hPoint = impoint(handles.axesMRI, coordinates);
                set(hPoint, 'UserData', nP);
                % Construct boundary constraint function
                fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
                % Enforce boundary constraint function using setPositionConstraintFcn
                setPositionConstraintFcn(hPoint, fcn);
                setColor(hPoint, 'g');
                % TO-DO: DECIDE WHETHER WE WANT TO LET THE POINTS TO BE DRAGGED OR
                % NOT
                %addNewPositionCallback(hPoint, @pointDragMRI);
                %hTxt = 
                text(coordinates(1), coordinates(2), sprintf('%d',k), ...
                    'Color', 'g', 'FontSize',8, ...
                    'HorizontalAlignment','left', 'VerticalAlignment','top');
            end
        end

        % Set again the properties lost in axesMRI
        set(handles.axesMRI, 'Tag', tag);
        %set(handles.axesMRI, 'ButtonDownFcn', @axesMRI_ButtonDownFcn);
        set(hFLAIR, 'ButtonDownFcn', {@getClicksMRI, handles});
        
        % ACTIVATE THE CORRESPONDING RADIOBUTTON
        set(handles.rbFLAIR, 'Value', 1);
        
    catch
        if(~isempty(fileImage))
            errordlg('Choose a DICOM or a JPEG image', 'Error reading the image');
        end
    end

end
% *************************************************************************

% --- Executes on button press in pushbuttonT1.
function pushbuttonT1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonT1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    valueT1w = 2;
    
    % GET WORKING DIRECTORY
    workDir = get(handles.editWorkDir, 'String');
    
    % Get MRI file name
    [fileImage, pathImage] = uigetfile({'*.dcm;*.jpg;*.jpeg',...
    'Image Files (*.dcm,*.jpg,*.jpeg)';
    '*.dcm',  'DICOM (*.dcm)'; ...
    '*.jpg;*.jpeg','JPEG image (*.jpg;*.jpeg)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Choose the T1w image', ...
    workDir);
        
    % READ IMAGE FILE 
    try
        [~,~,ext] = fileparts(fileImage);
        
        if(strcmp(ext,'.dcm'))
            image = dicomread([pathImage, fileImage]);
        else
            image = imread([pathImage, fileImage]);
        end
        
        % ADD TO THE AXIS. IF IT IS NOT THE FIRST ONE, CONTROL THE SLIDER
        % Store in the cell array inside axis1
        data = guidata(handles.axesMRI);
        
        % Get the mask from the cell array inside axis1
        mask = data.data.MRIImage{3,1};
        
        % Segment the T1W image using the mask (in a external function)
        [imgCropped, maskCrop] = segmentMR_Block(image, 'T1W', mask);
        
        % Save images
        data.data.MRIImage{1, valueT1w} = image;
        data.data.MRIImage{2, valueT1w} = imgCropped;
        data.data.MRIImage{3, valueT1w} = maskCrop;
        guidata(handles.axesMRI, data);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% SHOW IN THE AXIS %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%        
        % Next, set your image display axes as the current axes:
        axes(handles.axesMRI);
        
        % Get axesMRI tag
        tag = get(handles.axesMRI, 'Tag');
        
        % Display image in axesMRI using bicubic interpolation
        hT1w = imshow(data.data.MRIImage{2,valueT1w},[]);
%         resizePos = get(handles.axesMRI, 'Position');
%         scaleRes = resizePos(3)/max(size(data.data.MRIImage{2,valueT1w})); % Axis is squared
%         %imageRes = imresize(data.data.MRIImage{2,valueFLAIR}, [resizePos(4) resizePos(3)]);
%         imageRes = imresize(data.data.MRIImage{2,valueT1w}, scaleRes, 'bicubic');
%         hT1w = imshow(imageRes);

        % CHECK IF THERE ARE POINTS ALREADY SELECTED, AND SHOW THEM IF
        % NECESSARY
        selPointsMRI = data.data.tableCoordinatesMRI;
        dataTableMRI = get(handles.uitablePointsMRI, 'Data'); % Gets current data of the table
        numPointsMRI = size(dataTableMRI,1);
        if(numPointsMRI>0)
            axes(handles.axesMRI);
            for nP=1:size(selPointsMRI,1)
                coordinates = [selPointsMRI(nP,2), selPointsMRI(nP,1)]; %(x,y)
                % Put a marker on the point
                hold on;
                hPoint = plot(coordinates(1), coordinates(2), 'r*', 'MarkerSize', 9);
%                 hPoint = impoint(handles.axesMRI, coordinates);
                set(hPoint, 'UserData', nP);
% %                 % Construct boundary constraint function
% %                 fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
% %                 % Enforce boundary constraint function using setPositionConstraintFcn
% %                 setPositionConstraintFcn(hPoint, fcn);
%                 setColor(hPoint, 'r');
% %                 addNewPositionCallback(hPoint, @pointDragMRI);
%                 %setString(hPoint, num2str(k));
                %hTxt = 
                text(coordinates(1), coordinates(2), sprintf('%d',nP), ...
                    'Color', 'r', 'FontSize',8, ...
                    'HorizontalAlignment','left', 'VerticalAlignment','top');
            end
        end
        
        selPointsMRILab = data.data.tableCoordsPointsLab{1};
        dataTableMRILab = get(handles.uitablePointsLab, 'Data');
        numPointsMRILab = size(dataTableMRILab,1);
        if(numPointsMRILab>0)
            for nP=1:size(selPointsMRILab, 1)
                coordinates = [selPointsMRILab(nP,2),selPointsMRILab(nP,1)];
                hold on;
                % Put a marker on the point
                hPoint = impoint(handles.axesMRI, coordinates);
                set(hPoint, 'UserData', nP);
                % Construct boundary constraint function
                fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
                % Enforce boundary constraint function using setPositionConstraintFcn
                setPositionConstraintFcn(hPoint, fcn);
                setColor(hPoint, 'g');
                % TO-DO: DECIDE WHETHER WE WANT TO LET THE POINTS TO BE DRAGGED OR
                % NOT
                %addNewPositionCallback(hPoint, @pointDragMRI);
                %hTxt = 
                text(coordinates(1), coordinates(2), sprintf('%d',k), ...
                    'Color', 'g', 'FontSize',8, ...
                    'HorizontalAlignment','left', 'VerticalAlignment','top');
            end
        end
        
        % Set again the properties lost in axesMRI
        set(handles.axesMRI, 'Tag', tag);
        %set(handles.axesMRI, 'ButtonDownFcn', @axesMRI_ButtonDownFcn);
        set(hT1w, 'ButtonDownFcn', {@getClicksMRI, handles});
        
        % ACTIVATE THE CORRESPONDING RADIOBUTTON
        set(handles.rbT1w, 'Value', 1);
        
        % IF THE IMAGE IS WELL READ, MAKE THIS BUTTON INVISIBLE
        %set(handles.buttonT1w, 'Visible', 'off');
        
    catch
        if(~isempty(fileImage))
            errordlg('Choose a DICOM or a JPEG image', 'Error reading the image');
        end
    end

end
% *************************************************************************

% --- Executes on button press in pushbuttonT2.
function pushbuttonT2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonT2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    valueT2w = 3;
    
    % GET WORKING DIRECTORY
    workDir = get(handles.editWorkDir, 'String');
    
    % Get MRI file name
    [fileImage, pathImage] = uigetfile({'*.dcm;*.jpg;*.jpeg',...
    'Image Files (*.dcm,*.jpg,*.jpeg)';
    '*.dcm',  'DICOM (*.dcm)'; ...
    '*.jpg;*.jpeg','JPEG image (*.jpg;*.jpeg)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Choose the T2w image', ...
    workDir);
        
    % READ IMAGE FILE 
    try
        [~,~,ext] = fileparts(fileImage);
        
        if(strcmp(ext,'.dcm'))
            image = dicomread([pathImage, fileImage]);
        else
            image = imread([pathImage, fileImage]);
        end
        
        % ADD TO THE AXIS. IF IT IS NOT THE FIRST ONE, CONTROL THE SLIDER
        % Store in the cell array inside axis1
        data = guidata(handles.axesMRI);
        
        % Get the mask from the cell array inside axis1
        mask = data.data.MRIImage{3,1};
        
        % Segment the T1W image using the mask (in a external function)
        [imgCropped, maskCrop] = segmentMR_Block(image, 'T2W', mask);
        
        % Save images
        data.data.MRIImage{1, valueT2w} = image;
        data.data.MRIImage{2, valueT2w} = imgCropped;
        data.data.MRIImage{3, valueT2w} = maskCrop;
        guidata(handles.axesMRI, data);
        
     
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% SHOW IN THE AXIS %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%  
        % Next, set your image display axes as the current axes:
        axes(handles.axesMRI);
        
        % Get axesMRI tag
        tag = get(handles.axesMRI, 'Tag');
        
        % Display image in axesMRI using bicubic interpolation
        hT2w = imshow(data.data.MRIImage{2,valueT2w},[]);
%         resizePos = get(handles.axesMRI, 'Position');
%         scaleRes = resizePos(3)/max(size(data.data.MRIImage{2,valueT2w})); % Axis is squared
%         %imageRes = imresize(data.data.MRIImage{2,valueFLAIR}, [resizePos(4) resizePos(3)]);
%         imageRes = imresize(data.data.MRIImage{2,valueT2w}, scaleRes, 'bicubic');
%         hT2w = imshow(imageRes);
        
        % CHECK IF THERE ARE POINTS ALREADY SELECTED, AND SHOW THEM IF
        % NECESSARY
        selPointsMRI = data.data.tableCoordinatesMRI;
        dataTableMRI = get(handles.uitablePointsMRI, 'Data'); % Gets current data of the table
        numPointsMRI = size(dataTableMRI,1);
        if(numPointsMRI>0)
            axes(handles.axesMRI);
            for nP=1:size(selPointsMRI,1)
                coordinates = [selPointsMRI(nP,2), selPointsMRI(nP,1)]; %(x,y)
                % Put a marker on the point
                hold on;
                hPoint = plot(coordinates(1), coordinates(2), 'r*', 'MarkerSize', 9);
%                 hPoint = impoint(handles.axesMRI, coordinates);
                set(hPoint, 'UserData', nP);
%                 % Construct boundary constraint function
% %                 fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
% %                 % Enforce boundary constraint function using setPositionConstraintFcn
% %                 setPositionConstraintFcn(hPoint, fcn);
%                 setColor(hPoint, 'r');
% %                 addNewPositionCallback(hPoint, @pointDragMRI);
%                 %setString(hPoint, num2str(k));
                %hTxt = 
                text(coordinates(1), coordinates(2), sprintf('%d',nP), ...
                    'Color', 'r', 'FontSize',8, ...
                    'HorizontalAlignment','left', 'VerticalAlignment','top');
            end
        end
        
        selPointsMRILab = data.data.tableCoordsPointsLab{1};
        dataTableMRILab = get(handles.uitablePointsLab, 'Data');
        numPointsMRILab = size(dataTableMRILab,1);
        if(numPointsMRILab>0)
            for nP=1:size(selPointsMRILab, 1)
                coordinates = [selPointsMRILab(nP,2),selPointsMRILab(nP,1)];
                hold on;
                % Put a marker on the point
                hPoint = impoint(handles.axesMRI, coordinates);
                set(hPoint, 'UserData', nP);
                % Construct boundary constraint function
                fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
                % Enforce boundary constraint function using setPositionConstraintFcn
                setPositionConstraintFcn(hPoint, fcn);
                setColor(hPoint, 'g');
                % TO-DO: DECIDE WHETHER WE WANT TO LET THE POINTS TO BE DRAGGED OR
                % NOT
                %addNewPositionCallback(hPoint, @pointDragMRI);
                %hTxt = 
                text(coordinates(1), coordinates(2), sprintf('%d',k), ...
                    'Color', 'g', 'FontSize',8, ...
                    'HorizontalAlignment','left', 'VerticalAlignment','top');
            end
        end
        
        % Set again the properties lost in axesMRI
        set(handles.axesMRI, 'Tag', tag);
        %set(handles.axesMRI, 'ButtonDownFcn', @axesMRI_ButtonDownFcn);
        set(hT2w, 'ButtonDownFcn', {@getClicksMRI, handles});
        
        % ACTIVATE THE CORRESPONDING RADIOBUTTON
        set(handles.rbT2w, 'Value', 1);
        
        % IF THE IMAGE IS WELL READ, MAKE THIS BUTTON INVISIBLE
        %set(handles.buttonT2w, 'Visible', 'off');
        
    catch
        if(~isempty(fileImage))
            errordlg('Choose a DICOM or a JPEG image', 'Error reading the image');
        end
    end

end
% *************************************************************************

% --- Executes on button press in pushbuttonT2star.
function pushbuttonT2star_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonT2star (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    valueT2star = 4;
    
    % GET WORKING DIRECTORY
    workDir = get(handles.editWorkDir, 'String');
    
    % Get MRI file name
    [fileImage, pathImage] = uigetfile({'*.dcm;*.jpg;*.jpeg',...
    'Image Files (*.dcm,*.jpg,*.jpeg)';
    '*.dcm',  'DICOM (*.dcm)'; ...
    '*.jpg;*.jpeg','JPEG image (*.jpg;*.jpeg)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Choose the T2star image', ...
    workDir);
        
    % READ IMAGE FILE 
    try
        [~,~,ext] = fileparts(fileImage);
        
        if(strcmp(ext,'.dcm'))
            image = dicomread([pathImage, fileImage]);
        else
            image = imread([pathImage, fileImage]);
        end
        
        % ADD TO THE AXIS. IF IT IS NOT THE FIRST ONE, CONTROL THE SLIDER
        % Store in the cell array inside axis1
        data = guidata(handles.axesMRI);
        
        % Get the mask from the cell array inside axis1
        mask = data.data.MRIImage{3,1};
        
        % Segment the T1W image using the mask (in a external function)
        [imgCropped, maskCrop] = segmentMR_Block(image, 'T2star', mask);
        
        % Save images
        data.data.MRIImage{1, valueT2star} = image;
        data.data.MRIImage{2, valueT2star} = imgCropped;
        data.data.MRIImage{3, valueT2star} = maskCrop;
        guidata(handles.axesMRI, data);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% SHOW IN THE AXIS %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Next, set your image display axes as the current axes:
        axes(handles.axesMRI);
        
        % Get axesMRI tag
        tag = get(handles.axesMRI, 'Tag');
        
        % Display image in axesMRI
        hT2star = imshow(data.data.MRIImage{2,valueT2star},[]);
%         resizePos = get(handles.axesMRI, 'Position');
%         scaleRes = resizePos(3)/max(size(data.data.MRIImage{2,valueT2star})); % Axis is squared
%         %imageRes = imresize(data.data.MRIImage{2,valueFLAIR}, [resizePos(4) resizePos(3)]);
%         imageRes = imresize(data.data.MRIImage{2,valueT2star}, scaleRes, 'bicubic');
%         hT2star = imshow(imageRes);
        
        % CHECK IF THERE ARE POINTS ALREADY SELECTED, AND SHOW THEM IF
        % NECESSARY
        selPointsMRI = data.data.tableCoordinatesMRI;
        dataTableMRI = get(handles.uitablePointsMRI, 'Data'); % Gets current data of the table
        numPointsMRI = size(dataTableMRI,1);
        if(numPointsMRI>0)
            axes(handles.axesMRI);
            for nP=1:size(selPointsMRI,1)
                coordinates = [selPointsMRI(nP,2), selPointsMRI(nP,1)]; %(x,y)
                % Put a marker on the point
                hold on;
                hPoint = plot(coordinates(1), coordinates(2), 'r*', 'MarkerSize', 9);
%                 hPoint = impoint(handles.axesMRI, coordinates);
                set(hPoint, 'UserData', nP);
% %                 % Construct boundary constraint function
% %                 fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
% %                 % Enforce boundary constraint function using setPositionConstraintFcn
% %                 setPositionConstraintFcn(hPoint, fcn);
%                 setColor(hPoint, 'r');
% %                 addNewPositionCallback(hPoint, @pointDragMRI);
%                 %setString(hPoint, num2str(k));
                %hTxt = 
                text(coordinates(1), coordinates(2), sprintf('%d',nP), ...
                    'Color', 'r', 'FontSize',8, ...
                    'HorizontalAlignment','left', 'VerticalAlignment','top');
            end
        end
        
        selPointsMRILab = data.data.tableCoordsPointsLab{1};
        dataTableMRILab = get(handles.uitablePointsLab, 'Data');
        numPointsMRILab = size(dataTableMRILab,1);
        if(numPointsMRILab>0)
            for nP=1:size(selPointsMRILab, 1)
                coordinates = [selPointsMRILab(nP,2),selPointsMRILab(nP,1)];
                hold on;
                % Put a marker on the point
                hPoint = impoint(handles.axesMRI, coordinates);
                set(hPoint, 'UserData', nP);
                % Construct boundary constraint function
                fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
                % Enforce boundary constraint function using setPositionConstraintFcn
                setPositionConstraintFcn(hPoint, fcn);
                setColor(hPoint, 'g');
                % TO-DO: DECIDE WHETHER WE WANT TO LET THE POINTS TO BE DRAGGED OR
                % NOT
                %addNewPositionCallback(hPoint, @pointDragMRI);
                %hTxt = 
                text(coordinates(1), coordinates(2), sprintf('%d',k), ...
                    'Color', 'g', 'FontSize',8, ...
                    'HorizontalAlignment','left', 'VerticalAlignment','top');
            end
        end
        
        % Set again the properties lost in axesMRI
        set(handles.axesMRI, 'Tag', tag);
        %set(handles.axesMRI, 'ButtonDownFcn', @axesMRI_ButtonDownFcn);
        set(hT2star, 'ButtonDownFcn', {@getClicksMRI, handles});
        
        % ACTIVATE THE CORRESPONDING RADIOBUTTON
        set(handles.rbT2star, 'Value', 1);
        
        % IF THE IMAGE IS WELL READ, MAKE THIS BUTTON INVISIBLE
        %set(handles.buttonT2star, 'Visible', 'off');
        
    catch
        if(~isempty(fileImage))
            errordlg('Choose a DICOM or a JPEG image', 'Error reading the image');
        end
    end
end
% *************************************************************************

% --- Executes when selected object is changed in uibuttongroupMRIMods.
function uibuttongroupMRIMods_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroupMRIMods 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    %tagNewVal = get(eventdata.NewValue, 'Tag'); % Get Tag of selected object.
    switch get(eventdata.NewValue, 'Tag') % Get Tag of selected object.
        case 'rbFLAIR'
            newValMRI = 1;
            newModality = 'FLAIR';
        case 'rbT1w'
            newValMRI = 2;
            newModality = 'T1w';
        case 'rbT2w'
            newValMRI = 3;
            newModality = 'T2w';
        case 'rbT2star'
            newValMRI = 4;
            newModality = 'T2star';
    end
    
    % Loaded MRI images
    dataHandles = guidata(handles.axesMRI);
    
    % Look if the modality marked by newVal has already been loaded
    %if(exist('data.MRIImage','var'))
        indexesMRI = find(~cellfun(@isempty,dataHandles.data.MRIImage(:)));
    %else
    %    indexesMRI = 0;
    %end
    
    if(sum(indexesMRI==newValMRI)~=0) % Has been loaded
        %Show in the axis
        % Next, set your image display axes as the current axes:
        axes(handles.axesMRI);
        
        % Get axes1 tag
        tag = get(handles.axesMRI, 'Tag');
        
        % Finally, display the image using bicubic interpolation
        hMRI = imshow(dataHandles.data.MRIImage{2,newValMRI},[]);
%         resizePos = get(handles.axes1, 'Position');
%         scaleRes = resizePos(3)/max(size(dataHandles.data.MRIImage{2,newValMRI})); % Axis is squared
%         %imageRes = imresize(data.data.MRIImage{2,valueFLAIR}, [resizePos(4) resizePos(3)]);
%         imageRes = imresize(dataHandles.data.MRIImage{2,newValMRI}, scaleRes, 'bicubic');
%         hMRI = imshow(imageRes);
        
        % CHECK IF THERE ARE REFERENCE POINTS ALREADY LOADED, AND SHOW 
        % THEM IF NECESSARY
        refPointsMRI = dataHandles.data.tableCoordinatesMRI;
        dataRefPointsTableMRI = get(handles.uitablePointsMRI, 'Data'); % Gets current data of the table
        numRefPointsMRI = size(dataRefPointsTableMRI,1);
        if(numRefPointsMRI>0)
            axes(handles.axesMRI);
            for nP=1:size(refPointsMRI,1)
                coordinates = [refPointsMRI(nP,2), refPointsMRI(nP,1)]; %(x,y)
                % Put a marker on the point
                hPoint = plot(coordinates(1), coordinates(2), 'r*', 'MarkerSize', 9);
%                 hPoint = impoint(handles.axesMRI, coordinates);
                set(hPoint, 'UserData', nP);
% %                 % Construct boundary constraint function
% %                 fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
% %                 % Enforce boundary constraint function using setPositionConstraintFcn
% %                 setPositionConstraintFcn(hPoint, fcn);
%                 setColor(hPoint, 'r');
% %                 addNewPositionCallback(hPoint, @pointDragMRI);
%                 %setString(hPoint, num2str(k));
                %hTxt = 
                text(coordinates(1), coordinates(2), sprintf('%d',nP), ...
                    'Color', 'r', 'FontSize',8, ...
                    'HorizontalAlignment','left', 'VerticalAlignment','top');
            end
        end
        
        % CHECK IF THERE ARE SELECTED POINTS AND SHOW THEM IF NECESSARY
        %dataSelPointsTableMRI = get(handles.uitablePointsLab, 'Data'); % Gets current data of the table
        %selPointsMRI = handles.uitablePointsLab;
        selPointsMRI = get(handles.uitablePointsLab, 'Data'); % Gets current data of the table
        numSelPointsMRI = size(selPointsMRI,1);
        if(numSelPointsMRI>0)
            axes(handles.axesMRI);
            for nP=1:numSelPointsMRI
                coordinates = [selPointsMRI{nP,3}, selPointsMRI{nP,2}]; %(x,y)
                
                hPointMRI = impoint(handles.axesMRI, coordinates);
                set(hPointMRI, 'UserData', nP);
                % Construct boundary constraint function
                fcn = makeConstrainToRectFcn('impoint', get(handles.axesMRI,'XLim'), get(handles.axesMRI,'YLim'));
                % Enforce boundary constraint function using setPositionConstraintFcn
                setPositionConstraintFcn(hPointMRI, fcn);
                setColor(hPointMRI, 'g');
                % TO-DO: DECIDE WHETHER WE WANT TO LET THE POINTS TO BE DRAGGED OR
                % NOT
                %addNewPositionCallback(hPoint, @pointDragMRI);
                %hTxt = 
                text(coordinates(1), coordinates(2), sprintf('%d',nP), ...
                'Color', 'g', 'FontSize',8, ...
                'HorizontalAlignment','left', 'VerticalAlignment','top');
            end
        end
                
        % Set again the properties lost in axesMRI
        set(handles.axesMRI, 'Tag', tag);
        %set(handles.axes1, 'ButtonDownFcn', @axes1_ButtonDownFcn);
        set(hMRI, 'ButtonDownFcn', {@getClicksMRI, handles});
        
    else % Has not been loaded
        switch get(eventdata.OldValue, 'Tag') % Get Tag of previously selected object
            case 'rbFLAIR' 
                set(handles.rbFLAIR, 'Value', 1);
            case 'rbT1w'
                set(handles.rbT1w, 'Value', 1);
            case 'rbT2w'
                set(handles.rbT2w, 'Value', 1);
            case 'rbT2star'
                set(handles.rbT2star, 'Value', 1);
        end
        errordlg(['The modality ', newModality, ' has not been loaded'], 'Modality not loaded');
    end

end
% *************************************************************************

% --- Executes on button press in pushbuttonRefPoints.
function pushbuttonRefPoints_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRefPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % GET WORKING DIRECTORY
    workDir = get(handles.editWorkDir, 'String');
    
    % GET THE MAT FILE WITH THE POINTS AND LOAD IT
    [filePoints, pathPoints] = uigetfile({'*.mat', 'Mat Files (*.mat)';
    fullfile(workDir,'*.*'),  'All Files (*.*)'}, ...
    'Select a file to save the points', ...
    workDir);
    
    if(filePoints==0)
        errordlg('File not found','File Error');
        return;
    end
    load([pathPoints, filePoints]);

    % CHECK IF THE VARIABLES 'pointsMRI' AND 'pointsHistol' EXIST IN THE
    % WORKSPACE. IF THEY DONT, THE FILE HAS NOT BEEN WELL LOADED
    if(~(exist('pointsMRI','var') && exist('pointsHistol','var')))
        errordlg('This file does not contain points', 'File Error');
        return;
    end
    
    % GET THE TABLE WITH THE POINTS
    dataHandles = guidata(handles.axesMRI);
    
    % CHECK IF ANY MRI IMAGE HAS BEEN LOADED
    indexesMRI = find(~cellfun(@isempty,dataHandles.data.MRIImage(:)));
    if(isempty(indexesMRI) || isempty(dataHandles.data.histologyImage{1}))
        errordlg('You must load the MRI and the Histology image before loading the points', 'Error images loaded');
        return;
    end
    
    % COPY THE POINTS INTO THE VARIABLES dataHandles.data.tableCoordinatesMRI
    % AND dataHandles.data.tableCoordinatesHistol
    dataHandles.data.tableCoordinatesMRI = pointsMRI;
    dataHandles.data.tableCoordinatesHistol = pointsHistol;
    
    % SAVE THE POINTS INTO THE AXIS
    guidata(handles.axesMRI, dataHandles);
    
    % COPY THE POINTS INTO THE GUI TABLES: 
    newDataTableMRI = num2cell(round(pointsMRI));
    newDataTableHistol = num2cell(round(pointsHistol));
    set(handles.uitablePointsMRI, 'Data', newDataTableMRI);
    set(handles.uitablePointsHistol, 'Data', newDataTableHistol);
        
    % IF THERE ARE POINTS IN THE IMAGE, REMOVE THEM AND PAINT THEM AGAIN
    % axesMRI
    axes(handles.axesMRI);
    axesMRI_chil = get(handles.axesMRI, 'Children');
    class_axesMRI = arrayfun(@class, axesMRI_chil, 'UniformOutput', false);
    isImage = strcmpi('matlab.graphics.primitive.Image', class_axesMRI);
    delete(axesMRI_chil(~isImage));
    
    axes(handles.axesMRI);
    for nP_MRI=1:size(pointsMRI,1)
        coord_MRI = [pointsMRI(nP_MRI,2), pointsMRI(nP_MRI,1)]; %(x,y)
        % Next, set your image display axes as the current axes:
        hold on;
        hPoint = plot(coord_MRI(1), coord_MRI(2), 'r*', 'MarkerSize', 9);
        % Put a marker on the point
%         hPoint = impoint(handles.axesMRI, coord_MRI);
        set(hPoint, 'UserData', nP_MRI);
% %         % Construct boundary constraint function
% %         fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
% %         % Enforce boundary constraint function using setPositionConstraintFcn
% %         setPositionConstraintFcn(hPoint, fcn);
%         setColor(hPoint, 'r');
% %         addNewPositionCallback(hPoint, @pointDragMRI);
%         %setString(hPoint, num2str(k));
        %hTxt = 
        text(coord_MRI(1), coord_MRI(2), sprintf('%d',nP_MRI), ...
            'Color', 'r', 'FontSize',8, ...
            'HorizontalAlignment','left', 'VerticalAlignment','top');
    end
    
    % axesHist
    axes(handles.axesHist);
    axesHist_chil = get(handles.axesHist, 'Children');
    class_axesHist = arrayfun(@class, axesHist_chil, 'UniformOutput', false);
    isImage = strcmpi('matlab.graphics.primitive.Image', class_axesHist);
    delete(axesHist_chil(~isImage));
    
    %axes(handles.axesHist);
    for nP_Hist=1:size(pointsHistol,1)
        coord_Hist = [pointsHistol(nP_Hist,2), pointsHistol(nP_Hist,1)]; %(x,y)
        % Put a marker on the point
        hold on;
        hPoint = plot(coord_Hist(1), coord_Hist(2), 'b*', 'MarkerSize', 9);
%         hPoint = impoint(handles.axesHist, coord_Hist);
        set(hPoint, 'UserData', nP_Hist);
% %         % Construct boundary constraint function
% %         fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
% %         % Enforce boundary constraint function using setPositionConstraintFcn
% %         setPositionConstraintFcn(hPoint, fcn);
%         setColor(hPoint, 'b');
% %         addNewPositionCallback(hPoint, @pointDragHistol);
%         %setString(hPoint, num2str(k));
        %hTxt = 
        text(coord_Hist(1), coord_Hist(2), sprintf('%d',nP_Hist), ...
            'Color', 'b', 'FontSize',8, ...
            'HorizontalAlignment','left', 'VerticalAlignment','top');
    end

    % TO-DO: DELETE THE SELECTED POINTS THE REGIONS OF INTEREST ON ITS TABLE
    
    % TO-DO: DELETE THE SELECTED POINTS THE REGIONS OF INTEREST ON THE IMAGES
    
end
% *************************************************************************

% --- Executes on button press in pushbuttonSaveTable.
function pushbuttonSaveTable_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Gets current data of the selected points table
    dataTablePointsLab = get(handles.uitablePointsLab, 'Data');
    dataToSave = [{'MR label','MR row','MR column','Histology label','Histology row','Histology column'};...
        dataTablePointsLab];
        
    % Select path to save the excel file
    [fileName, pathName, filterIndex] = uiputfile({'*.xlsx','Excel Workbook';
        '*.xls', 'Excel 97-2003 Workbook';...
        '*.*',  'All Files (*.*)'}, 'Save selected points file');
    
    % Write file into the selected path
    if(filterIndex~=0)
        try
            if(filterIndex==1)
                xlswrite(fullfile(pathName, fileName), dataToSave);
            elseif(filterIndex==2)
                xlswrite(fullfile(pathName, fileName), dataToSave);
            else
                xlswrite([fullfile(pathName, fileName), '.xls'], dataToSave);
            end
        catch
            [~,name,~] = fileparts(fullfile(pathName, fileName));
            % Save in CSV format
            cell2csv([fullfile(pathName, name),'.csv'], dataToSave);
        end
    end
    
    
end
% *************************************************************************

% --- Executes on mouse press over MRI axes background.
%function getClicks(objectHandle, eventData)
function getClicksMRI(objectHandle, ~, hGlobal)
% objectHandle: Handle to the image itself
% The second input argument is eventData, but it is not used in this case

    axesMRIHandles = get(objectHandle, 'Parent');
    %uiTableMRIHandles = findobj('Tag', 'uitablePointsMRI');
    axesHistHandles = findobj('Tag', 'axesHist');
        
    dataHandles = guidata(axesMRIHandles);
    
    refPointsMRI = dataHandles.data.tableCoordinatesMRI;
    refPointsHistol = dataHandles.data.tableCoordinatesHistol;

    if(~isempty(dataHandles.data.MRIImage{1,1}))
        coordsAux = get(axesMRIHandles, 'CurrentPoint'); 
        pointMRI = coordsAux(1,1:2);

        % Add the coordinates to the uitablePointsMRI
        dataTablePointsLab = get(hGlobal.uitablePointsLab, 'Data'); % Gets current data of the table
        
        k = size(dataTablePointsLab,1)+1; % Row number of the new row of the data
        if(k==1) % It is the first point to insert in the data
            newDataTablePointsLab = num2cell(dataTablePointsLab);
%             dataHandles.data.tableCoordinatesMRI(k,:) = [coordinates(2), coordinates(1)];
        else
            newDataTablePointsLab = dataTablePointsLab;
%             dataHandles.data.tableCoordinatesMRI = [dataHandles.data.tableCoordinatesMRI; [coordinates(2), coordinates(1)]];
        end
        newDataTablePointsLab{k,2} = round(pointMRI(2)); % Row of the MRI point
        newDataTablePointsLab{k,3} = round(pointMRI(1)); % Column of the MRI point
        
        guidata(axesMRIHandles, dataHandles);
        
        % Put a marker on the MRI image
        axes(axesMRIHandles); % Set axes of the MRI the current axis
        hPointMRI = impoint(axesMRIHandles, pointMRI);
        set(hPointMRI, 'UserData', k);
        % Construct boundary constraint function TODO CHECK IF IT IS
        % CORRECT
        fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
        % Enforce boundary constraint function using setPositionConstraintFcn
        setPositionConstraintFcn(hPointMRI, fcn);
        setColor(hPointMRI, 'g');
        % TO-DO: DECIDE WHETHER WE WANT TO LET THE POINTS TO BE DRAGGED OR
        % NOT
        %addNewPositionCallback(hPoint, @pointDragMRI);
        %hTxt = 
        text(pointMRI(1), pointMRI(2), sprintf('%d',k), ...
            'Color', 'g', 'FontSize',8, ...
            'HorizontalAlignment','left', 'VerticalAlignment','top');
        
        % SHOW HELP DIALOG
        hDlg = helpdlg('Estimating point in Histology image.', 'Estimating point');
        
        % FUNCTION TO GET THE SAME POINT IN THE HISTOLOGY IMAGE
        pointHistol = fInferPointInHistol([pointMRI(2), pointMRI(1)], refPointsMRI, refPointsHistol, dataHandles.data.histologyImage{1,3});
        
        % Add coordinaes to the uitablePointsMRI
        newDataTablePointsLab{k,5} = round(pointHistol(2)); % Row of the Histology point
        newDataTablePointsLab{k,6} = round(pointHistol(1)); % Column of the Histology point
        
        % Set the data into the table
        %set(uiTableMRIHandles, 'Data', newDataTableMRI);
        set(hGlobal.uitablePointsLab, 'Data', newDataTablePointsLab);
        
        % Put a marker on the Histology image
        axes(axesHistHandles); % Set axes of the Histoloty the current axis
        hPointHistol = impoint(axesHistHandles, [pointHistol(2) pointHistol(1)]);
        set(hPointHistol, 'UserData', k);
        % Construct boundary constraint function
        fcn = makeConstrainToRectFcn('impoint', get(axesHistHandles,'XLim'), get(axesHistHandles,'YLim'));
        % Enforce boundary constraint function using setPositionConstraintFcn
        setPositionConstraintFcn(hPointHistol, fcn);
        setColor(hPointHistol, 'y');
        % TO-DO: DECIDE WHETHER WE WANT TO LET THE POINTS TO BE DRAGGED OR
        % NOT
        %addNewPositionCallback(hPoint, @pointDragMRI);
        %hTxt = 
        text(pointHistol(2), pointHistol(1), sprintf('%d',k), ...
            'Color', 'y', 'FontSize',8, ...
            'HorizontalAlignment','left', 'VerticalAlignment','top');
        
        % CLOSE HELP DIALOG
        close(hDlg);
    end

end
% *************************************************************************

% --- Executes on mouse press over MRI axes background.
%function getClicks(objectHandle, eventData)
function getClicksHistology(objectHandle, ~, hGlobal)
% objectHandle: Handle to the image itself
% The second input argument is eventData, but it is not used in this case

    axesHistHandles = get(objectHandle, 'Parent');
	axesMRIHandles = findobj('Tag', 'axesMRI');
    dataHandles = guidata(hGlobal.axesMRI);
    
    refPointsMRI = dataHandles.data.tableCoordinatesMRI;
    refPointsHistol = dataHandles.data.tableCoordinatesHistol;

    if(~isempty(dataHandles.data.histologyImage{1,1}))
        coordsAux = get(axesHistHandles, 'CurrentPoint'); 
        pointHistol = coordsAux(1,1:2);
        %disp(['You clicked in point of the Histology axes: (', num2str(coordinates(2)), ',', num2str(coordinates(1)),')']);
        
        % Add the coordinates to the uitablePointsMRI
        %dataTableMRI = get(uiTableMRIHandles, 'Data'); % Gets current data of the table
        dataTablePointsLab = get(hGlobal.uitablePointsLab, 'Data'); % Gets current data of the table
        
        k = size(dataTablePointsLab,1)+1; % Row number of the new row of the data
        
        if(k==1) % It is the first point to insert in the data
            newDataTablePointsLab = num2cell(dataTablePointsLab);
            %dataHandles.data.tableCoordinatesHistol(k,:) = [coordinates(2), coordinates(1)];
        else
            newDataTablePointsLab = dataTablePointsLab;
            %dataHandles.data.tableCoordinatesHistol = [dataHandles.data.tableCoordinatesHistol; [coordinates(2), coordinates(1)]];
        end
        newDataTablePointsLab{k,5} = round(pointHistol(2)); % Row of the point
        newDataTablePointsLab{k,6} = round(pointHistol(1)); % Column of the point
        
        % FUNCTION TO GET THE SAME POINT IN THE MR IMAGE
        pointMRI = fInferPointInMRI([pointHistol(2), pointHistol(1)], refPointsMRI, refPointsHistol, size(dataHandles.data.MRIImage{2,1}));
        
        % Add coordinaes to the uitablePointsMRI
        newDataTablePointsLab{k,2} = round(pointMRI(1)); % Row of the MRI point
        newDataTablePointsLab{k,3} = round(pointMRI(2)); % Column of the MRI point
        
        % Set the data into the table
        %set(uiTableMRIHandles, 'Data', newDataTableMRI);
        set(hGlobal.uitablePointsLab, 'Data', newDataTablePointsLab);
        
        % Save the updated data into axes1
        %guidata(axes1Handles, dataHandles);
        guidata(hGlobal.axesMRI, dataHandles);
        
        % Put a marker on the MR Image
        axes(axesMRIHandles); % Set axes of the MRI the current axis
        hPointMRI = impoint(axesMRIHandles, [pointMRI(2), pointMRI(1)]);
        set(hPointMRI, 'UserData', k);
        % Construct boundary constraint function TODO CHECK IF IT IS
        % CORRECT
        fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
        % Enforce boundary constraint function using setPositionConstraintFcn
        setPositionConstraintFcn(hPointMRI, fcn);
        setColor(hPointMRI, 'g');
        % TO-DO: DECIDE WHETHER WE WANT TO LET THE POINTS TO BE DRAGGED OR
        % NOT
        %addNewPositionCallback(hPoint, @pointDragMRI);
        %hTxt = 
        text(pointMRI(2), pointMRI(1), sprintf('%d',k), ...
            'Color', 'g', 'FontSize',8, ...
            'HorizontalAlignment','left', 'VerticalAlignment','top');
        
        % Put a marker on the Histology image
        axes(axesHistHandles); % Set axes of the Histology the current axis
        hPointHistol = impoint(axesHistHandles, pointHistol);
        set(hPointHistol, 'UserData', k);
        % Construct boundary constraint function
        fcn = makeConstrainToRectFcn('impoint', get(axesHistHandles,'XLim'), get(axesHistHandles,'YLim'));
        % Enforce boundary constraint function using setPositionConstraintFcn
        setPositionConstraintFcn(hPointHistol, fcn);
        setColor(hPointHistol, 'y');
        % TO-DO: DECIDE WHETHER WE WANT TO LET THE POINTS TO BE DRAGGED OR
        % NOT
        %addNewPositionCallback(hPoint, @pointDragMRI);
        %hTxt = 
        text(pointHistol(1), pointHistol(2), sprintf('%d',k), ...
            'Color', 'y', 'FontSize',8, ...
            'HorizontalAlignment','left', 'VerticalAlignment','top');
    end

end
% *************************************************************************

% --- Executes when selected cell(s) is changed in uitablePointsLab.
function uitablePointsLab_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitablePointsLab (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

    if(~isempty(eventdata.Indices) && eventdata.Indices(2)~=1 && eventdata.Indices(2)~=4)

%         tagMRI = get(handles.axesMRI, 'Tag');
%         tagHistol = get(handles.axesHist, 'Tag');
        
        % Get the selected row
        rowToDelete = eventdata.Indices(1);
        %disp(['Selected cell: ', num2str(rowToDelete)]);

        % Get the selected points table from the handles
        dataHandles = guidata(handles.axesMRI);
        dataTablePointsLab = get(handles.uitablePointsLab, 'Data'); % Gets current data of the table
        %pointsMRI = dataHandles.data.tableCoordinatesMRI;
        %pointsHistology = dataHandles.data.tableCoordinatesHistol;

        % Remove points from tables
        dataTablePointsLab(rowToDelete,:) = [];
        %pointsMRI(rowToDelete,:) = [];
        %pointsHistology(rowToDelete,:) = [];

        % Save data into the selected points' table
        %dataHandles.data.tableCoordinatesMRI = pointsMRI;
        %dataHandles.data.tableCoordinatesHistol = pointsHistology;
        %guidata(handles.axes1, dataHandles);

        % Update data in gui table
        set(handles.uitablePointsLab, 'Data', dataTablePointsLab);
        %newDataTableMRI = num2cell(round(pointsMRI));
        %newDataTableHistol = num2cell(round(pointsHistology));
        %set(handles.uitablePointsMRI, 'Data', newDataTableMRI);
        %set(handles.uitablePointsHistol, 'Data', newDataTableHistol);

        % IF THERE ARE REFERENCE POINTS LOADED, REMOVE THEM AND PAINT THEM AGAIN
        % Axes MRI
        axes(handles.axesMRI);
        axes1_chil = get(handles.axesMRI, 'Children');
        class_axes1 = arrayfun(@class, axes1_chil, 'UniformOutput', false);
        isImage = strcmpi('matlab.graphics.primitive.Image', class_axes1);
        delete(axes1_chil(~isImage));

        pointsMRI = dataHandles.data.tableCoordinatesMRI;
        for nP_MRI=1:size(pointsMRI,1)
            coord_MRI = [pointsMRI(nP_MRI,2), pointsMRI(nP_MRI,1)]; %(x,y)
            % Put a marker on the point
            hPoint = plot(coord_MRI(1), coord_MRI(2), 'r*', 'MarkerSize', 9);
            set(hPoint, 'UserData', nP_MRI);
            text(coord_MRI(1), coord_MRI(2), sprintf('%d',nP_MRI), ...
                'Color', 'r', 'FontSize',8, ...
                'HorizontalAlignment','left', 'VerticalAlignment','top');
        end

        % Axes Histology
        axes(handles.axesHist);
        axes2_chil = get(handles.axesHist, 'Children');
        class_axes2 = arrayfun(@class, axes2_chil, 'UniformOutput', false);
        isImage = strcmpi('matlab.graphics.primitive.Image', class_axes2);
        delete(axes2_chil(~isImage));

        pointsHistology = dataHandles.data.tableCoordinatesHistol;
        for nP_Hist=1:size(pointsHistology,1)
            coord_Hist = [pointsHistology(nP_Hist,2), pointsHistology(nP_Hist,1)]; %(x,y)
            % Put a marker on the point
            hPoint = plot(coord_Hist(1), coord_Hist(2), 'b*', 'MarkerSize', 9);
            set(hPoint, 'UserData', nP_Hist);
            text(coord_Hist(1), coord_Hist(2), sprintf('%d',nP_Hist), ...
                'Color', 'b', 'FontSize',8, ...
                'HorizontalAlignment','left', 'VerticalAlignment','top');
        end

        % IF THERE ARE REFERENCE POINTS SELECTED, REMOVE THEM AND PAINT THEM AGAIN
        selPoints = get(handles.uitablePointsLab, 'Data'); % Gets current data of the table
        numSelPointsMRI = size(selPoints,1);
        if(numSelPointsMRI>0)
            axes(handles.axesMRI);
            %coordinatesHistol = [selPoints{nP,6}, selPoints{nP,5}]; %(x,y)
            for nP=1:numSelPointsMRI
                coordinatesMRI = [selPoints{nP,3}, selPoints{nP,2}]; %(x,y)
                
                hPointMRI = impoint(handles.axesMRI, coordinatesMRI);
                set(hPointMRI, 'UserData', nP);
                % Construct boundary constraint function
                fcn = makeConstrainToRectFcn('impoint', get(handles.axesMRI,'XLim'), get(handles.axesMRI,'YLim'));
                % Enforce boundary constraint function using setPositionConstraintFcn
                setPositionConstraintFcn(hPointMRI, fcn);
                setColor(hPointMRI, 'g');
                % TO-DO: DECIDE WHETHER WE WANT TO LET THE POINTS TO BE DRAGGED OR
                % NOT
                %addNewPositionCallback(hPoint, @pointDragMRI);
                text(coordinatesMRI(1), coordinatesMRI(2), sprintf('%d',nP), ...
                'Color', 'g', 'FontSize',8, ...
                'HorizontalAlignment','left', 'VerticalAlignment','top');
            end
            
            axes(handles.axesHist);
            %coordinatesHistol = [selPoints{nP,6}, selPoints{nP,5}]; %(x,y)
            for nP=1:numSelPointsMRI
                coordinatesHistol = [selPoints{nP,6}, selPoints{nP,5}]; %(x,y)
                
                hPointHistol = impoint(handles.axesHist, coordinatesHistol);
                set(hPointHistol, 'UserData', nP);
                % Construct boundary constraint function
                fcn = makeConstrainToRectFcn('impoint', get(handles.axesHist,'XLim'), get(handles.axesHist,'YLim'));
                % Enforce boundary constraint function using setPositionConstraintFcn
                setPositionConstraintFcn(hPointHistol, fcn);
                setColor(hPointHistol, 'y');
                % TO-DO: DECIDE WHETHER WE WANT TO LET THE POINTS TO BE DRAGGED OR
                % NOT
                %addNewPositionCallback(hPoint, @pointDragMRI); 
                text(coordinatesHistol(1), coordinatesHistol(2), sprintf('%d',nP), ...
                'Color', 'y', 'FontSize',8, ...
                'HorizontalAlignment','left', 'VerticalAlignment','top');
            end
            
        end
                
        % Set again the properties lost in axesMRI and axesHistol
%         set(handles.axesMRI, 'Tag', 'axesMRI');
%         set(handles.axesHistol, 'Tag', 'axesHist');
        
        %set(handles.axes1, 'ButtonDownFcn', @axes1_ButtonDownFcn);
%         set(hMRI, 'ButtonDownFcn', {@getClicksMRI, handles});
        
    end
end
