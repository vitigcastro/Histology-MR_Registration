function varargout = Register_PM_MRI(varargin)
% REGISTER_PM_MRI MATLAB code for Register_PM_MRI.fig
%      REGISTER_PM_MRI, by itself, creates a new REGISTER_PM_MRI or raises the existing
%      singleton*.
%
%      H = REGISTER_PM_MRI returns the handle to a new REGISTER_PM_MRI or the handle to
%      the existing singleton*.
%
%      REGISTER_PM_MRI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REGISTER_PM_MRI.M with the given input arguments.
%
%      REGISTER_PM_MRI('Property','Value',...) creates a new REGISTER_PM_MRI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Register_PM_MRI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Register_PM_MRI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% Created by:	Victor Gonzalez Castro
% Funded by:    Row Fogo Charitable Trust
% 

    % Edit the above text to modify the response to help Register_PM_MRI

    % Last Modified by GUIDE v2.5 07-Jun-2016 10:57:37

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @Register_PM_MRI_OpeningFcn, ...
                       'gui_OutputFcn',  @Register_PM_MRI_OutputFcn, ...
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

% --- Executes just before Register_PM_MRI is made visible.
function Register_PM_MRI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Register_PM_MRI (see VARARGIN)

    % ADD B-SPLINE FUNCTIONS TO PATH
    addpath('./B-spline/');
    addpath('./B-spline/functions');
    addpath('./B-spline/functions_affine');
    addpath('./B-spline/functions_nonrigid');
    
    maxMRIImages = 4;
    
    % Choose default command line output for Register_PM_MRI
    handles.output = hObject;
    
    % Add additional data as a new field inside axes1. Specifically:
    % - A cell array to store the MRI images
    data = guidata(handles.axes1);
    data.data.MRIImage = cell(3,maxMRIImages); % Cell array for the MR images
    data.data.tableCoordinatesMRI = zeros(1,2);
    data.data.tableCoordinatesHistol = zeros(1,2);
    %data.indexMRILoaded = zeros(1, maxMRIImages);
    
    % Add additional data as a new field inside axes2. Specifically:
    % - A cell array to store the histology images
    data.data.histologyImage = cell(1,3); % Cell array for the MR images
    
    guidata(handles.axes1, data);
%     guidata(handles.axes2, dataHistology);
    
    % Update handles structure
    %guidata(hObject, handles);

    % UIWAIT makes Register_PM_MRI wait for user response (see UIRESUME)
    % uiwait(handles.figure1);
    
    % Set all buttons of the MRI except for the FLAIR to invisible
    set(handles.buttonT1w, 'Visible', 'off');
    set(handles.buttonT2w, 'Visible', 'off');
    set(handles.buttonT2star, 'Visible', 'off');
    
    % Set string of the text panel above the histology panel
    set(handles.text5, 'String', '');
    
    % Set properties of uitables
    set(handles.uitablePointsMRI, 'ColumnEditable', [false, false]);
    set(handles.uitablePointsHistol, 'ColumnEditable', [false, false]);
    
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
        
end
% *************************************************************************

% --- Outputs from this function are returned to the command line.
function varargout = Register_PM_MRI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    %varargout{1} = handles.output;
end
% *************************************************************************

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end
% *************************************************************************

% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup1 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

end
% *************************************************************************

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of edit2 as text
    %        str2double(get(hObject,'String')) returns contents of edit2 as a double
    
    
end
% *************************************************************************

% --- Executes on button press in buttonFLAIR.
function buttonFLAIR_Callback(hObject, eventdata, handles)
% hObject    handle to buttonFLAIR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    valueFLAIR = 1;

    % Get MRI file name
    [fileImage, pathImage] = uigetfile({'*.dcm;*.jpg;*.jpeg',...
    'Image Files (*.dcm,*.jpg,*.jpeg)';
    '*.dcm',  'DICOM (*.dcm)'; ...
    '*.jpg;*.jpeg','JPEG image (*.jpg;*.jpeg)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Choose the FLAIR image');
        
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
        data = guidata(handles.axes1);
        %data.MRIImage{valueFLAIR} = image;
        data.data.MRIImage{1,valueFLAIR} = image;
        data.data.MRIImage{2,valueFLAIR} = imgCropped;
        data.data.MRIImage{3,valueFLAIR} = maskCrop;
        guidata(handles.axes1, data);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% SHOW IN THE AXIS %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Next, set your image display axes as the current axes:
        axes(handles.axes1);
        
        % Get axes1 tag
        tag = get(handles.axes1, 'Tag');
        
%         % Display image in axes1 using bicubic interpolation
        hFLAIR = imshow(data.data.MRIImage{2,valueFLAIR},[]);
%         resizePos = get(handles.axes1, 'Position');
%         scaleRes = resizePos(3)/max(size(data.data.MRIImage{2,valueFLAIR})); % Axis is squared
%         %imageRes = imresize(data.data.MRIImage{2,valueFLAIR}, [resizePos(4) resizePos(3)]);
%         imageRes = imresize(data.data.MRIImage{2,valueFLAIR}, scaleRes, 'bicubic');
%         hFLAIR = imshow(imageRes);
        
        % Set again the properties lost in axes1
        set(handles.axes1, 'Tag', tag);
        %set(handles.axes1, 'ButtonDownFcn', @axes1_ButtonDownFcn);
        set(hFLAIR, 'ButtonDownFcn', {@getClicksMRI, handles});
        
        % ACTIVATE THE CORRESPONDING RADIOBUTTON
        set(handles.rbFLAIR, 'Value', 1);
        
        % IF THE IMAGE IS WELL READ, SET THE BUTTONS FOR THE OTHER
        % MODALITIES TO VISIBLE
        set(handles.buttonT1w, 'Visible', 'on');
        set(handles.buttonT2w, 'Visible', 'on');
        set(handles.buttonT2star, 'Visible', 'on');
                
    catch
        if(~isempty(fileImage))
            errordlg('Choose a DICOM or a JPEG image', 'Error reading the image');
        end
    end

end
% *************************************************************************

% --- Executes on button press in buttonT1w.
function buttonT1w_Callback(hObject, eventdata, handles)
% hObject    handle to buttonT1w (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    valueT1w = 2;
    
    % Get MRI file name
    [fileImage, pathImage] = uigetfile({'*.dcm;*.jpg;*.jpeg',...
    'Image Files (*.dcm,*.jpg,*.jpeg)';
    '*.dcm',  'DICOM (*.dcm)'; ...
    '*.jpg;*.jpeg','JPEG image (*.jpg;*.jpeg)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Choose the T1w image');
        
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
        data = guidata(handles.axes1);
        
        % Get the mask from the cell array inside axis1
        mask = data.data.MRIImage{3,1};
        
        % Segment the T1W image using the mask (in a external function)
        [imgCropped, maskCrop] = segmentMR_Block(image, 'T1W', mask);
        
        % Save images
        data.data.MRIImage{1, valueT1w} = image;
        data.data.MRIImage{2, valueT1w} = imgCropped;
        data.data.MRIImage{3, valueT1w} = maskCrop;
        guidata(handles.axes1, data);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% SHOW IN THE AXIS %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%        
        % Next, set your image display axes as the current axes:
        axes(handles.axes1);
        
        % Get axes1 tag
        tag = get(handles.axes1, 'Tag');
        
        % Display image in axes1 using bicubic interpolation
        hT1w = imshow(data.data.MRIImage{2,valueT1w},[]);
%         resizePos = get(handles.axes1, 'Position');
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
            for nP=1:size(selPointsMRI,1)
                coordinates = [selPointsMRI(nP,2), selPointsMRI(nP,1)]; %(x,y)
                % Put a marker on the point
                hPoint = impoint(handles.axes1, coordinates);
                set(hPoint, 'UserData', nP);
                % Construct boundary constraint function
                fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
                % Enforce boundary constraint function using setPositionConstraintFcn
                setPositionConstraintFcn(hPoint, fcn);
                setColor(hPoint, 'r');
                addNewPositionCallback(hPoint, @pointDragMRI);
                %setString(hPoint, num2str(k));
                hTxt = text(coordinates(1), coordinates(2), sprintf('%d',nP), ...
                    'Color', 'r', 'FontSize',8, ...
                    'HorizontalAlignment','left', 'VerticalAlignment','top');
            end
        end
        
        % Set again the properties lost in axes1
        set(handles.axes1, 'Tag', tag);
        %set(handles.axes1, 'ButtonDownFcn', @axes1_ButtonDownFcn);
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

% --- Executes on button press in buttonT2w.
function buttonT2w_Callback(hObject, eventdata, handles)
% hObject    handle to buttonT2w (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    valueT2w = 3;
    
    % Get MRI file name
    [fileImage, pathImage] = uigetfile({'*.dcm;*.jpg;*.jpeg',...
    'Image Files (*.dcm,*.jpg,*.jpeg)';
    '*.dcm',  'DICOM (*.dcm)'; ...
    '*.jpg;*.jpeg','JPEG image (*.jpg;*.jpeg)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Choose the T2w image');
        
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
        data = guidata(handles.axes1);
        
        % Get the mask from the cell array inside axis1
        mask = data.data.MRIImage{3,1};
        
        % Segment the T1W image using the mask (in a external function)
        [imgCropped, maskCrop] = segmentMR_Block(image, 'T2W', mask);
        
        % Save images
        data.data.MRIImage{1, valueT2w} = image;
        data.data.MRIImage{2, valueT2w} = imgCropped;
        data.data.MRIImage{3, valueT2w} = maskCrop;
        guidata(handles.axes1, data);
        
     
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% SHOW IN THE AXIS %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%  
        % Next, set your image display axes as the current axes:
        axes(handles.axes1);
        
        % Get axes1 tag
        tag = get(handles.axes1, 'Tag');
        
        % Display image in axes1 using bicubic interpolation
        hT2w = imshow(data.data.MRIImage{2,valueT2w},[]);
%         resizePos = get(handles.axes1, 'Position');
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
            for nP=1:size(selPointsMRI,1)
                coordinates = [selPointsMRI(nP,2), selPointsMRI(nP,1)]; %(x,y)
                % Put a marker on the point
                hPoint = impoint(handles.axes1, coordinates);
                set(hPoint, 'UserData', nP);
                % Construct boundary constraint function
                fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
                % Enforce boundary constraint function using setPositionConstraintFcn
                setPositionConstraintFcn(hPoint, fcn);
                setColor(hPoint, 'r');
                addNewPositionCallback(hPoint, @pointDragMRI);
                %setString(hPoint, num2str(k));
                hTxt = text(coordinates(1), coordinates(2), sprintf('%d',nP), ...
                    'Color', 'r', 'FontSize',8, ...
                    'HorizontalAlignment','left', 'VerticalAlignment','top');
            end
        end
        
        % Set again the properties lost in axes1
        set(handles.axes1, 'Tag', tag);
        %set(handles.axes1, 'ButtonDownFcn', @axes1_ButtonDownFcn);
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

% --- Executes on button press in buttonT2star.
function buttonT2star_Callback(hObject, eventdata, handles)
% hObject    handle to buttonT2star (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    valueT2star = 4;
    
    % Get MRI file name
    [fileImage, pathImage] = uigetfile({'*.dcm;*.jpg;*.jpeg',...
    'Image Files (*.dcm,*.jpg,*.jpeg)';
    '*.dcm',  'DICOM (*.dcm)'; ...
    '*.jpg;*.jpeg','JPEG image (*.jpg;*.jpeg)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Choose the T2star image');
        
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
        data = guidata(handles.axes1);
        
        % Get the mask from the cell array inside axis1
        mask = data.data.MRIImage{3,1};
        
        % Segment the T1W image using the mask (in a external function)
        [imgCropped, maskCrop] = segmentMR_Block(image, 'T2star', mask);
        
        % Save images
        data.data.MRIImage{1, valueT2star} = image;
        data.data.MRIImage{2, valueT2star} = imgCropped;
        data.data.MRIImage{3, valueT2star} = maskCrop;
        guidata(handles.axes1, data);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% SHOW IN THE AXIS %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Next, set your image display axes as the current axes:
        axes(handles.axes1);
        
        % Get axes1 tag
        tag = get(handles.axes1, 'Tag');
        
        % Display image in axes1
        hT2star = imshow(data.data.MRIImage{2,valueT2star},[]);
%         resizePos = get(handles.axes1, 'Position');
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
            for nP=1:size(selPointsMRI,1)
                coordinates = [selPointsMRI(nP,2), selPointsMRI(nP,1)]; %(x,y)
                % Put a marker on the point
                hPoint = impoint(handles.axes1, coordinates);
                set(hPoint, 'UserData', nP);
                % Construct boundary constraint function
                fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
                % Enforce boundary constraint function using setPositionConstraintFcn
                setPositionConstraintFcn(hPoint, fcn);
                setColor(hPoint, 'r');
                addNewPositionCallback(hPoint, @pointDragMRI);
                %setString(hPoint, num2str(k));
                hTxt = text(coordinates(1), coordinates(2), sprintf('%d',nP), ...
                    'Color', 'r', 'FontSize',8, ...
                    'HorizontalAlignment','left', 'VerticalAlignment','top');
            end
        end
        
        % Set again the properties lost in axes1
        set(handles.axes1, 'Tag', tag);
        %set(handles.axes1, 'ButtonDownFcn', @axes1_ButtonDownFcn);
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

% --- Executes on button press in pushbuttonHE.
function pushbuttonHE_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonHE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    valueHistology = 1;
    
    % Get MRI file name
    [fileImage, pathImage] = uigetfile({'*.jpg;*.jpeg;*.png',...
    'Image Files (*.jpg,*.jpeg,*.png)';
    '*.jpg;*.jpeg','JPEG image (*.jpg;*.jpeg)'; ...
    '*.png',  'Portable Network Graphics (*.png)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Choose the H&E histology image');

    try
        % Load image
        imageHistology = imread([pathImage, fileImage]);
        
        % Resize image (TEST: Set scale to 0.25 (i.e. 1/4 of the original
        % size))
        imageHistologyResize = imresize(imageHistology, 0.25);
        sizeHistol = size(imageHistology);
        clear imageHistology;
        
        % ADD TO THE AXIS. IF IT IS NOT THE FIRST ONE, CONTROL THE SLIDER
        % Store in the cell array inside axis1
        data = guidata(handles.axes1);
        %data.histologyImage{valueHistology} = imageHistology;
        data.data.histologyImage{valueHistology,1} = imageHistologyResize;
        % Store also path of the original histology image
        data.data.histologyImage{valueHistology, 2} = [pathImage, fileImage];
        data.data.histologyImage{valueHistology, 3} = sizeHistol;
        guidata(handles.axes1, data);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% SHOW IN THE AXIS %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Next, set your image display axes as the current axes:
        axes(handles.axes2);
        
        % Get axes2 tag
        tag = get(handles.axes2, 'Tag');
        
        % Display image in axes1, using bicubic interpolation
        hHistology = imshow(data.data.histologyImage{valueHistology},[]);
%         resizePos = get(handles.axes2, 'Position');
%         scaleRes = resizePos(3)/max(size(data.data.histologyImage{valueHistology})); % Axis is squared
%         %imageRes = imresize(data.data.MRIImage{2,valueFLAIR}, [resizePos(4) resizePos(3)]);
%         imageRes = imresize(data.data.histologyImage{valueHistology}, scaleRes, 'bicubic');
%         hHistology = imshow(imageRes);
        
        % Set again the properties lost in axes1
        set(handles.axes2, 'Tag', tag);
        %set(handles.axes1, 'ButtonDownFcn', @axes1_ButtonDownFcn);
        set(hHistology, 'ButtonDownFcn', {@getClicksHistology, handles});
        
        % ACTIVATE THE CORRESPONDING RADIOBUTTON
        set(handles.rbHE, 'Value', 1);
        
        % Change text above the histology image panel, indicating the user
        % to select a few points (needed to make an affine registration to 
        % the image)
        set(handles.text5, 'String', 'Select 4-5 points in both images');
        
        % IF THE IMAGE IS WELL READ, MAKE THIS BUTTON INVISIBLE
%         set(handles.buttonMTIse02, 'Visible', 'off');

        set(handles.pushbuttonAffRegistration, 'Enable', 'on');
        set(handles.pushbuttonRegistration, 'Enable', 'off');
        set(handles.radiobutton18, 'Enable', 'on');
        set(handles.radiobutton19, 'Enable', 'on');
        
    catch
        errordlg('Choose a JPEG image', 'Error reading the image');
    end

end
% *************************************************************************

% --- Executes on button press in pushbuttonSavePoints.
function pushbuttonSavePoints_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSavePoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    % Get the table with the points
    dataHandles = guidata(handles.axes1);
    
    pointsMRI = dataHandles.data.tableCoordinatesMRI;
    pointsHistol = dataHandles.data.tableCoordinatesHistol;
    
    % Select file name to save the points
    uisave({'pointsMRI','pointsHistol'}, 'Points.mat');
    
end
% *************************************************************************

% --- Executes on button press in pushbuttonLoadPoints.
function pushbuttonLoadPoints_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLoadPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%     % Disable the cellselectioncallback of the uitables
%     cellSelectionCallbackMRIAux = get(handles.uitablePointsMRI, 'CellSelectionCallback');
%     cellSelectionCallbackHistAux = get(handles.uitablePointsHistol, 'CellSelectionCallback');
%     set(handles.uitablePointsMRI, 'CellSelectionCallback', '');
%     set(handles.uitablePointsHistol, 'CellSelectionCallback', '');
    
    % GET THE MAT FILE WITH THE POINTS AND LOAD IT
    [filePoints, pathPoints] = uigetfile({'*.mat',...
    'Mat Files (*.mat)';
    '*.*',  'All Files (*.*)'}, ...
    'Select a file to save the points');
    
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
    dataHandles = guidata(handles.axes1);
    
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
    guidata(handles.axes1, dataHandles);
    
    % COPY THE POINTS INTO THE GUI TABLES: 
    newDataTableMRI = num2cell(round(pointsMRI));
    newDataTableHistol = num2cell(round(pointsHistol));
    set(handles.uitablePointsMRI, 'Data', newDataTableMRI);
    set(handles.uitablePointsHistol, 'Data', newDataTableHistol);
        
    % IF THERE ARE POINTS IN THE IMAGE, REMOVE THEM AND PAINT THEM AGAIN
    % Axes1
    axes(handles.axes1);
    axes1_chil = get(handles.axes1, 'Children');
    class_axes1 = arrayfun(@class, axes1_chil, 'UniformOutput', false);
    isImage = strcmpi('matlab.graphics.primitive.Image', class_axes1);
    delete(axes1_chil(~isImage));
    
    for nP_MRI=1:size(pointsMRI,1)
        coord_MRI = [pointsMRI(nP_MRI,2), pointsMRI(nP_MRI,1)]; %(x,y)
        % Put a marker on the point
        hPoint = impoint(handles.axes1, coord_MRI);
        set(hPoint, 'UserData', nP_MRI);
        % Construct boundary constraint function
        fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
        % Enforce boundary constraint function using setPositionConstraintFcn
        setPositionConstraintFcn(hPoint, fcn);
        setColor(hPoint, 'r');
        addNewPositionCallback(hPoint, @pointDragMRI);
        %setString(hPoint, num2str(k));
        hTxt = text(coord_MRI(1), coord_MRI(2), sprintf('%d',nP_MRI), ...
            'Color', 'r', 'FontSize',8, ...
            'HorizontalAlignment','left', 'VerticalAlignment','top');
    end
    
    % Axes2
    axes(handles.axes2);
    axes2_chil = get(handles.axes2, 'Children');
    class_axes2 = arrayfun(@class, axes2_chil, 'UniformOutput', false);
    isImage = strcmpi('matlab.graphics.primitive.Image', class_axes2);
    delete(axes2_chil(~isImage));
    
    for nP_Hist=1:size(pointsHistol,1)
        coord_Hist = [pointsHistol(nP_Hist,2), pointsHistol(nP_Hist,1)]; %(x,y)
        % Put a marker on the point
        hPoint = impoint(handles.axes2, coord_Hist);
        set(hPoint, 'UserData', nP_Hist);
        % Construct boundary constraint function
        fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
        % Enforce boundary constraint function using setPositionConstraintFcn
        setPositionConstraintFcn(hPoint, fcn);
        setColor(hPoint, 'b');
        addNewPositionCallback(hPoint, @pointDragHistol);
        %setString(hPoint, num2str(k));
        hTxt = text(coord_Hist(1), coord_Hist(2), sprintf('%d',nP_Hist), ...
            'Color', 'b', 'FontSize',8, ...
            'HorizontalAlignment','left', 'VerticalAlignment','top');
    end
    
%     % Enable again the cellselectioncallback of the uitables
%     set(handles.uitablePointsMRI, 'CellSelectionCallback', cellSelectionCallbackMRIAux);
%     set(handles.uitablePointsHistol, 'CellSelectionCallback', cellSelectionCallbackHistAux)
end
% *************************************************************************

% --- Executes on selection change in listboxPointsMRI.
function listboxPointsMRI_Callback(hObject, eventdata, handles)
% hObject    handle to listboxPointsMRI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxPointsMRI contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxPointsMRI
end
% *************************************************************************

% --- Executes on selection change in listboxPointsHistol.
function listboxPointsHistol_Callback(hObject, eventdata, handles)
% hObject    handle to listboxPointsHistol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxPointsHistol contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxPointsHistol
end
% *************************************************************************

% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup1 
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
    dataHandles = guidata(handles.axes1);
    
    % Look if the modality marked by newVal has already been loaded
    %if(exist('data.MRIImage','var'))
        indexesMRI = find(~cellfun(@isempty,dataHandles.data.MRIImage(:)));
    %else
    %    indexesMRI = 0;
    %end
    
    if(sum(indexesMRI==newValMRI)~=0) % Has been loaded
        %Show in the axis
        % Next, set your image display axes as the current axes:
        axes(handles.axes1);
        
        % Get axes1 tag
        tag = get(handles.axes1, 'Tag');
        
        % Finally, display the image using bicubic interpolation
        hMRI = imshow(dataHandles.data.MRIImage{2,newValMRI},[]);
%         resizePos = get(handles.axes1, 'Position');
%         scaleRes = resizePos(3)/max(size(dataHandles.data.MRIImage{2,newValMRI})); % Axis is squared
%         %imageRes = imresize(data.data.MRIImage{2,valueFLAIR}, [resizePos(4) resizePos(3)]);
%         imageRes = imresize(dataHandles.data.MRIImage{2,newValMRI}, scaleRes, 'bicubic');
%         hMRI = imshow(imageRes);
        
        % CHECK IF THERE ARE POINTS ALREADY SELECTED, AND SHOW THEM IF
        % NECESSARY
        selPointsMRI = dataHandles.data.tableCoordinatesMRI;
        dataTableMRI = get(handles.uitablePointsMRI, 'Data'); % Gets current data of the table
        numPointsMRI = size(dataTableMRI,1);
        if(numPointsMRI>0)
            for nP=1:size(selPointsMRI,1)
                coordinates = [selPointsMRI(nP,2), selPointsMRI(nP,1)]; %(x,y)
                % Put a marker on the point
                hPoint = impoint(handles.axes1, coordinates);
                set(hPoint, 'UserData', nP);
                % Construct boundary constraint function
                fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
                % Enforce boundary constraint function using setPositionConstraintFcn
                setPositionConstraintFcn(hPoint, fcn);
                setColor(hPoint, 'r');
                addNewPositionCallback(hPoint, @pointDragMRI);
                %setString(hPoint, num2str(k));
                hTxt = text(coordinates(1), coordinates(2), sprintf('%d',nP), ...
                    'Color', 'r', 'FontSize',8, ...
                    'HorizontalAlignment','left', 'VerticalAlignment','top');
            end
        end
                
        % Set again the properties lost in axes1
        set(handles.axes1, 'Tag', tag);
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

% --- Executes on mouse press over MRI axes background.
%function getClicks(objectHandle, eventData)
function getClicksMRI(objectHandle, ~, hGlobal)
% objectHandle: Handle to the image itself
% The second input argument is eventData, but it is not used in this case

    axes1Handles = get(objectHandle, 'Parent');
    %uiTableMRIHandles = findobj('Tag', 'uitablePointsMRI');
        
    dataHandles = guidata(axes1Handles);

    if(~isempty(dataHandles.data.MRIImage{1,1}))
        coordsAux = get(axes1Handles, 'CurrentPoint'); 
        coordinates = coordsAux(1,1:2);
%         disp(['You clicked in point of the MRI axes: (', num2str(coordinates(2)), ',', num2str(coordinates(1)),')']);

        % Add the coordinates to the uitablePointsMRI
        %dataTableMRI = get(uiTableMRIHandles, 'Data'); % Gets current data of the table
        dataTableMRI = get(hGlobal.uitablePointsMRI, 'Data'); % Gets current data of the table
        
        k = size(dataTableMRI,1)+1; % Row number of the new row of the data
        if(k==1) % It is the first point to insert in the data
            newDataTableMRI = num2cell(dataTableMRI);
            dataHandles.data.tableCoordinatesMRI(k,:) = [coordinates(2), coordinates(1)];
        else
            newDataTableMRI = dataTableMRI;
            dataHandles.data.tableCoordinatesMRI = [dataHandles.data.tableCoordinatesMRI; [coordinates(2), coordinates(1)]];
        end
        newDataTableMRI{k,1} = round(coordinates(2)); % Row of the point
        newDataTableMRI{k,2} = round(coordinates(1)); % Column of the point
        
        % Set the data into the table
        %set(uiTableMRIHandles, 'Data', newDataTableMRI);
        set(hGlobal.uitablePointsMRI, 'Data', newDataTableMRI);
        
        guidata(axes1Handles, dataHandles);
        
        % Put a marker on the point
        hPoint = impoint(axes1Handles, coordinates);
        set(hPoint, 'UserData', k);
        % Construct boundary constraint function
        fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
        % Enforce boundary constraint function using setPositionConstraintFcn
        setPositionConstraintFcn(hPoint, fcn);
        setColor(hPoint, 'r');
        addNewPositionCallback(hPoint, @pointDragMRI);
        %setString(hPoint, num2str(k));
        hTxt = text(coordinates(1), coordinates(2), sprintf('%d',k), ...
            'Color', 'r', 'FontSize',8, ...
            'HorizontalAlignment','left', 'VerticalAlignment','top');
    end

end
% *************************************************************************

% --- Executes when a draggable point is moved in the MRI axes
function pointDragMRI(pos)
% pos: New position

    % UPDATE TABLE VALUE
    %   Get number of the object
    %k = get(objectHandle, '')
    h = get(gco);
    k = h.UserData;
    % LOOK FOR THE OLD STRING 
    hAxes = get(gco, 'Parent');
    hObjectsOnAxes = get(hAxes, 'Children');
    stringFound = false;
    i=1;
    %for i=1:length(hObjectsOnAxes)
    while(~stringFound)
        if(strcmpi(class(hObjectsOnAxes(i)),'matlab.graphics.primitive.Text') &&...
           strcmpi(get(hObjectsOnAxes(i),'String'),num2str(k)))
           % CHANGE THE POSITION OF THE STRING
           set(hObjectsOnAxes(i),'Position',[pos(1),pos(2)]);
           stringFound = true;
        end
        i = i + 1;
    end
    
    % CHANGE VALUE OF THE K-TH ROW OF THE TABLE
    %%%%% findobj YIELDED A PROBLEM. CHANGED %%%%%
%     uiTableMRIHandles = findobj('Tag', 'uitablePointsMRI');
%     dataTableMRI = get(uiTableMRIHandles, 'Data'); % Gets current data of the table
%     dataTableMRI{k,1} = round(pos(2)); % Row of the point
%     dataTableMRI{k,2} = round(pos(1)); % Column of the point
%     set(uiTableMRIHandles, 'Data', dataTableMRI);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hFigure = get(hAxes, 'Parent');
    uiTableMRIHandles = findObjInObjChildren(hFigure, 'uitablePointsMRI');
    dataTableMRI = get(uiTableMRIHandles, 'Data'); % Gets current data of the table
    dataTableMRI{k,1} = round(pos(2)); % Row of the point
    dataTableMRI{k,2} = round(pos(1)); % Column of the point
    set(uiTableMRIHandles, 'Data', dataTableMRI);
    
    % CHANGE VALUE IN THE DATA STORED IN THE AXES
    dataHandles = guidata(hAxes);
    dataHandles.data.tableCoordinatesMRI(k,:) = [pos(2), pos(1)];
    guidata(hAxes, dataHandles);
    
end
% *************************************************************************

% --- Executes on mouse press over MRI axes background.
function getClicksHistology(objectHandle, ~, hGlobal)
% objectHandle: Handle to the image itself
% The second input argument is eventData, but it is not used in this case

    axes2Handles = get(objectHandle, 'Parent');
%     axes1Handles = findobj('Tag', 'axes1');
%     dataHandles = guidata(axes1Handles);
    dataHandles = guidata(hGlobal.axes1);

    if(~isempty(dataHandles.data.histologyImage{1,1}))
        coordsAux = get(axes2Handles, 'CurrentPoint'); 
        coordinates = coordsAux(1,1:2);
        %disp(['You clicked in point of the Histology axes: (', num2str(coordinates(2)), ',', num2str(coordinates(1)),')']);
        
        % Add the coordinates to the uitablePointsMRI
        %dataTableMRI = get(uiTableMRIHandles, 'Data'); % Gets current data of the table
        dataTableHistol = get(hGlobal.uitablePointsHistol, 'Data'); % Gets current data of the table
        
        k = size(dataTableHistol,1)+1; % Row number of the new row of the data
        
        if(k==1) % It is the first point to insert in the data
            newDataTableHistol = num2cell(dataTableHistol);
            dataHandles.data.tableCoordinatesHistol(k,:) = [coordinates(2), coordinates(1)];
        else
            newDataTableHistol = dataTableHistol;
            dataHandles.data.tableCoordinatesHistol = [dataHandles.data.tableCoordinatesHistol; [coordinates(2), coordinates(1)]];
        end
        newDataTableHistol{k,1} = round(coordinates(2)); % Row of the point
        newDataTableHistol{k,2} = round(coordinates(1)); % Column of the point
        
        % Set the data into the table
        %set(uiTableMRIHandles, 'Data', newDataTableMRI);
        set(hGlobal.uitablePointsHistol, 'Data', newDataTableHistol);
        
        % Save the updated data into axes1
        %guidata(axes1Handles, dataHandles);
        guidata(hGlobal.axes1, dataHandles);
        
        % Put a marker on the point
        hPoint = impoint(axes2Handles, coordinates);
        set(hPoint, 'UserData', k);
        % Construct boundary constraint function
        fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
        % Enforce boundary constraint function using setPositionConstraintFcn
        setPositionConstraintFcn(hPoint, fcn);
        setColor(hPoint, 'b');
        addNewPositionCallback(hPoint, @pointDragHistol);
        hTxt = text(coordinates(1), coordinates(2), sprintf('%d',k), ...
            'Color', 'b', 'FontSize',8, ...
            'HorizontalAlignment','left', 'VerticalAlignment','top');
    end

end
% *************************************************************************

% --- Executes when a draggable point is moved in the MRI axes
function pointDragHistol(pos)
% pos: New position

    % UPDATE TABLE VALUE
    %   Get number of the point
    h = get(gco);
    k = h.UserData;
    % LOOK FOR THE OLD STRING 
    hAxes = get(gco, 'Parent');
    hObjectsOnAxes = get(hAxes, 'Children');
    
    
    stringFound = false;
    i=1;
    %for i=1:length(hObjectsOnAxes)
    while(~stringFound)
        if(strcmpi(class(hObjectsOnAxes(i)),'matlab.graphics.primitive.Text') &&...
           strcmpi(get(hObjectsOnAxes(i),'String'),num2str(k)))
           % CHANGE THE POSITION OF THE STRING
           set(hObjectsOnAxes(i),'Position',[pos(1),pos(2)]);
           stringFound = true;
        end
        i = i + 1;
    end
    
    % CHANGE VALUE OF THE K-TH ROW OF THE TABLE
    %%%%% findobj YIELDED A PROBLEM. CHANGED %%%%%
%     uiTableHistolHandles = findobj('Tag', 'uitablePointsHistol');
%     dataTableHistol = get(uiTableHistolHandles, 'Data'); % Gets current data of the table
%     dataTableHistol{k,1} = round(pos(2)); % Row of the point
%     dataTableHistol{k,2} = round(pos(1)); % Column of the point
%     set(uiTableHistolHandles, 'Data', dataTableHistol);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hFigure = get(hAxes, 'Parent');
    uiTableHistolHandles = findObjInObjChildren(hFigure, 'uitablePointsHistol');
    dataTableHistol = get(uiTableHistolHandles, 'Data'); % Gets current data of the table
    dataTableHistol{k,1} = round(pos(2)); % Row of the point
    dataTableHistol{k,2} = round(pos(1)); % Column of the point
    set(uiTableHistolHandles, 'Data', dataTableHistol);
    
    % CHANGE VALUE IN THE DATA STORED IN THE AXES
    %%%%% findobj YIELDED A PROBLEM. CHANGED %%%%%
%     axes1Handles = findobj('Tag', 'axes1');
%     dataHandles = guidata(axes1Handles);
%     dataHandles.data.tableCoordinatesHistol(k,:) = [pos(2), pos(1)];
%     guidata(hAxes, dataHandles);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    axes1Handles = findObjInObjChildren(hFigure, 'axes1');
    dataHandles = guidata(axes1Handles);
    dataHandles.data.tableCoordinatesHistol(k,:) = [pos(2), pos(1)];
    guidata(hAxes, dataHandles);
    
end
% *************************************************************************

% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%     data = guidata(handles.axes1);
% 
%     if(~isempty(data.MRIImage{1,1}))
%         disp('Hola mundo. Ya hay una imagen cargada');
%     end

end
% *************************************************************************

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbuttonHE.
function pushbuttonHE_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbuttonHE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end
% *************************************************************************

% --- Executes on button press in pushbuttonAffRegistration.
function pushbuttonAffRegistration_Callback(hObject, eventdata, handles)
% When the push button pushbuttonAffRegistration is pressed, the image
% (histology or MRI) is registered to the other with an affine
% registration. Then, the new image is shown in the axes. Finally, the
% user is given the opportunity to save the image. The points in the
% corresponding table are also transformed and changed accordingly
%
% hObject    handle to pushbuttonAffRegistration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    
    % Get data stored in axes1
    dataHandles = guidata(handles.axes1);
    valueHistology = 1;
    
    % NEED TO CHECK IF THERE ARE POINTS LOADED
    if(isempty(dataHandles.data.tableCoordinatesMRI) || isempty(dataHandles.data.tableCoordinatesHistol))
        errordlg('You must have some points loaded before affine registration ', 'Error Affine registration');
        return;
    end
    
    % NEED TO CHECK IF THERE ARE AT LEAST AN MRI IMAGE LOADED
    indexesMRI = find(~cellfun(@isempty,dataHandles.data.MRIImage(:)));
    if(isempty(indexesMRI) || isempty(dataHandles.data.histologyImage{1}))
        errordlg('You must load the MRI and the Histology image before affine registation', 'Error Affine registration');
        return;
    end
    
    % NEED TO CHECK IF THE HISTOLOGY IMAGE HAS BEEN LOADAD
    indexesHistol = find(~cellfun(@isempty,dataHandles.data.histologyImage(:)));
    if(isempty(indexesHistol) || isempty(dataHandles.data.histologyImage{1}))
        errordlg('You must load the MRI and the Histology image before affine registation', 'Error Affine registration');
        return;
    end
        
    % Get what kind of registration to carry out: Histology to MRI or MRI
    % to Histology
    hist2MRI = false;
    hSelectedReg = get(handles.uibuttongroupRegistration, 'SelectedObject');
    switch get(hSelectedReg, 'Tag')
        case 'radiobutton18'
            hist2MRI = true;
        case 'radiobutton19'
            hist2MRI = false;
    end
    
    % Call the registration function
    pointsMRI = [dataHandles.data.tableCoordinatesMRI(:,1), dataHandles.data.tableCoordinatesMRI(:,2)]; %(row,col)
    pointsHistology = [dataHandles.data.tableCoordinatesHistol(:,1), dataHandles.data.tableCoordinatesHistol(:,2)]; %(row,col)
%     MRILoaded = ~cellfun(@isempty,dataHandles.data.MRIImage);
    imagesMRI = dataHandles.data.MRIImage(2,:);
    imagesMRI(cellfun(@isempty,imagesMRI)) = [];
    sizeHistolImageGUI = size(dataHandles.data.histologyImage{1});
%     sizeHistolImageOri = dataHandles.data.histologyImage{3};
%     pathHistolImage = dataHandles.data.histologyImage{2};
    histolImage = imread(dataHandles.data.histologyImage{2});
    % TRY WITH THE ORIGINAL IMAGE SIZE
    histolImage = imresize(histolImage, 0.5);
    clear histolImageOri;
    
    [imRegAffine, pointsReg] = fRegAffine(pointsMRI, pointsHistology, imagesMRI, sizeHistolImageGUI, histolImage, hist2MRI);
    
    if(hist2MRI) % IF THE HISTOLOGY IS REGISTERED TO THE MRI
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% 0 - SAVE THE REGISTERED IMAGE(S) IN THE HARD DISK %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [fileName, pathName, ~] = uiputfile( ...
        {'*.jpg', 'JPG Images (*.jpg)';...
         '*.png', 'PNG Images (*.png)';...
         '*.bmp', 'Bitmap files (*.bmp)';...
         '*.jpg;*.png;*.bmp', 'Image Files (*.jpg,*.png,*.bmp)'},...
         'Save Affine registration image');

        % In case the user does not select a filename.
        if(isequal(fileName,0) && isequal(pathName,0))
            fileName = 'AffineImage.jpg';
            pathName = '.';
        end
        
        fileNameNoExt = fileName(1:end-4);
        fileExt = fileName(end-2:end);
        fileNameAffine = [fileNameNoExt, '_Affine', '.', fileExt];
        imwrite(imRegAffine, [pathName, fileNameAffine]);
        
        dataHandles.data.histologyImage{valueHistology,1} = imRegAffine;
        % Store also the path of the affine registered histology image
        dataHandles.data.histologyImage{valueHistology, 2} = [pathName, fileNameAffine];
        dataHandles.data.histologyImage{valueHistology, 3} = size(imRegAffine);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% 1- CHANGE IMAGE IN AXIS %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes2);
        
        % Get axes2 tag
        tag = get(handles.axes2, 'Tag');
        
        % Display image in axes2 using bicubic interpolation
        hHistology = imshow(imRegAffine,[]);
%         resizePos = get(handles.axes2, 'Position');
%         scaleRes = resizePos(3)/max(size(imRegAffine)); % Axis is squared
%         %imageRes = imresize(data.data.MRIImage{2,valueFLAIR}, [resizePos(4) resizePos(3)]);
%         imageRes = imresize(imRegAffine, scaleRes, 'bicubic');
%         hHistology = imshow(imageRes);
        
        % Set again the properties lost in axes1
        set(handles.axes2, 'Tag', tag);
        %set(handles.axes1, 'ButtonDownFcn', @axes1_ButtonDownFcn);
        set(hHistology, 'ButtonDownFcn', {@getClicksHistology, handles});
        
        % Change text above the histology image panel, indicating the user
        % to select a few points (needed to make an affine registration to 
        % the image)
        set(handles.text5, 'String', 'Now, select points uniformly throughout the image for the final registration.');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% 2- CHANGE LIST OF POINTS AND DRAW POINTS AGAIN %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % COPY THE POINTS INTO THE VARIABLE dataHandles.data.tableCoordinatesHistol
        dataHandles.data.tableCoordinatesHistol = pointsReg;
        
        % COPY THE POINTS INTO THE GUI TABLES: 
        set(handles.uitablePointsHistol, 'Data', num2cell(round(pointsReg)));

        % IF THERE ARE POINTS IN THE IMAGE, REMOVE THEM AND PAINT THEM AGAIN
        % Axes2
        axes(handles.axes2);
        axes2_chil = get(handles.axes2, 'Children');
        class_axes2 = arrayfun(@class, axes2_chil, 'UniformOutput', false);
        isImage = strcmpi('matlab.graphics.primitive.Image', class_axes2);
        delete(axes2_chil(~isImage));

        for nP_Hist=1:size(pointsReg,1)
            coord_Hist = [pointsReg(nP_Hist,2), pointsReg(nP_Hist,1)]; %(x,y)
            % Put a marker on the point
            hPoint = impoint(handles.axes2, coord_Hist);
            set(hPoint, 'UserData', nP_Hist);
            % Construct boundary constraint function
            fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
            % Enforce boundary constraint function using setPositionConstraintFcn
            setPositionConstraintFcn(hPoint, fcn);
            setColor(hPoint, 'b');
            addNewPositionCallback(hPoint, @pointDragHistol);
            %setString(hPoint, num2str(k));
            hTxt = text(coord_Hist(1), coord_Hist(2), sprintf('%d',nP_Hist), ...
                'Color', 'b', 'FontSize',8, ...
                'HorizontalAlignment','left', 'VerticalAlignment','top');
        end
        
        guidata(handles.axes1, dataHandles);
        
    else % IF THE MRI IS REGISTERED TO THE HISTOLOGY
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% 0 - SAVE THE REGISTERED IMAGE(S) IN THE HARD DISK %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [fileName, pathName, ~] = uiputfile( ...
        {'*.jpg', 'JPG Images (*.jpg)';...
         '*.png', 'PNG Images (*.png)';...
         '*.bmp', 'Bitmap files (*.bmp)';...
         '*.jpg;*.png;*.bmp', 'Image Files (*.jpg,*.png,*.bmp)'},...
         'Save Affine registration image');

        % In case the user does not select a filename.
         if(isequal(fileName,0) && isequal(pathName,0))
            fileName = 'AffineImage.jpeg';
            pathName = '.';
         end
        
        fileNameNoExt = fileName(1:end-4);
        fileExt = fileName(end-2:end);
        
        % Save and store only if they were previously loaded
        if(~isempty(dataHandles.data.MRIImage{2,1}))

            imwrite(imRegAffine{1}, [pathName, fileNameNoExt, '_FLAIR_Affine', '.', fileExt]);
            dataHandles.data.MRIImage{2,1} = imRegAffine{1};
        end
        if(~isempty(dataHandles.data.MRIImage{2,2}))
            imwrite(imRegAffine{2}, [pathName, fileNameNoExt, '_T1_Affine', '.', fileExt]);
            dataHandles.data.MRIImage{2,2} = imRegAffine{2};
        end
        if(~isempty(dataHandles.data.MRIImage{2,3}))
            imwrite(imRegAffine{3}, [pathName, fileNameNoExt, '_T2_Affine', '.', fileExt]);
            dataHandles.data.MRIImage{2,3} = imRegAffine{3};
        end
        if(~isempty(dataHandles.data.MRIImage{2,4}))
            imwrite(imRegAffine{4}, [pathName, fileNameNoExt, '_T2star_Affine', '.', fileExt]);
            dataHandles.data.MRIImage{2,4} = imRegAffine{4};
        end
           
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% 1- CHANGE IMAGE IN AXIS %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % LOOK WHICH IS THE MRI IMAGE THAT IS CURRENTLY SHOWN IN AXES 1
        % (LOOK AT THE RADIOBUTTONS)
        valFlair = get(handles.rbFLAIR, 'Value');
        valT1w = get(handles.rbT1w, 'Value');
        valT2w = get(handles.rbT2w, 'Value');
        valT2star = get(handles.rbT2star, 'Value');
        
        idMRI = 0;
        if(valFlair)
            idMRI = 1;
        elseif(valT1w)
            idMRI = 2;
        elseif(valT2w)
            idMRI = 3;
        elseif(valT2star)
            idMRI = 4;
        end
        
        axes(handles.axes1);
        
        % Get axes2 tag
        tag = get(handles.axes1, 'Tag');
        
        % Display image in axes2
        hMRI = imshow(imRegAffine{idMRI},[]);
%         resizePos = get(handles.axes1, 'Position');
%         scaleRes = resizePos(3)/max(size(imRegAffine{idMRI})); % Axis is squared
%         %imageRes = imresize(data.data.MRIImage{2,valueFLAIR}, [resizePos(4) resizePos(3)]);
%         imageRes = imresize(imRegAffine{idMRI}, scaleRes, 'bicubic');
%         hMRI = imshow(imageRes);
        
        % Set again the properties lost in axes1
        set(handles.axes1, 'Tag', tag);
        %set(handles.axes1, 'ButtonDownFcn', @axes1_ButtonDownFcn);
        set(hMRI, 'ButtonDownFcn', {@getClicksMRI, handles});
        
        % Change text above the histology image panel, indicating the user
        % to select a few points (needed to make an affine registration to 
        % the image)
        set(handles.text5, 'String', 'Now, select points uniformly throughout the image for the final registration.');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% 2- CHANGE LIST OF POINTS AND DRAW POINTS AGAIN %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % COPY THE POINTS INTO THE VARIABLES dataHandles.data.tableCoordinatesMRI
        dataHandles.data.tableCoordinatesMRI = pointsReg;
        
        % COPY THE POINTS INTO THE GUI TABLES: 
        set(handles.uitablePointsMRI, 'Data', num2cell(round(pointsReg)));


        % IF THERE ARE POINTS IN THE IMAGE, REMOVE THEM AND PAINT THEM AGAIN
        axes(handles.axes1);
        axes1_chil = get(handles.axes1, 'Children');
        class_axes1 = arrayfun(@class, axes1_chil, 'UniformOutput', false);
        isImage = strcmpi('matlab.graphics.primitive.Image', class_axes1);
        delete(axes1_chil(~isImage));

        for nP_MRI=1:size(pointsReg,1)
            coord_MRI = [pointsReg(nP_MRI,2), pointsReg(nP_MRI,1)]; %(x,y)
            % Put a marker on the point
            hPoint = impoint(handles.axes1, coord_MRI);
            set(hPoint, 'UserData', nP_MRI);
            % Construct boundary constraint function
            fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
            % Enforce boundary constraint function using setPositionConstraintFcn
            setPositionConstraintFcn(hPoint, fcn);
            setColor(hPoint, 'r');
            addNewPositionCallback(hPoint, @pointDragMRI);
            %setString(hPoint, num2str(k));
            hTxt = text(coord_MRI(1), coord_MRI(2), sprintf('%d',nP_MRI), ...
                'Color', 'r', 'FontSize',8, ...
                'HorizontalAlignment','left', 'VerticalAlignment','top');
        end

        guidata(handles.axes1, dataHandles);
    end
    
    set(handles.pushbuttonAffRegistration, 'Enable', 'off');
    set(handles.pushbuttonRegistration, 'Enable', 'on');
    set(handles.radiobutton18, 'Enable', 'off');
    set(handles.radiobutton19, 'Enable', 'off');
    
end
% *************************************************************************

% --- Executes on button press in pushbuttonRegistration.
function pushbuttonRegistration_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRegistration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Get data stored in axes1
    dataHandles = guidata(handles.axes1);

    % NEED TO CHECK IF THERE ARE POINTS LOADED
    if(isempty(dataHandles.data.tableCoordinatesMRI) || isempty(dataHandles.data.tableCoordinatesHistol))
        errordlg('You must have some points loaded before registration ', 'Error Affine registration');
        return;
    end
    
    % NEED TO CHECK IF THERE ARE AT LEAST AN MRI IMAGE LOADED
    indexesMRI = find(~cellfun(@isempty,dataHandles.data.MRIImage(:)));
    if(isempty(indexesMRI) || isempty(dataHandles.data.histologyImage{1}))
        errordlg('You must load the MRI and the Histology image before registation', 'Error Affine registration');
        return;
    end
    
    % NEED TO CHECK IF THE HISTOLOGY IMAGE HAS BEEN LOADAD
    indexesHistol = find(~cellfun(@isempty,dataHandles.data.histologyImage(:)));
    if(isempty(indexesHistol) || isempty(dataHandles.data.histologyImage{1}))
        errordlg('You must load the MRI and the Histology image before registation', 'Error Affine registration');
        return;
    end
        
    % Get what kind of registration to carry out: Histology to MRI or MRI
    % to Histology
    hist2MRI = false;
    hSelectedReg = get(handles.uibuttongroupRegistration, 'SelectedObject');
    switch get(hSelectedReg, 'Tag')
        case 'radiobutton18'
            hist2MRI = true;
        case 'radiobutton19'
            hist2MRI = false;
    end
    
    % Call the registration function
    pointsMRI = [dataHandles.data.tableCoordinatesMRI(:,1), dataHandles.data.tableCoordinatesMRI(:,2)]; %(row,col)
    pointsHistology = [dataHandles.data.tableCoordinatesHistol(:,1), dataHandles.data.tableCoordinatesHistol(:,2)]; %(row,col)
%     MRILoaded = ~cellfun(@isempty,dataHandles.data.MRIImage);
    imagesMRI = dataHandles.data.MRIImage(2,:);
    imagesMRI(cellfun(@isempty,imagesMRI)) = [];
    sizeHistolImageGUI = size(dataHandles.data.histologyImage{1});
%     sizeHistolImageOri = dataHandles.data.histologyImage{3};
%     pathHistolImage = dataHandles.data.histologyImage{2};
    % IF THE HISTOLOGY WAS AFFINE REGISTERED, USE THE AFFINE-REGISTERED 
    % IMAGE STORED IN THE HANDLES
    if(hist2MRI)
        histolImage = dataHandles.data.histologyImage{1};
    else
        % IF THE HISTOLOGY IMAGE WAS NOT AFFINE REGISTERED, READ THE
        % ORIGINAL IMAGE FROM THE DISK.
        % TRY WITH THE ORIGINAL SIZE IMAGE
        histolImage = imread(dataHandles.data.histologyImage{2});
        histolImage = imresize(histolImage, 0.5);
    end
        
    
    %imReg = fRegBSpline(pointsMRI, pointsHistology, imagesMRI, sizeHistolImage, sizeHistolImageOri, pathHistolImage, hist2MRI);
    imRegBSpline = fRegBSpline(pointsMRI, pointsHistology, imagesMRI, sizeHistolImageGUI, histolImage, hist2MRI);
    
%     %Show both registered images
%     if(hist2MRI)
%         figure; 
%         subplot(1,2,1); imshow(imRegBSpline); title('B-Spline registration');
%         subplot(1,2,2); imshow(imRegAffine); title('Affine registration');
%     else
%         figure; 
%         subplot(4,2,1); imshow(imRegBSpline{1}); title('B-Spline registration FLAIR');
%         subplot(4,2,2); imshow(imRegAffine{1}); title('Affine registration FLAIR');
%         subplot(4,2,3); imshow(imRegBSpline{2}); title('B-Spline registration T1');
%         subplot(4,2,4); imshow(imRegAffine{2}); title('Affine registration T1');
%         subplot(4,2,5); imshow(imRegBSpline{3}); title('B-Spline registration T2');
%         subplot(4,2,6); imshow(imRegAffine{3}); title('Affine registration T2');
%         subplot(4,2,7); imshow(imRegBSpline{4}); title('B-Spline registration T2star');
%         subplot(4,2,8); imshow(imRegAffine{4}); title('Affine registration T2star');
%     end
    
%     choice = questdlg('Which image do you want to save?', ...
% 	'Save the image', ...
% 	'B-Spline', 'Affine', 'Both', 'Both');
    
    % Save the image
    [fileName, pathName, ~] = uiputfile( ...
    {'*.jpg', 'JPG Images (*.jpg)';...
     '*.png', 'PNG Images (*.png)';...
     '*.bmp', 'Bitmap files (*.bmp)';...
     '*.jpg;*.png;*.bmp', 'Image Files (*.jpg,*.png,*.bmp)'},...
     'Save registered images');
    
    fileNameNoExt = fileName(1:end-4);
    fileExt = fileName(end-2:end);
    
    if(~(isequal(fileName,0) || isequal(pathName,0)))
%         if(strcmpi(choice,'B-Spline') || strcmpi(choice, 'Both'))
        if(hist2MRI)
            imwrite(imRegBSpline, [pathName, fileNameNoExt, '_BSpline', '.', fileExt]);
        else
            imwrite(imRegBSpline{1}, [pathName, fileNameNoExt, '_FLAIR_BSpline', '.', fileExt]);
            imwrite(imRegBSpline{2}, [pathName, fileNameNoExt, '_T1_BSpline', '.', fileExt]);
            imwrite(imRegBSpline{3}, [pathName, fileNameNoExt, '_T2_BSpline', '.', fileExt]);
            imwrite(imRegBSpline{4}, [pathName, fileNameNoExt, '_T2star_BSpline', '.', fileExt]);
        end
%         end
%         if(strcmpi(choice,'Affine') || strcmpi(choice, 'Both'))
%             if(hist2MRI)
%                 imwrite(imRegAffine, [pathName, fileNameNoExt, '_Affine', '.', fileExt]);
%             else
%                 imwrite(imRegAffine{1}, [pathName, fileNameNoExt, '_FLAIR_Affine', '.', fileExt]);
%                 imwrite(imRegAffine{2}, [pathName, fileNameNoExt, '_T1_Affine', '.', fileExt]);
%                 imwrite(imRegAffine{3}, [pathName, fileNameNoExt, '_T2_Affine', '.', fileExt]);
%                 imwrite(imRegAffine{4}, [pathName, fileNameNoExt, '_T2star_Affine', '.', fileExt]);
%             end
%         end
    end
    
end
% *************************************************************************

% --- Executes when selected cell(s) is changed in uitablePointsMRI.
function uitablePointsMRI_CellSelectionCallback(hObject, eventdata, handles)
% When a cell is clicked in a table, that row is deleted in that table and
% in the other table. Then the tables are updated and the points are drawn
% again.
%
% hObject    handle to uitablePointsMRI (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)


    if(~isempty(eventdata.Indices))
%     cellSelectionCallbackMRIAux = get(handles.uitablePointsMRI, 'CellSelectionCallback');
%     cellSelectionCallbackHistAux = get(handles.uitablePointsHistol, 'CellSelectionCallback');
%     set(handles.uitablePointsMRI, 'CellSelectionCallback', '');
%     set(handles.uitablePointsHistol, 'CellSelectionCallback', '');

        % Get the selected row
        rowToDelete = eventdata.Indices(1);
        %disp(['Selected cell: ', num2str(rowToDelete)]);

        % Get the tables from the handles
        dataHandles = guidata(handles.axes1);
        pointsMRI = dataHandles.data.tableCoordinatesMRI;
        pointsHistology = dataHandles.data.tableCoordinatesHistol;

        % Remove points from tables
        pointsMRI(rowToDelete,:) = [];
        pointsHistology(rowToDelete,:) = [];

        % Save tables in handles
        dataHandles.data.tableCoordinatesMRI = pointsMRI;
        dataHandles.data.tableCoordinatesHistol = pointsHistology;
        guidata(handles.axes1, dataHandles);

        % Update data in gui tables
        newDataTableMRI = num2cell(round(pointsMRI));
        newDataTableHistol = num2cell(round(pointsHistology));
        set(handles.uitablePointsMRI, 'Data', newDataTableMRI);
        set(handles.uitablePointsHistol, 'Data', newDataTableHistol);

        % IF THERE ARE POINTS IN THE IMAGE, REMOVE THEM AND PAINT THEM AGAIN
        % Axes1
        axes(handles.axes1);
        axes1_chil = get(handles.axes1, 'Children');
        class_axes1 = arrayfun(@class, axes1_chil, 'UniformOutput', false);
        isImage = strcmpi('matlab.graphics.primitive.Image', class_axes1);
        delete(axes1_chil(~isImage));

        for nP_MRI=1:size(pointsMRI,1)
            coord_MRI = [pointsMRI(nP_MRI,2), pointsMRI(nP_MRI,1)]; %(x,y)
            % Put a marker on the point
            hPoint = impoint(handles.axes1, coord_MRI);
            set(hPoint, 'UserData', nP_MRI);
            % Construct boundary constraint function
            fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
            % Enforce boundary constraint function using setPositionConstraintFcn
            setPositionConstraintFcn(hPoint, fcn);
            setColor(hPoint, 'r');
            addNewPositionCallback(hPoint, @pointDragMRI);
            %setString(hPoint, num2str(k));
            hTxt = text(coord_MRI(1), coord_MRI(2), sprintf('%d',nP_MRI), ...
                'Color', 'r', 'FontSize',8, ...
                'HorizontalAlignment','left', 'VerticalAlignment','top');
        end

        % Axes2
        axes(handles.axes2);
        axes2_chil = get(handles.axes2, 'Children');
        class_axes2 = arrayfun(@class, axes2_chil, 'UniformOutput', false);
        isImage = strcmpi('matlab.graphics.primitive.Image', class_axes2);
        delete(axes2_chil(~isImage));

        for nP_Hist=1:size(pointsHistology,1)
            coord_Hist = [pointsHistology(nP_Hist,2), pointsHistology(nP_Hist,1)]; %(x,y)
            % Put a marker on the point
            hPoint = impoint(handles.axes2, coord_Hist);
            set(hPoint, 'UserData', nP_Hist);
            % Construct boundary constraint function
            fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
            % Enforce boundary constraint function using setPositionConstraintFcn
            setPositionConstraintFcn(hPoint, fcn);
            setColor(hPoint, 'b');
            addNewPositionCallback(hPoint, @pointDragHistol);
            %setString(hPoint, num2str(k));
            hTxt = text(coord_Hist(1), coord_Hist(2), sprintf('%d',nP_Hist), ...
                'Color', 'b', 'FontSize',8, ...
                'HorizontalAlignment','left', 'VerticalAlignment','top');
        end
    
%     set(handles.uitablePointsMRI, 'CellSelectionCallback', cellSelectionCallbackMRIAux);
%     set(handles.uitablePointsHistol, 'CellSelectionCallback', cellSelectionCallbackHistAux);
    end
end
% *************************************************************************

% --- Executes when selected cell(s) is changed in uitablePointsHistol.
function uitablePointsHistol_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitablePointsHistol (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

    if(~isempty(eventdata.Indices))
%     cellSelectionCallbackMRIAux = get(handles.uitablePointsMRI, 'CellSelectionCallback');
%     cellSelectionCallbackHistAux = get(handles.uitablePointsHistol, 'CellSelectionCallback');
%     set(handles.uitablePointsMRI, 'CellSelectionCallback', '');
%     set(handles.uitablePointsHistol, 'CellSelectionCallback', '');
    
        % Get the selected row
        rowToDelete = eventdata.Indices(1);
        %disp(['Selected cell: ', num2str(rowToDelete)]);

        % Get the tables from the handles
        dataHandles = guidata(handles.axes1);
        pointsMRI = dataHandles.data.tableCoordinatesMRI;
        pointsHistology = dataHandles.data.tableCoordinatesHistol;

        % Remove points from tables
        pointsMRI(rowToDelete,:) = [];
        pointsHistology(rowToDelete,:) = [];

        % Save tables in handles
        dataHandles.data.tableCoordinatesMRI = pointsMRI;
        dataHandles.data.tableCoordinatesHistol = pointsHistology;
        guidata(handles.axes1, dataHandles);

        % Update data in gui tables
        newDataTableMRI = num2cell(round(pointsMRI));
        newDataTableHistol = num2cell(round(pointsHistology));
        set(handles.uitablePointsMRI, 'Data', newDataTableMRI);
        set(handles.uitablePointsHistol, 'Data', newDataTableHistol);

        % IF THERE ARE POINTS IN THE IMAGE, REMOVE THEM AND PAINT THEM AGAIN
        % Axes1
        axes(handles.axes1);
        axes1_chil = get(handles.axes1, 'Children');
        class_axes1 = arrayfun(@class, axes1_chil, 'UniformOutput', false);
        isImage = strcmpi('matlab.graphics.primitive.Image', class_axes1);
        delete(axes1_chil(~isImage));

        for nP_MRI=1:size(pointsMRI,1)
            coord_MRI = [pointsMRI(nP_MRI,2), pointsMRI(nP_MRI,1)]; %(x,y)
            % Put a marker on the point
            hPoint = impoint(handles.axes1, coord_MRI);
            set(hPoint, 'UserData', nP_MRI);
            % Construct boundary constraint function
            fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
            % Enforce boundary constraint function using setPositionConstraintFcn
            setPositionConstraintFcn(hPoint, fcn);
            setColor(hPoint, 'r');
            addNewPositionCallback(hPoint, @pointDragMRI);
            %setString(hPoint, num2str(k));
            hTxt = text(coord_MRI(1), coord_MRI(2), sprintf('%d',nP_MRI), ...
                'Color', 'r', 'FontSize',8, ...
                'HorizontalAlignment','left', 'VerticalAlignment','top');
        end

        % Axes2
        axes(handles.axes2);
        axes2_chil = get(handles.axes2, 'Children');
        class_axes2 = arrayfun(@class, axes2_chil, 'UniformOutput', false);
        isImage = strcmpi('matlab.graphics.primitive.Image', class_axes2);
        delete(axes2_chil(~isImage));

        for nP_Hist=1:size(pointsHistology,1)
            coord_Hist = [pointsHistology(nP_Hist,2), pointsHistology(nP_Hist,1)]; %(x,y)
            % Put a marker on the point
            hPoint = impoint(handles.axes2, coord_Hist);
            set(hPoint, 'UserData', nP_Hist);
            % Construct boundary constraint function
            fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim'), get(gca,'YLim'));
            % Enforce boundary constraint function using setPositionConstraintFcn
            setPositionConstraintFcn(hPoint, fcn);
            setColor(hPoint, 'b');
            addNewPositionCallback(hPoint, @pointDragHistol);
            %setString(hPoint, num2str(k));
            hTxt = text(coord_Hist(1), coord_Hist(2), sprintf('%d',nP_Hist), ...
                'Color', 'b', 'FontSize',8, ...
                'HorizontalAlignment','left', 'VerticalAlignment','top');
        end

%     set(handles.uitablePointsMRI, 'CellSelectionCallback', cellSelectionCallbackMRIAux);
%     set(handles.uitablePointsHistol, 'CellSelectionCallback', cellSelectionCallbackHistAux);
    end
end
% *************************************************************************

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbuttonLoadPoints.
function pushbuttonLoadPoints_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbuttonLoadPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end
% *************************************************************************

