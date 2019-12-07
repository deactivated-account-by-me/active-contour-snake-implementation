function varargout = Snakes_GUIDE(varargin)
% SNAKES_GUIDE MATLAB code for Snakes_GUIDE.fig
%      SNAKES_GUIDE, by itself, creates a new SNAKES_GUIDE or raises the existing
%      singleton*.
%
%      H = SNAKES_GUIDE returns the handle to a new SNAKES_GUIDE or the handle to
%      the existing singleton*.
%
%      SNAKES_GUIDE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SNAKES_GUIDE.M with the given input arguments.
%
%      SNAKES_GUIDE('Property','Value',...) creates a new SNAKES_GUIDE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Snakes_GUIDE_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Snakes_GUIDE_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Snakes_GUIDE

% Last Modified by GUIDE v2.5 03-Dec-2019 22:16:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Snakes_GUIDE_OpeningFcn, ...
                   'gui_OutputFcn',  @Snakes_GUIDE_OutputFcn, ...
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

function UploadImage(hObject, eventdata, handles,imageFile)
    try
        handles.Image = imread(imageFile);
    catch ME
        % If problem reading image, display error message
        %uialert(app.UIFigure, ME.message, 'Image Error');
        return;
    end
    % axes(handles.Image_Out)
    imshow(handles.Image, 'Parent', handles.Image_Out)
    %hold(app.UIAxes, "on");

    if (ndims(handles.Image) == 3)
        handles.Image = rgb2gray(handles.Image);
    end

    [handles.x,handles.y,handles.x2, handles.y2] = InitializeUISnake(hObject, eventdata, handles);
    
    guidata(hObject, handles);
end

function [x, y, x2, y2] = InitializeUISnake(hObject, eventdata, handles)
    
    snake = handles.snake_type;
    
    if get(handles.closed,'Value') == 1
      handles.snake_type = "Closed";
      snake = "Closed";
    else
      handles.snake_type = "Open";
      snake = "Open";
    end  
    I = handles.Image;
    disp(snake)
    max_len = max(size(I)) - 1;

    [x, y] = getpts();

    hold(handles.Image_Out, "on");

    x = transpose(x);
    y = transpose(y);

    x2 = x;
    y2 = y+5;
    
    y= y-5;
    
    if(snake == "Closed")
        x = [x, x(1)];
        y = [y, y(1)];
        knots = [x ; y];
        number_of_points = length(x);
        distance_points = 1:number_of_points;
        final_distance_points = 1:0.05:number_of_points;
        closed_curve = spline(distance_points, knots, final_distance_points);
        closed_curve(closed_curve < 1) = 1;
        closed_curve(closed_curve > max_len) = max_len;
        x_new = closed_curve(1,:);
        y_new = closed_curve(2,:);
        plot(handles.Image_Out, x ,y, 'o', x_new, y_new, '--');
        x = x_new;
        y = y_new;
    else
        knots = [x ; y];
        knots2 = [x2; y2];
        
        number_of_points = length(x);
        distance_points = 1:number_of_points;
        final_distance_points = 1:0.09:number_of_points;
        
        Open_Curve = spline(distance_points, knots, final_distance_points);
        Open_Curve(Open_Curve < 1) = 1;
        Open_Curve(Open_Curve > max_len) = max_len;
        
        x_new = Open_Curve(1,:);
        y_new = Open_Curve(2,:);
        
        plot(handles.Image_Out,x ,y, 'o', x_new, y_new, '--');
        
        x = x_new;
        y = y_new;
        
        Open_Curve2 = spline(distance_points, knots2, final_distance_points);
        Open_Curve2(Open_Curve2 < 1) = 1;
        Open_Curve2(Open_Curve2 > max_len) = max_len;
        
        x_new2 = Open_Curve2(1,:);
        y_new2 = Open_Curve2(2,:);
        
        plot(handles.Image_Out,x2 ,y2, 'o', x_new2, y_new2, '--');
        
        x2 = x_new2;
        y2 = y_new2;
        
        
        
    end
%     
%     handles.x = x;
%     handles.y = y;
%     
%     handles.x2 = x2;
%     handles.y2 = y2;
    
    hold(handles.Image_Out, "on");
    
end

function Snake(hObject, eventdata, handles)

    N = handles.iterations;
    snakeType = handles.snake_type;
    alpha = handles.alpha;
    beta = handles.beta;
    sigma = handles.sigma;
    gamma = handles.gamma;
    wLine = handles.wLine;
    wEdge = handles.wEdge;
    wTerm = handles.wTerm;

    x1 = handles.x;
    y1 = handles.y;
    
    x2 = handles.x2;
    y2 = handles.y2;
    
    I = handles.Image;

    edgeFunction = handles.edge_function;

    I_after_gaussian_filter =  double(imgaussfilt(I, sigma));
    external_energy = ExternalEnergyCal(I_after_gaussian_filter, wLine, wEdge, wTerm);

    % calcualte the inverse matrix of A
    disp("points "+size(x1))
    a_inverse = InternalEnergyCal(size(x1, 2), alpha, beta, gamma, snakeType);
    x1 = x1';
    y1 = y1';
    
    x2 = x2';
    y2 = y2';

    edge_x = [1 0 -1;2 0 -2; 1 0 -1];
    edge_y = [1 2 1; 0 0 0; -1 -2 -1];

    if(edgeFunction == "Sobel")
        edge_x = [1 0 -1;2 0 -2; 1 0 -1];
        %edge_x = [0 0 0; 0 0 0; 0 0 0];
        edge_y = [1 2 1; 0 0 0; -1 -2 -1];
        %edge_y = [0 0 0; 0 0 0; 0 0 0];
    end

    if(edgeFunction == "Prewitt")
        edge_x = [1 0 -1;1 0 -1; 1 0 -1];
        edge_y = [1 1 1; 0 0 0; -1 -1 -1];
    end
     if(edgeFunction == "Roberts")
        edge_x = [1 0 0; 0 -1 0; 0 0 0];
        edge_y = [0 1 0; -1 0 0; 0 0 0];
     end

    fx = conv2(external_energy, edge_x, 'same');
    fy = conv2(external_energy, edge_y, 'same');

    steps = floor(N/30);

    set(gcf,'WindowButtonDownFcn',{@ButttonDownFcn});
    global click;
    global click_x;
    global click_y;
    global click_x2;
    global click_y2;
    global index_new_point;
    global constraint_type;
    global snake_constraint;
    global index_new_point2;

    disp(constraint_type)
    for i = 1:N
        if(click == 1)
          if (snake_constraint == 1)
            distance_measurement_vector = euclidean_distance([x1 y1], [click_x click_y]);
          else
            distance_measurement_vector = euclidean_distance([x2 y2], [click_x click_y]);
          end
          
          close_point = min(distance_measurement_vector);
          if(constraint_type == "Soft")
              if (snake_constraint == 1)
                index = find(distance_measurement_vector == close_point);
                force_x = click_x - x1(index);
                force_y = click_y - y1(index);
                x1(index) = handles.k * force_x + x1(index);
                y1(index) = handles.k * force_y + y1(index);
              else
                index = find(distance_measurement_vector == close_point);
                force_x = click_x2 - x2(index);
                force_y = click_y2 - y2(index);
                x2(index) = handles.k * force_x + x2(index);
                y2(index) = handles.k * force_y + y2(index);
                  
                  
              end
%             index = find(distance_measurement_vector == close_point);
%             force_x = click_x - x1(index);
%             force_y = click_y - y1(index);
%             x1(index) = handles.k * force_x + x1(index);
%             y1(index) = handles.k * force_y + y1(index);
          else
              if (snake_constraint == 1)
                index_new_point = find(distance_measurement_vector == close_point);
                force_x = click_x - x1(index_new_point);
                force_y = click_y - y1(index_new_point);
                x1(index_new_point) = handles.k * force_x + x1(index_new_point);
                y1(index_new_point) = handles.k * force_y + y1(index_new_point);
              else
                index_new_point2 = find(distance_measurement_vector == close_point);
                force_x = click_x2 - x2(index_new_point2);
                force_y = click_y2 - y2(index_new_point2);
                x2(index_new_point2) = handles.k * force_x + x2(index_new_point2);
                y2(index_new_point2) = handles.k * force_y + y2(index_new_point2);
                  
              end
%             index_new_point = find(distance_measurement_vector == close_point);
%             force_x = click_x - x1(index_new_point);
%             force_y = click_y - y1(index_new_point);
%             x1(index_new_point) = handles.k * force_x + x1(index_new_point);
%             y1(index_new_point) = handles.k * force_y + y1(index_new_point);
          end
        end
        [x1, y1] = iteration(a_inverse, x1, y1, external_energy, gamma, fx, fy);
        [x2, y2] = iteration(a_inverse, x2, y2, external_energy, gamma, fx, fy);
        
        y2 = max(y2, y1+5);
        
        imshow(handles.Image,'parent', handles.Image_Out);
        hold(handles.Image_Out, "on");
        plot(handles.Image_Out,x1, y1, 'r');
        plot(handles.Image_Out,x2, y2, 'g');
        if(mod(i, steps) == 0)
            fprintf('%d/%d interations\n', i, N);
        end
        click = 0;
        pause(0.1);
    end

    if(mod(i, steps) == 0)
        fprintf('%d/%d interations\n', N, N);
    end
end

function distance_array = euclidean_distance(vector_1, vector_2)
    d = abs(vector_1 - vector_2);
    distance_array = [];
    for k1 = 1:length(d)
        distance = norm(d(k1));
        distance_array = [distance_array, distance];
    end
end

function ButttonDownFcn(src,event)
    global click;
    global click_x;
    global click_y;
    global click_x2;
    global click_y2;
    global first_click;
    global snake_constraint;
    first_click =1;
    pt = get(gca,'CurrentPoint');
    if(snake_constraint == 1)
        click_x = pt(1,1);
        click_y = pt(1,2);
        click = 1;
end

% --- Executes just before Snakes_GUIDE is made visible.
function Snakes_GUIDE_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Snakes_GUIDE (see VARARGIN)

% Choose default command line output for Snakes_GUIDE
handles.output = hObject;

handles.iterations = str2num(get(handles.Iteration_Val,'String'));
handles.alpha = str2num(get(handles.Alpha_Val,'String'));
handles.beta = str2num(get(handles.Beta_Val,'String'));
handles.sigma = str2num(get(handles.Sigma_Val,'String'));
handles.gamma = str2num(get(handles.Gamma_Val,'String'));
handles.wLine = str2num(get(handles.WLine_Val,'String'));
handles.wEdge = str2num(get(handles.WEdge_Val,'String'));
handles.wTerm = str2num(get(handles.WTerm_Val,'String'));


if get(handles.sobel,'Value') == 1
  handles.edge_function = "Sobel";
elseif get(handles.prewitt,'Value') == 1
  handles.edge_function = "Prewitt";
elseif get(handels.roberts,'Value') == 1
  handles.edge_function = "Roberts";
end


if get(handles.closed,'Value') == 1
  handles.snake_type = "Closed";
else
  handles.snake_type = "Open";
end

if get(handles.soft,'Value') == 1
  handles.constraint_type = "Soft";
else
  handles.constraint_type = "Hard";
end

if get(handles.snake1,'Value') == 1
  handles.snake_constraint = 1;
else
  handles.snake_constraint = 2;
end

guidata(hObject, handles);
end

% UIWAIT makes Snakes_GUIDE wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = Snakes_GUIDE_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

% --- Executes on button press in Cancel_Button.
function Cancel_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --- Executes on button press in Start_Button.
function Start_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Start_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)\

handles.iterations = str2num(get(handles.Iteration_Val,'String'));
handles.alpha = str2num(get(handles.Alpha_Val,'String'));
handles.beta = str2num(get(handles.Beta_Val,'String'));
handles.sigma = str2num(get(handles.Sigma_Val,'String'));
handles.gamma = str2num(get(handles.Gamma_Val,'String'));
handles.wLine = str2num(get(handles.WLine_Val,'String'));
handles.wEdge = str2num(get(handles.WEdge_Val,'String'));
handles.wTerm = str2num(get(handles.WTerm_Val,'String'));


if get(handles.sobel,'Value') == 1
  handles.edge_function = "Sobel";
elseif get(handles.prewitt,'Value') == 1
  handles.edge_function = "Prewitt";
else
  handles.edge_function = "Roberts";
end

if get(handles.closed,'Value') == 1
  handles.snake_type = "Closed";
else
  handles.snake_type = "Open";
end

if get(handles.soft,'Value') == 1
  handles.constraint_type = "Soft";
else
  handles.constraint_type = "Hard";
end

if get(handles.snake1,'Value') == 1
  handles.snake_constraint = 1;
else
  handles.snake_constraint = 2;
end




global click;
global fix_point;
global first_click;
global constraint_type;
global snake_constraint;

click = 0;
fix_point = double.empty(3,0);
first_click = 0;
constraint_type = handles.constraint_type;
snake_constraint = handles.snake_constraint;

handles.k =1.25;

Snake(hObject, eventdata, handles);
guidata(hObject, handles);
end

% --- Executes on button press in Select_Image.
function Select_Image_Callback(hObject, eventdata, handles)
% hObject    handle to Select_Image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end


function Iteration_Val_Callback(hObject, eventdata, handles)
% hObject    handle to Iteration_Val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Iteration_Val as text
%        str2double(get(hObject,'String')) returns contents of Iteration_Val as a double
end

% --- Executes during object creation, after setting all properties.
function Iteration_Val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Iteration_Val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function Alpha_Val_Callback(hObject, eventdata, handles)
% hObject    handle to Alpha_Val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Alpha_Val as text
%        str2double(get(hObject,'String')) returns contents of Alpha_Val as a double
end

% --- Executes during object creation, after setting all properties.
function Alpha_Val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Alpha_Val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function Beta_Val_Callback(hObject, eventdata, handles)
% hObject    handle to Beta_Val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Beta_Val as text
%        str2double(get(hObject,'String')) returns contents of Beta_Val as a double
end

% --- Executes during object creation, after setting all properties.
function Beta_Val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Beta_Val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function Gamma_Val_Callback(hObject, eventdata, handles)
% hObject    handle to Gamma_Val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Gamma_Val as text
%        str2double(get(hObject,'String')) returns contents of Gamma_Val as a double
end

% --- Executes during object creation, after setting all properties.
function Gamma_Val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gamma_Val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function Sigma_Val_Callback(hObject, eventdata, handles)
% hObject    handle to Sigma_Val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Sigma_Val as text
%        str2double(get(hObject,'String')) returns contents of Sigma_Val as a double
end

% --- Executes during object creation, after setting all properties.
function Sigma_Val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Sigma_Val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function WEdge_Val_Callback(hObject, eventdata, handles)
% hObject    handle to WEdge_Val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WEdge_Val as text
%        str2double(get(hObject,'String')) returns contents of WEdge_Val as a double

% --- Executes during object creation, after setting all properties.
end


function WLine_Val_Callback(hObject, eventdata, handles)
% hObject    handle to WLine_Val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WLine_Val as text
%        str2double(get(hObject,'String')) returns contents of WLine_Val as a double
end


function WLine_Val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WLine_Val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes during object creation, after setting all properties.
function WEdge_Val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WEdge_Val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function WTerm_Val_Callback(hObject, eventdata, handles)
% hObject    handle to WTerm_Val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WTerm_Val as text
%        str2double(get(hObject,'String')) returns contents of WTerm_Val as a double
end

% --- Executes during object creation, after setting all properties.
function WTerm_Val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WTerm_Val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on key press with focus on Select_Image and none of its controls.
function Select_Image_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to Select_Image (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
     % Display uigetfile dialog
    filterspec = {'*.jpg;*.tif;*.png;*.gif','All Image Files'};
    [f, p] = uigetfile(filterspec);

    % Make sure user didn't cancel uigetfile dialog
    if (ischar(p))
       fname = [p f];
       UploadImage(hObject,eventdata, handles,fname);
    end
end



% --- Executes when selected object is changed in Constraint_Type.
function Constraint_Type_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in Constraint_Type
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  global constraint_type;
  global click;
  if get(handles.soft,'Value') == 1
    handles.constraint_type = "Soft";
  else
    handles.constraint_type = "Hard";
  end
  constraint_type = handles.constraint_type;
  click = 0;
  guidata(hObject, handles);
end


% --- Executes on button press in closed.
function closed_Callback(hObject, eventdata, handles)
% hObject    handle to closed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of closed
end


% --- Executes when selected object is changed in constraintSnake.
function constraintSnake_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in constraintSnake 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    global snake_constraint;
    global click;
    if get(handles.snake1, 'Value') == 1
        handles.snake_constraint = 1;
    else
        handles.snake_constraint = 2;
    end
    
    disp(handles.snake_constraint);
  
    snake_constraint = handles.snake_constraint;
    click = 0;
    guidata(hObject, handles);
end
