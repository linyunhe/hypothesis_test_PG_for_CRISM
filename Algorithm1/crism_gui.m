function varargout = crism_gui(varargin)
%CRISM_GUI M-file for crism_gui.fig
%      CRISM_GUI, by itself, creates a new CRISM_GUI or raises the existing
%      singleton*.
%
%      H = CRISM_GUI returns the handle to a new CRISM_GUI or the handle to
%      the existing singleton*.
%
%      CRISM_GUI('Property','Value',...) creates a new CRISM_GUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to crism_gui_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      CRISM_GUI('CALLBACK') and CRISM_GUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in CRISM_GUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help crism_gui

% Last Modified by GUIDE v2.5 07-Apr-2017 11:27:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @crism_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @crism_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before crism_gui is made visible.
function crism_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for crism_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes crism_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = crism_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in start_button.
function start_button_Callback(hObject, eventdata, handles)
% hObject    handle to start_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_react ( hObject, handles )

% --- Executes on button press in load_button.
function load_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_react ( hObject, handles )

% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_react ( hObject, handles )


function radius_name_Callback(hObject, eventdata, handles)
% hObject    handle to radius_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of radius_name as text
%        str2double(get(hObject,'String')) returns contents of radius_name as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function radius_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radius_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radius_browse.
function radius_browse_Callback(hObject, eventdata, handles)
% hObject    handle to radius_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_react ( hObject, handles )


function spice_name_Callback(hObject, eventdata, handles)
% hObject    handle to spice_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spice_name as text
%        str2double(get(hObject,'String')) returns contents of spice_name as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function spice_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spice_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in spice_browse.
function spice_browse_Callback(hObject, eventdata, handles)
% hObject    handle to spice_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_react ( hObject, handles )


function topo_name_Callback(hObject, eventdata, handles)
% hObject    handle to topo_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of topo_name as text
%        str2double(get(hObject,'String')) returns contents of topo_name as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function topo_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to topo_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in topo_browse.
function topo_browse_Callback(hObject, eventdata, handles)
% hObject    handle to topo_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_react ( hObject, handles )

% --- Executes on button press in spat_pen_flag.
function spat_pen_flag_Callback(hObject, eventdata, handles)
% hObject    handle to spat_pen_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of spat_pen_flag
gui_react ( hObject, handles )

% --- Executes on button press in spec_pen_flag.
function spec_pen_flag_Callback(hObject, eventdata, handles)
% hObject    handle to spec_pen_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of spec_pen_flag
gui_react ( hObject, handles )


function beta_spat_Callback(hObject, eventdata, handles)
% hObject    handle to beta_spat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of beta_spat as text
%        str2double(get(hObject,'String')) returns contents of beta_spat as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function beta_spat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beta_spat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function beta_spec_Callback(hObject, eventdata, handles)
% hObject    handle to beta_spec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of beta_spec as text
%        str2double(get(hObject,'String')) returns contents of beta_spec as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function beta_spec_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beta_spec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function delta_spat_Callback(hObject, eventdata, handles)
% hObject    handle to delta_spat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delta_spat as text
%        str2double(get(hObject,'String')) returns contents of delta_spat as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function delta_spat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delta_spat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function delta_spec_Callback(hObject, eventdata, handles)
% hObject    handle to delta_spec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delta_spec as text
%        str2double(get(hObject,'String')) returns contents of delta_spec as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function delta_spec_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delta_spec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function spec_span_Callback(hObject, eventdata, handles)
% hObject    handle to spec_span (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spec_span as text
%        str2double(get(hObject,'String')) returns contents of spec_span as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function spec_span_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spec_span (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function spec_range_Callback(hObject, eventdata, handles)
% hObject    handle to spec_span (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spec_span as text
%        str2double(get(hObject,'String')) returns contents of spec_span as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function spec_range_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spec_span (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pen_internal_iters_Callback(hObject, eventdata, handles)
% hObject    handle to pen_internal_iters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pen_internal_iters as text
%        str2double(get(hObject,'String')) returns contents of pen_internal_iters as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function pen_internal_iters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pen_internal_iters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function out_file_prefix_Callback(hObject, eventdata, handles)
% hObject    handle to out_file_prefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of out_file_prefix as text
%        str2double(get(hObject,'String')) returns contents of out_file_prefix as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function out_file_prefix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to out_file_prefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_c_flag.
function save_c_flag_Callback(hObject, eventdata, handles)
% hObject    handle to save_c_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_c_flag
gui_react ( hObject, handles )

% --- Executes on button press in save_mu_flag.
function save_mu_flag_Callback(hObject, eventdata, handles)
% hObject    handle to save_mu_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_mu_flag
gui_react ( hObject, handles )

function save_proj_flag_Callback(hObject, eventdata, handles)
% hObject    handle to save_mu_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_mu_flag
gui_react ( hObject, handles )

function save_input_flag_Callback(hObject, eventdata, handles)
% hObject    handle to save_mu_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_mu_flag
gui_react ( hObject, handles )

function save_spectrogram_flag_Callback(hObject, eventdata, handles)
% hObject    handle to save_mu_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_mu_flag
gui_react ( hObject, handles )

function save_iters_Callback(hObject, eventdata, handles)
% hObject    handle to save_iters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of save_iters as text
%        str2double(get(hObject,'String')) returns contents of save_iters as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function save_iters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to save_iters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_iters_Callback(hObject, eventdata, handles)
% hObject    handle to num_iters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_iters as text
%        str2double(get(hObject,'String')) returns contents of num_iters as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function num_iters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_iters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in switchover_flag.
function switchover_flag_Callback(hObject, eventdata, handles)
% hObject    handle to switchover_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of switchover_flag
gui_react ( hObject, handles )


function switchover_at_iter_Callback(hObject, eventdata, handles)
% hObject    handle to switchover_at_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of switchover_at_iter as text
%        str2double(get(hObject,'String')) returns contents of switchover_at_iter as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function switchover_at_iter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to switchover_at_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function switchover_new_vals_Callback(hObject, eventdata, handles)
% hObject    handle to switchover_new_vals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of switchover_new_vals as text
%        str2double(get(hObject,'String')) returns contents of switchover_new_vals as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function switchover_new_vals_CreateFcn(hObject, eventdata, handles)
% hObject    handle to switchover_new_vals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kernel_size_Callback(hObject, eventdata, handles)
% hObject    handle to kernel_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kernel_size as text
%        str2double(get(hObject,'String')) returns contents of kernel_size as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function kernel_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kernel_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in projected_start_flag.
function projected_start_flag_Callback(hObject, eventdata, handles)
% hObject    handle to projected_start_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of projected_start_flag
gui_react ( hObject, handles )


function proj_scene_name_Callback(hObject, eventdata, handles)
% hObject    handle to proj_scene_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of proj_scene_name as text
%        str2double(get(hObject,'String')) returns contents of proj_scene_name as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function proj_scene_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to proj_scene_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in proj_scene_browse.
function proj_scene_browse_Callback(hObject, eventdata, handles)
% hObject    handle to proj_scene_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_react ( hObject, handles )


function init_ssa_guess_Callback(hObject, eventdata, handles)
% hObject    handle to init_ssa_guess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of init_ssa_guess as text
%        str2double(get(hObject,'String')) returns contents of init_ssa_guess as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function init_ssa_guess_CreateFcn(hObject, eventdata, handles)
% hObject    handle to init_ssa_guess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pix_spacing_Callback(hObject, eventdata, handles)
% hObject    handle to pix_spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pix_spacing as text
%        str2double(get(hObject,'String')) returns contents of pix_spacing as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function pix_spacing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pix_spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ddr_name_Callback(hObject, eventdata, handles)
% hObject    handle to ddr_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ddr_name as text
%        str2double(get(hObject,'String')) returns contents of ddr_name as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function ddr_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ddr_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sb_name_Callback(hObject, eventdata, handles)
% hObject    handle to sb_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sb_name as text
%        str2double(get(hObject,'String')) returns contents of sb_name as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function sb_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sb_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wa_name_Callback(hObject, eventdata, handles)
% hObject    handle to wa_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wa_name as text
%        str2double(get(hObject,'String')) returns contents of wa_name as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function wa_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wa_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ssa_browse.
function ssa_browse_Callback(hObject, eventdata, handles)
% hObject    handle to ssa_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_react ( hObject, handles )

% --- Executes on button press in ddr_browse.
function ddr_browse_Callback(hObject, eventdata, handles)
% hObject    handle to ddr_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_react ( hObject, handles )

% --- Executes on button press in sb_browse.
function sb_browse_Callback(hObject, eventdata, handles)
% hObject    handle to sb_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_react ( hObject, handles )

% --- Executes on button press in wa_browse.
function wa_browse_Callback(hObject, eventdata, handles)
% hObject    handle to wa_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_react ( hObject, handles )


function ssa_row_min_Callback(hObject, eventdata, handles)
% hObject    handle to ssa_row_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ssa_row_min as text
%        str2double(get(hObject,'String')) returns contents of ssa_row_min as a double


% --- Executes during object creation, after setting all properties.
function ssa_row_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ssa_row_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ssa_row_max_Callback(hObject, eventdata, handles)
% hObject    handle to ssa_row_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ssa_row_max as text
%        str2double(get(hObject,'String')) returns contents of ssa_row_max as a double


% --- Executes during object creation, after setting all properties.
function ssa_row_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ssa_row_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ssa_row_all.
function ssa_row_all_Callback(hObject, eventdata, handles)
% hObject    handle to ssa_row_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ssa_row_all
gui_react ( hObject, handles )


function ssa_col_min_Callback(hObject, eventdata, handles)
% hObject    handle to ssa_col_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ssa_col_min as text
%        str2double(get(hObject,'String')) returns contents of ssa_col_min as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function ssa_col_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ssa_col_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ssa_col_max_Callback(hObject, eventdata, handles)
% hObject    handle to ssa_col_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ssa_col_max as text
%        str2double(get(hObject,'String')) returns contents of ssa_col_max as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function ssa_col_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ssa_col_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ssa_col_all.
function ssa_col_all_Callback(hObject, eventdata, handles)
% hObject    handle to ssa_col_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ssa_col_all
gui_react ( hObject, handles )


function ddr_col_min_Callback(hObject, eventdata, handles)
% hObject    handle to ddr_col_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ddr_col_min as text
%        str2double(get(hObject,'String')) returns contents of ddr_col_min as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function ddr_col_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ddr_col_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ddr_col_max_Callback(hObject, eventdata, handles)
% hObject    handle to ddr_col_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ddr_col_max as text
%        str2double(get(hObject,'String')) returns contents of ddr_col_max as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function ddr_col_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ddr_col_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ddr_col_all.
function ddr_col_all_Callback(hObject, eventdata, handles)
% hObject    handle to ddr_col_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ddr_col_all
gui_react ( hObject, handles )


function ddr_row_min_Callback(hObject, eventdata, handles)
% hObject    handle to ddr_row_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ddr_row_min as text
%        str2double(get(hObject,'String')) returns contents of ddr_row_min as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function ddr_row_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ddr_row_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ddr_row_max_Callback(hObject, eventdata, handles)
% hObject    handle to ddr_row_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ddr_row_max as text
%        str2double(get(hObject,'String')) returns contents of ddr_row_max as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function ddr_row_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ddr_row_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ddr_row_all.
function ddr_row_all_Callback(hObject, eventdata, handles)
% hObject    handle to ddr_row_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ddr_row_all
gui_react ( hObject, handles )


function sb_col_min_Callback(hObject, eventdata, handles)
% hObject    handle to sb_col_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sb_col_min as text
%        str2double(get(hObject,'String')) returns contents of sb_col_min as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function sb_col_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sb_col_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sb_col_max_Callback(hObject, eventdata, handles)
% hObject    handle to sb_col_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sb_col_max as text
%        str2double(get(hObject,'String')) returns contents of sb_col_max as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function sb_col_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sb_col_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in sb_col_all.
function sb_col_all_Callback(hObject, eventdata, handles)
% hObject    handle to sb_col_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sb_col_all
gui_react ( hObject, handles )


function wa_col_min_Callback(hObject, eventdata, handles)
% hObject    handle to wa_col_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wa_col_min as text
%        str2double(get(hObject,'String')) returns contents of wa_col_min as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function wa_col_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wa_col_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wa_col_max_Callback(hObject, eventdata, handles)
% hObject    handle to wa_col_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wa_col_max as text
%        str2double(get(hObject,'String')) returns contents of wa_col_max as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function wa_col_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wa_col_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in wa_col_all.
function wa_col_all_Callback(hObject, eventdata, handles)
% hObject    handle to wa_col_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of wa_col_all
gui_react ( hObject, handles )


function ssa_band_min_Callback(hObject, eventdata, handles)
% hObject    handle to ssa_band_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ssa_band_min as text
%        str2double(get(hObject,'String')) returns contents of ssa_band_min as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function ssa_band_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ssa_band_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ssa_band_max_Callback(hObject, eventdata, handles)
% hObject    handle to ssa_band_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ssa_band_max as text
%        str2double(get(hObject,'String')) returns contents of ssa_band_max as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function ssa_band_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ssa_band_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ssa_band_all.
function ssa_band_all_Callback(hObject, eventdata, handles)
% hObject    handle to ssa_band_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ssa_band_all
gui_react ( hObject, handles )


function sb_band_min_Callback(hObject, eventdata, handles)
% hObject    handle to sb_band_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sb_band_min as text
%        str2double(get(hObject,'String')) returns contents of sb_band_min as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function sb_band_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sb_band_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sb_band_max_Callback(hObject, eventdata, handles)
% hObject    handle to sb_band_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sb_band_max as text
%        str2double(get(hObject,'String')) returns contents of sb_band_max as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function sb_band_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sb_band_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in sb_band_all.
function sb_band_all_Callback(hObject, eventdata, handles)
% hObject    handle to sb_band_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sb_band_all
gui_react ( hObject, handles )


function wa_band_min_Callback(hObject, eventdata, handles)
% hObject    handle to wa_band_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wa_band_min as text
%        str2double(get(hObject,'String')) returns contents of wa_band_min as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function wa_band_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wa_band_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wa_band_max_Callback(hObject, eventdata, handles)
% hObject    handle to wa_band_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wa_band_max as text
%        str2double(get(hObject,'String')) returns contents of wa_band_max as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function wa_band_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wa_band_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in wa_band_all.
function wa_band_all_Callback(hObject, eventdata, handles)
% hObject    handle to wa_band_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of wa_band_all
gui_react ( hObject, handles )


% --- Executes on button press in nadir_geo_button.
function nadir_geo_button_Callback(hObject, eventdata, handles)
% hObject    handle to nadir_geo_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nadir_geo_button
gui_react ( hObject, handles )

% --- Executes on button press in off_nadir_geo_button.
function off_nadir_flat_geo_button_Callback(hObject, eventdata, handles)
% hObject    handle to off_nadir_geo_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of off_nadir_geo_button
gui_react ( hObject, handles )

% --- Executes on button press in off_nadir_geo_button.
function off_nadir_rolling_geo_button_Callback(hObject, eventdata, handles)
% hObject    handle to off_nadir_geo_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of off_nadir_geo_button
gui_react ( hObject, handles )


% --- Executes on button press in g1_checkbox.
function g1_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to g1_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of g1_checkbox



function output_dir_Callback(hObject, eventdata, handles)
% hObject    handle to output_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of output_dir as text
%        str2double(get(hObject,'String')) returns contents of output_dir as a double


% --- Executes during object creation, after setting all properties.
function output_dir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in output_dir_browse.
function output_dir_browse_Callback(hObject, eventdata, handles)
% hObject    handle to output_dir_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_react ( hObject, handles )


% --- Executes on button press in l_band_button.
function l_band_button_Callback(hObject, eventdata, handles)
% hObject    handle to l_band_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of l_band_button
gui_react ( hObject, handles )


% --- Executes on button press in s_band_button.
function s_band_button_Callback(hObject, eventdata, handles)
% hObject    handle to s_band_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of s_band_button
gui_react ( hObject, handles )



function spec_offset_Callback(hObject, eventdata, handles)
% hObject    handle to spec_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spec_offset as text
%        str2double(get(hObject,'String')) returns contents of spec_offset as a double


% --- Executes during object creation, after setting all properties.
function spec_offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spec_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in no_switch_radio.
function no_switch_radio_Callback(hObject, eventdata, handles)
% hObject    handle to no_switch_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of no_switch_radio
gui_react ( hObject, handles )

% --- Executes on button press in switch_penalty_radio.
function switch_penalty_radio_Callback(hObject, eventdata, handles)
% hObject    handle to switch_penalty_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of switch_penalty_radio
gui_react ( hObject, handles )


% --- Executes on button press in switch_penalty_spatial.
function switch_penalty_spatial_Callback(hObject, eventdata, handles)
% hObject    handle to switch_penalty_spatial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of switch_penalty_spatial
gui_react ( hObject, handles )

% --- Executes on button press in switch_penalty_spectral.
function switch_penalty_spectral_Callback(hObject, eventdata, handles)
% hObject    handle to switch_penalty_spectral (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of switch_penalty_spectral
gui_react ( hObject, handles )


% --- Executes on button press in output_dir_auto_button.
function output_dir_auto_button_Callback(hObject, eventdata, handles)
% hObject    handle to output_dir_auto_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_react ( hObject, handles )


% --- Executes on button press in auto_find_sb_wa_button.
function auto_find_sb_wa_button_Callback(hObject, eventdata, handles)
% hObject    handle to auto_find_sb_wa_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spec_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function version_num_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to version_num_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
try
    [ver,algo] = get_version_number();
catch e
    ver = '<VER>';
    algo = '<ALGO>';
end
set(hObject,'String',['CRISM ', algo, ' GUI version ',ver]);




function pix_spacing_output_Callback(hObject, eventdata, handles)
% hObject    handle to pix_spacing_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pix_spacing_output as text
%        str2double(get(hObject,'String')) returns contents of pix_spacing_output as a double
gui_react ( hObject, handles )

% --- Executes during object creation, after setting all properties.
function pix_spacing_output_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pix_spacing_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pre_end_Callback(hObject, eventdata, handles)
% hObject    handle to load_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_react ( hObject, handles )

function post_start_Callback(hObject, eventdata, handles)
% hObject    handle to load_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_react ( hObject, handles )

function pre_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pix_spacing_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function post_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pix_spacing_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
