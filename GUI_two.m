function varargout = GUI_two(varargin)
% GUI_TWO MATLAB code for GUI_two.fig
%      GUI_TWO, by itself, creates a new GUI_TWO or raises the existing
%      singleton*.
%
%      H = GUI_TWO returns the handle to a new GUI_TWO or the handle to
%      the existing singleton*.
%
%      GUI_TWO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_TWO.M with the given input arguments.
%
%      GUI_TWO('Property','Value',...) creates a new GUI_TWO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_two_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_two_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_two

% Last Modified by GUIDE v2.5 20-Aug-2018 21:20:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_two_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_two_OutputFcn, ...
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


% --- Executes just before GUI_two is made visible.
function GUI_two_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_two (see VARARGIN)

% Choose default command line output for GUI_two
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_two wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% ���ñ���ͼ��
backgroundlmage2 = importdata('4.jpg');
axes(handles.axes1);
image(backgroundlmage2);
axis off
% ������ʾ����ͼ��
backgroundlmage = importdata('C:\Users\jj\Desktop\QQͼƬ20180819165938.jpg');
axes(handles.axes3);
image(backgroundlmage);
axis off

% ���ð�ťͼ��
% anniu=imread('C:\Users\jj\Desktop\3.jpg');
% set(handles.pushbutton1,'cdata',anniu);
% set(handles.pushbutton2,'cdata',anniu);
% set(handles.pushbutton3,'cdata',anniu);
% set(handles.pushbutton4,'cdata',anniu);

% --- Outputs from this function are returned to the command line.
function varargout = GUI_two_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global b
% ѡ��ͼƬ
[filename, pathname] = uigetfile({'*.fts'; '*.jpg'; '*.bmp'; '*.gif'; '*.png'}, 'ѡ��ͼƬ');
% �ϳ�·��+�ļ���
name = [pathname filename];
%% ��ȡͼƬ
a=fitsinfo(name);
% �õ�һ��1*1�Ľṹ����a����������˸���fts�ļ���������Ϣ
b=fitsread(name);
% �õ�fts�ļ������ݾ���
% ʹ�õ�һ��axes
axes(handles.axes3);
% ��ʾͼƬ
imagesc(b,[min(min(b)),max(max(b))]);
% ��min��maxΪ�续���Ҷ�ͼ��
colormap(gray);
%ɫ��
axis equal
axis off % ȥ��������

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName,PathName,filespec] = uiputfile({'*.jpg','JPEG(*.jpg)';...
                                 '*.bmp','Bitmap(*.bmp)';...
                                 '*.*',  'All Files (*.*)'},...
                                 'Save Picture','Untitled');
% if FileName==0
%   return;
% else
h = getframe(handles.axes3);
imwrite(h.cdata,fullfile(PathName,FileName));
% name2 = [filename(1:end-4),'.jpg'];
% imwrite(h.cdata,['D:\������\����������\ͼƬ\',name2]);
hs = msgbox('ͼƬ����ɹ���','��ʾ��Ϣ','help','modal');
ht = findobj(hs,'Type','text');
set(ht,'FontSize',15,'Unit','normal');
% �ı�Ի����С
set(hs,'Resize','on');% �Զ��ı�
% set(hs,'Position',[400 200 200 200])


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*.fts'; '*.jpg'; '*.bmp'; '*.gif'; '*.png'}, 'ѡ��ͼƬ');
listing = dir([pathname,'*.fts']);
[m,n]=size(listing);
for i=1:m
    name=listing(i).name;
    a=fitsinfo(name);
    %�õ�һ��1*1�Ľṹ����a����������˸���fts�ļ���������Ϣ
    b=fitsread(name);
    %�õ�fts�ļ������ݾ���
    axes(handles.axes3);
    imagesc(b,[min(min(b)),max(max(b))]);
    % ��min��maxΪ�续���Ҷ�ͼ��
    colormap(gray);
    %ɫ��
    axis equal 
    %���������
    % pause(1);
    % �� fts ��ʽ  ת��Ϊ  jpg  ��ʽ
    name2 = [name(1:end-4),'.jpg'];
    h = getframe(handles.axes3);
    imwrite(h.cdata,[pathname,name2]);
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(gcf);

% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup1 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global b
str = get(hObject,'string');
axes(handles.axes3);

switch str
    case 'ԭͼ'
        imagesc(b,[min(min(b)),max(max(b))]);
        colormap(gray);
        axis equal
        axis off % ȥ��������
    case '����ͼ'
        bw = mat2gray(b);
        mm = grayslice(bw,16);
        imshow(mm,jet(16));
end
