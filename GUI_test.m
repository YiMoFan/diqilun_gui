function varargout = GUI_test(varargin)
% GUI_TEST MATLAB code for GUI_test.fig
%      GUI_TEST, by itself, creates a new GUI_TEST or raises the existing
%      singleton*.
%
%      H = GUI_TEST returns the handle to a new GUI_TEST or the handle to
%      the existing singleton*.
%
%      GUI_TEST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_TEST.M with the given input arguments.
%
%      GUI_TEST('Property','Value',...) creates a new GUI_TEST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_test_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_test_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_test

% Last Modified by GUIDE v2.5 21-Aug-2018 15:29:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_test_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_test_OutputFcn, ...
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


% --- Executes just before GUI_test is made visible.
function GUI_test_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_test (see VARARGIN)

% Choose default command line output for GUI_test
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_test wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% 设置背景图案
backgroundlmage2 = importdata('4.jpg');
axes(handles.axes1);
image(backgroundlmage2);
axis off
% 设置显示窗口图案
backgroundlmage = importdata('C:\Users\jj\Desktop\QQ图片20180819165938.jpg');
axes(handles.axes2);
image(backgroundlmage);
axis off

axes(handles.axes3);
image(backgroundlmage);
axis off

% 设置按钮图案
% anniu=imread('C:\Users\jj\Desktop\3.jpg');
% set(handles.pushbutton1,'cdata',anniu);
% set(handles.pushbutton2,'cdata',anniu);
% set(handles.pushbutton3,'cdata',anniu);
% set(handles.pushbutton4,'cdata',anniu);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_test_OutputFcn(hObject, eventdata, handles) 
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


global b filename2 pathname2 flag
% 选择图片
[filename2, pathname2] = uigetfile({'*.fts'; '*.jpg'; '*.bmp'; '*.gif'; '*.png'}, '选择图片');
% 合成路径+文件名
name = [pathname2 filename2];
flag = 0;
%% 读取图片
a=fitsinfo(name);
% 得到一个1*1的结构变量a，里面记载了各种fts文件的总体信息
b=fitsread(name);
% 得到fts文件的数据矩阵
% 使用第一个axes
axes(handles.axes2);
% 显示图片
imagesc(b,[min(min(b)),max(max(b))]);
% 以min和max为界画出灰度图像
colormap(gray);
%色调
axis equal
axis off % 去掉坐标轴

axes(handles.axes3);
% 显示图片
imagesc(b,[min(min(b)),max(max(b))]);
% 以min和max为界画出灰度图像
colormap(gray);
%色调
axis equal
axis off % 去掉坐标轴

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
h = getframe(handles.axes2);
imwrite(h.cdata,fullfile(PathName,FileName));
% name2 = [filename(1:end-4),'.jpg'];
% imwrite(h.cdata,['D:\第七轮\第七轮数据\图片\',name2]);
hs = msgbox('图片保存成功！','提示信息','help','modal');
ht = findobj(hs,'Type','text');
set(ht,'FontSize',15,'Unit','normal');
% 改变对话框大小
set(hs,'Resize','on');% 自动改变
% set(hs,'Position',[400 200 200 200])


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile({'*.fts'; '*.jpg'; '*.bmp'; '*.gif'; '*.png'}, '选择图片');
listing = dir([pathname,'*.fts']);
[m,n]=size(listing);
for i=1:m
    name=listing(i).name;
    a=fitsinfo(name);
    %得到一个1*1的结构变量a，里面记载了各种fts文件的总体信息
    b=fitsread(name);
    %得到fts文件的数据矩阵
    axes(handles.axes3);
    imagesc(b,[min(min(b)),max(max(b))]);
    % 以min和max为界画出灰度图像
    colormap(gray);
    %色调
    axis equal 
    %坐标轴比例
    % pause(1);
    % 把 fts 格式  转换为  jpg  格式
    name2 = [name(1:end-4),'.jpg'];
    h = getframe(handles.axes2);
    imwrite(h.cdata,[pathname,name2]);
end



% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(gcf);

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
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
% imwrite(h.cdata,['D:\第七轮\第七轮数据\图片\',name2]);
hs = msgbox('图片保存成功！','提示信息','help','modal');
ht = findobj(hs,'Type','text');
set(ht,'FontSize',15,'Unit','normal');
% 改变对话框大小
set(hs,'Resize','on');% 自动改变
% set(hs,'Position',[400 200 200 200])


% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup1 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global b
str = get(hObject,'string');
axes(handles.axes2);

switch str
    case '原图'
        imagesc(b,[min(min(b)),max(max(b))]);
        colormap(gray);
        axis equal
        axis off % 去掉坐标轴
    case '索引图'
        bw = mat2gray(b);
        mm = grayslice(bw,16);
        imshow(mm,jet(16));
end


% --- Executes when selected object is changed in uibuttongroup2.
function uibuttongroup2_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup2 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global b filename2 pathname2 flag name
str2 = get(hObject,'string');
axes(handles.axes3);


switch str2
    case '灰度图'
        imagesc(b,[min(min(b)),max(max(b))]);
        colormap(gray);
        axis equal
        axis off % 去掉坐标轴
    case '同态滤波'
        if flag == 0
            name = [filename2(1:end-4),'.bmp'];
            h = getframe(handles.axes2);
            imwrite(h.cdata,[pathname2,name]);
        end
        flag = flag + 1;
       %% 同台滤波测试函数
        img=imread(name);
        I = img; 
        % I=double(rgb2gray(I)); 
        I=double(I);
        [M,N,channals]=size(I);    
        rL=0.5;    
        rH=5;%可根据需要效果调整参数    
        c=3;    
        d0=9;  
        I3=I;
        for ch=1:channals
            I1=log(I(:,:,ch)+1);%取对数    
            FI=fft2(I1);%傅里叶变换    
            n1=floor(M/2);    
            n2=floor(N/2);    
            for i=1:M    
                for j=1:N    
                    D(i,j)=((i-n1).^2+(j-n2).^2);    
                    H(i,j)=(rH-rL).*(exp(c*(-D(i,j)./(d0^2))))+rL;%高斯同态滤波    
                end    
            end    
            I2=ifft2(H.*FI);%傅里叶逆变换    
            I3(:,:,ch)=real(exp(I2));
        end
        minV=min(min(min(I3)));
        maxV=max(max(max(I3)));
        for ch=1:channals
            for i=1:M   
                for j=1:N    
                    I3(i,j,ch)=255* (I3(i,j,ch)-minV)./(maxV-minV);
                end
            end
        end
        imshow(uint8(I3))
    case '参数同态'
        if flag == 0
            name = [filename2(1:end-4),'.bmp'];
            h = getframe(handles.axes2);
            imwrite(h.cdata,[pathname2,name]);
        end
        flag = flag + 1;
        prompt = {'参数a取值：','参数b的取值'};%设置提示字符串
        title = '同态滤波参数确定';%设置标题
        numlines = 1;%指定输入数据的行数
        defAns = {'0.5','1.2'};%设定默认值
        Resize = 'on';%设定对话框尺寸可调节
        answer = inputdlg(prompt,title,numlines,defAns,'on');%创建输入对话框
        aa = answer(1,1);
        bb = answer(2,1);
        aa = str2num(char(aa));
        bb = str2num(char(bb));
        
        img = imread(name);
        img = rgb2gray(img);
        [height, width] = size(img);
        
        f = double(img);
        % f = Centralize(f);
        [height, width] = size(f);
        for i = 1 : height
            for j = 1 : width
                if mod(i + j, 2) == 1
                    f(i, j) = -f(i, j);

                end
            end
        end
        F = fft2(f);
        D0 = 1;
        % H = BHPF(D0, height, width);

        % aa = 0.5; bb = 1.2;

        for i = 1 : height
            x = i - (height / 2);
            for j = 1 : width
                y = j - (width / 2);
                H(i, j) = (1 / (1 + (D0 ^ 2) / (x ^ 2 + y ^ 2)))*bb+aa;
            end
        end

        g = real(ifft2(H .* F));

        % g = Centralize(g);

        [height, width] = size(g);
        for i = 1 : height
            for j = 1 : width
                if mod(i + j, 2) == 1
                    g(i, j) = -g(i, j);

                end
            end
        end
        imshow(uint8(g));
    case '优化算法'
        if flag == 0
            name = [filename2(1:end-4),'.bmp'];
            h = getframe(handles.axes2);
            imwrite(h.cdata,[pathname2,name]);
        end
        flag = flag + 1;
        
        I = imread(name);

        R = I(:, :, 1);
        [N1, M1] = size(R);
        R0 = double(R);
        Rlog = log(R0+1);
        Rfft2 = fft2(R0);

        sigma = 250;
        F = fspecial('gaussian', [N1,M1], sigma);
        Efft = fft2(double(F));

        DR0 = Rfft2.* Efft;
        DR = ifft2(DR0);

        DRlog = log(DR +1);
        Rr = Rlog - DRlog;
        EXPRr = exp(Rr);
        MIN = min(min(EXPRr));
        MAX = max(max(EXPRr));
        EXPRr = (EXPRr - MIN)/(MAX - MIN);
        EXPRr = adapthisteq(EXPRr);

        G = I(:, :, 2);

        G0 = double(G);
        Glog = log(G0+1);
        Gfft2 = fft2(G0);

        DG0 = Gfft2.* Efft;
        DG = ifft2(DG0);
        
        DGlog = log(DG +1);
        Gg = Glog - DGlog;
        EXPGg = exp(Gg);
        MIN = min(min(EXPGg));
        MAX = max(max(EXPGg));
        EXPGg = (EXPGg - MIN)/(MAX - MIN);
        EXPGg = adapthisteq(EXPGg);

        B = I(:, :, 3);

        B0 = double(B); 
        Blog = log(B0+1);
        Bfft2 = fft2(B0);

        DB0 = Bfft2.* Efft;
        DB = ifft2(DB0);

        DBlog = log(DB+1);
        Bb = Blog - DBlog;
        EXPBb = exp(Bb);
        MIN = min(min(EXPBb));
        MAX = max(max(EXPBb));
        EXPBb = (EXPBb - MIN)/(MAX - MIN);
        EXPBb = adapthisteq(EXPBb);

        result = cat(3, EXPRr, EXPGg, EXPBb);
        imshow(result);
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global b  filename2 pathname2 flag name

if flag == 0
    name = [filename2(1:end-4),'.bmp'];
    h = getframe(handles.axes2);
    imwrite(h.cdata,[pathname2,name]);
end
flag = flag + 1;

I = imread(name);
% rng('default');
[m, n, p] = size(I);
X = reshape(double(I), m*n, p);
k = 7; b = 2;
iter = 0;
[N, p] = size(X);
P = randn(N, k);
P = P./(sum(P, 2)*ones(1, k));
J_prev = inf; J = [];
while true,
    iter = iter + 1;
    t = P.^b;
    C = (X'*t)'./(sum(t)'*ones(1, p));
    dist = sum(X.*X, 2)*ones(1, k) + (sum(C.*C, 2)*ones(1, N))'-2*X*C';
    t2 = (1./dist).^(1/(b-1));
    P = t2./(sum(t2, 2)*ones(1, k));
    J_cur = sum(sum((P.^b).*dist))/N;
    J = [J J_cur];
    if norm(J_cur-J_prev, 'fro') < 1e-3,
        break;
    end
%   display(sprintf('#iteration: %03d, objective function: %f', iter, J_cur));
   J_prev = J_cur; 
end
[~, label] = min(dist, [], 2);
axes(handles.axes3);
imshow(uint8(reshape(C(label, :), m, n, p)))