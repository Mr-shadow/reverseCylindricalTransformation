
clear;
%%% srcpath: 需要进行曲面反变换到平面的图像
%%% Pup.mat: 使用getInput_my 生成的文件，里面是对应曲面的网格纸使用鼠标点击取的坐标。
srcpath='ceshiTu.jpg'; 
load  'calibration.mat'; 

img = imread(srcpath);
[length,width,d] = size(img); 
if d==3
    img=rgb2gray(img);
end
imgDouble = double(img);
imagesc(imgDouble);

Px=B;
Py=A;
mZeta=size(Px,1);
mXi=size(Px,2);
PZeta=zeros(length,width);
PXi=zeros(length,width);

%%% choose ROI 
mini=1;
maxi=length;
minj=1;
maxj=width;

for i=mini:maxi
    for j=minj:maxj             
        [s,JI]=min(((Px-i).*(Px-i)+(Py-j).*(Py-j)));
        [~,Xi]=min(s);
        Zeta=JI(Xi);
        [PZeta(i,j),PXi(i,j),imgDouble(i,j)]=Pscale(i,j,Zeta,Xi,Px,Py,mZeta,mXi,imgDouble); 
    end
 end
PZeta=PZeta(mini:maxi,minj:maxj);
PXi=PXi(mini:maxi,minj:maxj);
imgDouble=imgDouble(mini:maxi,minj:maxj);

delta=0.025;
a=fix(min(min(PZeta)));b=fix(max(max(PZeta)));
c=fix(min(min(PXi)));d=fix(max(max(PXi)));

[X, Y] = meshgrid(a:delta:b, c:delta:d);
dstImg = griddata(PZeta',PXi',imgDouble', X, Y);
dstImg=dstImg';
dstImg=dstImg((1-a)/delta+1:(mZeta-a)/delta+1,(1-c)/delta+1:(mXi-c)/delta+1);
dstImg = int8(dstImg);
imshow(dstImg);
 imwrite()

function [PZeta,PXi,H]=Pscale(i,j,Zeta,Xi,Px,Py,mZeta,mXi,Hpic)
    %Test the 4 nearby elements
    PZeta_right=0;
    PXi_right=0;
    for m=-1:0
        for n=-1:0
            if(((Zeta+m>0)&&(Zeta+m<mZeta))&&((Xi+n>0)&&(Xi+n<mXi)))    %for boundarys
                [Zeta_Local,Xi_Local]=Plocal(i,j,Zeta+m,Xi+n,Px,Py);
                if((abs(Zeta_Local)<=0.5)&&(abs(Xi_Local)<=0.5))
                    PZeta_right=Zeta+m+0.5+Zeta_Local;
                    PXi_right=Xi+n+0.5+Xi_Local;
                end
                    PZeta=Zeta+m+0.5+Zeta_Local;
                    PXi=Xi+n+0.5+Xi_Local;  
                    PZeta=max(-2,PZeta);PZeta=min(mZeta+2,PZeta);
                    PXi=max(-2,PXi);PXi=min(mXi+2,PXi);
                    if(PZeta>=1&&PZeta<=mZeta&&PXi>=1&&PXi<=mXi)
                        PZeta=-2;
                        PXi=-2;
                    end
                
            end            
        end
    end    
    if((PZeta_right~=0)&&(PXi_right~=0))
        PZeta=PZeta_right;
        PXi=PXi_right;
        H=Hpic(i,j);
    else
        H=Hpic(i,j);
    end
end

function [Zeta_Local,Xi_Local]=Plocal(x,y,Zeta,Xi,Px,Py)
    %Corner data
    Zeta_Corner=[-1,1,1,-1];
    Xi_Corner  =[-1,-1,1,1];   
    x_Corner=[Px(Zeta,Xi),Px(Zeta+1,Xi),Px(Zeta+1,Xi+1),Px(Zeta,Xi+1)];
    y_Corner=[Py(Zeta,Xi),Py(Zeta+1,Xi),Py(Zeta+1,Xi+1),Py(Zeta,Xi+1)];    
    %Iteration    
    ZX_Local=[0;0];    
    dZX=[0;0];
    err=1;
    k=0;
    while(err>1e-10&&k<50)         
        %new Zeta, Fi & J
        ZX_Local=ZX_Local+dZX;         
        Fi=(0.5+ZX_Local(1)*Zeta_Corner).*(0.5+ZX_Local(2)*Xi_Corner);
        J=[x_Corner;y_Corner]*[Zeta_Corner.*(0.5+ZX_Local(2)*Xi_Corner);Xi_Corner.*(0.5+ZX_Local(1)*Zeta_Corner)]'; %care!
        %xy=ZX*Fi
        xy=[x_Corner;y_Corner]*Fi';
        dxy=[x;y]-xy;
        %dZX=[J-]*dx;  
        dZX=J\dxy; %for a  matrix smaller than 100*100, 'inv' is much faster than '\'
        %error
        err=norm(dZX,2);
        k=k+1;
    end 
    Zeta_Local=ZX_Local(1);
    Xi_Local=ZX_Local(2);
end
    
    