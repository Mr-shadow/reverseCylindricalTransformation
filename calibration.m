picPath = 'wangGe.jpg';
pic=imread(picPath);
imshow(pic);
try
    load 'calibration.mat';
catch 
    A=zeros(10,13);
    B=zeros(10,13);
end
a=ginput;
row=input('\nµÚ¼¸ÁÐ?');
A(:,row)=a(:,1);
B(:,row)=a(:,2);
save 'calibration.mat' A B;




