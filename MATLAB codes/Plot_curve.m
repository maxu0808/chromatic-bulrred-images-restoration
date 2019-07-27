%%% This MATLAB script is to compare spectral curves between spectral
%%% images and restored images.
clear all;
close all;
clc;
lambda=400:5:700;
spectralimage=load('spectralimage.mat');
spectralimage=spectralimage.spectralimage;
restoredimage=load('restoredimage.mat');
restoredimage=restoredimage.restoredimage_cut;
% picking two position randomly
green_position=[43,84];
blue_position=[15,5];
% set region size is 3*3
region_size=3;

%% green region
spec=zeros(1,61);
res=zeros(1,61);
for i=1:61
    spec(i)=mean2(spectralimage(green_position(1):green_position(1)+region_size,...
        green_position(2):green_position(2)+region_size,i));
    res(i)=mean2(restoredimage(green_position(1):green_position(1)+region_size,...
        green_position(2):green_position(2)+region_size,i));
end
figure(1)
plot(lambda,spec,'r-*')
hold on;
plot(lambda,res,'b-o')
grid on;
title('Green Region')
legend('ground truth','restored')

%% blue region
spec=zeros(1,61);
res=zeros(1,61);
for i=1:61
    spec(i)=mean2(spectralimage(blue_position(1):blue_position(1)+region_size,...
        blue_position(2):blue_position(2)+region_size,i));
    res(i)=mean2(restoredimage(blue_position(1):blue_position(1)+region_size,...
        blue_position(2):blue_position(2)+region_size,i));
end
figure(2)
plot(lambda,spec,'r-*')
hold on;
plot(lambda,res,'b-o')
grid on;
title('Blue Region')
legend('ground truth','restored')