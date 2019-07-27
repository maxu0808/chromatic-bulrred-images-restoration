close all;clear all;clc;
spectralimage=load('spectralimage.mat');
spectralimage=spectralimage.spectralimage;
restoredimage=load('restoredimage.mat');
restoredimage=restoredimage.restoredimage_cut;
[testrows,testcols,num]=size(spectralimage);
%% compute MSE
MSE=0;
for i=1:num
    i
   for j=1:testrows
       for k=1:testcols
       MSE=MSE+(255*restoredimage(j,k,i)-255*spectralimage(j,k,i))^2/(testrows*testcols*num);
       end
   end
end
%% compute PSNR
PSNR=10*log10(255^2/MSE);
