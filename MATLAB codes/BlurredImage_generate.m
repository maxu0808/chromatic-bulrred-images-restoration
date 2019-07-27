%%% This script is to generate chromatic blurred images based on
%%% spectral images and PSFs.
clear all;
close all;
clc;

lambda=400:5:700;
%% Load PSFs
PSF=cell(61,61);
path='.\All_PSF';
for i=1:61
    psf=load([path '\' num2str(lambda(i)) 'nm.mat']);
    psf=psf.psf;
    for j=1:61
        psf(:,:,j)=psf(:,:,j)./sum(sum(psf(:,:,j)));
        PSF{i,j}=psf(:,:,j);
    end
end
%% Load spectral images and generate chromatic blurred images
spectralimage=load('spectralimage.mat');
spectralimage=spectralimage.spectralimage;
[testrows,testcols,~]=size(spectralimage);
[psf_rows,psf_cols,~]=size(psf);
fft_rows=testrows+psf_rows-1;
fft_cols=testcols+psf_cols-1;
for i=1:61
    for j=1:61
        PSF_fft{i,j}=fftshift(fft2(PSF{i,j},fft_rows,fft_cols));
    end
end
for i=1:61
    spectralimage_fft(:,:,i)=fftshift(fft2(spectralimage(:,:,i),fft_rows,fft_cols));
end
blurredimage_fft=zeros(fft_rows,fft_cols,61);
for i=1:61 
    for j=1:61  
    blurredimage_fft(:,:,i)=blurredimage_fft(:,:,i)+spectralimage_fft(:,:,j).* PSF_fft{j,i};
    end
end
for i=1:61
    blurredimage(:,:,i)=abs(ifft2(ifftshift(blurredimage_fft(:,:,i))));
end
% save('blurredimage.mat','blurredimage');
value_max=max(max(max(blurredimage)));
blurredimage=blurredimage./value_max;
for i=1:61
    imshow(blurredimage(:,:,i));
    shg
    title(['The retinal position is at' num2str(lambda(i)) 'nm focused'])
    pause(0.5);
end
blurredimage=blurredimage.*value_max;
save('blurredimage.mat','blurredimage');