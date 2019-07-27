%%% This MATLAB script is to recover the spectral images using SVD
%%% algorithm
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
%% Load blurred images
blurredimage=load('blurredimage.mat');
blurredimage=blurredimage.blurredimage;
[fft_rows,fft_cols,~]=size(blurredimage);

%% Take the Fourier transform of the blurred images and PSFs
[psf_rows,psf_cols,~]=size(psf);
testrows=fft_rows-psf_rows+1;
testcols=fft_cols-psf_cols+1;
for i=1:61
    for j=1:61
        PSF_fft{i,j}=fftshift(fft2(PSF{i,j},fft_rows,fft_cols));
    end
end
for i=1:61
    blurredimage_fft(:,:,i)=fftshift(fft2(blurredimage(:,:,i),fft_rows,fft_cols));
end
%% SVD
psf_Matrix=zeros(61,61);
blurred_vector=zeros(1,61);
restoredimage_fft=zeros(fft_rows,fft_cols,61);
restoredimage=zeros(fft_rows,fft_cols,61);
alpha=1e-7;
for i=1:fft_rows
    rows=i
    for j=1:fft_cols
        for t=1:61
            for k=1:61  
        psf_Matrix(k,t)=PSF_fft{k,t}(i,j);
            end
        end
        for t=1:61
        blurred_vector(1,t)= blurredimage_fft(i,j,t);
        end
        [U S V]=svd(psf_Matrix);
        T=S./(S.^2+alpha);
        res_vector=blurred_vector*V*T*U';
        for k=1:61
        restoredimage_fft(i,j,k)=res_vector(k);
        end
    end
end
for i=1:61
    restoredimage(:,:,i)=abs(ifft2(ifftshift(restoredimage_fft(:,:,i))));
end
%% cut it and show
restoredimage_cut=zeros(testrows,testcols,61);
for i=1:61
    restoredimage_cut(:,:,i)=restoredimage(1:testrows,1:testcols,i);
end
figure(1)
for i=1:61
    imshow(restoredimage_cut(:,:,i));
    title([num2str(lambda(i)) 'nm spectral image'])
    pause(0.2);
end
save('restoredimage.mat','restoredimage_cut')
