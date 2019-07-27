%%% This is a MATLAB script for calculate the PSFs of color blind animals'
%%% lens model.The method is:we firstly set the incident light as a certain wavelength,
%%% we then calculate the three dimensional field distribution around the
%%% focus by Lomel function,the specific formula can refer to "Born, M.
%%% Principles of Optics,chapter 4:10".
%%% It should be stressed that the PSFs' computing is very time-consuming. 
%%% it can be mitigated by using parallel computing,however,we don't
%%% contribute the parallel computing program.We have been calculate
%%% the PSFs in "All_PSF" directory.

clear all;
close all;
clc;

a=4;% The radius of the pupil£¬mm
pixel=5e-3; %the detector pixel size 5um*5um
wavelength=400:5:700; %Wavelengths need to be calculate PSFs,nm
focallength=12;  %Take 550nm focused position as zero,the focallength of 550nm is 12mm. 
focusshift=-5.4676e-05*wavelength.^2 + 0.0794*wavelength -27.1047; % percent focus shift
focusshiftmm=1E-2*focusshift*focallength; 
wavelength2=wavelength.*1e-6;  %convert from nm to mm
psf_rows=101;   %set each psf image size is 101*101
psf_cols=101;
PSF_all=cell(61,61); 

%% Choose the incident light wavelength as a certain wavelength,
%%% and we hold it invarant when we calculate diffrent positions' PSF.
%%% The diffrent positions is determined by other wavelengths focused positions.
%%% we then cycle computing all PSFs by choose another incident light.
for foucuslambda=2:length(wavelength2)-1
foucuslambda
lambda=wavelength2(foucuslambda);
%the choosed incident wavelength's focallength
f=12+focusshiftmm(foucuslambda);
%% compute the PSF on the focal plane
    for maskrow=1:psf_rows   %
%     making_pupil_mask_at_row=maskrow
         for maskcol=1:psf_rows
            masky=maskrow-(psf_rows+1)/2; 
            maskx=maskcol-(psf_cols+1)/2;
            rpix=sqrt((maskx).^2+(masky).^2).*pixel;  % distance from center of pupil
            v=2*pi*a/(lambda*f)*rpix;                
            intensity2=(2*besselj(1,v)./v).^2;
            pupilmask(maskrow,maskcol)=intensity2;
         end 
    end  
pupilmask((psf_rows+1)/2,(psf_cols+1)/2)=1;
psf(:,:,foucuslambda)=pupilmask;
%% compute the PSFs when the defocusing direction is positive
    for i=foucuslambda+1:61
       u=(focusshiftmm(i)-focusshiftmm(foucuslambda))*(2*pi*a*a)/(lambda*f*f);
       for maskrow=1:psf_rows
    %     making_pupil_mask_at_row=maskrow
         for maskcol=1:psf_cols
            masky=maskrow-(psf_rows+1)/2; 
            maskx=maskcol-(psf_cols+1)/2;
            rpix=sqrt((maskx).^2+(masky).^2).*pixel;  % distance from center of pupil
            v=2*pi*a/(lambda*f)*rpix;
            V0=0;
            V1=0;
            U1=0;
            U2=0;
            for m=0:1:50
            V0=V0+(-v.^2/u.^2).^m.*besselj(2*m,v);
            V1=V1+(-1).^m.*(v./u).^(2*m+1).*besselj(2*m+1,v);
            U1=U1+(-1).^m.*(u./v).^(2*m+1).*besselj(2*m+1,v);
            U2=U2+(-1).^m.*(u./v).^(2*m+2).*besselj(2*m+2,v);
            end
            if v<u
            intensity2=(2/u).^2.*(1+V0.^2+V1.^2-2*V0*cos(0.5*(u+v.^2/u))-2*V1*sin(0.5*(u+v.^2/u)));
            else
            intensity2=(2/u).^2.*(U1.^2+U2.^2);
            end
            pupilmask(maskrow,maskcol)=intensity2;
        end 
       end  
    psf(:,:,i)=pupilmask;
    end

%% compute the PSFs when the defocusing direction is negetive
    for i=1:foucuslambda-1
       u=abs((focusshiftmm(i)-focusshiftmm(foucuslambda)))*(2*pi*a*a)/(lambda*f*f);
       for maskrow=1:psf_rows
    %     making_pupil_mask_at_row=maskrow
         for maskcol=1:psf_cols
            masky=maskrow-(psf_rows+1)/2; 
            maskx=maskcol-(psf_cols+1)/2;
            rpix=sqrt((maskx).^2+(masky).^2).*pixel;  % distance from center of pupil
            v=2*pi*a/(lambda*f)*rpix;
            V0=0;
            V1=0;
            U1=0;
            U2=0;
            for m=0:1:50
            V0=V0+(-v.^2/u.^2).^m.*besselj(2*m,v);
            V1=V1+(-1).^m.*(v./u).^(2*m+1).*besselj(2*m+1,v);
            U1=U1+(-1).^m.*(u./v).^(2*m+1).*besselj(2*m+1,v);
            U2=U2+(-1).^m.*(u./v).^(2*m+2).*besselj(2*m+2,v);
            end
            if v<u
            intensity2=(2/u).^2.*(1+V0.^2+V1.^2-2*V0*cos(0.5*(u+v.^2/u))-2*V1*sin(0.5*(u+v.^2/u)));
            else
            intensity2=(2/u).^2.*(U1.^2+U2.^2);
            end
            pupilmask(maskrow,maskcol)=intensity2;
        end 
       end  
    psf(:,:,i)=pupilmask;
end
path ='.\All_PSF2';
save([path num2str(wavelength(foucuslambda)) 'nm.mat'],'psf')
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% Head and Tail
for foucuslambda=1:1
foucuslambda
lambda=wavelength2(foucuslambda);
f=12+focusshiftmm(foucuslambda);
%% compute the PSF on the focal plane
    for maskrow=1:psf_rows
    %     making_pupil_mask_at_row=maskrow
         for maskcol=1:psf_cols
            masky=maskrow-(psf_rows+1)/2; 
            maskx=maskcol-(psf_cols+1)/2;
            rpix=sqrt((maskx).^2+(masky).^2).*pixel;  % distance from center of pupil
            v=2*pi*a/(lambda*f)*rpix;
            intensity2=(2*besselj(1,v)./v).^2;
            pupilmask(maskrow,maskcol)=intensity2;
        end 
    end  
pupilmask((psf_rows+1)/2,(psf_cols+1)/2)=1;
psf(:,:,foucuslambda)=pupilmask;
%% compute the PSFs when the defocusing direction is positive
    for i=foucuslambda+1:61
       u=(focusshiftmm(i)-focusshiftmm(foucuslambda))*(2*pi*a*a)/(lambda*f*f);
       for maskrow=1:101
    %     making_pupil_mask_at_row=maskrow
         for maskcol=1:101
            masky=maskrow-102/2; 
            maskx=maskcol-102/2;
            rpix=sqrt((maskx).^2+(masky).^2).*pixel;  % distance from center of pupil
            v=2*pi*a/(lambda*f)*rpix;
            V0=0;
            V1=0;
            U1=0;
            U2=0;
            for m=0:1:50
            V0=V0+(-v.^2/u.^2).^m.*besselj(2*m,v);
            V1=V1+(-1).^m.*(v./u).^(2*m+1).*besselj(2*m+1,v);
            U1=U1+(-1).^m.*(u./v).^(2*m+1).*besselj(2*m+1,v);
            U2=U2+(-1).^m.*(u./v).^(2*m+2).*besselj(2*m+2,v);
            end
            if v<u
            intensity2=(2/u).^2.*(1+V0.^2+V1.^2-2*V0*cos(0.5*(u+v.^2/u))-2*V1*sin(0.5*(u+v.^2/u)));
            else
            intensity2=(2/u).^2.*(U1.^2+U2.^2);
            end
            pupilmask(maskrow,maskcol)=intensity2;
        end 
       end  
psf(:,:,i)=pupilmask;
end
path ='.\All_PSF2\';
save([path num2str(wavelength(foucuslambda)) 'nm.mat'],'psf')
end

for foucuslambda=61:61
foucuslambda
lambda=wavelength2(foucuslambda);
f=12+focusshiftmm(foucuslambda);
%% compute the PSF on the focal plane
    for maskrow=1:psf_rows
    %     making_pupil_mask_at_row=maskrow
         for maskcol=1:psf_cols
            masky=maskrow-(psf_rows+1)/2; 
            maskx=maskcol-(psf_cols+1)/2;
            rpix=sqrt((maskx).^2+(masky).^2).*pixel;  % distance from center of pupil
            v=2*pi*a/(lambda*f)*rpix;
            intensity2=(2*besselj(1,v)./v).^2;
            pupilmask(maskrow,maskcol)=intensity2;
        end 
    end  
pupilmask((psf_rows+1)/2,(psf_cols+1)/2)=1;
psf(:,:,foucuslambda)=pupilmask;
%% compute the PSFs when the defocusing direction is negetive
    for i=1:foucuslambda-1
       u=abs((focusshiftmm(i)-focusshiftmm(foucuslambda)))*(2*pi*a*a)/(lambda*f*f);
       for maskrow=1:psf_rows
    %     making_pupil_mask_at_row=maskrow
         for maskcol=1:psf_cols
            masky=maskrow-(psf_rows+1)/2; 
            maskx=maskcol-(psf_cols+1)/2;
            rpix=sqrt((maskx).^2+(masky).^2).*pixel;  % distance from center of pupil
            v=2*pi*a/(lambda*f)*rpix;
            V0=0;
            V1=0;
            U1=0;
            U2=0;
            for m=0:1:50
            V0=V0+(-v.^2/u.^2).^m.*besselj(2*m,v);
            V1=V1+(-1).^m.*(v./u).^(2*m+1).*besselj(2*m+1,v);
            U1=U1+(-1).^m.*(u./v).^(2*m+1).*besselj(2*m+1,v);
            U2=U2+(-1).^m.*(u./v).^(2*m+2).*besselj(2*m+2,v);
            end
            if v<u
            intensity2=(2/u).^2.*(1+V0.^2+V1.^2-2*V0*cos(0.5*(u+v.^2/u))-2*V1*sin(0.5*(u+v.^2/u)));
            else
            intensity2=(2/u).^2.*(U1.^2+U2.^2);
            end
            pupilmask(maskrow,maskcol)=intensity2;
        end 
       end  
psf(:,:,i)=pupilmask;
end
path ='.\All_PSF2\';
save([path num2str(wavelength(foucuslambda)) 'nm.mat'],'psf')
end