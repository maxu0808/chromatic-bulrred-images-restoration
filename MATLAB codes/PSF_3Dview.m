clear all;close all;clc;
PSF_450nm=load('./All_PSF/450nm.mat');
PSF_550nm=load('./All_PSF/550nm.mat');
PSF_650nm=load('./All_PSF/650nm.mat');
PSF_450nm=PSF_450nm.psf;
PSF_550nm=PSF_550nm.psf;
PSF_650nm=PSF_650nm.psf;
a=PSF_450nm(:,:,11);
b=PSF_450nm(:,:,31);
c=PSF_450nm(:,:,51);
d=PSF_550nm(:,:,11);
e=PSF_550nm(:,:,31);
f=PSF_550nm(:,:,51);
g=PSF_650nm(:,:,11);
h=PSF_650nm(:,:,31);
i=PSF_650nm(:,:,51);
x_axis=0:100;
y_axis=0:100;
subplot(3,3,1)
mesh(x_axis,y_axis,a)
xlabel('x-axis','rotation',-40);ylabel('y-axis','rotation',20);zlabel('Intensity')
title({'wavelength = 450 nm','accommadation = -0.32 mm'});
view(55,30)
colorbar
subplot(3,3,2)
mesh(x_axis,y_axis,b)
xlabel('x-axis','rotation',-40);ylabel('y-axis','rotation',20);zlabel('Intensity')
title({'wavelength = 450 nm','accommadation = 0.00 mm'});
view(55,30)
colorbar
subplot(3,3,3)
mesh(x_axis,y_axis,c)
xlabel('x-axis','rotation',-40);ylabel('y-axis','rotation',20);zlabel('Intensity')
title({'wavelength = 450 nm','accommadation = 0.16 mm'});
view(55,30)
colorbar
subplot(3,3,4)
mesh(x_axis,y_axis,d)
xlabel('x-axis','rotation',-40);ylabel('y-axis','rotation',20);zlabel('Intensity')
title({'wavelength = 550 nm','accommadation = -0.32 mm'});
view(55,30)
colorbar
subplot(3,3,5)
mesh(x_axis,y_axis,e)
xlabel('x-axis','rotation',-40);ylabel('y-axis','rotation',20);zlabel('Intensity')
title({'wavelength = 550 nm','accommadation = 0.00 mm'});
view(55,30)
colorbar
subplot(3,3,6)
mesh(x_axis,y_axis,f)
colorbar
xlabel('x-axis','rotation',-40);ylabel('y-axis','rotation',20);zlabel('Intensity')
title({'wavelength = 550 nm','accommadation = 0.16 mm'});
view(55,30)
subplot(3,3,7)
mesh(x_axis,y_axis,g)
colorbar
xlabel('x-axis','rotation',-40);ylabel('y-axis','rotation',20);zlabel('Intensity')
title({'wavelength = 650 nm','accommadation = -0.32 mm'});
view(55,30)
subplot(3,3,8)
mesh(x_axis,y_axis,h)
xlabel('x-axis','rotation',-40);ylabel('y-axis','rotation',20);zlabel('Intensity')
title({'wavelength = 650 nm','accommadation = 0.00 mm'});
view(55,30)
colorbar
subplot(3,3,9)
mesh(x_axis,y_axis,i)
xlabel('x-axis','rotation',-40);ylabel('y-axis','rotation',20);zlabel('Intensity')
title({'wavelength = 650 nm','accommadation = 0.16 mm'});
view(55,30)
colorbar

