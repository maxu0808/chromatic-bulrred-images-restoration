# chromatic-bulrred-images-restoration
This repository is relevent to the article "Hyperspectral imaging bioinspired by chromatic blur vision in color blind animals".
Here,we provide necessary MATLAB codes for readers to validate our ideas.
The funtion of these scripts are described below:
1.Stripe_target_generate.m: this script is to generate green stripe taget or green-blue stripe target,it also generates their spectral
images correspondly.The spectralimages will save as   "spectralimage.mat".
2.ZJU_target_generate.m:this script is to generate ZJU target with 6 colors and it's spectral images.
3.All_PSFs_calculate:this script is to calculate the PSFs of the blind animals' eye optical model.Run it will take much time,so we have 
been provide all the PSFs in the directory "All_PSF".
4.BlurredImage_generate.m:this is to generate chromatic blurred images under different accommadations,the bulrredimages are
saved as "blurredimage.mat".
5.SVD.m:this is to use SVD algorithm to restore the spectral images,the restored images will be saved as "restoredimage.mat".
6.plot_curve.m:it is to plot spectral cuves to compare ground truth spectral and restored spectral.
7.calculate_MSE_PSNR.m:it is to calculate MSE and PSNR.
8.PSF—3Dview.m:it is to plot 3D psf images.
The readers should run the scripts in the following order:target_generat,blurredimage_generate,SVD,plot_curve,MSE.
