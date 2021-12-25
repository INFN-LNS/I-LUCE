%COLOR HISTOGRAM EQUALIZATION

 disp('Open one RCF image.')
 [filename, path] = imgetfile;

%READ THE INPUT IMAGE
I = imread(filename);

%CONVERT THE RGB IMAGE INTO HSV IMAGE FORMAT
HSV = rgb2hsv(I);
%https://www.imageeprocessing.com/2013/05/converting-rgb-image-to-hsi.html


%PERFORM HISTOGRAM EQUALIZATION ON INTENSITY COMPONENT
Heq = histeq(HSV(:,:,3));
%https://www.imageeprocessing.com/2011/04/matlab-code-histogram-equalization.html

HSV_mod = HSV;
HSV_mod(:,:,3) = Heq;

RGB = hsv2rgb(HSV_mod);
%https://www.imageeprocessing.com/2013/06/convert-hsi-image-to-rgb-image.html

figure
imshow(I);title('Before Histogram Equalization');
figure
imshow(RGB);title('After Histogram Equalization');

I_gray = rgb2gray(I);
imshow(I_gray)
colormap(gca, jet(100));
colorbar(gca);

imhist(I_gray)

bw = imbinarize(I_gray);

imshow(bw)

