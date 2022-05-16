clear;clc;
imdata = imread("ferrari.jpg");
imdata = rgb2gray(imdata);
imdata = im2double(imdata);

[u, s, v] = svd(imdata);
approx = 40;
imgapprox = u(:,1:approx)*s(1:approx,1:approx)*v(:,1:approx)';
figure
imshow(imdata)