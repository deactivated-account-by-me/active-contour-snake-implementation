
clear all;
snakeType = "Open";
N = 500;
alpha = 0.08;
beta = 0.5;
gamma = 0.8;
w_line = 2.5;
w_edge = 4.5;
w_term = 3.5;
sigma = 0.5;

% load image
I = imread('special_guassians.gif');

if (ndims(I) == 3)
    I = rgb2gray(I);
end

[x, y] = InitializeSnake(I,snakeType);

I_after_gaussian_filter = double(imgaussfilt(I, sigma));
%figure,imtool(I_after_gaussian_filter)
external_energy = ExternalEnergyCal(I_after_gaussian_filter, w_line, w_edge, w_term);
%figure,imtool(external_energy)

% calcualte the inverse matrix of A
a_inverse = InternalEnergyCal(size(x, 2), alpha, beta, gamma, snakeType);
x = x';
y = y';
sobel_x = [1 0 -1;2 0 -2; 1 0 -1];
sobel_y = [1 2 1; 0 0 0; -1 -2 -1];

fx = conv2(external_energy, sobel_x, 'same');
fy = conv2(external_energy, sobel_y, 'same');

steps = floor(N/30);



for i = 1:N
    [x, y] = iteration(a_inverse, x, y, external_energy, gamma, fx, fy);
    imshow(I);
    hold on;
    plot(x, y, 'r');
    if(mod(i, steps) == 0)
        fprintf('%d/%d interations\n', i, N);
    end
    pause(0.0001);
end

if(mod(i, steps) == 0)
    fprintf('%d/%d interations\n', N, N);
end
