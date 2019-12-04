clear all;

N = 500;
alpha = 0.2;
beta = 0.05;
gamma = 0.2;
w_line = 2.5;
w_edge = 4;
w_term = 3.5;
sigma = 0.5;
k = 1.25;
% load image
I = imread('circle.jpg');
if (ndims(I) == 3)
    I = rgb2gray(I);
end

[x, y] = initialize_snake(I);

I_after_gaussian_filter = double(imgaussfilt(I, sigma));
external_energy = external_energy_calculation(I_after_gaussian_filter, w_line, w_edge, w_term);

% calcualte the inverse matrix of A
a_inverse = internal_energy_calculation(size(x, 2), alpha, beta, gamma);
x = x';
y = y';
sobel_x = [1 0 -1;2 0 -2; 1 0 -1];
sobel_y = [1 2 1; 0 0 0; -1 -2 -1];

fx = conv2(external_energy, sobel_x, 'same');
fy = conv2(external_energy, sobel_y, 'same');

steps = floor(N/30);
set(gcf,'WindowButtonDownFcn',{@ButttonDownFcn}, 'WindowButtonUpFcn',{@ButttonUpFcn});

for i = 1:N
    global click;
    global click_x;
    global click_y;
    if (click == 1)
        distance_measurement_vector = euclidean_distance([x y], [click_x click_y]);
        close_point = min(distance_measurement_vector);
        index = find(distance_measurement_vector == close_point);
        force_x = click_x - x(index);
        force_y = click_y - y(index);
        x(index) = k * force_x + x(index);
        y(index) = k * force_y + y(index);
        click = 0;
    end
    [x, y] = iteration(a_inverse, x, y, external_energy, gamma, fx, fy);
    imshow(I);
    hold on;
    plot(x, y, 'r');
    if(mod(i, steps) == 0)
        fprintf('%d/%d interations\n', i, N);
    end
    pause(0.1);
end

if(mod(i, steps) == 0)
    fprintf('%d/%d interations\n', N, N);
end

function [x, y] = initialize_snake(I)
    % here is to show the snake in the space
    fig = figure;
    imshow(I);
    max_len = max(size(I)) - 1;
    axis([0 1 0 1]);
    imshow(I);
    [x, y] = getpts();
    hold on;
    x = transpose(x);
    y = transpose(y);
    x = [x, x(1)];
    y = [y, y(1)];
    knots = [x ; y];
    number_of_points = length(x);
    distance_points = 1:number_of_points;
    final_distance_points = 1:0.05:number_of_points;
    closed_curve = spline(distance_points, knots, final_distance_points);
    closed_curve(closed_curve < 1) = 1;
    closed_curve(closed_curve > max_len) = max_len;
    x_new = closed_curve(1,:);
    y_new = closed_curve(2,:);
    plot(x ,y, 'o', x_new, y_new, '--');
    x = x_new;
    y = y_new;
    hold on;
end

function ButttonDownFcn(src,event)
    set(gcf,'WindowButtonMotionFcn',{@ButtonMoveFcn});
end

function ButttonUpFcn(src,event)
    set(gcf,'WindowButtonMotionFcn',{});
end

function ButtonMoveFcn(src, event)
    global click;
    global click_x;
    global click_y;
    pt = get(gca,'CurrentPoint'); 
    click_x = pt(1,1);
    click_y = pt(1,2);
    click = 1;
end

function [E_external] = external_energy_calculation(I, w_line, w_edge, w_term)
    E_line = double(I);
    
    [gradient_x, gradient_y] = gradient(I);
    E_edge = -1 * (gradient_x.^2 + gradient_y.^2);
    
    %sobel_x = [1 0 -1; 2 0 -2; 1 0 -1];
    %sobel_y = [1 2 1; 0 0 0; -1 -2 -1];
    
    %c_x = conv2(I, sobel_x, 'same');
    %c_xx = conv2(c_x, sobel_x, 'same');
    %c_y = conv2(I, sobel_y, 'same');
    %c_yy = conv2(c_y, sobel_y, 'same');
    %c_xy = conv2(c_x, sobel_y, 'same');
    c_x = imgaussfilt(I,[4, 1]);
    c_xx = imgaussfilt(c_x, [4, 1]);
    c_y = imgaussfilt(I,[1, 4]);
    c_yy = imgaussfilt(c_y, [1, 4]);
    c_xy = imgaussfilt(c_x, [1, 4]);
    
   
    
    E_term = ((c_yy.*(c_x.*c_x)) - (2.*c_xy.*c_xy.*(c_x.*c_y)) + (c_xx.*(c_y.*c_y)))./((c_x.*c_x)+(c_y.*c_y)).^ 1.5;
    
    E_external = w_line * E_line + w_term * E_term + w_edge + E_edge;
    
    E_external = E_external / max(E_external(:));
    
    
end

function [a_inverse] = internal_energy_calculation(points, alpha, beta, gamma)
    a = zeros(points, points);
    a_row = zeros(1, points);
    a_row(1) = 2 * alpha + 6 * beta;
    a_row(2) = -alpha - 4 * beta;
    a_row(3) = beta;
    a_row(points-1) = beta;
    a_row(points) = -alpha - 4 * beta;
    
    
    for i = 1: points
        a(i,:) = a_row;
        a_row = [a_row(points) a_row(1:points-1)];
    end
    
    a_inverse = inv(a + gamma * eye(points, points));
end

function [new_x, new_y]= iteration(a_inverse, x, y, external_energy, gamma, fx, fy)
    new_x = a_inverse * (gamma * x - 0.15*interp2(fx, x, y));
    new_y = a_inverse * (gamma * y - 0.15*interp2(fy, x, y));
    
    max_x = max(x);
    max_y = max(y);
    new_y(new_y < 1) = 1;
    new_y(new_y > max_y) = max_y;
    new_x(new_x < 1) = 1;
    new_x(new_x > max_x) = max_x;
end

function distance_array = euclidean_distance(vector_1, vector_2)
    d = abs(vector_1 - vector_2);
    distance_array = [];
    for k1 = 1:length(d)
        distance = norm(d(k1));
        distance_array = [distance_array, distance];
    end
end

