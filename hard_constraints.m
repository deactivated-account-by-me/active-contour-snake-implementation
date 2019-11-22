clear all

N = 1000;
alpha = 0.002;
beta = 0.2;
gamma = 0.1;
w_line = 1;
w_edge = 0.1;
w_term = 0.1;
sigma = 0.5;
kappa = 0.15;
global fix_point
fix_point = double.empty(3,0);
global first_click;
first_click = 0;
k = 1.25;
% load image
I = imread('simple_image_g_xy.gif');
if (ndims(I) == 3)
    I = rgb2gray(I);
end

[x, y] = initialize_open_snake(I); % open snake initialization

Image_after_gaussian_filter = double(imgaussfilt(I, sigma));
external_energy = external_energy_calculation(Image_after_gaussian_filter, w_line, w_edge, w_term);

% calculate the inverse matrix of A
a_inverse = internal_energy_calculation(size(x, 2), alpha, beta, gamma);

x = x';
y = y';
sobel_x = [1 0 -1;2 0 -2; 1 0 -1];
sobel_y = [1 2 1; 0 0 0; -1 -2 -1];

fx = conv2(external_energy, sobel_x, 'same');
fy = conv2(external_energy, sobel_y, 'same');

steps = floor(N/30);
set(gcf,'WindowButtonDownFcn',{@ButttonDownFcn});

for i=1:N
    global click;
    global click_x;
    global click_y;
    global index_new_point;
    if (click == 1)
        distance_measurement_vector = euclidean_distance([x y], [click_x click_y]);
        close_point = min(distance_measurement_vector);
        index_new_point = find(distance_measurement_vector == close_point);
        force_x = click_x - x(index_new_point);
        force_y = click_y - y(index_new_point);
        x(index_new_point) = k * force_x + x(index_new_point);
        y(index_new_point) = k * force_y + y(index_new_point);
        [x, y] = iteration(a_inverse, x, y, external_energy, gamma, fx, fy, kappa);
         %Displaying the snake in its new position
        imshow(I,[]); 
        hold on;
    
        %plot([xs; xs(1)], [ys; ys(1)], 'b-'); for closed
        plot(x, y, 'r');%for open
    
        if(mod(i, steps) == 0)
            fprintf('%d/%d interations\n', i, N);
        end
        pause(0.1);    
    else
        [x, y] = iteration(a_inverse, x, y, external_energy, gamma, fx, fy, kappa);
         %Displaying the snake in its new position
        imshow(I,[]); 
        hold on;
    
        %plot([xs; xs(1)], [ys; ys(1)], 'b-'); for closed
        plot(x, y, 'r');%for open
    
        if(mod(i, steps) == 0)
            fprintf('%d/%d interations\n', i, N);
        end
        pause(0.1);    
    end
end


function [x, y] = initialize_open_snake(I)
    global knots;
    imshow(I);
    max_len = max(size(I)) - 1;
    axis([0 1 0 1]);
    imshow(I);
    [x, y] = getpts();
    x = transpose(x);
    y = transpose(y);
    hold on;
    % x=x';y=y';
    % temp=[x(1);y(1)];
    knots = [x ; y];
    %xy=[xy,temp];%for closed 
    number_of_points = length(x);
    distance_points = 1:number_of_points;
    final_distance_points = 1: 0.45: number_of_points;
    open_curve = spline(distance_points,knots,final_distance_points);
    open_curve(open_curve < 1) = 1;
    open_curve(open_curve > max_len) = max_len;
    x_new = open_curve(1,:);
    y_new = open_curve(2,:);
    plot(x ,y, 'o', x_new, y_new, '--');
    x = x_new;
    y = y_new;
    hold on ;
    
    
end
function ButttonDownFcn(src,event)
    global click;
    global click_x;
    global click_y;
    global first_click
    first_click =1;
    pt = get(gca,'CurrentPoint'); 
    click_x = pt(1,1);
    click_y = pt(1,2);
    click = 1;
end



    
function [E_external] = external_energy_calculation(I, w_line, w_edge, w_term)
    E_line = double(I);
    
    [gradient_x, gradient_y] = gradient(I);
    E_edge = -1 * (gradient_x.^2 + gradient_y.^2);
    
    sobel_x = [1 0 -1; 2 0 -2; 1 0 -1];
    sobel_y = [1 2 1; 0 0 0; -1 -2 -1];
    
    c_x = conv2(I, sobel_x, 'same');
    c_xx = conv2(c_x, sobel_x, 'same');
    c_y = conv2(I, sobel_y, 'same');
    c_yy = conv2(c_y, sobel_y, 'same');
    c_xy = conv2(c_x, sobel_y, 'same');
    %c_x = imgaussfilt(I,[4, 1]);
    %c_xx = imgaussfilt(c_x, [4, 1]);
    %c_y = imgaussfilt(I,[1, 4]);
    %c_yy = imgaussfilt(c_y, [1, 4]);
    %c_xy = imgaussfilt(c_x, [1, 4]);
    
   
    
    E_term = ((c_yy.*(c_x.*c_x)) - (2.*c_xy.*c_xy.*(c_x.*c_y)) + (c_xx.*(c_y.*c_y)))./((c_x.*c_x)+(c_y.*c_y)).^ 1.5;
    
    E_external = w_line * E_line + w_term * E_term + w_edge + E_edge;
    
    E_external = E_external / max(E_external(:));
    
    
end

function Ainv = internal_energy_calculation(points,alpha ,beta, gamma)

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
    a(1,1) = alpha + beta;
    a(1,2) = -alpha - 2 * beta;
    a(2,1) = -alpha - 2 * beta;
    a(2,2) = 2*alpha + 5 *beta;
        
    a(points - 1 , 1) = 0;
    a(points, 1) = 0;
    a(points-1 , 2) = 0;
    a(points, 2) = 0;
        
    a(points,points) = alpha + beta;
    a(points-1, points) = -alpha - 2 * beta;
    a(points,points - 1) = -alpha - 2 * beta;
    a(points-1, points-1) = 2*alpha + 5 * beta;
        
    a(1, points - 1) = 0;
    a(1, points) = 0;
    a(2, points-1) = 0;
    a(2, points) = 0;

    Ainv=inv(a + gamma.* eye(points));
end

function [new_x, new_y]= iteration(a_inverse, x, y, external_energy, gamma, fx, fy, kappa)
    global a_x
    global a_y
    global fix_point
    global first_click
    global index_new_point
    global click_x
    global click_y
    global click
    if(first_click==1)
        new_x = a_inverse * (gamma * x - kappa*interp2(fx, x, y));
        new_y = a_inverse * (gamma * y - kappa*interp2(fy, x, y));
    
   
        max_x = max(x);
        max_y = max(y);
        new_y(new_y < 1) = 1;
        new_y(new_y > max_y) = max_y;
        new_x(new_x < 1) = 1;
        new_x(new_x > max_x) = max_x;
    
        a_x = new_x;
        a_y = new_y;
        
        a_x(index_new_point) = click_x;
        a_y(index_new_point) = click_y;
        single_fix_point = [a_x(index_new_point); a_y(index_new_point); index_new_point];
        fix_point = [fix_point single_fix_point];
   
    
        for k = 1:size(fix_point,2)
            i = fix_point(3, k);
            t_x = fix_point(1, k);
            t_y = fix_point(2, k);
            new_x(i) = t_x;
            new_y(i) = t_y;
        end
    else
        new_x = a_inverse * (gamma * x - kappa*interp2(fx, x, y));
        new_y = a_inverse * (gamma * y - kappa*interp2(fy, x, y));
    
   
        max_x = max(x);
        max_y = max(y);
        new_y(new_y < 1) = 1;
        new_y(new_y > max_y) = max_y;
        new_x(new_x < 1) = 1;
        new_x(new_x > max_x) = max_x;
    
    end
    
    
   click = 0;
    


end

function distance_array = euclidean_distance(vector_1, vector_2)
    d = abs(vector_1 - vector_2);
    distance_array = [];
    for k1 = 1:length(d)
        distance = norm(d(k1));
        distance_array = [distance_array, distance];
    end
end


