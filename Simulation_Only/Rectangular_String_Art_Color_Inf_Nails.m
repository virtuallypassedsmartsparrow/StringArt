clc; clear all; close all;

% Created by Matthew James

% Description:
% Generates PNG image from string art (not nail list!). It assumes an
% infinite amount of nails. 

% Instructions:
% Save this .m file in a folder
% Save an image into the same folder
% Enter the name of the image in 'filename' as shown below
% (optional) edit User Inputs
% Click run

%% User Inputs:
filename = 'Euler.png'
ideal_image_size = 400;
Num_Nails = 300;

% Line properties
d = 0.036; %darkness of line
p_min = 0.00016; % minimum of p_line
tstart = 0.0014; % 0<= tstart <0 tend 
tend = 0.0161; % Thickness of web

% Greedy Algorithm properties
p_theshold = 0.0019; % if p(a,s) during any iteration gets smaller then this, break.
num_max_lines = 12000;

% Plotting
transparency = 0.06; %Each line darkness. Just for visualizng plot.

%% Image f
A = imread(filename);
BW_raw = rgb2gray(A); %returns a matrix of the intensity. 255 = white, 0 = black
[h,b] = size(BW_raw);
ratio = h/b;
b = 2; % base width of rectangle
h = ratio * b;
A = imresize(A,[round(ideal_image_size*ratio), ideal_image_size]);
A = double(A)/255;

% RGB to CMYK conversion
C = 1 - A(:,:,1);
M = 1 - A(:,:,2);
Y = 1 - A(:,:,3);
K = min(min(C, M), Y);
% Avoid division by zero for pure black
K(K == 1) = 0.999;
% Update CMY values based on K
C = (C - K) ./ (1 - K);
M = (M - K) ./ (1 - K);
Y = (Y - K) ./ (1 - K);
% Fill the CMYK channels
f_c = C;
f_m = M;
f_y = Y;
f_k = K;


%% Plot CMYK components of image
figure;
imagesc(A)
title('Raw Image')
axis square;

% Plot for the Cyan component
figure;
ax(1) = subplot(2,2,1);
imagesc(f_c);
colormap(ax(1),flipud([linspace(0,1,256)' linspace(1,1,256)' linspace(1,1,256)'])); % Smooth gradient from white to cyan
colorbar;
title('Cyan Component');
axis square;

% Plot for the Magenta component
ax(2) = subplot(2,2,2);
imagesc(f_m);
colormap(ax(2),flipud([linspace(1,1,256)' linspace(0,1,256)' linspace(1,1,256)'])); % Smooth gradient from white to magenta
colorbar;
title('Magenta Component');
axis square;

% Plot for the Yellow component
ax(3) = subplot(2,2,3);
imagesc(f_y);
colormap(ax(3),flipud([linspace(1,1,256)' linspace(1,1,256)' linspace(0,1,256)'])); % Smooth gradient from white to yellow
colorbar;
title('Yellow Component');
axis square;

% Plot for the Black (Key) component
ax(4) = subplot(2,2,4);
imagesc(f_k);
colormap(ax(4),[linspace(1,0,256)' linspace(1,0,256)' linspace(1,0,256)']); % Smooth gradient from white to black
colorbar;
title('Black Component');
axis square;

%% Calculations for CMYK
% Cyan
[Plot_Lines_x_c, Plot_Lines_y_c] = StringArtCalc(f_c,ideal_image_size,Num_Nails,d,p_min,tstart,tend,p_theshold,num_max_lines,b,h);
% Magenta
[Plot_Lines_x_m, Plot_Lines_y_m] = StringArtCalc(f_m,ideal_image_size,Num_Nails,d,p_min,tstart,tend,p_theshold,num_max_lines,b,h);
% Yellow
[Plot_Lines_x_y, Plot_Lines_y_y] = StringArtCalc(f_y,ideal_image_size,Num_Nails,d,p_min,tstart,tend,p_theshold,num_max_lines,b,h);
% Black
[Plot_Lines_x_k, Plot_Lines_y_k] = StringArtCalc(f_k,ideal_image_size,Num_Nails,d,p_min,tstart,tend,p_theshold,num_max_lines,b,h);

%% Plot
width = 800;
height = 800 * ratio;
figure('Name', 'String Art','color', 'w','Position', [50 10 width height]);
plot(Plot_Lines_x_c,Plot_Lines_y_c,'Color',[0,1,1,transparency])
hold on
plot(Plot_Lines_x_m,Plot_Lines_y_m,'Color',[1,0,1,transparency])
hold on
plot(Plot_Lines_x_y,Plot_Lines_y_y,'Color',[1,1,0,transparency])
hold on
plot(Plot_Lines_x_k,Plot_Lines_y_k,'Color',[0,0,0,transparency])
hold on
xlim([-b/2,b/2])
ylim([-h/2,h/2])
axis off;
box off;

%% Functions
function [Plot_Lines_x, Plot_Lines_y] = StringArtCalc(f,ideal_image_size,Num_Nails,d,p_min,tstart,tend,p_theshold,num_max_lines,b,h)
    %% Radon Transform image p
    alpha = linspace(0,180,Num_Nails);
    [p,s] = radon_fun(f,alpha,b,ideal_image_size);
    
    % Remove values of |s| >= max(L)/2
    ind_keep = abs(s) < sqrt(b^2+h^2)/2;
    s = s(ind_keep);
    p = p(ind_keep,:);
    
    alpha = deg2rad(alpha);
    [ALPHA, S] = meshgrid(alpha,s);
    
    L = zeros(size(ALPHA));
    for i = 1:length(alpha)
        for j = 1:length(s)
            L(j, i) = Length_Line_Rectangle_Calc(alpha(i),s(j),b,h);
        end
    end
    
    p = p./L;
    
    % Remove values of alpha and s that don't intersect the rectangle
    % deliberally make mask a bit larger to avoid algorithm selecting close to the edge
    b_crop = 0.99*b;
    h_crop = 0.99*h;
    Union1 = abs(S) >= (b_crop/2)*cos(ALPHA) + (h_crop/2)*sin(ALPHA) & 0<=ALPHA & ALPHA<=pi/2;
    Union2 = abs(S) >= (-b_crop/2)*cos(ALPHA) + (h_crop/2)*sin(ALPHA) & pi/2<=ALPHA & ALPHA<=pi;
    Not_in_rectangle = Union1 | Union2;
    p(Not_in_rectangle) = -inf; %artifically make "super white" to make algorithm not suggest it
    p_orig = p;
    
    figure
    imagesc(alpha,s,p)
    xlabel('\alpha [degrees]')
    ylabel('s [m]')
    title('Radon Transform of image p = R(f)')
    colormap(gca,hot), colorbar
    
    %% Calculations (Greedy Algorithm)
    
    % Add one line at a time by looking at brihtest parts of p
    size_p = size(p);
    for i = 1:num_max_lines
        % Find max p
        [p_max, ind] = max(p(:));
    
        if p_max < p_theshold
            disp('Threshold reached')
            break
        end
    
        [row, col] = ind2sub(size_p, ind);
        
        % Find corresponding alpha and s values at max p. Store these and
        % generate lines from them later.
        s0(i) = s(row);
        alpha0(i) = alpha(col);
        
        % Subtract p_line from p
        p_line_theory = p_line_fun(alpha0(i),s0(i),ALPHA,S,tstart,tend,d,p_min,b,h,L);
        p = p - p_line_theory;
    end
    disp(['Number of lines used = ', num2str(i)])
    
    
    %% Convert s,alpha --> psi_1, psi_2
    s = s0;
    alpha = alpha0;
    
    [Plot_Lines_x, Plot_Lines_y] = Plot_Lines_Rectangle_Calc(alpha,s,b,h);
end

function [p,s] = radon_fun(f,alpha,b,ideal_image_size)
    % Calculate radon transform of f. MATLAB assumes distance between
    % pixels is 1 so rescale.
    [p,s] = radon(f,alpha);
    p = p * b/ ideal_image_size;
    s = s * b/ ideal_image_size;
end

function p_line = p_line_fun(alpha0,s0,ALPHA,S,tstart,tend,d,p_min,b,h,L)
    % Calculates p(alpha,s) from theoretical model of line
    p_line = d*p_min ./((d.*L-p_min).*abs(sin(ALPHA-alpha0)) + p_min);
    mask = maskfun(alpha0,s0,ALPHA,S,b,h,tstart,tend);
    p_line = p_line.*mask;
end 

function mask = maskfun(alpha0,s0,ALPHA,S,b,h,tstart,tend)
    % Define where p_line = 0 and where it doesn't 
    [x_region, y_region] = p_regionfun(alpha0,s0,ALPHA,S,b,h,tstart);
    mask = ones(size(ALPHA));
    Too_Wide = x_region > 0;
    Too_Tall = y_region > 0;
    Union = Too_Wide | Too_Tall;
    mask(Union) = 0;
    n = 4;
    t = linspace(tstart,tend,n);
    for i = 1:n-1
        [x_region, y_region] = p_regionfun(alpha0,s0,ALPHA,S,b,h,t(i));
        Too_Wide = x_region > 0;
        Too_Tall = y_region > 0;
        Union1 = Too_Wide | Too_Tall;

        [x_region, y_region] = p_regionfun(alpha0,s0,ALPHA,S,b,h,t(i+1));
        Too_Wide = x_region < 0;
        Too_Tall = y_region < 0;
        Union2 = Too_Wide & Too_Tall;

        mask(Union1 & Union2) = (tend-t(i))/(tend-tstart);
    end
end

function [x_region, y_region] = p_regionfun(alpha0,s0,ALPHA,S,b,h,t)
    % x_region > 0 is where p_line = 0
    % y_region > 0 is where p_line = 0
    x_region = abs(S.*sin(alpha0)-s0*sin(ALPHA)) - (b/2).*abs(sin(ALPHA-alpha0)) - (t/2)*abs(sin(alpha0));
    y_region = abs(S.*cos(alpha0)-s0*cos(ALPHA)) - (h/2).*abs(sin(ALPHA-alpha0)) - (t/2)*abs(cos(alpha0));
end

function [Plot_Lines_x, Plot_Lines_y] = Plot_Lines_Rectangle_Calc(alpha,s,b,h)
    % Preallocation for speed
    num_lines = length(alpha);
    Start_Line_x = zeros(1, num_lines);
    Start_Line_y = zeros(1, num_lines);
    End_Line_x = zeros(1, num_lines);
    End_Line_y = zeros(1, num_lines);

    % Calculate (x,y) of start and end of line along rectangle edge
    for i = 1:num_lines
        alpha_local = alpha(i);
        s_local = s(i);
        [x1,y1,x2,y2] = Start_End_Line_Rectangle_Calc(alpha_local,s_local,b,h);
        
        % Store all lines
        Start_Line_x(1,i) = x1;
        Start_Line_y(1,i) = y1;
        End_Line_x(1,i) = x2;
        End_Line_y(1,i) = y2;
    end
    
    Plot_Lines_x = [Start_Line_x; End_Line_x];
    Plot_Lines_y = [Start_Line_y; End_Line_y];
end

function L = Length_Line_Rectangle_Calc(alpha_local,s_local,b,h)
    [x1,y1,x2,y2] = Start_End_Line_Rectangle_Calc(alpha_local,s_local,b,h);
    L = sqrt((x2-x1)^2+(y2-y1)^2);
end

function [x1,y1,x2,y2] = Start_End_Line_Rectangle_Calc(alpha_local,s_local,b,h)
        side_count = 0;

        % Initialize
        x1 = NaN;
        y1 = NaN;
        x2 = NaN;
        y2 = NaN;
    
        % Right side
        x = b/2;
        y = (s_local - 0.5*b*cos(alpha_local))./sin(alpha_local);
        if -0.5*h <= y && y <= 0.5*h
            side_count = side_count + 1;
            if side_count == 1
                x1 = x;
                y1 = y;
            elseif side_count == 2
                x2 = x;
                y2 = y;
            end
        end
        
        % Top side
        x = (s_local - 0.5*h*sin(alpha_local))./cos(alpha_local);
        y = h/2;
        if -0.5*b <= x && x <= 0.5*b
            side_count = side_count + 1;
            if side_count == 1
                x1 = x;
                y1 = y;
            elseif side_count == 2
                x2 = x;
                y2 = y;
            end
        end
        
        % Left side
        x = -b/2;
        y = (s_local + 0.5*b*cos(alpha_local))./sin(alpha_local);
        if -0.5*h <= y && y <= 0.5*h
            side_count = side_count + 1;
            if side_count == 1
                x1 = x;
                y1 = y;
            elseif side_count == 2
                x2 = x;
                y2 = y;
            end
        end
    
        % Bot side
        x = (s_local + 0.5*h*sin(alpha_local))./cos(alpha_local);
        y = -h/2;
        if -0.5*b <= x && x <= 0.5*b
            side_count = side_count + 1;
            if side_count == 1
                x1 = x;
                y1 = y;
            elseif side_count == 2
                x2 = x;
                y2 = y;
            end
        end
end