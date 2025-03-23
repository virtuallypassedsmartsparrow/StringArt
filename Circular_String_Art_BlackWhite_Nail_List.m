clc; clear all; close all;

% Created by Matthew James

% Instructions:
% Save this .m file in a folder
% Save an image in the same folder
% Enter the name of the image in 'filename' in User Inputs
% Enter in number of nails in User Inputs
% (optional) edit Advanced User Inputs
% Click run
% Note this will save a .txt file which contains the nail list order


%% User Inputs:
filename = 'mona_lisa.png'      % Write the image filename here
Num_Nails = 400;                % Enter the number of nails around the Canvas

%% Advanced User Inputs:
% Image properties
ideal_image_size = 400;         % This compresses the image size. Larger is better but slower.

% Line properties
d = 0.036;                      % Darkness of line. The line is modelled with thickness with negligible darkness at the edge, and maximum darkness, d, in the center
p_min = 0.00016;                % Minimum value of radon transform of a line. Effectively the 'area' of darknes under the line with a perpendicular slice
tstart = 0.0014;                % 0<= tstart < tend. tstart is the thickness value of the line at which the darkness begins to decay as it approaches the edge of the line
tend = 0.0161;                  % tend is the thickness of the line. At this thickness, the line effectively has no darkness.

% Greedy Algorithm properties
p_theshold = 0.0037;            % If p(a,s) during any iteration gets smaller then this, break. In theory we want p_threshold = 0, but in practice this is not ideal
num_max_lines = 9000;           % Algorithm will halt if number of lines > num_max_lines

% Plotting
transparency = 0.06;            % Each line darkness. Just for visualizng plot. If you are using thick string, increase this value to match how dark it appears 

%% Image f
R = 1; % Radius

A = imread(filename);
BW_raw = rgb2gray(A); %returns a matrix of the intensity. 255 = white, 0 = black
BW = imresize(BW_raw,[ideal_image_size,ideal_image_size]);
BW = 1-double(BW)/255;
f = BW;

%Image crop
x = linspace(-R,R,ideal_image_size);
y = linspace(-R,R,ideal_image_size);
[X,Y] = meshgrid(x,y);
f(X.^2+Y.^2>R^2) = 0;

figure
imagesc(f)
title('Image f')
colormap(flipud(gray(256)));
colorbar
axis square

%% Radon Transform image p
alpha = linspace(0,180,3*Num_Nails); %Increase resolution of alpha for R(image) > 2*N for better interpolation
[p,s] = radon_fun(f,alpha,R,ideal_image_size);

% Remove values of |s| >= R
ind_keep = abs(s) < R;
s = s(ind_keep);
p = p(ind_keep,:);

% Technically I only care about p/L = R[f] per unit length of line
alpha = deg2rad(alpha);
[ALPHA, S] = meshgrid(alpha,s);
L = 2*sqrt(R^2-S.^2); %length of lines
p = p./L;

% Create ideal uniform psi_1 and psi_2 grid from Num_Nails
psi_1 = linspace(-pi, pi, Num_Nails+1);       % from -pi to pi
psi_2 = linspace(0, 2*pi, Num_Nails+1);       % from 0 to 2pi
[PSI_1, PSI_2] = meshgrid(psi_1, psi_2);
L = 2*R*sin(abs(PSI_2-PSI_1)/2);

% Convert radon of image to psi_1 psi_2 by interpolation
p = AlphaS2Phi(ALPHA,S,PSI_1,PSI_2,p,R);

figure
h = pcolor(PSI_1,PSI_2,p);
xlabel('\psi_1 [rad]')
ylabel('\psi_2 [rad]')
title('Radon Transform of image p = R(f)')
set(h, 'EdgeColor', 'none');
colormap(gca,hot), colorbar

%% Calculations (Greedy Algorithm)

% Add one line at a time by looking at brihtest parts of p
size_p = size(p);
psi_10 = NaN;  % Initial placeholder
psi_20 = NaN;  % Initial placeholder

for i = 1:num_max_lines
    if mod(i, 2) == 1  % Odd iteration: Keep psi_1 constant and search for psi_2
        if i == 1
            % Full search for the first iteration
            [p_max, ind] = max(p(:));
            [row, col] = ind2sub(size_p, ind);
        else
            % Partial search by fixing psi_1 at col
            [p_max, row] = max(p(:, col));  % Search only in the fixed column
        end
        psi_10(i) = psi_1(col);
        psi_20(i) = psi_2(row);
        
    else  % Even iteration: Keep psi_2 constant and search for psi_1
        % Partial search by fixing psi_2 at row
        [p_max, col] = max(p(row, :));  % Search only in the fixed row
        psi_10(i) = psi_1(col);
        psi_20(i) = psi_2(row);
    end
    
    % Check threshold condition
    if p_max < p_theshold
        disp('Threshold reached')
        break
    end
    
    % Convert to alpha0 and s0 values for the current iteration
    s0 = R * cos((psi_20(i) - psi_10(i)) / 2);
    alpha0 = (psi_10(i) + psi_20(i)) / 2;
    
    % Subtract p_line from p
    p_line_theory = p_line_fun(alpha0, s0, PSI_1, PSI_2, R, L, tstart, tend, d, p_min);
    p = p - p_line_theory;
end
disp(['Number of lines used = ', num2str(i)])

% p_line plot
figure
h = pcolor(PSI_1,PSI_2,p_line_theory);
xlabel('\psi_1 [rad]')
ylabel('\psi_2 [rad]')
title('Theoretical Radon Transform of line p_{theory}')
colormap(gca,hot), colorbar
set(h, 'EdgeColor', 'none');

% p plot
figure
h = pcolor(PSI_1,PSI_2,p);
xlabel('\psi_1 [rad]')
ylabel('\psi_2 [rad]')
title('Final Error: p - \Sigma p_{line}')
set(h, 'EdgeColor', 'none');
colormap(gca,hot), colorbar

%% Lines + Plot
psi_1 = mod(psi_10 + 2*pi, 2*pi); % converts to [0,2pi] instead of [-pi,pi]
psi_2 = psi_20;

Psi_1_plot = psi_1;
Psi_2_plot = psi_2;

Start_Line_x = R*cos(Psi_1_plot);
Start_Line_y = R*sin(Psi_1_plot);
End_Line_x = R*cos(Psi_2_plot);
End_Line_y = R*sin(Psi_2_plot);

Plot_Lines_x = [Start_Line_x; End_Line_x];
Plot_Lines_y = [Start_Line_y; End_Line_y];

%circle
circle_theta = linspace(0,2*pi,1000);
x_circ = R*cos(circle_theta);
y_circ = R*sin(circle_theta);

% Plot string art
figure('Name','String Art','Position',[50 100 800 600])
plot(Plot_Lines_x,Plot_Lines_y,'Color',[0,0,0,transparency])
axis square
hold on
plot(x_circ,y_circ,'k')
set(gcf,'color','w');
title(['Nails used = ',num2str(Num_Nails)])
set(gca,'XTick',[], 'YTick', [])

%% Generate .txt file
psi_1_number = round(psi_1 * Num_Nails/(2*pi));
psi_2_number = round(psi_2 * Num_Nails/(2*pi));

for i = 1:length(psi_1_number)
    if mod(i, 2) == 0 % even
        List(i) = psi_1_number(i);
    else
        List(i) = psi_2_number(i);
    end
end
List = List';
% writematrix(List, sprintf('Nail_Order_List_Num_Nails_%d.txt', Num_Nails));

WriteTextFile(filename, Num_Nails, List);

%% Functions
function [p,s] = radon_fun(f,alpha,R,ideal_image_size)
    % Calculate radon transform of f. MATLAB assumes distance between
    % pixels is 1 so rescale.
    [p,s] = radon(f,alpha);
    p = p * R*2/ideal_image_size;
    s = s * R*2/ideal_image_size;
end

function p_line = p_line_fun(alpha0,s0,PSI_1,PSI_2,R,L,tstart,tend,d,p_min)
    % Calculates radon transform of line p_line(alpha,s) using theoretical
    % model
    ALPHA = (PSI_1+PSI_2)/2;
    S = R*cos((PSI_2-PSI_1)/2);
    mask = maskfun(alpha0,s0,ALPHA,S,R,tstart,tend);
    p_line = d*p_min ./((d.*L-p_min).*abs(sin(ALPHA-alpha0)) + p_min);
    p_line = p_line.*mask;
end

function mask = maskfun(alpha0,s0,ALPHA,S,R,tstart,tend)
    % Creates a mask for calculating p_line. 1 within p_region, 0 outside.
    % Also fades from tstart to tend
    mask = zeros(size(ALPHA));
    mask(p_regionfun(alpha0,s0,ALPHA,S,R,tstart) < 0) = 1;
    n = 4;
    t = linspace(tstart,tend,n);
    for i = 1:n-1
        mask(p_regionfun(alpha0,s0,ALPHA,S,R,t(i)) > 0 &  p_regionfun(alpha0,s0,ALPHA,S,R,t(i+1)) < 0) = (tend-t(i))/(tend-tstart);
    end
end

function p_region = p_regionfun(alpha0,s0,ALPHA,S,R,t)
    % Used for generating mask for p_line
    % p_region > 0 is where p_line = 0
    p_region = (S.^2+s0.^2-2.*S.*s0.*cos(ALPHA-alpha0)) ./ ((sin(ALPHA - alpha0)).^2 + (t/(2*R))^2) - R^2;
end

function p_q = AlphaS2Phi(ALPHA,S,PSI_1,PSI_2,p,R)
    % Convert to psi
    PSI_1_convert = ALPHA - acos(S/R);
    PSI_2_convert = ALPHA + acos(S/R);
    
    % Perform 2D interpolation on p values using griddata
    p_q = griddata(PSI_1_convert(:), PSI_2_convert(:), p(:), PSI_1, PSI_2, 'linear');  % Linear interpolation
end

function WriteTextFile(filename, Num_Nails, Lines_List)
    % WriteTextFile creates a .txt file with:
    %   - A descriptive preamble
    %   - A list of lines (nail indices)
    %
    % The name of the .txt file is automatically generated as:
    %   CircleStringArtNailList_<imageFileNameWithoutExt>_<Num_Nails>.txt
    %
    % Parameters:
    %   filename   - The name/path of the image file, e.g. 'my_image.png'
    %   Num_Nails  - Total number of nails around the circumference
    %   Lines_List - Nx1 vector of nail indices representing the lines

    % Extract the file name (without path and extension)
    [~, imageFileNameWithoutExt, ~] = fileparts(filename);
    
    % Construct the output file name
    outFileName = sprintf('CircleStringArtNailList_%s_%d.txt', ...
                          imageFileNameWithoutExt, Num_Nails);

    % Open the file for writing
    fid = fopen(outFileName, 'w');
    if fid == -1
        error('Could not open file %s for writing.', outFileName);
    end
    
    % Preamble
    fprintf(fid, 'Image file is: filename = %s\n', filename);
    fprintf(fid, 'Number of lines of string is: Lines_Drawn = %d\n', length(Lines_List));
    fprintf(fid, 'Number of nails is: Num_Nails = %d\n', Num_Nails);
    fprintf(fid, ...
        'It''s up to you what radius you choose for the Canvas. For default parameters, I''d recommend R = 0.3 meters.\n\n');

    % Instructions
    fprintf(fid, ...
        'Instructions: Insert the nails (equally spaced) around the circumference of the circle.\n');
    fprintf(fid, ...
        'Put the nails around Counter Clockwise and label them 1 to %d.\n\n', Num_Nails);

    % Heading before the nail list
    fprintf(fid, 'The following nail list is:\n');
    
    % Print the nail list
    for i = 1:length(Lines_List)
        fprintf(fid, '%d\n', Lines_List(i));
    end

    % Close the file
    fclose(fid);

    % Display a message indicating where the file was saved
    fprintf('Nail list successfully written to: %s\n', outFileName);
end
