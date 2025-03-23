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
filename = 'Euler.png'
Num_Nails = 450;                % Enter the number of nails around the Canvas (must be even number)

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
num_max_lines = 12000;          % Algorithm will halt if number of lines > num_max_lines

% Plotting
transparency = 0.06;            % Each line darkness. Just for visualizng plot. If you are using thick string, increase this value to match how dark it appears 

%% Image f
A = imread(filename);
BW_raw = rgb2gray(A); % 255 = white, 0 = black
[origH, origW] = size(BW_raw);
ratio = origW / origH;
b = 2; 
h = b / ratio;   % initial h
[ratio, spacing, N_h, N_b] = Sanity_Check(Num_Nails, b, h);
h = b / ratio;   % updated h

finalW = origW;
finalH = round(finalW / ratio);
if finalH > origH
    finalH = origH;
    finalW = round(finalH * ratio);
end
if finalW > origW
    finalW = origW;
    finalH = round(finalW / ratio);
end
if (finalH > origH) || (finalW > origW)
    error('Desired crop size (%d x %d) exceeds original image size (%d x %d).', ...
          finalH, finalW, origH, origW);
end
% crop image
rowStart = floor((origH - finalH) / 2) + 1;
rowEnd   = rowStart + finalH - 1;
colStart = floor((origW - finalW) / 2) + 1;
colEnd   = colStart + finalW - 1;
A = A(rowStart:rowEnd, colStart:colEnd, :);
A = imresize(A,[round(ideal_image_size / ratio), ideal_image_size]);
BW = 1-double(rgb2gray(A))/255;
f = BW;


%% Plot CMYK components of image
width = 800;
height = width / ratio;
figure('Name', 'String Art','color', 'w','Position', [50 10 width height]);
imagesc(f)
title('Image f')
colormap(flipud(gray(256)));
axis off;
box off;


%% Calculations for CMYK

% Cyan
[Plot_Lines_x, Plot_Lines_y, Lines_list] = StringArtCalc(f,ideal_image_size,Num_Nails,d,p_min,tstart,tend,p_theshold,num_max_lines,spacing,b,h);

%% Plot
width = 800;
height = width / ratio;
figure('Name', 'String Art','color', 'w','Position', [50 10 width height]);
plot(Plot_Lines_x,Plot_Lines_y,'Color',[0,0,0,transparency])
xlim([-b/2,b/2])
ylim([-h/2,h/2])
axis off;
box off;

%% Write .txt file
WriteTextFile(filename, Num_Nails, N_h, N_b, ratio, Lines_list)

%% Functions
function [Plot_Lines_x, Plot_Lines_y,Lines_list] = StringArtCalc(f,ideal_image_size,Num_Nails,d,p_min,tstart,tend,p_theshold,num_max_lines,spacing,b,h)
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

        % Convert (alpha,s) to (L1,L2)
        [L1_local,L2_local] = AlphaS2L1L2_local(alpha0(i),s0(i),b,h);
        L1_convert(i) = L1_local;
        L2_convert(i) = L2_local;

        % Subtract p_line from p
        p_line_theory = p_line_fun(alpha0(i),s0(i),ALPHA,S,tstart,tend,d,p_min,b,h,L);
        p = p - p_line_theory;
    end
    Num_Lines_Used = i;
    % disp(['Number of lines used = ', num2str(Num_Lines_Used)])  % Only for continuous case

    %% Convert s,alpha --> L1,L2
    [L1, L2] = meshgrid(0:spacing:(2*b+2*h),0:spacing:(2*b+2*h));
    % Create matrix p. p = 1 where line exists and p = 0 else
    p = zeros(size(L1));
    for k = 1 : length(L1_convert)
        distMat = hypot(L1 - L1_convert(k), L2 - L2_convert(k));
        [~, idxMin] = min(distMat(:));
        p(idxMin) = 1;
    end

    % Note L starts at bottom right of rectangular canvas and goes CCW
    % Lines denotes index of L
    [Lines, p] = Interp2LinesCalc(p,Num_Lines_Used);
    Lines(Lines == Num_Nails + 1) = 1; % Wrap around since last nail = first nail

    % Plot
    [Plot_Lines_x, Plot_Lines_y, Lines] = Plot_Lines_Rectangle_Calc(Lines,spacing,b,h);

    % Create .txt file. Lines constains is Nx2 and Lines_list is Nx1.
    % Contain same info but Lines_list is more userfriendly
    Num_Lines_Drawn = length(Lines);    
    Lines_list = zeros(Num_Lines_Drawn,1);
    for i = 1:Num_Lines_Drawn
        Lines_list(i) = Lines(i,1+mod(i,2));
    end

    Num_Lines_Used_Discrete = length(Lines_list);
    disp(['Number of lines used = ', num2str(Num_Lines_Used_Discrete)])



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

function [r,spacing,N_h,N_b] = Sanity_Check(N,b,h)
    % N = total number of nails
    % r = b/h that ensures equal spacing
    % spacing = distance between adjacent nails
    % N_h = number of nails on vertical (h) side (including ends)
    % N_b = number of nails on horizontal (b) side (including ends)
    % N = 2(N_h + N_b) - 4
    % spacing = h/(N_h-1) = b/(N_b-1)

    % Check if N/2 is an integer
    if mod(N, 2) ~= 0
        N = N + 1;
        disp(['Nail Number must be even, so made Num_Nails = ', num2str(N)]);
    end
    
    B = (1:1:N/2 - 1)';
    H = 0.5*N - B;
    
    % User defined ratio r0 (might not correspond to integer nails on each
    % side)
    r0 = b/h;
    B0 = r0*N/(2+2*r0);
    H0 = N/(2+2*r0);
    
    % Find closest integer solution
    D = (B-B0).^2 + (H-H0).^2;
    [~,ind] = min(D);
    B = B(ind);
    H = H(ind);
    r = B/H;
    
    N_h = H + 1;
    N_b = B + 1;

    h = b/r;
    spacing = h/(N_h-1);
    
    disp(['Ratio of canvas needed to ensure equal spacing is: r = b/h = ', num2str(r)]);
    % disp(['Spacing between nails should be: spacing = ', num2str(spacing)]);    % This is only if b = 2
    disp(['Number of nails on side (including ends) should be: N_h = ', num2str(N_h)]);
    disp(['Number of nails on top (including ends) should be: N_b = ', num2str(N_b)]);
end

function [L1_local,L2_local] = AlphaS2L1L2_local(alpha_local,s_local,b,h)
        side_count = 0;

        % Initialize
        L1_local = NaN;
        L2_local = NaN;

        % Right side
        y = (s_local - 0.5*b*cos(alpha_local))./sin(alpha_local);
        if -0.5*h <= y && y <= 0.5*h
            side_count = side_count + 1;
            if side_count == 1
                L1_local = h/2 + y;
            elseif side_count == 2
                L2_local = h/2 + y;
            end
        end

        % Top side
        x = (s_local - 0.5*h*sin(alpha_local))./cos(alpha_local);
        if -0.5*b <= x && x <= 0.5*b
            side_count = side_count + 1;
            if side_count == 1
                L1_local = h + (0.5*b -x);
            elseif side_count == 2
                L2_local = h + (0.5*b -x);
            end
        end

        % Left side
        y = (s_local + 0.5*b*cos(alpha_local))./sin(alpha_local);
        if -0.5*h <= y && y <= 0.5*h
            side_count = side_count + 1;
            if side_count == 1
                L1_local = h+b + (0.5*h-y);
            elseif side_count == 2
                L2_local = h+b + (0.5*h-y);
            end
        end

        % Bot side
        x = (s_local + 0.5*h*sin(alpha_local))./cos(alpha_local);
        if -0.5*b <= x && x <= 0.5*b
            side_count = side_count + 1;
            if side_count == 1
                L1_local = 2*h+b + (0.5*b+x);
            elseif side_count == 2
                L2_local = 2*h+b + (0.5*b+x);
            end
        end
end



function [Lines, Grid] = Interp2LinesCalc(Grid, Max_Lines)
    % Lines = Coordinates of starting index and ending index of each line (L1,L2). Stored as Nx2 
    % Grid = Final grid after lines drawn (should ideally be close to zero matrix)
    
    Halt_Threshold = 2;     % Specify max number of times no 1 is found in a row before halting
    Lines = [];             % Will store moves as Nx2
    Line_Num = 1;
    
    % ---- FIRST MOVE (STARTING MOVE) ----
    c = 1;
    r = find(Grid(:,c) == 1, 1, 'first');
    if isempty(r)
        r = 1; 
    end
    
    %% CHANGED/NEW: If the first pick is (1,1), pick a different row
    if r == c
        rowCounts = Grid(:, c);          % how many 1's in each row (or just the row's value)
        rowCounts(c) = -Inf;            % exclude the diagonal index
        [~, bestRow] = max(rowCounts);
        r = bestRow;
    end
    
    Lines(Line_Num,:) = [r, c];
    Grid(r,c) = Grid(r,c) - 1;
    Line_Num = Line_Num + 1;
    currentPos = [r, c];
    
    % Keep track of consecutive times no '1' is found
    noOneFoundCount = 0;
    
    % ---- SUBSEQUENT MOVES ----
    % Even moves: horizontal (h), Odd moves: vertical (v)
    for i = 2 : Max_Lines
        if mod(i, 2) == 0
            direction = 'h';  % even move number -> horizontal
        else
            direction = 'v';  % odd move number -> vertical
        end
        
        [currentPos, foundOne] = findNextPosition(Grid, currentPos, direction);
        Lines(Line_Num,:) = currentPos;
        Grid(currentPos(1), currentPos(2)) = Grid(currentPos(1), currentPos(2)) - 1;
        Line_Num = Line_Num + 1;
        
        if ~foundOne
            noOneFoundCount = noOneFoundCount + 1;
        else
            noOneFoundCount = 0;
        end
        
        % If we fail to find a '1' twice in a row, halt
        if noOneFoundCount >= Halt_Threshold
            disp('Halt loop because no 1s found too many times in a row')
            break;
        end
    end
end


function [newPos, foundOne] = findNextPosition(Grid, currentPos, direction)
    r = currentPos(1);
    c = currentPos(2);

    switch direction
        case 'h'
            % 1) Find all columns that have a '1' in the current row r
            nextC = find(Grid(r, :) == 1);

            %% CHANGED/NEW: Remove columns == r to avoid picking (r,r)
            nextC = nextC(nextC ~= r);

            if isempty(nextC)
                % -> No '1' in the current row, so do the "max-ones" logic globally
                foundOne = false;

                % Sum(Grid,1) is a 1xN vector giving number of ones in each column
                colCounts = sum(Grid, 1);

                %% CHANGED/NEW: Exclude column == r
                if r <= length(colCounts)
                    colCounts(r) = -Inf;
                end

                [~, bestCol] = max(colCounts);
                nextC = bestCol;  
            else
                % -> We did find columns with a '1', so pick the column
                %    that leads to the largest sum of ones (in that column).
                foundOne = true;
                
                % Among the columns that have a 1 in row r, pick the column
                % with the largest total sum of ones in that column.
                [~, idx] = max(sum(Grid(:, nextC), 1));  
                nextC = nextC(idx);
            end
            newPos = [r, nextC];

        case 'v'
            % 1) Find all rows that have a '1' in the current column c
            nextR = find(Grid(:, c) == 1);

            %% CHANGED/NEW: Remove rows == c to avoid picking (c,c)
            nextR = nextR(nextR ~= c);

            if isempty(nextR)
                % -> No '1' in the current column, so do the "max-ones" logic globally
                foundOne = false;

                % Sum(Grid,2) is an Nx1 vector giving number of ones in each row
                rowCounts = sum(Grid, 2);

                %% CHANGED/NEW: Exclude row == c
                if c <= length(rowCounts)
                    rowCounts(c) = -Inf;
                end

                [~, bestRow] = max(rowCounts);
                nextR = bestRow; 
            else
                % -> We did find rows with a '1', so pick the row
                %    that leads to the largest sum of ones in that row.
                foundOne = true;
                
                % Among the rows that have a 1 in column c, pick the row
                % with the largest total sum of ones in that row.
                [~, idx] = max(sum(Grid(nextR, :), 2));
                nextR = nextR(idx);
            end
            newPos = [nextR, c];
    end
end
  

function [Plot_Lines_x, Plot_Lines_y, Lines] = Plot_Lines_Rectangle_Calc(Lines,spacing,b,h)
    
    % Lines = index, L = value
    L = (Lines - 1) * spacing;

    % Preallocation for speed
    num_lines = length(Lines);
    Start_Line_x = zeros(1, num_lines);
    Start_Line_y = zeros(1, num_lines);
    End_Line_x = zeros(1, num_lines);
    End_Line_y = zeros(1, num_lines);

    % Calculate (x,y) of start and end of line along rectangle edge
    % indx_remove = zeros(num_lines,1); % Initialize
    for i = 1:num_lines
        L1 = L(i,2);
        L2 = L(i,1);
        [x1,y1] = L2xy(L1,b,h);
        [x2,y2] = L2xy(L2,b,h);

        % % Store index of lines that go along edge (actaully better to
        % keep this so one string can be used.
        % tol = 1e-6;  % Tollerance
        % if abs(x1 - b/2) < tol && abs(x2 - b/2) < tol
        %     % Right side
        %     indx_remove(i) = 1;
        % elseif abs(y1 - h/2) < tol && abs(y2 - h/2) < tol
        %     % Top side
        %     indx_remove(i) = 1;
        % elseif abs(x1 + b/2) < tol && abs(x2 + b/2) < tol
        %     % Left side
        %     indx_remove(i) = 1;
        % elseif abs(y1 + h/2) < tol && abs(y2 + h/2) < tol
        %     % Bottom side
        %     indx_remove(i) = 1;
        % end

        % Store all lines
        Start_Line_x(1,i) = x1;
        Start_Line_y(1,i) = y1;
        End_Line_x(1,i) = x2;
        End_Line_y(1,i) = y2;
    end

    % % Manually delete lines that go along edge
    % indx_remove = find(indx_remove);
    % Lines(indx_remove,:) = [];
    % Start_Line_x(indx_remove) = [];
    % Start_Line_y(indx_remove) = [];
    % End_Line_x(indx_remove) = [];
    % End_Line_y(indx_remove) = [];

    Plot_Lines_x = [Start_Line_x; End_Line_x];
    Plot_Lines_y = [Start_Line_y; End_Line_y];
end


function [x,y] = L2xy(L,b,h)
    if 0 <= L && L <= h
        % Right side
        x = b/2;
        y = L - h/2;
    elseif h <= L && L <= h+b
        % Top side
        x = h + b/2 - L;
        y = h/2;
    elseif h+b <= L && L <= 2*h+b
        % Left side
        x = - b/2;
        y = b + 1.5*h - L;
    else % 2*h+b <= L && L <= 2*h+2*b
        % Bot side
        x = L - 2*h - 1.5*b;
        y = - h/2;
    end
end

function WriteTextFile(filename, Num_Nails, N_h, N_b, r, Lines_list)
    % WriteTextFile creates a .txt file with:
    %   - A descriptive preamble
    %   - A list of lines for a particular Color string (Nx1 vector)
    %
    % The name of the .txt file is automatically generated as:
    %   ColorStringArtNailList_<Color>_<imageFileNameWithoutExt>_<Num_Nails>.txt
    %
    % Parameters:
    %   Lines_list - Nx1 vector of nail indices
    %   Color      - The color of the string, e.g. 'Cyan'
    %   imageFile  - The name/path of the image file, e.g. 'my_image.png'
    %   N_h        - Number of nails on the "height" side (including corners)
    %   N_b        - Number of nails on the "base" side (including corners)
    %   Num_Nails  - Total number of nails around the perimeter
    %   r          - Ratio b/h
    
    % Extract the file name (without path and extension) from imageFile
    %   e.g. 'my_image.png' -> 'my_image'
    [~, imageFileNameWithoutExt, ~] = fileparts(filename);
    
    % Construct the output file name
    %   e.g. ColorStringArtNailList_Cyan_my_image_24.txt
    outFileName = sprintf('RectangularColorStringArtNailList_%s_%d.txt', ...
                          imageFileNameWithoutExt, Num_Nails);

    % Open the file for writing
    fid = fopen(outFileName, 'w');
    if fid == -1
        error('Could not open file %s for writing.', outFileName);
    end
    
    % Preamble
    fprintf(fid, 'Image file is: filename = %s\n', filename);
    fprintf(fid, 'Number of lines of string is: Lines_Drawn = %d\n', length(Lines_list));
    fprintf(fid, 'Number of total nails is: Num_Nails = %d\n', Num_Nails);
    fprintf(fid, 'Number of nails on the side (including ends) is: N_h = %d\n', N_h);
    fprintf(fid, 'Number of nails on the top (including ends) is: N_b = %d\n', N_b);
    fprintf(fid, ...
        'Ratio (b=base to h=height) of canvas needed to ensure equal spacing between nails is: r = b/h = %.3g\n', ...
        r);
    fprintf(fid, ...
        'The spacing between nails is calculated by the formula: spacing = h/(N_h-1) = b/(N_b-1)\n');
    fprintf(fid, ...
        'It''s up to you what base width, b, you choose. For default parameters, I''d recommend b = 0.6 meters.\n\n');
    
    % Instructions
    fprintf(fid, 'Instructions: Insert the nails (equally spaced) around the perimeter of the rectangle.\n');
    fprintf(fid, 'The first nail should start in the bottom right corner and go around Counter Clockwise.\n');
    fprintf(fid, 'Label these nails 1 to %d. \n\n', Num_Nails);
    
    % Heading before the nail list
    fprintf(fid, 'The following nail list is\n');

    % Print the nail list (assuming Nx1)
    for i = 1:length(Lines_list)
        fprintf(fid, '%d\n', Lines_list(i));
    end

    % Close the file
    fclose(fid);

    % Display a message indicating where the file was saved
    fprintf('Nail list successfully written to: %s\n', outFileName);
end