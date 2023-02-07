%**************************************************************************
% Update Log -- editted by Z.Wan  
% V2 Detect whether the ROI is a spot.
% No contrast is needed.
% V3 Auto choose the turning points of the oscillation.
% V4 Auto ROIs choosing. Does not scan the whole frame & auto choose the
% possible oscillators by evaluate the correlation coefficient at linear
% part and the standard deviation in plateau part.
% V5 Keep all the points in res even they are unused in calculation.
% V6 Change the order of evaluating the oscillation curves. Use the slope
% ratio between linear regime and plateau to distinguish the turning point.
% Then use the standard deviation of plateau and linearity to select the
% good one.
% V62 Fine tune thresholds are added for multipeaks selection. Larger signals
% have larger numerical diversion from their mean value. Use the different
% parameters can have a better dynamic detection.

%% Load the images and get the image number and their size 
clear
clc

file_path = 'D:\Data\20201228_2 single IgG (anti-IgA) 10K\res31-40\fft39';
file_list = dir(file_path);
file_list = file_list(3:end);
count = length(file_list);
image0 = imread([file_path, '\', file_list(end).name]);
[M, N] = size(image0);
% Delete a rectangle area of degraded image part.
% image0(32:192, 64:192) = 0;    
% image0(1:140, 1:256) = 0;

%% All pre-defined parameters used in following sections
% Spots selection parameters
bgpct_thres = 0.75;
max_pt = max(image0, [], 'all');
[N0, ~] = histcounts(image0, max_pt + 1);

% Oscillation curve selection parameters
% vlt = [100 150 200 250 300 350 400 450 500];    % Voltage applied mV
vlt = [0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5];    % Voltage applied V
% vlt = wfm(1,:);
s = 1;                              % Starting point of the images
pp = 3;                             % Earliest plateau point
is_fine = true;
fine_thres = 8;                    % Threshold for fine tuning

%************************************
corr_thres = 0.85;                  % Correlation coefficient, linearity
fb_thres = 8;                       % Slope ratio value
std_pct = 0.1;                     % Percentage of plateau variation 
%************************************
corr_fine = 0.88;                   % Fine tuning parameters
fb_fine = 16;
std_pct_fine = 0.05;
%************************************

gf_thres = 0.88;                    % Spot detection by gaussian filter
I_high_thres = 40;                  % Upper intensity limit
I_low_thres = 12;                    % Lower intensity

% ROI size and mask filters
rs = 7;                             % ROI_size
rr = floor(rs/2);                   % ROI_radius

% Round shape roi initialization
mask = [0, 0, 1, 1, 1, 0, 0;0, 1, 1, 1, 1, 1, 0;1, 1, 1, 1, 1, 1, 1;
        1, 1, 1, 1, 1, 1, 1;1, 1, 1, 1, 1, 1, 1;0, 1, 1, 1, 1, 1, 0;
        0, 0, 1, 1, 1, 0, 0];
mask_area = sum(mask, 'all');
%Gaussian filter, sgima = 1
gf = [
0.000036	0.000363	0.001446	0.002291	0.001446	0.000363	0.000036;
0.000363	0.003676	0.014662	0.023226	0.014662	0.003676	0.000363;
0.001446	0.014662	0.058488	0.092651	0.058488	0.014662	0.001446;
0.002291	0.023226	0.092651	0.146768	0.092651	0.023226	0.002291;
0.001446	0.014662	0.058488	0.092651	0.058488	0.014662	0.001446;
0.000363	0.003676	0.014662	0.023226	0.014662	0.003676	0.000363;
0.000036	0.000363	0.001446	0.002291	0.001446	0.000363	0.000036];
    
%% Find possible oscillator ROIs
% Dynamic threshold for binarization
bin_thres = 0;
temp_sum = N0(1);
i = 1;
while(temp_sum/(M*N) < bgpct_thres)
    bin_thres = i;
    temp_sum = temp_sum + N0(i+1);
    i = i + 1;
end
disp(i);

% bin_thres = 5;            
bin_mask = image0;                              % Binarization
bin_mask(bin_mask <= bin_thres) = 0;
bin_mask(bin_mask > bin_thres) = 1;
nhood = ones(3);
bin_mask = imopen(bin_mask, nhood);             % Opening

[label, num_a] = bwlabel(bin_mask,4);           % Find connected area
sub_img = cell(1, num_a);
location = zeros(2, 2, num_a);
corner = zeros(2);
for i = 1 : num_a
    [row, col] = find(label == i);
    corner(1,1) = min(row);
    corner(1,2) = min(col);
    corner(2,1) = max(row);
    corner(2,2) = max(col);
    sub_img{i} = label(corner(1,1): corner(2,1), corner(1,2): corner(2,2));
    location(:,:,i) = corner;
end

c_row = [];
c_col = [];
rest_row = [];
rest_col = [];
subtotal = 0;
temp_step = 4;

for i = 1 : num_a
    temp = sub_img{i}';
    [temp_m, temp_n] = size(temp);
    if temp_m*temp_n <= rs*rs
        c_row(i) = round((location(1,1,i) + location(2,1,i))/2);
        c_col(i) = round((location(1,2,i) + location(2,2,i))/2);
    else
        c_row(i) = NaN;
        c_col(i) = NaN;
        seq_row = location(1,1,i): temp_step: location(2,1,i);
        seq_col = location(1,2,i): temp_step: location(2,2,i);
        n_seq_row = length(seq_row);
        n_seq_col = length(seq_col);
        temp_rest = zeros(n_seq_row * n_seq_col, 2);
        for p1 = 1 : n_seq_row
            for q1 = 1 : n_seq_col
                if seq_row(p1) > M - rr
                    seq_row(p1) = M - rr;
                end
                if seq_row(p1) > N - rr
                    seq_row(p1) = N - rr;
                end
                if seq_row(p1) <= rr
                    seq_row(p1) = rr + 1;
                end
                if seq_row(p1) <= rr
                    seq_row(p1) = rr + 1;
                end
                if seq_col(q1) > M - rr
                    seq_col(q1) = M - rr;
                end
                if seq_col(q1) > N - rr
                    seq_col(q1) = N - rr;
                end
                if seq_col(q1) <= rr
                    seq_col(q1) = rr + 1;
                end
                if seq_col(q1) <= rr
                    seq_col(q1) = rr + 1;
                end
                temp_rest( (p1-1)*n_seq_col + q1, 1) = seq_row(p1);
                temp_rest( (p1-1)*n_seq_col + q1, 2) = seq_col(q1);
            end
        end
        for j = 1 : size(temp_rest,1)
            roi = label(temp_rest(j,1) - rr: temp_rest(j,1) + rr, ...
                temp_rest(j,2) - rr: temp_rest(j,2) + rr );
            roi(roi ~= i) = 0;
            roi_area = sum(roi, 'all')/i;           
            if roi_area < sum(mask_area)/2
                temp_rest(j,:) = NaN;
            end
        end
        temp_rest = (rmmissing(temp_rest));
        rest_row = [rest_row; temp_rest(:,1)];
        rest_col = [rest_col; temp_rest(:,2)]; 
    end
end

c_row = [c_row, rest_row'];
c_col = [c_col, rest_col'];
c_row = rmmissing(c_row);
c_col = rmmissing(c_col);

%% Get the intensity sequence 
na = size(c_row,2);             % number of ROIs
res = zeros(na, count + 3);     % 3 more columns for turning point, slope and plateau
% meanv = zeros(1,count);         % Mean value of imgs

% Intensity summary
for i = 1 : count
    img = double(imread([file_path, '\', file_list(i).name]));
%     meanv(i) = mean(img,'all');
    for k = 1 : na              % Possible area index
        
        if c_row(k) > M - rr
           c_row(k) = M - rr;  
        end
        if c_col(k) > N - rr
           c_col(k) = N - rr; 
        end
        if c_row(k) <= rr
           c_row(k) = rr + 1;  
        end
        if c_col(k) <= rr
           c_col(k) = rr + 1; 
        end
        tempi = img(c_row(k)-rr:c_row(k)+rr, c_col(k)-rr: c_col(k)+rr) .* mask;
        
        % Center check: sigma = 1 Gaussian filter
        if i == count
            bf = tempi;
            bf = bf .* gf;        % Buffer
            bf1 = bf(2,4)+ sum(bf(3,3:5)) + sum(bf(4,2:6)) + sum(bf(5,3:5)) + bf(6,4);
            bf2 = sum(bf, 'all');
            if bf1/bf2 < gf_thres
                res(k,:) = NaN;
                c_row(k) = NaN;
                c_col(k) = NaN;
                continue
            end
        end
        
        res(k, i) = sum(tempi, 'all');      
    end
end

res = rmmissing(res);
c_row = rmmissing(c_row);
c_col = rmmissing(c_col);
% res = res(:,start_pos: end);
% vlt = vlt(start_pos: end);
pts = length(vlt);
 
%% Oscillation criteria
% 1. First point should be in low voltage region
osc_n = size(res, 1);        % number of oscillators 
for i = 1 : osc_n
    [~, rmin_pos] = min(res(i, s:pts) ); 
    if rmin_pos > s + 1
        res(i,:) = NaN;
        c_row(i) = NaN;
        c_col(i) = NaN;
    end
end
res = rmmissing(res);
c_row = rmmissing(c_row);
c_col = rmmissing(c_col);
res(:, 1:pts) = res(:, 1:pts) / sum(mask,'all');

% 2. Choose the best linearity & std
osc_n = size(res, 1);           % number of oscillators
for i = 1: osc_n
    front = zeros(1, pts-pp-s);
    back = zeros(1, pts-pp-s);
    for j = s+pp-1 : pts-2
        temp_f = polyfit(vlt(s:j), res(i, s:j), 1);
        temp_b = polyfit(vlt(j:end), res(i, j:pts), 1);
        front(1, j-pp-s+2) = abs(temp_f(1));
        back(1, j-pp-s+2) = abs(temp_b(1));
    end
    
    fb_ratio  = front ./ back;
    [fb_max, tp] = max(fb_ratio);
    slope = front(tp);
    tp = tp + pp + s - 2;
    
    plateau = mean(res(i, tp:pts));
    std_thres = std_pct * plateau;
    
    r2_linear = corr2(res(i, s:tp), vlt(s:tp));
    std_plateau = std(res(i, tp:pts));
    
    if fb_max < fb_thres || r2_linear < corr_thres || std_plateau > std_thres || slope < 0
        res(i,:) = NaN;
        c_row(i) = NaN;
        c_col(i) = NaN;
    else
        res(i, end-2) = tp;  
        res(i, end-1) = slope;
        res(i, end) = plateau;
    end
end
res = rmmissing(res);
c_row = rmmissing(c_row);
c_col = rmmissing(c_col);

% 3. Threshold
osc_n = size(res, 1);           % number of oscillators
for i = 1: osc_n
    
    if res(i, end) > I_high_thres || res(i, end) < I_low_thres
        res(i,:) = NaN;
        c_row(i) = NaN;
        c_col(i) = NaN;
    end
end
res = rmmissing(res);
c_row = rmmissing(c_row);
c_col = rmmissing(c_col);

% 4. Fine tuning
if is_fine == true
    osc_n = size(res, 1);           % number of oscillators
    for i = 1: osc_n
        if res(i, end) < fine_thres 
            tp_f = res(i,end-2);
            fit_f = polyfit(vlt(s:tp_f), res(i,s:tp_f),1);
            fit_b = polyfit(vlt(tp_f:pts), res(i, tp_f:pts), 1);
            slope_f = fit_f(1);
            slope_b = fit_b(1);
            fb_f = abs(slope_f / slope_b);
            r2_f = corr2(vlt(s:tp_f), res(i, s:tp_f));
            mp_f = mean(res(i, tp_f:pts));          % mean_plateau_fine
            std_thres_f = mp_f * std_pct_fine;
            std_f  = std(res(i, tp_f:pts));
            if r2_f < corr_fine || std_f > std_thres_f || fb_f < fb_fine
                res(i,:) = NaN;
                c_row(i) = NaN;
                c_col(i) = NaN;
            end
        else
            continue
        end
    end
    res = rmmissing(res);
    c_row = rmmissing(c_row);
    c_col = rmmissing(c_col);
end

%% Highlight selected ROIs
osc_n = size(res,1);
figure(1)
imshow(image0, [0 max_pt])
% imshow(image0, [bin_thres max_pt])
for i = 1 :  osc_n 
    rectangle('Position', [c_col(i)-4,c_row(i)-4,7,7], 'Curvature', [1,1], 'Edgecolor', 'r');
end


figure(2)
imshow(bin_mask, [0 1])
for i = 1 :  osc_n 
    rectangle('Position', [c_col(i)-4,c_row(i)-4,7,7], 'Curvature', [1,1], 'Edgecolor', 'r');
end