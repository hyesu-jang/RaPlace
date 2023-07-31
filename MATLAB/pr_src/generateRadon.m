function [sinoffts, rowkeys, xy_poses] = generateRadon(data_dir,down_shape)

%%
radar_data_dir = fullfile(data_dir, 'sensor_data/radar/backward_cart/');
data_names = osdir(radar_data_dir);

%% gps to xyz
% for mulran dataset
gtpose = csvread(strcat(data_dir, 'global_pose.csv'));
gtpose_time = gtpose(:, 1);
gtpose_xy = gtpose(:, [5,9]);

% for the oxford dataset
% radar_data_dir = fullfile(data_dir, 'backward_cart/');
% data_names = osdir(radar_data_dir);
% gtpose = csvread(strcat(data_dir, 'gps/gps_edit.csv'));
% gtpose_time = gtpose(:, 1);
% gtpose_xy = gtpose(:, [9,10]);

% figure(1); hold on;
% plot(traj_x, traj_y);

%%
num_data = length(data_names);
 theta=0:179;
%theta=0:359;
% sinograms = cell(1, num_data);
sinoffts = cell(1, num_data);
rowkeys = zeros(num_data, down_shape *180);
xy_poses = zeros(num_data, 2);
tradon = tic;
radon_time = zeros(1,num_data);
for data_idx = 1:num_data
    file_name = data_names{data_idx};
    data_time = str2double(file_name(1:end-4));
    data_path = strcat(radar_data_dir, file_name);
    
    % get
    tmp = imread(data_path);
    [R,xp] = radon(tmp,theta);
    R = R/max(max(R));
    R = double(R);
    R = imresize(R,down_shape);
    sinofft = abs(fft(R));
    sinofft_rows = size(sinofft,1);
    sinofft = sinofft(1:fix(sinofft_rows/2),:);
    
    [nearest_time_gap, nearest_idx] = min(abs(repmat(data_time, length(gtpose_time), 1) - gtpose_time));
    xy_pose = gtpose_xy(nearest_idx, :);
    radon_time(data_idx) = toc(tradon);
    % save 
%     sinograms{data_idx} = R;
    sinoffts{data_idx} = sinofft;
    rowkeys(data_idx, :) = rowkey(R);
    xy_poses(data_idx, :) = xy_pose;
   
    % log
    if(rem(data_idx, 100) == 0)
        message = strcat(num2str(data_idx), " / ", num2str(num_data), " processed.");
        disp(message); 
    end
end

end

