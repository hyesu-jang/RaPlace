clear; clc;
%%%%
%%%% The base of this PR code was written by Giseop Kim, from Scan Context(T-RO 2021).
%%%% Base Codes can be found in "https://github.com/irapkaist/scancontext"
%%%%
%%%% This is the sample MATLAB code for operating RaPlace, so it might have some errors.
%%%% Please feel free to notice me if the error exists in the code.
%%%% Maintainer : Hyesu Jang (dortz@snu.ac.kr)
%%%
%%%

addpath(genpath('pr_src'));
addpath(genpath('../data'));

%% data preparation 
global data_path; 
data_path = '/somewhere/ICRA20_MULRAN/KAIST/20190823/';
disp(strcat('Processing :  ',data_path(end-12:end-10), data_path(end-4:end-1)," within ", num2str(criteriaaa(crit)), "m/"));
% ### NOTE: Use this sequence directory structure
% example:
% /your/MulRan/sequence/dir/Riverside02/
%   L sensor_data/ 
%       L radar/
%           L polar/
%               L {unix_times}.png
%   L global_pose.csv

down_shape = 0.1;
[data_sinofft,data_rowkeys, data_poses] = loadData(down_shape);

%% main - global recognizer
revisit_criteria = 20; % in meter (recommend test for 5, 10, 20 meters)
keyframe_gap = 1; % for_fast_eval (if 1, no skip)

global num_candidates; num_candidates = 10; %% elements of the hierarchical group
global num_node_enough_apart; num_node_enough_apart = 50; 

% policy (top N)
num_top_n = 25;
top_n = linspace(1, num_top_n, num_top_n);

% Entropy thresholds 
middle_thres = 0.001;
thresholds1 = linspace(0, middle_thres, 50); 
thresholds2 = linspace(middle_thres, 0.01, 50);
thresholds = [thresholds1, thresholds2];
num_thresholds = length(thresholds);

% Main variables to store the result for drawing PR curve 
num_hits = zeros(num_top_n, num_thresholds); 
num_false_alarms = zeros(num_top_n, num_thresholds); 
num_correct_rejections = zeros(num_top_n, num_thresholds); 
num_misses = zeros(num_top_n, num_thresholds);

% Check the threshold for max value
case_for_hit = zeros(length(data_poses),3);
case_for_fa = [];

% main 
loop_log = [];

exp_poses = [];
exp_rowkeys = [];
exp_sinofft = {};
exp_sinograms = {};

num_queries = length(data_poses);
for query_idx = 1:num_queries - 1
    query_sinofft = data_sinofft{query_idx};
    query_sinofft = (query_sinofft-mean(query_sinofft(:)))/std(query_sinofft(:));
    query_rowkey = data_rowkeys(query_idx,:);
    query_pose = data_poses(query_idx,:);

    exp_sinofft{end+1} = query_sinofft;
    exp_rowkeys = [exp_rowkeys; query_rowkey];
    exp_poses = [exp_poses; query_pose];
    
    if(rem(query_idx, keyframe_gap) ~= 0)
       continue;
    end

    if( length(exp_sinofft) < num_node_enough_apart )
       continue;
    end
    
    can_sinofft = exp_sinofft(1:end-(num_node_enough_apart-1));
    
%%%%%% Single process
    tmpval = 0;
    maxval = 0;
    rotval = 0;
    candnum = 0;
    for cands = 1:size(can_sinofft,2)
        tmp_sinofft = can_sinofft{1,cands};
          [fftresult,tmpval] = fast_dft(query_sinofft, tmp_sinofft);
        if (maxval < tmpval)
            maxval = tmpval;
            candnum = cands;
        end
    end
    nearest_idx = candnum;
    [fftresult,tmpval] = fast_dft(query_sinofft, query_sinofft);
    min_dist = (tmpval-maxval)/1000;
    real_dist = dist_btn_pose(query_pose, exp_poses(nearest_idx, :));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Hierarchical Process
%%%%%%  --We don't need this if the number of images are small
%%%%%%  --All the radar datasets can be operated in realtime without this
%     group_size = num_candidates;
%     groupnum = 0;
%     groups = ceil(size(can_sinofft,2) / group_size); % Determine the number of groups
%     for g = 1:groups
%         group_representative = can_sinofft{1,(g-1)*group_size+1};
%         [fftresult,tmpval] = fast_dft(query_sinofft, group_representative);
%         if (maxval < tmpval)
%             maxval = tmpval;
%             groupnum = g-1;
%         end
%     end
%     for cands = 1:group_size
%         tmp_sinofft = can_sinofft{1,groupnum*group_size+cands};
%         [fftresult,tmpval] = fast_dft(query_sinofft, tmp_sinofft);
%         if (maxval < tmpval)
%             maxval = tmpval;
%             candnum = groupnum*group_size+cands;
%         end
%     end
%     nearest_idx = candnum;
%     [fftresult,tmpval] = fast_dft(query_sinofft, query_sinofft);
%     min_dist = (tmpval-maxval)/1000;
%     real_dist = dist_btn_pose(query_pose, exp_poses(nearest_idx, :));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

%% prcurve analysis 
    for topk = 1:num_top_n
        for thres_idx = 1:num_thresholds
            threshold = thresholds(thres_idx);
            
            reject = 0;
            if( min_dist > threshold)
                reject = 1; 
            end

            if(reject == 1) 
                if(revisitness == 0)
                    % TN: Correct Rejection
                    num_correct_rejections(topk, thres_idx) = num_correct_rejections(topk, thres_idx) + 1;
                else            
                    % FN: MISS
                    num_misses(topk, thres_idx) = num_misses(topk, thres_idx) + 1; 
                end
            else
                % if under the theshold, it is considered seen.
                % and then check the correctness
                if( real_dist < revisit_criteria)
                    % TP: Hit
                    num_hits(topk, thres_idx) = num_hits(topk, thres_idx) + 1;
                    case_for_hit(query_idx,1) = 1;
                    case_for_hit(query_idx,2) = nearest_idx;
                    case_for_hit(query_idx,3) = dist_btn_pose(query_pose, exp_poses(nearest_idx, :));
                else
                    % FP: False Alarm 
                    num_false_alarms(topk, thres_idx) = num_false_alarms(topk, thres_idx) + 1;            
                end
            end
            
        end
    end

    if( rem(query_idx, 100) == 0)
        disp( strcat(num2str(query_idx/num_queries * 100), ' % processed') );
    end
    
end


%% save the log 
savePath = strcat("pr_result/pp ", data_path(end-12:end-10), data_path(end-4:end-1)," within ", num2str(revisit_criteria), "m/");
if((~7==exist(savePath,'dir')))
    mkdir(savePath);
end
save(strcat(savePath, 'nCorrectRejections.mat'), 'num_correct_rejections');
save(strcat(savePath, 'nMisses.mat'), 'num_misses');
save(strcat(savePath, 'nHits.mat'), 'num_hits');
save(strcat(savePath, 'nFalseAlarms.mat'), 'num_false_alarms');

hittingPath = strcat("/path_for_saving_pr/",data_path(end-12:end));
if((~7==exist(hittingPath,'dir')))
    mkdir(hittingPath);
end
save(strcat(hittingPath, num2str(revisit_criteria), "m ", 'pp hit.mat'), 'case_for_hit');


