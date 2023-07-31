function [sinofft, rowkeys, poses] = loadData(down_shape)

%%
global data_path;
data_save_path = fullfile('somewhere',data_path(end-12:end)); 

%%
% newly make
is_already_made_data_exist = exist(data_save_path);
if is_already_made_data_exist == 0 
    % make 
    [sinofft,rowkeys, poses] = generateRadon(data_path,down_shape);    
    
    % save
    mkdir(data_save_path);

%     filename = strcat(data_save_path, 'sinograms', '.mat');
%     save(filename, 'sinograms');
    filename = strcat(data_save_path, 'sinofft', '.mat');
    save(filename, 'sinofft','-v7.3');
    filename = strcat(data_save_path, 'rowkeys', '.mat');
    save(filename, 'rowkeys');
    filename = strcat(data_save_path, 'poses', '.mat');
    save(filename, 'poses');

% or load 
else
%     filename = strcat(data_save_path, 'sinograms', '.mat');
%     load(filename);
%     % fix
%     for iii = 1:length(sinograms)
%         sinograms{iii} = double(sinograms{iii});
%     end
    
    filename = strcat(data_save_path, 'sinofft', '.mat');
    load(filename);
    filename = strcat(data_save_path, 'rowkeys', '.mat');
    load(filename);
    filename = strcat(data_save_path, 'poses', '.mat');
    load(filename);
    
    disp('- successfully loaded.');
end

%%
disp(' ');

end

