function [ row_key ] = rowkey(sino)

% row_key = sum(ting);
% row_key = normalize(row_key);
[row,col] = size(sino);
row_key = sino(fix((row+1)/2),:);


end

