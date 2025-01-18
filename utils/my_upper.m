function [upper] = my_upper(matrix)
idx_upper = tril(true(size(matrix)), -1);
upper = matrix(idx_upper);
endma