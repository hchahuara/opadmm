function params = read_params_table(filename)

    opts = spreadsheetImportOptions('NumVariables',2,'DataRange','A:B');
    data = readcell(filename,opts);

    params = struct();

    for i = 2:height(data)
        param_name = data{i,1};
        param_val  = data{i,2};
        params.(param_name) = param_val;
    end

end