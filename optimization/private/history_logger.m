%% history_logger
% Logs history.

function history_logger(filename, get_epsilon, index, x, z, p)

    N = length(x); % Number of modes.
    addnum = @(str, num) [str, num2str(num)];

    % Get the unvectorization functions.
    for k = 1 : N
        get_eps{k} = @(z) my_get_eps(get_epsilon{k}, z);
        eps = get_eps{k}(z);
        dims = size(eps);
        unvec{k} = @(u) reshape(u, dims);
    end

    if index == 1 % Assume file and datasets does not yet exist, so create it.

        % Delete the file if it currently exists.
        try
            delete(filename);
        end

        my_create_dataset(filename, 'z', length(z));
        my_create_dataset(filename, 'p', length(p));

        for k = 1 : N
            my_create_dataset(filename, addnum('x', k), length(x{k}));
            my_create_dataset(filename, addnum('E', k), dims);
            my_create_dataset(filename, addnum('eps', k), dims);
        end
    end

    % Write the datasets.
    my_write_dataset(filename, 'z', z, index);
    my_write_dataset(filename, 'p', p, index);

    for k = 1 : N
        my_write_dataset(filename, addnum('x', k), x{k}, index);
        my_write_dataset(filename, addnum('E', k), unvec{k}(x{k}), index);
        my_write_dataset(filename, addnum('eps', k), get_eps{k}(z), index);
    end
end

function my_create_dataset(filename, datasetname, dims)
% Creates datasets for complex arrays.
% Includes real, imaginary, and absolute value datasets.

    type = {'_real', '_imag', '_abs'};
    for k = 1 : length(type)
        h5create(filename, ['/', datasetname, type{k}], [dims Inf], ...
                    'ChunkSize', [dims 1]);
    end
end

function my_write_dataset(filename, datasetname, data, index)
% Writes to datasets for complex arrays.
% Includes real, imaginary, and absolute value datasets.

    dims = size(data);
    if length(dims) == 2 % Detect one-dimensional arrays.
        dims = dims(1);
    end

    type = {'_real', '_imag', '_abs'};
    f = {@real, @imag, @abs};
    for k = 1 : length(type)
        h5write(filename, ['/', datasetname, type{k}], f{k}(data), ...
                [ones(1, length(dims)), index], [dims, 1]);
    end
end


function [eps] = my_get_eps(get_epsilon, z)
% Gets epsilon and puts it into a 4D array.
    epsilon = get_epsilon(z);
    eps = cat(4, epsilon{:});
end
