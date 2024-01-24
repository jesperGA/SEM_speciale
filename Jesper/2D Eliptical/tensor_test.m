% Assuming the dimensions of the tensor A
I = 3;  % Example value for I
J = 3;  % Example value for J
M = 3;  % Example value for M
N = 3;  % Example value for N

% Initialize a cell array to store the strings
A_strings = cell(I*J, M*N);

% Nested for loops to populate A_strings
for i = 1:I
    for j = 1:J
        for m = 1:M
            for n = 1:N
                % Compute row and column indices for A_strings
                row = (i-1)*J + j;
                col = (m-1)*N + n;

                % Assign the string representing the indices to the cell array
                A_strings{row, col} = sprintf('A_{%d%d%d%d}', i, j, m, n);
            end
        end
    end
end

% Display the resulting array of strings
disp(A_strings);

% Initialize a cell array to store the strings for U
U_strings = cell(M*N, 1);

% Loop to populate U_strings
for m = 1:M
    for n = 1:N
        % Compute the index for U_strings
        idx = (m-1)*N + n;

        % Assign the string representing the indices to the cell array
        U_strings{idx} = sprintf('U_{%d%d}', m, n);
    end
end

% Display the resulting vector of strings
disp('U tensor (string representation):');
disp(U_strings);