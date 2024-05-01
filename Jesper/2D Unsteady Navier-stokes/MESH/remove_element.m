function [X, IX, newBoundary] = remove_element(X, IX, n_remove)
    %[X, IX, newBoundary] = remove_element(X, IX, n_remove)
    % Remove element #n_remove
    remove_dofs = IX(:,:,n_remove);
    IX = IX(:,:,[1:n_remove-1 n_remove+1:end]);
    newBoundary = intersect(remove_dofs(:), IX(:));
    remove_dofs = setdiff(remove_dofs(:), IX(:));

    % Create logical index for rows to keep
    rows_to_keep = true(length(X), 1);
    rows_to_keep(remove_dofs) = false;

    % Index the array using logical indexing
    X = X(rows_to_keep, :);
    mapping_from = X(:,1);
    mapping_to   = (1:length(X))';

    IX = replace_numbers(IX, mapping_from, mapping_to);
    X(:,1) = 1:length(X);

    % Create a mapping from replace_from to replace_to
    [~, idx] = ismember(newBoundary, mapping_from);
    % Replace elements using the mapping
    newBoundary = mapping_to(idx);    
    
end

function IX_new = replace_numbers(IX, mapping_from, mapping_to)
    % Convert 3D array to 1D array
    IX_1d = IX(:);
    
    % Create a mapping from replace_from to replace_to
    [~, idx] = ismember(IX_1d, mapping_from);
    
    % Replace elements using the mapping
    IX_new = mapping_to(idx);    

    % Reshape back to 3D array
    IX_new = reshape(IX_new, size(IX));
end