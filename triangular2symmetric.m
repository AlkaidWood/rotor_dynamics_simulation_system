%% triangular2symmetric
% transfer a triangular matrix to symmetric matrix
%% Syntax
% B = triangular2symmetric(A)
%% Description
% A is a triangular matrix (upper or lower)
%
% B is a symmetric matrix with elements in A


function B = triangular2symmetric(A)

% Determine whether the matrix is upper triangular or lower triangular
isUpper = istriu(A);
isLower = istril(A);
if ~isUpper && ~isLower
    error('please input a triangular matrix')
end

%%

% transfer
diagElement = diag(A);
B = A - diag(diagElement) +A';

end