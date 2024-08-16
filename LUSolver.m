function [x] = LUSolver(A,b)
%LUSolver Solves system Ax=b using LU decomposition with partial pivoting
%   A is a square NxN matrix, b is Nx1 vector, x is Nx1 solution vector.
%   The function performs partial pivoting of A and then LU decomposition
%   to calculate x.

%% Error check and preparation

%  Check that A is a square matrix
if size(A,1) ~= size(A,2)
    error('Error: A must be a square matrix\n')
end % End square matrix check

%  Check that b is a vector and the same length as a dimension of A
if size(b,2) ~= 1
    error('Error: b must be a vector\n')
end % End vector check

if length(b) ~= size(A,1)
    error('Error: Length of b must be equal to N of NxN matrix A\n')
end % End length check

%  Check that determinant of A is not zero to have a unique solution
if det(A) == 0
    error('Error: Determinant of A is zero. No unique solution.\n')
end % End determinant check

%  Assign dimension of NxN matrix
N = size(A,1);

%  Preallocate NxN L and U matrices as well as Nx1 vectors y and x
L = zeros(N,N);
U = zeros(N,N);
y = zeros(N,1);
x = zeros(N,1);

%% Partial Pivoting

%  For each column
for col = 1:1:N
    %  Find the index of the maximum absolute value in the column
    [~, rowIndx] = max(abs(A(:,col)));

    %  Replace the row of the current submatrix with the row that has the
    %  largest absolute value in the column
    A([col rowIndx],:) = A([rowIndx col],:);

    %  Switch corresponding entries of vector b
    b([col rowIndx],:) = b([rowIndx col],:);

end % End pivoting

%% LU Decomposition

%  Assign diagonal of L to be filled with 1s
for i=1:1:N
    L(i,i) = 1;
end

%  Assign first row of U to be the first row of A
U(1,:) = A(1,:);
%  Compute the first column of L
L(:,1) = A(:,1)/U(1,1);

%  For the rest of indicies
for i=2:1:N
    %  From column i to N
    for j=i:1:N
        %  Compute U entries
        U(i,j) = A(i,j) - L(i,1:i-1)*U(1:i-1,j);
    end % End computation of U
    %  From row i+1 to N
    for k=i+1:1:N
        %  Compute L entries
        L(k,i) = (A(k,i) - L(k,1:i-1)*U(1:i-1,i))/U(i,i);
    end % End computation of L
end % End computation of L and U

%  Forward substitution
%  Compute first element of y
y(1) = b(1)/L(1,1);

%  For rest of y elements
for k=2:N
    %  Compue all elements of y vector
    y(k) = (b(k) - L(k,1:k-1)*y(1:k-1))/L(k,k);
end %  End y vector computation

%  Back substitution
%  Compute last element of vector x
x(N) = y(N)/U(N,N);

%  For rest of x elements
for k=N-1:-1:1
    %  Compute all elements of x
    x(k) = (y(k) - U(k,k+1:N)*x(k+1:N))/U(k,k);
end % End computation of vector x

end % End of LUSolver function