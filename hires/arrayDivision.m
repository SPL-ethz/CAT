function [Astar] = arrayDivision(A,dim,order)
% divides array A element by element along dimension dim, so that 
% arrayDivision(A,1) gives A(2,:,:)./A(1,:,:) etc. for a 3D array    

%shift the dimensions of the matrix in such a way that the dimension to
%be divised is the last one
if length(size(A))>=3

    A = shiftdim(A,dim);

    %calculate the indeces of elements to be divided in the different
    %layers. variable blub contains the pairs of indices to be divided in
    %two adjacent columns
    s = size(A);
    inc = prod(s(1:end-1));
    basis = linspace(1,inc,inc)';
    increments = (1:inc:numel(A))-1;
    blub = basis*ones(1,length(increments))+ones(inc,1)*increments;

    %do the division
    if order==2
        Astar = A(blub(:,2:end))./A(blub(:,1:end-1));
    elseif order==1
        Astar = A(blub(:,1:end-1))./A(blub(:,2:end));
    end

    %this is the new size of the divided matrix
    s(end) = s(end)-1;

    %reshaping
    Astar = reshape(Astar, s);

    %re-shifting to reflect the initial dimensions

    Astar = shiftdim(Astar,length(s)-dim);

elseif length(size(A))<3
    % do the same also for matrices and vectors
    A=shiftdim(A,dim-1);
    if order == 2
        Astar=A(2:end,:)./A(1:end-1,:);
    elseif order == 1
        Astar=A(1:end-1,:)./A(2:end,:);
    end

    Astar=shiftdim(Astar,dim-1);
end