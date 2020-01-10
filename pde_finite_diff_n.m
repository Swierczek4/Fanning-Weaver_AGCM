function num_diff = pde_finite_diff_n(U,variable,derivative,order,h)
%% Finite Difference calculation for array U
% Calculates the first or second derivatives in the x and y direction
% (x if you input 1 for variable, y if you input 2)
% Select order out of 2 or 4.
% 
% NOTE! Assumes / Enforces Neumann zero flux boundary conditions!
% https://en.wikipedia.org/wiki/Finite_difference_coefficient
%%

if variable==1
    n = size(U,2);
    if derivative==1
        if order==2
            num_diff = ([U(:,2:n),U(:,n-1)]-[U(:,2),U(:,1:n-1)])./(2*h);
        elseif order==4
            num_diff = (-0.25.*[U(:,3:n),U(:,(n-1):-1:(n-2))]+2.*[U(:,2:n),U(:,n-1)]-...
                2.*[U(:,2),U(:,1:n-1)]+0.25.*[U(:,3:-1:2),U(:,1:n-2)])./(3*h);
        end
    elseif derivative==2
        if order==2
            num_diff = ([U(:,2:n),U(:,n-1)]-2.*U+[U(:,2),U(:,1:n-1)])./(2*h.^2);
        elseif order==4
            num_diff = (-0.25.*[U(:,3:n),U(:,(n-1):-1:(n-2))]+4.*[U(:,2:n),U(:,n-1)]-7.5.*U+...
                4.*[U(:,2),U(:,1:n-1)]-0.25.*[U(:,3:-1:2),U(:,1:n-2)])./(3*h.^2);
        end
    end
elseif variable==2
    m = size(U,1);
    if derivative==1
        if order==2
            num_diff = ([U(2:m,:);U(m-1,:)]-[U(2,:);U(1:m-1,:)])./(2*h);
        elseif order==4
            num_diff = (-0.25.*[U(3:m,:);U((m-1):-1:(m-2),:)]+2.*[U(2:m,:);U(m-1,:)]-...
                2.*[U(2,:);U(1:m-1,:)]+0.25.*[U(3:-1:2,:);U(1:(m-2),:)])./(3*h);
        end        
    elseif derivative==2
        if order==2
            num_diff = ([U(2:m,:);U(m-1,:)]-2.*U+[U(2,:);U(1:m-1,:)])./(h.^2);
        elseif order==4
            num_diff = (-0.25.*[U(3:m,:);U((m-1):-1:(m-2),:)]+4.*[U(2:m,:);U(m-1,:)]-7.5.*U+...
                4.*[U(2,:);U(1:m-1,:)]+0.25.*[U(3:-1:2,:);U(1:(m-2),:)])./(3*h.^2);
        end
    end
end
end

