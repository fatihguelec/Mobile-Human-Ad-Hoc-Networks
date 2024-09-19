function  N_inf = infection2(no_of_nodes, no_of_steps, d_inf, Nx, Ny, N_inf)
% Trace the infection state of the nodes
for k = 1:no_of_steps 
    for ix = 1:no_of_nodes
        if N_inf(ix,k) == 1
            for i = 1:no_of_nodes
                r_d = sqrt( (Nx(ix,k)-Nx(i,k))^2 + (Ny(ix,k)-Ny(i,k))^2 ); %Euclidean distance between infected node 1 and others
                if r_d <= d_inf && i ~= ix
                    N_inf(i,k:end) = 1;
                end
            end
        end
    end
end



end