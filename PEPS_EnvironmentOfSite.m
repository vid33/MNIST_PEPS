function [ env_out ] = PEPS_EnvironmentOfSite(A, B, location, MPO)

    %funciton either takes left and right env. MPS, and calculates the
    %enviroment on either the left or the right-hand boundaris, or it takes
    %and left and right environment MPS and an MPO sandwitche by these, 
    %and calculates the envirnment of a
    %site away from the left- or right- boundaries.

    [N, ~] = size(A);   
    
    % for left boundary specify site as [k,0], for right
    % as [0, k], and when away from left/right boundaries,
    % as [k, k]
    if location(1) ~= location(2) && min(abs(location)) ~= 0
        error('Location vector not specified correctly.');
    elseif location(1) == location(2) && nargin < 4
        error('For env. away from left- and right- boundaries an MPO needs to be specified.');
    end
    
    if location(1) == 0 ||  location(2) == 0
    
        position_vert = max(location);

        side = 'left';
        if location(1) == 0
            side = 'right';
        elseif location(2) == 0
            side = 'left';
        end

        if position_vert > 1
            rho_down = Contract({ A{1}, B{1} }, {[-1, 1], [-2, 1]});
        end
        if position_vert > 2
            for kk=2:position_vert-1
                rho_down = Contract({rho_down, A{kk},B{kk}}, {[1, 2], [-1, 3, 1], [-2, 2, 3]});            
            end
        end

        if position_vert < N
            rho_up = Contract({A{N}, B{N} }, {[1, -1], [-2, 1]});
        end

        if position_vert < N-1
            for kk=N-1:-1:position_vert+1
                rho_up = Contract({rho_up, A{kk},B{kk}}, {[1, 2], [1, 3, -1], [2, -2, 3]});            
            end
        end 

        if position_vert == N && strcmp(side, 'left') == true
            env_out = Contract({rho_down, B{N}}, {[-2, 1], [1, -1]});
        elseif position_vert == N && strcmp(side, 'right') == true
            env_out = Contract({rho_down, A{N}}, {[1, -1], [-2, 1]});
        elseif position_vert == 1 && strcmp(side, 'left') == true
            env_out = Contract({rho_up, B{1}}, {[-1, 1], [1, -2]});
        elseif position_vert == 1 && strcmp(side, 'right') == true
            env_out = Contract({rho_up, A{1}}, {[1, -1], [1, -2]});
        elseif strcmp(side, 'left') == true
            env_out = Contract({rho_up, B{position_vert}, rho_down}, {[-1, 1], [1, 2, -2], [-3, 2]});
        elseif strcmp(side, 'right') == true
            env_out = Contract({rho_up, A{position_vert}, rho_down}, {[1, -1], [1, -3, 2], [2, -2]});   
        end
    
    elseif location(1) == location(2)
        position_vert = location(1);
        if position_vert > 1
            rho_down = Contract({ A{1}, MPO{1}, B{1} }, ...
                    {[-1, 1],[-2, 2, 1], [-3, 2]});
        end
        if position_vert > 2
            for kk=2:position_vert-1
                rho_down = Contract({rho_down, A{kk}, MPO{kk}, B{kk}},...
                    {[1, 2, 3], [-1, 4, 1], [-2, 5, 2, 4], [-3, 3, 5]});  
                %%%
                % NB v1= rand(2, 3) v2 = rand(2,3); Contract({v1, v2}, {[1,
                % 2], [2,1]}) compiles, calculating trace(v1*reshape(v2,
                % 3,2)). Potentially a nasty source of bugs.
            end
        end
        if position_vert < N
            rho_up = Contract({A{N}, MPO{N}, B{N} }, {[1, -1], [2, -2, 1], [-3, 2]});
        end

        if position_vert < N-1
            for kk=N-1:-1:position_vert+1
                rho_up = Contract({rho_up, A{kk}, MPO{kk}, B{kk}}, ...
                    {[1, 2, 3], [1, 4, -1], [2, 5, -2, 4] [3, -3, 5]});            
            end
        end 
        if position_vert == N
            env_out = Contract({A{position_vert}, B{position_vert}, rho_down},...
                    {[-3, 1], [2, -1], [1, -2, 2]}); 
        elseif position_vert == 1
            env_out = Contract({rho_up, A{position_vert}, B{position_vert}},...
                    {[1, -1, 2], [1, -3], [2, -2]});            
        else
            env_out = Contract({rho_up, A{position_vert}, B{position_vert}, rho_down},...
                    {[1, -1, 2], [1, -4, 3], [2, 4, -2], [3, -3, 4]}); 
        end
    end
    
end

