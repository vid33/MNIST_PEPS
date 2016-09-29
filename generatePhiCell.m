function [ Phi_cell ] = generatePhiCell( N, images_cell )


    Phi_cell = cell(1, size(images_cell, 2));
    
    
    for nn=1:size(images_cell,2)
        nn
        Phi_cell{nn} = zeros(N, N, 2);

                Phi_cell{nn}(:, :, 1) = cos(pi/2*images_cell{nn});
                Phi_cell{nn}(:, :, 2) = sin(pi/2*images_cell{nn});
    end
    
end

