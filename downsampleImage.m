function [ image_out ] = downsampleImage( image_in, factor )

    [Ninit, ~] = size(image_in);
    
    Nnew = Ninit/factor;
  
    if rem(Nnew, 1) ~= 0
        fprintf('Original dim not divisible by factor.\n');
        return;
    else
        image_out = zeros(Nnew, Nnew);
    end
    
    for kk=1:Nnew
        for mm=1:Nnew
            image_out(kk, mm) = (1/4)*( image_in(kk*2, mm*2)+ image_in(kk*2-1, mm*2) ...
                                        + image_in(kk*2, mm*2-1) + image_in(kk*2-1, mm*2-1));         
        end
    end
    image_out = image_out/max(reshape(image_out, 1, Nnew^2)); %larges value 1, pure white
    

end

