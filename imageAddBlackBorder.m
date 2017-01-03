function [ image_out ] = imageAddBlackBorder( image_in, border_width )

    [N, ~] = size(image_in);
    
    image_out = zeros(N+2*border_width, N+2*border_width);
    
    image_out( (1+border_width):(end-border_width), (1+border_width):(end-border_width)) = image_in;

end

