clear;

d=2;
N_large = 28;
downsample_factor = 4;

%Load training data
test_images = loadMNISTImages('data/t10k-images-idx3-ubyte');
test_labels = loadMNISTLabels('data/t10k-labels-idx1-ubyte');

SAMPLE_NO = size(test_images,2);

%Store images in cell, original and downsampled
test_images_cell = cell(1, SAMPLE_NO);
border_width=2;
for kk=1:SAMPLE_NO
    kk
    test_images_cell{kk} = reshape(test_images(:,kk), N_large, N_large);
    test_images_cell{kk} = imageAddBlackBorder( test_images_cell{kk}, border_width );
    test_images_cell{kk} = imageDownsample(test_images_cell{kk}, downsample_factor);
end
N_large = N_large+2*border_width; %black border
N_small = N_large/downsample_factor;

%work with small for now
clearvars  test_images
test_images = test_images_cell;
clearvars test_images_cell;

%draw a random image from training set
tmp = randi([1, SAMPLE_NO ]);
fprintf('Image %i, digit %i\n', tmp, test_labels(tmp));
%figure; colormap(gray(256)); image(training_images_cell{tmp}*256);
figure; colormap(gray(256)); image(test_images{tmp}*256)


% Display_network from the autoencoder code
%display_network(images(:,1:100)); % Show the first 100 images
%disp(labels(1:10));


%NxN grid of 2-dim vectors, as a cell, encoding the image data
[Phi] = generatePhiCell(N_small, test_images);
N = N_small;

fOut = sprintf('data/test_data_N=%i_d=%i_samples=%i', N, d, SAMPLE_NO);

save(fOut, 'test_images', 'test_labels', 'Phi', 'N', 'd', 'SAMPLE_NO');

clear;