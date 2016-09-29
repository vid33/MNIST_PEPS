clear;

d=2;
N_large = 28;
N_small = 14;

%Load training data
training_images = loadMNISTImages('data/train-images-idx3-ubyte');
training_labels = loadMNISTLabels('data/train-labels-idx1-ubyte');

SAMPLE_NO = size(training_images,2);

%Store images in cell, original and downsampled
training_images_cell = cell(1, SAMPLE_NO);
for kk=1:SAMPLE_NO
    kk
    training_images_cell{kk} = reshape(training_images(:,kk), N_large, N_large);
    training_images_cell{kk} = downsampleImage(training_images_cell{kk}, 2);
end


%work with small for now
clearvars  training_images
training_images = training_images_cell;
clearvars training_images_cell;

%draw a random image from training set
tmp = randi([1, SAMPLE_NO ]);
fprintf('Image %i, digit %i\n', tmp, training_labels(tmp));
%figure; colormap(gray(256)); image(training_images_cell{tmp}*256);
figure; colormap(gray(256)); image(training_images{tmp}*256)


% Display_network from the autoencoder code
%display_network(images(:,1:100)); % Show the first 100 images
%disp(labels(1:10));


%NxN grid of 2-dim vectors, as a cell, encoding the image data
[Phi] = generatePhiCell(N_small, training_images);
N = N_small;

fOut = sprintf('data/training_data_N=%i_d=%i_samples=%i', N, d, SAMPLE_NO);

save(fOut, 'training_images', 'training_labels', 'Phi', 'N', 'd', 'SAMPLE_NO');

clear;