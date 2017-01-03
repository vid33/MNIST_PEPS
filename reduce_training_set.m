clear;

%generate toy data files with fewer training images, for testing

d=2;
N = 8;
SAMPLE_NO = 60000;

fIn = sprintf('data/training_data_N=%i_d=%i_samples=%i', N, d, SAMPLE_NO);
load(fIn); clearvars fIn;

SAMPLE_NO = 100;

training_images(SAMPLE_NO+1:end) = [];
training_labels(SAMPLE_NO+1:end) =[];
Phi(SAMPLE_NO+1:end) = [];

fOut = sprintf('data/training_data_N=%i_d=%i_samples=%i', N, d, SAMPLE_NO);

save(fOut);

clear;