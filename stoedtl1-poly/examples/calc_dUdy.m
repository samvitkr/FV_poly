clear;  clc;  close all;
addpath 'C:\Users\samvi\Dropbox\SimonsProject\Finite_Vol_data\dmsuite'
% automatically add code modules to matlab path
functionCallStack = dbstack;
[scriptFolder, ~, ~] = fileparts(which(functionCallStack(1).file));
[sourceFolder, ~, ~] = fileparts(scriptFolder);
addpath(fullfile(sourceFolder, 'data_readers'), fullfile(sourceFolder, 'utilities'));

% user input: specify run parameters. This assumes your data is stored in a directory with full path dataFolder/runName
%runName = 'canonical_4pi_2pi_5';
%runFolder = fullfile(dataFolder, runName);
runFolder='C:\Users\samvi\Dropbox\SimonsProject\Finite_Vol_data\Wi6p67';
m=matfile('C:\Users\samvi\Dropbox\SimonsProject\Finite_Vol_data\Wi6p67\mean_profiles.mat',...
    'Writable',true)
nCheb=220;
[yCheb, ~] = chebdif(nCheb, 1);
dy = FirstDerivativeYChebyshev(nCheb);
m.dUdy=dy.compute_derivative(m.Um);
