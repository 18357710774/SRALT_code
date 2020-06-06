clear;
clc;
close all ;

% addpath
addpath SRALT_toolbox ;    

%% define images' path
currentPath = cd ;    
imagePath = fullfile(currentPath,'data') ;
pointPath = fullfile(currentPath,'data') ; % path to files containing initial feature coordinates
SaltNoiseRatio = 0.4;
GauNoiseSig = 25;
userName = ['Windows' '_SaltNoise' num2str(SaltNoiseRatio) '_GauNoise' num2str(GauNoiseSig)];
para.optmet = 'ADMM'; % the used method for optimization

%% define the path for saving results
destRoot = fullfile(currentPath,['results_real\L1+'  para.optmet]) ;
destDir = fullfile(destRoot,userName) ;
if ~exist(destDir,'dir')
    mkdir(destRoot,userName) ;
end

%% define parameters
lambdac = 1; 
alpha1 = 0.1;
layout.xI = 4 ;
layout.yI = 4 ;
layout.gap = 4 ;
layout.gap2 = 2 ;
para.DISPLAY = 0 ;
para.saveStart = 1 ;
para.saveEnd = 1 ;
para.saveIntermedia = 0 ;
para.tau = [0.15 0.95] ; 

% for windows images
para.canonicalImageSize = [ 200 200  ];
para.canonicalCoords = [ 26  176  100; ...
                             24  24   184 ];
% parametric tranformation model
para.transformType = 'HOMOGRAPHY'; 
% one of 'TRANSLATION', 'EUCLIDEAN', 'SIMILARITY', 'AFFINE','HOMOGRAPHY'

para.numScales = 1 ; % if numScales > 1, we use multiscales;

% main loop
para.stoppingDelta = .01; % stopping condition of main loop
para.maxIter = 50; % maximum iteration number of main loops

% inner loop
para.inner_tol = 1e-6; 
para.muBound = 1e+30;
para.inner_maxIter = 500 ;        
para.continuationFlag = 1 ;  
para.rho0 = 1.25;

%% Get training images
% get initial transformation
transformationInit = 'AFFINE';

[fileNames, transformations, numImages] = get_training_images( imagePath, pointPath, userName, para.canonicalCoords, transformationInit) ;

para.p = 1;
para.lambdac = lambdac ;    
para.alpha = [alpha1 alpha1 1-2*alpha1];  
 %% Tensor main loop: do robust batch image alignment
[Do, A, E, xi] = SRALT_main(fileNames, transformations, numImages, para, destDir);
sralt_plot(destDir, numImages, para.canonicalImageSize, layout);