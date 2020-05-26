%% Demo getGrating
%This is a demonstration of the getGrating function of the analyzePattern
%   toolkit. To get accurate measurements, you must begin with good data.
%   For grating analysis, good data requires that the SEM image of the
%   grating is aligned with the horizontal or vertical edge of the image.
%   You can read how the algorithm works to understand why.
%
%execute "help analyzePatterns" for more information
clear all;
clc;
scale = 200/134;                        %scale [nm/pixels]
c = {'b' 'g' 'r' 'm' 'c' 'k' 'b--' 'g--' 'r--' 'm--' 'c--' 'k--'};

%% Determine Crop Coordinates
%   Run this cell to determine the crop coordinates

%index of file to load for the crop program
fileNumber = 1;

%This is a list of all the files that will be analyzed
fn = {'Data/S1_200kX.tif' ...
    'Data/S2_200kX.tif' ...
    'Data/S3_200kX.tif'};

%Loads an image file and waits for the user to click on two crop position
%   Comment out this section if it is not needed
% g = analyzePatterns();
% g = g.setScale(scale);
% g = g.loadImage(fn{fileNumber});
% g = g.setShowPlot(true);
% g.cropImage();

%% Record the crop positions
cropPos = [[66 121 948 789]; ...
    [73 129 941 806]; ...
    [73 157 934 826]];

%% Automatically analyzes the grating pattern on the images
g = analyzePatterns();
g = g.setScale(scale);

p = zeros([1 length(fn)]);
w1 = p;
w2 = p;
I = [];
for i = 1:length(fn)
    c = g.loadImage(fn{i});
    c = c.cropImage(cropPos(i,:));
    c = c.rotateImage();
    c = c.gaussianFilter();
    c = c.setShowPlot(true);
    [p(i), w1(i), w2(i), tmp] = c.getGrating(2, 0.55);
end