%% Demo getEllipse()
%This is a demonstration of the getEllipse function of the analyzePattern
%   toolkit. To get accurate measurements, you must begin with good data.
%   For ellipse analysis, good data requires that the SEM image of the
%   ellipse/circles looks like rings and are arranged in an array.
%   You can read how the algorithm works to understand why.
%
%execute "help analyzePatterns" for more information
clear all;
close all;
warning('off');
clc;

%% Initialization

%Image Parameters
fn = 'Data/Square_100nm_100kX.tif';     %filename in current directory
r = 50;                                 %expected radius [nm]
pitch = 200;                            %expected pitch [nm]
scale = 500/168;                        %[nm/pixel]

c = analyzePatterns();                  %initialize analyzePattern object
c = c.setScale(scale);                  %sets the scale
c = c.setRadius(r);                     %sets anticipated radius
c = c.setPitch(pitch);                  %sets anticipated pitch
c = c.setShowPlot(true);                %Enable plots to show for each function

c = c.loadImage(fn);                    %Load the image file
c = c.cropImage();                      %Crop the image
c = c.atfFilter(6);                     %Apply an adaptive threshold filter
c = c.findCenter();                     %Find the center of each ellipse
c = c.isolateDots();                    %Isolate each ellipse
c.getEllipse(2);                        %Fit ellipse to the image data
                                        %The results are in units of nm