%% Demo getEllipse()
%This is a demonstration of the getEllipse function of the analyzePattern
%   toolkit. To get accurate measurements, you must begin with good data.
%   For ellipse analysis, good data requires that the SEM image of the
%   grating is aligned with the horizontal or vertical edge of the image.
%   You can read how the algorithm works to understand why.
%
%execute "help analyzePatterns" for more information
clear all;
close all;
warning('off');
clc;

%% Initialization

%Image Parameters
fn = 'Data/Mask_Circles.png';               %filename
r = 10;                                     %expected radius [um]
pitch = 20;                                 %expected pitch [um]
scale = 8*40/1107;                          %scale [um/pixels]

%Load Image
g = analyzePatterns();
g = g.setScale(scale);
g = g.setRadius(r);
g = g.setPitch(pitch);
g = g.setShowPlot(true);
%g = g.cropImage();                         %Use to crop by hand

crop = [267 111 626 455];                   %Acquired using line 20

g = g.loadImage(fn);                        %Loads the image
g = g.cropImage(crop);                      %Crop the image
g = g.atfFilter(6);                         %Apply adaptive threshold filter
g = g.findCenter(r/3,.5,6,.5);              %Find ellipse centers
g = g.inverseFilter();                      %Invert the image
g = g.isolateDots(3,.7);                    %Isolate each ellipse
g.getEllipse(2);                            %Measure each ellipse
                                            %Results are in units of [um]
