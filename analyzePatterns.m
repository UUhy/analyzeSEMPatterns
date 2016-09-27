%analyzePatterns Class
%
%This class is designed to help measure simple patterns from a grayscale
%image.  It can currently measure:
%   1)  Ellipse/circles that are sparsely spaced
%   2)  Gratings
%
%Long Chang, UH, March 17 2016
%
%See also getEllipse, getGrating, getGratingLER

classdef analyzePatterns
    
    properties(GetAccess=public, SetAccess=protected, Hidden=false)
        %The get and set privileges are as shown above.
        
        %image scale in [nm/pixel]
        scale;
        
        %raw Image
        Iraw;
        
        %processed Image
        I;
        
        %dimension of processed image
        dim;
        
        %radius of circle in units of pixel
        r;
        
        %pitch of pattern in units of pixel
        pitch;
        
        %half pitch of the pattern in units of pixel
        halfPitch;
        
        %show plot enable
        showPlot;
        
        %center points of each pattern
        cx;
        cy;
        
        %an array of individual patterns
        dot;
        
        %an array of filtered individual patterns
        gdot;
        
        %a matrix of fitted ellipse parameters
        pEllipse;
        
        %image processing toolkit
        IPK;
    end
    
    properties(GetAccess=public, SetAccess=protected, Hidden=true)
        %No Private Variables
    end
    
    properties(GetAccess=private, Constant=true, Hidden=true)
        %No Constants
    end
    
    methods (Access=public, Sealed=true)
        
        function obj = analyzePatterns()
        %obj = analyzePatterns()
        %
        %   Output:
        %       obj     	:   class instance
        %
        %   Description:
        %   Constructor for the analyzePatterns class.
     
            obj.scale = 1;
            obj.Iraw = [];
            obj.I = [];
            obj.r = [];
            obj.pitch = [];
            obj.halfPitch = [];
            obj.showPlot = false;
            obj.cx = [];
            obj.cy = [];
            obj.dot = [];
            obj.gdot = [];
            obj.pEllipse = [];
            obj.IPK = imageProcessingKit();
        end
        
        function obj = setScale(obj, scale)
        %obj = setScale(obj)
        %
        %   Output:
        %       obj         :   class instance
        %
        %   Input:
        %       obj         :   class instance
        %       scale       :   double  =   scale in [units/pixel]
        %
        %   Description:
        %   The default scale is pixels.  Setting the scale parameter will
        %   cause all outputs to be reported in [units].
            
            obj.scale = scale;
        end
        
        function obj = setRadius(obj, radius)
        %obj = setRadius(obj, radius)
        %
        %   Output:
        %       obj         :   class instance
        %
        %   Input:
        %       obj         :   class instance
        %       radius      :   integer  =   expected radius of circle
        %
        %   Description:
        %   Set the radius parameter.  The radius parameter is the
        %   expected radius of circles.
            
            if ~isempty(obj.scale)
                obj.r = round(radius/obj.scale);
            else
                error('Use setScale to set a scale');
            end
        end
        
        function obj = setPitch(obj, pitch)
        %obj = setPitch(obj, pitch)
        %
        %   Output:
        %       obj         :   class instance
        %
        %   Input:
        %       obj         :   class instance
        %       pitch       :   integer  =   expected pitch of the array
        %
        %   Description:
        %   Set the pitch parameter.  The pitch parameter is the expected
        %   pitch of an array.
            
            if ~isempty(obj.scale)
                obj.pitch = round(pitch/obj.scale);
                obj.halfPitch = round(obj.pitch/2);
            else
                error('Use setScale to set a scale');
            end
        end
        
        function obj = setShowPlot(obj, show)
        %obj = setShowPlot(obj, show)
        %
        %   Output:
        %       obj         :   class instance
        %
        %   Input:
        %       obj         :   class instance
        %       show        :   boolean     =   show plot when available
        %
        %   Description:
        %   Set the showPlot parameter.  If showPlot is set to true, then
        %   some functions will show plots to help understand the result of
        %   the function.
            
            obj.showPlot = show;
            
        end
        
        function obj = loadImage(obj, filename)
        %obj = loadImage(obj, filename)
        %
        %   Output:
        %       obj         :   class instance
        %
        %   Input:
        %       obj         :   class instance
        %       filename    :   string  =   name of image file
        %
        %   Description:
        %   Loads and image from a file.  Only grayscale images are
        %   acceptable.
            
            obj.Iraw = imread(filename);
            obj.I = obj.Iraw;
            if size(obj.Iraw,3) > 1
                error('loadImage():  The image needs to be a grayscale image');
            end
        end
        
        function obj = setImage(obj, img)
        %obj = setImage(obj, img)
        %
        %   Output:
        %       obj         :   class instance
        %
        %   Input:
        %       obj         :   class instance
        %       img         :   NxM array   =   image matrix
        %
        %   Description:
        %   Set the image parameters.  An rgb image can be processed into
        %   grayscale before using this function to set the image for
        %   further processing.
            
            obj.Iraw = img;
            obj.I = obj.Iraw;

        end
        
        function obj = cropImage(obj, varargin)
        %obj = cropImage(obj, corners = [])
        %
        %   Output:
        %       obj         :   class instance
        %
        %   Input:
        %       obj         :   class instance
        %       corners*    :   [x1 y1 x2 y2]   =   corners of a box to crop
        %                          
        %   Description:
        %   Crops the image to remove unwanted data.  If the corners
        %   parameter is not specified, then the user will need to specify
        %   the crop boundary by clicking on the figure.
        
            if nargin > 1
                tmp = varargin{1};
                x = [tmp(1) tmp(3)];
                y = [tmp(2) tmp(4)];
            else
                figure
                colormap(gray)
                imagesc(obj.I)
                title('Click top left corner and bottom right corner to crop image')
                [x, y] = ginput(2);
            end
            x = round(x);
            y = round(y);
            obj.I = double(obj.I(y(1):y(2),x(1):x(2)));   %crop image
            obj.dim = size(obj.I);
            
            if obj.showPlot
                figure
                colormap(gray)
                imagesc(obj.I)
                title(['xy1(' num2str(x(1)) ',' num2str(y(1)) ')' '   xy2(' num2str(x(2)) ',' num2str(y(2)) ')'])
            end
        end
        
        function obj = rotateImage(obj, varargin)
        %obj = rotateImage(obj, ang = 90)
        %
        %   Output:
        %       obj         :   class instance
        %
        %   Input:
        %       obj         :   class instance
        %       ang*        :   double  =   degrees to rotate
        %                                   counterclockwise
        %
        %   Description:
        %   Rotates the image counterclockwise about its center at the
        %   specified angle in degrees.
        %
        %   showPlot:
        %   Shows image after rotation
        
            ang = 90;
            if nargin > 1
                ang = varargin{1};
            end
            obj.I = obj.IPK.rotateImage(obj.I,ang);
            
            if obj.showPlot
                figure
                colormap(gray)
                imagesc(obj.I)
                title('Image processed by rotateImage')
            end
        end
        
        function obj = alignX(obj)
        %obj = alignX(obj)
        %
        %   Output:
        %       obj         :   class instance
        %
        %   Input:
        %       obj         :   class instance
        %
        %   Description:
        %   The function is used to align a line on the image along the x
        %   axis.  The line is determined by 2 clicks that the user
        %   performs on the image.
        %
        %   showPlot:
        %   Shows image after aligned

            figure
            colormap(gray)
            imagesc(obj.I)
            title('Click on 2 points from left to right to align with X axis')
            [x, y] = ginput(2);
            ang = atan(diff(y)/diff(x))*180/pi;
            obj.I = obj.IPK.rotateImage(obj.I,ang);
            
            if obj.showPlot
                figure
                colormap(gray)
                imagesc(obj.I)
                title(['Angle(' num2str(ang) ')']);
            end
        end
        
        function obj = alignY(obj)
        %obj = alignY(obj)
        %
        %   Output:
        %       obj         :   class instance
        %
        %   Input:
        %       obj         :   class instance
        %
        %   Description:
        %   The function is used to align a line on the image along the y
        %   axis.  The line is determined by 2 clicks that the user
        %   performs on the image.
        %
        %   showPlot:
        %   Shows image after aligned

            figure
            colormap(gray)
            imagesc(obj.I)
            title('Click on 2 points from top to bottom to align with Y axis')
            [x, y] = ginput(2);
            ang = atan(diff(y)/diff(x))*180/pi;
            obj.I = obj.IPK.rotateImage(obj.I,ang+90);
            
            if obj.showPlot
                figure
                colormap(gray)
                imagesc(obj.I)
                title(['Angle(' num2str(ang) ')']);
            end
        end
        
        function obj = gaussianFilter(obj, varargin)
        %obj = gaussianFilter(obj, blur = 3)
        %
        %   Output:
        %       obj         :   class instance
        %
        %   Input:
        %       obj         :   class instance
        %       blur*       :   integer     =   amount of blur
        %
        %   Description:
        %   Filters the image by convolving it with a 2D gaussian function.
        %   The blur parameter is used as sigma when generating the
        %   gaussian function.
        %
        %   showPlot:
        %   Shows the filtered image
        
            blur = 3;
            if nargin > 1
                blur = varargin{1};
            end
            obj.I = obj.IPK.gaussianFilter(obj.I, blur);
            
            if obj.showPlot
                figure
                colormap(gray)
                imagesc(obj.I)
                title('Image processed by gaussianFilter');
            end
        end
        
        function obj = atfFilter(obj, varargin)
        %obj = atfFilter(obj, radius = 10)
        %
        %   Output:
        %       obj         :   class instance
        %
        %   Input:
        %       obj         :   class instance
        %       radius*     :   integer     =   radius of disk filter
        %
        %   Description:
        %   The adaptive threshold filter is very useful for removing
        %   background noise.  A background image is acquired by processing
        %   the image with a disk filter using a large blur.  Then the
        %   background image is subtracted from the raw image.
        %
        %   showPlot:
        %   Shows the filtered image
        
            radius = 10;
            if nargin > 1
                radius = varargin{1};
            end
            obj.I = obj.IPK.atf(obj.I, radius);
            
            if obj.showPlot
                figure
                colormap(gray)
                imagesc(obj.I)
                title('Image processed by atfFilter');
            end
        end
        
        function obj = findCenter(obj, varargin)
        %obj = findCenter(obj, blur = radius, threshold = 0.5, noiseBlur = 6, noiseThreshold = 0.5)
        %
        %   Output:
        %       obj         :   class instance
        %
        %   Input:
        %       obj             :   class instance
        %       blur*           :   integer         =   blur amount
        %       threshold*      :   double [0 to 1] =   threshold after blur
        %       noiseBlur*      :   integer         =   blur amount to reduce
        %                                               noise
        %       noiseThreshold* :   double [0 to 1] =   threshold after a
        %                                               noise blur
        %
        %   Description:
        %   The centers for a set of circles can be found by:
        %       1)  noiseBlur to reduce noise
        %       2)  noiseThreshold to remove background
        %       3)  blur the circles into a peak
        %       4)  threshold to remove background
        %       5)  find the peaks
        %   This algorithm assumes:
        %       1)  circles are approximately the same size
        %       2)  circles are at least 2*diameter apart
        %
        %   showPlot:
        %   Shows the center points overlayed on the image
        
            blur  = obj.r;
            threshold = 0.5;
            noiseBlur = 6;
            noiseThreshold = 0.5;
            
            if nargin > 1
                blur = varargin{1};
            end
            if nargin > 2
                threshold = varargin{2};
            end
            if nargin > 3
                noiseBlur = varargin{3};
            end
            if nargin > 4
                noiseThreshold = varargin{4};
            end
            
            gimg = obj.IPK.gaussianFilter(obj.I,noiseBlur);
            gimg = obj.IPK.threshold(gimg, noiseThreshold);
            gimg = obj.IPK.gaussianFilter(gimg,blur);
            gimg = obj.IPK.threshold(gimg, threshold);
            fimg = obj.IPK.peakImg(gimg,1);

            [obj.cy, obj.cx] = find(fimg==1);

            if obj.showPlot
                figure
                colormap(gray)
                imagesc(obj.I)
                hold on
                plot(obj.cx,obj.cy,'m*')
                hold off
            end
        end
        
        function obj = isolateDots(obj, varargin)
        %obj = isolateDots(obj, blur = 3, threshold = 0.5)
        %
        %   Output:
        %       obj         :   class instance
        %
        %   Input:
        %       obj         :   class instance
        %       blur*           :   integer         =   blur amount
        %       threshold*      :   double [0 to 1] =   threshold after blur
        %
        %   Description:
        %   Each circle will need to be isolated from the image before
        %   the fitEllipse function can be used.  This function can isolate
        %   each circle if the following is known:
        %       1)  centers of each circle (use findCenter)
        %       2)  approximate distance between circles (use setPitch)
        %
        %   showPlot:
        %   Shows the datapoints to be fitted overlayed on the image
        
            blur = 3;
            threshold = 0.5;
            
            if nargin > 1
                blur = varargin{1};
            end
            if nargin > 2
                threshold = varargin{2};
            end
        
            %Preallocate memory
            obj.dot = cell([1 length(obj.cx)]);
            obj.gdot = obj.dot;

            %Crop dot from image
            for i = 1:length(obj.dot)
                if(obj.cx(i)-obj.halfPitch < 1)
                    xl = 1;
                else xl = obj.cx(i)-obj.halfPitch;
                end
                if(obj.cx(i)+obj.halfPitch > obj.dim(2))
                    xu = obj.dim(2);
                else xu = obj.cx(i)+obj.halfPitch;
                end
                if(obj.cy(i)-obj.halfPitch < 1)
                    yl = 1;
                else yl = obj.cy(i)-obj.halfPitch;
                end
                if(obj.cy(i)+obj.halfPitch > obj.dim(1))
                    yu = obj.dim(1);
                else yu = obj.cy(i)+obj.halfPitch;
                end    
                obj.dot{i} = obj.I(yl:yu,xl:xu);
                tmp = obj.IPK.gaussianFilter(obj.dot{i},blur);
                tmp = obj.IPK.threshold(tmp,threshold);
                obj.gdot{i} = tmp > 0;    
            end
            
            if obj.showPlot
                figure
                colormap(gray)
                imagesc(obj.I)
                %Fit ellipse to dot
                hold on
                for i = 1:length(obj.dot)
                    [tY, tX] = ind2sub(size(obj.gdot{i}),find(obj.gdot{i}>0));
                    plot(obj.cx(i)+tX-obj.halfPitch, obj.cy(i)+tY-obj.halfPitch,'m*');
                end
                hold off
            end
        end
        
        function obj = getEllipse(obj, varargin)
        %obj = getEllipse(obj, type = 2)
        %
        %   Output:
        %       obj         :   class instance
        %
        %   Input:
        %       obj         :   class instance
        %       type        :   integer [1, 2, 3]   =   1   = fmin
        %                                           =   2   = fmean
        %                                           =   3   = fmax
        %
        %   Description:
        %   This function will fit an ellipse to all the circles that are
        %   found in the image.  Make sure you set showPlot to true before
        %   running this function in order to see the results.
        %
        %   This function assumes that the processed isolated circles looks
        %   like a ring.  An ellipse will be fit to the inner diameter of
        %   the ring, the mean diameter of the ring and the outer diameter
        %   of the ring.  The showPlot from the isolateDots function plots
        %   the data that is to be fitted.
        %
        %   Example:
        %   fn = 'Square_100nm_100kX.tif';      %filename in current directory
        %   r = 50;                             %expected radius [nm]
        %   pitch = 200;                        %expected pitch [nm]
        % 
        %   c = analyzePatterns();
        %   c = c.setScale(scale);
        %   c = c.setRadius(r);
        %   c = c.setPitch(pitch);
        %   c = c.setShowPlot(true);
        % 
        %   c = c.loadImage(fn);
        %   c = c.cropImage();
        %   c = c.atfFilter(6);
        %   c = c.findCenter();
        %   c = c.isolateDots();
        %   c.getEllipse(1);
        %
        %   See also setShowPlot, cropImage, atfFilter, findCenter,
        %   isolateDots, rotateImage, alignX, setImage, removeData
        
            type = 2;
            a = obj.r;
            b = obj.r;
            nseg = 128;
            
            if nargin > 1
                type = varargin{1};
            end
            
            %Preallocate memory
            obj.pEllipse = zeros([length(obj.dot) 8 3]);

            %Fit ellipse to dot
            for i = 1:length(obj.dot)
                [tY, tX] = ind2sub(size(obj.gdot{i}),find(obj.gdot{i}>0));
                [fmean, fmin, fmax] = obj.fitEllipse(tX,tY,a,b,nseg);

                %Store only the fit to the perimeter of dot
                obj.pEllipse(i,:,1) = fmin + [obj.cx(i)-obj.halfPitch obj.cy(i)-obj.halfPitch 0 0 0 0 0 0]; 
                obj.pEllipse(i,:,2) = fmean + [obj.cx(i)-obj.halfPitch obj.cy(i)-obj.halfPitch 0 0 0 0 0 0]; 
                obj.pEllipse(i,:,3) = fmax + [obj.cx(i)-obj.halfPitch obj.cy(i)-obj.halfPitch 0 0 0 0 0 0]; 
            end
            
            if obj.showPlot
                figure
                colormap(gray)
                imagesc(obj.I)
                hold on
                for i = 1:size(obj.pEllipse,1)
                    obj.drawEllipse(obj.pEllipse(i,1,type), obj.pEllipse(i,2,type), obj.pEllipse(i,3,type), obj.pEllipse(i,4,type), 0, 32, 'r-');
                    plot(obj.pEllipse(i,1,type), obj.pEllipse(i,2,type),'*')    
                end
                axis tight
                hold off
                title(['Mean('  num2str(mean(obj.pEllipse(:,8,type))*obj.scale) ')  Spread(' num2str(std(obj.pEllipse(:,8,type))*obj.scale) ')']);
                
                disp(['Mean Diameter = ' num2str(mean(obj.pEllipse(:,8,type))*obj.scale) ' [nm]'])
                disp(['Std Diameter = ' num2str(std(obj.pEllipse(:,8,type))*obj.scale) ' [nm]'])
            end
        end
        
        function [fmean, fmin, fmax] = fitEllipse(obj, x,y,a,b,varargin)
        %[fmean, fmin, fmax] = fitEllipse(obj, x, y, a, b, nseg = 128, debug = false)
        %
        %   Output:
        %       fmean   :   [x y a b phi ellipticity area diameter] =   fit to mean diameter of a ring
        %       fmin    :   [x y a b phi ellipticity area diameter] =   fit to inner diameter of a ring
        %       fmax    :   [x y a b phi ellipticity area diameter] =   fit to outer diameter of a ring       
        %
        %   Input:
        %       obj     :   class instance
        %       x       :   array   =   x data
        %       y       :   array   =   y data
        %       a       :   double  =   major axis length
        %       b       :   double  =   minor axis length
        %       nseg*   :   integer =   number of segments for averaging
        %       debug*  :   boolean =   debug mode enable 
        %
        %   Description:
        %   Automatically fits an ellipse to the data set using least
        %   square minimization.  Rotation in the ellipse is acceptable.
        %
        %   Future Development:
        %   Perhaps add a line edge rougness calculation to this script
        %   since we have all the information to calculate it
        
            nseg = 128;
            debug = false;

            if nargin>5
                nseg = varargin{1};
            end
            if nargin>6
                debug = varargin{2};
            end

            xc = mean(x);
            yc = mean(y);
            x  = x - xc;
            y  = y - yc;

            %Calculate tilt of the ellipse major axis
            [theta,rho] = cart2pol(x,y);

            ang = linspace(-pi,pi,nseg);
            rmean = zeros(size(ang));
            rmin = rmean;
            rmax = rmean;
            for i = 1:nseg
                ind = find(theta<ang(i));
                rmean(i) = mean(rho(ind));
                if isnan(rmean(i))
                    rmin(i) = NaN;
                    rmax(i) = NaN;
                else
                    rmin(i) = min(rho(ind));
                    rmax(i) = max(rho(ind));        
                end
                theta(ind) = [];
                rho(ind) = [];
            end
            idel = isnan(rmean);
            rmean(idel) = [];
            rmin(idel) = [];
            rmax(idel) = [];
            ang(idel) = [];

            %Anonymous Function
            options = optimset('Display','off','TolFun',1e-6);
            %Optimize Fit
            fc = lsqnonlin(@(g) (g(1)*cos(2*(ang+g(2)))+g(3))-rmean, [max(rmean)-min(rmean) 0 mean(rmean)], [], [],options);

            phi = fc(2);
            if phi > pi/2
                phi = pi-phi;
            elseif phi < -pi/2
                phi = pi+phi;
            end

            [xmean, ymean] = pol2cart(ang,rmean);
            [xmin, ymin] = pol2cart(ang,rmin);
            [xmax, ymax] = pol2cart(ang,rmax);

            %Rotate ellipse such that the major axis is aligned with the x axis
            xrmean = xmean*cos(-phi)+ymean*sin(-phi);
            yrmean = ymean*cos(-phi)-xmean*sin(-phi);
            xrmin = xmin*cos(-phi)+ymin*sin(-phi);
            yrmin = ymin*cos(-phi)-xmin*sin(-phi);
            xrmax = xmax*cos(-phi)+ymax*sin(-phi);
            yrmax = ymax*cos(-phi)-xmax*sin(-phi);
            
            %Fit Ellipse to Data
            %Initial Guess
            guess = [0 0 a b];   %[xc a yc b]

            %Anonymous Function
            % f = @(a) ((xr-a(1)).^2)/a(3).^2 + ((yr-a(2)).^2)/a(4).^2 -1;
            options = optimset('Display','off','TolFun',1e-6);

            %Optimize Fit
            fmean = lsqnonlin(@(g) ((xrmean-g(1)).^2)/g(3).^2 + ((yrmean-g(2)).^2)/g(4).^2 -1, guess, [], [],options);
            fmin = lsqnonlin(@(g) ((xrmin-g(1)).^2)/g(3).^2 + ((yrmin-g(2)).^2)/g(4).^2 -1, guess, [], [],options);
            fmax = lsqnonlin(@(g) ((xrmax-g(1)).^2)/g(3).^2 + ((yrmax-g(2)).^2)/g(4).^2 -1, guess, [], [],options);

            %Reconstruct the center of ellipse
            fmean(1) = fmean(1)*cos(phi)+fmean(2)*sin(phi);
            fmean(2) = fmean(2)*cos(phi)-fmean(1)*sin(phi);
            fmean(1) = fmean(1)+xc;
            fmean(2) = fmean(2)+yc;
            fmean(5) = phi;
            if fmean(3)>fmean(4)
                fmean(6) = fmean(3)/fmean(4);
            else
                fmean(6) = fmean(4)/fmean(3);
            end
            fmean(7) = fmean(3)*fmean(4)*pi;
            fmean(8) = 2*sqrt(fmean(7)/pi);

            fmin(1) = fmin(1)*cos(phi)+fmin(2)*sin(phi);
            fmin(2) = fmin(2)*cos(phi)-fmin(1)*sin(phi);
            fmin(1) = fmin(1)+xc;
            fmin(2) = fmin(2)+yc;
            fmin(5) = phi;
            if fmin(3)>fmin(4)
                fmin(6) = fmin(3)/fmin(4);
            else
                fmin(6) = fmin(4)/fmin(3);
            end
            fmin(7) = fmin(3)*fmin(4)*pi;
            fmin(8) = 2*sqrt(fmin(7)/pi);

            fmax(1) = fmax(1)*cos(phi)+fmax(2)*sin(phi);
            fmax(2) = fmax(2)*cos(phi)-fmax(1)*sin(phi);
            fmax(1) = fmax(1)+xc;
            fmax(2) = fmax(2)+yc;
            fmax(5) = phi;
            if fmax(3)>fmax(4)
                fmax(6) = fmax(3)/fmax(4);
            else
                fmax(6) = fmax(4)/fmax(3);
            end
            fmax(7) = fmax(3)*fmax(4)*pi;
            fmax(8) = 2*sqrt(fmax(7)/pi);

            if debug
                figure
                hold on
                plot(x+xc,y+yc,'b*');
                plot(xrmean+xc,yrmean+yc,'g*');
                plot(xrmin+xc,yrmin+yc,'g*');
                plot(xrmax+xc,yrmax+yc,'g*');
                obj.drawEllipse(fmean(1), fmean(2), fmean(3), fmean(4), 0, 32, 'r-');
                obj.drawEllipse(fmin(1), fmin(2), fmin(3), fmin(4), 0, 32, 'r-');
                obj.drawEllipse(fmax(1), fmax(2), fmax(3), fmax(4), 0, 32, 'r-');
                axis tight
                hold off  
            end
        end
        
        function [px, py] = drawEllipse(~, x, y, a, b, phi, nseg, S)
        %[px, py] = drawEllipse(obj, x, y, a, b, phi, nseg, S)
        %
        %   Output:
        %       px          :   array 	=   x points for the ellipse
        %       py          :   array 	=   y points for the ellipse
        %
        %   Input:
        %       obj         :   class instance
        %       x           :   double  =   x center of circle
        %       y           :   double  =   y center of circle
        %       a           :   double  =   radius of major axis
        %       b           :   double  =   radius of minor axis
        %       phi         :   double  =   angle of ellipse
        %       nseg        :   integer =   number of segments
        %       S           :   string  =   plot symbols such as 'r-'
        %
        %   Description:
        %   Plots an ellipse

            theta = 0 : (2 * pi / nseg) : (2 * pi);
            px = a*cos(theta)*cos(phi) - b*sin(theta)*sin(phi) + x;
            py = a*cos(theta)*sin(phi) + b*sin(theta)*cos(phi) + y;

            plot(px, py, S);
        end
        
        function [px, py] = drawCircle(~, x, y, r, nseg, S)
        %[px, py] = drawCircle(obj, x, y, r, nseg, S)
        %
        %   Output:
        %       px          :   array 	=   x points for the circle
        %       py          :   array 	=   y points for the circle
        %
        %   Input:
        %       obj         :   class instance
        %       x           :   double  =   x center of circle
        %       y           :   double  =   y center of circle
        %       r           :   double  =   radius of circle
        %       nseg        :   integer =   number of segments
        %       S           :   string  =   plot symbols such as 'r-'
        %
        %   Description:
        %   Plots a circle

            theta = 0 : (2 * pi / nseg) : (2 * pi);
            px = r * cos(theta) + x;
            py = r * sin(theta) + y;

            plot(px, py, S);
        end
        
        function [px, py] = drawRectangle(~, x, y, a, b, phi, n, nseg, S)
        %[px, py] = drawRectangle(obj, x, y, w, h, phi, n, nseg, S)
        %
        %   Output:
        %       px          :   array 	=   x points for the circle
        %       py          :   array 	=   y points for the circle
        %
        %   Input:
        %       obj         :   class instance
        %       x           :   double  =   x center of circle
        %       y           :   double  =   y center of circle
        %       a           :   double  =   length of major axis
        %       b           :   double  =   height of major axis
        %       phi         :   double  =   angle of major axis
        %       n           :   double  =   squareness factor
        %       nseg        :   integer =   number of segments
        %       S           :   string  =   plot symbols such as 'r-'
        %
        %   Description:
        %   Plots a circle

            theta = 0 : (2 * pi / nseg) : (2 * pi);
            px = a*cos(theta).^(2/n)*cos(phi) - b*sin(theta).^(2/n)*sin(phi) + x;
            py = a*cos(theta).^(2/n)*sin(phi) + b*sin(theta).^(2/n)*cos(phi) + y;

            plot(px, py, S);
        end
        
        function obj = removeData(obj, varargin)
        %obj = removeData(obj, nPointsMax = 20)
        %
        %   Output:
        %       obj         :   class instance
        %
        %   Input:
        %       obj         :   class instance
        %       nPointsMax  :   integer     =   maximum number of points to
        %                                       remove
        %
        %   Description:
        %   Remove unwanted getEllipse data by clicking on the centers of
        %   each fitted circle and then pressing Enter.  This function is
        %   used to remove erroneous fits which adversely affects the
        %   statistics.
        
            nPointsMax = 20;
        
            type = 3;
            if nargin > 1
                type = varargin{1};
            end
            
            figure
            colormap(gray)
            imagesc(obj.I)
            hold on
            for i = 1:size(obj.pEllipse,1)
                obj.drawEllipse(obj.pEllipse(i,1,type), obj.pEllipse(i,2,type), obj.pEllipse(i,3,type), obj.pEllipse(i,4,type), 0, 32, 'r-');
                plot(obj.pEllipse(i,1,type), obj.pEllipse(i,2,type),'*')  
            end
            axis tight
            hold off
            title('Click at centers of circles to remove.  Press Enter when done.');

            [x, y] = ginput(nPointsMax);
            pIndex = obj.IPK.findCluster([obj.pEllipse(:,1,type) ; x],[obj.pEllipse(:,2,type) ; y],30);

            index = 1;
            rindex = zeros([1 length(pIndex)- size(obj.pEllipse,1)]);
            for i = size(obj.pEllipse,1)+1:length(pIndex)
                rindex(index) = pIndex{i}(2);
                index = index+1;
            end
            obj.pEllipse(rindex,:,:) = [];

            if obj.showPlot
                figure
                colormap(gray)
                imagesc(obj.I)
                for i = 1:size(obj.pEllipse,1)
                    hold on
                    obj.drawEllipse(obj.pEllipse(i,1,type), obj.pEllipse(i,2,type), obj.pEllipse(i,3,type), obj.pEllipse(i,4,type), 0, 32, 'r-');
                    plot(obj.pEllipse(i,1,type), obj.pEllipse(i,2,type),'*')
                    axis tight
                    hold off    
                end    
            end

            disp(['Mean Diameter = ' num2str(mean(obj.pEllipse(:,8,type))*obj.scale) ' [nm]'])
            disp(['Std Diameter = ' num2str(std(obj.pEllipse(:,8,type))*obj.scale) ' [nm]'])    
        end
        
        function [pitchMean, width1, width2, varargout] = getGrating(obj, varargin)
        %[pitch width1 width2 data] = getGrating(obj, blur = 0, threshold = 0.8, type = 1)
        %
        %   Output:
        %       pitch       :   double                  =   mean pitch
        %       width1      :   double                  =   smaller width
        %       width2      :   double                  =   larger width
        %       data        :   {p w1 w2 edge}          =   data (type = 1)
        %                   :   {p w1 w2 w3 w4 edge}    =   data (type = 2)
        %
        %   Input:
        %       obj         :   class instance
        %       blur        :   integer         =   blur amount
        %       threshold   :   double [0 to 1] =   threshold value
        %       type        :   integer [1, 2]  =   peaks per edge
        %                   :   1   =   1 peak per edge
        %                   :   2   =   2 peaks per edge
        %
        %   Description:
        %   Measures the pitch and duty cycle of a grating.  The lines of
        %   the grating needs to be near parallel to the y axis.
        %
        %   The pitch measured with this function is consistent as long as
        %   the blue is at least 1 to suppress noise.  The duty cycle
        %   measured depends on the amount of blur.  See example for
        %   specifics.
        %
        %   Due to the interaction of the electron beam with the sample,
        %   the derivative on an edge oten results in a negative peak and a
        %   positive peak.  Sometimes the magnitude of the negative peak
        %   and the positive peaks are similar and it becomes difficult to
        %   identify a single peak to assign as the edge.  In such case,
        %   set the type parameter to 2 and the algorithm will assume 2
        %   peaks per edge and calculate 4 widths.  At this point, it is up
        %   to the user to decides how to assign these widths.  The outputs
        %   width3 and width4 is only available when the type is set to 2.
        %
        %   Example:
        %     fn = 'Grating_200nm_100kX.tif';
        %     scale = 500/167;
        %     c = {'b' 'g' 'r' 'm' 'c' 'k' 'b--' 'g--' 'r--' 'm--' 'c--' 'k--'};
        %     g = analyzePatterns();
        %     g = g.setScale(scale);
        %     g = g.loadImage(fn);
        %     g = g.cropImage();
        % 
        %     x = 0:9;
        %     index = 1;
        % 
        %     figure
        %     for j = -3:3
        %         pitch = zeros([1 10]);
        %         dutyCycle = pitch;
        % 
        %         for i = 1:9
        %             [p, w1] = g.getGrating(i);
        %             pitch(i+1) = p;
        %             dutyCycle(i+1) = w1/p;
        %         end
        %         h = g.rotateImage(j);
        %         [p, w1] = h.getGrating();
        %         pitch(1) = p;
        %         dutyCycle(1) = w1/p;
        % 
        %         subplot(1,2,1)
        %         hold on
        %         plot(x,pitch,c{index})
        %         xlabel('Blur [pixels]');
        %         ylabel('Pitch [nm]');
        %         subplot(1,2,2)
        %         hold on
        %         plot(x,dutyCycle,c{index})
        %         xlabel('Blur [pixels]');
        %         ylabel('Duty Cycle [%]');
        %         index = index + 1;
        %     end
        %     legend('-3 degrees','-2 degrees','-1 degrees',' 0 degrees', ...
        %         ' 1 degrees',' 2 degrees',' 3 degrees');
        %     subplot(1,2,1)
        %     hold on
        %     legend('-3 degrees','-2 degrees','-1 degrees',' 0 degrees', ...
        %         ' 1 degrees',' 2 degrees',' 3 degrees'); 
        %
        %   See also setShowPlot, cropImage, atfFilter, findCenter,
        %   isolateDots, rotateImage, alignX, setImage, removeData,
        %   getGratingLER

            blur = 0;
            threshold = 0.8;
            type = 1;
            
            if nargin > 1
                blur = varargin{1};
            end
            if nargin > 2
                threshold = varargin{2};
            end
            if nargin > 3
                type = varargin{3};
            end
            
            if blur ~= 0
                img = obj.IPK.gaussianFilter(obj.I, blur);
            else
                img = obj.I;
            end
            
            a = mean(img);
            da = diff(a);
            pa = abs(da);
            
            iPeak = obj.IPK.findPeaks(pa);
            iEdge = iPeak(pa(iPeak) > max(pa)*threshold);
            try
                if type == 1
                    p = [diff(iEdge(1:2:end)) diff(iEdge(2:2:end))]*obj.scale;
                    tmp = diff(iEdge)*obj.scale;
                    w1 = tmp(1:2:end);
                    w2 = tmp(2:2:end);

                    pitchMean = mean(p);
                    width1 = mean(w1);
                    width2 = mean(w2);

                    dutyCycle = mean(width1)/pitchMean;
                    varargout{1} = {p w1 w2 iEdge};
                elseif type == 2
                    p = [diff(iEdge(1:4:end)) diff(iEdge(3:4:end))]*obj.scale;
                    tmp = diff(iEdge)*obj.scale;
                    w1 = tmp(1:4:end);
                    w2 = tmp(2:4:end);
                    w3 = tmp(3:4:end);
                    w4 = tmp(4:4:end);

                    pitchMean = mean(p);
                    width1 = mean(w1)+mean(w2)+mean(w3);
                    width2 = mean(w4);

                    dutyCycle = mean(width1)/pitchMean;
                    varargout{1} = {p w1 w2 w3 w4 iEdge};
                end
            catch
                error('There should be multiples of 2 or 4 edge peaks for type 1 or 2 respectively');
            end
            
            if obj.showPlot
                y = da/range(da)-mean(da/range(da));        %normalize and center
                ya = abs(y);                                %absolute value
                t = ones(size(y))*max(ya)*threshold;
                
                y = y*size(obj.I,1)/3+size(obj.I,1)/3;
                ya = ya*size(obj.I,1)/3+size(obj.I,1)/3*2;
                t = t*size(obj.I,1)/3+size(obj.I,1)/3*2;
                figure
                colormap(gray)
                imagesc(obj.I)
                hold on
                plot(y,'r')
                plot(ya,'m');
                plot(t, 'c');
                hold off
                title(['pitch(' num2str(pitchMean) ')  dutyCylce(' num2str(dutyCycle) ')']);
            end
        end
        
        function [sdev, varargout] = getGratingLER(obj, varargin)
        %[sdev, sdevc] = getGratingLER(obj, peakThreshold = 0.8, binWidth = 10, binMinThreshold = 0.5)
        %
        %   Output:
        %       sdev    :   array of double     =   std of edge data
        %       sdevc 	:   array of double     =   tilt corrected std
        %       rang    :   array of double     =   range of data
        %       rangc   :   array of double     =   tilt corrected range
        %
        %   Input:
        %       obj             :   class instance
        %       peakThreshold   :   double          =   accept peaks larger
        %                                               than the threshold
        %       binWidth        :   integer         =   width of bin for
        %                                               clustering data points
        %       binMinThreshold :   double          =   threshold for
        %                                               acceptable number
        %                                               of counts in bin
        %
        %   Description:
        %   Determines the line edge roughness of a grating pattern,
        %   reported as the standard deviation of points that deviates from
        %   a straight line drawn along the edge of the grating.
        %
        %   If the edge of the grating isn't exactly vertical or aligned
        %   with the y access, then the reported LER will be wrong.  The
        %   corrected LER fits a line through the dataset and removes the
        %   slant from the data before calculating the LER.
        %
        %   Example:
        %   g = analyzePatterns();
        %   g = g.setShowPlot(true);
        %   g = g.loadImage('G200_F0_250kX.TIF');
        %   g = g.cropImage([214 71 896 811]);
        %   g = g.gaussianFilter(1);
        %   [a, b, c, d] = g.getGratingLER();
        %
        %   See also setShowPlot, cropImage, atfFilter, findCenter,
        %   isolateDots, rotateImage, alignX, setImage, removeData,
        %   getGrating
        
            peakThreshold = 0.8;
            binWidth = 10;
            binMinThreshold = 0.5;
            
            if nargin > 1
                peakThreshold = varargin{1};
            end
            if nargin > 2
                binWidth = varargin{2};
            end
            if nargin > 3
                binMinThreshold = varargin{3};
            end
            
            dimI = size(obj.I);
            
            peakMap = zeros(dimI);
            for i = 1:dimI(1)
                iP = obj.IPK.findPeaks(obj.I(i,:));
                peakMax = max(obj.I(i,iP));
                ind = iP(obj.I(i,iP)>peakMax*peakThreshold);
                peakMap(i,ind) = 1;
            end

            %Determine position of peaks
            [iy, ix] = ind2sub(size(obj.I),find(peakMap==1));

            %Find points that are along same column
            [n, xb] = hist(ix, max(ix)/binWidth);
            ibin = find(n>max(n)*binMinThreshold);

            cluster = cell(1,length(ibin));
            sdev = zeros([1, length(ibin)]);
            sdevc = sdev;
            rang = sdev;
            rangc = sdev;

            options = optimset('Display','off','TolFun',1e-6);
            for i = 1:length(ibin)
                a = ix > xb(ibin(i))-binWidth;
                b = ix < xb(ibin(i))+binWidth;
                cluster{i} = find(and(a,b));

                %Fit a line through the data
                x = ix(cluster{i});
                y = iy(cluster{i});
                f = lsqnonlin(@(h) (h(1)*y+h(2)-x), [1000 y(1)], [], [],options);
                
                %Remove tilt from data
                xc = x-f(1)*(y-mean(y));
                
                %Determine std of data and tilt corrected data
                sdev(i) = std(x)*obj.scale;
                sdevc(i) = std(xc)*obj.scale;
                rang(i) = range(x)*obj.scale;
                rangc(i) = range(xc)*obj.scale;
            end
            
            varargout = {sdevc rang rangc};
            
            if obj.showPlot
                figure
                colormap(gray)
                imagesc(obj.I)
                hold on
                for i = 1:length(cluster)
                    plot(ix(cluster{i}),iy(cluster{i}),'.')
                end
                hold off
                title(['LER(' num2str(mean(sdev)*obj.scale) ')  cLER(' num2str(mean(sdevc)*obj.scale) ')']);
            end
        end
        
        function [tang, bang] = getGratingAngle(obj, varargin)
        %[tang, bang] = getGratingAngle(obj, peakThreshold = 0.8, binWidth = 10, binMinThreshold = 0.5)
        %
        %   Output:
        %       sdev    :   array of double     =   std of edge data
        %       sdevc 	:   array of double     =   tilt corrected std
        %       rang    :   array of double     =   range of data
        %       rangc   :   array of double     =   tilt corrected range
        %
        %   Input:
        %       obj             :   class instance
        %       peakThreshold   :   double          =   accept peaks larger
        %                                               than the threshold
        %       binWidth        :   integer         =   width of bin for
        %                                               clustering data points
        %       binMinThreshold :   double          =   threshold for
        %                                               acceptable number
        %                                               of counts in bin
        %
        %   Description:
        %   Determines the line edge roughness of a grating pattern,
        %   reported as the standard deviation of points that deviates from
        %   a straight line drawn along the edge of the grating.
        %
        %   If the edge of the grating isn't exactly vertical or aligned
        %   with the y access, then the reported LER will be wrong.  The
        %   corrected LER fits a line through the dataset and removes the
        %   slant from the data before calculating the LER.
        %
        %   Example:
        %   g = analyzePatterns();
        %   g = g.setShowPlot(true);
        %   g = g.loadImage('G200_F0_250kX.TIF');
        %   g = g.cropImage([214 71 896 811]);
        %   g = g.gaussianFilter(1);
        %   [a, b, c, d] = g.getGratingLER();
        %
        %   See also setShowPlot, cropImage, atfFilter, findCenter,
        %   isolateDots, rotateImage, alignX, setImage, removeData,
        %   getGrating
        
            peakThreshold = 0.8;
            binWidth = 10;
            binMinThreshold = 0.5;
            deadZoneRegion = 0.2;
            
            if nargin > 1
                peakThreshold = varargin{1};
            end
            if nargin > 2
                binWidth = varargin{2};
            end
            if nargin > 3
                binMinThreshold = varargin{3};
            end
            
            dimI = size(obj.I);
            
            peakMap = zeros(dimI);
            for i = 1:dimI(1)
                iP = obj.IPK.findPeaks(obj.I(i,:));
                peakMax = max(obj.I(i,iP));
                ind = iP(obj.I(i,iP)>peakMax*peakThreshold);
                peakMap(i,ind) = 1;
            end

            %Determine position of peaks
            [iy, ix] = ind2sub(size(obj.I),find(peakMap==1));

            tix = ix(iy < dimI(1)/2-dimI(1)*deadZoneRegion);
            tiy = iy(iy < dimI(1)/2-dimI(1)*deadZoneRegion);
            bix = ix(iy > dimI(1)/2+dimI(1)*deadZoneRegion);
            biy = iy(iy > dimI(1)/2+dimI(1)*deadZoneRegion);
            
            %Find points that are along same column for top section
            [n, xb] = hist(tix, max(tix)/binWidth);
            ibin = find(n>max(n)*binMinThreshold);

            tcluster = cell(1,length(ibin));
            options = optimset('Display','off','TolFun',1e-6);
            tf = zeros([length(ibin) 2]);
            for i = 1:length(ibin)
                a = tix > xb(ibin(i))-binWidth;
                b = tix < xb(ibin(i))+binWidth;
                tcluster{i} = find(and(a,b));

                %Fit a line through the data
                tx = tix(tcluster{i});
                ty = tiy(tcluster{i});
                tf(i,:) = lsqnonlin(@(h) (h(1)*ty+h(2)-tx), [1000 ty(1)], [], [],options);
            end

            %Find points that are along same column for bottom section
            [n, xb] = hist(bix, max(bix)/binWidth);
            ibin = find(n>max(n)*binMinThreshold);

            bcluster = cell(1,length(ibin));
            options = optimset('Display','off','TolFun',1e-6);
            bf = zeros([length(ibin) 2]);
            for i = 1:length(ibin)
                a = bix > xb(ibin(i))-binWidth;
                b = bix < xb(ibin(i))+binWidth;
                bcluster{i} = find(and(a,b));

                %Fit a line through the data
                bx = bix(bcluster{i});
                by = biy(bcluster{i});
                bf(i,:) = lsqnonlin(@(h) (h(1)*by+h(2)-bx), [1000 by(1)], [], [],options);
            end

            tang = atan(1/mean(tf(:,1)));
            bang = atan(1/mean(bf(:,1)));
            
            if obj.showPlot
                figure
                colormap(gray)
                imagesc(obj.I)
                hold on
                for i = 1:length(tcluster)
                    plot(tix(tcluster{i}),tiy(tcluster{i}),'.')
                end
                for i = 1:length(bcluster)
                    plot(bix(bcluster{i}),biy(bcluster{i}),'.')
                end
                title(['TopAng(' num2str(tang,'%2.2f') ')  BotAng(' num2str(bang,'%2.2f') ')']);
            end
        end
    end
    
end