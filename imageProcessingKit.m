%imageProcessingKit Class
%
%This class contains a set of imaging processing functions that was
%developed for analyzing bit patterned media.  Bit patterned media is an
%array of magnetic disks.  SEM images of BPM or contact tester image of BPM
%can be analyzed semi-automatically using the functions from this kit.
%
%Long Chang, UH, March 17 2016
%
%See also gaussianFilter, atf, threshold, findPeaks, peakImg, findCluster,
%removeClusters, getNeighbors, neightborWindow, indexWindow, sortPolygon,
%insidePoly, rotateImage, findShift

classdef imageProcessingKit
    
    properties(GetAccess=public, SetAccess=protected, Hidden=false)
        %The get and set privileges are as shown above.
        
        %Number of turns
        hexIndex;
        
    end
    
    properties(GetAccess=public, SetAccess=protected, Hidden=true)
        %No Private Variables
    end
    
    properties(GetAccess=private, Constant=true, Hidden=true)
        %No Constants
    end
    
    methods (Access=public, Sealed=true)
        
        function obj = imageProcessingKit()
        %obj = imageProcessingKit()
        %
        %   Output:
        %       obj     	:   class instance
        %
        %   Description:
        %   Constructor for imageProcessingKit class.
     
            obj.hexIndex = [];
            
        end
        
        function gimg = gaussianFilter(~, img, sigma)
        %gimg = gaussianFilter(obj, img, sigma)
        %
        %   Output:
        %       gimg        :   2D matrix   =   filtered image
        %
        %   Input:
        %       obj         :   class instance
        %       img         :   2D matrix       =   image
        %       sigma       :   integer         =   amount of blur
        %
        %   Description:
        %   Filters the image by convolving it with a 2D gaussian function.
        %   The blur parameter is used as sigma when generating the
        %   gaussian function.
        %
        %   See also atf, threshold
            
            dim = size(img');
            x = (1:dim(1))-mean(1:dim(1));
            y = (1:dim(2))-mean(1:dim(2));
            [X, Y] = meshgrid(x,y);
            g = 1/(2*pi*sigma^2)*exp((-X.^2-Y.^2)/(2*sigma^2));
            fg = fft2(g);
            fimg = fft2(double(img));
            gimg = fftshift(ifft2(fg.*fimg));
        end

        function fimg = atf(~, img, radius)
        %fimg = atf(obj, img, radius)
        %
        %   Output:
        %       fimg        :   2D matrix   =   filtered image
        %
        %   Input:
        %       obj         :   class instance
        %       img         :   2D matrix       =   image
        %       radius      :   integer         =   radius of disk filter
        %
        %   Description:
        %   The adaptive threshold filter is very useful for removing
        %   background noise.  A background image is acquired by processing
        %   the image with a disk filter using a large blur.  Then the
        %   background image is subtracted from the raw image.
        %
        %   See also gaussianFilter, threshold

            imgf = imfilter(img,fspecial('disk',radius),'replicate');
            fimg = img-imgf;
        end
        
        function img = threshold(~, img, threshold)
        %img = threshold(~, img, threshold)
        %
        %   Output:
        %       img         :   2D matrix       =   processed image
        %
        %   Input:
        %       obj         :   class instance
        %       img         :   2D matrix       =   image
        %       threshold   :   double [0 to 1] =   threshold value
        %
        %   Description:
        %   Sets any pixel in the image that is below the threshold to
        %   zero.
        %       threshold = (maxPixelValue - minPixelValue)*threshold
        %
        %   See also gaussianFilter, atf
        
            img(img<(max(max(img))-min(min(img)))*threshold + min(min(img))) = 0;
        end
        
        function [peaks, troughs] = findPeaks(~, M)
        %[peaks, troughs] = findPeaks(obj, M)
        %
        %   Output:
        %       peaks       :   array   =   index of peaks
        %       troughs     :   array   =   index of trough
        %
        %   Input:
        %       obj         :   class instance
        %       img         :   2D matrix       =   image
        %                   :   array           =   data
        %
        %   Description:
        %   Find the peaks in an image or a curve according to the
        %   dimensions of the data that is passed.
        %
        %   See also peakImg

            dim = size(M);
            if min(dim) ~= 1

                %Obtain the differential between rows and colums
                dx = M(:,2:end) - M(:,1:end-1);
                dy = M(2:end,:) - M(1:end-1,:);

                %Convert the differential matrix to a binary matrix
                Ax = dx>0;
                Ay = dy>0;

                %Peaks are located where the gradient changes sign from + to -
                %   The matrix is padded with zeros because each differential
                %   reduces the matrix by a row/col
                Bx = [zeros(dim(1),1) zeros(dim(1),1) Ax(:,1:end-1)-Ax(:,2:end)];
                By = [zeros(1, dim(2)); zeros(1, dim(2)); Ay(1:end-1,:)-Ay(2:end,:)];

                %Locate peaks
                peaks = find(Bx+By==2);

                %Locate troughs
                troughs = find(Bx+By == -2);
            else
                M = reshape(M,1,numel(M));
                dx = M(2:end) - M(1:end-1);
                A = dx>0;
                B = [0 A(1:end-1)-A(2:end)];

                peaks = find(B == 1);
                troughs = find(B == -1);
            end
        end
    
        function img = peakImg(obj, M, type)
        %img = peakImg(obj, M, type)
        %
        %   Output:
        %       img         :   2D matrix       =   image where each point
        %                                           is a peak or a trough
        %
        %   Input:
        %       M           :   2D matrix           =   image
        %       type        :   integer [1,2,3,4]   =   1   =   peak only
        %                                               2   =   trough only
        %                                               3   =   peaks and trough are both 1
        %                                               4   =   peaks are 1 and troughs are -1
        %                                               
        %
        %
        %   Description:
        %   This function finds all peaks and troughs in an image and then
        %   marks them in the output image
        %
        %   See also findPeaks
        
            [peaks, troughs] = obj.findPeaks(M);

            dim = size(M);
            pimg = zeros(dim);
            pimg(peaks) = 1;
            timg = zeros(dim);
            timg(troughs) = 1;

            switch type
                case 1
                    img = pimg;
                case 2
                    img = timg;
                case 3
                    img = pimg + timg;
                case 4
                    img = pimg - timg;
                otherwise
                    img = pimg + timg;
            end
        end
        
        function neigh = findCluster(obj, px , py, r)
        %neigh = findClusters(obj, px, py, r)
        %
        %   Output:
        %       neigh       :   cell    =   a list of neighboring points
        %
        %   Input:
        %       obj         :   class instance
        %       px          :   array   =   x points
        %       py          :   array   =   y points
        %       r           :   double  =   radius of cluster
        %
        %   Description:
        %   Find the peaks in an image or a curve according to the
        %   dimensions of the data that is passed.
        %
        %   See also findCluster, removeClusters, getNeighbors,
        %   neighborWindow, indexWindow, sortPolygon, insidePoly,
        %   rotateImage, findShift
        
            %Remove points that are within r distance
            polygon = obj.neighborWindow(1,1);
            polygon(1,:) = [];
            polygon = obj.sortPolygon(polygon);
            polygon = polygon*r;

            neigh = obj.getNeighbors(px,py,polygon);
            
        end

        function [px, py] = removeClusters(~, px,py,r)
        %[px py] = removeClusters(obj, px,py,r)
        %
        %   Output:    
        %       px  :   array   =   x points with clusters removed 
        %     	py  :   array   =   y points with clusters removed
        %
        %   Input:     
        %       obj         :   class instance
        %       px  :   array   =   x points
        %      	py  :   array   =   y points
        %    	r   :   double  =   radius of cluster
        %
        %   Description:  
        %   Assumes a square cluster where 2*r is the length of each edge
        %
        %   Note:  
        %   The results of this function is dependent on the sequence
        %   of the data points.  For example:
        %   x = [1 5 9] and r = 5
        %       if we start with x(1)
        %           result = [3 9]
        %       if we start with x(2)
        %           result = [5]
        %       if we start with x(3)
        %           result = [1 7]
        %
        %   See also findCluster, removeClusters, getNeighbors,
        %   neighborWindow, indexWindow, sortPolygon, insidePoly,
        %   rotateImage, findShift
        
            i = 1;
            while i <= length(px)
                    %find points around reference
                    ind = find(px>=px(i)-r & px<=px(i)+r & py>=py(i)-r & py<=py(i)+r);

                    %find points that are within r of the reference
                    z = abs(px(ind)-px(i)+(py(ind)-py(i))*1i);
                    ind = ind(z<=r);

                    if(~isempty(ind))
                        %find centroid of cluster
                        px(ind(1)) = mean(px(ind));
                        py(ind(1)) = mean(py(ind));
                        %remove all points in cluster
                        if(length(ind)>1)
                            px(ind(2:end)) = [];
                            py(ind(2:end)) = [];
                        end
                    end
                i = i+1;
            end
        end
        
        function neighbor = getNeighbors(obj, px, py, polygon)
        %neighbor = getNeighbors(obj, px,py,polygon)
        %
        %   Output:
        %       neighbor    :   cell    =   cell containing indices of all neigbors
        %
        %   Input:
        %       obj         :   class instance
        %       px          :   array           =   x position
        %       py          :   array           =   y position
        %       polygon     :   array           =   polygon to detect neighbors
        %
        %   Description:
        %   Draws a polygon around every point in [px,py] and determines
        %   neighboring points be deciding if that point is inside the
        %   polygon.
        %
        %   Example:
        %   %create polygon
        %   polygon = neighborWindow(1,1);
        %   polygon = polygon(find(~(polygon(:,1)==0 & polygon(:,2)==0)),:);
        %   polygon = sortPolygon(polygon);
        %   polygon = polygon*vm(2)*1.5;        %scale
        %   [polygon(:,1) polygon(:,2)] = rotatePoints(polygon(:,1) ...
        %   	,polygon(:,2),vm(1));           %rotate
        % 
        %   neigh = getNeighbors(px,py,polygon);
        %
        %   See also findCluster, removeClusters, getNeighbors,
        %   neighborWindow, indexWindow, sortPolygon, insidePoly,
        %   rotateImage, findShift
        
            %preallocate memory
            neighbor = cell([length(px) 1]);

            for a=1:length(px);
                nw(:,1) = polygon(:,1)+px(a);
                nw(:,2) = polygon(:,2)+py(a);

                %find local points
                bound = max(abs(polygon(:,1)+polygon(:,2)*1i));
                pi = find(px<px(a)+bound & px>px(a)-bound & py<py(a)+bound & py>py(a)-bound);

                %determine if points are inside the polygon
                nindex = 0;
                for k = 1:length(pi)
                    if(obj.insidePoly([px(pi(k)) py(pi(k))],nw))
                        nindex = [nindex pi(k)];
                    end
                end
                nindex = nindex(2:end);
                neighbor{a} = [nindex(nindex==a) nindex(nindex~=a)];
            end
        end
        
        function nindex = neighborWindow(obj, n, varargin)
        %nindex = neighborWindow(obj, n, type = 0)
        %
        %   Output:
        %       nindex  :   array   =   neighbor index window
        %
        %   Input:
        %       obj     :   class instance
        %       n       :   integer     =   n closest neighbors
        %       *type   =   1   =   hexagonal (n <= 248)
        %               =   0   =   square
        %
        %   Description:
        %   Generates an array of neighboring points according to the type
        %   of lattice specified
        %
        %   See also findCluster, removeClusters, getNeighbors,
        %   neighborWindow, indexWindow, sortPolygon, insidePoly,
        %   rotateImage, findShift
        
            type = 0;

            if(nargin>2)
                type = varargin{1};
            end

            if (type == 1)
                if isempty(obj.hexIndex)
                    obj = obj.genHexIndex();
                end
                index = obj.hexIndex;
            else
                index = obj.indexWindow(n);
            end

            z = abs(index(:,1)+index(:,2)*1i);
            d = sort(obj.removeClusters(z,zeros(size(z)),.01));
            r = d(n+1);
            dr = (d(n+2)-d(n+1))/2;
            ind = find(z<r+dr);
            nindex = [index(ind,1) index(ind,2)];
        end
        
        function obj = genHexIndex(obj)
        %obj = genHexIndex(obj)
        %
        %   Output:
        %       obj         :   class instance
        %
        %   Input:
        %       obj         :   class instance
        %
        %   Description:
        %   Generates the hex indices for neighborWindow.  A larger hex
        %   window than what is currently provided can be done using the
        %   same code, but it would take a while.  All you need to do is to
        %   change k = 1:n where n is larger than 4.  This algorithm is
        %   very bad.
        
            h = [0 -.5 .5 -.5 .5 -1 1; 0 -.866 -.866 .866 .866 0 0]';

            for k = 1:4
                H = [0 0];
                for i = 1:length(h)
                    H = [H; repmat(h(i,:),length(h),1)+h];
                end

                [hx, hy] = obj.removeClusters(H(:,1),H(:,2),.01);
                h = [hx hy];    
            end
            obj.hexIndex = h;
        end
        
        function iWindow = indexWindow(~, max_shift)
        %iWindow = indexWindow(obj, max_shift)
        %
        %	Ouput:  
        %       iWindow     : 	[row col]   =   row and col of window 
        %
        %   Input:  
        %       obj         :   class instance
        %       max_shift	:   integer     =   maximum shift
        %
        %   Description:
        %  	This function produces a set of integer coordinates
        %  	of a square with side length of 2*max_shift
        %
        %   Example:
        %  	iWindow = indexWindow(1)
        %  	iWindow = [-1 -1; -1 0; -1 1; 0 -1; 0 0; 0 1; 1 -1; 1 0; 1 1]
        %
        %   See also findCluster, removeClusters, getNeighbors,
        %   neighborWindow, indexWindow, sortPolygon, insidePoly,
        %   rotateImage, findShift

            shift_map_row = repmat((-max_shift:max_shift),max_shift*2+1,1);
            shift_map_row = reshape(shift_map_row,numel(shift_map_row),1);
            shift_map_col = repmat((-max_shift:max_shift)',max_shift*2+1,1);
            iWindow = [shift_map_row shift_map_col];
        end       
        
        function points = sortPolygon(~, points)
        %points = sortPolygon(obj, points)
        %
        %   Output:
        %       sortedPoints	:   array   =  sorted points
        %
        %   Input:
        %       obj         :   class instance
        %       points  =   [x y]
        %
        %   Description:
        %   The points are sorted such that the distance between each point
        %   is minimized
        %
        %   See also findCluster, removeClusters, getNeighbors,
        %   neighborWindow, indexWindow, sortPolygon, insidePoly,
        %   rotateImage, findShift

            for i = 2:length(points)
                d = abs((points(i-1,1)-points(i:end,1))+(points(i-1,2)-points(i:end,2))*1i);
                ind = find(d == min(d));
                temp = points(i,:);
                points(i,:) = points(ind(1)+(i-1),:);
                points(ind(1)+(i-1),:) = temp;
            end
        end

        function result = insidePoly(~, point, polygon, varargin)
        %result = insidePoly(obj, point, polygon, debug = false)
        %
        %   Ouput:  
        %       result      :   integer     =   1   = point is inside
        %                                   =   0 	= point is outside
        %   Input:  
        %       obj         :   class instance
        %       point       :   [x y]   =   test point
        %    	polygon 	:   [x y]   =   coordinates of polygon
        %   	debug       :   boolean	=   plots point and polygon for inspection
        %
        %   Procedure:
        %       1)  Shift the space such that the test point is at the origin
        %       2)  Determine the number of edges of the polygon cross the positive
        %       x axis
        %       3)  If crosses is odd, the test point is inside the polygon.
        %
        %   See also findCluster, removeClusters, getNeighbors,
        %   neighborWindow, indexWindow, sortPolygon, insidePoly,
        %   rotateImage, findShift
        
            debug = false;
            if(nargin>3)
                debug = varargin{1};
            end

            %reshapes the input matrices for processing
            dim = size(point);
            if(dim(2)<dim(1))
                point = point';
            end
            point = repmat(point,length(polygon),1);
            dim = size(polygon);
            if(dim(1)<dim(2))
                polygon = polygon';
            end

            %move the coordinate space about the point
            polygon = polygon - point;

            %insures that no vertex lies on the x-axis, see Note
            if(find(polygon(:,2)==0))
                polygon(:,2) = polygon(:,2)+.0001;
            end

            %shift the polygon
            polygon = [polygon;polygon(1,:)];    

            ap = atan2(polygon(:,2),polygon(:,1));
            %count number of times an edge crosses the positive x axis
            counter = 0;
            for i = 1:length(ap)-1
                if((ap(i)*ap(i+1)<0) && (abs(ap(i)-ap(i+1))<pi))
                    counter = counter+1;
                end
            end

            %An odd number of crossing means that the point is inside the polygon
            result = mod(counter,2);

            if debug
                patch([polygon(:,1) polygon(:,1)/100],[polygon(:,2) polygon(:,2)/100],2);
            end
        end
        
        function I = rotateImage(~, M, varargin)
        %I = rotateImage(obj, M, ang = 90)
        %
        %   Output:
        %       I       :   2D matrix   =   rotated image
        %
        %   Input:
        %       M       :   2D matrix   =   image
        %       ang     :   double      =   angle in degrees
        %
        %   Description:
        %   Rotates an image M by ang degrees in a counterclockwise 
        %   direction around its center point.  See documentation for 
        %   imrotate for more information.
        %
        %   See also findCluster, removeClusters, getNeighbors,
        %   neighborWindow, indexWindow, sortPolygon, insidePoly,
        %   rotateImage, findShift

            if nargin > 2
                ang = varargin{1};
            end
            I = imrotate(M,ang,'bilinear');

        end
        
        function I = flatten(~, M, varargin)
        %I = flatten(obj, M, type = 'x', order = 1)
        %
        %   Output:
        %       I       :   2D matrix   =   flattened image
        %
        %   Input:
        %       M       :   2D matrix   =   image
        %       type    :   string      =   x   =   flatten along x
        %                                   y   =   flatten along y
        %                                   xy  =   flatten x then y
        %                                   yx  =   flatten y then x
        %       order   :   integer     =   order of polynomial fit
        %                               =   0, 1 or 2
        %
        %   Description:
        %   Flatten the image by subtracting a fitted polynomial one line
        %   at a time
        %
        %   See also findCluster, removeClusters, getNeighbors,
        %   neighborWindow, indexWindow, sortPolygon, insidePoly,
        %   rotateImage, findShift

            type = 'x';
            order = 1;
            if nargin > 2
                type = varargin{1};
            end
            if nargin > 3
                order = varargin{2};
                if order < 0
                    order = 0;
                end
                if order > 2
                    order = 2;
                end
            end
            
            %flatten along x
            if (strcmp(type,'x') || strcmp(type,'xy'))
                x = (1:numel(M(1,:)));
                for i = 1:numel(M(:,1))
                    parm = polyfit(x,M(i,:),order);
                    y = polyval(parm,x);
                    M(i,:) = M(i,:)-y;
                end
            end
            
            %flatten along y
            if (strcmp(type,'y') || strcmp(type,'xy') || strcmp(type,'yx'))
                x = (1:numel(M(:,1)))';
                for i = 1:numel(M(1,:))
                    parm = polyfit(x,M(:,i),order);
                    y = polyval(parm,x);
                    M(:,i) = M(:,i)-y;
                end
            end
            
            %flatten along x
            if strcmp(type,'yx')
                x = (1:numel(M(1,:)));
                for i = 1:numel(M(:,1))
                    parm = polyfit(x,M(i,:),order);
                    y = polyval(parm,x);
                    M(i,:) = M(i,:)-y;
                end
            end
            
            I = M;
            
        end

        function [shift, varargout] = findShift(obj, img1, img2, max_shift)
        %[shift corr*] = findDrift(img1, img2, max_shift)
        %   Output: 
        %       shift   :   [x y]   =   The shift between img1 and img2
        %                               img1+shift = img2
        %       corr*   :   array   =   correlation data
        %
        %   Input:  
        %       img1        :   2D matrix   =   image 1
        %    	img2        :   2D matrix   =   image 2
        %     	max_drift   :   integer     =   max expected shift
        %
        %   Description:
        %   This function assumes the two input images are the identical,
        %   but are only different because of a translation/shift.  The
        %   amount of shift is determined by:
        %   	1)  Determine largest window possible for calculation
        %      	2)  Make all dots the same polarity
        %     	3)  Determine correlation between the 2 images
        %      	4)  Shift second image and repeat step 3
        %      	5)  Find the shifted value with largest correlation
        %
        %   This function works for 1D array data as well.
        %
        %   Example:
        %   shift = findDrift(img1,img2,5);
        %       img1+shift = img2;
        %
        %   See also findCluster, removeClusters, getNeighbors,
        %   neighborWindow, indexWindow, sortPolygon, insidePoly,
        %   rotateImage, findShift
            
            if ~any(size(img1) == [1 1])
                %transpose so that x = rows and y = columns
                img1 = img1';
                img2 = img2';

                dim = size(img1);
                %index of the center of image
                center = ceil(dim/2);

                %estimate largest window size. window size is smaller if the expected
                %   drift is larger.
                window = floor((dim-max_shift*2-1)/2);

                %Create the index matrix for a window centered on the image
                i_x = repmat(center(1)-window(1):center(1)+window(1),window(2)*2+1,1);
                i_x = reshape(i_x,numel(i_x),1);
                i_y = repmat((center(2)-window(2):center(2)+window(2))',window(1)*2+1,1);
                i_center = sub2ind(dim,i_x,i_y);

                %Create a map that outlines all the possible shift within the expected
                %   limit
                shift_map = obj.indexWindow(max_shift);

                %Find the correlation between the 2 images
                r = zeros(1,length(shift_map));

                %Correlate images
                for i = 1:length(shift_map)
                    i_shift = sub2ind(dim,i_x+shift_map(i,1),i_y+shift_map(i,2));
                    r(i) = sum(img1(i_center).*img2(i_shift));
                end

                %Find the maximum correlation and related shift
                shift = shift_map(r==max(r),:);
                shift = shift(1,:);

                varargout{1} = reshape(r,max_shift*2+1,max_shift*2+1);
            else
                %generate indices
                index = max_shift+1:numel(img1)-max_shift;
                iShift = -max_shift:max_shift;
                
                r = zeros(1, length(iShift));
                for i = 1:length(iShift);
                    r(i) = sum(img1(index).*img2(index+iShift(i)));
                end
                
                %Find the maximum correlation
                shift = iShift(r == max(r));
                
                varargout{1} = r;
            end
        end
    end
end