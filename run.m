%{
Recursion Erosion Program
Name: Saurabh Goel
Roll No.: 10CH10061

Include imoverlay.m in the root directory
Type test in matlab to run the program
Images are present in the /IMAGES folder
Output is present in the /RESULTS folder
%}

% ------------------------------------------------------------------------------------------------------------
% FUNCTION to display the menu
function [] = test()
    result = 1;
    disp('************************************************************');
    disp('PARTICLE SIZE DISTRIBUTION USING RECURSION EROSION ALGORITHM');
    disp('************************************************************');
    disp('Select the folder containing the images:')

    global images magnification final_distribution;
    magnification = 50;
    final_distribution = [];
    
    field1 = 'original_image';      value1 = {};
    field2 = 'process_image';       value2 = {};
    field3 = 'final_image';         value3 = {};
    images = struct(field1, value1, field2, value2, field3, value3);

    warning('off','all');
    while result~=0
        read();
        for k=1:length(images)
            intialize(images(k).original_image, k);
        end
        figure, hist(final_distribution, 16);
        xlabel('Size of Particles (nm)');
        ylabel('Number of Particles');
        title('Combined Particle Size Distribution');
        
        prompt = '\nDo you want to continue (1=yes, 0=no): ';
        result = input(prompt);
    end
end


% ------------------------------------------------------------------------------------------------------------
% FUNCTION to read images using the image browser and return the image and magnification
function read()

    % Accessing the global variables
    global images magnification;

    % Selecting the folder to read images from
    pathname = uigetdir('','Choose a folder');
    if (pathname)
        
        dirData = dir(pathname);
        dirIndex = [dirData.isdir];
        fileList = {dirData(~dirIndex).name};
        
        % Reading the images and storing them in cell array
        for k=1:length(fileList)
            image = imread([strcat(pathname,'\'), fileList{k}]);
            images(length(images) + 1).original_image = image;
        end
        
        % Asking the Magnification from the user
        prompt = 'Choose the magnification (50/100): ';
        magnification = input(prompt);
    
    else
        error('Folder not selected.');
    end
end


% ------------------------------------------------------------------------------------------------------------
% FUNCTION to initialize (set up the variables and call functions)
function intialize(image, n)


    % Initializing the Global Variables
    global images process_image magnification erosion_radius bwarea_radius image_size final_distribution;

    field1 = 'segment';     value1 = {};
    field2 = 'diameter';    value2 = {};
    field3 = 'centroid_x';  value3 = {};
    field4 = 'centroid_y';  value4 = {};
    process_image = struct(field1, value1, field2, value2, field3, value3, field4 , value4);

    % Printing status bars
    fprintf(strcat('\nProcessing (Image ',num2str(n),') '));

     % Computing the size of the image
    image_size = size(image);

    % Display the Original Image
    figure, imshow(image);

    % Intializing the erosion radius and bwarea depending on magnification
    if (magnification==50)
        erosion_radius = 50;
        bwarea_radius = 60;
    else
        erosion_radius = 20;
        bwarea_radius = 20;
    end
    
    % Initializing the Final Image and Boundary Image as empty images of the same size as actual image
    final_image = zeros(image_size(1),image_size(2));
    final_image = double(cat(3, final_image, final_image, final_image));
    boundary_image = zeros(image_size(1),image_size(2));
    
    % Calling the watershed algorithm on the original image
    s = watershed_algo(image);
    
    % Traversing the watershed segments 
    for k=2:length(s)
        % Pre Erosion Dilation - Dilating and eroding to remove irregular elements
        se = strel('disk',erosion_radius);                      
        s(k).Image = imerode(s(k).Image,se);
        s(k).Image = imdilate(s(k).Image, se);
        
        % Calling the recursion function
        recursion(build_image(s(k)), 0);

    end
    
    % Calling the Module - Removing Big Objects
    module_remove_big_objects();

    % Calling the Module - Remove Repeating Objects
    module_adjust_repeating();

    % Traversing the processed segments
    for k=1:length(process_image)
        
        % Creating the Boundaries
        se = strel('disk', 2);
        boundary_mask = imdilate(process_image(k).segment, se) - process_image(k).segment; 

        % Overlaying the individual segments to get the final image
        final_image = imoverlay(final_image, process_image(k).segment, [rand(1,1) rand(1,1) rand(1,1)]);
      
        % Overlaying the boundaries on final image for a better view
        final_image = imoverlay(final_image, boundary_mask, [1 1 1]);
        
    end

    % Storing the segments
    images(length(images)).process_image = process_image;

    % Storing the final image
    images(length(images)).final_image = final_image;

    % Overlaying the results on original image
    final_image = imfuse(final_image, image,'blend','Scaling','joint');
    
    % Diplsaying the final image
    figure, imshow(final_image);
    
    % Overlaying the Centroids on the final image for a better view
    hold on
    for k = 1:length(process_image)
        plot([process_image(k).centroid_x],[process_image(k).centroid_y],'ro');
    end
    hold off;
    
    % Computing the actual diameter(nm) using the pixel ratio depending on the magnification
    if (magnification==50)
        distribution = ([process_image(1,:).diameter])*50/400;
    else
        distribution = ([process_image(1,:).diameter])*100/360;
    end

    final_distribution = [final_distribution distribution];

    % Plotting the Particle Size Distribution using the array of diameters
    figure, hist(distribution, 10);
    xlabel('Size of Particles (nm)');
    ylabel('Number of Particles');
    title(strcat('Particle Size Distribution (Image ',num2str(n),')'));
end


% ------------------------------------------------------------------------------------------------------------
% FUNCTION that inputs an image and computes the watershed transform and returns the labelled image variable s
function s = watershed_algo(image)

    global erosion_radius bwarea_radius;

    % Converting image to Grayscale
    gray_image = rgb2gray(image);
    
    % Computing the graythreshold
    thresh = graythresh(gray_image);
    
    % Converting the image to Black and White
    bw_image = im2bw(gray_image,thresh);
    
    % Remove small objects from the image
    bw_image = bwareaopen(bw_image, bwarea_radius);    
        
    % Computed the Extended-maxima transform
    mask = imextendedmax(imcomplement(gray_image), 60);
    
    % Morphologically Closing the image
    mask = imclose(mask, ones(5,5));
    
    % Filling Image regions and holes
    mask = imfill(mask, 'holes');
    
    % Removing small areas from the image
    mask = bwareaopen(mask, 40);
    
    % Morphological Reconstruction to impose regional minima
    I_mod = imimposemin(gray_image, bw_image | mask);
   
    % Computing the Watershed Transform
    L = watershed(I_mod,8);
    
    % Compute the region properties of the watershed segments
    s = regionprops(L, 'Image', 'Centroid','BoundingBox', 'EquivDiameter');

    % Displaying the watershed image
    figure, imshow(label2rgb(L));
end


% ------------------------------------------------------------------------------------------------------------
% FUNCTION that imputs a labelled image variable (s) and returns the complete image
function image = build_image(s)
    
    global image_size;

    fprintf('.');
    % Computing the bounding coordinates of the segment of the image
    a = double(uint32(s.BoundingBox(2)));                                                       
    b = double(uint32(s.BoundingBox(1)));
    c = double(uint32(s.BoundingBox(4)));
    d = double(uint32(s.BoundingBox(3)));
    
    % Case when the segment touches the boundary
    if(a==1 || b==1 || (a+c)>image_size(1) || (b+d)>image_size(1))
        
        % Removing the segment by making it blank (zero)
        image = zeros(image_size(1), image_size(2));                               
    
    % Case when the segment does not touch the boundary
    else
        
        % Build the image using the coordinates
        image = padarray(s.Image, [a b], 'pre');                
        a = image_size(1) - a - c;
        b = image_size(1) - b - d;
        image = padarray(image, [a b], 'post');
    end
end


% ------------------------------------------------------------------------------------------------------------
% FUNCTION that removes outlier objects or objects that are oversized (using outlier detection method)
function module_remove_big_objects()
    
    global process_image;
    
    x = [process_image(1,:).diameter];
    index = [];
    
    % Storing the size of the array
    Nx = size(x,1);

    % compute the median
    medianx = median(x);

    % STEP 1 - rank the data
    y = sort(x);

    % compute 25th percentile (first quartile)
    Q(1) = median(y(find(y<median(y))));

    % compute 50th percentile (second quartile)
    Q(2) = median(y);

    % compute 75th percentile (third quartile)
    Q(3) = median(y(find(y>median(y))));

    % compute Interquartile Range (IQR)
    IQR = Q(3)-Q(1);

    % Get the index which have outliers
    for k=1:length(x)
        if(x(k)<(Q(1)-3*IQR) | x(k)>(Q(1)+5*IQR))
            index = [index k];
        end
    end

    % Remove the segments that correspond to outliers
    if(~isempty(index))
        for k=1:length(index)
            index(k) = index(k) - k + 1;
            process_image(index(k)) = [];
        end
    end
end


% ------------------------------------------------------------------------------------------------------------
% FUNCTION that combines the repeating units into single particle
function module_adjust_repeating()

    global process_image image_size;
    
    arr1 = [];
    arr2 = [];
    final = {};
    
    % Identify particles that are repeating
    for k=1:(length(process_image)-1)
        
        for l=(k+1):length(process_image)
            
            % Define the position of the centroids as coordinates
            X = [process_image(k).centroid_x process_image(k).centroid_y ; process_image(l).centroid_x process_image(l).centroid_y];
            
            % Calculate the avergage equivalent distance between 2 particles using diameter
            equ_dia = (process_image(k).diameter + process_image(l).diameter)/2;
        
            % Check how close are the particles 
            if(pdist(X,'euclidean') < equ_dia)
        
                % Compute the areas of the two particles
                a1 = bwarea(process_image(k).segment);
                a2 = bwarea(process_image(l).segment);
                
                % Compute the amount of overlapping        
                diff = immultiply(process_image(k).segment,process_image(l).segment);
                a3 = bwarea(diff);

                % Identify particles that are heavily overlapped        
                if((a3/a1)>0.75 && (a3/a2)>0.75)
        
                    arr1 = [arr1 k];
                    arr2 = [arr2 l];
        
                end
            end
        end
    end

    % Form sets of segments to club together into a single particle
    k = 1;
    while(k <= length(arr1))
        
        current = [arr1(k) arr2(k)];
        if(k==length(arr1))
            l = k;
        else
            l = k + 1;
        end

        while (l <= length(arr1))
            if(ismember(arr1(l), current))
                current = [current arr2(l)];
                arr1(l) = [];
                arr2(l) = [];
            else
                if (ismember(arr2(l), current))
                    current = [current arr1(l)];
                    arr1(l) = [];
                    arr2(l) = [];
                else
                    l = l + 1;
                end 
            end
        end
        final{length(final) + 1} = unique(current);
        k = k + 1;
    end

    % Clubbing the specified elements together(from final cell)
    for k=1:length(final)
        
        for l=1:length(final{k})
            if (l==1)
                image = process_image(final{k}(l)).segment; 
            else 
                % Adding the segments
                image = process_image(final{k}(l)).segment + image;
            end
        end
        
        image = logical(image);

        % Labelling the image
        label_image = bwlabel(image);
    
        % Computing the region properties of the segment
        s = regionprops(label_image, 'Centroid', 'EquivDiameter');

        current = length(process_image) + 1;
 
        % Storing the segment
        process_image(current).segment = image; 

        % Storing the centroids of the segment
        process_image(current).centroid_x = s(1).Centroid(1); 
        process_image(current).centroid_y = s(1).Centroid(2); 
        
        % Storing the diameter of the segment
        process_image(current).diameter = s(1).EquivDiameter; 
    end

    % Remove the segments from the image
    counter = 0;
    for k=1:length(final)
        
        for l=1:length(final{k})
            final{k}(l) = final{k}(l) - counter;
            process_image(final{k}(l)) = [];
            counter = counter + 1;
        end
    end
end


% ------------------------------------------------------------------------------------------------------------
% FUNCTION to traverse the tree using recursion algorithm
function recursion(image, rad_o)
    
    % Accesing the global variables in this function
    global process_image;
    
    % Initialize the erosion radius variables
    rad_min = 10; 
    rad_max = 250;
    rad_interval = 10;
    rad = rad_min;
    
    % Initializing the flag variable
    flag = 1;
    
    % Building the structuring element
    se = strel('disk',rad);
    
    % Initialize eroded image with the original image
    eroded_image = image;
    
    while (rad < rad_max && flag==1)                                 
        
        % Eroding the image
        eroded_image = imerode(eroded_image,se);
        
        % Labelling the image
        label_image = bwlabel(eroded_image);
        
        % Computing the region properties of the segment
        s = regionprops(label_image, 'Image', 'Centroid','BoundingBox', 'EquivDiameter');
        
        % If multiple segments present then apply recursion on all segments
        if(length(s)>1)
            for k = 1:length(s)
                recursion(build_image(s(k)), rad_o+rad);
                flag = 0;
            end
        end
        
        % If no segments are present then the last biggest single element is the result
        if(length(s)==0)
            
            % Stopping the iteration
            flag= 0;
            
            % Structuring element equal to the amount of total erosion
            se = strel('disk', rad_o);
            
            % Dilating the image by the starting radius
            image = imdilate(image, se);
            
            s = regionprops(image, 'Centroid', 'EquivDiameter');
            
            if(length(s)>0)

                current = length(process_image) + 1;

                % Storing the segment
                process_image(current).segment = image; 

                % Storing the centroids of the segment
                process_image(current).centroid_x = s(1).Centroid(1); 
                process_image(current).centroid_y = s(1).Centroid(2); 
                
                % Storing the diameter of the segment
                process_image(current).diameter = s(1).EquivDiameter; 
               
            end
        end
        
        % If regionprops contains only one segment then continue erosion using a bigger radius
        rad = rad + rad_interval;
    end
end


%{
Variable Directory - 
s           : Labelling Variable containing the different segments of an image as
              an array (with regionprops properties).
mag         : Magnification of image (50 or 100)
diameter    : Array containing the diameters of particles detected
centroid_x  : Array containing the x coordinates of the centroids
centroid_y  : Array containing the y coordinated of the centroids
rad_o       : The total erosion that has been applied on a segment

%}
