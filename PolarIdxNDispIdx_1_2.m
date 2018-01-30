function PolarIdxNDispIdx_1_2()

%% this function calculates polarization index and dispersion index
%
%  http://www.cell.com/action/showMethods?pii=S2211-1247%2811%2900019-2
%  see reference for details
%
%  PI or polarization index is calculated as:
%  sqrt( centroid_of_molecule - centroid_of_cytoplasm )^2 / Radius_of_gyration
%
%  Rg or Radius of gyration is caclculated by root-mean-square distance of all
%  pixels of the cell from the centroids (centroid as axis) to normalize
%  cell shape
%
%  DI or dispersion index is calculated as:
%  variance_of_molecule / second_moment_of_cytoplasm
%  (i.e., variance of molecule / variance of cytoplasm)
%
%  color of perim has changed to magenta, now is RGB tiff image

clear; close all force; clc;

folderName = uigetdir('F:\matlab');
cd(folderName);
fileList = dir(fullfile(folderName, '*.tif'));
numFiles = length(fileList);
chkDir = 1;
se = strel('disk', 3);
xlsLine = 0;

for cntFiles = 1 : numFiles % all files within a folder will be analyzed
    
    cntCell = 1;
    multipleCell = true;
    initMask = true;
    
    % Read Images ---------------------
    FileName = fileList(cntFiles).name; % first file starts with 'GFP' = cluster
    cellImage = imread(FileName);

    imgInfo = imfinfo(FileName);
    imgDepth = imgInfo.BitDepth;

    if imgDepth ~= 8
        cellImage = uint8(cellImage / 16); % NIS file is in 12-bit - need to be divided by 2^4 = 16 to be an 8-bit file
    else
        cellImage = uint8(cellImage); % if is 8 bit then use as is
    end

    disp(strcat('analyzing file name = ', FileName)); % check console for appropriate files

    % Overall cell area is detected with actin
    filtImage = medfilt2(cellImage);
    filtImage = imgaussfilt(filtImage);
    
    % mask for exclusion
    baseMask = true(size(filtImage));
    
    % analysis starts
    while multipleCell
        
        if initMask
            adjImg = imadjust(uint8(filtImage));
            initMask = false;
        else
            adjImg = adjImg .* uint8(baseMask);
        end
        
        % isolate each cell manually
        
        
        CellMask = roipoly(adjImg); 

        close all;

        % select cell area
        ImageSelected = adjImg .* uint8(CellMask); 
        seg_I = imquantize(ImageSelected, 5);
        CellMask = seg_I > 1;
        CellMask = imdilate(CellMask, se);
        CellMask = imfill(CellMask, 'holes');
        CellMask = imerode(CellMask, se);
        CellMask = bwareafilt(CellMask, 1);
        
        baseMask = baseMask & ~CellMask;
        
        % size of cell in pixel
        numPix = sum(sum(CellMask)); 
        perimImg = bwperim(CellMask);

        % obtain centroid and distance information
        imgInfo = regionprops(CellMask, ImageSelected, 'Centroid', 'WeightedCentroid');
        centroidCell = imgInfo.Centroid;
        centroidImg = imgInfo.WeightedCentroid;
        cellCenter = false(size(CellMask));
        cellCenter(round(imgInfo.Centroid(2)), round(imgInfo.Centroid(1))) = true;
        distCenter = bwdist(cellCenter) .* double(CellMask);
        
        wcellCenter = false(size(CellMask));
        wcellCenter(round(imgInfo.WeightedCentroid(2)), round(imgInfo.WeightedCentroid(1))) = true;
        wdistCenter = bwdist(wcellCenter) .* double(CellMask);
        
        
        % calculate root mean square distance = RG = standard deviation?!
        RG = sqrt(sum(sum(distCenter.^2)) / numPix);
        % calculate distance of signals from center
        distCent = sqrt(sum((centroidCell - centroidImg).^2));
        % calculate PI
        PI = distCent / RG;
        
        % calculate variance of pixels within cell
        mu2p = sum(sum(wdistCenter.^2)) / numPix;
        % calculate variance of molecules: intensity based (see page 181 of reference)
        mu2 = sum(sum( (wdistCenter.^2) .* double(ImageSelected) / double(sum(sum(ImageSelected))) ));
        % calculate DI
        DI = mu2 / mu2p;
         
        % concat info
        PIDI = {PI, DI, FileName, cntCell};
        
        % save results - create folder
        if chkDir == 1
            mkdir(folderName, 'isolated');
            chkDir = 0;
        end
        cd(strcat(folderName, '\', 'isolated'));
        
        % save results - write image
        rgbImg = uint8(zeros(size(cellImage, 1), size(cellImage, 2), 3));
        tempRGB1 = imadjust(cellImage);
        tempRGB1(perimImg) = 255;
        tempRGB2 = imadjust(cellImage);
        tempRGB2(perimImg) = 0;
        rgbImg(:, :, 1) = tempRGB1;
        rgbImg(:, :, 2) = tempRGB2;
        rgbImg(:, :, 3) = tempRGB1;
        
        imwrite(rgbImg, strcat('AnalyzedCell', num2str(cntFiles), '_', num2str(cntCell), '.tif'));
        
        % save results - write excel
        cd(folderName);
        xlsName = strcat(folderName, '.xlsx');
        xlswrite(xlsName, PIDI, 'Sheet1', strcat('A', num2str( 2 + xlsLine )));
        xlsLine = xlsLine + 1;
        
        % continue to next cell
        cntCell = cntCell + 1;
        baseMask(CellMask) = 0;
        choice = questdlg('analyze the same image?', 'Reanalysis', 'Yes', 'No', 'Yes');
        
        switch choice
            case 'Yes'
                multipleCell = true;
            case 'No'
                multipleCell = false;
        end

    end

end

% label excel datasheet
dataLabel = {'PI', 'DI', 'FileName', 'cell'};
xlswrite(xlsName, dataLabel, 'Sheet1', strcat('A', num2str( 1 )));

disp(strcat('finished analyss, check... ', xlsName, ' in your folder'));
close all;
 


end