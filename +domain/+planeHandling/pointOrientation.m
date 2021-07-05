function pointOrient = pointOrientation(relPointToCenterOfPlane,textForPoint)

nonZeroElement = find(relPointToCenterOfPlane);
if nonZeroElement == 1
    if relPointToCenterOfPlane(nonZeroElement) > 0
        pointOrient = 'L';
    else
        pointOrient = 'R';
    end
    %disp(['The ',textForPoint,' is along the X axis - ',pointOrient])
elseif nonZeroElement == 2
    if relPointToCenterOfPlane(nonZeroElement) > 0
        pointOrient = 'P';
    else
        pointOrient = 'A';
    end
    %disp(['The ',textForPoint,' is along the Y axis - ',pointOrient])
elseif nonZeroElement == 3
    if relPointToCenterOfPlane(nonZeroElement) > 0
        pointOrient = 'H';
    else
        pointOrient = 'F';
    end
    %disp(['The ',textForPoint,' is along the Z axis - ',pointOrient])
end