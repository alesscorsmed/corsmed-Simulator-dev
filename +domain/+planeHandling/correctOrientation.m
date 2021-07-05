function [imageCorrOrient,point1New,point2New,point3New,point4New,...
    topNew,bottomNew,rightNew,leftNew,fovXNew,fovYNew,...
    InPlanePhaseEncodingDirection] = ...
    correctOrientation(initialImage,point1,point2,point3,point4,...
    top,bottom,right,left,fovX,fovY,foldoverdir)

dimensions = strcat(bottom,top,right,left);
dimensions = sort(dimensions);

switch dimensions
    case 'AFHP'
        topNew = 'H';
        bottomNew = 'F';
        leftNew = 'A';
        rightNew = 'P';
        
        if(top=='A') 
            if(right=='H')
                point1New = point2;
                point2New = point4;
                point3New = point1;
                point4New = point3;
                imageCorrOrient = imrotate(initialImage,90);
            else
                % mirror + rotate counterclock
                point1New = point1;
                point2New = point3;
                point3New = point2;
                point4New = point4;
                initialImage = flip(initialImage,2);
                imageCorrOrient = imrotate(initialImage,90);
            end            
            fovXNew = fovY;
            fovYNew = fovX;
        elseif(top=='F')
                if(right=='A')
                    % rotate counterclock 2
                    point1New = point4;
                    point2New = point3;
                	point3New = point2;
                	point4New = point1;
                	imageCorrOrient = imrotate(initialImage,180);
                else
                    % mirror + rotate counterclock 2
                    point1New = point3;
                    point2New = point4;
                    point3New = point1;
                    point4New = point2;
                    initialImage = flip(initialImage,2);
                    imageCorrOrient = imrotate(initialImage,180);
                end
            fovYNew = fovY;
            fovXNew = fovX;
        elseif(top=='P')
                if(right=='F')
                    % rotate clock wise 1
                    point1New = point3;
                    point2New = point1;
                    point3New = point4;
                    point4New = point2;
                    imageCorrOrient = imrotate(initialImage,-90);
                else
                    % mirror + clock wise 1
                    point1New = point4;
                    point2New = point2;
                    point3New = point3;
                    point4New = point1;
                    initialImage = flip(initialImage,2);
                    imageCorrOrient = imrotate(initialImage,-90);
                end
            fovXNew = fovY;
            fovYNew = fovX;
        else
            if(right=='P')
                % Proper Orientation
                point1New = point1;
                point2New = point2;
                point3New = point3;
                point4New = point4;
                imageCorrOrient = initialImage;
            else
                % mirror
                point1New = point2;
                point2New = point1;
                point3New = point4;
                point4New = point3; 
                imageCorrOrient = flip(initialImage,2);
            end
            fovYNew = fovY;
            fovXNew = fovX;
        end
        
        
    case 'ALPR'
        topNew = 'A';
        bottomNew = 'P';
        leftNew = 'R';
        rightNew = 'L';       
               
        if(top=='A') 
            if(right=='L')
                % Proper Orientation
                point1New = point1;
                point2New = point2;
                point3New = point3;
                point4New = point4;
                imageCorrOrient = initialImage;
            else
                % mirror  
                point1New = point2;
                point2New = point1;
                point3New = point4;
                point4New = point3;
                imageCorrOrient = flip(initialImage,2);          
            end
            fovYNew = fovY;
            fovXNew = fovX;
        elseif(top=='L')
             if(right=='P')
                 % rotate clockwise 1
                 point1New = point3;
                 point2New = point1;
                 point3New = point4;
                 point4New = point2;
                 imageCorrOrient = imrotate(initialImage,-90);
             else
                 % mirror + clockwise 1
                 point1New = point4;
                 point2New = point2;
                 point3New = point3;
                 point4New = point1;
                 initialImage = flip(initialImage,2);
                 imageCorrOrient = imrotate(initialImage,-90);
             end
            fovXNew = fovY;
            fovYNew = fovX;
        elseif(top=='R')
            if(right=='A')
                % rotate counterclock 1
                point1New = point2;
                point2New = point4;
                point3New = point1;
                point4New = point3;
                imageCorrOrient = imrotate(initialImage,90);
            else
                % mirror + counterclock 1
                point1New  = point1;
                point2New  = point3;
                point3New  = point2;
                point4New  = point4;
                initialImage      = flip(initialImage,2);
                imageCorrOrient  = imrotate(initialImage,90);
            end
            fovXNew = fovY;
            fovYNew = fovX;
        else
            if(right=='R')
                % rotate clockwise 2
                point1New = point4;
                point2New = point3;
                point3New = point2;
                point4New = point1;
                imageCorrOrient = imrotate(initialImage,180);
            else
                % mirror + rotate 
                point1New  = point3;
                point2New  = point4;
                point3New  = point1;
                point4New  = point2;
                initialImage      = flip(initialImage,2);
                imageCorrOrient  = imrotate(initialImage,180);
            end
            fovYNew = fovY;
            fovXNew = fovX;
        end
        
        
    case 'FHLR'
        topNew = 'H';
        bottomNew = 'F';
        leftNew = 'R';
        rightNew = 'L';         
                      
        if(top=='F')
            if(right=='R')
                % rotate clockwise 2
                point1New = point4;
                point2New = point3;
                point3New = point2;
                point4New = point1;
                imageCorrOrient = imrotate(initialImage,180);
            else
                % mirror + rotate clockwise 2
                point1New = point3;
                point2New = point4;
                point3New = point1;
                point4New = point2;
                initialImage = flip(initialImage,2);
                imageCorrOrient = imrotate(initialImage,180);
            end
            fovYNew = fovY;
            fovXNew = fovX;
        elseif(top=='R')
            if(right=='F')
                % mirror + rotate counterclock
                point1New  = point3;
                point2New  = point4;
                point3New  = point1;
                point4New  = point2;
                initialImage      = flip(initialImage,2);
                imageCorrOrient  = imrotate(initialImage,90);
            else
                % rotate counterclock
                point1New = point2;
                point2New = point4;
                point3New = point1;
                point4New = point3;
                imageCorrOrient = imrotate(initialImage,90);
            end
            fovXNew = fovY;
            fovYNew = fovX;
        elseif(top=='L')
            if(right=='F')
                % rotate clockwise 1
                point1New = point3;
                point2New = point1;
                point3New = point4;
                point4New = point2;
                imageCorrOrient = imrotate(initialImage,-90);
            else
                % mirror + rotate clockwise 1
                point1New  = point4;
                point2New  = point2;
                point3New  = point3;
                point4New  = point1;
                initialImage      = flip(initialImage,2);
                imageCorrOrient  = imrotate(initialImage,-90);
            end
            fovXNew = fovY;
            fovYNew = fovX;
        else 
            if(right=='R')
                % mirror
                point1New = point2;
                point2New = point1;
                point3New = point4;
                point4New = point3; 
                imageCorrOrient = flip(initialImage,2);
            else
                % proper orientation
                point1New = point1;
                point2New = point2;
                point3New = point3;
                point4New = point4;
                imageCorrOrient = initialImage;
            end
            fovYNew = fovY;
            fovXNew = fovX;
        end
end

if strcmp(foldoverdir(1,1),topNew) || strcmp(foldoverdir(1,1),bottomNew)
    InPlanePhaseEncodingDirection = 'ROW';
else
    InPlanePhaseEncodingDirection = 'COL';
end