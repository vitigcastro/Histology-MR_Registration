function hObjFound = findObjInObjChildren(hObject, tag)
% This function returns the handles of the object which have a specific tag
% among the children of a given object handles
%
% INPUT
%   - hObject: handles to the object whose children we want to search.
%   - tag: string with the name of the object
%
% OUTPUT
%   Handles to the searched object, or 0 if it is not found
%

    hObjFound = 0;
    
    hObject_chil = get(hObject, 'Children');
    
    i=1;
    found = false;
    while(~found)
        if(strcmp(hObject_chil(i).Tag, tag))
            hObjFound = hObject_chil(i);
            found = true;
        end
        i = i + 1;
    end

end