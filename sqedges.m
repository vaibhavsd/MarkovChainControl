function [otpt ] = sqedges( shape,k )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% shape= [3 3]; k=1;

x= shape(1); y= shape(2);

allelem= 1:1:x*y;
mat= vec2mat(allelem,x);

mat1= mat(:,1:x-1);
mat2= mat(:,2:x);
mat3= mat(1:y-1,:);

mat4= mat(2:y,:);

grid= [mat1(:),mat2(:);mat3(:),mat4(:)];

switch k
    case 1
otpt= [grid(:,1)',grid(:,2)']; %s

    case 2
otpt= [grid(:,2)',grid(:,1)']; %t

end


end

