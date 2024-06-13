%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHILIPPE SCHUCHERT            %
% SCI-STI-AK, EPFL              %
% philippe.schuchert@epfl.ch    %
% March 2021                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% cone from YALMIP is not serialized -> use this instead 
%

function rcone = rcone_serialized(x,y,z)
if ~isempty(z)
    % ||z||^2<2xy, x+y>0
    v = 1/sqrt(2);
    T = blkdiag([v,v;v,-v], eye(2*size(z,2)));
    
    Z = T*[x.';y.';sqrt(2)*real(z).';sqrt(2)*imag(z).'];
    
    rcone = cone(Z);
    % rcone = cone(Z(2:4),Z(1));
    % rcone = cone(Z(1:3));
    % rcone = [Z(1)<=Z(3),Z(2)<=Z(4),Z(1)>=0,Z(2)>=0];
    % rcone1 = rcone(z,x,y);
else
    rcone = [];
end
end