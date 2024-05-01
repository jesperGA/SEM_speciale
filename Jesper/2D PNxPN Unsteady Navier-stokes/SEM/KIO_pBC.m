function [omega] = KIO_pBC(Dx,Dy,ux,uy,free)
%Function to calculate double curl globally, and then return boundary terms
%for the second split calculation. 
temp = [Dx*(Dy*uy)-Dy*(Dy*ux);
          -Dx*(Dx*uy)+Dy*(Dx*ux)];
omega = zeros(size(temp));
omega(~free) = temp(~free);




end