function [interp_data] = oneD_element_interpolator(mesh,solution)
% 
IX = mesh.IX; %Element Topology matrix
X = mesh.X; %Position data matrix
N = size(IX(:,2:end-1),2)-1;
interp_data = [];
for e = 1:size(IX,1)
     
    nen = IX(e,2:end-1);  
    edof = nen; %Only true for 1D.
    x_data = X(nen,2).';
    x_fit = linspace(x_data(1),x_data(end),1000);

    y_data = solution(edof).';
    y_fit = LagrangeFormInterpolation(x_data,y_data,x_fit);
    interp_data = [interp_data; x_fit.',y_fit.'];
end

end

function fit = LagrangeFormInterpolation(knots,ydata,t)
%LagrangeFormInterpolation: CalculateS the values of the interpolating 
%                           polynomial in Lagrange form
%
%   fit = LagrangeFormInterpolation(knots,ydata,t)
%
% Input:
%   knots=[x0 x1 ... xn]   is a row of n+1 knot-values
%   ydata=[y0 y1 ... yn]   is a row of the corresponding n+1 y-values
%   t=[t1 ... tm]          is a row of all the m values that the inter-
%                          polating polynomial is to be evaluated in
% Output:
%   fit=[P(t1) ... P(tm)]  a row with the m function values of the
%                          interpolating polynomial

m=length(knots);
if (m~=length(ydata) ) 
    disp('Der skal være lige mange knuder og ydata'); 
    fit=NaN; return; 
end

cardinals=zeros(m,length(t));
for i=1:m
    cardinals(i,:)=cardinalpoly(knots,i,t)';
end
Y=repmat(ydata',1,length(t));
fit=sum(cardinals.*Y);

% --------- Subfunction ----------------------------------
end
function l = cardinalpoly(knots,i,t)
%cardinalpoly: Calculates the function values of i'th cardinal 
%              polynomial in all t-values based on the knots 
%
% Input:
%    i      the index (between 1 and n+1) of the cardinal polynomial
%    knots: the knots [x_1 ... x_(n+1)]
%    t      a row [t1 ... tm] of all the values that the 
%           i'th cardinal polynomial is to be evaluated in  
% Output:
%    l      a row with the m function values of the
%           i'th cardinal polynomial

l=ones(size(t)); 
m=length(knots);
for j=[1:i-1 i+1:m]
    l=l.*(t-knots(j))./(knots(i)-knots(j));
end
end