function [mesh] = liddriven(xNV,yNV,iglobV)
xNV = xNV ;yNV = yNV ;
topdofs = find(yNV == 1);
lower_dofs = find(yNV == 0);

left_dofs = find(xNV == 0);
right_dofs = find(xNV == 1);


mesh.bound = [topdofs,ones(size(topdofs)),-ones(size(topdofs))
    lower_dofs,ones(size(lower_dofs)),zeros(size(lower_dofs))
    left_dofs,ones(size(left_dofs)),zeros(size(left_dofs))
    right_dofs,ones(size(right_dofs)),zeros(size(right_dofs))
    topdofs,ones(size(topdofs))*2,zeros(size(topdofs))
    lower_dofs,ones(size(lower_dofs))*2,zeros(size(lower_dofs))
    left_dofs,ones(size(left_dofs))*2,zeros(size(left_dofs))
    right_dofs,ones(size(right_dofs))*2,zeros(size(right_dofs))];
mesh.IXv = iglobV;
mesh.Xv = [(1:length(xNV))',xNV,yNV]; %Save grid to mesh-struct

end