function mesh = modify_to_roenquist_mesh(xN,yN,iglob)

    for e = 3:4
        glob_mod = iglob(:,:,e);
        for k = 1:size(glob_mod,1)
            if e==3
                x_spec = xN(glob_mod(k,end));
                yN(glob_mod(k,end)) = 1+1/4*sin(pi*x_spec);
                ratio(k) = (1/2+1/4*sin(pi*x_spec))/(1/2);
                yN(glob_mod(k,2:end-1)) = (yN(glob_mod(k,2:end-1))-1/2)*(ratio(k))+1/2;
            elseif e==4 && k==1
            else
                x_spec = xN(glob_mod(k,end));
                yN(glob_mod(k,end)) = 1+1/4*sin(pi*x_spec);
                % ratio = (1/2-1/4*x_spec)/(1/2);
                yN(glob_mod(k,2:end-1)) = (yN(glob_mod(k,2:end-1))-1/2)*(ratio(end-k+1))+1/2;
            end

        end

    end

    %Boundaries on west
    west_bound = iglob(1,:,[1,3]);
    west_bound_mat = [west_bound(:),ones(numel(west_bound),1),zeros(numel(west_bound),1)];

    south_bound = iglob(:,1,[1,2]);
    south_bound_mat = [south_bound(:),ones(numel(south_bound),1),sin(xN(south_bound(:))).*exp(-yN(south_bound(:)))];

    east_bound = iglob(end,:,[2,4]);
    east_bound_mat = [east_bound(:),ones(numel(east_bound),1),sin(xN(east_bound(:))).*exp(-yN(east_bound(:)))];

    north_bound = iglob(:,end,[3,4]);
    north_bound_mat = [north_bound(:),ones(numel(north_bound),1),sin(xN(north_bound(:))).*exp(-yN(north_bound(:)))];

    mesh.bound = [west_bound_mat;south_bound_mat;east_bound_mat;north_bound_mat];
    % figure()
    % plotmesh(iglob,xN,yN)
    mesh.IX = iglob;
    mesh.X = [(1:length(xN))',xN,yN]; %Save grid to mesh-struct

end