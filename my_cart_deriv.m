function [dgdx, dvolum] = my_cart_deriv(npoin,nelem,nnode,nstre,ndime,ndofn, ngaus,ntype,lnods,coord,posgp,weigp)
ngaus2 = ngaus;
if(nnode == 3)          % Gaussian integration point at every direction is 1 for the three-node triangular element
    ngaus2 = 1;
end
mgaus  = ngaus * ngaus2; %
dvolum = zeros(nelem, mgaus);
dgdx   = zeros(nelem, mgaus, ndime, nnode);

elcod  = zeros(ndime, nnode); % element coordinates
for ielem = 1 : nelem
    for inode = 1 : nnode
        lnode = lnods(ielem,inode);
        for idime = 1 : ndime
            elcod(idime, inode) = coord(lnode, idime);
        end
    end
    
    % gauss points:
    kgasp = 0;
    for igaus = 1 : ngaus         % first loop along the x axis
        exisp = posgp(igaus);     % first local coordinate
        for jgaus=1 : ngaus2      % first loop along the y axis
            etasp = posgp(jgaus); % second local coordinate
            if(nnode ==3)         % for three-node triangular element
                etasp = posgp(ngaus + igaus);
            end
            
            kgasp = kgasp + 1;
            mgaus = mgaus + 1;
            
            % The shape function values and derivative values are obtained based on the element information and integration point coordinates.
            [shape,deriv] =my_shape_fun(exisp, etasp, nnode);
            
            % Jacobian matrix, Jacobian matrix determinant, and derivatives of shape functions of Cartesian coordinates
            [cartd,djacb,gpcod]=my_jacobian(ielem,elcod,kgasp,shape,deriv,nnode,ndime);
            
            % thickness ?
            dvolu = djacb * weigp(igaus) * weigp(jgaus);
            if(nnode == 3)
                dvolu = djacb * weigp(igaus); % one more integration derivative
            end
            dvolum(ielem, kgasp) = dvolu;
            
            % matrix merge for storing
            for idime = 1 : ndime
                for inode = 1 : nnode
                    dgdx(ielem, kgasp, idime, inode) = cartd(idime, inode);    % derivatives[2 * 4]
                end
            end
            
        end %igaus
    end %jgaus
end % ielem

end

%endfunction




