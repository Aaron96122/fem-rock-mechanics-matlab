function [cartd,djacb,gpcod] =my_jacobian(ielem,elcod,kgasp,shape,deriv,nnode,ndime)
format long;
%% gauss point coordinates =========================================
cg              = shape * elcod';  % use shape functions to calculate the global coordinates at the Gaussian integration points
gpcod(1, kgasp) = cg(1);           % matrix expansion
gpcod(2, kgasp) = cg(2);

%% jacobian matrix =================================================
xjacm = deriv * elcod';   % [2 * 4] * [4 * 2]

%% Determinate of Jacobian ========================================
djacb = xjacm(1,1) * xjacm(2,2) - xjacm(1,2) * xjacm(2,1);
if(djacb <= 0.0)
    fprintf('Element No: %5d\n',ielem);
    error('Program terminated zero or negative area');
end

%% cartesion derivatives ==========================================
xjaci(1,1) = xjacm(2,2)/djacb;
xjaci(2,2) = xjacm(1,1)/djacb;
xjaci(1,2) = -xjacm(1,2)/djacb;
xjaci(2,1) = -xjacm(2,1)/djacb;

cartd = xjaci * deriv; % derivatives of shape functions with respect to global Cartesian coordinates

end %endfunction
