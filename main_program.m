clc, clear, close all
time0 = clock();                      % initial time

%% ==== File Reading
% global in;
inp   = fopen('fract_1ca.inp','r');     % FEM Mesh    [fract_1hca.inp]
out1  = fopen('result_1.out','w');      % Result File for data output
out2  = fopen('force-disp.out','w');    % Force-disp curve

%% ==== time integration parameters
nstep    =    4000;    % total time step
nprnt    =    01;      % time interval for data output
dtime    =    1.0;     % time step for calculation
miter    =    5;       % maximum iteration steps for N-R solution
toler    =    5.0e-3;  % error tolerance

tfacto   =    0.0;
dfacto   =    0.0005;
isolve   =    2;

%% ==== Material specific parameters
constk    = 1.0e-6;    % Small value to avoid overflow for cracked elements
cenerg    = 0.001;     % Critical strain energy for fracture
constl    = 0.125;     % Interface control parameter
constn    = 2;         % Power term for the penalty parameter
coneta    = 20.0e5*2;  % Magnitude of the penalty term

%% ==== input data
[npoin, nelem, nvfix, ntype, nnode, ndofn, ndime, ngaus, nstre, ...
    nmats, nprop, lnods, matno, coord, props, nofix, iffix, fixed] = ...
    my_input_data(inp, out1);

ntotv    = npoin * ndofn;
ndofn2   = 2;
ntotv2   = npoin * ndofn2;

%% ===== Initialization for the initial fracture
icase    = 1;          % for plate with nothing
[tdisp, stres, stran] = my_initiallize(nelem, npoin, nnode, ngaus, ndofn, nstre, icase);

%% ==== Gauss Integration Information
[posgp, weigp] = my_gauss(ngaus, nnode); % Gaussian integration points determined by the distribution of integration points and the number of nodes within the element

%% ==== Calculate the derivative matrix of each Gaussian node with respect to the global Cartesian coordinates within each element
[dgdx, dvolum] = my_cart_deriv(npoin,nelem,nnode,nstre,ndime,ndofn,ngaus,ntype,lnods,coord,posgp,weigp);

%% ==== time integration
for istep = 1 : nstep
    tfacto     = tfacto + dfacto;

    tdisp_old  = tdisp;             % Update [displacement, phase field], using the previous step [displacement, phase field] as old data.

    gforce     = zeros(ntotv, 1);   % Global force vector (ntotv ¼ npoin  ndofn).

    gstif      = my_stiff(npoin,nelem,nnode,nstre,ndime,ndofn,ngaus,ntype,lnods,coord,props,posgp,weigp,dgdx,dvolum,dtime,constk,cenerg,constl,constn,coneta,stres,stran,tdisp,tdisp_old);

    % ==== newton iteration
    for iter = 1 : miter
        % ==== boundary conditions
        [gstif,gforce,treac] = my_BC(npoin,nvfix,nofix,iffix,fixed,ndofn,tfacto,gstif,gforce,tdisp);

        % ==== Results Solving
        asdis = gstif\gforce;

        % ==== Data Updating
        tdisp = tdisp + asdis;

        % ==== adjust small deviations
        dummy   = tdisp(ntotv2+1: ntotv);
        inrange        = (dummy > 0.999);
        dummy(inrange) = 1.0;
        inrange        = ( dummy < 0.0);
        dummy(inrange) = 0.0;
        tdisp(ntotv2+1: ntotv) = dummy;

        % ==== calculate stress & strain increments

        [stran,stres] = my_stress_fract(asdis,nelem,npoin,nnode,ngaus,nstre,props,ntype,ndofn,ndime,lnods,matno,coord,posgp,weigp,dgdx,dvolum,tdisp);

        % ==== check norm for convergence
        normF = norm(gforce, 2);
        if(normF <= toler); break;  end

        % ==== calculate residual force vector
        [gforce]= my_residual(npoin,nelem,nnode,nstre,ndime,ndofn,ngaus,ntype,lnods,coord,props,posgp, ...
            weigp,dgdx,dvolum,dtime,constk,cenerg,constl,constn,coneta,stres,stran,tdisp,tdisp_old);

    end % end of Newton

    %% ==== print data for force-disp curves 
    lnode  = nofix(nvfix);
    nvfix2 = nvfix/2;       % extract the crack tip nodes.
    sumr   = 0.0;
    for ivfix = 1 : nvfix2
        sumr = sumr + treac(ivfix, 2);
    end
    fprintf(out2,'%14.6e %14.6e\n', tdisp((lnode-1)*2+2), sumr);  % Output the displacement in the y-direction at a specific point.

    %% ==== print results
    if(mod(istep,nprnt) == 0 )
        fprintf('Done step: %5d\n',istep);

        %--- write to vtk file with updated mesh
        for ipoin = 1 : npoin
            for idofn = 1 : ndofn2
                itotv               = (ipoin - 1) * ndofn2 + idofn;
                cord2(ipoin, idofn) = coord(ipoin,idofn) + 10.0 * tdisp(itotv);
                jtotv               = ntotv2 + ipoin;
                cont1(ipoin)        = tdisp(jtotv);
            end
        end
        my_vtk(npoin, nelem, nnode, lnods, cord2, istep, cont1);
    end %if

end % istep

%% ==== Computing Time
compute_time = etime(clock(), time0); %calculate the computation time

fprintf(out1,'compute time: %10.4f\n',compute_time);

fprintf('compute time: %10.4f\n',compute_time);

fclose('all');