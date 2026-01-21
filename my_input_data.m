function [npoin,nelem,nvfix,ntype,nnode,ndofn,ndime,ngaus,nstre,nmats,nprop,lnods,matno,coord,props,nofix,iffix,fixed] = my_input_data(inp, out1)

%% read the Control data: ==========================================
npoin=fscanf(inp,'%d',1);  % number of total nodes  
nelem=fscanf(inp,'%d',1);  % element number
nvfix=fscanf(inp,'%d',1);  % Number of nodes having prescribed displacements.
ntype=fscanf(inp,'%d',1);  % Solution type, ntype = 1 for planestress, ntype = 2 for plane-strain.
nnode=fscanf(inp,'%d',1);  % Number of nodes per element.
ndofn=fscanf(inp,'%d',1);  % Number of DOF per node.
ndime=fscanf(inp,'%d',1);  % Number of Cartesian component dimension.
ngaus=fscanf(inp,'%d',1);  % Order of numerical integration
nstre=fscanf(inp,'%d',1);  % Number of stress components
nmats=fscanf(inp,'%d',1);  % Total number of different materials in the solution.
nprop=fscanf(inp,'%d',1);  % Number of material properties

%% Element node numbers & material property number=================
matno     = zeros(nelem, 1);         % Material types for the elements.
lnods     = zeros(nelem, nnode);     % Element nodal connectivity list.  
for ielem = 1 : nelem
    jelem = fscanf(inp, '%d', 1);
    dummy = fscanf(inp, '%d', [nnode+1,1]);
    for inode = 1 : nnode
        lnods(jelem, inode) = dummy(inode);
    end
    matno(jelem) = dummy(nnode + 1); 
end

%% Nodal coordinates:==============================================
coord     = zeros(npoin, ndime);    % Cartesian coordinates of each node
for ipoin = 1 : npoin
    jpoin = fscanf(inp, '%d', 1);
    dummy = fscanf(inp, '%lf %lf', [2,1]);
    for idime=1 : ndime
        coord(ipoin, idime)=dummy(idime);
    end
end

%% Constraint nodes and their values ==============================
nofix = zeros(nvfix, 1);          % Node numbers at which one or more DOFs are constrained
iffix = zeros(nvfix, 2);          % List of constrained DOFs
fixed = zeros(nvfix, 2);          % Prescribed value of any constrained DOFs.

for ivfix = 1 : nvfix
    nofix(ivfix) = fscanf(inp,'%d',1);
    dummy1 = fscanf(inp, '%d %d', [2,1]);
    dummy2 = fscanf(inp, '%lf %lf', [2,1]);
    for idime=1 : ndime
        iffix(ivfix, idime)=dummy1(idime);
        fixed(ivfix, idime)=dummy2(idime);
    end
end

%% Material properties ============================================
props = zeros(nmats, nprop);
for imats = 1 : nmats
    jmats = fscanf(inp, '%d',1);
    dummy = fscanf(inp, '%lf %lf',[2,1]);
    for iprop = 1 : nprop
        props(jmats, iprop) = dummy(iprop);
    end
end

%% printout =======================================================
fprintf(out1,'*****************************\n');
fprintf(out1,'*      FEM input data       *\n');
fprintf(out1,'*****************************\n');
fprintf(out1,'\n');

fprintf(out1,'Number of Elements           : %5d\n',nelem);
fprintf(out1,'Number of Node               : %5d\n',npoin);
fprintf(out1,'Number of Fixed nodes        : %5d\n',nvfix);
fprintf(out1,'Number of Nodes per element  : %5d\n',nnode);
fprintf(out1,'Number of Integration points : %5d\n',ngaus);
fprintf(out1,'Number of Materials          : %5d\n',nmats);
fprintf(out1,'Number of properties         : %5d\n',nprop);
fprintf(out1,'\n');

%% ================================================================
if(ntype == 1)
    fprintf(out1,'Plane-stress elasticity solution\n');
end

if(ntype == 2)
    fprintf(out1,'Plane-strain elasticity solution\n');
end

end %endfunction
