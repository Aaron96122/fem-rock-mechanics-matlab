function [gstif]= my_stiff(npoin,nelem,nnode,nstre,ndime,ndofn,ngaus,ntype,lnods,coord,props,posgp,weigp,dgdx,dvolum,dtime,constk,cenerg,constl,constn,coneta,stres,stran,tdisp,tdisp_old)
format long;
%% order of integration 
ngaus2 = ngaus;
if(nnode == 3)
    ngaus2 = 1;
end
mgaus = ngaus * ngaus2;

%% initialize global and local stiffness 
ntotv  = npoin * ndofn;   % global df
nevab  = nnode * ndofn;   % df for nodes at FEM element
ndofn2 = ndofn - 1;       % df for displcement [2]
ndofn3 = ndofn - 2;       % df for field [1]
ntotv2 = npoin * ndofn2;  % index
nevab2 = nnode * ndofn2;  % df for displacement at each FEM element [4*2]
nevab3 = nnode * ndofn3;  % df for field at each FEM element [1*4]

%% global stiffness 
gstif = sparse(ntotv,ntotv);

%% element stiffnesses 
estif  = zeros(nelem ,nevab,   nevab); %K_all
estif1 = zeros(nelem, nevab2, nevab2); %K_dis_dis
estif2 = zeros(nelem, nevab2, nevab3); %K_dis_phi
estif3 = zeros(nelem, nevab3, nevab2); %K_phi_dis
estif4 = zeros(nelem, nevab3, nevab3); %K_phi_phi

eload  = zeros(nelem, nevab);   %F_all
eload1 = zeros(nelem, nevab2);  %F_dis
eload2 = zeros(nelem, nevab3);  %F_phi

%% element nodal values 
ephi  = zeros(nelem, nnode);  % phase field value
ephir = zeros(nelem, nnode);  % phase field differential value  (d + dt) - (t)

for inode =1 : nnode
    lnode = lnods(:, inode);
    itotv = ntotv2 + lnode;
    ephi(:, inode)  = tdisp(itotv);                     % marix based on the element list
    ephir(:, inode) = tdisp(itotv) - tdisp_old(itotv);  % same
end

%% integrate element stiffness
kgasp = 0;
for igaus = 1 : ngaus          % first gaussian integral loop 
    exisp = posgp(igaus);         % x gauss point coordinate
    for jgaus =  1 : ngaus2    % second gaussian integral loop
        etasp  = posgp(jgaus);    % y gauss point coordinate
        if(nnode == 3)
            etasp = posgp(ngaus + igaus);
        end
        
        kgasp = kgasp+1;  % point counter  
        [shape, deriv] = my_shape_fun(exisp, etasp, nnode);
        
        phigp  = zeros(nelem, 1);
        phirgp = zeros(nelem, 1);
        for inode = 1 : nnode
            phigp  = phigp  +  ephi(:, inode)  * shape(inode);
            phirgp = phirgp + ephir(:, inode)  * shape(inode);
        end
        
        %matrix computation of D and C
        mtype   = 1;
        [dmatx] = my_D_Mat(mtype, ntype, nstre, props);
        
        %matrix computation of B
        [bmatx1] = my_bmat_v1(dgdx, nelem, nnode, nstre, nevab2, kgasp);
        [bmatx2] = my_bmat_v2(dgdx, nelem, nnode, nstre, nevab2, kgasp);
        [dbmat]  = my_dbe2(nelem,nevab2,nstre,bmatx2,dmatx);
        
        %-matrix construction of element stiffness
        % element stiffness ==== estif1:
        for ievab = 1 : nevab2
            for jevab = 1 : nevab2
                for istre = 1 : nstre
                    estif1( :,ievab,jevab) = estif1( :,ievab,jevab) + ...
                        ((1.0 - phigp).^2 + constk) .* bmatx2( :,istre,ievab) .* dbmat( :,istre,jevab) .* dvolum(:,kgasp);
                end
            end
        end
        
        % element stiffness ==== estif2:
        dummy = zeros(nelem,nevab2);
        for istre =1:nstre
            for inode =1:nnode
                dummy( :,istre,inode) =stres( :,kgasp,istre) * shape(inode);
            end
        end
        for ievab = 1 : nevab2
            for jevab = 1 : nevab3
                for istre = 1 : nstre
                    estif2( :,ievab,jevab) = estif2( :,ievab,jevab) ...
                        -2.0 * (1.0-phigp) .* bmatx2( :,istre,ievab) .* dummy( :,istre,jevab) .* dvolum( :,kgasp);
                end
            end
        end
        
        % element stiffness ==== estif4 --- part1
        for ievab = 1 : nevab3
            for jevab = 1 : nevab3
                for istre = 1 : ndime
                    estif4( :, ievab, jevab) = estif4( :, ievab, jevab) + ...
                        cenerg * constl * bmatx1( :,istre,ievab) .* bmatx1( :,istre,jevab) .* dvolum( :, kgasp);
                end
            end
        end
        
        % element stiffness ==== estif4 --- part2 strain energ
        senerg = zeros(nelem,1);
        for istre = 1 : nstre
            senerg = senerg + 0.5 * stres( :, kgasp, istre) .* stran( :, kgasp, istre); %弹性应变能累加
        end
        
        for inode = 1 : nnode
            for jnode = 1 : nnode
                estif4( :, inode, jnode) = estif4( :, inode, jnode) + ...
                    ((cenerg/constl) + 2.0 * senerg) .* shape(inode) * shape(jnode) .* dvolum( :,kgasp);
            end
        end
        
        % element stiffness ==== estif4 --- part3 penalty term:
        constx = zeros(nelem,1);
        inrange = (phirgp < 0 );
        constx(inrange) = -phirgp(inrange);
        for inode = 1 : nnode
            for jnode = 1 : nnode
                estif4( :, inode, jnode) = estif4( :, inode, jnode) + ...
                    (coneta/dtime) * constx.^(constn-1) * shape(inode) * shape(jnode) .* dvolum( :, kgasp);
            end
        end
    end %jgaus
end %igaus

%% assemble global stiffness matrix 
%% --- assemble estif1:
for inode = 1 : nnode               %% one loop
    lnode = lnods( :, inode);
    
    for idofn = 1 : ndofn2       %index: (2-1)*2+1=3, (2-1)*2+2=4
        ievab = (inode - 1) * ndofn2 + idofn;           
        itotv = (lnode - 1) * ndofn2 + idofn;          
        %
        for jnode = 1 : nnode       %% one loop
            knode = lnods( :, jnode);
            
            for jdofn =1 : ndofn2 %same
                jevab = (jnode-1) * ndofn2 + jdofn;     
                jtotv = (knode-1) * ndofn2 + jdofn;     
                
                gstif = gstif +sparse(itotv,jtotv,estif1( :,ievab,jevab),ntotv,ntotv);
            end
        end
    end
end

%% assemble estif2 :
for inode = 1 : nnode
    lnode = lnods( :,inode);
    for idofn = 1 : ndofn2
        ievab = (inode - 1) * ndofn2 + idofn;
        itotv = (lnode - 1) * ndofn2 + idofn;
        %
        for jnode = 1 : nnode
            knode = lnods( :, jnode);
            for jdofn = 1 : ndofn3
                jevab = (jnode - 1) * ndofn3 + jdofn;
                jtotv = (knode - 1) * ndofn3 + jdofn + ntotv2;  
                
                gstif = gstif +sparse(itotv,jtotv,estif2( :,ievab,jevab),ntotv,ntotv);
                
            end
        end
    end
end

%% assembe estif 3 as transpose of estif2:
for inode = 1 : nnode
    lnode = lnods( :, inode);
    for idofn = 1 : ndofn3
        ievab = (inode - 1) * ndofn3 + idofn;
        itotv = (lnode - 1) * ndofn3 + idofn + ntotv2;  

        for jnode = 1 : nnode
            knode = lnods( :, jnode);
            for jdofn = 1 : ndofn2
                jevab = (jnode - 1) * ndofn2+jdofn;
                jtotv = (knode - 1) * ndofn2+jdofn;
                gstif = gstif + sparse(itotv, jtotv, estif2( :,jevab,ievab), ntotv, ntotv); 
            end
        end
    end
end

%% assemble estif4 =================================================
for inode = 1 : nnode
    lnode = lnods( :, inode);
    for idofn = 1 : ndofn3
        ievab = (inode - 1) * ndofn3 + idofn;
        itotv = (lnode - 1) * ndofn3 + idofn + ntotv2;
        %
        for jnode = 1 : nnode
            knode = lnods( :, jnode);
            for jdofn = 1 : ndofn3
                jevab = (jnode - 1) * ndofn3 + jdofn;
                jtotv = (knode - 1) * ndofn3 + jdofn + ntotv2;
                gstif = gstif + sparse(itotv, jtotv, estif4( :,ievab,jevab), ntotv, ntotv);
            end
        end
    end
end

end %endfunction





