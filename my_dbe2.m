function [dbmat] = my_dbe2(nelem,nevab,nstre,bmatx,dmatx)
%multiply bmatx with dmatx

format long;
dbmat = zeros(nelem,nstre,nevab);

for istre=1 : nstre                  % row>3
    for ievab = 1 : nevab            % colume>8
        %dbmat( :,istre,ievab)=0.0;
        for jstre = 1 : nstre        % matrix multiplication
            dbmat( :, istre, ievab) = dbmat( :, istre, ievab) + dmatx(istre, jstre) * bmatx( :,jstre,ievab);
        end
    end
end
end %endfunction
