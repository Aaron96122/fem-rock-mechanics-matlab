function [bmatx1]=my_bmat_v1(dgdx,nelem,nnode,nstre,nevab,kgasp)
format long;
bmatx1 = zeros(nelem, nstre, nevab);
for inode = 1 : nnode
    bmatx1(:, 1, inode) = dgdx(:, kgasp, 1, inode);
    bmatx1(:, 2, inode) = dgdx(:, kgasp, 2, inode);
end

end %endfunction

