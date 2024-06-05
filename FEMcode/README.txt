===================================================

    File directory
===================================================


1). pstress.m 
    
    Main file to solve a 2D plane stress FEM problem
    for a rectangular slab.

2). elementF.m
    
    Compute the element force vector

3). elementK.m

    Compute the element stiffness matrix

4). globalF.m/ globalK.m

    Stitch together the global stiffness matrix
    and force vector.

5). list_dofs.m/list_nodes.m

    Enumerate the indices for a given element.