C++ codes for calculating the infinitesimal degree of freedom in Miura-ori structure

Codes for the paper:
Chen, S., & Mahadevan, L. (2019). Rigidity percolation and geometric information
in floppy origami. Proceedings of the National Academy of Sciences, 116(17), 8119-8124.


Installation of the SuiteSparse is required. The codes ran with SuiteSparse 5.2.0.
origami::gen_DoF is the main function calculating the rank of the rigidity matrix (thus the DoF)

The codes for calculating the DoF in a main file can be:

    origami Origami(L_quad, L_quad, n_cst_coplanar,0,0);
    Origami.gen_random_cst_list();
    long dof = Origami.gen_DoF(Origami.constraint_all);
    
where L_quad is the "L" in the paper, n_cst_coplanar is the number of rigid quads. The constraint pattern is generated randomly in this case.
It is also possible to specify which quads are rigid by modifying the input vector of gen_DoF( )


Please contact siheng_chen@g.harvard.edu for any questions.
