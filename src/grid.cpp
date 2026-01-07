#include <iostream>
#include <petscsnes.h>
#include <petscdmda.h>


// following tutorials at 
// - https://petsc.org/release/manual/vec/
// - page 43: https://drive.google.com/drive/u/0/folders/18FI9Vwx9tLODL8zyXvJLeR_C8hAPYoPd 


PetscErrorCode formMatrix (DM da, Mat A) {
    DMDALocalInfo info;
    MatStencil row, col[5];
    PetscReal hx, hy, v[5];
    PetscInt i, j, ncols;
    DMDAGetLocalInfo (da, &info);
    hx = 1.0 / (info.mx-1); 
    hy = 1.0 / (info.my-1);

    for ( j=info.ys; j<info.ys+info.ym; j++ ) 
    {
        for ( i=info.xs; i<info.xs+info.xm; i++ ) 
        {
            // do stuff, then
            // MatSetValuesStencil(A, 1, &row, ncols, col, v, INSERT_VALUES);
        }
    }
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY) ;
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY) ;
    return 0;
}


int main(int argc, char** argv)
{
    DM da; // structured grid

    Vec u,  // displacement vector
        f,  // external wind pressure
        p;  // internal pressure

    SNES snes;

    char        file[256] = "path/to/file/containing/options";
    char        help[256] = "thing to print when used with -help";

    PetscFunctionBeginUser;
    // PetscCall(PetscInitialize(&argc, &argv, file, help) );
    PetscCall(PetscInitialize(&argc, &argv, NULL, help) );

    PetscCall(SNESCreate(PETSC_COMM_WORLD, &snes));

    DMDACreate2d(
        PETSC_COMM_WORLD,
        DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,  // boundary types (if periodic or not)
        DMDA_STENCIL_STAR, // type of stencil?
        9, 9,  // M, N the global dimensions of the grid, change default 9x9 size using −da_grid_x M −da_grid_y N
        PETSC_DECIDE, PETSC_DECIDE,  // nb of communicators for each axis
        1,  // nb of degrees of freedom per node
        1,  // stencil width: the number of grid points away from the center of the stencil
        NULL, NULL,  // arrays containing the number of nodes in each process’ portion of the grid, or NULL; if non-null, these
                    // must be of length m and n, respectively, and the sums of lx[] and ly[] must be M and N, respectively
        &da // output grid
    );

    DMCreateGlobalVector(da, &u);
    // VecDuplicate(u, &b);
    

    // SNESSetDM(SNES snes,DM dm);

}