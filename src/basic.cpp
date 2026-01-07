#include <iostream>
#include <petscsnes.h>
#include <petscdmplex.h>


/* DONNEES DANS EXEMPLE GENEPI 

Membrane:
    primary E for membrane material = 4000
    Poisson = 0.3
    shear modulus G12 = 700
    surface density = 0.05

Damping:
    Axial edge damping coefficient sigma = 0
    Point-wise damping coefficient = -2.0e-03
    Relative velocity damping coefficient = 0

*/


/* https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
    Lamé parameter lambda = E*nu / ( (1+nu)*(1-2*nu) )
*/


/*
    Membrane mechanics with wrinkling behavior based on principal stress criterion:
    - TAUT: sigma_2 > 0 (both principal stresses positive)
    - WRINKLED: sigma_1 > 0 && sigma_2 <= 0 (only one principal stress positive)
    - SLACK: sigma_1 <= 0 (both principal stresses negative)
*/



// // https://petsc.gitlab.io/annual-meetings/2023/tutorials/petsc_annual_meeting_2023_tutorial.pdf
// Different flags:
//      -snes_ksp_ew            : set adaptive forcing in Newton-Krylov
//      -ksp_type / -pc_type    : set Krylov method, for example lu,mg for pc_type and bcgs for ksp
//      -snes_monitor
//      -snes_converged_reason  : why did non-linear solver halt
//      -lidvelocity 100
//      -da_grid_x 16
//      -da_grid_y 16
//      -ksp_converged_reason   : why did linear solver halt
//      -log_view:log.txt       : log PETSc performance in log.txt
//      -ksp_rtol               : linear solve tolerance


class PhysicalProperties
{
private:
    PetscScalar E_Young,
                nu_Poisson,
                G12_shear,
                surface_density;
public:
};

class FieldNode
{
private:
    // needs to be adapted
    PetscScalar u,
                sigma,
                eps;
public:
};

class Context
{
private:
    PhysicalProperties  properties;
    Vec                 x_0[];
public:
};


/*
    Function to take in some custom command line options
*/
static PetscErrorCode ProcessOptions(MPI_Comm comm, Context *context)
{
    PetscInt sol;

    PetscFunctionBeginUser;
    //   context->sol = SOL_QUADRATIC;  
    PetscOptionsBegin(comm, "", "Problem Options", "DMPLEX");
    //   sol = context->sol;
    //   PetscCall(PetscOptionsEList("-sol", "The MMS solution", "basic.cpp", SolTypes, PETSC_STATIC_ARRAY_LENGTH(SolTypes) - 3, SolTypes[context->sol], &sol, NULL));
    //   context->sol = (SolType)sol;
    PetscOptionsEnd();
    PetscFunctionReturn(PETSC_SUCCESS);
}


/*
    Basic function to create a mesh
*/
static PetscErrorCode CreateMesh(
    MPI_Comm comm, Context* context, 
    DM *dm)
{
    PetscFunctionBeginUser;
    PetscCall(DMCreate(comm, dm));
    PetscCall(DMSetType(*dm, DMPLEX)); // THERE IS A CHOICE HERE 
    PetscCall(DMSetFromOptions(*dm));
    PetscCall(DMViewFromOptions(*dm, NULL, "-dm_view"));
    PetscFunctionReturn(PETSC_SUCCESS);
}


// void f_function(
//     PetscInt dim, PetscInt Nf, PetscInt NfAux,                                  // the coordinate dimension, the number of fields and the number of auxilliary fields
//     const PetscInt uOff[], const PetscInt uOff_x[],                             // the offset into u[] and u_t[] for each field and the offset into u_x[] for each field
//     const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],    // each field, time derivative and gradient evaluated at the current point
//     const PetscInt aOff[], const PetscInt aOff_x[],                             // the offset into a[] and a_t[] for each field and the offset into a_x[] for each field
//     const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],    // each auxiliary field evaluated at the current point ie not the main one we're solving for (p in Stokes) ?
//     PetscReal t, const PetscReal x[],                                           // time t and 3D position x
//     PetscInt numConstants, const PetscScalar constants[],                       // constants to be used in the function ie nu, E, ...
//     PetscScalar f0[])                                                           // output of the function
// {
//     for (int i=0; i<dim; i++) f0[i] = 0;
// }

// void f_jacobian(
//     PetscInt dim, PetscInt Nf, PetscInt NfAux,                                  // the coordinate dimension, the number of fields and the number of auxilliary fields
//     const PetscInt uOff[], const PetscInt uOff_x[],                             // the offset into u[] and u_t[] for each field and the offset into u_x[] for each field
//     const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],    // each field, time derivative and gradient evaluated at the current point
//     const PetscInt aOff[], const PetscInt aOff_x[],                             // the offset into a[] and a_t[] for each field and the offset into a_x[] for each field
//     const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],    // each auxiliary field evaluated at the current point ie not the main one we're solving for (p in Stokes) ?
//     PetscReal t, const PetscReal x[],                                           // time t and 3D position x
//     PetscInt numConstants, const PetscScalar constants[],                       // constants to be used in the function ie nu, E, ...
//     PetscScalar j0[])                                                           // output of the function
// {
//     for (int i=0; i<dim; i++) j0[i] = 0;
// }


PetscErrorCode local_function(
    const Vec x_1, const Vec x_2, const Vec x_3,
    Vec f,
    void *ctx)
{
    const PetscScalar*  uu_1;
    const PetscScalar*  uu_2;
    const PetscScalar*  uu_3;

    PetscScalar*        ff;

    PetscFunctionBeginUser;
    /*
    Get pointers to vector data.
        - For default PETSc vectors, VecGetArray() returns a pointer to
            the data array.  Otherwise, the routine is implementation dependent.
        - You MUST call VecRestoreArray() when you no longer need access to
            the array.
    */
    PetscCall(VecGetArrayRead(x_1, &uu_1));
    PetscCall(VecGetArrayRead(x_2, &uu_2));
    PetscCall(VecGetArrayRead(x_3, &uu_3));

    PetscCall(VecGetArray(f, &ff));





    PetscCall(VecRestoreArray(f, &ff));

    PetscFunctionReturn(PETSC_SUCCESS);
}



// PetscErrorCode thing()
// {
//     PetscErrorCode ierr;
//     const PetscReal hx = 1.0 / ( info->mx-1), 
//                     hy = 1.0 / ( info->my-1);
//     const Quad1D q = gausslegendre[ user->quadpts-1 ];
//     PetscReal x, y, lobj = 0.0;
//     PetscInt i, j, r, s;
//     MPI_Comm com;

//     // loop over all elements
//     for (j = info->ys; j < info->ys + info->ym; j++) 
//     {
//         if (j == 0) continue;
//         y = j*hy;
        
//         for (i = info->xs; i < info->xs + info->xm; i++) 
//         {
//             if (i == 0) continue;
//             x = i*hx;
//             const PetscReal ff[4] = { 
//                 user->f( x, y, user->p, user->eps ),
//                 user->f( x-hx, y, user->p, user->eps ),
//                 user->f( x-hx, y-hy, user->p, user->eps ),
//                 user->f( x, y-hy, user->p, user->eps ) 
//             } ;
            
//             const PetscReal uu[4] = {
//                 au[j] [i], au[j] [i-1],
//                 au[j-1][i-1], au[j-1][i] 
//             };

//             // loop over quadrature points on this element
//             for (r = 0; r < q.n; r++) for (s = 0; s < q.n; s++) 
//             {
//                 lobj += q.w[r] * q.w[s] * ObjIntegrandRef( info, ff, uu,q.xi[r], q.xi[s], user );
//             }
//         }
//     }

//     lobj *= hx * hy / 4.0; // from change of variables formula

//     PetscObjectGetComm ( (PetscObject) (info->da) ,&com);
//     MPI_Allreduce(&lobj, obj, 1, MPIU_REAL, MPIU_SUM, com);
//     PetscLogFlops(129*info->xm*info->ym);

//     PetscFunctionReturn(PETSC_SUCCESS);
// }



/*
    The F in F(u)=0 that the SNES is solving where
    Inputs:
        - snes: the SNES context
        - u: the input vector to F
        - ctx: optional user context
    Outputs:
        - F: the output value of the function for this u
*/
PetscErrorCode global_function(
    SNES snes, 
    Vec u, Vec f,
    void *ctx)
{
    
    PetscFunctionReturn(PETSC_SUCCESS);
}


/*
    The jacobian of F in F(u)=0 that the SNES is solving where
    Inputs:
        - snes: the SNES context
        - u: the input vector to F
        - ctx: optional user context
    Outputs:
        - J: Jacobian matrix output
        - B: optionally different matrix used to construct the preconditioner
*/
PetscErrorCode global_jacobian(
    SNES snes, 
    Vec u, Mat J, Mat B,
    void *ctx)
{
    PetscFunctionReturn(PETSC_SUCCESS);
}



int main(int argc, char** argv)
{
    Vec         x,      // approx solution vector
                r;      // residual vector
    Mat         J;      // Jacobian matrix

    SNES        snes;   // nonlinear solver context
    KSP         ksp;    // linear solver context
    PC          pc;     // preconditioner context

    PetscScalar *xx;    // util to set the initial guess x

    Context     contxt; // context that contains some info on physical properties etc
    //  -> should be made to contain user defined inputs!

    char        file[256] = "path/to/file/containing/options";
    char        help[256] = "thing to print when used with -help";

    PetscFunctionBeginUser;
    // PetscCall(PetscInitialize(&argc, &argv, file, help) );
    PetscCall(PetscInitialize(&argc, &argv, NULL, help) );


    /*
        Get some parameter from command line
    */
    // PetscCall(PetscOptionsGetInt(NULL, NULL, "-n", &n, NULL));

    

    /*
        Create abstract non-linear solver and its type and name(?)
    */
    PetscCall(SNESCreate(PETSC_COMM_WORLD, &snes));
    PetscCall(SNESSetType(snes, SNESNRICHARDSON));         // type of solver
    PetscCall(SNESSetOptionsPrefix(snes, "mysolver_"));


    /*
        Create an abstract vector
    */
    PetscCall(VecCreate(PETSC_COMM_SELF, &x));    // create a vector from scratch
    PetscCall(PetscObjectSetName((PetscObject)x,  // set its name
                                 "Approximate solution"));


    /*
        Set its local size (arg 2) and global size (arg 3)
        I'm letting PETSc decide the local size from
            context and size
    */
    PetscCall(VecSetSizes(x, PETSC_DECIDE, 3));
    /*
        Finalise creation
    */
    PetscCall(VecSetFromOptions(x));

    /*
        Duplicate a vector to another one to not
            create from scratch again
    */
    PetscCall(VecDuplicate(x, &r));


    /*
        Matrix is the same with an additional call of
            SetUp at the end
    */
    PetscCall(MatCreate(PETSC_COMM_SELF, &J));
    PetscCall(MatSetSizes(J,
                PETSC_DECIDE, PETSC_DECIDE,
                3, 3));
    MatSetType(J, MATAIJ);      // Set the type of matrix you want ie here compressed sparse row storage format
    PetscCall(MatSetFromOptions(J));
    PetscCall(MatSetUp(J));


    /*
        Set function evaluation routine and vector, with the user defined context
    */
    PetscCall(SNESSetFunction(snes, r, global_function, &contxt));

    /*
        Set Jacobian matrix data structure and Jacobian evaluation routine
    */
    PetscCall(SNESSetJacobian(snes, J, J, global_jacobian, &contxt));


    /*
        Set linear solver defaults for this problem. By extracting the
        KSP and PC contexts from the SNES context, we can then
        directly call any KSP and PC routines to set various options.
    */
    PetscCall(SNESGetKSP(snes, &ksp));
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCSetType(pc, PCNONE));       // set the type of pc (like LU, ...)
    PetscCall(KSPSetTolerances(ksp, 1.e-4, PETSC_CURRENT, PETSC_CURRENT, 20));  // set the tolerance we want for convergence


    /*
        Set SNES/KSP/KSP/PC runtime options, e.g.,
            -snes_view -snes_monitor -ksp_type <ksp> -pc_type <pc>
        These options will override those specified above as long as
        SNESSetFromOptions() is called _after_ any other customization
        routines.
    */
    PetscCall(SNESSetFromOptions(snes));


    /*
        Set the initial guess
    */
    PetscCall(VecGetArray(x, &xx));
    xx[0] = 2.0;
    xx[1] = 3.0;
    xx[2] = 4.0;
    PetscCall(VecRestoreArray(x, &xx));


    /*
        Solve the non-linear system
    */
    PetscCall(SNESSolve(snes, NULL, x));



    /*
        View the results
    */
    Vec f;
    PetscCall(VecView(x, PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(SNESGetFunction(snes, &f, 0, 0));
    PetscCall(VecView(r, PETSC_VIEWER_STDOUT_WORLD));



    /*
        Free everything we properly created (with VecCreate, ...)
    */
    PetscCall(VecDestroy(&x));
    PetscCall(VecDestroy(&r));
    PetscCall(MatDestroy(&J));
    PetscCall(SNESDestroy(&snes));


    PetscCall(PetscFinalize());


    return 0;
}








