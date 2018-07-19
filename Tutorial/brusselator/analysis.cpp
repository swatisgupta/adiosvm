/*
 * Analysis code for the Brusselator application.
 * Reads variable u_real, u_imag, v_real, v_imag, and computes the norm of U and V.
 * Writes the variables and their computed norm to an ADIOS2 file.
 *
 * Kshitij Mehta
 */
#include <iostream>
#include <stdexcept>
#include <cstdint>
#include "adios2.h"

double compute_norm (double *** array, int d1, int d2, int d3) {
    norm = 0.0;
    int i,j,k;

    for(i=0; i<d1; i++) {
        for(j=0; j<d2; j++) {
            for(k=0; k<d3; k++) {
                norm = 
            }
        }
    }
}

/*
 * Helper function to dynamically create a fully contiguous 3D array
 */
double *** create3DArray(int x, int y, int z) {
    double *** array = new double ** [x];
    array[0] = new double * [x*y];
    array[0][0] = new double [x*y*z];

    int i,j,k;
    for( i = 0; i < x; i++) {
        if (i < x -1 ) {
            array[0][(i+1)*y] = &(array[0][0][(i+1)*z*y]);
            array[i+1] = &(array[0][(i+1)*y]);
        }

        for( j = 0; j < y; j++) {     
            if (j > 0) array[i][j] = array[i][j-1] + z;
        }
    }

    //cout << endl;
    return array;
}

/*
 * Delete (de-allocate) a 3D complex array created using the function 'create3DArray'
 */
void delete3DArray (double *** array) {
    delete[] array[0][0];
    delete[] array[0];
    delete[] array;
}

/*
 * MAIN
 */
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank, comm_size;
    bool firstStep = true;

    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &comm_size);

    double *** u_real_data = NULL;
    double *** u_imag_data = NULL;
    double *** v_real_data = NULL;
    double *** v_imag_data = NULL;

    std::vector<long unsigned int> shape_u_real;
    std::vector<long unsigned int> shape_u_imag;
    
    // adios variable declarations
    adios2::Variable<double> var_u_real, var_u_imag, var_v_real, var_v_imag;

    adios2::ADIOS ad ("analysis_adios2.xml", MPI_COMM_WORLD, adios2::DebugOFF);
    adios2::IO ad_io = ad.DeclareIO("analysis");
    adios2::Engine ad_engine = ad_io.Open("./data/brusselator.bp", adios2::Mode::Read, MPI_COMM_WORLD);

    // read data per timestep
    while(true) {
        adios2::StepStatus status = ad_engine.BeginStep (adios2::StepMode::NextAvailable, 0.0f);
        if (status != adios2::StepStatus::OK)
            break;

        var_u_real = ad_io.InquireVariable<double>("u_real");
        var_u_imag = ad_io.InquireVariable<double>("u_imag");

        if (!var_u_real) {
            std::cout << "ERROR: u_real not found. Exiting.." << std::endl;
            break;
        }

        shape_u_real = var_u_real.Shape();
        shape_u_imag = var_u_imag.Shape();

        //std::cout << "u_real shape: " << std::endl;
        //for (int i=0; i<shape_u.size(); i++)
        //    std::cout << shape_u[i] << std::endl;

        // Set decomposition -- use 2decomp library?
        
        // Create data arrays
        if (NULL == u_real_data) {
            u_real_data = create3DArray (shape_u_real[0]/comm_size, shape_u_real[1], shape_u_real[2]);
            u_imag_data = create3DArray (shape_u_imag[0]/comm_size, shape_u_imag[1], shape_u_imag[2]);
        }
        var_u_real.SetSelection(adios2::Box<adios2::Dims>(
                    {shape_u_real[0]/comm_size*rank,0,0},{shape_u_real[0], shape_u_real[1], shape_u_real[2]}));
        var_u_imag.SetSelection(adios2::Box<adios2::Dims>(
                    {shape_u_imag[0]/comm_size*rank,0,0},{shape_u_imag[0], shape_u_imag[1], shape_u_imag[2]}));

        ad_engine.Get<double>(var_u_real, **u_real_data, adios2::Mode::Deferred);
        ad_engine.Get<double>(var_u_imag, **u_imag_data, adios2::Mode::Deferred);
        ad_engine.PerformGets();
        ad_engine.EndStep();

        std::cout << u_real_data[0][0][0] << std::endl;

        norm_u = compute_norm(u_real_data, shape_u_real[0]/comm_size, shape_u_real[1], shape_u_real[2]);
        norm_v = compute_norm(v_real_data, shape_v_real[0]/comm_size, shape_v_real[1], shape_v_real[2]);

        // write them norms out motherfucker
    }

    // cleanup
    ad_engine.Close();
    MPI_Finalize();
    return 0;
}

