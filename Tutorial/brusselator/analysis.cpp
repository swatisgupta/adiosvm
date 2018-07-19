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
#include <cmath>
#include "adios2.h"

std::vector<double> compute_norm (std::vector<double> real_part, std::vector<double> imag_part, int dims) {
    int i;
    std::vector<double> norm;
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    norm.resize(dims);
    for (i=0; i<dims; i++)
        norm[i] = sqrt( real_part[i]*real_part[i] + imag_part[i]*imag_part[i] );

    return norm;
}

/*
 * MAIN
 */
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank, comm_size;
    
    int u_global_size, v_global_size;
    int u_local_size, v_local_size;
    
    bool firstStep = true;

    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &comm_size);

    std::vector<long unsigned int> shape_u_real;
    std::vector<long unsigned int> shape_u_imag;
    std::vector<long unsigned int> shape_v_real;
    std::vector<long unsigned int> shape_v_imag;
    
    std::vector<double> u_real_data;
    std::vector<double> u_imag_data;
    std::vector<double> v_real_data;
    std::vector<double> v_imag_data;

    std::vector<double> norm_u;
    std::vector<double> norm_v;
    
    // adios2 variable declarations
    adios2::Variable<double> var_u_real, var_u_imag, var_v_real, var_v_imag;

    // adios2 io object and engine init
    adios2::ADIOS ad ("analysis_adios2.xml", MPI_COMM_WORLD, adios2::DebugOFF);
    adios2::IO ad_io = ad.DeclareIO("analysis");
    adios2::Engine ad_engine = ad_io.Open("./data/brusselator.bp", adios2::Mode::Read, MPI_COMM_WORLD);

    // read data per timestep
    while(true) {

        // Begin step
        adios2::StepStatus status = ad_engine.BeginStep (adios2::StepMode::NextAvailable, 0.0f);
        if (status != adios2::StepStatus::OK)
            break;

        // Inquire variable
        var_u_real = ad_io.InquireVariable<double>("u_real");
        var_u_imag = ad_io.InquireVariable<double>("u_imag");
        var_v_real = ad_io.InquireVariable<double>("v_real");
        var_v_imag = ad_io.InquireVariable<double>("v_imag");

        if (!var_u_real) {
            std::cout << "ERROR: u_real not found. Exiting.." << std::endl;
            break;
        }

        shape_u_real = var_u_real.Shape();
        shape_u_imag = var_u_imag.Shape();
        shape_v_real = var_v_real.Shape();
        shape_v_imag = var_v_imag.Shape();

        // Calculate global and local sizes of U and V
        u_global_size = shape_u_real[0] * shape_u_real[1] * shape_u_real[2];
        u_local_size  = u_global_size/comm_size;
        v_global_size = shape_v_real[0] * shape_v_real[1] * shape_v_real[2];
        v_local_size  = v_global_size/comm_size;

        //std::cout << "u_real shape: " << std::endl;
        //for (int i=0; i<shape_u.size(); i++)
        //    std::cout << shape_u[i] << std::endl;
        
        // Set selection
        var_u_real.SetSelection(adios2::Box<adios2::Dims>(
                    {shape_u_real[0]/comm_size*rank,0,0},{shape_u_real[0], shape_u_real[1], shape_u_real[2]}));
        var_u_imag.SetSelection(adios2::Box<adios2::Dims>(
                    {shape_u_imag[0]/comm_size*rank,0,0},{shape_u_imag[0], shape_u_imag[1], shape_u_imag[2]}));

        // Allocate vectors to hold U and V data and their norms
        u_real_data.resize(u_local_size);
        u_imag_data.resize(u_local_size);
        v_real_data.resize(v_local_size);
        v_imag_data.resize(v_local_size);
        
        // Read adios2 data
        ad_engine.Get<double>(var_u_real, u_real_data.data(), adios2::Mode::Deferred);
        ad_engine.Get<double>(var_u_imag, u_imag_data.data(), adios2::Mode::Deferred);
        ad_engine.Get<double>(var_v_real, v_real_data.data(), adios2::Mode::Deferred);
        ad_engine.Get<double>(var_v_imag, v_imag_data.data(), adios2::Mode::Deferred);
        ad_engine.PerformGets();

        // End adios2 step
        ad_engine.EndStep();

        std::cout << u_real_data[0] << std::endl;

        // Compute norms
        norm_u = compute_norm(u_real_data, u_imag_data, u_local_size);
        norm_v = compute_norm(v_real_data, v_imag_data, v_local_size);

        std::cout << "Printing norm" << std::endl;
        std::cout << norm_u[1] << std::endl;
        std::cout << norm_v[1] << std::endl;

        // write U, V, and their norms out
    }

    // cleanup
    ad_engine.Close();
    MPI_Finalize();
    return 0;
}

