/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IO_ADIOS2.cpp
 *
 *  Created on: Feb 2017
 *      Author: Norbert Podhorszki
 */

#include "IO.h"

#include <iostream>
#include <sstream>
#include <string>

#include <adios2.h>

adios2::ADIOS *ad = nullptr;
adios2::ADIOS *ad2 = nullptr;
adios2::Engine writer;
adios2::Engine writer2;
adios2::Variable<double> varT;
adios2::Variable<double> varT1;
adios2::Variable<unsigned int> varGndx;

IO::IO(const Settings &s, MPI_Comm comm)
{
    //Settings s1 = s;
    ad = new adios2::ADIOS(s.configfile, comm, adios2::DebugON);
    //ad2 = new adios2::ADIOS(s1.configfile, comm, adios2::DebugON);

    adios2::IO io = ad->DeclareIO("SimulationOutput");
    adios2::IO io2 = ad->DeclareIO("SimulationOutput2");
    if (!io.InConfigFile())
    {
        // if not defined by user, we can change the default settings
        // BPFile is the default writer
        io.SetEngine("BPFile");
        io.SetParameters({{"num_threads", "1"}});

        // ISO-POSIX file output is the default transport (called "File")
        // Passing parameters to the transport
        io.AddTransport("File", {{"Library", "POSIX"}});
    }
    // if not defined by user, we can change the default settings
    // BPFile is the default writer
    io2.SetEngine("BPFile");
    io2.SetParameters({{"num_threads", "1"}});

        // ISO-POSIX file output is the default transport (called "File")
        // Passing parameters to the transport
    io2.AddTransport("File", {{"Library", "POSIX"}});

    if (!s.rank)
    {
//        std::cout << "Using " << io.m_EngineType << " engine for output" << std::endl;
    }

    // define T as 2D global array
    varT = io.DefineVariable<double>(
        "T",
        // Global dimensions
        {s.gndx, s.gndy},
        // starting offset of the local array in the global space
        {s.offsx, s.offsy},
        // local size, could be defined later using SetSelection()
        {s.ndx, s.ndy});
        
    varT1 = io2.DefineVariable<double>(
        "T",
        // Global dimensions
        {s.gndx, s.gndy},
        // starting offset of the local array in the global space
        {s.offsx, s.offsy},
        // local size, could be defined later using SetSelection()
        {s.ndx, s.ndy});

    std::stringstream ss;
    ss << s.outputfile << ".tmp";
    writer = io.Open(s.outputfile, adios2::Mode::Write, comm);
    writer2 = io2.Open(ss.str(), adios2::Mode::Write, comm);

    // Some optimization:
    // we promise here that we don't change the variables over steps
    // (the list of variables, their dimensions, and their selections)
    writer.FixedSchedule();
    writer2.FixedSchedule();
}

IO::~IO()
{
    writer.Close();
    writer2.Close();
    delete ad;
}

void IO::write(int step, const HeatTransfer &ht, const Settings &s,
               MPI_Comm comm)
{
    int control_msg;

    if (s.rank == 0) {
      std::cin >> control_msg;
    }

    MPI_Bcast(&control_msg, 1, MPI_INT, 0, comm);
    
    if ( control_msg == 0) { 
    	writer.BeginStep();
    // using Put() you promise the pointer to the data will be intact
    // until the end of the output step.
    // We need to have the vector object here not to destruct here until the end
    // of function.
    	std::vector<double> v = ht.data_noghost();
    	writer.Put<double>(varT, v.data());
    	writer.EndStep();
    } else {
       std::cout<<"Switching to second engine";
       writer2.BeginStep();
    // using Put() you promise the pointer to the data will be intact
    // until the end of the output step.
    // We need to have the vector object here not to destruct here until the end
    // of function.
        std::vector<double> v = ht.data_noghost();
        writer2.Put<double>(varT1, v.data());
        writer2.EndStep();
    }
    
}
