
#include <iostream>
#include <mpi.h>
#include <vector>
#include <string>
#include <adios2.h>
#include <vtkm/cont/DataSetBuilderExplicit.h>
#include <vtkm/rendering/Camera.h>
#include <vtkm/rendering/CanvasRayTracer.h>
#include <vtkm/rendering/MapperRayTracer.h>
#include <vtkm/rendering/Scene.h>
#include <vtkm/rendering/View3D.h>
//#include "vtkh.hpp"

typedef struct {
    adios2::IO *inIO;
    adios2::Engine *reader;
    int nprocs;
    int rank;
    adios2::Engine *writer;
} Context;

MPI_Comm comm;

void render(const vtkm::cont::DataSet& outputData, std::string outfile) 
{


// compute the bounds and extends of the output data
  vtkm::Bounds coordsBounds = outputData.GetCoordinateSystem().GetBounds();
  vtkm::Vec<vtkm::FloatDefault,3> totalExtent( coordsBounds.X.Length(),
                                          coordsBounds.Y.Length(),
                                          coordsBounds.Z.Length() );
  vtkm::FloatDefault mag = vtkm::Magnitude(totalExtent);
  vtkm::Normalize(totalExtent);

  // setup a camera and point it to towards the center of the input data
  vtkm::rendering::Camera camera;
  camera.ResetToBounds(coordsBounds);

  camera.SetLookAt(totalExtent*(mag * .5f));
  camera.SetViewUp(vtkm::make_Vec(0.f, 1.f, 0.f));
  camera.SetClippingRange(1.f, 100.f);
  camera.SetFieldOfView(60.f);
  camera.SetPosition(totalExtent*(mag * 2.f));
  vtkm::cont::ColorTable colorTable("inferno");

 
  // Create a mapper, canvas and view that will be used to render the scene
  vtkm::rendering::Scene scene;
  vtkm::rendering::MapperRayTracer mapper;
  vtkm::rendering::CanvasRayTracer canvas(512, 512);
  vtkm::rendering::Color bg(0.2f, 0.2f, 0.2f, 1.0f);

  // Render an image of the output isosurface
  scene.AddActor(vtkm::rendering::Actor(outputData.GetCellSet(),
                                        outputData.GetCoordinateSystem(),
                                        outputData.GetField("pointvar"),
                                        colorTable));
  vtkm::rendering::View3D view(scene, mapper, canvas, camera, bg);
  view.Initialize();
  view.Paint();
  view.SaveAs(outfile);
  MPI_Barrier(comm);
}

vtkm::cont::DataSet read_mesh(const std::vector<double> &bufPoints,
                                       const std::vector<int> &bufCells,
                                       const std::vector<double> &bufNormals,
                                       const std::vector<double> &bufScalars)
{
    int nPoints = bufPoints.size() / 3;
    int nCells = bufCells.size() / 3;
    
    
    vtkm::cont::DataSetBuilderExplicitIterative  dsBuild;

    for( vtkm::Id index = 0; index < nPoints; index ++) {
          dsBuild.AddPoint((float)bufPoints[index*3 + 0], 
                  (float)bufPoints[index*3 + 1], (float)bufPoints[index*3 + 2]);
    }  


    for( vtkm::Id index = 0; index < nCells; index ++) {
          dsBuild.AddCell(vtkm::CELL_SHAPE_TRIANGLE);
          dsBuild.AddCellPoint(bufCells[index*3 + 0]);
          dsBuild.AddCellPoint(bufCells[index*3 + 1]);
          dsBuild.AddCellPoint(bufCells[index*3 + 2]);
    }  

    vtkm::cont::DataSet  dataSet = dsBuild.Create();

    vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::FloatDefault,3>> normalsArray;
    normalsArray.Allocate(nPoints);
    auto portal1 = normalsArray.GetPortalControl();

    for( vtkm::Id index = 0; index < nPoints; index ++) {
          portal1.Set(index, vtkm::Vec<vtkm::FloatDefault,3>(bufNormals[index*3 + 0], 
                             bufNormals[index*3 + 1], bufNormals[index*3 + 2]));
    }
    
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> scalarsArray;
    scalarsArray.Allocate(nPoints);
    auto portal2 = scalarsArray.GetPortalControl();

    for( vtkm::Id index = 0; index < nPoints; index ++) {
          portal2.Set(index, bufScalars[index]);
    }
    dataSet.AddField(vtkm::cont::Field("normals", 
                         vtkm::cont::Field::Association::POINTS, normalsArray));
    dataSet.AddField(vtkm::cont::Field("pointvar", 
                         vtkm::cont::Field::Association::POINTS, scalarsArray));

    return dataSet;
}

size_t get_size_to_read(int rank, size_t t_sz, int procs) {
    if (rank == procs -1) 
       return (t_sz%procs)?(t_sz/procs):(t_sz/procs+(t_sz%procs));
    else
       return t_sz/procs;
}

bool timer_func(Context *context)
{

    std::vector<double> points;
    std::vector<int> cells;
    std::vector<double> normals;
    std::vector<double> scalars;

    int step;

    adios2::StepStatus status = context->reader->BeginStep();

    if (status != adios2::StepStatus::OK) {
        return false;
    }

    auto varPoint = context->inIO->InquireVariable<double>("point");
    auto varCell = context->inIO->InquireVariable<int>("cell");
    auto varNormal = context->inIO->InquireVariable<double>("normal");
    auto varScalar = context->inIO->InquireVariable<double>("scalar");
    auto varStep = context->inIO->InquireVariable<int>("step");
    
    
    size_t npoints = get_size_to_read(context->rank, varPoint.Shape()[0], context->nprocs);
    size_t ncells = get_size_to_read(context->rank, varCell.Shape()[0], context->nprocs);
    size_t nnorms = get_size_to_read(context->rank, varNormal.Shape()[0], context->nprocs);
    size_t nscalars = get_size_to_read(context->rank, varScalar.Shape()[0], context->nprocs);
  
    size_t pstart, cstart, nstart, sstart;

    if (varPoint.Shape().size() > 0 || varCell.Shape().size() > 0) {
        if (context->rank < context-> nprocs - 1 ) {
            pstart = context->rank * npoints;
            cstart = context->rank * ncells;
            nstart = context->rank * nnorms;
            sstart = context->rank * nscalars;
        } else {
            pstart = context->rank * varPoint.Shape()[0]/context->nprocs;
            cstart = context->rank * varCell.Shape()[0]/context->nprocs;
            nstart = context->rank * varNormal.Shape()[0]/context->nprocs;
            sstart = context->rank * varScalar.Shape()[0]/context->nprocs;
        }    
        varPoint.SetSelection(
            {{pstart, 0}, {npoints, varPoint.Shape()[1]}});
        varCell.SetSelection(
            {{cstart, 0}, {ncells, varCell.Shape()[1]}});
        varNormal.SetSelection(
            {{nstart, 0}, {nnorms, varNormal.Shape()[1]}});
        varScalar.SetSelection(
            {{sstart}, {nscalars}});

        context->reader->Get<double>(varPoint, points);
        context->reader->Get<int>(varCell, cells);
        context->reader->Get<double>(varNormal, normals);
        context->reader->Get<double>(varScalar, scalars);
    }

    context->reader->Get<int>(varStep, &step);

    context->reader->EndStep();

    std::cout << "rendering_isosurface at step " << step << std::endl;

    vtkm::cont::DataSet dataSet = read_mesh(points, cells, normals, scalars);

    std::string outfile = "demo-" + std::to_string(step) + "-" + std::to_string(context->rank) + ".png";

    render(dataSet, outfile); 

    context->writer->BeginStep();

    context->writer->Put("step", step);

    context->writer->EndStep();

    std::ifstream ifile("kill_run");

    if (ifile) {
        // The file exists, and is open for input
        context->writer->Close();
        return false;
    }
       
    return true;
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int rank, procs, wrank;

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    const unsigned int color = 7;
    MPI_Comm_split(MPI_COMM_WORLD, color, wrank, &comm);

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &procs);

    if (argc < 2) {
        if (rank == 0) {
            std::cerr << "Too few arguments" << std::endl;
            std::cout << "Usage: render_isosurface input" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
/**
    if (procs != 1) {
        if (rank == 0) {
            std::cerr << "render_isosurface only supports serial execution"
                      << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
**/
    const std::string input_fname(argv[1]);
    const std::string output_fname("stat.bp");

//    std::cout << vtkh::AboutVTKH(rank);

    adios2::ADIOS adios("adios2.xml", comm, adios2::DebugON);

    adios2::IO inIO = adios.DeclareIO("IsosurfaceOutput");
    adios2::Engine reader = inIO.Open(input_fname, adios2::Mode::Read);

    adios2::IO outIO = adios.DeclareIO("StatisticOutput");
    adios2::Engine writer = outIO.Open(output_fname, adios2::Mode::Write);
    auto varOutStep = outIO.DefineVariable<int>("step");

    Context context = {
        .inIO = &inIO,
        .reader = &reader,
        .nprocs = procs,
        .rank = rank,
        .writer = &writer,
    };

    while(timer_func(&context) == true);
    writer.Close();
    reader.Close();
    MPI_Finalize();
}
