#include <chrono>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <adios2.h>
#include <vtkm/filter/MarchingCubes.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>
//#include "vtkh.hpp"

vtkm::cont::DataSet MakeIsoSurfaceDataSet(const adios2::Variable<double> &varField,
                   const std::vector<double> &field) {

  vtkm::cont::DataSet dataSet;
  
  // Convert field values to vtkm dataset
  std::string fieldName = "pointvar";

  int nfields = field.size();
  //vtkm::Id = nfields;
  vtkm::cont::ArrayHandle<vtkm::Float64> fieldArray = vtkm::cont::make_ArrayHandle(field.data(), nfields);

  // Global dimensions or local dimensions??
  vtkm::Id3 dims(varField.Count()[0], varField.Count()[1],
                            varField.Count()[2]);

  const vtkm::Id3 vdims(dims[0], dims[1], dims[2]);

  vtkm::Vec<vtkm::FloatDefault, 3> origin(varField.Start()[0], varField.Start()[1],
                                          varField.Start()[2]);

  vtkm::Vec<vtkm::FloatDefault, 3> spacing(1.0f, 1.0f, 1.0f);

  vtkm::cont::ArrayHandleUniformPointCoordinates coordinates(vdims, origin, spacing);

  dataSet.AddCoordinateSystem(vtkm::cont::CoordinateSystem("coordinates", coordinates));

  dataSet.AddField(
    vtkm::cont::Field("pointvar", vtkm::cont::Field::Association::POINTS, fieldArray));

  static const vtkm::IdComponent ndim = 3;
  vtkm::cont::CellSetStructured<ndim> cellSet("cells");
  cellSet.SetPointDimensions(vdims);
  dataSet.AddCellSet(cellSet);

  return dataSet;
}

vtkm::cont::DataSet
compute_isosurface(const adios2::Variable<double> &varField,
                   const std::vector<double> &field, double isovalue)
{   

    std::string fieldName = "pointvar";

    vtkm::cont::DataSet dataSet = MakeIsoSurfaceDataSet(varField, field);
    vtkm::Range range;

    dataSet.GetPointField(fieldName).GetRange(&range);
    //vtkm::Float64 isovalue = isovalue; //range.Center();

    // Create an isosurface filter
    vtkm::filter::MarchingCubes filter;
    filter.SetIsoValue(0, isovalue);
    filter.SetGenerateNormals(true);
    filter.SetActiveField(fieldName);
    filter.SetFieldsToPass({ "pointvar" });
    
    vtkm::cont::DataSet outputData = filter.Execute(dataSet);
    
//    vtkm::cont::DataSet outputData = dataSet; //filter.Execute(dataSet);
    return outputData;
}

void write_vtk(const std::string &fname,
               const vtkm::cont::DataSet data)
{
    vtkm::io::writer::VTKDataSetWriter writer(fname); 
    writer.WriteDataSet(data);
}

void write_adios(adios2::Engine &writer,
                 adios2::Engine &stat_writer,
                 const vtkm::cont::DataSet dataSet,
                 adios2::Variable<double> &varPoint,
                 adios2::Variable<int> &varCell,
                 adios2::Variable<double> &varNormal,
                 adios2::Variable<double> &varScalar,
                 adios2::Variable<int> &varOutStep, int step, MPI_Comm comm)
{
    
    int rank;

    MPI_Comm_rank(comm, &rank);

    vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float64,3>> verticesArray, normalsArray;
    vtkm::cont::DynamicCellSet cellsArray;
    vtkm::cont::ArrayHandle<vtkm::Float64> scalarsArray;

    using VertType = vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float64,3>>;
    vtkm::cont::CoordinateSystem coords = dataSet.GetCoordinateSystem();
    verticesArray = coords.GetData().Cast<VertType>();
    normalsArray = dataSet.GetField("normals").GetData().Cast<VertType>();
    scalarsArray = dataSet.GetField("pointvar").GetData().Cast<vtkm::cont::ArrayHandle<vtkm::Float64>>();

    cellsArray = dataSet.GetCellSet();
    //vtkm::cont::CellSetSingleType<> cellsA = cellsArray.Cast<vtkm::cont::CellSetSingleType<>>(); 
    auto cellsA = cellsArray.GetCellSetBase(); 

    auto readPortal = verticesArray.GetPortalConstControl(); 
    int numPoints = readPortal.GetNumberOfValues();
    int numCells = cellsArray.GetNumberOfCells();
    std::vector<double> points(numPoints * 3);
    std::vector<double> normals(numPoints * 3);
    std::vector<double> scalars(numPoints);
    std::vector<int> cells(numCells * 3); // Assumes that cells are triangles
 
    //Extract point coordinates
    for( vtkm::Id index = 0; index < readPortal.GetNumberOfValues(); index ++) {
      auto vector = readPortal.Get(index);
      points[index * 3 + 0] = vector[0];
      points[index * 3 + 1] = vector[1];
      points[index * 3 + 2] = vector[2];
    }

    auto readPortal1 = normalsArray.GetPortalConstControl();
    // Extract normals
    for( vtkm::Id index = 0; index < readPortal1.GetNumberOfValues(); index ++) {
        auto vector = readPortal1.Get(index);
        normals[index * 3 + 0] = vector[0];
        normals[index * 3 + 1] = vector[1];
        normals[index * 3 + 2] = vector[2];
    }
    
    auto readPortal2 = scalarsArray.GetPortalConstControl();
    for( vtkm::Id index = 0; index < readPortal2.GetNumberOfValues(); index ++) {
        scalars[index] = readPortal2.Get(index);
    }

    for( vtkm::Id index = 1; index < cellsA->GetNumberOfCells(); index ++) {
        
        int npts = cellsA->GetNumberOfPointsInCell(index);
        //vtkm::cont::ArrayHandle<vtkm::Id> ptids; 
        vtkm::Id* ptids = new vtkm::Id(npts);
        //cellsA.GetIndices(index, ptids);
        /* 
        auto readPortal3 = ptids.GetPortalConstControl();
        cells[index * 3 + 0 ] =  readPortal3.Get(0);
        cells[index * 3 + 1 ] =  readPortal3.Get(1);
        cells[index * 3 + 2 ] =  readPortal3.Get(2);
        */
        cellsA->GetCellPointIds(index, ptids);
        cells[index * 3 + 0 ] =  ptids[0];
        cells[index * 3 + 1 ] =  ptids[1];
        cells[index * 3 + 2 ] =  ptids[2];
    }

    int totalPoints, offsetPoints;
    MPI_Allreduce(&numPoints, &totalPoints, 1, MPI_INT, MPI_SUM, comm);
    MPI_Scan(&numPoints, &offsetPoints, 1, MPI_INT, MPI_SUM, comm);

    writer.BeginStep();

    varPoint.SetShape({static_cast<size_t>(totalPoints),
                       static_cast<size_t>(totalPoints > 0 ? 3 : 0)});
    varPoint.SetSelection({{static_cast<size_t>(offsetPoints - numPoints), 0},
                           {static_cast<size_t>(numPoints),
                            static_cast<size_t>(numPoints > 0 ? 3 : 0)}});

    varNormal.SetShape(varPoint.Shape());
    varNormal.SetSelection({varPoint.Start(), varPoint.Count()});

    varScalar.SetShape({static_cast<size_t>(totalPoints)});
    varScalar.SetSelection({{static_cast<size_t>(offsetPoints - numPoints)},
                           {static_cast<size_t>(numPoints)}});

    if (numPoints) {
        writer.Put(varPoint, points.data());
        writer.Put(varNormal, normals.data());
        writer.Put(varScalar, scalars.data());
    }

    int totalCells, offsetCells;
    MPI_Allreduce(&numCells, &totalCells, 1, MPI_INT, MPI_SUM, comm);
    MPI_Scan(&numCells, &offsetCells, 1, MPI_INT, MPI_SUM, comm);

    for (int i = 0; i < cells.size(); i++) {
        cells[i] += (offsetPoints - numPoints);
    }

    varCell.SetShape({static_cast<size_t>(totalCells),
                      static_cast<size_t>(totalCells > 0 ? 3 : 0)});
    varCell.SetSelection({{static_cast<size_t>(offsetCells - numCells), 0},
                          {static_cast<size_t>(numCells),
                           static_cast<size_t>(numCells > 0 ? 3 : 0)}});

    if (numCells) {
        writer.Put(varCell, cells.data());
    }

    if (!rank) {
        std::cout << "isosurface at step " << step << " writing out "
                  << totalCells << " cells and " << totalPoints << " points"
                  << std::endl;
    }

    writer.Put(varOutStep, step);

    writer.EndStep();

    stat_writer.BeginStep();

    stat_writer.Put("step", step);

    stat_writer.EndStep();

}

std::chrono::milliseconds
diff(const std::chrono::steady_clock::time_point &start,
     const std::chrono::steady_clock::time_point &end)
{
    auto diff = end - start;

    return std::chrono::duration_cast<std::chrono::milliseconds>(diff);
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int rank, procs, wrank;

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    const unsigned int color = 5;
    MPI_Comm comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, wrank, &comm);

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &procs);

    int dims[3] = {0};
    MPI_Dims_create(procs, 3, dims);
    size_t npx = dims[0];
    size_t npy = dims[1];
    size_t npz = dims[2];

    int coords[3] = {0};
    int periods[3] = {0};
    MPI_Comm cart_comm;
    MPI_Cart_create(comm, 3, dims, periods, 0, &cart_comm);
    MPI_Cart_coords(cart_comm, rank, 3, coords);
    size_t px = coords[0];
    size_t py = coords[1];
    size_t pz = coords[2];

    if (argc < 4) {
        if (rank == 0) {
            std::cerr << "Too few arguments" << std::endl;
            std::cout << "Usage: isosurface input output isovalue" << std::endl;
        }
	    MPI_Abort(MPI_COMM_WORLD, -1);
    }

    const std::string input_fname(argv[1]);
    const std::string output_fname(argv[2]);
    const double isovalue = std::stod(argv[3]);

    //std::cout << vtkh::AboutVTKH(rank); 

    adios2::ADIOS adios("adios2.xml", comm, adios2::DebugON);

    adios2::IO inIO = adios.DeclareIO("SimulationOutput");
    adios2::Engine reader = inIO.Open(input_fname, adios2::Mode::Read);

    adios2::IO outIO = adios.DeclareIO("IsosurfaceOutput");
    adios2::Engine writer = outIO.Open(output_fname, adios2::Mode::Write);
    auto varOutStep = outIO.DefineVariable<int>("step");

    adios2::IO statIO = adios.DeclareIO("StatisticOutput");
    adios2::Engine stat_writer = statIO.Open("stat.bp",  adios2::Mode::Write);
    auto varStatStep = statIO.DefineVariable<int>("step");

    auto varPoint =
        outIO.DefineVariable<double>("point", {1, 3}, {0, 0}, {1, 3});
    auto varCell = outIO.DefineVariable<int>("cell", {1, 3}, {0, 0}, {1, 3});
    auto varNormal =
        outIO.DefineVariable<double>("normal", {1, 3}, {0, 0}, {1, 3});
    auto varScalar =
        outIO.DefineVariable<double>("scalar", {1}, {0}, {1});

    std::vector<double> u;
    int step;

    auto start_total = std::chrono::steady_clock::now();

    std::ofstream log("isosurface.log");
    log << "step\tread_iso\tcompute_iso\twrite_iso" << std::endl;
    std::cout<<"starting iterations"<<std::endl;
    int i = 0;
    while (true) {
        auto start_step = std::chrono::steady_clock::now();
    
        adios2::StepStatus status = reader.BeginStep();

        if (status != adios2::StepStatus::OK) {
            break;
        }

        adios2::Variable<double> varU = inIO.InquireVariable<double>("U");
        const adios2::Variable<int> varStep = inIO.InquireVariable<int>("step");

        adios2::Dims shape = varU.Shape();

        size_t size_x = (shape[0] + npx - 1) / npx;
        size_t size_y = (shape[1] + npy - 1) / npy;
        size_t size_z = (shape[2] + npz - 1) / npz;

        size_t offset_x = size_x * px;
        size_t offset_y = size_y * py;
        size_t offset_z = size_z * pz;

        if (px == npx - 1) {
            size_x -= size_x * npx - shape[0];
        }
        if (py == npy - 1) {
            size_y -= size_y * npy - shape[1];
        }
        if (pz == npz - 1) {
            size_z -= size_z * npz - shape[2];
        }

        varU.SetSelection({{offset_x, offset_y, offset_z},
                           {size_x + (px != npx - 1 ? 1 : 0),
                            size_y + (py != npy - 1 ? 1 : 0),
                            size_z + (pz != npz - 1 ? 1 : 0)}});

        reader.Get<double>(varU, u);
        reader.Get<int>(varStep, step);
        reader.EndStep();

        auto end_read = std::chrono::steady_clock::now();

        auto dataSet = compute_isosurface(varU, u, isovalue);

        auto end_compute = std::chrono::steady_clock::now();
        std::cout<<"read step" << step <<std::endl;
        write_adios(writer, stat_writer, dataSet, varPoint, varCell, varNormal, varScalar, varOutStep,
                    step, comm);

        auto end_step = std::chrono::steady_clock::now();

        log << step << "\t" << diff(start_step, end_read).count() << "\t"
            << diff(end_read, end_compute).count() << "\t"
            << diff(end_compute, end_step).count() << std::endl;

        std::ifstream ifile("kill_run");
        if (ifile) {
            writer.Close();
            reader.Close();
            stat_writer.Close();
            MPI_Finalize();
            return 0;
        }
     }

    log.close();

    auto end_total = std::chrono::steady_clock::now();

    // std::cout << "Total runtime: " << diff(start_total, end_total).count()
    //           << " [ms]" << std::endl;

    writer.Close();
    reader.Close();

    MPI_Finalize();
}

