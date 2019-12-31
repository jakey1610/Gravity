// Translate this file with
//
// g++ -O3 --std=c++11 assignment-2019.c -o assignment
//
// Run it with
//
// ./assignment
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
// (C) 2018-2019 Tobias Weinzierl


//NEED TO SORT THIS PARALLELISATION OUT

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <limits>
#include <iomanip>

#include <time.h>


double t          = 0;
double tFinal     = 0;
double tPlot      = 0;
double tPlotDelta = 0;

int NumberOfBodies = 0;

/**
 * Pointer to pointers. Each pointer in turn points to three coordinates, i.e.
 * each pointer represents one molecule/particle/body.
 */
double** x;

/**
 * Equivalent to x storing the velocities.
 */
double** v;

/**
 * One mass entry per molecule/particle.
 */
double*  mass;

/**
 * Global time step size used.
 */
double   timeStepSize = 0.0;

/**
 * Maximum velocity of all particles.
 */
double   maxV;

/**
 * Minimum distance between two elements.
 */
double   minDx;


/**
 * Set up scenario from the command line.
 *
 * This operation is not to be changed in the assignment.
 */
void setUp(int argc, char** argv) {
  NumberOfBodies = (argc-3) / 7;

  x    = new double*[NumberOfBodies];
  v    = new double*[NumberOfBodies];
  mass = new double [NumberOfBodies];

  int readArgument = 1;

  tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
  tFinal       = std::stof(argv[readArgument]); readArgument++;
  timeStepSize = std::stof(argv[readArgument]); readArgument++;

  for (int i=0; i<NumberOfBodies; i++) {
    x[i] = new double[3];
    v[i] = new double[3];

    x[i][0] = std::stof(argv[readArgument]); readArgument++;
    x[i][1] = std::stof(argv[readArgument]); readArgument++;
    x[i][2] = std::stof(argv[readArgument]); readArgument++;

    v[i][0] = std::stof(argv[readArgument]); readArgument++;
    v[i][1] = std::stof(argv[readArgument]); readArgument++;
    v[i][2] = std::stof(argv[readArgument]); readArgument++;

    mass[i] = std::stof(argv[readArgument]); readArgument++;

    if (mass[i]<=0.0 ) {
      std::cerr << "invalid mass for body " << i << std::endl;
      exit(-2);
    }
  }

  std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;

  if (tPlotDelta<=0.0) {
    std::cout << "plotting switched off" << std::endl;
    tPlot = tFinal + 1.0;
  }
  else {
    std::cout << "plot initial setup plus every " << tPlotDelta << " time units" << std::endl;
    tPlot = 0.0;
  }
}


std::ofstream videoFile;


/**
 * This operation is not to be changed in the assignment.
 */
void openParaviewVideoFile() {
  videoFile.open( "result.pvd" );
  videoFile << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
            << "<Collection>";
}





/**
 * This operation is not to be changed in the assignment.
 */
void closeParaviewVideoFile() {
  videoFile << "</Collection>"
            << "</VTKFile>" << std::endl;
}


/**
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 *
 * This operation is not to be changed in the assignment.
 */
void printParaviewSnapshot() {
  static int counter = -1;
  counter++;
  std::stringstream filename;
  filename << "result-" << counter <<  ".vtp";
  std::ofstream out( filename.str().c_str() );
  out << "<VTKFile type=\"PolyData\" >" << std::endl
      << "<PolyData>" << std::endl
      << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
      << "  <Points>" << std::endl
      << "   <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
//      << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

  for (int i=0; i<NumberOfBodies; i++) {
    out << x[i][0]
        << " "
        << x[i][1]
        << " "
        << x[i][2]
        << " ";
  }

  out << "   </DataArray>" << std::endl
      << "  </Points>" << std::endl
      << " </Piece>" << std::endl
      << "</PolyData>" << std::endl
      << "</VTKFile>"  << std::endl;

  videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
}



/**
 * This is the only operation you are allowed to change in the assignment.
 */
void updateBody() {
  maxV   = 0.0;
  minDx  = std::numeric_limits<double>::max();
  // double tPrev = tFinal;
  // tFinal = tFinal * 2;

  // force0 = force along x direction
  // force1 = force along y direction
  // force2 = force along z direction
  auto* force0 = new double[NumberOfBodies]();
  auto* force1 = new double[NumberOfBodies]();
  auto* force2 = new double[NumberOfBodies]();

  for(int j=0; j<NumberOfBodies; ++j){
    for (int i=j+1; i<NumberOfBodies; i++) {
        const double distance = sqrt(
          (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
          (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
          (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
        );

        // x,y,z forces acting on particles
        double c = mass[i]*mass[j] / distance / distance / distance ;

        force0[j] += (x[i][0]-x[j][0]) * c;
        force1[j] += (x[i][1]-x[j][1]) * c;
        force2[j] += (x[i][2]-x[j][2]) * c;
        force0[i] -= (x[i][0]-x[j][0]) * c;
        force1[i] -= (x[i][1]-x[j][1]) * c;
        force2[i] -= (x[i][2]-x[j][2]) * c;
        minDx = std::min(minDx, distance);
    }
  }

  const int numBuckets = 11;

  bool** buckets = new bool*[numBuckets]();
  for (int i=0; i<numBuckets; i++){
    buckets[i] = new bool[NumberOfBodies]();
  }

  int* timeStepDivisor = new int[NumberOfBodies]();
  const double vBucket = maxV / (numBuckets-1);
  double velocity = 0;
  #pragma omp parallel for
  for(int j=0; j<NumberOfBodies; ++j){
      //Sort into buckets
      #pragma omp parallel for
      for (int ii=0; ii<numBuckets-1; ii++){
        velocity = std::sqrt(v[j][0] * v[j][0] + v[j][1] * v[j][1] + v[j][2] * v[j][2]);
        if (velocity >= ii*vBucket && velocity < (ii+1)*vBucket){
          timeStepDivisor[j] = numBuckets-(ii+1);
          buckets[ii][j] = true;
        }
      }
      if (velocity >= (numBuckets+1)*vBucket){
        timeStepDivisor[j] = numBuckets-10;
        buckets[9][j] = true;
      }


  }
  #pragma omp parallel for collapse(2)
  for (int jj = numBuckets-2; jj>0; jj--){
    for (int ii = 0; ii<NumberOfBodies; ii++){
      if(buckets[jj][ii]){
        double timeStepAltered = timeStepSize / timeStepDivisor[ii];
        for (int num = 1; num<jj+1; num++){
          x[ii][0] = x[ii][0] + timeStepAltered * v[ii][0];
          x[ii][1] = x[ii][1] + timeStepAltered * v[ii][1];
          x[ii][2] = x[ii][2] + timeStepAltered * v[ii][2];
          v[ii][0] = v[ii][0] + timeStepAltered * force0[ii] / mass[ii];
          v[ii][1] = v[ii][1] + timeStepAltered * force1[ii] / mass[ii];
          v[ii][2] = v[ii][2] + timeStepAltered * force2[ii] / mass[ii];
        }
        maxV = std::max(maxV, std::sqrt(v[ii][0] * v[ii][0] + v[ii][1] * v[ii][1] + v[ii][2] * v[ii][2]));
      }
    }
  }

  //Update the velocity of each item if collision occurs

  for (int j=NumberOfBodies-1; j>0; --j){
    for(int i=0; i<j; ++i){
      const double distance = sqrt(
        (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
        (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
        (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
      );
      if (distance <= 0.01){
        const double comb = mass[i]+mass[j];
        v[i][0] = mass[i]*(1/comb)*v[i][0] + mass[j]*(1/comb)*v[j][0];
        v[i][1] = mass[i]*(1/comb)*v[i][1] + mass[j]*(1/comb)*v[j][1];
        v[i][2] = mass[i]*(1/comb)*v[i][2] + mass[j]*(1/comb)*v[j][2];

        mass[i] = comb;

        x[i][0] = (x[i][0] + x[j][0]) * 0.5;
        x[i][1] = (x[i][1] + x[j][1]) * 0.5;
        x[i][2] = (x[i][2] + x[j][2]) * 0.5;


        //Swap i for last item in each list then reduce the number of bodies.
        if(j!=NumberOfBodies-1){
          mass[j] = mass[NumberOfBodies-1];
          x[j][0] = x[NumberOfBodies-1][0];
          x[j][1] = x[NumberOfBodies-1][1];
          x[j][2] = x[NumberOfBodies-1][2];
          v[j][0] = v[NumberOfBodies-1][0];
          v[j][1] = v[NumberOfBodies-1][1];
          v[j][2] = v[NumberOfBodies-1][2];
        }
        NumberOfBodies -= 1;

      }
    }
  }
  if (NumberOfBodies<=1){
    std::cout << "Position: (" << x[0][0] << ", " << x[0][1] << ", " << x[0][2] << ")" << std::endl;
    t=tFinal+1;
  }

  t += timeStepSize;
  delete[] force0;
  delete[] force1;
  delete[] force2;
}


/**
 * Main routine.
 *
 * Not to be changed in assignment.
 */
int main(int argc, char** argv) {
  if (argc==1) {
    std::cerr << "usage: " + std::string(argv[0]) + " snapshot final-time dt objects" << std::endl
              << "  snapshot        interval after how many time units to plot. Use 0 to switch off plotting" << std::endl
              << "  final-time      simulated time (greater 0)" << std::endl
              << "  dt              time step size (greater 0)" << std::endl
              << std::endl
              << "Examples:" << std::endl
              << "0.01  100.0   0 0 0 1.0   0   0 1.0 \t One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
              << "0.01  100.0   0 0 0 1.0   0   0 1.0 0 1.0 0 1.0 0   0 1.0 \t One spiralling around the other one" << std::endl
              << "0.01  100.0 3.0 0 0   0 1.0   0 0.4 0   0 0   0 0   0 0.2 2.0 0 0 0 0 0 1.0 \t Three body setup from first lecture" << std::endl
              << std::endl
              << "In this naive code, only the first body moves" << std::endl;

    return -1;
  }
  else if ( (argc-4)%7!=0 ) {
    std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)" << std::endl;
    return -2;
  }

  std::cout << std::setprecision(15);

  setUp(argc,argv);

  openParaviewVideoFile();

  int snapshotCounter = 0;
  if (t > tPlot) {
    printParaviewSnapshot();
    std::cout << "plotted initial setup" << std::endl;
    tPlot = tPlotDelta;
  }

  int timeStepCounter = 0;
  while (t<=tFinal) {
    updateBody();
    timeStepCounter++;
    if (t >= tPlot) {
      printParaviewSnapshot();
      std::cout << "plot next snapshot"
    		    << ",\t time step=" << timeStepCounter
    		    << ",\t t="         << t
				<< ",\t dt="        << timeStepSize
				<< ",\t v_max="     << maxV
				<< ",\t dx_min="    << minDx
				<< std::endl;

      tPlot += tPlotDelta;
    }
  }

  closeParaviewVideoFile();

  return 0;
}
