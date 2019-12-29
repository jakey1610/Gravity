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

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <limits>
#include <iomanip>


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

  // std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;

  if (tPlotDelta<=0.0) {
    // std::cout << "plotting switched off" << std::endl;
    tPlot = tFinal + 1.0;
  }
  else {
    // std::cout << "plot initial setup plus every " << tPlotDelta << " time units" << std::endl;
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

  for(int j=0; j<NumberOfBodies; ++j){
      x[j][0] = x[j][0] + timeStepSize * v[j][0];
      x[j][1] = x[j][1] + timeStepSize * v[j][1];
      x[j][2] = x[j][2] + timeStepSize * v[j][2];
      v[j][0] = v[j][0] + timeStepSize * force0[j] / mass[j];
      v[j][1] = v[j][1] + timeStepSize * force1[j] / mass[j];
      v[j][2] = v[j][2] + timeStepSize * force2[j] / mass[j];
      maxV = std::max(maxV, std::sqrt(v[j][0] * v[j][0] + v[j][1] * v[j][1] + v[j][2] * v[j][2]));
  }

  //Update the velocity of each item if collision occurs
  for (int j=NumberOfBodies-1; j>0; --j){
    for(int i=0; i<j; ++i){
      const double distance = sqrt(
        (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
        (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
        (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
      );
      if (distance <= 0.01){ //timeStepSize*std::sqrt(v[j][0] * v[j][0] + v[j][1] * v[j][1] + v[j][2] * v[j][2])
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
  if (NumberOfBodies==1){
    //(2.00224819850178, 0.998171674255814)
    //5e-9: 20022481985.0178   9981716742.55814
    //1e-10: 20022482311.0063   9981716781.97356
    // std::cout << "Position: (" << x[0][0] << ", " << x[0][1] << ", " << x[0][2] << ")" << std::endl;
    std::cout << timeStepSize <<" "<<sqrt(pow((20022482277.631-pow(10,10)*x[0][0]),2) + pow((9981716775.90204-pow(10,10)*x[0][1]),2))/pow(10,10) << std::endl;
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
    // printParaviewSnapshot();
    // std::cout << "plotted initial setup" << std::endl;
    tPlot = tPlotDelta;
  }

  int timeStepCounter = 0;
  while (t<=tFinal) {
    updateBody();
    timeStepCounter++;
    if (t >= tPlot) {
      // printParaviewSnapshot();
      // std::cout << "plot next snapshot"
    	// 	    << ",\t time step=" << timeStepCounter
    	// 	    << ",\t t="         << t
			// 	<< ",\t dt="        << timeStepSize
			// 	<< ",\t v_max="     << maxV
			// 	<< ",\t dx_min="    << minDx
			// 	<< std::endl;

      tPlot += tPlotDelta;
    }
  }

  closeParaviewVideoFile();

  return 0;
}
