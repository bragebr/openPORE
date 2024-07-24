/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2011-2013 Mathias J. Krause, Thomas Henn, Tim Dornieden
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */

/* cylinder3d.cpp:
 * This example examines a steady flow past a cylinder placed in a channel.
 * The cylinder is offset somewhat from the center of the flow to make the
 * steady-state symmetrical flow unstable. At the inlet, a Poiseuille profile is
 * imposed on the velocity, whereas the outlet implements a Dirichlet pressure
 * condition set by p = 0.
 * Inspired by "Benchmark Computations of Laminar Flow Around
 * a Cylinder" by M.Schäfer and S.Turek. For high resolution, low
 * latticeU, and enough time to converge, the results for pressure drop, drag
 * and lift lie within the estimated intervals for the exact results.
 * An unsteady flow with Karman vortex street can be created by changing the
 * Reynolds number to Re=100.
 * It also shows the usage of the STL-reader and explains how
 * to set boundary conditions automatically.
 */


#include "olb3D.h"
#include "olb3D.hh"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = FLOATING_POINT_TYPE;
// 1. IMPLEMENT NAVIER - STOKES FIELD (NS) USING D3Q19
using DESCRIPTOR = D3Q7<VELOCITY,OMEGA>;

// 2. FIND GOOD WAY TO DO SEMI COUPLING OF FIELDS
// 3. TEMPLATES FOR SIMPLE REACTIONS
// 4. GEOMETRY UPDATES - RENAME FUNCTION?
using MOMENTA = momenta::AdvectionDiffusionBulkTuple;
using MomentaF = typename momenta::Tuple<
  momenta::SourcedDensity<typename MOMENTA::density>,
  momenta::FixedVelocityMomentum,
  typename MOMENTA::stress,
  typename MOMENTA::definition>;
using ReactionADBulkDynamics = dynamics::Tuple<T, DESCRIPTOR, MomentaF, equilibria::FirstOrder, collision::BGK>;

// Parameters for the simulation setup
const int N = 5;        // resolution of the model
const T Pe = 1.;       // Peclet number
const T physDiffusivity = 0.01; // diffusivity
const T physLength = 0.1;
const T maxPhysT = 8.; // max. simulation time in s, SI unit
const T relaxationTime = 0.53; // tau
const T Ceq = 50.; // eq. concentration at cylinder boundary

// Stores data from stl file in geometry in form of material numbers
void prepareGeometry( AdeUnitConverter<T,DESCRIPTOR> const& converter,
                      IndicatorF3D<T>& indicator,
                      STLreader<T>& stlReader,
                      SuperGeometry<T,3>& superGeometry )
{

  // === MATERIAL NUMBERS === //
  // solid: MN = 0            //
  // bulk: MN = 1             //
  // wall: MN = 2             //
  // inlet: MN = 3            //
  // outlet: MN = 4           //
  // cylinder: MN = 5         //
  // ======================== //

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0,2,indicator );
  superGeometry.rename( 2,1,stlReader );
  superGeometry.clean();

  Vector<T,3> origin = superGeometry.getStatistics().getMinPhysR( 2 );
  origin[1] += converter.getConversionFactorLength()/2.;
  origin[2] += converter.getConversionFactorLength()/2.;

  Vector<T,3> extend = superGeometry.getStatistics().getMaxPhysR( 2 );
  extend[1] = extend[1]-origin[1]-converter.getConversionFactorLength()/2.;
  extend[2] = extend[2]-origin[2]-converter.getConversionFactorLength()/2.;

  // Set material number for inflow
  origin[0] = superGeometry.getStatistics().getMinPhysR( 2 )[0]-converter.getConversionFactorLength();
  extend[0] = 2*converter.getConversionFactorLength();
  IndicatorCuboid3D<T> inflow( extend,origin );
  superGeometry.rename( 2,3,inflow );

  // Set material number for outflow
  origin[0] = superGeometry.getStatistics().getMaxPhysR( 2 )[0]-converter.getConversionFactorLength();
  extend[0] = 2*converter.getConversionFactorLength();
  IndicatorCuboid3D<T> outflow( extend,origin );
  superGeometry.rename( 2,4,outflow );

  // Set material number for cylinder
  origin[0] = superGeometry.getStatistics().getMinPhysR( 2 )[0]+converter.getConversionFactorLength();
  extend[0] = ( superGeometry.getStatistics().getMaxPhysR( 2 )[0]-superGeometry.getStatistics().getMinPhysR( 2 )[0] )/2.;
  std::shared_ptr<IndicatorF3D<T>> cylinder = std::make_shared<IndicatorCuboid3D<T>>( extend, origin );
  superGeometry.rename( 2,5, cylinder );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice( SuperLattice<T,DESCRIPTOR> &sLattice,
                     AdeUnitConverter<T,DESCRIPTOR> const& converter,
                     STLreader<T>& stlReader,
                     SuperGeometry<T,3>& superGeometry )
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1 -->bulk dynamics
  auto bulkIndicator = superGeometry.getMaterialIndicator({1,3,4,5});
  sLattice.template defineDynamics<NoDynamics>(superGeometry,0);
  sLattice.defineDynamics<AdvectionDiffusionBGKdynamics>(bulkIndicator);

  setBounceBackBoundary(sLattice, superGeometry, 2);

  // Material=5 -->reactive Robin
  setRobinBoundary<T,DESCRIPTOR>(sLattice, omega,
                                 superGeometry.getMaterialIndicator({5}));
  T reactionRate = converter.getCharLatticeVelocity();
  AnalyticalConst3D<T,T> coefficients(reactionRate,
                                      -converter.getLatticeDiffusivity(),
                                      reactionRate * Ceq);
  sLattice.template defineField<descriptors::G>(superGeometry.getMaterialIndicator({5}),
                                                (coefficients));
  // Material=4 -->zero field at outlet
  setRobinBoundary<T, DESCRIPTOR>(sLattice, omega,
                                  superGeometry.getMaterialIndicator({4}));
  AnalyticalConst3D<T,T> coefficients2(0., 1., 0.);
  sLattice.template defineField<descriptors::G>(superGeometry.getMaterialIndicator({4}),
                                                (coefficients2));

  setBounceBackBoundary(sLattice, superGeometry, 2);
  //setInterpolatedVelocityBoundary(sLattice, omega, superGeometry, 3);


  AnalyticalConst3D<T,T> u0(converter.getCharLatticeVelocity(),
                            converter.getCharLatticeVelocity(),
                            converter.getCharLatticeVelocity());
  AnalyticalConst3D<T,T> rho0(0.);

  sLattice.defineField<descriptors::VELOCITY>(bulkIndicator, u0);
  sLattice.defineRho(bulkIndicator, rho0);
  sLattice.iniEquilibrium(bulkIndicator, rho0, u0);

  sLattice.setParameter<descriptors::OMEGA>(omega);
  sLattice.initialize();


  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setBoundaryValues( SuperLattice<T, DESCRIPTOR>& sLattice,
                        AdeUnitConverter<T,DESCRIPTOR> const& converter, int iT,
                        SuperGeometry<T,3>& superGeometry )
{
  OstreamManager clout( std::cout,"setBoundaryValues" );

  // No of time steps for smooth start-up
  int iTmaxStart = converter.getLatticeTime( maxPhysT*0.4 );
  int iTupdate = 30;

  if ( iT%iTupdate == 0 && iT <= iTmaxStart ) {
    // Smooth start curve, sinus
    // SinusStartScale<T,int> StartScale(iTmaxStart, T(1));

    // Smooth start curve, polynomial
    PolynomialStartScale<T,int> StartScale( iTmaxStart, T( 1 ) );

    // Creates and sets the Poiseuille inflow profile using functors
    int iTvec[1] = {iT};
    T frac[1] = {};
    StartScale( frac,iTvec );
    std::vector<T> maxVelocity( 3,0 );
    maxVelocity[0] = 2.25*frac[0]*converter.getCharLatticeVelocity();

    T distance2Wall = converter.getConversionFactorLength()/2.;
    RectanglePoiseuille3D<T> poiseuilleU( superGeometry, 3, maxVelocity, distance2Wall, distance2Wall, distance2Wall );
    sLattice.defineU( superGeometry, 3, poiseuilleU );

    clout << "step=" << iT << "; maxVel=" << maxVelocity[0] << std::endl;

    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

// Computes the pressure drop between the voxels before and after the cylinder
void getResults( SuperLattice<T, DESCRIPTOR>& sLattice,
                 AdeUnitConverter<T,DESCRIPTOR> const& converter,
                 std::size_t iT,
                 SuperGeometry<T,3>& superGeometry,
                 util::Timer<T>& timer,
                 STLreader<T>& stlReader )
{

  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter3D<T> vtkWriter( "reactingCylinder" );
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
  SuperLatticeDensity3D<T, DESCRIPTOR> concentration(sLattice);
  vtkWriter.addFunctor( concentration );
  concentration.getName() = "concentration";

  vtkWriter.addFunctor( velocity );
  vtkWriter.addFunctor( pressure );


  T vtkIter  = 0.04;
  T statIter = 0.4;

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );
    vtkWriter.write( geometry );
    vtkWriter.write( cuboid );
    vtkWriter.write( rank );
    vtkWriter.addFunctor(geometry);

    vtkWriter.createMasterFile();
    vtkWriter.write(iT);

    timer.update(iT);
    timer.printStep();
  }

  // Writes output on the console
  else if ( iT % converter.getLatticeTime(statIter) == 0 ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    // Timer console output
    timer.update( iT );
    timer.printStep();
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
    vtkWriter.write(iT);


  }

  // Writes the vtk files
  else if ( iT% converter.getLatticeTime(vtkIter) == 0 ) {
    vtkWriter.write( iT );

    /*{
      SuperEuklidNorm3D<T> normVel( velocity );
      BlockReduction3D2D<T> planeReduction( normVel, Vector<T,3>({0, 0, 1}) );
      // write output as JPEG
      heatmap::write(planeReduction, iT);
    }*/
  }
}

int main( int argc, char* argv[] )
{

  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );
  // display messages from every single mpi process
  //clout.setMultiOutput(true);

  AdeUnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    N,              // resolution: number of voxels per charPhysL
    relaxationTime,           // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    physLength,            // charPhysLength: reference length of simulation geometry
    Pe*physDiffusivity/physLength,            // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    physDiffusivity, // physDiffusivity
    1.0             // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("reactingCylinder");


  // === 2nd Step: Prepare Geometry ===

  // Instantiation of the STLreader class
  // file name, voxel size in meter, stl unit in meter, outer voxel no., inner voxel no.
  STLreader<T> stlReader( "reactingCylinder.stl", converter.getConversionFactorLength(), 0.001 );
  IndicatorLayer3D<T> extendedDomain( stlReader, converter.getConversionFactorLength() );

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif
  CuboidGeometry3D<T> cuboidGeometry( extendedDomain, converter.getConversionFactorLength(), noOfCuboids );
  cuboidGeometry.setPeriodicity(false,true,true);
  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry<T,3> superGeometry( cuboidGeometry, loadBalancer );

  prepareGeometry( converter, extendedDomain, stlReader, superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

  //prepareLattice and set boundaryCondition
  prepareLattice( sLattice, converter, stlReader, superGeometry );

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for (std::size_t iT = 0; iT < converter.getLatticeTime( maxPhysT ); ++iT) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    //setBoundaryValues( sLattice, converter, iT, superGeometry );

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT, superGeometry, timer, stlReader );
  }

  timer.stop();
  timer.printSummary();
}
