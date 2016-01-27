/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011 FlowKit Sarl
 * Avenue de Chailly 23
 * 1012 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "palabos3D.h"
#include "palabos3D.hh"
#include <cstdlib>
#include <iostream>

#include "Zheng06Dynamics.cpp"
#include "Zheng2006Processor.cpp"
#include "Zheng06Lattice.h"

using namespace plb;
using namespace slilab;
using namespace std;

// Use double-precision arithmetics
typedef double T;
#define DESCRIPTOR   descriptors::ZhengD3Q19Descriptor
#define DESCRIPTOR_F descriptors::ZhengD3Q7Descriptor

/// Initial condition: heavy fluid at center, light fluid at ....
/** This functional is going to be used as an argument to the function "applyIndexed",
 *  to setup the initial condition. For efficiency reasons, this approach should
 *  always be preferred over explicit space loops in end-user codes.
 */
template<typename T, template<typename U> class Descriptor>
class DropletInitializer : public BoxProcessingFunctional3D_L<T,Descriptor>
{
public :
	DropletInitializer(plint nx_, plint ny_, plint nz_, T radius, T cn,
			T a, T k, T pc) :
				nx(nx_), ny(ny_), nz(nz_),
				r2(radius*radius),
				width(cn),
				aCoef(a),
				kappa(k),
				phiCoef(pc)
{ };
	virtual void process(Box3D domain,BlockLattice3D<T,Descriptor>& lattice)
	{
		Array<T,3> zeroVelocity(0,0,0);
		Dot3D Offset = lattice.getLocation();
		for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
			for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
				for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ)
				{
					const plint gX = iX+Offset.x;
					const plint gY = iY+Offset.y;
					const plint gZ = iZ+Offset.z;

					T rx = gX-nx/2;
					T ry = gY-ny/2;
					T rz = gZ-nz/2;

					T rSqr = sqrt(rx*rx + ry*ry + rz*rz);
					T phi = phiCoef*tanh( (T)2*(rSqr-sqrt(r2))/width );
					T lapPhi = 0;

					T mu = aCoef*( (T)4*phi*phi*phi - (T)4*phiCoef*phiCoef*phi ) - kappa * lapPhi;

					T phiBar = Descriptor<T>::rhoBar(phi);
					Array< T, SymmetricTensor< T, Descriptor >::n > PiNeq;

					lattice.get(iX,iY,iZ).regularize( phiBar, zeroVelocity, (T)0, PiNeq);

					lattice.get(iX,iY,iZ).setExternalField (
							Descriptor<T>::ExternalField::velocityBeginsAt,
							Descriptor<T>::ExternalField::sizeOfVelocity, &zeroVelocity[0] );
				}
			}
		}
	};
	virtual DropletInitializer<T,Descriptor>*
	clone() const
	{
		return new DropletInitializer<T,Descriptor>(*this);
	};
	virtual void getModificationPattern(std::vector<bool>& isWritten) const
	{
		isWritten[0] = true;
	};
	virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const
	{
		modified[0] = modif::staticVariables;
	}
	virtual bool includeEnvelope() const
	{
		return false;
	};

private :
	T width, r2;
	plint nx, ny, nz;
	T aCoef, kappa, phiCoef;
};


int main(int argc, char *argv[])
{
	plbInit(&argc, &argv);
	global::directories().setOutputDir("./tmp/");

	const T omegaF = 1.0/0.7;
	const T omegaP = 1.0/0.7;
	const plint nx   = 64;
	const plint ny   = 128;
	const plint nz   = 64;
	const plint radius = 10;
//	const plint maxIter  = 20000;
//	const plint saveIter = 500;
	const plint maxIter  = 5;
	const plint saveIter = 1;
	const plint statIter = 1000;

	const T rho_l = 1.0;
	const T rho_h = 1000;
	const T phiCoef = (rho_h-rho_l)/2.0;
	const T W = 3;
	const T sigma = 0.01;
	const T kappa = (T)3/(T)8*sigma*W/phiCoef/phiCoef;
	const T aCoef = 3./4*sigma/W/pow(phiCoef,4.0);
	//	const T aCoef = 0.001;
	const T Gamma = 100;
	const T qCoef = (T)1/((T)1/omegaP+(T)0.5);
	const T grav = -1e-6*phiCoef;

	Dynamics<T,DESCRIPTOR_F> *cellDynamic = new CahnHilliardDynamics<T, DESCRIPTOR_F>(omegaP);
	cellDynamic->setParameter( slilab::Zheng06::Gamma, Gamma );
	cellDynamic->setParameter( slilab::Zheng06::qCoef, qCoef );

	MultiBlockLattice3D<T, DESCRIPTOR> Fluid ( nx,ny,nz, new Zheng06BGKdynamics<T,DESCRIPTOR>(omegaF) );
	MultiBlockLattice3D<T, DESCRIPTOR_F> PhaseField ( nx,ny,nz, cellDynamic );

	// Checking the lattice type
	PLB_ASSERT(DESCRIPTOR<T>::ExternalField::numScalars);
	PLB_ASSERT(DESCRIPTOR_F<T>::ExternalField::numScalars);

//	OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
//		// = createInterpBoundaryCondition3D<T,DESCRIPTOR>();
//		// = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
//		  = createZouHeBoundaryCondition3D<T,DESCRIPTOR>();
//
//	boundaryCondition->addVelocityBoundary2N( Box3D(0,nx-1, 0,ny-1, 0,0), Fluid );
//	boundaryCondition->addVelocityBoundary2P( Box3D(0,nx-1, 0,ny-1, nz-1,nz-1), Fluid );
//	boundaryCondition->addVelocityBoundary0N( Box3D(0,0, 0,ny-1, 0,nz-1), Fluid );
//	boundaryCondition->addVelocityBoundary0N( Box3D(nx-1,nx-1, 0,ny-1, 0,nz-1), Fluid );

//	boundaryCondition->addExternalVelocityEdge1NN( Box3D(0,0, 0,ny-1, 0,0), Fluid);
//	boundaryCondition->addExternalVelocityEdge1NP( Box3D(0,0, 0,ny-1, nz-1,nz-1), Fluid);
//	boundaryCondition->addExternalVelocityEdge1PP( Box3D(nx-1,nx-1, 0,ny-1, nz-1,nz-1), Fluid);
//	boundaryCondition->addExternalVelocityEdge1PN( Box3D(nx-1,nx-1, 0,ny-1, 0,0), Fluid);

//	setBoundaryVelocity(Fluid, Fluid.getBoundingBox(), Array<T,3>(0.,0.,0.) );

	defineDynamics( Fluid, Box3D(0,nx-1, 0,ny-1, 0,0), new BounceBack<T,DESCRIPTOR>());
	defineDynamics( Fluid, Box3D(0,nx-1, 0,ny-1, nz-1,nz-1), new BounceBack<T,DESCRIPTOR>());
	defineDynamics( Fluid, Box3D(1,1, 0,ny-1, 0,nz-1), new BounceBack<T,DESCRIPTOR>());
	defineDynamics( Fluid, Box3D(nx-2,nx-2, 0,ny-1, 0,nz-1), new BounceBack<T,DESCRIPTOR>());

	// Make periodic.
	Fluid.periodicity().toggle(1, true);
	PhaseField.periodicity().toggleAll(true);

	initializeAtEquilibrium(Fluid, Fluid.getBoundingBox(), phiCoef+1, Array<T,3>(0.,0.,0.) );

	plint processorLevel = 1;
	integrateProcessingFunctional (
			new Zheng2006Processor3D<T,DESCRIPTOR,T,DESCRIPTOR_F>(aCoef, kappa, phiCoef, grav),
			Fluid.getBoundingBox(),
			Fluid,
			PhaseField,
			processorLevel );

	// laplaceSetup(Fluid, r_phase1, r_phase2, radius);
	applyProcessingFunctional(new DropletInitializer<T,DESCRIPTOR_F>(nx, ny, nz, radius, W, aCoef, kappa, phiCoef),
			PhaseField.getBoundingBox(), PhaseField);

	PhaseField.initialize();
	Fluid.initialize(); // .executeInternalProcessors(1);
	Fluid.executeInternalProcessors();

	pcout << "Starting simulation" << endl;
	// Main loop over time iterations.
	for (int iT=0; iT<maxIter; ++iT) {
		if (iT%saveIter==0) {
			VtkImageOutput3D<T> vtkOut(createFileName("data", iT, 6), 1.);
			vtkOut.writeData<float>(*computeDensity( PhaseField ), "phi", 1.);
			vtkOut.writeData<float>(*computeDensity( Fluid ), "density", 1.);
			vtkOut.writeData<3,float>(*computeVelocity( Fluid ), "velocity", 1.);
		}

		// Time iteration for the light fluid.
		PhaseField.collideAndStream();
		Fluid.collideAndStream();
	}
}

