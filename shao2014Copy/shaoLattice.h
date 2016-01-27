/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2012 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
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

/* The original version of this file was written by Minh Do-Quang
 */

/** \file
 *  Lattice descriptors for shao_shu_huang_chew (PRE 2014)  -- header file
 */
#ifndef SHAO_LATTICES_H
#define SHAO_LATTICES_H

#include "core/globalDefs.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"

namespace plb {

namespace descriptors {

/// Density, Momentum and Force as external scalars
struct ShaoNSExternals2D {
	static const int numScalars = 5;
	static const int numSpecies = 3;

	static const int densityBeginsAt  = 0;
	static const int sizeOfDensity    = 1;

	static const int forceBeginsAt    = 1;
	static const int sizeOfForce      = 2;

	static const int GradDensityBeginsAt  = 4;
	static const int sizeOfGradDensity    = 2;

  
};
struct ShaoNSExternals3D {
	static const int numScalars = 4;
	static const int numSpecies = 3;

	static const int densityBeginsAt  = 0;
	static const int sizeOfDensity    = 1;

	static const int forceBeginsAt    = 1;
	static const int sizeOfForce      = 3;

};

//////////////
struct ShaoPFExternals2D {
	static const int numScalars = 8; //D2D5 only
	static const int numSpecies = 2;

	static const int densityBeginsAt  = 0;
	static const int sizeOfDensity    = 1;

	static const int velocityBeginsAt = 1;
	static const int sizeOfVelocity   = 2;

	static const int fBeginsAt        = 3;

	static const int forceBeginsAt    = 0;
	static const int sizeOfForce      = 0;
};

//////////////
struct ShaoPFExternals3D {
	static const int numScalars = 11; //D3Q7 only
	static const int numSpecies = 2;

	static const int densityBeginsAt  = 0;
	static const int sizeOfDensity    = 1;

	static const int velocityBeginsAt = 1;
	static const int sizeOfVelocity   = 3;

	static const int fBeginsAt         = 4;

	static const int forceBeginsAt    = 0;
	static const int sizeOfForce      = 0;
};

///////////////////
struct ShaoNSBase2D {
	typedef ShaoNSExternals2D ExternalField;
};
struct ShaoNSBase3D {
	typedef ShaoNSExternals3D ExternalField;
};

///////////////////
struct ShaoPFBase2D {
	typedef ShaoPFExternals2D ExternalField;
};
struct ShaoPFBase3D {
	typedef ShaoPFExternals3D ExternalField;
};

////////////////////////////////////////////////////////////////////////////////////////
/// MRT D2Q5 based lattice. The numbering follows the one in "Viscous flow computations
/// with the method of lattice Boltzmann equation", D. Yu, L.-S. Luo, W. Shi,
/// Progress in Aerospace Sciences 39, (2003), p. 329-367
template <typename T>
struct MRTD2Q5DescriptorBase : public D2Q5DescriptorBase<T>
{
	typedef D2Q5DescriptorBase<T> BaseDescriptor;
	typedef MRTD2Q5DescriptorBase<T> SecondBaseDescriptor;
	enum { numPop=BaseDescriptor::q };

	static const T M[BaseDescriptor::q][BaseDescriptor::q];    // Matrix of base change between f and moments : moments=M.f
	static const T invM[BaseDescriptor::q][BaseDescriptor::q]; // inverse of base change matrix : f=invM.moments
	static const T S[BaseDescriptor::q];       // relaxation times
	enum {jIndexes = 2};
	static const int momentumIndexes[jIndexes]; // relevant indexes of r. t. for shear viscosity
	enum {shearIndexes = 2};
	static const int shearViscIndexes[shearIndexes]; // relevant indexes of r. t. for shear viscosity
	enum {qIndexes = 2};
	static const int qViscIndexes[qIndexes]; // relevant indexes of r. t. for q
	static const int bulkViscIndex  = 1; // relevant index of r. t. for bulk viscosity
	static const int epsilonIndex   = 2; // relevant index of r. t. for epsilon
};

template<typename T>
struct VariableRoundOffPolicy {
	static int SkordosFactor() {
		return mean;
	}
	static T rhoBar(T rho) {
		return rho - mean;
	}
	static T fullRho(T rhoBar) {
		return rhoBar + mean;
	}
	static T invRho(T rhoBar) {
		return (T)1 / (rhoBar + mean);
	}
	static T rhoMinus1(T rhoBar) {
		return rhoBar;
	}

public:
	static T mean;
};
template<typename T> T VariableRoundOffPolicy<T>::mean;

////////////////////////////////////////////////////////////////////////////////////////
/// D2Q5 lattice for Zheng model (PF)
template <typename T>
struct ShaoD2Q5Descriptor
: public D2Q5Constants<T>, public NoOptimizationRoundOffPolicy<T>, public ShaoPFBase2D
  {
	typedef ShaoD2Q5Descriptor<T> BaseDescriptor;
	enum { numPop=D2Q5Constants<T>::q };
	static const char name[];
  };

template<typename T>
const char ShaoD2Q5Descriptor<T>::name[] = "ShaoD2Q5";

/// D2Q9 lattice for Zheng model (NS)
template <typename T>
struct ShaoD2Q9Descriptor
: public D2Q9Constants<T>, public VariableRoundOffPolicy<T>, public ShaoNSBase2D
  {
	typedef ShaoD2Q9Descriptor<T> BaseDescriptor;
	enum { numPop=D2Q9Constants<T>::q };
	static const char name[];
  };

template<typename T>
const char ShaoD2Q9Descriptor<T>::name[] = "ShaoD2Q9";

///////////////////
/// D3Q7 lattice for Zheng model (PF)
template <typename T>
struct ShaoD3Q7Descriptor
: public D3Q7Constants<T>, public NoOptimizationRoundOffPolicy<T>, public ShaoPFBase3D
  {
	typedef ShaoD3Q7Descriptor<T> BaseDescriptor;
	enum { numPop=D3Q7Constants<T>::q };
	static const char name[];
  };

template<typename T>
const char ShaoD3Q7Descriptor<T>::name[] = "ShaoD3Q7";

/// D3Q19 lattice for Zheng model (NS)
template <typename T>
struct ShaoD3Q19Descriptor
: public D3Q19Constants<T>, public VariableRoundOffPolicy<T>, public ShaoNSBase3D
  {
	typedef ShaoD3Q19Descriptor<T> BaseDescriptor;
	enum { numPop=D3Q19Constants<T>::q };
	static const char name[];
  };

template<typename T>
const char ShaoD3Q19Descriptor<T>::name[] = "ShaoD3Q19";

}  // namespace descriptors

}  // namespace plb

#endif  // ZHENG06_LATTICES_H
