/*
 * ZouHeWettingBCDynamics.cpp
 *
 *  Created on: Mar 4, 2013
 *      Author: minh
 */

#include "ZouHeWettingBCDynamics.h"

namespace slilab {

template<typename T, template<typename U> class Descriptor,int direction, int orientation>
int ZouHeWettingBCDynamics<T,Descriptor,direction,orientation>::id =
		meta::registerGeneralDynamics<T,Descriptor, ZouHeWettingBCDynamics<T,Descriptor,direction,orientation> >
( std::string("Boundary_ZouHeWetting_")+util::val2str(direction) +
		std::string("_")+util::val2str(orientation) );

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
ZouHeWettingBCDynamics<T,Descriptor,direction,orientation>::ZouHeWettingBCDynamics (
		Dynamics<T,Descriptor>* baseDynamics, bool automaticPrepareCollision )
		: DensityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>(baseDynamics, automaticPrepareCollision)
		  { }

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
ZouHeWettingBCDynamics<T,Descriptor,direction,orientation>*
ZouHeWettingBCDynamics<T,Descriptor, direction, orientation>::clone() const
{
	return new ZouHeWettingBCDynamics<T,Descriptor,direction,orientation>(*this);
}

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
ZouHeWettingBCDynamics<T,Descriptor,direction,orientation>::
ZouHeWettingBCDynamics(HierarchicUnserializer& unserializer)
: DensityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>(0, false)
  {
	unserialize(unserializer);
  }

template<typename T, template<typename U> class Descriptor,
int direction, int orientation>
void ZouHeWettingBCDynamics<T,Descriptor,direction,orientation>::serialize(HierarchicSerializer& serializer) const
{
	DensityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>::serialize(serializer);
}

template<typename T, template<typename U> class Descriptor,
int direction, int orientation>
void ZouHeWettingBCDynamics<T,Descriptor,direction,orientation>::unserialize(HierarchicUnserializer& unserializer)
{
	DensityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>::unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor,
int direction, int orientation>
int ZouHeWettingBCDynamics<T,Descriptor,direction,orientation>::getId() const {
	return id;
}

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
void ZouHeWettingBCDynamics<T,Descriptor,direction,orientation>::completePopulations(Cell<T,Descriptor>& cell) const
{
	// First, compute feq
	T rhoBar;
	Array<T,Descriptor<T>::d> j;
	this -> computeRhoBarJ(cell, rhoBar, j);
	T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);

	for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
		cell[iPop] = this->getBaseDynamics().computeEquilibrium(iPop, rhoBar, j, jSqr);
	}

	// Then call to the orginal ZouHe function.
	ZouHeClosure<T,Descriptor,direction,orientation>(cell, *this);
}

} /* namespace slilab */
