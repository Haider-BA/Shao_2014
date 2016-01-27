/*
 * ShaoDynamics.h
 *
 *  Created on: Dec 04, 2015
 *      
 */

#ifndef SHAODYNAMICS_H_
#define SHAODYNAMICS_H_

using namespace plb;

namespace slilab {

namespace Shao {
const plint Gamma = 1000;
const plint qCoef = 1001;
}

template<typename T, template<typename U> class Descriptor>
class CahnHilliardDynamics : public AdvectionDiffusionDynamics<T,Descriptor> {
public:
	CahnHilliardDynamics(T omega_);
	CahnHilliardDynamics(HierarchicUnserializer& unserializer);
	CahnHilliardDynamics(CahnHilliardDynamics<T,Descriptor> const &rhs);

	/// Serialize the dynamics object. Including the private variables
	virtual void serialize ( HierarchicSerializer & serializer ) const;

	/// Unserialize the dynamics object. Including the private variables
	virtual void unserialize ( HierarchicUnserializer & unserializer);
	virtual void setParameter (plint whichParameter, T value);
	virtual T getParameter (plint whichParameter) const;

	virtual CahnHilliardDynamics<T,Descriptor>* clone() const;
	virtual int getId() const;
	virtual void collide(Cell<T,Descriptor>& cell, BlockStatistics& statistics );
	virtual T computeEquilibrium(Cell<T,Descriptor> const& cell, plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j, T jSqr, T thetaBar=T()) const;

private:
	static int id;
	T Gamma;
	T qCoef;
};

/////////////////////////////////////////
template<typename T, template<typename U> class Descriptor>
class ShaoBGKdynamics : public IsoThermalBulkDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    ShaoBGKdynamics(T omega_);

    /// Clone the object on its dynamic type.
    virtual ShaoBGKdynamics<T,Descriptor>* clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

/* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics& statistics_);

    /// Implementation of the collision step, with imposed macroscopic variables
    /*  virtual void collideExternal(Cell<T,Descriptor>& cell, T rhoBar,
	Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& stat);*/

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(Cell<T,Descriptor> const& cell, plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;

    /// Implementation of velocity computation.
    virtual void computeVelocity ( Cell< T, Descriptor > const & cell,
    		Array< T, Descriptor< T >::d > &u ) const;

    /// Implementation of Density computation.
    virtual T computeDensity ( Cell< T, Descriptor > const & cell ) const;

private:
    static int id;
};

} /* namespace slilab */

#endif /* SHAODYNAMICS_H_ */
