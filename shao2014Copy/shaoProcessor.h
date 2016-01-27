/*
 * ShaoProcessor.h
 *
 *  Created on: Dec. 04, 2015
 *     
 */

#ifndef SHAO_PROCESSOR_H_
#define SHAO_PROCESSOR_H_

using namespace plb;

namespace slilab {

/// Shao coupling between CH & NS
template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
class ShaoProcessor2D : public BoxProcessingFunctional2D {
public:
    ShaoProcessor2D( T2 A, T2 K, T2 Phi, T1 grav_, T1 rhoHigh, T1 rhoLow);
    virtual ~ShaoProcessor2D();
    ShaoProcessor2D(ShaoProcessor2D<T1,Descriptor1,T2,Descriptor2> const& rhs);
    ShaoProcessor2D& operator=(ShaoProcessor2D<T1,Descriptor1,T2,Descriptor2> const& rhs);
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> fields );
    // virtual BlockDomain::DomainT appliesTo() const;
    virtual ShaoProcessor2D<T1,Descriptor1,T2,Descriptor2>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    T2 aCoef, kappa, phiCoef;
    T1 grav, rho_h, rho_l;
};

// ______________________________________________________________________________
/*template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
class ShaoProcessor3D : public BoxProcessingFunctional3D_LL<T1,Descriptor1,T2,Descriptor2> {
public:
	ShaoProcessor3D( T2 A, T2 K, T2 Phi, T1 grav_);
    virtual ~ShaoProcessor3D();
    ShaoProcessor3D(ShaoProcessor3D<T1,Descriptor1,T2,Descriptor2> const& rhs);
    ShaoProcessor3D& operator=(ShaoProcessor3D<T1,Descriptor1,T2,Descriptor2> const& rhs);
    virtual void process(Box3D domain, BlockLattice3D<T1,Descriptor1>& lattice1, BlockLattice3D<T2,Descriptor2>& lattice2 );
    virtual ShaoProcessor3D<T1,Descriptor1,T2,Descriptor2>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    T2 aCoef, kappa, phiCoef2;
    T1 grav;
};
*/
}
#endif /* SHAO_PROCESSOR_H_ */
