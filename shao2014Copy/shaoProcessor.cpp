/*
 * ShaoProcessor.cpp
 *
 *  Created on: Dec 04, 2015
 *      
 */

#include "shaoProcessor.h"

using namespace plb;

namespace slilab {

/* *************** ShaoProcessor2D ***************** */
template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
ShaoProcessor2D<T1,Descriptor1,T2,Descriptor2>::ShaoProcessor2D
(T2 A, T2 K, T2 phi, T1 grav_, T1 rhoH_, T1 rhoL_) :
aCoef(A),
kappa(K),
phiCoef(phi),
grav(grav_),
rho_h(rhoH_),
rho_l(rhoL_)
{
}

template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
ShaoProcessor2D<T1,Descriptor1,T2,Descriptor2>::~ShaoProcessor2D() {
}

template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
ShaoProcessor2D<T1,Descriptor1,T2,Descriptor2>::ShaoProcessor2D (
		ShaoProcessor2D<T1,Descriptor1,T2,Descriptor2> const& rhs
)
{
	aCoef = rhs.aCoef;
	kappa = rhs.kappa;
	phiCoef = rhs.phiCoef;
	grav = rhs.grav;
	rho_h = rhs.rho_h;
	rho_l = rhs.rho_l;
}

template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
ShaoProcessor2D<T1,Descriptor1,T2,Descriptor2>&
ShaoProcessor2D<T1,Descriptor1,T2,Descriptor2>::operator= (
		ShaoProcessor2D<T1,Descriptor1,T2,Descriptor2> const& rhs )
{
	aCoef = rhs.aCoef;
	kappa = rhs.kappa;
	phiCoef = rhs.phiCoef;
	grav = rhs.grav;
	rho_h = rhs.rho_h;
	rho_l = rhs.rho_l;
	return *this;
}

template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
ShaoProcessor2D<T1,Descriptor1,T2,Descriptor2>*
ShaoProcessor2D<T1,Descriptor1,T2,Descriptor2>::clone() const
{
	return new ShaoProcessor2D<T1,Descriptor1,T2,Descriptor2>(*this);
}

  /*template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
BlockDomain::DomainT ShaoProcessor2D<T1,Descriptor1,T2,Descriptor2>::appliesTo() const
{
	return BlockDomain::bulk;
	}*/
template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
void ShaoProcessor2D<T1,Descriptor1,T2,Descriptor2>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
  //	modified[0] = modif::staticVariables;
  //	modified[1] = modif::staticVariables;
  modified[0] = modif::allVariables;
  modified[1] = modif::allVariables;

}

template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
void ShaoProcessor2D<T1,Descriptor1,T2,Descriptor2>::processGenericBlocks
(Box2D domain, std::vector<AtomicBlock2D*> blocks)
{
	// Short-hand notation for the lattice2 descriptor
	typedef Descriptor2<T2> D;
	// Handle to external scalars
	enum {
		densityOffset  = D::ExternalField::densityBeginsAt,
		velocityOffset = D::ExternalField::velocityBeginsAt,
		fOffset = D::ExternalField::fBeginsAt,
	};

	//	plint nx = domain.getNx() + 2;  // Include a one-cell boundary
	//	plint ny = domain.getNy() + 2;  // Include a one-cell boundary

	BlockLattice2D<T1,Descriptor1> *lattice1 = static_cast<BlockLattice2D<T1,Descriptor1>*>(blocks[0]);
	BlockLattice2D<T2,Descriptor2> *lattice2 = static_cast<BlockLattice2D<T2,Descriptor2>*>(blocks[1]);
	ScalarField2D<T2>* phiField = dynamic_cast<ScalarField2D<T2>*>(blocks[2]);
	//	phiField.setLocation( lattice2.getLocation() );

	
	// Get the offset of lattice2
	//	Dot2D Offset = lattice2.getLocation();
	Dot2D offset = plb::computeRelativeDisplacement(*blocks[1], *blocks[2]);
	for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
	  for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            phiField->get(iX+offset.x,iY+offset.y)
	      = lattice2->get(iX,iY).computeDensity();
	  }
	}

	// Compute the derivative of density, and store they in the external velocity field.
	for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
		for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
		  
			T2 phi = phiField->get(iX+offset.x,iY+offset.y);
			
			T2 phi3_iMinusOne = phiField->get(iX-1+offset.x,iY+offset.y)*phiField->get(iX-1+offset.x,iY+offset.y)*phiField->get(iX-1+offset.x,iY+offset.y);
			
			T2 phi3_iPlusOne = phiField->get(iX+1+offset.x,iY+offset.y)*phiField->get(iX+1+offset.x,iY+offset.y)*phiField->get(iX+1+offset.x,iY+offset.y);
			T2 phi3_jMinusOne = phiField->get(iX+offset.x,iY-1+offset.y)*phiField->get(iX+offset.x,iY-1+offset.y)*phiField->get(iX+offset.x,iY-1+offset.y);
			T2 phi3_jPlusOne = phiField->get(iX+offset.x,iY+1+offset.y)*phiField->get(iX+offset.x,iY+1+offset.y)*phiField->get(iX+offset.x,iY+1+offset.y);
			T2 phiCoef2 = phiCoef*phiCoef;
			// Computation and storage of the final momentum, including tho momentum
			//   difference due to interaction potential and the external force.
			Cell<T1,Descriptor1>& cell1 = lattice1->get(iX,iY);
			Cell<T2,Descriptor2>& cell2 = lattice2->get(iX,iY);

			// Laplace (phi(iX,iY))
			T2 lapPhi = (T2)-4*phi+
					phiField->get(iX-1+offset.x,iY+offset.y)+
					phiField->get(iX+1+offset.x,iY+offset.y)+
					phiField->get(iX+offset.x,iY-1+offset.y)+
					phiField->get(iX+offset.x,iY+1+offset.y);
			T2 lapPhi_iMinusOne =  (T2)-4*phiField->get(iX-1+offset.x,iY+offset.y)+
					phiField->get(iX-2+offset.x,iY+offset.y)+
					phiField->get(iX+offset.x,iY+offset.y)+
					phiField->get(iX-1+offset.x,iY-1+offset.y)+
					phiField->get(iX-1+offset.x,iY+1+offset.y);
			T2 lapPhi_iPlusOne = (T2)-4*phiField->get(iX+1+offset.x,iY+offset.y)+
					phiField->get(iX+offset.x,iY+offset.y)+
					phiField->get(iX+2+offset.x,iY+offset.y)+
					phiField->get(iX+1+offset.x,iY-1+offset.y)+
					phiField->get(iX+1+offset.x,iY+1+offset.y);
			T2 lapPhi_jMinusOne =  (T2)-4*phiField->get(iX+offset.x,iY-1+offset.y)+
					phiField->get(iX-1+offset.x,iY-1+offset.y)+
					phiField->get(iX+1+offset.x,iY-1+offset.y)+
					phiField->get(iX+offset.x,iY-2+offset.y)+
					phiField->get(iX+offset.x,iY+offset.y);
			T2 lapPhi_jPlusOne = (T2)-4*phiField->get(iX+offset.x,iY+1+offset.y)+
					phiField->get(iX-1+offset.x,iY+1+offset.y)+
					phiField->get(iX+1+offset.x,iY+1+offset.y)+
					phiField->get(iX+offset.x,iY+offset.y)+
					phiField->get(iX+offset.x,iY+2+offset.y);

			// Chemical potential
			T2 mu = aCoef*( (T2)4*phi*phi*phi - (T2)4*phiCoef2*phi ) - kappa * lapPhi;
			T2 mu_iMinusOne = aCoef*( (T2)4*phi3_iMinusOne - (T2)4*phiCoef2*phiField->get(iX-1+offset.x,iY+offset.y)) 
			  - kappa * lapPhi_iMinusOne;
			T2 mu_iPlusOne = aCoef*( (T2)4*phi3_iPlusOne - (T2)4*phiCoef2*phiField->get(iX+1+offset.x,iY+offset.y)) 
			  - kappa * lapPhi_iPlusOne;
			T2 mu_jMinusOne = aCoef*( (T2)4*phi3_jMinusOne - (T2)4*phiCoef2*phiField->get(iX+offset.x,iY-1+offset.y)) 
			  - kappa * lapPhi_jMinusOne;
			T2 mu_jPlusOne = aCoef*( (T2)4*phi3_jPlusOne - (T2)4*phiCoef2*phiField->get(iX+offset.x,iY+1+offset.y)) 
			  - kappa * lapPhi_jPlusOne;

			T2 Grad_mu_x = fd::ctl_diff(mu_iPlusOne, mu_iMinusOne);
			T2 Grad_mu_y = fd::ctl_diff(mu_jPlusOne, mu_jMinusOne);

			T2 Grad_phi_x = fd::ctl_diff(phiField->get(iX+1+offset.x,iY+offset.y), phiField->get(iX-1+offset.x,iY+offset.y));
			T2 Grad_phi_y = fd::ctl_diff(phiField->get(iX+offset.x,iY+1+offset.y), phiField->get(iX+offset.x,iY-1+offset.y));

			T1 rho = rho_l + ((phi + phiCoef)/(2 * phiCoef))*(rho_h - rho_l);
			/*	if(rho>1010){
			pcout<<"in processor rho is "<<rho<<" rho_l "<<rho_l<<" rho_h " <<rho_h << " phiCoef " <<phiCoef<<" phi "<< phi<<std::endl;
			getchar();}*/
			T1 Grad_rho_x = rho_l + ((Grad_phi_x + phiCoef)/(2 * phiCoef))*(rho_h - rho_l);
			T1 Grad_rho_y = rho_l + ((Grad_phi_y + phiCoef)/(2 * phiCoef))*(rho_h - rho_l);

			*cell1.getExternal(Descriptor1<T1>::ExternalField::GradDensityBeginsAt) = Grad_rho_x;
			*cell1.getExternal(Descriptor1<T1>::ExternalField::GradDensityBeginsAt+1) = Grad_rho_y;
			
			// Saving f
			for (plint iPop=0; iPop < Descriptor2<T2>::q; ++iPop) {
				*cell2.getExternal(fOffset+iPop) = lattice2->get(iX+D::c[iPop][0],iY+D::c[iPop][1])[iPop];
			}

			

			// Saving \mu to lattice2
			*cell2.getExternal(densityOffset) = mu;

			// Saving \rho (local density) to lattice1
			*cell1.getExternal(Descriptor1<T1>::ExternalField::densityBeginsAt) = rho;
			  

			//fix it \nab(\mu)\phi	// Saving \nab(\phi)\mu to lattice1
			*cell1.getExternal(Descriptor1<T1>::ExternalField::forceBeginsAt)   =
					phi*Grad_mu_x;
			*cell1.getExternal(Descriptor1<T1>::ExternalField::forceBeginsAt+1) =
					phi*Grad_mu_y+
					((phi>(T1)0)?grav:(T1)0);

		
			// velocity
			Array<T1,Descriptor1<T1>::d> u;
			cell1.computeVelocity(u);
			*cell2.getExternal(velocityOffset)   = u[0];
			*cell2.getExternal(velocityOffset+1) = u[1];
			

		}
	}
	//	getchar();
}

/*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * */
/* *************** ShaoProcessor3D ***************** */
/*template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
ShaoProcessor3D<T1,Descriptor1,T2,Descriptor2>::ShaoProcessor3D
(T2 A, T2 K, T2 phi, T1 grav_) :
aCoef(A),
kappa(K),
phiCoef2(phi*phi),
grav(grav_)
{
}

template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
ShaoProcessor3D<T1,Descriptor1,T2,Descriptor2>::~ShaoProcessor3D() {
}

template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
ShaoProcessor3D<T1,Descriptor1,T2,Descriptor2>::ShaoProcessor3D (
		ShaoProcessor3D<T1,Descriptor1,T2,Descriptor2> const& rhs
)
{
	aCoef = rhs.aCoef;
	kappa = rhs.kappa;
	phiCoef2 = rhs.phiCoef2;
	grav = rhs.grav;
}

template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
ShaoProcessor3D<T1,Descriptor1,T2,Descriptor2>&
ShaoProcessor3D<T1,Descriptor1,T2,Descriptor2>::operator= (
		ShaoProcessor3D<T1,Descriptor1,T2,Descriptor2> const& rhs )
{
	aCoef = rhs.aCoef;
	kappa = rhs.kappa;
	phiCoef2 = rhs.phiCoef2;
	grav = rhs.grav;
	return *this;
}

template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
ShaoProcessor3D<T1,Descriptor1,T2,Descriptor2>*
ShaoProcessor3D<T1,Descriptor1,T2,Descriptor2>::clone() const
{
	return new ShaoProcessor3D<T1,Descriptor1,T2,Descriptor2>(*this);
}

template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
void ShaoProcessor3D<T1,Descriptor1,T2,Descriptor2>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
	modified[0] = modif::staticVariables;
	modified[1] = modif::staticVariables;
}

template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
void ShaoProcessor3D<T1,Descriptor1,T2,Descriptor2>::process
(
		Box3D domain, BlockLattice3D<T1,Descriptor1>& lattice1, BlockLattice3D<T2,Descriptor2>& lattice2
)
{
	// Short-hand notation for the lattice2 descriptor
	typedef Descriptor2<T2> D;
	// Handle to external scalars
	enum {
		densityOffset  = D::ExternalField::densityBeginsAt,
		velocityOffset = D::ExternalField::velocityBeginsAt,
		fOffset = D::ExternalField::fBeginsAt,
	};

	plint nx = domain.getNx() + 2;  // Include a one-cell boundary
	plint ny = domain.getNy() + 2;  // Include a one-cell boundary
	plint nz = domain.getNz() + 2;  // Include a one-cell boundary
	ScalarField3D<T2> phiField(nx,ny,nz);
	phiField.setLocation( lattice2.getLocation() );

	Dot3D offset = computeRelativeDisplacement(lattice2, phiField);
	//	pcout << "offset="<<offset.x<<","<<offset.y<<","<<offset.z<<std::endl;

	// This function compute the phase filed on both "sub-domain" and the envelope
	computeDensity(lattice2,phiField);

	// Compute the derivative of density, and store they in the external velocity field.
	for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
		for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
			for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {

				T2 phi = phiField.get(iX,iY,iZ);

				// Computation and storage of the final momentum, including tho momentum
				//   difference due to interaction potential and the external force.
				Cell<T1,Descriptor1>& cell1 = lattice1.get(iX,iY,iZ);
				Cell<T2,Descriptor2>& cell2 = lattice2.get(iX,iY,iZ);

				// Laplace (phi)
				T2 lapPhi = (T2)-6*phi+
						phiField.get(iX-1,iY,iZ)+
						phiField.get(iX+1,iY,iZ)+
						phiField.get(iX,iY-1,iZ)+
						phiField.get(iX,iY+1,iZ)+
						phiField.get(iX,iY,iZ-1)+
						phiField.get(iX,iY,iZ+1);

				// Saving f
				for (plint iPop=0; iPop < Descriptor2<T2>::q; ++iPop) {
					*cell2.getExternal(fOffset+iPop) =
							lattice2.get(iX+D::c[iPop][0],iY+D::c[iPop][1],iZ+D::c[iPop][2])[iPop];
				}

				// Chemical potention
				T2 mu = aCoef*( (T2)4*phi*phi*phi - (T2)4*phiCoef2*phi ) - kappa * lapPhi;

				// Saving \mu to lattice2
				*cell2.getExternal(densityOffset) = mu;

				// Saving \phi\mu to lattice1
				*cell1.getExternal(Descriptor1<T1>::ExternalField::densityBeginsAt) = mu*phi;

				// Saving \nab(\phi)\mu to lattice1
				*cell1.getExternal(Descriptor1<T1>::ExternalField::forceBeginsAt)   = mu*fd::ctl_diff(phiField.get(iX+1,iY,iZ), phiField.get(iX-1,iY,iZ));
				*cell1.getExternal(Descriptor1<T1>::ExternalField::forceBeginsAt+1) =
						mu*fd::ctl_diff(phiField.get(iX,iY+1,iZ), phiField.get(iX,iY-1,iZ))+
						((phi<(T1)0)?grav:(T1)0);
				*cell1.getExternal(Descriptor1<T1>::ExternalField::forceBeginsAt+2) = mu*fd::ctl_diff(phiField.get(iX,iY,iZ+1), phiField.get(iX,iY,iZ-1));

				// velocity
				Array<T1,Descriptor1<T1>::d> u;
				cell1.computeVelocity(u);
				u.to_cArray(cell2.getExternal(velocityOffset));
			}
		}
	}
}
*/
} // namespace
