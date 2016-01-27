/*
 * ShaoDynamics.cpp
 *
 *  Created on: Dec 04, 2015
 *     
 */


#include "shaoDynamics.h"

namespace slilab {

template<typename T, template<typename U> class Descriptor>
int CahnHilliardDynamics<T,Descriptor>::id =
		meta::registerGeneralDynamics<T,Descriptor,CahnHilliardDynamics<T,Descriptor> >("CahnHilliardDynamics");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
CahnHilliardDynamics<T,Descriptor>::CahnHilliardDynamics (T omega_ )
: AdvectionDiffusionDynamics<T,Descriptor>(omega_)
  { }

template<typename T, template<typename U> class Descriptor>
CahnHilliardDynamics<T,Descriptor>* CahnHilliardDynamics<T,Descriptor>::clone() const {
	return new CahnHilliardDynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int CahnHilliardDynamics<T,Descriptor>::getId() const {
	return id;
}

template<typename T, template<typename U> class Descriptor>
CahnHilliardDynamics<T,Descriptor>::CahnHilliardDynamics(HierarchicUnserializer& unserializer)
:AdvectionDiffusionDynamics<T,Descriptor>((T)1)
 {
	unserialize( unserializer );
 }

template<typename T, template<typename U> class Descriptor>
CahnHilliardDynamics<T,Descriptor>::CahnHilliardDynamics(CahnHilliardDynamics<T,Descriptor> const &rhs)
: AdvectionDiffusionDynamics<T,Descriptor>(rhs.getOmega())
  {
	Gamma = rhs.Gamma;
	qCoef = rhs.qCoef;
  }

template<typename T, template<typename U> class Descriptor>
void CahnHilliardDynamics<T,Descriptor>::serialize ( HierarchicSerializer & serializer ) const
{
	AdvectionDiffusionDynamics<T,Descriptor>::serialize(serializer);
	serializer.addValue(Gamma);
	serializer.addValue(qCoef);
}

template<typename T, template<typename U> class Descriptor>
void CahnHilliardDynamics<T,Descriptor>::unserialize ( HierarchicUnserializer & unserializer)
{
	AdvectionDiffusionDynamics<T,Descriptor>::unserialize(unserializer);
	Gamma = unserializer.readValue<T>();
	qCoef = unserializer.readValue<T>();
}

template<typename T, template<typename U> class Descriptor>
void CahnHilliardDynamics<T,Descriptor>::setParameter (plint whichParameter, T value)
{
	switch(whichParameter) {
	case slilab::Shao::Gamma:
		Gamma = value;
		break;
	case slilab::Shao::qCoef:
		qCoef = value;
		break;
	default:
		AdvectionDiffusionDynamics<T,Descriptor>::setParameter(whichParameter, value);
		break;
	}
}

template<typename T, template<typename U> class Descriptor>
T CahnHilliardDynamics<T,Descriptor>::getParameter (plint whichParameter) const
{
	switch(whichParameter) {
	case slilab::Shao::Gamma:
		return(Gamma);
		break;
	case slilab::Shao::qCoef:
		return( qCoef );
		break;
	default:
		return(AdvectionDiffusionDynamics<T,Descriptor>::getParameter(whichParameter));
		break;
	}
}

// ___________
template<typename T, template<typename U> class Descriptor>
void CahnHilliardDynamics<T,Descriptor>::collide (
		Cell<T,Descriptor>& cell, BlockStatistics& statistics
)
{
	typedef Descriptor<T> D;
	// Handle to external scalars
	enum {
		scalarOffset  = D::ExternalField::densityBeginsAt,
		velocityBeginsAt = D::ExternalField::velocityBeginsAt,
		fBeginsAt = D::ExternalField::fBeginsAt,
	};

	T omega = this->getOmega();
	T one_m_omega = qCoef-omega;
	T one_m_q = (T)1-qCoef;
	T oneOver2q = (T)0.5/qCoef;

	// Compute phi from the distribution
	T phiBar = getRhoBar( cell );
	T phi = Descriptor<T>::fullRho( phiBar );

	// Get the nab(phi) & velocity from external variables
	T mu   = *cell.getExternal(scalarOffset);
	Array<T,Descriptor<T>::d> u;
	Array<T,Descriptor<T>::q> fs;
	u.from_cArray( cell.getExternal(velocityBeginsAt) );
	fs.from_cArray( cell.getExternal(fBeginsAt) );

	Array<T,D::q> feq;
	Array<T,D::q> &f = cell.getRawPopulations();

	const T B = Gamma*mu*D::invCs2;

	T edotu = (T)0;
	for(plint d=0; d < Descriptor<T>::d; ++d)
	{
		edotu += D::c[0][d]*u[d];
	}

	feq[0] = (D::t[0] - (T)1) * B +
			phi +
			D::t[0] * oneOver2q * phi * edotu;

	for (plint iPop=1; iPop < Descriptor<T>::q; ++iPop) {
		T edotu = (T)0;
		for(plint d=0; d < Descriptor<T>::d; ++d)
		{
			edotu += D::c[iPop][d]*u[d];
		}

		feq[iPop] = D::t[iPop] * (B + oneOver2q * phi * edotu);
	}

	for(plint i=0; i<Descriptor< T >::q; i++) {
		f[i] *= one_m_omega;
		f[i] += (feq[i] * omega);
		f[i] += (one_m_q * fs[i]);
	}

	//	if (cell.takesStatistics()) {
	//		gatherStatistics(statistics, rhoBar, uSqr);
	//	}

}

template<typename T, template<typename U> class Descriptor>
T CahnHilliardDynamics<T,Descriptor>::computeEquilibrium ( Cell<T,Descriptor> const& cell,
		plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j, T jSqr, T thetaBar
) const
{
	typedef Descriptor<T> D;
	enum {
		scalarOffset  = D::ExternalField::densityBeginsAt,
		velocityBeginsAt = D::ExternalField::velocityBeginsAt,
	};

	T oneOver2q = (T)0.5/qCoef;
	T phi = Descriptor<T>::fullRho( rhoBar );
	T mu   = *cell.getExternal(scalarOffset);
	Array<T,Descriptor<T>::d> u;
	u.from_cArray( cell.getExternal(velocityBeginsAt) );

	const T B = Gamma*mu*D::invCs2;

	T edotu = (T)0;
	for(plint d=0; d < Descriptor<T>::d; ++d)
	{
		edotu += D::c[iPop][d]*u[d];
	}
	if (iPop==0)
		return( (D::t[0] - (T)1) * B +
			phi +
			D::t[0] * oneOver2q * phi * edotu );
	else
		return( D::t[iPop] * (B + oneOver2q * phi * edotu) );
}
/* *************** Class ShaoBGKdynamics *********************************************** */

template<typename T, template<typename U> class Descriptor>
int ShaoBGKdynamics<T,Descriptor>::id =
		meta::registerOneParamDynamics<T,Descriptor,ShaoBGKdynamics<T,Descriptor> >("ShaoBGK");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
ShaoBGKdynamics<T,Descriptor>::ShaoBGKdynamics(T omega_ )
: IsoThermalBulkDynamics<T,Descriptor>(omega_)
  { }

template<typename T, template<typename U> class Descriptor>
ShaoBGKdynamics<T,Descriptor>* ShaoBGKdynamics<T,Descriptor>::clone() const {
	return new ShaoBGKdynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int ShaoBGKdynamics<T,Descriptor>::getId() const {
	return id;
}

template<typename T, template<typename U> class Descriptor>
void ShaoBGKdynamics<T,Descriptor>::collide (
		Cell<T,Descriptor>& cell,
		BlockStatistics& statistics
)
{
	typedef Descriptor<T> D;
        Array<T,D::d> force, Grad_Density;
	T omega = this->getOmega();
	T OneMinusHalfOmega = 1.-0.5*omega;

	// Compute f from the distribution
	Array< T, D::q > &f = cell.getRawPopulations();

	const T rho = (*cell.getExternal(D::ExternalField::densityBeginsAt));
	Array<T,D::d> u, j;
	T rhoBar, u_gradDensity ;
	momentTemplates<T,Descriptor>::get_rhoBar_j(cell,rhoBar, j);
	const T jSqr = VectorTemplateImpl<T,D::d>::normSqr(j);

	//	pcout<<"before bgk_ma2_collision"<<std::endl;
	//	T uSqr = dynamicsTemplates<T,Descriptor>::bgk_ma2_collision(cell, rho, j, this->getOmega());
	//computing collision instead of calling bgk_ma2_collision
	 T invRho = D::invRho(rhoBar);
	 //this will be used only for statistics
	 T uSqr = jSqr*invRho*invRho;

	 //pcout<<"after bgk_ma2_collision"<<std::endl;
	//	getchar();

	Grad_Density.from_cArray(cell.getExternal(D::ExternalField::GradDensityBeginsAt));
	// External forces
	force.from_cArray(cell.getExternal(D::ExternalField::forceBeginsAt));


	for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
	  {
	    u[iD] = (j[iD] - 0.5*force[iD])/rho;
    
	  }


	for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
	  {
	  
	    u_gradDensity +=u[iD]*Grad_Density[iD];
	  }
   
	T rho_0 = rhoBar + 0.5*u_gradDensity;
	

	for (plint iPop = 0; iPop < D::q; ++iPop)
	  {
	    
	    //computing force term
	   
	    T feq = dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rho, (T)1/rho, j, jSqr);
	   
	    T Gamma_alpha = feq / rho;

	    //to fix the missing term from equilibrium function	    
	    feq += D::t[iPop]*(rho_0 - rho);

	    T forceTerm = (T)0;
	    for (plint iD = 0; iD < D::d; ++iD){
	      T force_firstTerm = (D::c[iPop][iD]-u[iD]) * D::invCs2 * OneMinusHalfOmega;
	      T force_secondTerm = Grad_Density[iD] * D::cs2*(Gamma_alpha - D::t[iPop]);
	      force_secondTerm -= force[iD]*Gamma_alpha;
	      forceTerm +=  force_firstTerm * force_secondTerm;
	      
	    }

	    f[iPop] *= (T)1-omega;
	    f[iPop] += omega*feq;
	    f[iPop] += forceTerm;
	   
	
	  }

	

	if (cell.takesStatistics()) {
	  gatherStatistics(statistics, rhoBar, uSqr);
	}
}

  /*template<typename T, template<typename U> class Descriptor>
void ShaoBGKdynamics<T,Descriptor>::collideExternal (
		Cell<T,Descriptor>& cell, T rhoBar,
		Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& stat
)
{
	typedef Descriptor<T> D;
	T omega = this->getOmega();

	Array<T,D::d> u;
	T rho = D::fullRho( rhoBar );
	for (plint iD = 0; iD < D::d; ++iD)
	{
		u[iD] = j[iD] / rho;
	}

	T uSqr = dynamicsTemplates<T,Descriptor>::bgk_ma2_collision(cell, rhoBar, j, this->getOmega());

	const T A = omega*D::invCs2*(*cell.getExternal(D::ExternalField::densityBeginsAt));

	Array< T, D::q > &f = cell.getRawPopulations();
	f[0] += A*(D::t[0]-(T)1);
	for (plint iPop = 1; iPop < Descriptor<T>::q; ++iPop)
	{
		f[iPop] += D::t[iPop]*A;
	}

	// External forces
	externalForceTemplates<T,Descriptor>::addGuoForce(cell, u, this->getOmega(), (T)1);

	}*/

template<typename T, template<typename U> class Descriptor>
T ShaoBGKdynamics<T,Descriptor>::computeEquilibrium( Cell< T, Descriptor > const & cell,
		plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
		T jSqr, T thetaBar
) const
{
  
  
	typedef Descriptor<T> D;



	//compute rho_0

	Array<T,D::d> force,Grad_density,u;
	T u_gradDensity=0.0;

	//compute velocity 
	const T rho = (*cell.getExternal(D::ExternalField::densityBeginsAt));
	force.from_cArray(cell.getExternal(D::ExternalField::forceBeginsAt));
	for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
	  {
	    u[iD] = (j[iD] - 0.5*force[iD])/rho;
	   
	  }
	//end:compute velocity 
	
	//compute density
	Grad_density.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::GradDensityBeginsAt));
	
	for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
	  {
	  
	    u_gradDensity +=u[iD]*Grad_density[iD];
	   
	    
	  }
	
	T rho_0= rhoBar + 0.5*u_gradDensity;
	
	//end computeDensity (rho_0)


	T feq = dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rho, (T)1/rho, j, jSqr);

	return feq + D::t[iPop]*(rho_0 - rho);




}

template<typename T, template<typename U> class Descriptor>
void ShaoBGKdynamics<T,Descriptor>::computeVelocity ( Cell< T, Descriptor > const & cell,
		Array< T, Descriptor< T >::d > &u
) const
{
  
  typedef Descriptor<T> D;
  Array<T,D::d> force, j;
  momentTemplates<T,Descriptor>::get_j(cell, j);
  force.from_cArray(cell.getExternal(D::ExternalField::forceBeginsAt));
	
  const T rho = (*cell.getExternal(D::ExternalField::densityBeginsAt));

  for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
    {
      u[iD] = (j[iD] - 0.5*force[iD])/rho;
     
    }
 
}

template<typename T, template<typename U> class Descriptor>
T ShaoBGKdynamics<T,Descriptor>::computeDensity ( Cell< T, Descriptor > const & cell
) const
{
  typedef Descriptor<T> D;
  T rhoBar, rho_0, density, u_gradDensity=0.0;
  Array<T,D::d> force, u, Grad_density, j;

  momentTemplates<T,Descriptor>::get_rhoBar_j(cell,rhoBar,j);
  
  Grad_density.from_cArray(cell.getExternal(D::ExternalField::GradDensityBeginsAt));

  force.from_cArray(cell.getExternal(D::ExternalField::forceBeginsAt));
 
  const T rho = (*cell.getExternal(D::ExternalField::densityBeginsAt));
   
  for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
    {
      u[iD] = (j[iD] - 0.5*force[iD])/rho;
    
    }


  for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
    {
	  
      u_gradDensity +=u[iD]*Grad_density[iD];
    }
 
  rho_0= rhoBar + 0.5*u_gradDensity;
  /* if(rho_0<500 || rho_0>501){
    pcout<<" u "<<u[0]<<" ,"<<u[1]<<std::endl;
    pcout<<" Grad_Density "<<Grad_density[0]<<" ,"<<Grad_density[1]<<std::endl;
    pcout<<"rhoBar "<<rhoBar<<" u_gradDensity "<<u_gradDensity<<std::endl; 
    pcout<<"HEllo computeDensity rho_0"<<rho_0<<std::endl;}
  */
  return rho_0;
 
  
}

} /* namespace slilab */
