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
 *  Lattice descriptors for Shao-Shu-Huang-Chew (PRE 2014) -- header file
 */
#ifndef SHAO_LATTICES_HH
#define SHAO_LATTICES_HH

#include "shaoLattice.h"

namespace plb {

namespace descriptors {

// MRT D2Q5 ////////////////////////////////////////////////////////////
template<typename T>
const T MRTD2Q5DescriptorBase<T>::M[BaseDescriptor::q][BaseDescriptor::q] =
{
		{ (T)1, (T)1, (T)1, (T)1, (T)1 },
		{ (T), -(T)1, (T), (T)1, (T) },
		{ (T), (T), -(T)1, (T), (T)1 },
		{ -(T)4, (T)1, (T)1, (T)1, (T)1 },
		{ (T), (T)1, -(T)1, (T)1, -(T)1 }
};

template<typename T>
const T MRTD2Q5DescriptorBase<T>::invM[BaseDescriptor::q][BaseDescriptor::q] =
{
		{ (T)1/(T)5, (T), (T), -(T)1/(T)5, (T) },
		{ (T)1/(T)5, -(T)1/(T)2, (T), (T)1/(T)20, (T)1/(T)4 },
		{ (T)1/(T)5, (T), -(T)1/(T)2, (T)1/(T)20, -(T)1/(T)4 },
		{ (T)1/(T)5, (T)1/(T)2, (T), (T)1/(T)20, (T)1/(T)4 },
		{ (T)1/(T)5, (T), (T)1/(T)2, (T)1/T)20, -(T)1/(T)4 }
};

template<typename T>
const T MRTD2Q5DescriptorBase<T>::S[BaseDescriptor::q] =
{ (T), (T), (T), (T), (T) };

template<typename T>
const int MRTD2Q5DescriptorBase<T>::momentumIndexes[MRTD2Q5DescriptorBase<T>::jIndexes] = {3, 5};

template<typename T>
const int MRTD2Q5DescriptorBase<T>::shearViscIndexes[MRTD2Q5DescriptorBase<T>::shearIndexes] = {7, 8};

template<typename T>
const int MRTD2Q5DescriptorBase<T>::qViscIndexes[MRTD2Q5DescriptorBase<T>::qIndexes] = {4, 6};

template<typename T>
const char MRTD2Q5Descriptor<T>::name[] = "MRTD2Q5";

// EXTERNAL FORCE MRT D2Q5 ////////////////////////////////////////////////////////////

template<typename T>
const char ForcedMRTD2Q5Descriptor<T>::name[] = "ForcedMRTD2Q5";

}  // namespace descriptors

}  // namespace plb

#endif  // SHAO_LATTICES_HH
