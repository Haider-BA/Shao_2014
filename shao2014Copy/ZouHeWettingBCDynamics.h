/*
 * ZouHeWettingBCDynamics.h
 *
 *  Created on: Mar 4, 2013
 *      Author: minh
 */

#ifndef ZOUHEWETTINGBCDYNAMICS_H_
#define ZOUHEWETTINGBCDYNAMICS_H_

namespace slilab {

class ZouHeWettingBCDynamics {
public:
	ZouHeWettingBCDynamics(Dynamics<T,Descriptor>* baseDynamics, bool automaticPrepareCollision = true);
	ZouHeWettingBCDynamics(HierarchicUnserializer& unserializer);
    /// Clone the object on its dynamic type.
    virtual ZouHeWettingBCDynamics<T, Descriptor, direction, orientation>* clone() const;
    virtual void serialize(HierarchicSerializer& serializer) const;
    virtual void unserialize(HierarchicUnserializer& unserializer);
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T,Descriptor>& cell) const;
private:
    static int id;
};

} /* namespace slilab */
#endif /* ZOUHEWETTINGBCDYNAMICS_H_ */
