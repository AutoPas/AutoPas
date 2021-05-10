/**
 * @file MDFlexParticel.cpp
 * @author J. Körner
 * @date 04.05.2021
 */

#include "MDFlexParticle.h"

MDFlexParticle::MDFlexParticle() {
	_attributes = reinterpret_cast<Attributes*>(&_r[0]);
}

void MDFlexParticle::setAttributes(char* &attributeData) {
	_attributes = reinterpret_cast<Attributes*>(&attributeData[0]);
}

