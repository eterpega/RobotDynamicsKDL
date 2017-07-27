/*
 * Copyright (C) 2017 Fondazione Istituto Italiano di Tecnologia
 * Authors: Silvio Traversaro
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#ifndef IDYNTREE_IKIN_CONVERSIONS_H
#define IDYNTREE_IKIN_CONVERSIONS_H

#include <string>

namespace iCub {
namespace iKin {

class iKinChain;

}
}

namespace iDynTree

{

class Model;
class DHChain;

/**
 *
 * Load a iDynTree::DHChain object from a iCub::iKin::iKinChain .
 *
 * @return true if all went ok, false otherwise.
 */
bool DHChainFromiKinChain(iCub::iKin::iKinChain & ikinChain,
                          DHChain & out);

/**
 *
 * Load a iDynTree::Model object from a iCub::iKin::iKinChain .
 *
 * @return true if all went ok, false otherwise.
 */
bool modelFromiKinChain(iCub::iKin::iKinChain & ikinChain,
                        Model & output);

/**
 *
 * Extract an iCub::iKin::iKinChain from an iDynTree::Model .
 *
 * @return true if all went ok, false otherwise.
 */
bool iKinChainFromModel(const Model & model,
                        const std::string& baseFrame,
                        const std::string& distalFrame,
                        iCub::iKin::iKinChain & ikinChain);

/**
 *
 * Create a iCub::iKin::iKinChain from an iDynTree::DHChain
 */
bool iKinChainFromDHChain(const DHChain & dhChain,
                          iCub::iKin::iKinChain& ikinChain);

}

#endif
