/*
 * Copyright (C) 2016 Fondazione Istituto Italiano di Tecnologia
 * Authors: Silvio Traversaro
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 */

#ifndef IDYNTREE_IKIN_MODEL_IMPORT_H
#define IDYNTREE_IKIN_MODEL_IMPORT_H

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

bool DHChainFromiKinChain(iCub::iKin::iKinChain & ikinChain,
                           DHChain & out);

/**
 * \ingroup iDynTreeModelIO
 *
 * Load a iDynTree::Model object from a iKin chain.
 *
 *
 * @return true if all went ok, false otherwise.
 */
bool modelFromiKinChain(iCub::iKin::iKinChain & ikinChain,
                        iDynTree::Model & output);



}

#endif