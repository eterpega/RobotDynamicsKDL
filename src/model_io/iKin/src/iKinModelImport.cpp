/*
 * Copyright (C) 2016 Fondazione Istituto Italiano di Tecnologia
 * Authors: Silvio Traversaro
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 */

#include <iDynTree/iKinModelImport.h>
#include <iDynTree/Model/DenavitHartenberg.h>

#include <iDynTree/yarp/YARPConversions.h>

#include <iCub/iKin/iKinFwd.h>

#include <cassert>

namespace iDynTree
{

DHLink iKinLink2DHLink(const iCub::iKin::iKinLink & ikinlink)
{
    DHLink ret;

    ret.A = ikinlink.getA();
    ret.D = ikinlink.getD();
    ret.Alpha = ikinlink.getAlpha();
    ret.Offset = ikinlink.getOffset();

    return ret;
}

bool DHChainFromiKinChain(iCub::iKin::iKinChain& ikinChain,
                                DHChain& dhChain)
{
    assert(ikinChain.getN() == ikinChain.getDOF());

    iDynTree::Transform H0, HN;
    dhChain.setNrOfDOFs(ikinChain.getN());

    toiDynTree(ikinChain.getH0(),H0);
    dhChain.setH0(H0);

    for(int i=0; i < dhChain.getNrOfDOFs(); i++)
    {
        dhChain(i) = iKinLink2DHLink(ikinChain(i));
    }

    toiDynTree(ikinChain.getHN(),HN);
    dhChain.setHN(HN);

    return true;
}

bool modelFromiKinChain(iCub::iKin::iKinChain& ikinChain, Model& output)
{
    DHChain chain;
    bool ok = DHChainFromiKinChain(ikinChain,chain);

    if( !ok )
    {
        return false;
    }

    ok = CreateModelFromDHChain(chain,output);

    return ok;
}



}

