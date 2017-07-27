/*
 * Copyright (C) 2017 IIT - Istituto Italiano di Tecnologia - http://www.iit.it
 * Author: Silvio Traversaro
 * CopyPolicy: Released under the terms of the GNU LGPL v2.0 (or any later version)
 *
 * The development of this software was supported by the FP7 EU project
 * CoDyCo (No. 600716 ICT 2011.2.1 Cognitive Systems and Robotics (b))
 * http://www.codyco.eu
 */

#include <iDynTree/iCub/iKinConversions.h>
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
    ret.Min = ikinlink.getMin();
    ret.Max = ikinlink.getMax();

    return ret;
}

iCub::iKin::iKinLink DHLink2iKinLink(const DHLink & dhLink)
{
    return iCub::iKin::iKinLink(dhLink.A,dhLink.D,
                                dhLink.Alpha,dhLink.Offset,
                                dhLink.Min,dhLink.Max);
}

bool DHChainFromiKinChain(iCub::iKin::iKinChain& ikinChain,
                                DHChain& dhChain)
{
    assert(ikinChain.getN() == ikinChain.getDOF());

    iDynTree::Transform H0, HN;
    dhChain.setNrOfDOFs(ikinChain.getN());

    toiDynTree(ikinChain.getH0(),H0);
    dhChain.setH0(H0);

    for(size_t i=0; i < dhChain.getNrOfDOFs(); i++)
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
    bool ok = DHChainFromiKinChain(ikinChain, chain);

    if( !ok )
    {
        return false;
    }

    ok = CreateModelFromDHChain(chain,output);

    return ok;
}

bool iKinChainFromDHChain(const DHChain & dhChain,
                         iCub::iKin::iKinChain& ikinChain)
{
    ikinChain.clear();

    yarp::sig::Matrix yarpMatBuf;

    toYarp(dhChain.getH0().asHomogeneousTransform(),yarpMatBuf);
    ikinChain.setH0(yarpMatBuf);

    for(size_t i=0; i < dhChain.getNrOfDOFs(); i++)
    {
        // TODO : This is a memory leak, but is currently the only way of generating a iKinChain from
        // code without using the iKinLimb::fromLinksProperties method
        iCub::iKin::iKinLink* p_ikinLink = new iCub::iKin::iKinLink(DHLink2iKinLink(dhChain(i)));
        ikinChain.pushLink(*p_ikinLink);
    }

    toYarp(dhChain.getHN().asHomogeneousTransform(),yarpMatBuf);
    ikinChain.setHN(yarpMatBuf);

    return true;
}

bool iKinChainFromModel(const Model & model,
                        const std::string& baseFrame,
                        const std::string& distalFrame,
                        iCub::iKin::iKinChain & ikinChain)
{
    DHChain chain;
    bool conversionSuccessful = ExtractDHChainFromModel(model, baseFrame, distalFrame, chain);
    conversionSuccessful = conversionSuccessful && iKinChainFromDHChain(chain, ikinChain);
    return conversionSuccessful;
}

}