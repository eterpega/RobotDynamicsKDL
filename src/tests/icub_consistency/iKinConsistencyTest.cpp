/*
 * Copyright (C) 2015 Fondazione Istituto Italiano di Tecnologia
 * Authors: Silvio Traversaro
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <iDynTree/Core/TestUtils.h>

#include <iDynTree/iKinModelImport.h>

#include <iDynTree/KinDynComputations.h>

#include <iDynTree/Model/Model.h>
#include <iDynTree/Model/JointState.h>

#include <iDynTree/yarp/YARPConversions.h>

#include <iCub/iKin/iKinFwd.h>

#include <cassert>
#include <cstdio>
#include <cstdlib>

using namespace iDynTree;

void check_iKinImportFromChain(iCub::iKin::iKinChain & chain)
{
    // First unblock all DOF of the chain
    for(int i=0; i < chain.getN(); i++)
    {
        chain.releaseLink(i);
    }


    Model model;
    bool ok = modelFromiKinChain(chain,model);
    ASSERT_IS_TRUE(ok);

    ASSERT_EQUAL_DOUBLE(chain.getDOF(),model.getNrOfDOFs());

    KinDynComputations dynComp;
    ok = dynComp.loadRobotModel(model);
    ASSERT_IS_TRUE(ok);

    // Given how the model is built by modelFromiKinChain, we now that link 0 is the base
    // and the one with the highest index is the ee
    Vector3 grav;
    grav.zero();
    JointPosDoubleArray qj(model);
    JointDOFsDoubleArray dqj(model);
    dqj.zero();

    yarp::sig::Vector qj_yarp(model.getNrOfPosCoords());
    yarp::sig::Vector dqj_yarp(model.getNrOfDOFs());

    // Check for some various positions
    for(int i=0; i < 10; i++ )
    {
        getRandomVector(qj,-3.14,3.14);
        toYarp(qj,qj_yarp);
        toYarp(dqj,dqj_yarp);

        // First set the position in iKin to check any limit
        qj_yarp = chain.setAng(qj_yarp);
        toiDynTree(qj_yarp,qj);

        dynComp.setRobotState(qj,dqj,grav);

        // Get b_H_ee for both iKin and iDynTree
        LinkIndex baseIndex = 0;
        LinkIndex eeIndex   = model.getNrOfLinks()-1;
        Transform b_H_ee = dynComp.getRelativeTransform(baseIndex,eeIndex);
        yarp::sig::Matrix b_H_ee_ikin_yarp = chain.getH();
        Transform b_H_ee_ikin;
        toiDynTree(b_H_ee_ikin_yarp,b_H_ee_ikin);
        
        ASSERT_EQUAL_TRANSFORM(b_H_ee,b_H_ee_ikin);
    }
}

void check_iKinImport()
{
    iCub::iKin::iCubLeg leg;
    check_iKinImportFromChain(leg);
    iCub::iKin::iCubArm arm;
    check_iKinImportFromChain(arm);
}


int main()
{
    check_iKinImport();
    return EXIT_SUCCESS;
}