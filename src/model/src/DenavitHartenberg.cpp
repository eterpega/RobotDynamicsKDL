/*
 * Copyright (C) 2016 Fondazione Istituto Italiano di Tecnologia
 * Authors: Silvio Traversaro
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <iDynTree/Model/DenavitHartenberg.h>
#include <iDynTree/Model/ForwardKinematics.h>
#include <iDynTree/Model/Traversal.h>
#include <iDynTree/Model/RevoluteJoint.h>
#include <iDynTree/Model/FixedJoint.h>
#include <iDynTree/Model/LinkState.h>
#include <iDynTree/Model/FreeFloatingState.h>

#include <iDynTree/Core/EigenHelpers.h>


namespace iDynTree
{

void DHChain::setNrOfDOFs(size_t nDofs)
{
    dhParams.resize(nDofs);
}

size_t iDynTree::DHChain::getNrOfDOFs() const
{
    return dhParams.size();
}

void iDynTree::DHChain::setH0(const Transform & _H0)
{
    H0 = _H0;
}

const Transform & DHChain::getH0() const
{
    return H0;
}

void iDynTree::DHChain::setHN(const Transform & _HN)
{
    HN = _HN;
}

const Transform& DHChain::getHN() const
{
    return HN;
}

iDynTree::DHLink& iDynTree::DHChain::operator()(const size_t i)
{
    return dhParams[i];
}

const iDynTree::DHLink& iDynTree::DHChain::operator()(const size_t i) const
{
    return dhParams[i];
}

/**
 * Given two lines, find their closest points (i.e. the points belonging to the common normal)
 *
 * Using the algorithm in http://geomalgorithms.com/a07-_distance.html
 *
 * @param[in]
 * @param[in]

 * @return true if the all went well, false it the lines are parallel.
 */
bool closestPoints(const iDynTree::Axis line_A,
                   const iDynTree::Axis line_B,
                         iDynTree::Position & closest_point_line_A,
                         iDynTree::Position & closest_point_line_B,
                         double tol = 1e-6
                  )
{
    /*
      Using the notation in : http://geomalgorithms.com/a07-_distance.html
      direction_line_A is u
      origin_line_A is P_0
      direction_line_B is v
      origin_line_b is Q_0
      closest_point_line_A is P_C
      closest_point_line_B is Q_C
    */
    /*
    cerr << "origin_line_A: " << origin_line_A << std::endl;
    cerr << "origin_line_B: " << origin_line_B << std::endl;
    cerr << "direction_line_A: " << direction_line_A << std::endl;
    cerr << "direction_line_B: " << direction_line_B << std::endl;
    */
    Eigen::Vector3d origin_line_A = toEigen(line_A.getOrigin());
    Eigen::Vector3d origin_line_B = toEigen(line_B.getOrigin());
    Eigen::Vector3d direction_line_A = toEigen(line_A.getDirection());
    Eigen::Vector3d direction_line_B = toEigen(line_B.getDirection());

    Eigen::Vector3d w0 = origin_line_A-origin_line_B;
    double a = direction_line_A.dot(direction_line_A);
    double b = direction_line_A.dot(direction_line_B);
    double c = direction_line_B.dot(direction_line_B);
    double d = direction_line_A.dot(w0);
    double e = direction_line_B.dot(w0);

    double denominator = a*c-b*b;

    // test with fabs because sometimes the numerical zero is computed as -epsilon
    if( fabs(denominator) < tol )
    {
        return false;
    }

    //Denominator should be nonnegative
    assert(denominator >= 0.0);

    double s_C = (b*e-c*d)/denominator;
    double t_C = (a*e-b*d)/denominator;

    toEigen(closest_point_line_A) = origin_line_A + s_C*direction_line_A;
    toEigen(closest_point_line_B) = origin_line_B + t_C*direction_line_B;

    return true;
}

bool ExtractDHChainFromModel(const Model& model,
                             const std::string baseFrame,
                             const std::string eeFrame,
                                   DHChain& outputChain,
                                   double tolerance = 1e-5)
{
    // We check that the frame exist
    if( model.getFrameIndex(baseFrame) == iDynTree::FRAME_INVALID_INDEX ||
        model.getFrameIndex(eeFrame) == iDynTree::FRAME_INVALID_INDEX )
    {
        std::string err = "Frame " + baseFrame + " or frame " + eeFrame + " not found in model.";
        reportError("","ExtractDHChainFromModel",err.c_str());
        return false;
    }

    LinkIndex linkOfBaseFrameIdx = model.getFrameLink(model.getFrameIndex(baseFrame));
    LinkIndex linkOfEEFrameIdx = model.getFrameLink(model.getFrameIndex(eeFrame));

    // To extract the DH parameters, we need first to express all the
    // axis of the chain in a common frame, to simplify the computations.
    // As a first step, we compute the traversal that has as base link
    // the link at which the base frame is attached
    Traversal traversal;
    model.computeFullTreeTraversal(traversal,linkOfBaseFrameIdx);

    // We then loop from the eeLink to the baseLink, counting
    // the number of revolute joints in the middle, giving an error if
    // there is a non-revolute or non-fixed joint
    size_t nrOfDofs = 0;
    for(LinkIndex visitedLinkIdx = linkOfEEFrameIdx;
        visitedLinkIdx != linkOfBaseFrameIdx;
        visitedLinkIdx = traversal.getParentLinkFromLinkIndex(visitedLinkIdx)->getIndex())
    {
        IJointConstPtr joint = traversal.getParentJointFromLinkIndex(visitedLinkIdx);
        bool isRevoluteJoint =
            ( dynamic_cast<const RevoluteJoint*>(joint) != 0 );
        bool isFixedJoint    =
            ( dynamic_cast<const FixedJoint*>(joint) != 0 );

        if( !isFixedJoint && !isRevoluteJoint )
        {
            std::string err = "Joint  " + model.getJointName(joint->getIndex()) + " is not revolute neither fixed, but the DH converter only support this two joints.";
            reportError("","ExtractDHChainFromModel",err.c_str());
        }

        if( isRevoluteJoint )
        {
            nrOfDofs++;
        }
    }

    // Now we know that the output DH chain will have nrOfDofs dofs
    outputChain.setNrOfDOFs(nrOfDofs);

    // We need to write all the revolution axis in the same frame, to simplify computations.
    // The first step is to compute the link position with respect to the base frame, that
    // we can do with the usual fwd pos kinematics.

    FreeFloatingPos chainPos(model);
    // We use the baseFrame as "world", to easily compute the baseFrame_X_link transform
    // using the usual floating base forward kinematics.
    chainPos.worldBasePos() = model.getFrameTransform(model.getFrameIndex(baseFrame)).inverse();
    chainPos.jointPos().zero();
    LinkPositions baseFrame_X_link(model);
    bool ok = ForwardPositionKinematics(model,traversal,chainPos,baseFrame_X_link);

    if( !ok )
    {
        reportError("","ExtractDHChainFromModel","Error in computing ForwardKinematics.");
    }

    // Let's store all the axis in the baseFrame
    size_t nrOFDHFrames = nrOfDofs+1;
    std::vector<iDynTree::Axis> zAxes(nrOFDHFrames);
    std::vector<iDynTree::Direction> xDirection(nrOFDHFrames);
    std::vector<iDynTree::Direction> yDirection(nrOFDHFrames);

    size_t counter = 0;
    for(LinkIndex visitedLinkIdx = linkOfEEFrameIdx;
        visitedLinkIdx != linkOfBaseFrameIdx;
        visitedLinkIdx = traversal.getParentLinkFromLinkIndex(visitedLinkIdx)->getIndex())
    {
        IJointConstPtr joint = traversal.getParentJointFromLinkIndex(visitedLinkIdx);
        bool isRevoluteJoint =
            ( dynamic_cast<const RevoluteJoint*>(joint) != 0 );

        if( isRevoluteJoint )
        {
            const RevoluteJoint * revJoint = dynamic_cast<const RevoluteJoint*>(joint);
            zAxes[counter] = baseFrame_X_link(visitedLinkIdx)*revJoint->getAxis(visitedLinkIdx);
            counter++;
        }
    }

    // We have the H0 transformation to account for the
    // transform between the base frame and the first frame (0)
    // of the DH representation. Hence, we have to choose the
    // first frame for the DH representation. We already have its origin
    // and the z axis, we have then to choose an arbitrary direction for the
    // x and y axis.
    // The strategy that we use for getting this x and y axes is the following:
    //      * if the z axis of the base frame and the 0 DH frame are
    //        parallel, we can assign the x and y axes of the 0 DH frame
    //        to be the same x and y axes of the base frame.
    //      * if the z axis of the base frame and the 0 DH frame are
    //        not parallel, we can find a rotation around an axis that is
    //        transforming the z axis of the base frame to the 0 frame, and
    //        we apply the same transformation to the x and y axes of the base frame
    //        to obtain the x and y axes fo the DH frame.

    // We write the z axis in the base frame
    Axis z_in_base = Axis(Direction(0,0,1.0),Position::Zero());

    // We check if z_in_base and the joint axes of the first joint are parallel
    bool z_base_and_z_dh_0_are_parallel = z_in_base.isParallel(zAxes[0],tolerance);

    if( z_base_and_z_dh_0_are_parallel )
    {
        xDirection[0] = Direction(1,0,0);
        toEigen(yDirection[0]) = toEigen(zAxes[0].getDirection()).cross(toEigen(xDirection[0]));
        yDirection[0].Normalize();
    }
    else
    {
        iDynTree::Direction rotAxis;
        Eigen::Vector3d nonNormalizedRotAxis = toEigen(Direction(0,0,1.0)).cross(toEigen(zAxes[0].getDirection()));
        toEigen(rotAxis) = nonNormalizedRotAxis;
        rotAxis.Normalize();
        double rot_angle = atan2(toEigen(rotAxis).dot(nonNormalizedRotAxis),toEigen(Direction(0,0,1.0)).dot(toEigen(zAxes[0].getDirection())));
        iDynTree::Rotation base_R_0 = iDynTree::Rotation::RotAxis(rotAxis,rot_angle);
        xDirection[0] = base_R_0*Direction(1,0,0);
        yDirection[0] = base_R_0*Direction(0,1,0);
    }


    assert(xDirection[0].isPerpendicular(zAxes[0].getDirection(),tolerance));
    assert(yDirection[0].isPerpendicular(zAxes[0].getDirection(),tolerance));

    // The origin of each DH frame is assigned when computing the DH parameters
    std::vector<Position> dh_origin(nrOFDHFrames);

    // For the first one, we are free to choose it
    dh_origin[0] = zAxes[0].getOrigin();

    // We actually compute the DH parameters
    // and select the origin, x and y axis of the DH frames
    // (the z axis are already determined by the axes directions)



    return true;
}



}