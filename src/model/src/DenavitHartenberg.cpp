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

bool DHChain::fromModel(const Model& model,
                        const std::string baseFrame,
                        const std::string eeFrame)
{
    return ExtractDHChainFromModel(model,baseFrame,eeFrame,*this);
}

bool DHChain::toModel(Model& outputModel)
{
    return CreateModelFromDHChain(*this,outputModel);
}



Transform TransformFromDHCraig1989(double a,double alpha,double d,double theta)
{
    // returns Modified Denavit-Hartenberg parameters (According to Craig)
    double ct,st,ca,sa;
    ct = cos(theta);
    st = sin(theta);
    sa = sin(alpha);
    ca = cos(alpha);
    return Transform(Rotation(ct,       -st,     0,
                              st*ca,  ct*ca,   -sa,
                              st*sa,  ct*sa,    ca),
                     Position(a, -sa*d, ca*d));
}

Transform TransformFromDH(double a,double alpha,double d,double theta)
// returns Denavit-Hartenberg parameters (Non-Modified DH)
{
        double ct,st,ca,sa;
        ct = cos(theta);
        st = sin(theta);
        sa = sin(alpha);
        ca = cos(alpha);
        return Transform(Rotation(ct,    -st*ca,   st*sa,
                                  st,     ct*ca,  -ct*sa,
                                   0,        sa,      ca   ),
                         Position(a*ct,      a*st,  d));
}

/**
 * Check if two axes are incident (i.e. they share at least a point).
 */
bool checkIfAxesAreIncident(const iDynTree::Axis line_A,
                            const iDynTree::Axis link_B,
                                  double tol=1e-6)
{
//    TODO
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

/**
 * Given two lines, find the dh parameters that describe  the transformation between the two axis
 *                  and return the origin and the x,y axis of the DH reference frame
 *                  (run step 3,4,5,7 of section 3.2.3 of
 *                       http://www.cs.duke.edu/brd/Teaching/Bio/asmb/current/Papers/chap3-forward-kinematics.pdf )
 *                  The input lines and the output origin are expressed in the same frame (an inertial frame)
 *
 * The _hint variable are necessary when the DH convention allow freedom of choice in choosing dh_direction_axis_x_n or
 *   dh_origin_n, to preserve injectivity.
 *
 * @param[in] zAxis_i_minus_1 axis of the i-1 joint, given by the chain structure.
 * @param[in] xAxis_i_minus_1 x axis of the i-1 frame, given by a previous call to calculateDH.
 * @param[in] origin_i_minus_1 origin of the i-1 frame, given by a previous call to calculateDH.
 * @param[in] zAxis_i axis of the i joint, given by the chain structure.
 * @param[out] origin_i origin of the i frame, computed by this function.
 * @param[out] xAxis_i x axis of the i-th frame, output of this function.
 * @param[out] yAxis_i y axis of the i-th frame, output of this function.
 * @param[out] dhParams dh params of i-th frame, output of this function.
 */
/*
bool calculateDH(const iDynTree::Axis zAxis_i_minus_1,
                 const iDynTree::Axis xAxis_i_minus_1,
                 const iDynTree::Position origin_i_minus_1,
                 const iDynTree::Axis zAxis_i,
                 const iDynTree::Direction xAxis_n_direction_hint,
                       iDynTree::Position & origin_i,
                       iDynTree::Axis & xAxis_i,
                       iDynTree::Axis & yAxis_i,
                       DHLink & dhParams,
                       double tol = 1e-6,
                       int verbose = 0)
{

    // STEP 3 : Locate O_i (the origin of i-th frame)
    //          Locate the origin O_i where the common normal to z_i and z_{i-1} intersects
    //          z_i . If z_i intersects z_{i-1}, locate O_i at this intersection. If z_i and
    //          z_{i-1} are parallel, locate O_i to any convenenient location along z_i .
    bool zAxes_not_parallel;
    iDynTree::Position point_on_zAxis_i_minus_1_closest_to_zAxis_i;

    // Compute the closest point of the axis z_i-1 and z_i
    // During this computation we get as a byproduct if the two
    // axis are parallel. If the axis are not parallel, then the
    // origin of the frame i is given by the point on z_i closest
    // to z_i_minus_1
    zAxes_not_parallel = closestPoints(zAxis_i_minus_1,
                                 zAxis_i,
                                 point_on_zAxis_i_minus_1_closest_to_zAxis_i,
                                 origin_i,
                                 tol);

    if( !zAxes_not_parallel )
    {
        //if parallel, the origin is not specified and we resort to the original one of the chain
        origin_i = zAxis_i.getOrigin();
    }

    // At this point, we can already set the origin of the the axis of i-th frame

    // \todo check that origin_i actually lies on zAxis_i

    xAxis_i.setOrigin(origin_i);
    yAxis_i.setOrigin(origin_i);

    //
    // We check if z_i and z_i-1 are incident
    //
    bool zAxes_incident = checkIfAxesAreIncident(zAxis_i_minus_1,zAxis_i);

    bool zAxes_coincident = checkIfAxesAreCoincident(zAxis_i_minus_1,zAxis_i);


    //STEP 4 : Establish the direction of x axis of the i-th frame
    //         Establish x_i along the common normal between z_{i-1} and z_i throught O_i,
    //         or in the direction normal to the z_{i-1} - z_{i} plane if z_{i-1} and z_i intersect.

    // While the book distingush just two cases for STEP 4 , we need to actually distinguish three
    // different cases :
    // * z_{i-1} and z_{i} are not incident : in this case the direction of the common normal is always defined,
    //   and in particular we compute it as the difference between the two origins, projected on the spaces
    //   orthogonal to z_i.
    // * z_{i-1} and z_{i} are incident, but not coincident: then the two axis define a plane with a unique axis
    // * z_{i-1} and z_{i} are coincident (i.e. incident and parallel): in this case there are infinite normals,
    //   and so we choose x_i using the hint provided (tipically the x_i used in the previous representation).
    if( !zAxes_incident )
    {
        // not incident

        // If the axis are not incident, the x axis is the common normal of the two axis
        Eigen::Vector3d origin_diff = toEigen(origin_i-origin_i_minus_1);

        // To actually get the normal, we remove from the origin_diff the projection on z_i_minus_1
        Eigen::Vector3d x_i_candidate_eig = origin_diff-origin_diff.dot(toEigen(zAxis_i_minus_1.getDirection()))*toEigen(zAxis_i_minus_1.getDirection());

        x_i_candidate_eig.normalize();

        iDynTree::Direction x_i_candidate;
        iDynTree::toEigen(x_i_candidate) = x_i_candidate_eig;

        xAxis_i.setDirection(x_i_candidate);
    }
    else
    {
        if( !zAxes_coincident )
        {
            // incident but not coincident
            Eigen::Vector3d x_i_candidate_eig  = toEigen(zAxis_i_minus_1.getDirection()).cross(toEigen(zAxis_i.getDirection()));
            x_i_candidate_eig.normalize();

            // Check direction, use hint
            //The positive direction of the axis_x_n is arbitrary, however for dealing with limit case (link where
           // only alpha is different from zero)
           double dh_direction_axis_x_n_sign = x_i_candidate_eig.dot(toEigen(xAxis_n_direction_hint) );
           x_i_candidate_eig = dh_direction_axis_x_n_sign >= 0 ? x_i_candidate_eig : -x_i_candidate_eig;

           iDynTree::Direction x_i_candidate;
           iDynTree::toEigen(x_i_candidate) = x_i_candidate_eig;

           xAxis_i.setDirection(x_i_candidate);

           std::cerr << "ExtractDHChainFromModel : DH representation is not unique, made an arbitrary choice for x_i sign\n";
        }
        else
        {
            // coincident
            // if the two axis are coincident, the direction of dh_direction_axis_x_n is totally arbitrary
            // as long as it is perpendicular. We will take then the hint provided by : direction_axis_x_n_hint
            // to get the dh_direction_axis_x_n, we will project direction_axis_x_n_hint onto the plane
            // perpendicular to axis_z_n_minus_1 == axis_z_n, and will normalize the resulting vector
            iDynTree::Direction direction_axis_z_i = zAxis_i.getDirection();
            Eigen::Vector3d direction_axis_z_i_eig = toEigen(direction_axis_z_i);
            Eigen::Vector3d x_i_candidate_eig = toEigen(xAxis_n_direction_hint)- toEigen(xAxis_n_direction_hint).dot(direction_axis_z_i_eig)*direction_axis_z_i_eig;

            x_i_candidate_eig.normalize();

            iDynTree::Direction x_i_candidate;
            iDynTree::toEigen(x_i_candidate) = x_i_candidate_eig;

           xAxis_i.setDirection(x_i_candidate);

            std::cerr << "ExtractDHChainFromModel : DH representation is not unique, made an arbitrary choice for x_i direction\n";
        }
    }

    // STEP 4 ADDENDUM : get y_i

    // Once the direction of axis z_n and x_n has been determined,
    //  the direction of axis y_n is simply given by a cross product to
    //  ensure a right handed coordinates system
    iDynTree::Direction y_i_candidate;
    toEigen(y_i_candidate) = toEigen(zAxis_i.getDirection()).cross(toEigen(xAxis_i.getDirection()));

    yAxis_i.setDirection(y_i_candidate);


    ////////////////////////////////////////////////////////////////////
    // STEP 5: Computation of DH parameters
    ////////////////////////////////////////////////////////////////////

    //calculation of a_i
    //distance along x_i from O_i to the intersection of the x_i and z_{i-1} axes
    iDynTree::Position x_i_z_i_minus_1_intersection_A, x_i_z_i_minus_1_intersection_B;

    closestPoints(zAxis_i_minus_1
                  xAxis_i,
                  x_i_z_i_minus_1_intersection_A,
                  x_i_z_i_minus_1_intersection_B,
                  tol);

    //x_i and z_{i-1} should intersecate
    assert(checkIfAxesAreIncident(xAxis_i,zAxis_i_minus_1));

    dhParams.A = -(toEigen(x_i_z_i_minus_1_intersection_B)-toEigen(origin_i)).dot(toEigen(xAxis_i.getDirection());

    //calculation of d_i
    //distance along z_{i-1} from O_{i-1} to the intersection of the x_i and z_{i-1} axes
    dhParams.D = (toEigen(x_i_z_i_minus_1_intersection_A)-toEigen(origin_i_minus_1)).dot(toEigen(zAxis_i_minus_1.getDirection()));

    //calculation of alpha_i
    //angle between z_{i-1} and z_i measured about x_i
    double cos_alpha_i = toEigen(zAxis_i_minus_1.getDirection()).dot(toEigen(zAxis_i.getDirection()));
    assert(((direction_axis_z_n_minus_1*direction_axis_z_n)*dh_direction_axis_x_n).Norm() < tol);
    double sin_alpha_i = (toEigen(zAxis_i_minus_1.getDirection()).cross(toEigen(zAxis_i.getDirection()))).dot(toEigen(xAxis_i.getDirection());
    assert( fabs(cos_alpha_i*cos_alpha_i + sin_alpha_i*sin_alpha_i - 1) < tol);

    dhParams.Alpha = atan2(sin_alpha_i,cos_alpha_i);

    //calculation of theta_i
    //angle between x_{i-1} and x_i measure about z_{i-1}
    double cos_theta_i = toEigen(xAxis_i_minus_1.getDirection()).dot(toEigen(xAxis_i.getDirection()));
    double sin_theta_i = (toEigen(xAxis_i_minus_1.getDirection()).cross(toEigen(xAxis_i.getDirection()))).dot(toEigen(zAxis_i_minus_1.getDirection());
    dhParams.Offset = atan2(sin_theta_i,cos_theta_i);

    return true;
}
*/

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
    std::vector<iDynTree::Axis> xAxes(nrOFDHFrames);
    std::vector<iDynTree::Axis> yAxes(nrOFDHFrames);

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

    Direction xDirection0;
    Direction yDirection0;

    if( z_base_and_z_dh_0_are_parallel )
    {
        xDirection0  = Direction(1,0,0);
        toEigen(yDirection0) = toEigen(zAxes[0].getDirection()).cross(toEigen(xDirection0));
        yDirection0.Normalize();
    }
    else
    {
        iDynTree::Direction rotAxis;
        Eigen::Vector3d nonNormalizedRotAxis = toEigen(Direction(0,0,1.0)).cross(toEigen(zAxes[0].getDirection()));
        toEigen(rotAxis) = nonNormalizedRotAxis;
        rotAxis.Normalize();
        double rot_angle = atan2(toEigen(rotAxis).dot(nonNormalizedRotAxis),toEigen(Direction(0,0,1.0)).dot(toEigen(zAxes[0].getDirection())));
        iDynTree::Rotation base_R_0 = iDynTree::Rotation::RotAxis(rotAxis,rot_angle);
        xDirection0 = base_R_0*Direction(1,0,0);
        yDirection0 = base_R_0*Direction(0,1,0);
    }

    assert(xDirection0.isPerpendicular(zAxes[0].getDirection(),tolerance));
    assert(yDirection0.isPerpendicular(zAxes[0].getDirection(),tolerance));

    // The origin of each DH frame is assigned when computing the DH parameters
    std::vector<Position> dh_origin(nrOFDHFrames);

    // For the first one, we are free to choose it
    dh_origin[0] = zAxes[0].getOrigin();
    xAxes[0].setOrigin(dh_origin[0]);
    xAxes[0].setDirection(xDirection0);
    yAxes[0].setOrigin(dh_origin[0]);
    yAxes[0].setDirection(yDirection0);

    // We actually compute the DH parameters
    // and select the origin, x and y axis of the DH frames
    // (the z axis are already determined by the axes directions)
    for(size_t i=0; i < nrOfDofs; i++ )
    {
        assert(fabs(axis_z_i[i].Norm()-1) < tol);
        assert(fabs(axis_z_i[i+1].Norm()-1) < tol);
        /*
        calculateDH(zAxes[i],
                    xAxes[i],
                    dh_origin[i],
                    zAxes[i+1],
                    dh_origin[i+1],
                    xAxes[i+1],
                    yAxes[i+1],
                    outputChain(i+1));*/
    }

    // transform matrix from the 0th frame to the base frame
    iDynTree::Transform base_H_0;

    //transform from the end-effector to the N-th frame
    iDynTree::Transform N_H_ee;

    // H0_kdl = KDL::Frame(KDL::Rotation(axis_x_i[0],axis_y_i[0],axis_z_i[0]),origin_O_i[0]);


    return false;
}

std::string inline intToString(const int inInt)
{
    std::stringstream ss;
    ss << inInt;
    return ss.str();
}

bool CreateModelFromDHChain(const DHChain& inputChain,
                                  Model& outputModel)
{
    // First we clear the model
    outputModel = Model();

    // All the inertial will be of one kg in the origin
    iDynTree::SpatialInertia inertiaDefault(1.0,Position::Zero(),RotationalInertiaRaw::Zero());
    iDynTree::Link linkDefault;
    linkDefault.setInertia(inertiaDefault);

    // Then we create a base link
    outputModel.addLink("link0",linkDefault);

    std::string previousLinkName = "link0";

    for(int i=0; i<inputChain.getNrOfDOFs(); i++)
    {
        std::string addedJointName = "joint"+intToString(i);
        std::string addedLinkName = "link"+intToString(i+1);

        Transform parent_H_child;
        Axis rot_axis_in_parent_frame;
        Axis axisOnZThroughOrigin = Axis(Direction(0,0,1.0),Position::Zero());
        Transform dhTransform = TransformFromDH(inputChain(i).A,inputChain(i).Alpha,inputChain(i).D,inputChain(i).Offset);
        if( i == 0 )
        {
            // account for H0
            parent_H_child = inputChain.getH0()*dhTransform;
            rot_axis_in_parent_frame = inputChain.getH0()*axisOnZThroughOrigin;
        }
        else if( i == inputChain.getNrOfDOFs()-1 )
        {
            // Normal case
            parent_H_child = dhTransform*inputChain.getHN();
            rot_axis_in_parent_frame = axisOnZThroughOrigin;
        }
        else
        {
            // Normal case
            parent_H_child = dhTransform;
            rot_axis_in_parent_frame = axisOnZThroughOrigin;
        }

        RevoluteJoint addedJoint = RevoluteJoint(parent_H_child,rot_axis_in_parent_frame);
        outputModel.addJointAndLink(previousLinkName,
                                    addedJointName,&(addedJoint),
                                    addedLinkName,linkDefault);

        previousLinkName = addedLinkName;
    }

    return true;
}



}