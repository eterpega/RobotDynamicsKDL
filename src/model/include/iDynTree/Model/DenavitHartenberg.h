/*
 * Copyright (C) 2016 Fondazione Istituto Italiano di Tecnologia
 * Authors: Silvio Traversaro
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#ifndef IDYNTREE_DENAVIT_HARTENBERG_H
#define IDYNTREE_DENAVIT_HARTENBERG_H

#include <iDynTree/Model/Indeces.h>
#include <iDynTree/Model/Model.h>

#include <iDynTree/Core/Transform.h>

#include <vector>

namespace iDynTree
{

/**
 * Structure representing
 * the four DH parameters of
 * a link in a Chain.
 */
struct DHLink
{
   double A;
   double D;
   double Alpha;
   double Offset;
};

/**
 * Simple representation of a chain described with
 * Denavit-Hartenberg Parameters.
 * Directly inspiered by the iCub::iKin::iKinChain class.
 */
class DHChain
{
private:
    Transform H0;
    std::vector<DHLink> dhParams;
    Transform HN;

public:
    void   setNrOfDOFs(size_t nDofs);
    size_t getNrOfDOFs() const;

    void setH0(const iDynTree::Transform & _H0);
    const iDynTree::Transform & getH0() const;

    void setHN(const Transform & _HN);
    const Transform & getHN() const;

    /**
     * Return a reference to the i-th link of the chain.
     */
    DHLink & operator() (const size_t i);

    /**
     * Return a reference to the i-th link of the chain (const version).
     */
    const DHLink & operator() (const size_t i) const;
};

/**
 * Extract a Denavit Hartenberg Chain from a iDynTree::Model.
 * In particular you can specify the base frame and the end effector frame
 * of the chain. The chain is then extracted using the algorithm found
 * in:
 * Chapter 3, Spong, Mark W., Seth Hutchinson, and M. Vidyasagar. "Robot modeling and control". 2006.
 * (section 3.2.3 of http://www.cs.duke.edu/brd/Teaching/Bio/asmb/current/Papers/chap3-forward-kinematics.pdf)
 *
 *  \note The DH representation supports only revolute and translational
 *        1-DOF joints. In this implementation only revolute joints are supported.
 */
bool ExtractDHChainFromModel(const iDynTree::Model & model,
                             const std::string baseFrame,
                             const std::string eeFrame,
                                   DHChain & outputChain);

}


#endif