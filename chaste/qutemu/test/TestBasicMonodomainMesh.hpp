#ifndef TESTBASICMONODOMAINMESH_HPP_
#define TESTBASICMONODOMAINMESH_HPP_

#include <cxxtest/TestSuite.h>

#include "SlabSimulation.hpp"

class TestBasicMonodomainMesh : public CxxTest::TestSuite
{
public:
    void TestSimpleSimulation() throw(Exception)
    {
        SlabSimulation::BasicSquare(400, 10, 500, 100);

        HeartEventHandler::Headings();
        HeartEventHandler::Report();
    }
};

#endif /*TESTBASICMONODOMAINMESH_HPP_*/
