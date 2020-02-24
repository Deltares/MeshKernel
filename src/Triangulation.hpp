#pragma once

#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>
#include "Entities.hpp"



namespace GridGeom
{
    class Triangulation
    {
    public:
        Triangulation() : m_projection(Projections::cartesian)
        {
        };





        Projections m_projection;
    };

}
