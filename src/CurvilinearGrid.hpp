#pragma once

#include <vector>
#include <algorithm>
#include "Entities.hpp"

namespace GridGeom 
{
    class CurvilinearGrid
    {

    public:

        bool IncreaseGrid(const int m, const int n) 
        {

            int mMax = m + 1;
            int nMax = n + 1;
            int mnMax = std::max(m, n);

            m_grid.resize(mMax);
            for (int i = 0; i < m_grid.size(); ++i)
            {
                m_grid[i].resize(nMax);
            }

            return true;
        }

        bool Set(const std::vector<std::vector<Point>>& grid )
        {
            m_grid = grid;
            return true;
        }

        std::vector<std::vector<Point>> m_grid;

        int m_n;
        int m_m;

    };
}