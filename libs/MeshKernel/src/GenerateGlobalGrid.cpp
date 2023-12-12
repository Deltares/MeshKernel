#include "MeshKernel/GenerateGlobalGrid.hpp"
#include <cmath>
#include <iostream>

double meshkernel::GenerateGlobalGrid::getDeltaY(const double y, const double deltaX)
{

    double deltaY = deltaX * std::cos(constants::conversion::degToRad * y);

    for (int i = 0; i < 5; ++i)
    {
        double phi = constants::conversion::degToRad * (y + 0.5 * deltaY);
        double c = std::cos(phi);
        double s = std::sqrt(1.0 - c * c);
        double f = deltaY - deltaX * c;
        double df = 1.0 + 0.5 * constants::conversion::degToRad * deltaX * s;
        double yd = f / df;
        deltaY = deltaY - yd;

        // TODO is this a good convergence tolerance
        if (yd < 1.0e-14)
        {
            break;
        }
    }

    return deltaY;
}

void meshkernel::GenerateGlobalGrid::isNodeDB(const Mesh& mesh, const Point& x, UInt& kp)
{
    constexpr double tolerance = 1.0e-6;

    kp = constants::missing::uintValue;

    for (int i = static_cast<int>(mesh.m_nodes.size() - 1); i >= 0; --i)
    {
        if (IsEqual(x, mesh.m_nodes[i], tolerance))
        {
            kp = static_cast<UInt>(i);
            break;
        }
    }
}

void meshkernel::GenerateGlobalGrid::addMaze(Mesh& mesh, const std::array<Point, 8>& points, const double ySign, const UInt pointCount, const bool jafive)
{
    // std::array<Point, 8> newPoints;
    std::array<UInt, 8> kk;

    // UInt numk = 0;

    for (UInt i = 0; i < pointCount; ++i)
    {
        Point p = {points[i].x, ySign * points[i].y};
        isNodeDB(mesh, p, kk[i]);

        if (kk[i] == constants::missing::uintValue)
        {
            // newPoints[numk] = points[i];
            // kc[numk] = 1;
            kk[i] = mesh.InsertNode(p);
            // ++numk;
        }
    }

    for (UInt i = 0; i < pointCount; ++i)
    {
        UInt k2 = i + 1;

        if (i == pointCount - 1)
        {
            k2 = 0;
        }

        mesh.ConnectNodes(kk[i], kk[k2]);
    }

    if (jafive)
    {
        mesh.ConnectNodes(kk[1], kk[3]);
        mesh.ConnectNodes(kk[1], kk[4]);
    }
}

void meshkernel::GenerateGlobalGrid::Compute(const UInt nx, const UInt ny, Mesh2D& mesh)
{

    constexpr double xwleft = -180.0;

    constexpr double dxdouble = 30000.0;

    std::array<Point, 8> points;
    UInt numberOfPoints = 0;
    [[maybe_unused]] UInt n12;

    double dx0 = 360.0 / static_cast<double>(nx);
    double dy0 = dx0;
    double yy = 0.0;

    bool jafive = false;
    bool jaklaar = false;

    [[maybe_unused]] Projection projection = Projection::sphericalAccurate;
    mesh.m_projection = projection;

    for (UInt i = 0; i < ny; ++i)
    {
        dy0 = getDeltaY(yy, dx0);

        if (yy + 1.5 * dy0 > 90.0)
        {
            dy0 = 90.0 - yy;
            jaklaar = true;
            jafive = false;
        }
        else
        {
            if (((dy0 * constants::conversion::degToRad * constants::geometric::earth_radius) < dxdouble) && !jafive)
            {
                dx0 = 2.0 * dx0;
                jafive = true;
                dy0 = getDeltaY(yy, dx0);
            }
            else
            {
                jafive = false;
                // Not sure if this is a global variable and used elsewhere
                n12 = 0;
            }
            if (yy + 1.5 * dy0 > 90.0)
            {
                dy0 = 0.51 * (90.0 - yy);
            }
        }

        for (UInt j = 1; j <= nx; ++j)
        {
            double xx = static_cast<double>(j - 1) * dx0 + xwleft;

            points[0] = {xx, yy};

            if (!jafive)
            {
                points[1] = {xx + dx0, yy};
                points[2] = {xx + dx0, yy + dy0};
                points[3] = {xx, yy + dy0};
                numberOfPoints = 4;
            }
            else
            {
                points[1] = {xx + 0.5 * dx0, yy};
                points[2] = {xx + dx0, yy};
                points[3] = {xx + dx0, yy + dy0};
                points[4] = {xx, yy + dy0};
                numberOfPoints = 5;
            }

            bool isIn = true; // IsPointInPolygon(pnt, polygon);

            if (isIn && points[2].x <= xwleft + 360.0)
            {
                jafive = false;
                addMaze(mesh, points, 1.0, numberOfPoints, jafive);
                addMaze(mesh, points, -1.0, numberOfPoints, jafive);
            }
        }

        if (jaklaar)
        {
            break;
        }

        yy += dy0;
    }

    mergenodesinpolygon(mesh);

#if 0
    for (UInt l = 0; l < numL; ++l)
    {
        k1 = mesh.m_edges[l].first;
        k2 = mesh.m_edges[l].second;

        if (k1 != nullvalue && k2 != nullvalue) // jammer dan, nb na setnodadm nog zooi
        {
            if ((nmk[k1] == 5 || nmk[k1] == 6) && (nmk[k2] == 5 || nmk[k2] == 6))
            {
                if (yk[k1] == yk[k2])
                {
                    DELLINK(L);
                }
            }
        }
    }
#endif
}

void meshkernel::GenerateGlobalGrid::mergenodesinpolygon(Mesh2D& mesh [[maybe_unused]])
{
    mesh.AdministrateNodesEdges();

    UInt numberOfPoints = mesh.GetNumNodes();

    for (UInt i = 0; i < numberOfPoints; ++i)
    {
        // auto [isInPolygon, whichPolygon] = polygons.IsPointInPolygons(mesh.m_nodes[i]);
        bool isInPolygon = true;

        if (isInPolygon)
        {
            for (UInt j = 0; j < mesh.m_nodesNumEdges[i]; ++j)
            {

#if 0
                [[maybe_unused]] UInt ll = mesh.m_nodesEdges[i][j];
                // UInt ll = std::abs(mesh.m_nodesEdges[i][j]);

                if (mesh.m_edges[ll] == 1 || mesh.m_edges[ll] == 6)
                {
                    itp = 1;
                }
                else if (mesh.m_edges[ll] == 3 || mesh.m_edges[ll] == 4 || mesh.m_edges[ll] == 5 || mesh.m_edges[ll] == 7)
                {
                    itp = 1;
                }
                else if (mesh.m_edges[ll] == 2)
                {
                    itp = 2;
                }
                else
                {
                    itp = 0;
                }

                mesh.kc[i] = max(mesh.kc[i], itp);
#endif
            }
        }
    }

    // if (mesh.m_projection == Projection::sphericalAccurate)
    // {
    //     getMeshBounds();
    // }

    // [[maybe_unused]] UInt kint = std::max<UInt>(numberOfPoints / 100, 1U);
    double tooClose = 1.0e-3;

    if (tooClose > 0.0)
    {
        double searchRadius = tooClose * tooClose;
        [[maybe_unused]] bool jadone = true;
        UInt numberMerged = 0;

        // Do we have to rebuild the tree each time a point has been merged?
        mesh.BuildTree(Mesh::Location::Nodes);

        for (UInt i = 0; i < numberOfPoints; ++i)
        {
            // Do we have to rebuild the tree each time a point has been merged?
            // How does the tree handle invalid points?
            mesh.BuildTree(Mesh::Location::Nodes);
            mesh.SearchNearestLocation(mesh.m_nodes[i], searchRadius, Mesh::Location::Nodes);

            UInt count = mesh.GetNumLocations(Mesh::Location::Nodes);

            std::cout << "number of locations: " << count << std::endl;

            if (count > 1)
            {
                for (UInt j = 0; j < count; ++j)
                {
                    UInt nodeToMerge = mesh.GetLocationsIndices(j, Mesh::Location::Nodes);
                    std::cout << "merging nodes: " << i << "  " << j << "   " << mesh.GetLocationsIndices(j, Mesh::Location::Nodes) << "  "
                              << mesh.m_nodes[nodeToMerge].x << "  " << mesh.m_nodes[nodeToMerge].y << "  "
                              << std::endl;

                    if (nodeToMerge != i && mesh.m_nodes[nodeToMerge].IsValid())
                    {
                        mesh.MergeTwoNodes(nodeToMerge, i);
                        ++numberMerged;
                    }
                }
            }
        }

        std::cout << " numberMerged " << numberMerged << std::endl;
    }

    [[maybe_unused]] int dummy;

    dummy = 1;
}
