#pragma once
#include "Entities.hpp"

// include boost
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>


#include <vector>
#include <memory>

// kd-tree
// https://gist.github.com/alexcpn/1f187f2114976e748f4d3ad38dea17e8

// r-tree
// https://gist.github.com/logc/10272165

namespace GridGeom
{
    namespace SpatialTrees 
    {
        namespace bg = boost::geometry;
        namespace bgi = boost::geometry::index;

        class RTree
        {

            typedef bg::model::point<double, 2, bg::cs::cartesian> Point2D;
            typedef bg::model::box<Point2D> Box2D;
            typedef std::pair<Point2D, int> value2D;
            typedef bgi::rtree<value2D, bgi::quadratic<16>> RTree2D;

            typedef bg::model::point<double, 3, bg::cs::cartesian> Point3D;
            typedef bg::model::box<Point2D> Box3D;
            typedef std::pair<Point3D, int> value3D;
            typedef bgi::rtree<value3D, bgi::quadratic<16>> RTree3D;

        public:

            bool BuildTree(std::vector<Point>& nodes, const Projections projection)
            {
                m_points.resize(nodes.size());
                for (int n = 0; n < nodes.size(); ++n)
                {
                    m_points[n] = std::make_pair(Point2D{ nodes[n].x, nodes[n].y }, n);
                }
                m_rtree2D = RTree2D(m_points.begin(), m_points.end());
                return true;
            }

            std::vector<int> NearestNeighbours(Point node, const double searchRadius) const
            {
                  
                std::vector<value2D> queryResult;
                Box2D box(Point2D(node.x - searchRadius, node.y - searchRadius), Point2D(node.x + searchRadius, node.y + searchRadius));
                Point2D nodeSought = Point2D(node.x, node.y);
                m_rtree2D.query(
                    bgi::within(box) &&
                    bgi::satisfies([&](value2D const& v) {return bg::distance(v.first, nodeSought) < searchRadius; }),
                    std::back_inserter(queryResult));

                std::vector<int> result(queryResult.size());
                for (size_t i = 0; i < queryResult.size(); i++)
                {
                    result[i] = queryResult[i].second;
                }
                return std::move(result);
            }

            bool RemoveNode(int position) 
            {
                int numberRemoved = m_rtree2D.remove(m_points[position]);
                if (numberRemoved != 1) 
                {
                    return false;
                }
                return true;
            }

            int Size()
            {
                return m_rtree2D.size();
            }


        private:

            Projections m_projection;
            RTree2D m_rtree2D;
            std::vector<std::pair<Point2D, int>> m_points;

        };   

        class KDTree
        {

            class Node
            {
               public:
                   using NodePtr = std::shared_ptr< Node >;
                   size_t index;
                   std::vector< double > x;
                   NodePtr left;
                   NodePtr right;
            };



        public:
            bool BuildTree(std::vector<Point>& nodes, const Projections projection)
            {
                CreateKDTree(nodes);
                return true;
            }

            bool CreateKDTree(std::vector<Point>& nodes)
            {
                return true;
            }

            std::vector<int> NearestNeighbours(Point node, const double searchRadius)
            {
                return  std::move(std::vector<int>(0));
            } 
        
        };

    }
}


