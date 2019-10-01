#pragma once

#include "Gridgeom.hpp"
#include "Mesh.cpp"

 int ggeo_new_grid(int& gridStateId)
 {
     gridStateId = meshInstances.size() + 1;
     meshInstances.reserve(gridStateId);
     return 0;
 };

 int ggeo_deallocate_state(int& gridStateId)
 {
     meshInstances.erase(meshInstances.begin() + gridStateId);
     return 0;
 }

 int ggeo_set_state(int gridStateId, MeshGeometryDimensions& meshGeometryDimensions, MeshGeometry& meshGeometry, bool IsGeographic)
 {
     if (IsGeographic)
     {
         auto instance = std::make_unique<GridGeom::Mesh<GridGeom::cartesianPoint>>();

         std::vector<GridGeom::Edge> edges(meshGeometryDimensions.numedge);
         for (int e = 0; e < edges.size(); e++)
         {
             edges[e].first = 0;
             edges[e].second = 0;
         }

         std::vector<GridGeom::cartesianPoint> nodes(meshGeometryDimensions.numnode);
         for (int n = 0; n < edges.size(); n++)
         {
             nodes[n].x = meshGeometry.nodex[n];
             nodes[n].y = meshGeometry.nodey[n];
         }

         instance->setState(edges, nodes);

         meshInstances[gridStateId] = std::move(instance);

     }
     return 0;
 }
