#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Mesh2D.hpp"

namespace meshkernel
{

    const std::vector<int>& GetValidProjections()
    {
        static std::vector<int> validProjections{static_cast<int>(Projection::cartesian),
                                                 static_cast<int>(Projection::spherical),
                                                 static_cast<int>(Projection::sphericalAccurate)};
        return validProjections;
    }

    const std::vector<int>& GetValidDeletionOptions()
    {
        static std::vector validMesh2DDeletionOptions{static_cast<int>(Mesh2D::DeleteMeshOptions::InsideNotIntersected),
                                                      static_cast<int>(Mesh2D::DeleteMeshOptions::InsideAndIntersected),
                                                      static_cast<int>(Mesh2D::DeleteMeshOptions::FacesWithIncludedCircumcenters)};
        return validMesh2DDeletionOptions;
    }

    Projection GetProjectionValue(const int projection)
    {
        static const int Cartesian = static_cast<int>(Projection::cartesian);
        static const int Spherical = static_cast<int>(Projection::spherical);
        static const int SphericalAccurate = static_cast<int>(Projection::sphericalAccurate);

        if (projection == Cartesian)
        {
            return Projection::cartesian;
        }
        else if (projection == Spherical)
        {
            return Projection::spherical;
        }
        else if (projection == SphericalAccurate)
        {
            return Projection::sphericalAccurate;
        }
        else
        {
            throw ConstraintError("Cannot convert integer value, {}, for projection to projection enumeration value", projection);
        }
    }

    const std::string& ProjectionToString(const Projection projection)
    {
        static const std::string Cartesian = "Projection::Cartesian";
        static const std::string Spherical = "Projection::Spherical";
        static const std::string SphericalAccurate = "Projection::SphericalAccurate";
        static const std::string Unknown = "UNKNOWN";

        switch (projection)
        {
        case Projection::cartesian:
            return Cartesian;
        case Projection::spherical:
            return Spherical;
        case Projection::sphericalAccurate:
            return SphericalAccurate;
        default:
            return Unknown;
        }
    }

    CurvilinearDirection GetCurvilinearDirectionValue(int direction)
    {
        static const int mDirection = static_cast<int>(CurvilinearDirection::M);
        static const int nDirection = static_cast<int>(CurvilinearDirection::N);

        if (direction == mDirection)
        {
            return CurvilinearDirection::M;
        }
        else if (direction == nDirection)
        {
            return CurvilinearDirection::N;
        }
        else
        {
            throw ConstraintError("Cannot convert integer value, {}, for direction to curvilinear direction enumeration value", direction);
        }
    }

    const std::string& CurvilinearDirectionToString(CurvilinearDirection direction)
    {
        static const std::string mDirection = "CurvilinearDirection::M";
        static const std::string nDirection = "CurvilinearDirection::N";
        static const std::string Unknown = "UNKNOWN";

        switch (direction)
        {
        case CurvilinearDirection::M:
            return mDirection;
        case CurvilinearDirection::N:
            return nDirection;
        default:
            return Unknown;
        }
    }
} // namespace meshkernel
