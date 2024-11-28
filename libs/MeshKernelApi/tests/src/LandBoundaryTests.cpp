#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "MeshKernel/MeshTransformation.hpp"
#include "MeshKernel/Parameters.hpp"

#include "MeshKernelApi/BoundingBox.hpp"
#include "MeshKernelApi/GeometryList.hpp"
#include "MeshKernelApi/Mesh1D.hpp"
#include "MeshKernelApi/Mesh2D.hpp"
#include "MeshKernelApi/MeshKernel.hpp"
#include "Version/Version.hpp"

#include "TestUtils/Definitions.hpp"
#include "TestUtils/MakeCurvilinearGrids.hpp"
#include "TestUtils/MakeMeshes.hpp"
#include "TestUtils/SampleFileReader.hpp"

#include <memory>
#include <numeric>

#include "CartesianApiTestFixture.hpp"

TEST(LandBoundaryTests, MKernelSnapSplineToLandBoundary_ShouldSnap)
{
    const double tolerance = 1e-6;

    // Setup
    int meshKernelId = 0;
    int setProjectionType = 0; // Cartesian
    auto errorCode = meshkernelapi::mkernel_allocate_state(setProjectionType, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // The land boundary to which the spline is to be snapped.
    std::vector<double> landBoundaryPointsX{257.002197, 518.753845, 938.006470};
    std::vector<double> landBoundaryPointsY{442.130066, 301.128662, 416.629822};

    // The original spline points.
    std::vector<double> splinePointsX{281.0023, 367.2529, 461.7534, 517.2538, 614.0045, 720.5051, 827.7558, 923.7563};
    std::vector<double> splinePointsY{447.3801, 401.6296, 354.3792, 318.3788, 338.629, 377.6294, 417.3798, 424.1299};

    // The expected spline values after snapping to land boundary.
    std::vector<double> expectedSplinePointsX{273.5868719643935, 359.5998304717778, 451.5303458337523, 517.7962262926076,
                                              616.7325138813335, 725.7358644094627, 836.2627853156330, 923.5001778441060};

    std::vector<double> expectedSplinePointsY{434.2730022174478, 386.1712239047134, 338.3551703843473, 306.3259738916997,
                                              327.9627689164845, 358.0902879743862, 388.6415116416172, 412.5818685325169};

    meshkernelapi::GeometryList landBoundaryGeometry{};
    landBoundaryGeometry.geometry_separator = meshkernel::constants::missing::doubleValue;
    landBoundaryGeometry.coordinates_x = landBoundaryPointsX.data();
    landBoundaryGeometry.coordinates_y = landBoundaryPointsY.data();
    landBoundaryGeometry.num_coordinates = static_cast<int>(landBoundaryPointsX.size());

    meshkernelapi::GeometryList splineGeometry{};
    splineGeometry.geometry_separator = meshkernel::constants::missing::doubleValue;
    splineGeometry.coordinates_x = splinePointsX.data();
    splineGeometry.coordinates_y = splinePointsY.data();
    splineGeometry.num_coordinates = static_cast<int>(splinePointsX.size());

    errorCode = mkernel_splines_snap_to_landboundary(meshKernelId, landBoundaryGeometry, splineGeometry, 0, static_cast<int>(splinePointsX.size() - 1));
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    for (size_t i = 0; i < splinePointsX.size(); ++i)
    {
        EXPECT_NEAR(splineGeometry.coordinates_x[i], expectedSplinePointsX[i], tolerance);
    }

    for (size_t i = 0; i < splinePointsX.size(); ++i)
    {
        EXPECT_NEAR(splineGeometry.coordinates_y[i], expectedSplinePointsY[i], tolerance);
    }
}

TEST(LandBoundaryTests, MKernelSnapSplineToLandBoundary_ShouldThrowException)
{

    // Setup
    int meshKernelId = 0;
    int setProjectionType = 0; // Cartesian
    auto errorCode = meshkernelapi::mkernel_allocate_state(setProjectionType, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // The land boundary to which the spline is to be snapped.
    std::vector<double> landBoundaryPointsX{257.002197, 518.753845, 938.006470};
    std::vector<double> landBoundaryPointsY{442.130066, 301.128662, 416.629822};

    // The original spline points.
    std::vector<double> splinePointsX{281.0023, 367.2529, 461.7534, 517.2538, 614.0045, 720.5051, 827.7558, 923.7563};
    std::vector<double> splinePointsY{447.3801, 401.6296, 354.3792, 318.3788, 338.629, 377.6294, 417.3798, 424.1299};

    meshkernelapi::GeometryList landBoundaryGeometry{};
    landBoundaryGeometry.geometry_separator = meshkernel::constants::missing::doubleValue;

    meshkernelapi::GeometryList splineGeometry{};
    splineGeometry.geometry_separator = meshkernel::constants::missing::doubleValue;

    //--------------------------------
    // Start index is less than 0
    errorCode = mkernel_splines_snap_to_landboundary(meshKernelId,
                                                     landBoundaryGeometry,
                                                     splineGeometry,
                                                     -2,
                                                     1);
    EXPECT_EQ(meshkernel::ExitCode::ConstraintErrorCode, errorCode);

    //--------------------------------
    // Start index is greater than end index
    errorCode = mkernel_splines_snap_to_landboundary(meshKernelId,
                                                     landBoundaryGeometry,
                                                     splineGeometry,
                                                     2,
                                                     1);
    EXPECT_EQ(meshkernel::ExitCode::ConstraintErrorCode, errorCode);

    //--------------------------------
    // The land boundary is not set
    errorCode = mkernel_splines_snap_to_landboundary(meshKernelId,
                                                     landBoundaryGeometry,
                                                     splineGeometry,
                                                     0,
                                                     static_cast<int>(splinePointsX.size() - 1));
    EXPECT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    // First define the number of land boundary points
    landBoundaryGeometry.num_coordinates = static_cast<int>(landBoundaryPointsX.size());

    // The land boundary points are null
    errorCode = mkernel_splines_snap_to_landboundary(meshKernelId,
                                                     landBoundaryGeometry,
                                                     splineGeometry, 0,
                                                     static_cast<int>(splinePointsX.size() - 1));
    EXPECT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    //--------------------------------
    // Now define the land boundary
    landBoundaryGeometry.coordinates_x = landBoundaryPointsX.data();
    landBoundaryGeometry.coordinates_y = landBoundaryPointsY.data();

    // The number of spline points is 0
    errorCode = mkernel_splines_snap_to_landboundary(meshKernelId,
                                                     landBoundaryGeometry,
                                                     splineGeometry,
                                                     0,
                                                     static_cast<int>(splinePointsX.size() - 1));
    EXPECT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    // define the number of spline points
    splineGeometry.num_coordinates = static_cast<int>(splinePointsX.size());

    // The spline values are null
    errorCode = mkernel_splines_snap_to_landboundary(meshKernelId,
                                                     landBoundaryGeometry,
                                                     splineGeometry,
                                                     0,
                                                     static_cast<int>(splinePointsX.size() - 1));
    EXPECT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    splineGeometry.coordinates_x = splinePointsX.data();
    splineGeometry.coordinates_y = splinePointsY.data();

    // Start spline index is greater than the number of spline points
    errorCode = mkernel_splines_snap_to_landboundary(meshKernelId,
                                                     landBoundaryGeometry,
                                                     splineGeometry,
                                                     static_cast<int>(splinePointsX.size()) + 1,
                                                     static_cast<int>(splinePointsX.size()) + 2);
    EXPECT_EQ(meshkernel::ExitCode::ConstraintErrorCode, errorCode);

    // End spline index is greater than the number of spline points
    errorCode = mkernel_splines_snap_to_landboundary(meshKernelId,
                                                     landBoundaryGeometry,
                                                     splineGeometry,
                                                     0,
                                                     static_cast<int>(splinePointsX.size()));
    EXPECT_EQ(meshkernel::ExitCode::ConstraintErrorCode, errorCode);
}

TEST(LandBoundaryTests, PolygonSnapToLandboundary_ShouldSnapPolygonToLandBoundary)
{
    const double tolerance = 1e-6;

    // Setup
    int meshKernelId = 0;
    int setProjectionType = 0; // Cartesian
    auto errorCode = meshkernelapi::mkernel_allocate_state(setProjectionType, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> landBoundaryPointsX{139.251465, 527.753906, 580.254211, 194.001801};
    std::vector<double> landBoundaryPointsY{497.630615, 499.880676, 265.878296, 212.627762};

    std::vector<double> polygonPointsX{170.001648, 263.002228, 344.002747,
                                       458.753448, 515.753845, 524.753906,
                                       510.503754, 557.754089, 545.004028,
                                       446.003387, 340.252716, 242.752106,
                                       170.001648};
    std::vector<double> polygonPointsY{472.880371, 472.880371, 475.130432,
                                       482.630493, 487.130554, 434.630005,
                                       367.129333, 297.378601, 270.378357,
                                       259.128235, 244.128067, 226.877884,
                                       472.880371};

    // The expected polygon values after snapping to land boundary.
    std::vector<double> expectedSnappedPointX = {169.8572772242283, 262.8547378163090, 343.8655709877979,
                                                 458.6558591358565, 515.6804060372598, 541.5480568270806,
                                                 555.2836667233159, 572.4472626165707, 546.2703464583593,
                                                 447.5942143903486, 341.7865993173012, 243.7707524316129,
                                                 169.8572772242283};

    std::vector<double> expectedSnappedPointY = {497.8078724305628, 498.3464789799546, 498.8156634613377,
                                                 499.4804859264834, 499.8107507986815, 438.3979070214996,
                                                 377.1760644727631, 300.6751319852315, 261.1931241088368,
                                                 247.5891786326750, 233.0020541046851, 219.4891385810638,
                                                 497.8078724305628};

    meshkernelapi::GeometryList landBoundaryGeometry{};
    landBoundaryGeometry.geometry_separator = meshkernel::constants::missing::doubleValue;
    landBoundaryGeometry.coordinates_x = landBoundaryPointsX.data();
    landBoundaryGeometry.coordinates_y = landBoundaryPointsY.data();
    landBoundaryGeometry.num_coordinates = static_cast<int>(landBoundaryPointsX.size());

    meshkernelapi::GeometryList polygonGeometry{};
    polygonGeometry.geometry_separator = meshkernel::constants::missing::doubleValue;
    polygonGeometry.coordinates_x = polygonPointsX.data();
    polygonGeometry.coordinates_y = polygonPointsY.data();
    polygonGeometry.num_coordinates = static_cast<int>(polygonPointsX.size());

    errorCode = meshkernelapi::mkernel_polygon_snap_to_landboundary(meshKernelId,
                                                                    landBoundaryGeometry,
                                                                    polygonGeometry,
                                                                    0,
                                                                    static_cast<int>(polygonPointsX.size()) - 1);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    for (size_t i = 0; i < polygonPointsX.size(); ++i)
    {
        EXPECT_NEAR(polygonGeometry.coordinates_x[i], expectedSnappedPointX[i], tolerance);
    }

    for (size_t i = 0; i < polygonPointsX.size(); ++i)
    {
        EXPECT_NEAR(polygonGeometry.coordinates_y[i], expectedSnappedPointY[i], tolerance);
    }
}

TEST_F(CartesianApiTestFixture, Mesh2DSnapToLandboundary_ShouldSnapToLandBoundary)
{
    // Prepare
    MakeMesh(30, 30, 1);
    auto const meshKernelId = GetMeshKernelId();
    meshkernelapi::GeometryList landBoundaries{};
    std::vector landBoundariesX{-1.47, 9.57, 26.50, 36.86};
    std::vector landBoundariesY{31.64, 32.71, 32.13, 31.83};
    landBoundaries.coordinates_x = landBoundariesX.data();
    landBoundaries.coordinates_y = landBoundariesY.data();
    landBoundaries.num_coordinates = static_cast<int>(landBoundariesX.size());

    meshkernelapi::GeometryList selectingPolygon{};
    std::vector selectingPolygonX{16.7, -8.48, -7.07, 61.60, 42.23, 16.7};
    std::vector selectingPolygonY{41.59, 38.24, -8.06, -9.82, 43.53, 41.59};
    selectingPolygon.coordinates_x = selectingPolygonX.data();
    selectingPolygon.coordinates_y = selectingPolygonY.data();
    selectingPolygon.num_coordinates = static_cast<int>(selectingPolygonX.size());

    auto errorCode = meshkernelapi::mkernel_mesh2d_snap_to_landboundary(meshKernelId, selectingPolygon, landBoundaries);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_deallocate_state(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}
