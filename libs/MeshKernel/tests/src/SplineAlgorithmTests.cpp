#include <gtest/gtest.h>

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/LandBoundary.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/SplineAlgorithms.hpp>

TEST(SplineAlgorithms, SnapToLandBoundaryTest)
{
    // Test the algorithm for snapping splines to land boundaries.

    constexpr double tolerance = 1.0e-8;

    constexpr int iterationCountForTest = 5;

    // The land boundary to which the spline is to be snapped.
    std::vector<meshkernel::Point> landBoundaryPoints{{257.002197, 442.130066},
                                                      {518.753845, 301.128662},
                                                      {938.006470, 416.629822}};

    meshkernel::LandBoundary landBoundary(landBoundaryPoints);

    // The original spline points.
    std::vector<meshkernel::Point> splinePoints{{281.0023, 447.3801},
                                                {367.2529, 401.6296},
                                                {461.7534, 354.3792},
                                                {517.2538, 318.3788},
                                                {614.0045, 338.629},
                                                {720.5051, 377.6294},
                                                {827.7558, 417.3798},
                                                {923.7563, 424.1299}};

    // The expected spline values after snapping to land boundary.
    std::vector<meshkernel::Point> expectedSplinePoints{{273.5868719643935, 434.2730022174478},
                                                        {359.5998304717778, 386.1712239047134},
                                                        {451.5303458337523, 338.3551703843473},
                                                        {517.7962262926076, 306.3259738916997},
                                                        {616.7325138813335, 327.9627689164845},
                                                        {725.7358644094627, 358.0902879743862},
                                                        {836.262785315633, 388.6415116416172},
                                                        {923.500177844106, 412.5818685325169}};

    // Second derivative values of the spline at the spline points.
    std::vector<meshkernel::Point> splineDerivative = meshkernel::SplineAlgorithms::SecondOrderDerivative(splinePoints, 0, static_cast<meshkernel::UInt>(splinePoints.size()) - 1);

    // Snap the spline to the land boundary
    meshkernel::SplineAlgorithms::SnapSplineToBoundary(splinePoints, splineDerivative, landBoundary, meshkernel::Projection::Type::Cartesian, iterationCountForTest);

    // The number of points in the spline should be unchanged
    ASSERT_EQ(splinePoints.size(), expectedSplinePoints.size()) << ", expected the number of points to be unchanged.";

    // Now test the snapped spline points are close to the expected values
    for (size_t i = 0; i < splinePoints.size(); ++i)
    {
        EXPECT_TRUE(meshkernel::IsEqual(expectedSplinePoints[i].x, splinePoints[i].x, tolerance))
            << "Expected x-value: " << expectedSplinePoints[i].x << ", actual: " << splinePoints[i].x << ", relative tolerance: " << tolerance;
        EXPECT_TRUE(meshkernel::IsEqual(expectedSplinePoints[i].y, splinePoints[i].y, tolerance))
            << "Expected y-value: " << expectedSplinePoints[i].y << ", actual: " << splinePoints[i].y << ", relative tolerance: " << tolerance;
    }
}

TEST(SplineAlgorithms, SnapToLandBoundaryTestMoreComplex)
{
    // Test the algorithm for snapping splines to land boundaries.

    constexpr double tolerance = 1.0e-8;

    constexpr int iterationCountForTest = 5;

    // The land boundary to which the spline is to be snapped.
    std::vector<meshkernel::Point> landBoundaryPoints{{42.0, 423.0},
                                                      {257.002197, 442.130066},
                                                      {518.753845, 331.128662},
                                                      {938.006470, 416.629822},
                                                      {1321.0, 438.0}};

    meshkernel::LandBoundary landBoundary(landBoundaryPoints);

    // The original spline points.
    std::vector<meshkernel::Point> splinePoints{{39.0, 420.0},
                                                {149.5, 436.68},
                                                {241.0023, 447.3801},
                                                {367.2529, 401.6296},
                                                {461.7534, 354.3792},
                                                {517.2538, 318.3788},
                                                {614.0045, 338.629},
                                                {720.5051, 377.6294},
                                                {827.7558, 417.3798},
                                                {923.7563, 424.1299},
                                                {1092.6, 436.18},
                                                {1210.3, 437.14}};

    // The expected spline values after snapping to land boundary.
    std::vector<meshkernel::Point> expectedSplinePoints{{38.42192722461485, 424.8066056739462},
                                                        {149.4212294844463, 431.4279980793879},
                                                        {243.4366472927346, 439.4248664919662},
                                                        {363.3393876763865, 397.2884569516308},
                                                        {463.52960329766, 354.9022144214105},
                                                        {519.8125580315373, 334.4846251591024},
                                                        {610.0830875529675, 349.4756505685123},
                                                        {722.3285941717641, 372.9106676262226},
                                                        {832.0186256015455, 394.6007897085619},
                                                        {925.9583042672142, 413.2531703957142},
                                                        {1093.033003418322, 425.553458360068},
                                                        {1210.129613351779, 431.5895114864473}};

    // Second derivative values of the spline at the spline points.
    std::vector<meshkernel::Point> splineDerivative = meshkernel::SplineAlgorithms::SecondOrderDerivative(splinePoints, 0, static_cast<meshkernel::UInt>(splinePoints.size()) - 1);

    // Snap the spline to the land boundary
    meshkernel::SplineAlgorithms::SnapSplineToBoundary(splinePoints, splineDerivative, landBoundary, meshkernel::Projection::Type::Cartesian, iterationCountForTest);

    // The number of points in the spline should be unchanged
    ASSERT_EQ(splinePoints.size(), expectedSplinePoints.size()) << ", expected the number of points to be unchanged.";

    for (size_t i = 0; i < splinePoints.size(); ++i)
    {
        EXPECT_TRUE(meshkernel::IsEqual(expectedSplinePoints[i].x, splinePoints[i].x, tolerance))
            << "Expected x-value: " << expectedSplinePoints[i].x << ", actual: " << splinePoints[i].x << ", relative tolerance: " << tolerance;
        EXPECT_TRUE(meshkernel::IsEqual(expectedSplinePoints[i].y, splinePoints[i].y, tolerance))
            << "Expected y-value: " << expectedSplinePoints[i].y << ", actual: " << splinePoints[i].y << ", relative tolerance: " << tolerance;
    }
}

TEST(SplineAlgorithms, SnapToLandBoundarExceptionalCasesTest)
{
    // Test the algorithm for snapping splines to land boundaries
    // Should fail with no spline or derivative points or different length vectors

    // The values of the spline points and the land boundary do not really matter in this test,
    // as they will not be used.

    // The land boundary to which the spline is to be snapped.
    std::vector<meshkernel::Point> landBoundaryPoints{{0.0, 0.0}, {1.0, 1.0}};

    meshkernel::LandBoundary landBoundary(landBoundaryPoints);

    std::vector<meshkernel::Point> splinePoints{{281.0023, 447.3801},
                                                {367.2529, 401.6296},
                                                {461.7534, 354.3792}};

    std::vector<meshkernel::Point> splineDerivative;

    //--------------------------------
    // Should throw exception ConstraintError, splineDerivative is empty vector
    EXPECT_THROW(meshkernel::SplineAlgorithms::SnapSplineToBoundary(splinePoints, splineDerivative, landBoundary, meshkernel::Projection::Type::Cartesian),
                 meshkernel::ConstraintError);

    splineDerivative = meshkernel::SplineAlgorithms::SecondOrderDerivative(splinePoints, 0, static_cast<meshkernel::UInt>(splinePoints.size()) - 1);

    // Add point to end of spline derivative, so number of points is not the same as spline
    // The value of the point does not matter, the SnapSplineToBoundary should throw an exception
    splineDerivative.push_back({100.0, 100.0});

    //--------------------------------
    // Should throw exception ConstraintError, spline and splineDerivative are not same length
    EXPECT_THROW(meshkernel::SplineAlgorithms::SnapSplineToBoundary(splinePoints, splineDerivative, landBoundary, meshkernel::Projection::Type::Cartesian),
                 meshkernel::ConstraintError);

    splinePoints.clear();

    //--------------------------------
    // Should throw exception ConstraintError, spline is empty vector
    EXPECT_THROW(meshkernel::SplineAlgorithms::SnapSplineToBoundary(splinePoints, splineDerivative, landBoundary, meshkernel::Projection::Type::Cartesian),
                 meshkernel::ConstraintError);
}
