#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/RangeCheck.hpp"

#include <limits>
#include <string>
#include <utility>

#include <gtest/gtest.h>

using namespace meshkernel;

template <typename T>
class RangeCheckFixture : public testing::Test
{
public:
    void CheckEqual()
    {
        EXPECT_NO_THROW(range_check::CheckEqual<T>(1, 1, name));
        EXPECT_THROW(range_check::CheckEqual<T>(1, 2, name), RangeError);

        if (std::same_as<T, float> ||
            std::same_as<T, double> ||
            std::same_as<T, long double>)
        {
            double constexpr x = 1.0;
            double constexpr epsilon = std::numeric_limits<T>::epsilon();
            EXPECT_THROW(range_check::CheckEqual<T>(x, x - epsilon, name), RangeError);
            EXPECT_NO_THROW(range_check::CheckEqual<T>(x, x, name));
            EXPECT_THROW(range_check::CheckEqual<T>(x, x + epsilon, name), RangeError);
        }
    }

    void CheckNotEqual()
    {
        EXPECT_NO_THROW(range_check::CheckNotEqual<T>(1, 2, name));
        EXPECT_THROW(range_check::CheckNotEqual<T>(1, 1, name), RangeError);

        if (std::same_as<T, float> ||
            std::same_as<T, double> ||
            std::same_as<T, long double>)
        {
            double constexpr x = 1.0;
            double constexpr epsilon = std::numeric_limits<T>::epsilon();
            EXPECT_NO_THROW(range_check::CheckNotEqual<T>(x, x - epsilon, name));
            EXPECT_THROW(range_check::CheckNotEqual<T>(x, x, name), RangeError);
            EXPECT_NO_THROW(range_check::CheckNotEqual<T>(x, x + epsilon, name));
        }
    }

    void CheckGreaterTest()
    {
        EXPECT_NO_THROW(range_check::CheckGreater<T>(2, 1, name));
        EXPECT_THROW(range_check::CheckGreater<T>(2, 2, name), RangeError);
        EXPECT_THROW(range_check::CheckGreater<T>(3, 3, name), RangeError);
    }

    void CheckGreaterEqualTest()
    {
        EXPECT_NO_THROW(range_check::CheckGreaterEqual<T>(2, 1, name));
        EXPECT_NO_THROW(range_check::CheckGreaterEqual<T>(2, 2, name));
        EXPECT_THROW(range_check::CheckGreaterEqual<T>(2, 3, name), RangeError);
    }

    void CheckLessTest()
    {
        EXPECT_NO_THROW(range_check::CheckLess<T>(1, 2, name));
        EXPECT_THROW(range_check::CheckLess<T>(2, 2, name), RangeError);
        EXPECT_THROW(range_check::CheckLess<T>(3, 2, name), RangeError);
    }

    void CheckLessEqualTest()
    {
        EXPECT_NO_THROW(range_check::CheckLessEqual<T>(1, 2, name));
        EXPECT_NO_THROW(range_check::CheckLessEqual<T>(2, 2, name));
        EXPECT_THROW(range_check::CheckLessEqual<T>(3, 2, name), RangeError);
    }

    void CheckInClosedIntervalTest()
    {
        EXPECT_THROW(range_check::CheckInClosedInterval<T>(1, interval, name), RangeError);
        EXPECT_NO_THROW(range_check::CheckInClosedInterval<T>(3, interval, name));
        EXPECT_NO_THROW(range_check::CheckInClosedInterval<T>(4, interval, name));
        EXPECT_NO_THROW(range_check::CheckInClosedInterval<T>(5, interval, name));
        EXPECT_THROW(range_check::CheckInClosedInterval<T>(6, interval, name), RangeError);
    }

    void CheckInOpenIntervalTest()
    {
        EXPECT_THROW(range_check::CheckInOpenInterval<T>(1, interval, name), RangeError);
        EXPECT_THROW(range_check::CheckInOpenInterval<T>(3, interval, name), RangeError);
        EXPECT_NO_THROW(range_check::CheckInOpenInterval<T>(4, interval, name));
        EXPECT_THROW(range_check::CheckInOpenInterval<T>(5, interval, name), RangeError);
        EXPECT_THROW(range_check::CheckInOpenInterval<T>(6, interval, name), RangeError);
    }

    void CheckInSemiOpenFromAboveIntervalTest()
    {
        EXPECT_THROW(range_check::CheckInSemiOpenFromAboveInterval<T>(1, interval, name), RangeError);
        EXPECT_NO_THROW(range_check::CheckInSemiOpenFromAboveInterval<T>(3, interval, name));
        EXPECT_NO_THROW(range_check::CheckInSemiOpenFromAboveInterval<T>(4, interval, name));
        EXPECT_THROW(range_check::CheckInSemiOpenFromAboveInterval<T>(5, interval, name), RangeError);
        EXPECT_THROW(range_check::CheckInSemiOpenFromAboveInterval<T>(6, interval, name), RangeError);
    }

    void CheckInSemiOpenFromBelowIntervalTest()
    {

        EXPECT_THROW(range_check::CheckInSemiOpenFromBelowInterval<T>(1, interval, name), RangeError);
        EXPECT_THROW(range_check::CheckInSemiOpenFromBelowInterval<T>(3, interval, name), RangeError);
        EXPECT_NO_THROW(range_check::CheckInSemiOpenFromBelowInterval<T>(4, interval, name));
        EXPECT_NO_THROW(range_check::CheckInSemiOpenFromBelowInterval<T>(5, interval, name));
        EXPECT_THROW(range_check::CheckInSemiOpenFromBelowInterval<T>(6, interval, name), RangeError);
    }

    void CheckOneOfTest()
    {
        EXPECT_THROW(range_check::CheckOneOf<T>(1, values, name), RangeError);
        EXPECT_NO_THROW(range_check::CheckOneOf<T>(2, values, name));
        EXPECT_NO_THROW(range_check::CheckOneOf<T>(3, values, name));
        EXPECT_NO_THROW(range_check::CheckOneOf<T>(4, values, name));
        EXPECT_NO_THROW(range_check::CheckOneOf<T>(5, values, name));
        EXPECT_THROW(range_check::CheckOneOf<T>(6, values, name), RangeError);
    }

    void CheckNoneOfTest()
    {
        EXPECT_NO_THROW(range_check::CheckNoneOf<T>(1, values, name));
        EXPECT_THROW(range_check::CheckNoneOf<T>(2, values, name), RangeError);
        EXPECT_THROW(range_check::CheckNoneOf<T>(3, values, name), RangeError);
        EXPECT_THROW(range_check::CheckNoneOf<T>(4, values, name), RangeError);
        EXPECT_THROW(range_check::CheckNoneOf<T>(5, values, name), RangeError);
        EXPECT_NO_THROW(range_check::CheckNoneOf<T>(6, values, name));
    }

private:
    inline static std::string const name{"dummy"};
    std::pair<T, T> const interval{T{3}, T{5}};
    std::vector<T> const values{2, 3, 4, 5};
};

// Cannot include long when using MSVC 16 2019. It's a bug. See:
// https://github.com/microsoft/STL/issues/2765
// Fixed in MSVC 17 2022

using TestTypes = ::testing::Types<int,
                                   short,
                                   // long,
                                   int32_t, // same as int
                                   uint32_t,
                                   int64_t,
                                   uint64_t,
                                   float,
                                   double,
                                   long double>;

TYPED_TEST_SUITE(RangeCheckFixture, TestTypes);

TYPED_TEST(RangeCheckFixture, TestTypes)
{
    this->CheckGreaterTest();
    this->CheckGreaterEqualTest();
    this->CheckLessTest();
    this->CheckLessEqualTest();
    this->CheckInClosedIntervalTest();
    this->CheckInOpenIntervalTest();
    this->CheckInSemiOpenFromAboveIntervalTest();
    this->CheckInSemiOpenFromBelowIntervalTest();
    this->CheckOneOfTest();
    this->CheckNoneOfTest();
}