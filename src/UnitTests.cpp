#include "CppUnitTest.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTests
{
    TEST_CLASS(UnitTests)
    {
    public:

        TEST_METHOD(TestMethod2)
        {
            std::string string1 = "C++Rocks";
            std::string string2 = "C++Rocks";
            Assert::AreEqual(string1, string2);
        }
    };
}
