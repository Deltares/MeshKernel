using System.Collections.Generic;
using System.IO;
using DelftTools.TestUtils;
using DeltaShell.Plugins.GridEditor.Importers;
using GeoAPI.Extensions.Coverages;
using NUnit.Framework;

namespace DeltaShell.Plugins.GridEditor.Tests.Importers
{
    [TestFixture]
    public class AscReaderFastTest
    {
        [Test]
        public void TryingToReadWithInValidPath_ShouldReturnEmptyList()
        {
            // Act
            IList<IPointValue> points = null;
            TestHelper.AssertLogMessagesCount(() =>
            {
                points = AscReaderFast.Read("");
            },1);
            

            // Assert
            Assert.NotNull(points);
            Assert.AreEqual(0, points.Count);
        }

        [Test]
        public void TryingToReadALockedFile_ShouldReturnEmptyList()
        {
            //Arrange
            var path = TestHelper.GetTestFilePath("LockedFile.asc");

            // Act
            IList<IPointValue> points = null;

            using (var stream = new FileStream(path,FileMode.OpenOrCreate))
            {
                TestHelper.AssertLogMessagesCount(() =>
                {
                    points = AscReaderFast.Read(path);
                }, 1);
            }

            // Assert
            Assert.NotNull(points);
            Assert.AreEqual(0, points.Count);
        }

        [Test]
        public void TryingToReadFileWithInvalidValues_ShouldReturnEmptyList()
        {
            //Arrange
            var path = TestHelper.GetTestFilePath("SchematisatieInt9x9Invalid.asc");

            // Act
            IList<IPointValue> points = null;

            TestHelper.AssertLogMessagesCount(() =>
                {
                    points = AscReaderFast.Read(path);
                }, 1);

            // Assert
            Assert.NotNull(points);
            Assert.AreEqual(0, points.Count);
        }

        [Test]
        public void GivenAscFile_ReadValues_ShouldReturnCorrectPointValueList()
        {
            //Arrange
            var path = TestHelper.GetTestFilePath("SchematisatieInt9x9.asc");

            // Act
            var points = AscReaderFast.Read(path);

            // Assert
            Assert.AreEqual(81, points.Count);
        }

        [Test, Category(TestCategory.Performance)]
        public void ImportPerformance()
        {
            //Arrange
            var path = TestHelper.GetTestFilePath("PerformanceTest.asc");

            // Act
            IList<IPointValue> points = null;

            TestHelper.AssertIsFasterThan(400,
                () =>
                {
                    points = AscReaderFast.Read(path);
                });

            // Assert
            Assert.NotNull(points);
            Assert.AreEqual(137484, points.Count);
        }
    }
}