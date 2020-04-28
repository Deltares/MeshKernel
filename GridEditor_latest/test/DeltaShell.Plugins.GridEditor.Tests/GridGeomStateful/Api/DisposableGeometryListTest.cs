using System;
using DeltaShell.Plugins.GridEditor.GridGeomStateful.Api;
using NUnit.Framework;

namespace DeltaShell.Plugins.GridEditor.Tests.GridGeomStateful.Api
{
    [TestFixture]
    public class DisposableGeometryListTest
    {
        [Test]
        public void WhenCallingCreateGeometryListNative_AGeometryListNativeShouldBeCreated()
        {
            //Arrange
            var list = new DisposableGeometryList
            {
                XCoordinates = new[] {0.0, 1.0, 2.0, -999.0, 0.0, 1.0, 2.0 },
                YCoordinates = new[] {0.0, 1.0, 2.0, -999.0, 0.0, 1.0, 2.0 },
                ZCoordinates = new[] {0.0, 1.0, 2.0, -999.0, 0.0, 1.0, 2.0 },
                NumberOfCoordinates = 7,
                GeometrySeparator = -999.0
            };

            using (list)
            {
                // Act
                var nativeGeometryList = list.CreateGeometryListNative();

                // Assert

                Assert.AreEqual(list.NumberOfCoordinates, nativeGeometryList.numberOfCoordinates);
                Assert.AreNotEqual(IntPtr.Zero, nativeGeometryList.xCoordinates);
                Assert.AreNotEqual(IntPtr.Zero, nativeGeometryList.yCoordinates);
                Assert.AreNotEqual(IntPtr.Zero, nativeGeometryList.zCoordinates);
            }
        }
    }
}