using System.Drawing;
using DelftTools.TestUtils;
using DeltaShell.Plugins.GridEditor.Gui.MapLayers;
using NUnit.Framework;

namespace DeltaShell.Plugins.GridEditor.Gui.Tests.MapLayers
{
    [TestFixture]
    public class DrawingHelperTest
    {
        [Test]
        public void GivenDrawingHelper_DrawingUsingGenerateEllipseBitmap_ShouldGiveCorrectBitmap()
        {   
            //Arrange
            var pen = new Pen(Color.Blue, 1);
            var brush = new SolidBrush(Color.Coral);
            var path = TestHelper.GetTestFilePath("GenerateEllipseBitmapReference.png");

            // Act
            
            using (var bitmap = DrawingHelper.GenerateEllipseBitmap(pen, brush, 16, 16))
            {
                // Assert
                GuiTestHelper.CompareBitmapWithReferenceFile(bitmap, path);
            }
        }
    }
}   