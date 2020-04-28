using System.Collections.Generic;
using System.Linq;
using DelftTools.Shell.Core;
using DeltaShell.Plugins.GridEditor.Importers;
using GeoAPI.Geometries;
using NetTopologySuite.Extensions.Grids;
using NUnit.Framework;

namespace DeltaShell.Plugins.GridEditor.Tests.Importers
{
    [TestFixture]
    public class NetFileFileImporterTest
    {
        [Test]
        public void CheckNetFileImporterDefaultProperties()
        {
            //Arrange
            var importer = new NetFileFileImporter();
            
            // Assert
            Assert.AreEqual("Net file", importer.Name);
            Assert.AreEqual("Grid importers", importer.Category);
            Assert.AreEqual("Net file", importer.Description);
            Assert.Null(importer.Image);
            Assert.Contains(typeof(UnstructuredGrid), importer.SupportedItemTypes.ToArray());
            Assert.IsTrue(importer.CanImportOnRootLevel);
            Assert.IsTrue(importer.OpenViewAfterImport);
            Assert.AreEqual("Net file|*.nc", importer.FileFilter);
        }

        [Test]
        public void GivenType_Importer_ShouldReturnIfImportCanBeDone()
        {
            //Arrange
            var importer = new NetFileFileImporter();
            
            // Assert
            Assert.IsTrue(importer.CanImportOn(new UnstructuredGrid()));
            Assert.IsTrue(importer.CanImportOn(new Project()));
            Assert.IsFalse(importer.CanImportOn(""));
            Assert.IsFalse(importer.CanImportOn(null));
        }

        [Test]
        public void GivenAnUnstructuredGrid_Importing_ShouldOverrideValuesOfExistingGrid()
        {
            //Arrange
            var path = "TestPath";
            var grid = new UnstructuredGrid
            {
                Vertices = new List<Coordinate>
                {
                    new Coordinate(0,0),
                    new Coordinate(1,1),
                },
                Edges = new List<Edge>(new List<Edge>
                {
                    new Edge(0,1)
                })
            };

            var importer = new NetFileFileImporter
            {
                ImportGridFromFile = p =>
                {
                    Assert.AreEqual(path, p);
                    return new UnstructuredGrid();
                }
            };

            // Act
            
            importer.ImportItem(path, grid);

            // Assert
            Assert.AreEqual(0, grid.Vertices.Count);
            Assert.AreEqual(0, grid.Edges.Count);
        }

        [Test]
        public void GivenAProject_Importing_ShouldCreateANewGrid()
        {
            //Arrange
            var path = "TestPath";
            var grid = new UnstructuredGrid
            {
                Vertices = new List<Coordinate>
                {
                    new Coordinate(0,0),
                    new Coordinate(1,1),
                },
                Edges = new List<Edge>(new List<Edge>
                {
                    new Edge(0,1)
                })
            };

            var importer = new NetFileFileImporter
            {
                ImportGridFromFile = p =>
                {
                    Assert.AreEqual(path, p);
                    return grid;
                }
            };

            // Act
            var newGrid = importer.ImportItem(path, new Project());

            // Assert
            Assert.AreEqual(grid, newGrid);
        }
    }
}