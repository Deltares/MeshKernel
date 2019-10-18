using System.Diagnostics;
using System.Linq;
using DelftTools.Shell.Core;
using DelftTools.Shell.Core.Workflow.DataItems;
using DeltaShell.Plugins.GridEditor.Importers;
using NetTopologySuite.Extensions.Grids;
using NUnit.Framework;

namespace DeltaShell.Plugins.GridEditor.Tests
{
    [TestFixture]
    public class GridEditorApplicationPluginTest
    {
        [Test]
        public void TestDefaultSettings()
        {
            var plugin = new GridEditorApplicationPlugin();

            Assert.IsTrue(plugin?.Name?.Length > 0, "Plugin name should be specified");
            Assert.IsTrue(plugin?.DisplayName?.Length > 0, "Plugin display name should be specified");
            Assert.IsTrue(plugin?.Description?.Length > 0, "Plugin description should be specified");

            var versionInfo = FileVersionInfo.GetVersionInfo(typeof(GridEditorApplicationPlugin).Assembly.Location);
            Assert.AreEqual(versionInfo.ProductVersion, plugin.Version);
            Assert.AreEqual(versionInfo.FileVersion, plugin.FileFormatVersion);
        }

        [Test]
        public void WhenCallingGetFileImporters_ExpectedFileImportersShouldBeReturned()
        {
            //Arrange
            var plugin = new GridEditorApplicationPlugin();

            // Act
            var importers = plugin.GetFileImporters().ToList();

            // Assert
            Assert.AreEqual(1, importers.Count);
            Assert.IsInstanceOf<NetFileFileImporter>(importers[0]);
        }

        [Test]
        public void WhenCallingGetDataItemInfos_ExpectedDataItemInfosShouldBeReturned()
        {
            //Arrange
            var plugin = new GridEditorApplicationPlugin();

            // Act
            var dataItemInfos = plugin.GetDataItemInfos().ToList();

            // Assert
            Assert.AreEqual(1, dataItemInfos.Count);
            Assert.AreEqual("Unstructured grid", dataItemInfos[0].Name);
            Assert.AreEqual(typeof(UnstructuredGrid), dataItemInfos[0].ValueType);

            var grid = dataItemInfos[0].CreateData(null);

            Assert.IsInstanceOf<UnstructuredGrid>(grid);
        }

        [Test]
        public void WhenCallingProjectTemplates_ExpectedTemplatesShouldBeReturned()
        {
            //Arrange
            var project = new Project();
            var plugin = new GridEditorApplicationPlugin();

            // Act
            var templates = plugin.ProjectTemplates().ToArray();
            
            // Assert
            Assert.AreEqual(1, templates.Length);
            Assert.AreEqual("NewGrid", templates[0].Id);

            templates[0].ExecuteTemplate(project, null);

            var dataItem = project.RootFolder.Items[0] as IDataItem;

            Assert.NotNull(dataItem, "Added item to root folder should be a IDataItem");
            Assert.IsInstanceOf<UnstructuredGrid>(dataItem.Value, "DataItem value should be an unstructured grid");
        }
    }
}