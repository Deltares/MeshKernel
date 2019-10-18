using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using DelftTools.Shell.Core;
using DelftTools.Shell.Core.Workflow.DataItems;
using DeltaShell.Plugins.GridEditor.Importers;
using Mono.Addins;
using NetTopologySuite.Extensions.Grids;

namespace DeltaShell.Plugins.GridEditor
{
    [Extension(typeof(IPlugin))]
    public class GridEditorApplicationPlugin : ApplicationPlugin
    {
        public override string Name
        {
            get { return "GridEditor"; }
        }

        public override string DisplayName
        {
            get { return "Grid Editor Plugin"; }
        }

        public override string Description
        {
            get { return "Editor for creating and changing unstructured grids"; }
        }

        public override string Version
        {
            get { return GetVersionInfo().ProductVersion; }
        }

        public override string FileFormatVersion
        {
            get { return GetVersionInfo().FileVersion; }
        }

        public override IEnumerable<IFileImporter> GetFileImporters()
        {
            yield return new NetFileFileImporter();
        }

        public override IEnumerable<DataItemInfo> GetDataItemInfos()
        {
            yield return new DataItemInfo
            {
                Name = "Unstructured grid",
                Description = "",
                Image = new Bitmap(32,32),
                Category = "Grids",
                CreateData = p => new UnstructuredGrid(),
                ValueType = typeof(UnstructuredGrid)
            };
        }

        private FileVersionInfo GetVersionInfo()
        {
            return FileVersionInfo.GetVersionInfo(GetType().Assembly.Location);
        }

        public override IEnumerable<ProjectTemplate> ProjectTemplates()
        {
            yield return new ProjectTemplate
            {
                Id = "NewGrid",
                Name = "New unstructured grid",
                Description = "Creates a project that contains an empty grid.\r\n" + 
                              "This is meant for creating an unstructured grid from scratch.",
                Category = "Grid editing",
                ExecuteTemplate = (project, settings) =>
                {
                    project.RootFolder.Add(new DataItem(new UnstructuredGrid(), "New grid"));
                }
            };
        }
    }
}
