using System.ComponentModel;
using System.Drawing.Design;
using DelftTools.Shell.Gui;
using DelftTools.Utils;
using GeoAPI.Extensions.CoordinateSystems;
using NetTopologySuite.Extensions.Grids;
using SharpMap.UI.Forms;

namespace DeltaShell.Plugins.GridEditor.Gui
{
    public class UnstructuredGridProperties : ObjectProperties<UnstructuredGrid>
    {
        [Description("Status")]
        [Category("General")]
        [PropertyOrder(2)]
        [TypeConverter(typeof(CoordinateSystemStringTypeConverter))]
        [Editor(typeof(CoordinateSystemTypeEditor), typeof(UITypeEditor))]
        public ICoordinateSystem CoordinateSystem
        {
            get
            {
                return data.CoordinateSystem;
            }
            set
            {
                data.CoordinateSystem = value;
            }
        }
    }
}