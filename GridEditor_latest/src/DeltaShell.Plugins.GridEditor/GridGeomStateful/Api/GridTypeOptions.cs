using System.ComponentModel;
using DelftTools.Utils;

namespace DeltaShell.Plugins.GridEditor.GridGeomStateful.Api
{
    [TypeConverter(typeof(EnumDescriptionAttributeTypeConverter))]
    public enum GridTypeOptions
    {
        [Description("Square grid")]
        Square = 0,

        [Description("Wiber grid")]
        Wieber = 1,

        [Description("Hexagonal grid type one")]
        HexagonalTypeOne = 2,

        [Description("Hexagonal grid type two")]
        HexagonalTypeTwo = 3,

        [Description("Triangular grid")]
        Triangular = 4
    }
}
