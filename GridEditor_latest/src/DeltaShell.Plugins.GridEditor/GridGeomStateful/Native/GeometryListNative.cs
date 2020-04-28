using System;
using System.Runtime.InteropServices;

namespace DeltaShell.Plugins.GridEditor.GridGeomStateful.Native
{
    [StructLayout(LayoutKind.Sequential)]
    public struct GeometryListNative
    {
        public int type;

        public double geometrySeperator;

        public double innerOuterSeperator;

        public int numberOfCoordinates;

        public IntPtr xCoordinates;

        public IntPtr yCoordinates;

        public IntPtr zCoordinates;
    }
}
