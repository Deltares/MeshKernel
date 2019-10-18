using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using DeltaShell.Plugins.GridEditor.GridGeomStateful.Native;
using ProtoBuf;

namespace DeltaShell.Plugins.GridEditor.GridGeomStateful.Api
{
    [ProtoContract]
    public sealed class DisposableGeometryList : IDisposable
    {
        private readonly List<GCHandle> objectGarbageCollectHandles = new List<GCHandle>();

        [ProtoMember(1)] private double geometrySeparator;
        [ProtoMember(2)] private double innerOuterSeparator;
        [ProtoMember(3)] private int numberOfCoordinates;
        [ProtoMember(4)] private double[] xCoordinates;
        [ProtoMember(5)] private double[] yCoordinates;
        [ProtoMember(6)] private double[] zCoordinates;

        /// <summary>
        /// Separator used for separating multiple geometries
        /// </summary>
        
        public double GeometrySeparator
        {
            get { return geometrySeparator; }
            set { geometrySeparator = value; }
        }

        /// <summary>
        /// Separator used for separating the outer and inner rings of a polygon geometry
        /// </summary>
        
        public double InnerOuterSeparator
        {
            get { return innerOuterSeparator; }
            set { innerOuterSeparator = value; }
        }

        /// <summary>
        /// Total number of coordinates (including separators)
        /// </summary>
        
        public int NumberOfCoordinates
        {
            get { return numberOfCoordinates; }
            set { numberOfCoordinates = value; }
        }

        /// <summary>
        /// The x coordinates of the geometries
        /// </summary>
        
        public double[] XCoordinates
        {
            get { return xCoordinates; }
            set { xCoordinates = value; }
        }

        /// <summary>
        /// The y coordinates of the geometries
        /// </summary>
        
        public double[] YCoordinates
        {
            get { return yCoordinates; }
            set { yCoordinates = value; }
        }

        /// <summary>
        /// The z coordinates of the geometries
        /// </summary>
        
        public double[] ZCoordinates
        {
            get { return zCoordinates; }
            set { zCoordinates = value; }
        }

        private bool IsMemoryPinned
        {
            get { return objectGarbageCollectHandles.Count > 0; }
        }

        /// <summary>
        /// Creates a <see cref="GeometryListNative"/> for the specified geometries
        /// </summary>
        /// <returns></returns>
        public GeometryListNative CreateGeometryListNative()
        {
            if (!IsMemoryPinned)
            {
                PinMemory();
            }

            var lookup = objectGarbageCollectHandles.ToDictionary(h => h.Target, h => h);

            return new GeometryListNative
            {
                xCoordinates = lookup[XCoordinates].AddrOfPinnedObject(),
                yCoordinates = lookup[YCoordinates].AddrOfPinnedObject(),
                zCoordinates = lookup[ZCoordinates].AddrOfPinnedObject(),
                numberOfCoordinates = NumberOfCoordinates
            };
        }

        public void Dispose()
        {
            UnPinMemory();
        }

        private void UnPinMemory()
        {
            foreach (var handle in objectGarbageCollectHandles)
            {
                handle.Free();
            }

            objectGarbageCollectHandles.Clear();
        }

        private void PinMemory()
        {
            // compensate for null arrays
            XCoordinates = GetArray(XCoordinates);
            YCoordinates = GetArray(YCoordinates);
            ZCoordinates = GetArray(ZCoordinates);

            objectGarbageCollectHandles.Add(GCHandle.Alloc(XCoordinates, GCHandleType.Pinned));
            objectGarbageCollectHandles.Add(GCHandle.Alloc(YCoordinates, GCHandleType.Pinned));
            objectGarbageCollectHandles.Add(GCHandle.Alloc(ZCoordinates, GCHandleType.Pinned));
        }

        private static T[] GetArray<T>(T[] array)
        {
            return array ?? new T[0];
        }
    }
}