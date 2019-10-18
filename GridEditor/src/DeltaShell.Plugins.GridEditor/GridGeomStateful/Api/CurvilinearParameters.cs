using ProtoBuf;

namespace DeltaShell.Plugins.GridEditor.GridGeomStateful.Api
{
    [ProtoContract(AsReferenceDefault = true)]
    public class CurvilinearParameters
    {
        public static CurvilinearParameters CreateDefault()
        {
            return new CurvilinearParameters
            {
                MRefinement = 2000,
                NRefinement = 40,
                SmoothingIterations = 10,
                SmoothingParameter = 0.5,
                AttractionParameter = 0.0
            };
        }

        /// <summary>
        /// *M-refinement factor for regular grid generation (2000) 
        /// </summary>
        [ProtoMember(1)]
        public int MRefinement { get; set; }

        /// <summary>
        /// *N-refinement factor for regular grid generation (40) 
        /// </summary>
        [ProtoMember(2)]
        public int NRefinement { get; set; }

        /// <summary>
        /// *Nr. of inner iterations in regular grid smoothing (10).
        /// </summary>
        [ProtoMember(3)]
        public int SmoothingIterations { get; set; }

        /// <summary>
        /// * Smoothing parameter (0.5).
        /// </summary>
        [ProtoMember(4)]
        public double SmoothingParameter { get; set; }

        /// <summary>
        /// *Attraction/repulsion parameter  (0.0).
        /// </summary>
        [ProtoMember(5)]
        public double AttractionParameter { get; set; }

    }
}
