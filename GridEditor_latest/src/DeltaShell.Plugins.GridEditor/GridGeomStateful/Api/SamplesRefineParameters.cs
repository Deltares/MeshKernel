using ProtoBuf;

namespace DeltaShell.Plugins.GridEditor.GridGeomStateful.Api
{
    [ProtoContract(AsReferenceDefault = true)]
    public class SamplesRefineParameters
    {
        public static SamplesRefineParameters CreateDefault()
        {
            return new SamplesRefineParameters
            {
                SampleVectorDimension = 1,
                MinimumCellSize = 50000.0,
                DirectionalRefinement = 0,
                RefinementType = 2,
                ConnectHangingNodes = 1,
                MaximumTimeStepInCourantGrid = 0.0,
                AccountForSamplesOutside = 1
            };
        }

        /// <summary>
        /// Sample vector dimension (1) 
        /// </summary>
        [ProtoMember(1)]
        public int SampleVectorDimension { get; set; }

        /// <summary>
        /// *Minimum cell size (50000.0)
        /// </summary>
        [ProtoMember(2)]
        public double MinimumCellSize { get; set; }

        /// <summary>
        /// *Directional refinement, 1 yes 0 no (0)
        /// </summary>
        [ProtoMember(3)]
        public int DirectionalRefinement { get; set; }

        /// <summary>
        /// *Refinement criterion type, 1 ridge, 2 wave courant, 3 meshwidth (2)
        /// </summary>
        [ProtoMember(4)]
        public int RefinementType { get; set; }

        /// <summary>
        /// *Connect hanging nodes at the end of the iteration, 1 yes 0 no (1) 
        /// </summary>
        [ProtoMember(5)]
        public int ConnectHangingNodes { get; set; }

        /// <summary>
        /// *Maximum time-step in courant grid (0.0)
        /// </summary>
        [ProtoMember(6)]
        public double MaximumTimeStepInCourantGrid { get; set; }

        /// <summary>
        /// *Take samples outside cell into account , 1 yes 0 no (1)
        /// </summary>
        [ProtoMember(7)]
        public int AccountForSamplesOutside { get; set; }
    }
}