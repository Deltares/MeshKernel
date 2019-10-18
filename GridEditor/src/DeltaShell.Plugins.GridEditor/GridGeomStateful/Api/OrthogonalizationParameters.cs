using ProtoBuf;

namespace DeltaShell.Plugins.GridEditor.GridGeomStateful.Api
{
    [ProtoContract(AsReferenceDefault = true)]
    public class OrthogonalizationParameters
    {
        public static OrthogonalizationParameters CreateDefault()
        {
            return new OrthogonalizationParameters
            {
                OuterIterations = 1,
                BoundaryIterations = 25,
                InnerIterations = 25,
                OrthogonalizationToSmoothingFactor = 0.975
            };
        }

        /// <summary>
        /// *Number of outer iterations in orthogonalization. Increase this parameter for complex grids (2) 
        /// </summary>
        [ProtoMember(1)]
        public int OuterIterations { get; set; }
        
        /// <summary>   
        /// *Number of boundary iterations in grid/net orthogonalization within itatp (25) 
        /// </summary>
        [ProtoMember(2)]
        public int BoundaryIterations { get; set; }
        
        /// <summary>
        /// *Number of inner iterations in grid/net orthogonalization within itbnd (25)
        /// </summary>
        [ProtoMember(3)]
        public int InnerIterations { get; set; }

        /// <summary>
        /// *Factor from 0 to 1. between grid smoothing and grid orthogonality (0.975)	 
        /// </summary>
        [ProtoMember(4)]
        public double OrthogonalizationToSmoothingFactor { get; set; }
    }
}