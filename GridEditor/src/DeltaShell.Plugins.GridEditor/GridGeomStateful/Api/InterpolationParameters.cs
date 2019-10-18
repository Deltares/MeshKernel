using ProtoBuf;

namespace DeltaShell.Plugins.GridEditor.GridGeomStateful.Api
{
    [ProtoContract(AsReferenceDefault = true)]
    public class InterpolationParameters
    {
        public static InterpolationParameters CreateDefault()
        {
            return new InterpolationParameters
            {
                InterpolationType = 1,
                DisplayInterpolationProcess = 0,
                MaxNumberOfRefinementIterations = 1,
                AveragingMethod = 1,
                MinimumNumberOfPoints = 1,
                RelativeSearchRadius = 1.01,
                InterpolateTo = 2
            };
        }

        /// <summary>
        /// *Actual interpolation type (1)
        /// </summary>
        [ProtoMember(1)]
        public int InterpolationType { get; set; }

        /// <summary>
        /// Variable related to interactor behaviour (0) 
        /// </summary>
        [ProtoMember(2)]
        public int DisplayInterpolationProcess { get; set; }

        /// <summary>
        /// *Maximum number of refinement iterations, set to 1 if only one refinement is wanted (10) 
        /// </summary>
        [ProtoMember(3)]
        public int MaxNumberOfRefinementIterations { get; set; }

        /// <summary>
        /// *Averaging method : 1 = simple averaging, 2 = closest point, 3 = max, 4 = min, 5 = inverse weighted distance, 6 = minabs, 7 = kdtree (1)
        /// </summary>
        [ProtoMember(4)]
        public int AveragingMethod { get; set; }

        /// <summary>
        /// *Minimum number of points needed inside cell to handle the cell (1)
        /// </summary>
        [ProtoMember(5)]
        public int MinimumNumberOfPoints { get; set; }

        /// <summary>
        /// *Relative search cell size, default 1= actual cell size, 2= twice as large, search radius can be larger than cell so more sample are included. (1.01)
        /// </summary>
        [ProtoMember(6)]
        public double RelativeSearchRadius { get; set; }

        /// <summary>
        /// *Interpolation settings, 1=bathy, 2=zk, 3=s1, 4=Zc (2)
        /// </summary>
        [ProtoMember(7)]
        public int InterpolateTo { get; set; }

    }
}
