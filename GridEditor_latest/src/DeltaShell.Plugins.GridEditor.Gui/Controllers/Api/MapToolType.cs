namespace DeltaShell.Plugins.GridEditor.Gui.Controllers.Api
{
    public enum MapToolType
    {
        /// <summary>
        /// No map tool
        /// </summary>
        None,
        /// <summary>
        /// Polygon map tool
        /// </summary>
        Polygon,
        /// <summary>
        /// Spline map tool
        /// </summary>
        Spline,
        /// <summary>
        /// Land boundaries map tool
        /// </summary>
        LandBoundaries,
        /// <summary>
        /// Insert edges map tool
        /// </summary>
        InsertVertices,
        /// <summary>
        /// Delete vertex map tool
        /// </summary>
        DeleteVertices,
        /// <summary>
        /// Delete edges map tool
        /// </summary>
        DeleteEdges,
        /// <summary>
        /// Move vertices map tool
        /// </summary>
        MoveVertices
    }
}