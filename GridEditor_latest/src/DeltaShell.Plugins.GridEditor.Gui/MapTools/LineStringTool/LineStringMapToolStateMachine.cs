using System.Collections.Generic;

namespace DeltaShell.Plugins.GridEditor.Gui.MapTools.LineStringTool
{
    /// <summary>
    /// State machine for handling state transitions for <see cref="LineStringMapTool"/>
    /// </summary>
    internal class LineStringMapToolStateMachine
    {
        private class StateTransition
        {
            private readonly LineStringMapToolState currentState;
            private readonly LineStringMapToolStateCommand command;

            public StateTransition(LineStringMapToolState currentState, LineStringMapToolStateCommand command)
            {
                this.currentState = currentState;
                this.command = command;
            }

            public override int GetHashCode()
            {
                return 17 + 31 * currentState.GetHashCode() + 31 * command.GetHashCode();
            }

            public override bool Equals(object obj)
            {
                var other = obj as StateTransition;
                return other != null && currentState == other.currentState && command == other.command;
            }
        }

        // defines all transitions between states (current state + command => new state)
        private readonly Dictionary<StateTransition, LineStringMapToolState> transitions = new Dictionary<StateTransition, LineStringMapToolState>
        {
            { new StateTransition(LineStringMapToolState.Selecting, LineStringMapToolStateCommand.MouseDownNotOnFeatureNoSelection), LineStringMapToolState.AddingNew},
            { new StateTransition(LineStringMapToolState.Selecting, LineStringMapToolStateCommand.AfterMouseDownFeatureSelected), LineStringMapToolState.Editing},

            { new StateTransition(LineStringMapToolState.AddingNew, LineStringMapToolStateCommand.AfterMouseDoubleClick), LineStringMapToolState.Selecting},

            { new StateTransition(LineStringMapToolState.Editing, LineStringMapToolStateCommand.MouseDownOnFeatureOnTracker), LineStringMapToolState.Dragging},
            { new StateTransition(LineStringMapToolState.Editing, LineStringMapToolStateCommand.MouseDownNotOnFeatureHasSelection), LineStringMapToolState.Selecting},
            { new StateTransition(LineStringMapToolState.Editing, LineStringMapToolStateCommand.MouseDownNotOnFeatureNoSelection), LineStringMapToolState.Selecting},
            { new StateTransition(LineStringMapToolState.Editing, LineStringMapToolStateCommand.MouseDownOnOtherFeature), LineStringMapToolState.Selecting},

            { new StateTransition(LineStringMapToolState.Dragging, LineStringMapToolStateCommand.AfterMouseUp), LineStringMapToolState.Editing}
        };

        /// <summary>
        /// Current <see cref="LineStringMapToolState"/>
        /// </summary>
        public LineStringMapToolState CurrentState { get; private set; } = LineStringMapToolState.Selecting;

        /// <summary>
        /// This function will return the state that should be transitioned to
        /// </summary>
        /// <param name="command"></param>
        /// <returns></returns>
        public LineStringMapToolState MoveNext(LineStringMapToolStateCommand command)
        {
            CurrentState = GetNextState(command);
            return CurrentState;
        }

        private LineStringMapToolState GetNextState(LineStringMapToolStateCommand command)
        {
            var stateTransition = new StateTransition(CurrentState, command);
            return !transitions.TryGetValue(stateTransition, out var nextState) ? CurrentState : nextState;
        }
    }
}