using System;
using System.Diagnostics.CodeAnalysis;
using System.Windows;
using System.Windows.Media;
using DelftTools.Shell.Core;
using DelftTools.Shell.Gui;

namespace DeltaShell.Plugins.GridEditor.Gui
{
    [ExcludeFromCodeCoverage]
    internal sealed class GridGraphicsProvider : IGraphicsProvider
    {
        private readonly ResourceDictionary resources = new ResourceDictionary
        {
            Source = new Uri("pack://application:,,,/DeltaShell.Plugins.GridEditor.Gui;component/Controls/Dictionaries/IconBrushes.xaml")
        };

        public bool CanProvideDrawingGroupFor(object item)
        {
            return item is GridEditorGuiPlugin ||
                   item is GridEditorApplicationPlugin ||
                   (item is ProjectTemplate template && template.Id == "NewGrid");

        }

        public DrawingGroup CreateDrawingGroupFor(object item)
        {
            if (item is GridEditorGuiPlugin || 
                item is GridEditorApplicationPlugin ||
                item is ProjectTemplate)
            {
                return (DrawingGroup)resources["GridDrawing"];
            }

            return null;
        }
    }
}