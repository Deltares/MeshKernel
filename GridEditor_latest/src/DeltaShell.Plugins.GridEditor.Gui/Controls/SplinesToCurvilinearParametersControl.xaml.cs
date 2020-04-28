using System.Diagnostics.CodeAnalysis;
using System.Windows;
using System.Windows.Controls;
using DeltaShell.Plugins.GridEditor.GridGeomStateful.Api;

namespace DeltaShell.Plugins.GridEditor.Gui.Controls
{
    /// <summary>
    /// Interaction logic for SplinesToCurvilinearParametersControl.xaml
    /// </summary>
    public partial class SplinesToCurvilinearParametersControl : UserControl
    {
        public static readonly DependencyProperty SplinesToCurvilinearParametersProperty = DependencyProperty.Register(
            "SplinesToCurvilinearParameters", typeof(SplinesToCurvilinearParameters), typeof(SplinesToCurvilinearParametersControl), new PropertyMetadata(default(SplinesToCurvilinearParameters), PropertyChangedCallback));

        public static readonly DependencyProperty CurvilinearParametersProperty = DependencyProperty.Register(
            "CurvilinearParameters", typeof(CurvilinearParameters), typeof(SplinesToCurvilinearParametersControl), new PropertyMetadata(default(CurvilinearParameters), PropertyChangedCallback));

        public SplinesToCurvilinearParametersControl()
        {
            InitializeComponent();
        }

        [ExcludeFromCodeCoverage]
        public CurvilinearParameters CurvilinearParameters
        {
            get { return (CurvilinearParameters)GetValue(CurvilinearParametersProperty); }
            set { SetValue(CurvilinearParametersProperty, value); }
        }

        [ExcludeFromCodeCoverage]
        public SplinesToCurvilinearParameters SplinesToCurvilinearParameters
        {
            get { return (SplinesToCurvilinearParameters)GetValue(SplinesToCurvilinearParametersProperty); }
            set { SetValue(SplinesToCurvilinearParametersProperty, value); }
        }

        private static void PropertyChangedCallback(DependencyObject d, DependencyPropertyChangedEventArgs e)
        {
            if (!(d is SplinesToCurvilinearParametersControl control)) return;

            if (e.NewValue is SplinesToCurvilinearParameters splinesToCurvilinearParameters)
            {
                control.MainGrid.DataContext = splinesToCurvilinearParameters;
            }

            if (e.NewValue is CurvilinearParameters interpolationParameters)
            {
                control.NRefinementTextBox.DataContext = interpolationParameters;
                control.MRefinementTextBox.DataContext = interpolationParameters;
            }
        }
    }
}
