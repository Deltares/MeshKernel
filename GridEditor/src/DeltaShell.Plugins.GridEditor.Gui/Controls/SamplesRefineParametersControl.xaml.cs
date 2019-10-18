using System.Diagnostics.CodeAnalysis;
using System.Windows;
using System.Windows.Controls;
using DeltaShell.Plugins.GridEditor.GridGeomStateful.Api;

namespace DeltaShell.Plugins.GridEditor.Gui.Controls
{
    /// <summary>
    /// Interaction logic for SamplesRefineParametersControl.xaml
    /// </summary>
    public partial class SamplesRefineParametersControl : UserControl
    {
        public static readonly DependencyProperty SamplesRefineParametersProperty = DependencyProperty.Register(
            "SamplesRefineParameters", typeof(SamplesRefineParameters), typeof(SamplesRefineParametersControl), new PropertyMetadata(default(SamplesRefineParameters), PropertyChangedCallback));

        public static readonly DependencyProperty InterpolationParametersProperty = DependencyProperty.Register(
            "InterpolationParameters", typeof(InterpolationParameters), typeof(SamplesRefineParametersControl), new PropertyMetadata(default(InterpolationParameters), PropertyChangedCallback));

        public SamplesRefineParametersControl()
        {
            InitializeComponent();
        }

        [ExcludeFromCodeCoverage]
        public SamplesRefineParameters SamplesRefineParameters
        {
            get { return (SamplesRefineParameters)GetValue(SamplesRefineParametersProperty); }
            set { SetValue(SamplesRefineParametersProperty, value); }
        }

        [ExcludeFromCodeCoverage]
        public InterpolationParameters InterpolationParameters
        {
            get { return (InterpolationParameters)GetValue(InterpolationParametersProperty); }
            set { SetValue(InterpolationParametersProperty, value); }
        }

        private static void PropertyChangedCallback(DependencyObject d, DependencyPropertyChangedEventArgs e)
        {
            var control = d as SamplesRefineParametersControl;
            if (control == null) return;

            if (e.NewValue is SamplesRefineParameters samplesRefineParameters)
            {
                control.MainGrid.DataContext = samplesRefineParameters;
            }

            if (e.NewValue is InterpolationParameters interpolationParameters)
            {
                control.MaxNumberOfRefinementIterationsTextBox.DataContext = interpolationParameters;
            }
        }
    }
}
