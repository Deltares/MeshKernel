using System.Diagnostics.CodeAnalysis;
using System.Windows;
using System.Windows.Controls;
using DeltaShell.Plugins.GridEditor.GridGeomStateful.Api;

namespace DeltaShell.Plugins.GridEditor.Gui.Controls
{
    /// <summary>
    /// Interaction logic for OrthogonalizationParametersControl.xaml
    /// </summary>
    public partial class OrthogonalizationParametersControl : UserControl
    {
        public static readonly DependencyProperty OrthogonalizationParametersProperty = DependencyProperty.Register(
            "OrthogonalizationParameters", typeof(OrthogonalizationParameters), typeof(OrthogonalizationParametersControl), new PropertyMetadata(default(OrthogonalizationParameters), PropertyChangedCallback));

        public OrthogonalizationParametersControl()
        {
            InitializeComponent();
        }

        [ExcludeFromCodeCoverage]
        public OrthogonalizationParameters OrthogonalizationParameters
        {
            get { return (OrthogonalizationParameters)GetValue(OrthogonalizationParametersProperty); }
            set { SetValue(OrthogonalizationParametersProperty, value); }
        }

        private static void PropertyChangedCallback(DependencyObject d, DependencyPropertyChangedEventArgs e)
        {
            var control = d as OrthogonalizationParametersControl;
            if (control == null) return;

            control.DataContext = e.NewValue;
        }
    }
}
