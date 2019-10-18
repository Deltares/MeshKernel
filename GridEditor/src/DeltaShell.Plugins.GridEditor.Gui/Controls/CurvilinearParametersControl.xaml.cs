using System.Diagnostics.CodeAnalysis;
using System.Windows;
using System.Windows.Controls;
using DeltaShell.Plugins.GridEditor.GridGeomStateful.Api;

namespace DeltaShell.Plugins.GridEditor.Gui.Controls
{
    /// <summary>
    /// Interaction logic for CurvilinearParametersControl.xaml
    /// </summary>
    public partial class CurvilinearParametersControl : UserControl
    {
        public static readonly DependencyProperty CurvilinearParametersProperty = DependencyProperty.Register(
            "CurvilinearParameters", typeof(CurvilinearParameters), typeof(CurvilinearParametersControl), new PropertyMetadata(default(CurvilinearParameters), PropertyChangedCallback));

        public CurvilinearParametersControl()
        {
            InitializeComponent();
        }

        [ExcludeFromCodeCoverage]
        public CurvilinearParameters CurvilinearParameters
        {
            get { return (CurvilinearParameters)GetValue(CurvilinearParametersProperty); }
            set { SetValue(CurvilinearParametersProperty, value); }
        }

        private static void PropertyChangedCallback(DependencyObject d, DependencyPropertyChangedEventArgs e)
        {
            var control = d as CurvilinearParametersControl;
            if (control == null) return;

            if (e.NewValue is CurvilinearParameters samplesRefineParameters)
            {
                control.DataContext = samplesRefineParameters;
            }
        }
    }
}
