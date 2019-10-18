using System.Windows;
using System.Windows.Forms;
using DeltaShell.Plugins.GridEditor.GridGeomStateful.Api;
using UserControl = System.Windows.Controls.UserControl;

namespace DeltaShell.Plugins.GridEditor.Gui.Controls
{
    /// <summary>
    /// Interaction logic for MakeGridParametersControl.xaml
    /// </summary>
    public partial class MakeGridParametersControl : UserControl
    {
        public static readonly DependencyProperty MakeGridParametersProperty = DependencyProperty.Register(
            "MakeGridParameters", typeof(MakeGridParameters), typeof(MakeGridParametersControl), new PropertyMetadata(default(MakeGridParameters), PropertyChangedCallback));

        public MakeGridParametersControl()
        {
            InitializeComponent();

            ViewModel.GetFilePath = () =>
            {
                var dialog = new OpenFileDialog {Filter = "ascii file|*.asc"};

                return dialog.ShowDialog() == DialogResult.OK
                    ? dialog.FileName
                    : null;
            };
        }

        public MakeGridParameters MakeGridParameters
        {
            get { return (MakeGridParameters)GetValue(MakeGridParametersProperty); }
            set { SetValue(MakeGridParametersProperty, value); }
        }

        private static void PropertyChangedCallback(DependencyObject d, DependencyPropertyChangedEventArgs e)
        {
            var viewModel = (d as MakeGridParametersControl)?.ViewModel;
            if (viewModel == null) return;

            viewModel.GridParameters = e.NewValue as MakeGridParameters;
        }
    }
}
