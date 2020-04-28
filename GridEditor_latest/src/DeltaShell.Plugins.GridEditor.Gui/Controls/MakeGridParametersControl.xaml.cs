using System.Windows;
using System.Windows.Controls;
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

        public static readonly DependencyProperty IsBasedOnPolygonProperty = DependencyProperty.Register(
            "IsBasedOnPolygon", typeof(bool), typeof(MakeGridParametersControl), new PropertyMetadata(default(bool), PropertyChangedCallback));
        
        public MakeGridParametersControl()
        {
            InitializeComponent();
        }

        public MakeGridParameters MakeGridParameters
        {
            get { return (MakeGridParameters)GetValue(MakeGridParametersProperty); }
            set { SetValue(MakeGridParametersProperty, value); }
        }

        public bool IsBasedOnPolygon
        {
            get { return (bool)GetValue(IsBasedOnPolygonProperty); }
            set { SetValue(IsBasedOnPolygonProperty, value); }
        }

        private static void PropertyChangedCallback(DependencyObject d, DependencyPropertyChangedEventArgs e)
        {
            var viewModel = (d as MakeGridParametersControl)?.ViewModel;
            if (viewModel == null) return;

            if (e.Property == MakeGridParametersProperty)
            {
                viewModel.GridParameters = e.NewValue as MakeGridParameters;
            }

            if (e.Property == IsBasedOnPolygonProperty)
            {
                viewModel.IsBasedOnPolygon = (bool)e.NewValue;
            }
        }

        private void ButtonBase_OnClick(object sender, RoutedEventArgs e)
        {
            var flyout = FindFlyoutToCollapse();
            flyout.Visibility = Visibility.Collapsed;

            var dialog = new OpenFileDialog { Filter = "ascii file|*.asc|tiff file|*.tif" };

            var result = dialog.ShowDialog() == DialogResult.OK ? dialog.FileName : null;

            flyout.Visibility = Visibility.Visible;

            if (!string.IsNullOrEmpty(result))
            {
                ViewModel.AsciiFilePath = result;
            }
        }

        private Grid FindFlyoutToCollapse()
        {
            FrameworkElement target = this;
            Grid grid;
            do
            {
                target = FindLogicalTreeParent<Grid>(target);
                grid = (Grid) target;
            } while (grid != null && !grid.Name.Equals("GenerateMenuButtonContent"));

            return grid;
        }

        private static T FindLogicalTreeParent<T>(DependencyObject dependencyObject) where T : class
        {
            if (dependencyObject == null)
            {
                return null;
            }
            
            var target = dependencyObject;
            do
            {
                target = LogicalTreeHelper.GetParent(target);
            } while (target != null && !(target is T));

            return target as T;
        }
    }
}
