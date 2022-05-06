using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using Syncfusion.UI.Xaml.Charts;
using Syncfusion.Windows.Tools.Controls;

namespace ONDAS_GUI
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : RibbonWindow
    {
        public MainWindow()
        {
            InitializeComponent();
            MainRibbon.LoadRibbonState();
            //_ = MainDock.LoadDockState();

            OutputBlock.Inlines.Add("=================\nSimulation status\n=================\n Iteration no.1821 COMPLETE\n Total time			=	0.100015 s\n===========================================================\nJob complete\n============\n End time			Mon Nov  1 11:50:10 2021\n CPU time accumulated		1.000000 s (895.000000 ms)\n Real time elapsed		11.000000 s (11713.000000 ms)");

            List<Change> changes = new();
            changes.Add(new() { Model = "Model1", Parameter = "Temperature", Previous = 250, New = 300 });
            changes.Add(new() { Model = "Model1", Parameter = "Pressure", Previous = 0, New = 1.23 });
            changes.Add(new() { Model = "Model1", Parameter = "Radius", Previous = 0, New = 0.05 });
            changes.Add(new() { Model = "Model1", Parameter = "Length", Previous = 100, New = 1000 });
            lvVersionChanges.ItemsSource = changes;

            var velocity = new double[] { 0.000000, -0.007574, -0.082428, -0.156407, -0.162508, -0.505783, -0.591622, -0.726569, -0.909640, -0.864997, -1.398059, -1.162415, -1.845595, -1.421567, -1.396582, -1.052465, -0.927407, -0.756535, -0.498979, -0.169457, 0.000000 };
            var flow = new double[] { 0.000000, -0.000041, -0.000450, -0.000856, -0.000891, -0.002791, -0.003282, -0.004068, -0.005136, -0.004949, -0.008221, -0.007080, -0.011702, -0.009232, -0.009210, -0.006986, -0.006186, -0.005062, -0.003349, -0.001138, 0.000000 };
            var pressure = new double[] { 3.169640, 3.169586, 3.170714, 3.172817, 3.175852, 3.179763, 3.183299, 3.189583, 3.194671, 3.201934, 3.209730, 3.220615, 3.223213, 3.222565, 3.224821, 3.226803, 3.229929, 3.232723, 3.235013, 3.235529, 3.234152 };
            var temperature = new double[] { 930.284836, 930.278509, 929.569686, 928.266961, 927.150860, 922.533276, 918.874655, 912.095486, 905.969210, 896.044493, 874.042245, 846.629066, 813.981565, 794.509404, 782.954273, 778.439463, 775.395684, 773.534802, 771.866882, 771.632411, 771.558949 };
            
            List<Result> results = new();
            for (int i = 0; i < velocity.Length; i++)
            {
                results.Add(new Result() { Velocity = velocity[i], Flow = flow[i], Pressure = pressure[i], Temperature = temperature[i] });
            }
            dgResults.ItemsSource = results;

            NumericalAxis xAxis = new() { Header = "Velocity in m/s" };
            Chart1.PrimaryAxis = xAxis;
            NumericalAxis yAxis = new() { Header = "Pressure in bar" };
            Chart1.SecondaryAxis = yAxis;
            LineSeries series = new();
            series.ItemsSource = results;
            series.XBindingPath = "Velocity";
            series.YBindingPath = "Pressure";
            Chart1.Series.Add(series);

        }


        public class Change
        {
            public string Model { get; set; }
            public string Parameter { get; set; }
            public double Previous { get; set; }
            public double New { get; set; }
        }

        public class Result
        {
            public double Velocity { get; set; }
            public double Flow { get; set; }
            public double Pressure { get; set; }
            public double Temperature { get; set; }
        }
    }
}
