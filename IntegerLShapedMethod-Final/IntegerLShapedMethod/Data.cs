using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ILOG.Concert;
using ILOG.CPLEX;

namespace IntegerLShapedMethod
{
    class Data
    {
        public class InputDataReaderException : System.Exception
        {
            internal InputDataReaderException(string file)
                : base("'" + file + "' contains bad data format") { }
        }

        public static void getInputData_GenearalCSV(string DataFolder)
        {
            string InputFileName;
            string line;

            // Read A_Matrix
            InputFileName = "A_Matrix.csv";
            System.IO.StreamReader file = new System.IO.StreamReader(DataFolder + InputFileName);
            while ((line = file.ReadLine()) != null)
            {
                string[] arr_line = line.Split(',');
                Program._ILShaped.A_Matrix.Add(Array.ConvertAll(arr_line, double.Parse));
            }
            file.Close();

            // Read c_Vector
            InputFileName = "c_Vector.csv";
            file = new System.IO.StreamReader(DataFolder + InputFileName);
            while ((line = file.ReadLine()) != null)
            {
                string[] arr_line = line.Split(',');
                Program._ILShaped.c_Vector.AddRange(Array.ConvertAll(arr_line, double.Parse));
            }
            file.Close();

            // Read b_vector
            InputFileName = "b_Vector.csv";
            file = new System.IO.StreamReader(DataFolder + InputFileName);
            while ((line = file.ReadLine()) != null)
            {
                string[] arr_line = line.Split(',');
                Program._ILShaped.b_Vector.AddRange(Array.ConvertAll(arr_line, double.Parse));
            }
            file.Close();

            // Read Variables
            InputFileName = "FirstStageVariables.csv";
            file = new System.IO.StreamReader(DataFolder + InputFileName);
            while ((line = file.ReadLine()) != null)
            {
                string[] arr_line = line.Split(',');
                Program._ILShaped.FirStageVars.AddRange(arr_line);
            }
            file.Close();

            InputFileName = "SecondStageVariables.csv";
            file = new System.IO.StreamReader(DataFolder + InputFileName);
            while ((line = file.ReadLine()) != null)
            {
                string[] arr_line = line.Split(',');
                Program._ILShaped.SecStageVars.AddRange(arr_line);
            }
            file.Close();

            // Read Probability of each scenario
            InputFileName = "Scenario_Probability.csv";
            file = new System.IO.StreamReader(DataFolder + InputFileName);
            while ((line = file.ReadLine()) != null)
            {
                string[] arr_line = line.Split(',');
                Program._ILShaped.ScenarioProbability.AddRange(Array.ConvertAll(arr_line, double.Parse));
            }
            file.Close();
        }

        public static void getInputData_forScenarioCSV(string DataFolder, string myScenario)
        {
            string InputFileName;
            string line;

            // Read W_Martrix
            Program._ILShaped.W_Matrix.Clear();
            InputFileName = myScenario + "_W_Matrix.csv";
            System.IO.StreamReader file = new System.IO.StreamReader(DataFolder + InputFileName);
            while ((line = file.ReadLine()) != null)
            {
                string[] arr_line = line.Split(',');
                Program._ILShaped.W_Matrix.Add(Array.ConvertAll(arr_line, double.Parse));
            }
            file.Close();

            // Read T_Matrix
            Program._ILShaped.T_Matrix.Clear();
            InputFileName = myScenario + "_T_Matrix.csv";
            file = new System.IO.StreamReader(DataFolder + InputFileName);
            while ((line = file.ReadLine()) != null)
            {
                string[] arr_line = line.Split(',');
                Program._ILShaped.T_Matrix.Add(Array.ConvertAll(arr_line, double.Parse));
            }
            file.Close();

            // Read d_Matrix
            Program._ILShaped.d_Matrix.Clear();
            InputFileName = myScenario + "_d_Matrix.csv";
            file = new System.IO.StreamReader(DataFolder + InputFileName);
            while ((line = file.ReadLine()) != null)
            {
                string[] arr_line = line.Split(',');
                Program._ILShaped.d_Matrix.Add(Array.ConvertAll(arr_line, double.Parse));
            }
            file.Close();

            // Read q_Vector
            Program._ILShaped.q_Vector.Clear();
            InputFileName = myScenario + "_q2_Vector.csv";
            file = new System.IO.StreamReader(DataFolder + InputFileName);
            while ((line = file.ReadLine()) != null)
            {
                string[] arr_line = line.Split(',');
                Program._ILShaped.q_Vector.AddRange(Array.ConvertAll(arr_line, double.Parse));
            }
            file.Close();

            // Read q_Vector
            Program._ILShaped.q_Matrix.Clear();
            InputFileName = myScenario + "_q_Matrix.csv";
            file = new System.IO.StreamReader(DataFolder + InputFileName);
            while ((line = file.ReadLine()) != null)
            {
                string[] arr_line = line.Split(',');
                Program._ILShaped.q_Matrix.Add(Array.ConvertAll(arr_line, double.Parse));
            }
            file.Close();

            // Read h_Vector
            Program._ILShaped.h_Vector.Clear();
            InputFileName = myScenario + "_h_Vector.csv";
            file = new System.IO.StreamReader(DataFolder + InputFileName);
            while ((line = file.ReadLine()) != null)
            {
                string[] arr_line = line.Split(',');
                Program._ILShaped.h_Vector.AddRange(Array.ConvertAll(arr_line, double.Parse));
            }
            file.Close();
        }
        public static void getInputData_General(string DataFolder)
        {
            string InputFileName;
            string line;

            // Read A_Matrix
            InputFileName = "A_Matrix.txt";
            System.IO.StreamReader file = new System.IO.StreamReader(DataFolder + InputFileName);
            while ((line = file.ReadLine()) != null)
            {
                string[] arr_line = line.Split(Program._ILShaped.delimiterChars);
                Program._ILShaped.A_Matrix.Add(Array.ConvertAll(arr_line, double.Parse));
            }
            file.Close();

            // Read c_vector
            InputFileName = "c_Vector.txt";
            file = new System.IO.StreamReader(DataFolder + InputFileName);
            if ((line = file.ReadLine()) != null)
            {
                string[] arr_line = line.Split(Program._ILShaped.delimiterChars);
                Program._ILShaped.c_Vector.AddRange(Array.ConvertAll(arr_line, double.Parse));
            }
            file.Close();

            // Read b_vector
            InputFileName = "b_Vector.txt";
            file = new System.IO.StreamReader(DataFolder + InputFileName);
            if ((line = file.ReadLine())!=null)
            {
                string[] arr_line = line.Split(Program._ILShaped.delimiterChars);
                Program._ILShaped.b_Vector.AddRange(Array.ConvertAll(arr_line, double.Parse));
            }
            file.Close();

            // Read Variables
            InputFileName = "FirstStageVariables.txt";
            file = new System.IO.StreamReader(DataFolder + InputFileName);
            if ((line = file.ReadLine())!=null)
            {
                string[] arr_line = line.Split(Program._ILShaped.delimiterChars);
                Program._ILShaped.FirStageVars.AddRange(arr_line);
            }
            file.Close();

            InputFileName = "SecondStageVariables.txt";
            file = new System.IO.StreamReader(DataFolder + InputFileName);
            if ((line = file.ReadLine())!=null)
            {
                string[] arr_line = line.Split(Program._ILShaped.delimiterChars);
                Program._ILShaped.SecStageVars.AddRange(arr_line);
            }
            file.Close();

            // Read Probability of each scenario
            InputFileName = "Scenario_Probability.txt";
            file = new System.IO.StreamReader(DataFolder + InputFileName);
            if((line = file.ReadLine())!=null)
            {
                string[] arr_line = line.Split(Program._ILShaped.delimiterChars);
                Program._ILShaped.ScenarioProbability.AddRange(Array.ConvertAll(arr_line, double.Parse));
            }
            file.Close();
        }

        public static void getInputData_forScenario(string DataFolder, string myScenario)
        {
            string InputFileName;
            string line;

            // Read W_Martrix
            Program._ILShaped.W_Matrix.Clear();
            InputFileName = myScenario + "_W_Matrix.txt";
            System.IO.StreamReader file = new System.IO.StreamReader(DataFolder + InputFileName);
            while ((line = file.ReadLine()) != null)
            {
                string[] arr_line = line.Split(Program._ILShaped.delimiterChars);
                Program._ILShaped.W_Matrix.Add(Array.ConvertAll(arr_line, double.Parse));
            }
            file.Close();

            // Read T_Matrix
            Program._ILShaped.T_Matrix.Clear();
            InputFileName = myScenario + "_T_Matrix.txt";
            file = new System.IO.StreamReader(DataFolder + InputFileName);
            while ((line = file.ReadLine()) != null)
            {
                string[] arr_line = line.Split(Program._ILShaped.delimiterChars);
                Program._ILShaped.T_Matrix.Add(Array.ConvertAll(arr_line, double.Parse));
            }
            file.Close();

            // Read q_Vector
            Program._ILShaped.q_Vector.Clear();
            InputFileName = myScenario + "_q_Vector.txt";
            file = new System.IO.StreamReader(DataFolder + InputFileName);
            if ((line = file.ReadLine()) != null)
            {
                string[] arr_line = line.Split(Program._ILShaped.delimiterChars);
                Program._ILShaped.q_Vector.AddRange(Array.ConvertAll(arr_line, double.Parse));
            }
            file.Close();

            // Read h_Vector
            Program._ILShaped.h_Vector.Clear();
            InputFileName = myScenario + "_h_Vector.txt";
            file = new System.IO.StreamReader(DataFolder + InputFileName);
            while ((line = file.ReadLine()) != null)
            {
                string[] arr_line = line.Split(Program._ILShaped.delimiterChars);
                Program._ILShaped.h_Vector.AddRange(Array.ConvertAll(arr_line, double.Parse));
            }
            file.Close();

            // Read R_Vector
            Program._ILShaped.R_Vector.Clear();
            InputFileName = myScenario + "_R_Vector.txt";
            file = new System.IO.StreamReader(DataFolder + InputFileName);
            while ((line = file.ReadLine()) != null)
            {
                string[] arr_line = line.Split(Program._ILShaped.delimiterChars);
                Program._ILShaped.R_Vector.AddRange(Array.ConvertAll(arr_line, double.Parse));
            }
            file.Close();
        }

        public static int ReadCSV(string filename, string NullName)
        {
            System.IO.StreamReader sr = new System.IO.StreamReader(filename);
            string read;
            int param = 0;

            if ((read = sr.ReadLine()) != null)
                param = Convert.ToInt32(read);
            return param;
        }

        public static List<double> ReadCSV(string filename)
        {
            System.IO.StreamReader sr = new System.IO.StreamReader(filename);
            string read;
            string[] arr_read;
            List<double> param = new List<double>();

            while ((read = sr.ReadLine()) != null)
            {
                arr_read = read.Split(',');
                param = Array.ConvertAll(arr_read, double.Parse).ToList();
            }
            return param;
        }

        public static List<string> ReadVar(string filename)
        {
            System.IO.StreamReader sr = new System.IO.StreamReader(filename);
            string read;
            string[] arr_read;
            List<string> param = new List<string>();

            while ((read = sr.ReadLine()) != null)
            {
                arr_read = read.Split(',');
                param = arr_read.ToList();
            }
            return param;
        }
        
        public static double[] ReadCSV(string filename, int numRows)
        {
            System.IO.StreamReader sr = new System.IO.StreamReader(filename);
            string read;
            string[] arr_read;
            double[] param = new double[numRows];

            while ((read = sr.ReadLine()) != null)
            {
                try
                {
                    arr_read = read.Split(',');
                    param = Array.ConvertAll(arr_read, double.Parse);
                }
                catch { continue; }
            }
            return param;
        }

        public static List<double[]> ReadCSV(string filename, int numVars, string NULL)
        {
            System.IO.StreamReader sr = new System.IO.StreamReader(filename);
            string read;
            string[] arr_read;
            List<double[]> list_param = new List<double[]>();
            double[] param = new double[numVars];

            while ((read = sr.ReadLine()) != null)
            {
                try
                {
                    arr_read = read.Split(',');
                    param = Array.ConvertAll(arr_read, double.Parse);
                    list_param.Add(param);
                }
                catch { continue; }
            }
            return list_param;
        }

        public static void WriteCSV(string filename, INumVar[] _INumVar)
        {
            System.IO.StreamWriter sw = new System.IO.StreamWriter(filename);
            
        }

        public static void WriteCSV(string filename, INumVar[][] _INumVar)
        {
            System.IO.StreamWriter sw = new System.IO.StreamWriter(filename);

        }
    }
}
