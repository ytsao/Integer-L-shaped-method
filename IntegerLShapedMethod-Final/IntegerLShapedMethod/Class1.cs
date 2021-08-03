using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace IntegerLShapedMethod
{
    class Class1
    {
        //_ILShapedMethod = new ILShapedMethod();
        //Data.getInputData_General(".//olddata/", _ILShapedMethod);

        //_ILShapedMethod.K = _ILShapedMethod.ScenarioProbability.Count();
        //_ILShapedMethod.numFirstStageVars = _ILShapedMethod.FirstStageVars.Count();
        //_ILShapedMethod.numSecondStageVars = _ILShapedMethod.SecondStageVars.Count();

        //try
        //{
        //    double EV_Function = 0;
        //    int numFeasibilityCuts = 0;
        //    int numOptimalityCuts = 0;
        //    double theta;
        //    int v = 1;
        //    _ILShapedMethod.BuildMasterProb();
        //    //Sampling
        //    //Sampling();
        //    // Only Build Master once at the beginning

        //    double w = 0;
        //    while (v <= 100)
        //    {
        //        //*****************************
        //        // Step1. Solve the Master Problem
        //        //*****************************
        //        _ILShapedMethod._FisrtStageSolver.MasterProb.SetOut(null);
        //        if (_ILShapedMethod._FisrtStageSolver.MasterProb.Solve())
        //            theta = _ILShapedMethod._FisrtStageSolver.MasterProb.GetValue(_ILShapedMethod._FisrtStageSolver.theta);
        //        else
        //        {
        //            System.Console.WriteLine("Infeasible Master Problem !");
        //            break;
        //        }

        //        double myObjValue = _ILShapedMethod._FisrtStageSolver.MasterProb.ObjValue - theta + w;
        //        System.Console.WriteLine(string.Format("Iteration {0}, Solution Status = {1}, Objective Value = {2}"
        //            , v, _ILShapedMethod._FisrtStageSolver.MasterProb.GetStatus(), myObjValue));
        //        _ILShapedMethod._FisrtStageSolver.MasterProb.ExportModel("Master_" + v + ".lp");
        //        _ILShapedMethod.x_v[0] = _ILShapedMethod._FisrtStageSolver.MasterProb.GetValues(_ILShapedMethod._FisrtStageSolver.myFirstStageVars);

        //        //***************************
        //        // Step2. Add Feasibility Cut
        //        //***************************
        //        bool isSecondStageProbFeasible = true;
        //        _ILShapedMethod._SecondStageSolver = new SecondStageSolver(_ILShapedMethod.numSecondStageVars, _ILShapedMethod.numSecondStageCons);
        //        for (int k = 1; k < _ILShapedMethod.K; k++)
        //        {
        //            double negative_d = 0;
        //            double[] negative_D = new double[_ILShapedMethod.numFirstStageVars + 1];      // D_vector size = number of first-stage variables + theta

        //            // Read General Input Data Files
        //            Data.getInputData_forScenario(".//data/", k.ToString(), _ILShapedMethod);

        //            _ILShapedMethod.numSecondStageCons = _ILShapedMethod.T_Matrix.Count();

        //             //Build SubProblem
        //            _ILShapedMethod._SecondStageSolver.Feasibility.ClearModel();
        //            _ILShapedMethod.BuildFeasibilityProb(_ILShapedMethod.x_v);
        //            _ILShapedMethod._SecondStageSolver.Feasibility.ExportModel("Feasiblity.lp");
        //            _ILShapedMethod._SecondStageSolver.Feasibility.SetOut(null);

        //            // Solve SubModel
        //            if (_ILShapedMethod._SecondStageSolver.Feasibility.Solve())
        //            {
        //                if (_ILShapedMethod._SecondStageSolver.Feasibility.ObjValue > 0) // Add Feasibility Cut, if obj > 0
        //                {
        //                    isSecondStageProbFeasible = false;

        //                    // Get Dual Values
        //                    double[] sigma = _ILShapedMethod._SecondStageSolver.Feasibility.GetDuals(_ILShapedMethod._SecondStageSolver.mySecondStageCons);

        //                    // Compute d
        //                    for (int i = 0; i < _ILShapedMethod._SecondStageSolver.mySecondStageCons.Count(); i++)
        //                        negative_d = negative_d - sigma[i] * _ILShapedMethod.h_Vector[i];

        //                    // Compute D
        //                    for (int j = 0; j < _ILShapedMethod.numFirstStageVars; j++)
        //                    {
        //                        for (int i = 0; i < _ILShapedMethod._SecondStageSolver.mySecondStageCons.Count(); i++)
        //                            negative_D[j] = negative_D[j] - sigma[i] * _ILShapedMethod.T_Matrix[i][j];
        //                    }

        //                    numFeasibilityCuts++;
        //                    _ILShapedMethod._FisrtStageSolver.MasterProb.AddGe(negative_d, _ILShapedMethod._FisrtStageSolver.MasterProb.ScalProd(negative_D
        //                        , _ILShapedMethod._FisrtStageSolver.myFirstStageVars), "FeaCut" + numFeasibilityCuts);
        //                }
        //            }
        //        }

        //        //****************************
        //        // Step3. Add Optimality Cut
        //        //****************************
        //        if (isSecondStageProbFeasible)
        //        {
        //            EV_Function = 0;
        //            double negative_es = 0;
        //            double[] negative_ES = new double[_ILShapedMethod.numFirstStageVars + 1]; // The ES_Vector size = num of first-stage variables + 1

        //            for (int k = 1; k <= _ILShapedMethod.K; k++)
        //            {
        //                // Read General input data files
        //                Data.getInputData_forScenario(".//olddata/", k.ToString(), _ILShapedMethod);

        //                _ILShapedMethod.numSecondStageCons = _ILShapedMethod.T_Matrix.Count();

        //                // Build SubProblem for generating Optimality Cut
        //                _ILShapedMethod._SecondStageSolver.Optimality.ClearModel();
        //                _ILShapedMethod.BuildOptimalityProb();
        //                _ILShapedMethod._SecondStageSolver.Optimality.SetOut(null);

        //                // Solve SubModel
        //                if (_ILShapedMethod._SecondStageSolver.Optimality.Solve())
        //                {
        //                    // Get Expectation of Objective Value
        //                    EV_Function = EV_Function + _ILShapedMethod.ScenarioProbability[k - 1] * _ILShapedMethod._SecondStageSolver.Optimality.ObjValue;

        //                    // Get Dual Values
        //                    double[] pi = _ILShapedMethod._SecondStageSolver.Optimality.GetDuals(_ILShapedMethod._SecondStageSolver.mySecondStageCons);

        //                    // Compute es
        //                    for (int i = 0; i < _ILShapedMethod._SecondStageSolver.mySecondStageCons.Count(); i++)
        //                        negative_es = negative_es - _ILShapedMethod.ScenarioProbability[k - 1] * pi[i] * _ILShapedMethod.h_Vector[i];

        //                    // Compute ES
        //                    for (int j = 0; j < _ILShapedMethod.numFirstStageVars; j++)
        //                    {
        //                        for (int i = 0; i < _ILShapedMethod._SecondStageSolver.mySecondStageCons.Count(); i++)
        //                            negative_ES[j] = negative_ES[j] - _ILShapedMethod.ScenarioProbability[k - 1] * pi[i] * _ILShapedMethod.T_Matrix[i][j];
        //                    }
        //                    _ILShapedMethod._SecondStageSolver.Optimality.ExportModel("Optimality_v" + v + "_k" + k + ".lp");
        //                }
        //            }

        //            // Compute w = e - Ex
        //            w = -negative_es;
        //            for (int j = 0; j < _ILShapedMethod.numFirstStageVars; j++)
        //                w = w + negative_ES[j] * _ILShapedMethod.x_v[0][j];

        //            if (w <= theta * (1 - ILShapedMethod.error))
        //                break;
        //            else
        //            {
        //                negative_ES[_ILShapedMethod.numFirstStageVars] = -1;
        //                numOptimalityCuts++;
        //                _ILShapedMethod._FisrtStageSolver.MasterProb.AddGe(negative_es
        //                    , _ILShapedMethod._FisrtStageSolver.MasterProb.ScalProd(negative_ES, _ILShapedMethod._FisrtStageSolver.myFirstStageVars), "OptCut" + numOptimalityCuts);
        //            }
        //        }
        //        v++;
        //    }
        //    System.Console.WriteLine("Finished!");
        //    System.Console.WriteLine(_ILShapedMethod._FisrtStageSolver.MasterProb.ObjValue);
        //    System.Console.ReadLine();

        //    _ILShapedMethod._FisrtStageSolver.MasterProb.End();
        //    //ILShapedMethod._Solver.Feasibility.End();
        //    _ILShapedMethod._SecondStageSolver.Optimality.End();
        //}
        //catch (ILOG.Concert.Exception ex)
        //{
        //    System.Console.WriteLine("Concert Error: " + ex);
        //    System.Console.ReadLine();
        //}
        //catch (Data.InputDataReaderException ex)
        //{
        //    System.Console.WriteLine("Data Error: " + ex);
        //    System.Console.ReadLine();
        //}
        //catch (System.IO.IOException ex)
        //{
        //    System.Console.WriteLine("IO Error: " + ex);
        //    System.Console.ReadLine();
        //}
    }
}
