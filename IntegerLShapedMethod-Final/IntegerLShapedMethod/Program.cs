using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.Integration;

namespace IntegerLShapedMethod
{
    class Program
    {
        #region Declare global variables
        public static BranchAndBound BB;
        public static DisjunctiveCuts _DisjunctiveCuts;
        public static TreeNode _TreeNode;
        public static ILShaped _ILShaped;
        #endregion

        [STAThread]
        static void Main(string[] args)
        {
            #region Choose Algorithm
            //MainBranchAndBound();
            GomoryCutsPaper();
            #endregion
        }

        public static void MainBranchAndBound()
        {
            #region Initialize 
            System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();

            int algorithm = 0;
            sw.Reset();
            sw.Start();
            int iteration = 1;
            BB = new BranchAndBound("ISVM_15data.lp");
            _DisjunctiveCuts = new DisjunctiveCuts();
            _TreeNode = new TreeNode(BB);

            BB.CreateRootNode(_TreeNode);
            string SelectedNodeID = "0";
            #endregion

            #region  Termination condition ex.Node set is empty
            do
            {
                try
                {
                    SelectedNodeID = BB.NodeSelection(_TreeNode, NODESELECTIONTYPE.DFS);
                    if (SelectedNodeID == "-1") break;
                    BB.CreateNode(SelectedNodeID);
                    BB.SolveNode(SelectedNodeID, algorithm);
                    BB.CheckIPFeasible();
                    BB.IdentifyNodeStatus(SelectedNodeID, VARIABLESELECTION.MOSTFRAC);
                    if (BB.NodeStatus == 3 || BB.NodeStatus == 5)
                        _DisjunctiveCuts.D_CutGeneration(BB.BranchingID, SelectedNodeID, BB.NodeStatus);
                }
                catch (ILOG.Concert.Exception e)
                {
                    BB.UpdateLB(SelectedNodeID);
                    int status = 5;
                    BB.myTree[status].Add(SelectedNodeID, _TreeNode);
                    BB.myTree[0].Remove(SelectedNodeID);
                    System.Console.WriteLine(string.Format("Node {0} No Solution", SelectedNodeID));
                }
                if (BB._LBCandidate.First().Value.LPR == BB.BestUB) break;
                iteration++;
            }//while (iteration <= 10);
            while (BB.myTree[0].Count != 0);
            //while (sw.Elapsed.TotalMinutes <= 100);
            sw.Stop();
            #endregion

            #region Show results
            System.Console.WriteLine("*******************************" + sw.Elapsed.TotalSeconds);
            System.Console.WriteLine("Best LB = " + BB.BestLB);
            foreach (var i in BB.myTree[1])
                System.Console.WriteLine("Best IP = " + i.Value.LPR);

            System.Console.ReadLine(); 
            #endregion
        }

        /// <summary>
        /// Integer L-Shaped Method (Branching + Bender's)
        /// </summary>
        public static void GomoryCutsPaper()
        {
            #region Initialize (Read data set & Create root node & set off preprocessing and cuts of CPLE)
            int algorithm = 0;
            string SelectedNodeID;
            double Best_theta = 0;

            _ILShaped = new ILShaped();
            Data.getInputData_General(".//GomoryCutsPaper_data/");
            ILShaped.Model.BuildMasterProb();
            _ILShaped = new ILShaped(_ILShaped, "Master_0.lp");

            int numScenarioProbability = _ILShaped.ScenarioProbability.Count;
            _TreeNode = new TreeNode(_ILShaped);

            _ILShaped.node_id = "0";
            _ILShaped.BestUB = 99999999;
            _ILShaped.BestLB = -99999999;

            _ILShaped.CreateRootNode(_TreeNode);

            SetAlgorithm.Cuts_off(_ILShaped.RootCplex);
            SetAlgorithm.Preprocessing_off(_ILShaped.RootCplex);

            SetAlgorithm.Cuts_off(_ILShaped._SecStageProb.OptimalityProb);
            SetAlgorithm.Preprocessing_off(_ILShaped._SecStageProb.OptimalityProb);
            #endregion

            do
            {
                #region Step 1. Node Selection & Create it
                SelectedNodeID = _ILShaped.NodeSelection(_TreeNode, NODESELECTIONTYPE.BESTBOUND);
                if (SelectedNodeID == "-1") break;
                _ILShaped.CreateNode(SelectedNodeID);
                #endregion

                try
                {
                    #region Step 2. Solve current problem & get solution
                    _ILShaped.SolveNode(SelectedNodeID, algorithm);

                    _ILShaped.x_v_Sln[0] = _ILShaped.RootCplex.GetValues(_ILShaped._INumVar);
                    _ILShaped.theta_v_Sln = _ILShaped.RootCplex.GetValue(_ILShaped.theta);
                    #endregion

                    #region Step 3. & 4. if current solution > UB then fathom the current node [Check for integrality restrictions (Create Child Node) => IdentifyNodeStatus]
                    _ILShaped.IdentifyNodeStatus(SelectedNodeID, VARIABLESELECTION.MOSTFRAC);     // Including Create Child Node
                    #endregion

                    #region Step 5. Compute qs = E[Q(x_v)] and Update LB & UB
                    while (_ILShaped.NodeStatus == 1)
                    {
                        #region Compute Q(x_v) = qs
                        ILShaped.OptCutGeneration.qs = 0;
                        for (int k = 1; k <= numScenarioProbability; k++)
                        {
                            Data.getInputData_forScenario(".//GomoryCutsPaper_data/", k.ToString());
                            _ILShaped.numSecStageVars = _ILShaped.SecStageVars.Count;
                            _ILShaped.numSecStageCons = _ILShaped.W_Matrix.Count;
                            _ILShaped._SecStageProb = new SecStageProb(_ILShaped.numSecStageVars, _ILShaped.numSecStageCons);
                            ILShaped.Model.BuildOptimalityProb(_ILShaped);
                            ILShaped.Model.Compute_param(k - 1, _ILShaped, SelectedNodeID, ".//GomoryCutsPaper_data/", OPTIMALITYCUTS.LAPROT);
                        }
                        #endregion

                        #region Compute z_v = cTx_v + Q(x_v)
                        double z_v = 0;
                        for (int i = 0; i < _ILShaped.c_Vector.Count; i++) z_v = z_v + _ILShaped.c_Vector[i] * _ILShaped.x_v_Sln[0][i];
                        z_v = z_v + ILShaped.OptCutGeneration.qs;

                        // Update  UB
                        if (z_v < _ILShaped.BestUB) _ILShaped.BestUB = z_v;
                        #endregion

                        #region Step 6. Impose Optimaliyt Cut
                        if (Math.Abs(ILShaped.OptCutGeneration.qs - _ILShaped.theta_v_Sln) <= 0.000000001)
                        {
                            TreeNode _oldTreeNode = _ILShaped.myTree[0][SelectedNodeID];
                            _ILShaped.myTree[1].Add(SelectedNodeID, _oldTreeNode);
                            _ILShaped.myTree[0].Remove(SelectedNodeID);
                            _ILShaped.UpdateLB(SelectedNodeID);
                            Best_theta = _ILShaped.theta_v_Sln;
                            break;
                        }
                        else if (ILShaped.OptCutGeneration.qs < _ILShaped.theta_v_Sln)
                        {
                            TreeNode _oldTreeNode = _ILShaped.myTree[0][SelectedNodeID];
                            _ILShaped.myTree[4].Add(SelectedNodeID, _oldTreeNode);
                            _ILShaped.myTree[0].Remove(SelectedNodeID);
                            _ILShaped.UpdateLB(SelectedNodeID);
                            break;
                        }
                        //else ILShaped.OptCutGeneration.Generate_OptimalityCut(_ILShaped, SelectedNodeID, ADDITIONALCUTTYPE.LAPROTCUT);
                        else if (SelectedNodeID.Length - 1 != _ILShaped._INumVar.Length)
                        {
                            ILShaped.OptCutGeneration.Lambda_s(_ILShaped, SelectedNodeID, ".//GomoryCutsPaper_data/", ADDITIONALCUTTYPE.PROPOSITION5);
                            ILShaped.OptCutGeneration.Generate_OptimalityCut(_ILShaped, SelectedNodeID, ADDITIONALCUTTYPE.PROPOSITION5);
                        }
                        else
                        {
                            TreeNode _oldTreeNode = _ILShaped.myTree[0][SelectedNodeID];
                            _ILShaped.myTree[1].Add(SelectedNodeID, _oldTreeNode);
                            _ILShaped.myTree[0].Remove(SelectedNodeID);
                            _ILShaped._LBCandidate.Remove(SelectedNodeID);
                            break;
                        }

                        //else if (SelectedNodeID.Length - 1 != _ILShaped._INumVar.Length)
                        //{
                        //    ILShaped.CutGeneration.Lambda_s(_ILShaped, SelectedNodeID, ADDITIONALCUTTYPE.PROPOSITION6);
                        //    ILShaped.CutGeneration.Generate_OptimalityCut(_ILShaped, SelectedNodeID, ADDITIONALCUTTYPE.PROPOSITION6);
                        //}
                        //else
                        //{
                        //    TreeNode _oldTreeNode = _ILShaped.myTree[0][SelectedNodeID];
                        //    _ILShaped.myTree[1].Add(SelectedNodeID, _oldTreeNode);
                        //    _ILShaped.myTree[0].Remove(SelectedNodeID);
                        //    _ILShaped._LBCandidate.Remove(SelectedNodeID);
                        //    break;
                        //}

                        // Solve current node again
                        _ILShaped.CreateNode(SelectedNodeID);
                        _ILShaped.SolveNode(SelectedNodeID, algorithm);
                        _ILShaped.x_v_Sln[0] = _ILShaped.RootCplex.GetValues(_ILShaped._INumVar);
                        _ILShaped.theta_v_Sln = _ILShaped.RootCplex.GetValue(_ILShaped.theta);

                        // Check Node Status again
                        _ILShaped.IdentifyNodeStatus(SelectedNodeID, VARIABLESELECTION.MOSTFRAC);     // Including Create Child Node
                        #endregion
                    }
                    #endregion
                }
                catch (Exception error)
                {
                    #region Infeasible node 
                    _ILShaped.UpdateLB(SelectedNodeID);
                    int status = 5;
                    _ILShaped.myTree[status].Add(SelectedNodeID, _TreeNode);
                    _ILShaped.myTree[0].Remove(SelectedNodeID);
                    System.Console.WriteLine(string.Format("{0} No Solution", SelectedNodeID));

                    _ILShaped._LBCandidate.Remove(SelectedNodeID);
                    #endregion
                }
            } while (_ILShaped.myTree[0].Count != 0);

            #region Show the Results
            System.Console.WriteLine("Best LB :" + _ILShaped._LBCandidate.First().Value.LPR);
            System.Console.WriteLine("Best UB : " + _ILShaped.BestUB);

            System.Console.WriteLine("\nFirst Stage ObjValue: " + _ILShaped.BestUB);
            System.Console.WriteLine("Seconde Stage ObjValue: " + Best_theta);

            System.Console.WriteLine("\nBest Node ID: " + _ILShaped._LBCandidate.First().Key);
            System.Console.WriteLine("Best theta Sln: " + Best_theta);
            System.Console.ReadLine();
            #endregion
        }

        public static void Scenario_5_25_50()
        {
            //*********************************
            // Integer L-Shaped Method
            // Branching + Bender's 
            //
            // Step 0. Initialization
            // Create RootNode
            //*********************************
            _ILShaped = new ILShaped();
            Data.getInputData_GenearalCSV(".//5_25_50data/");
            ILShaped.Model.BuildSIPLIBMaster(_ILShaped);
            _ILShaped = new ILShaped(_ILShaped, "Master_0.lp");

            int numScenarioProbability = _ILShaped.ScenarioProbability.Count;
            _TreeNode = new TreeNode(_ILShaped);

            _ILShaped.node_id = "0";
            _ILShaped.BestUB = 9999;
            _ILShaped.BestLB = -99999999;

            _ILShaped.CreateRootNode(_TreeNode);

            int algorithm = 0;
            string SelectedNodeID;

            double Best_theta = 0;

            SetAlgorithm.Cuts_off(_ILShaped.RootCplex);
            SetAlgorithm.Preprocessing_off(_ILShaped.RootCplex);

            SetAlgorithm.Cuts_off(_ILShaped._SecStageProb.OptimalityProb);
            SetAlgorithm.Preprocessing_off(_ILShaped._SecStageProb.OptimalityProb);

            //*****************************************************************************************
            // Step 1. Node Selection
            //      SELECTIONTYPE : DFS, BFS, BESTBOUND
            //
            // Step 2. Solve Current Problem
            //
            // Step 3 & 4. if current solution > UB then fathom the current node
            //             Check for integrality restrictions (Create Child Node) => IdentifyNodeStatus
            //
            // Step 5. Compute qs = E[Q(x_v)] and Update LB & UB
            //
            // Step 6. Add Optimality Cuts (Improved Optimality Cuts) 
            //*****************************************************************************************
            do
            {
                // Step 1.
                SelectedNodeID = _ILShaped.NodeSelection(_TreeNode, NODESELECTIONTYPE.BESTBOUND);
                if (SelectedNodeID == "-1") break;
                _ILShaped.CreateNode(SelectedNodeID);

                try
                {
                    // Step 2. Solve current problem
                    _ILShaped.SolveNode(SelectedNodeID, algorithm);

                    _ILShaped.x_v_Sln[0] = _ILShaped.RootCplex.GetValues(_ILShaped._INumVar);
                    _ILShaped.theta_v_Sln = _ILShaped.RootCplex.GetValue(_ILShaped.theta);

                    // Step 3. & 4.
                    _ILShaped.IdentifyNodeStatus(SelectedNodeID, VARIABLESELECTION.MOSTFRAC);     // Including Create Child Node

                    // Step 5.
                    while (_ILShaped.NodeStatus == 1)
                    {
                        // Compute Q(x_v) = qs
                        ILShaped.OptCutGeneration.qs = 0;
                        for (int k = 1; k <= numScenarioProbability; k++)
                        {
                            Data.getInputData_forScenarioCSV(".//5_25_50data/", k.ToString());
                            _ILShaped.numSecStageVars = _ILShaped.SecStageVars.Count;
                            _ILShaped.numSecStageCons = _ILShaped.W_Matrix.Count;
                            _ILShaped._SecStageProb = new SecStageProb(_ILShaped.numSecStageCons, _ILShaped.numFirStageVars, _ILShaped.numSecStageCons);
                            ILShaped.Model.BuildSIPLIBSubProb(_ILShaped);
                            ILShaped.Model.Compute_param(k - 1, _ILShaped, SelectedNodeID, ".//5_25_50data/", OPTIMALITYCUTS.LAPROT);
                        }

                        //Compute z_v = cTx_v + Q(x_v)
                        double z_v = 0;
                        for (int i = 0; i < _ILShaped.c_Vector.Count; i++) z_v = z_v + _ILShaped.c_Vector[i] * _ILShaped.x_v_Sln[0][i];
                        z_v = z_v + ILShaped.OptCutGeneration.qs;

                        // Update  UB
                        if (z_v < _ILShaped.BestUB) _ILShaped.BestUB = z_v;

                        // Step 6. Impose Optimaliyt Cut
                        if (ILShaped.OptCutGeneration.qs == _ILShaped.theta_v_Sln)
                        {
                            TreeNode _oldTreeNode = _ILShaped.myTree[0][SelectedNodeID];
                            _ILShaped.myTree[1].Add(SelectedNodeID, _oldTreeNode);
                            _ILShaped.myTree[0].Remove(SelectedNodeID);
                            _ILShaped.UpdateLB(SelectedNodeID);
                            Best_theta = _ILShaped.theta_v_Sln;
                            _ILShaped.BestNodeID = SelectedNodeID;
                            break;
                        }
                        else if (ILShaped.OptCutGeneration.qs < _ILShaped.theta_v_Sln)
                        {
                            TreeNode _oldTreeNode = _ILShaped.myTree[0][SelectedNodeID];
                            _ILShaped.myTree[4].Add(SelectedNodeID, _oldTreeNode);
                            _ILShaped.myTree[0].Remove(SelectedNodeID);
                            _ILShaped.UpdateLB(SelectedNodeID);
                            break;
                        }
                        //else ILShaped.OptCutGeneration.Generate_OptimalityCut(_ILShaped, SelectedNodeID, ADDITIONALCUTTYPE.LAPROTCUT);
                        else if (SelectedNodeID.Length - 1 != _ILShaped._INumVar.Length)
                        {
                            ILShaped.OptCutGeneration.Lambda_s(_ILShaped, SelectedNodeID, ".//5_25_50data/", ADDITIONALCUTTYPE.PROPOSITION5);
                            ILShaped.OptCutGeneration.Generate_OptimalityCut(_ILShaped, SelectedNodeID, ADDITIONALCUTTYPE.PROPOSITION5);
                        }
                        else
                        {
                            TreeNode _oldTreeNode = _ILShaped.myTree[0][SelectedNodeID];
                            _ILShaped.myTree[1].Add(SelectedNodeID, _oldTreeNode);
                            _ILShaped.myTree[0].Remove(SelectedNodeID);
                            _ILShaped._LBCandidate.Remove(SelectedNodeID);
                            break;
                        }
                        //else if (SelectedNodeID.Length - 1 != _ILShaped._INumVar.Length)
                        //{
                        //    ILShaped.CutGeneration.Lambda_s(_ILShaped, SelectedNodeID, ADDITIONALCUTTYPE.PROPOSITION6);
                        //    ILShaped.CutGeneration.Generate_OptimalityCut(_ILShaped, SelectedNodeID, ADDITIONALCUTTYPE.PROPOSITION6);
                        //}
                        //else
                        //{
                        //    TreeNode _oldTreeNode = _ILShaped.myTree[0][SelectedNodeID];
                        //    _ILShaped.myTree[1].Add(SelectedNodeID, _oldTreeNode);
                        //    _ILShaped.myTree[0].Remove(SelectedNodeID);
                        //    _ILShaped._LBCandidate.Remove(SelectedNodeID);
                        //    break;
                        //}


                        // Solve current node again
                        _ILShaped.CreateNode(SelectedNodeID);
                        _ILShaped.SolveNode(SelectedNodeID, algorithm);
                        _ILShaped.x_v_Sln[0] = _ILShaped.RootCplex.GetValues(_ILShaped._INumVar);
                        _ILShaped.theta_v_Sln = _ILShaped.RootCplex.GetValue(_ILShaped.theta);

                        // Check Node Status again
                        _ILShaped.IdentifyNodeStatus(SelectedNodeID, VARIABLESELECTION.MOSTFRAC);     // Including Create Child Node
                    }
                }
                catch (Exception error)
                {
                    _ILShaped.UpdateLB(SelectedNodeID);
                    int status = 5;
                    _ILShaped.myTree[status].Add(SelectedNodeID, _TreeNode);
                    _ILShaped.myTree[0].Remove(SelectedNodeID);
                    System.Console.WriteLine(string.Format("{0} No Solution", SelectedNodeID));

                    _ILShaped._LBCandidate.Remove(SelectedNodeID);
                }
            } while (_ILShaped.myTree[0].Count != 0);


            System.Console.WriteLine("Best LB :" + _ILShaped._LBCandidate.First().Value.LPR);
            System.Console.WriteLine("Best UB : " + _ILShaped.BestUB);

            System.Console.WriteLine("\nFirst Stage ObjValue: " + _ILShaped.BestUB);
            System.Console.WriteLine("Seconde Stage ObjValue: " + Best_theta);

            System.Console.WriteLine("\nBest Node ID: " + _ILShaped._LBCandidate.First().Key);
            System.Console.WriteLine("Best theta Sln: " + Best_theta);
            System.Console.ReadLine();
        }
    }
}
