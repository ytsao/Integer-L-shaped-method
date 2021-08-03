using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ILOG.Concert;
using ILOG.CPLEX;
using MathNet.Numerics.Integration.Algorithms;
using MathNet.Numerics;

namespace IntegerLShapedMethod
{
    public enum OPTIMALITYCUT { NORMALCUT, IMPROVEDCUT}

    public struct FisrtStageSolver
    {
        public Cplex MasterProb;
        public INumVar[] myFirstStageVars;
        public INumVar theta;
        
        public FisrtStageSolver(int numFirstStageVar)
        {
            MasterProb = new Cplex();

            myFirstStageVars = MasterProb.NumVarArray(numFirstStageVar, 0, 1, NumVarType.Bool);
            theta = MasterProb.NumVar(-ILShapedMethod.infinity, ILShapedMethod.infinity, NumVarType.Float,"theta");
        }
    }

    public struct SecondStageSolver
    {
        public Cplex Optimality;
        public Cplex Feasibility;
        public INumVar y_0;
        public INumVar R;
        public INumVar[] mySecondStageVars;
        public IRange[] mySecondStageCons;

        public SecondStageSolver(int numSecondStageVars, int numSecondStageCons)
        {
            Optimality = new Cplex();
            Feasibility = new Cplex();

            y_0 = Optimality.NumVar(-ILShapedMethod.infinity, ILShapedMethod.infinity, NumVarType.Float);
            R = Optimality.NumVar(0, ILShapedMethod.infinity, NumVarType.Float);
            mySecondStageVars = Optimality.NumVarArray(numSecondStageVars, 0, 1, NumVarType.Bool);
            mySecondStageCons = new IRange[numSecondStageCons];
        }
    }

    class ILShapedMethod : BranchAndBound
    {
        public const double error = 0.000000001;
        public const double infinity = 99;
        public char[] delimiterChars = { ' ', '\t' };

        public int w = 0;
        public string v = "0";
        public const int V = 100;
        public int numFirstStageVars = 0;
        public int numSecondStageVars = 0;
        public int numSecondStageCons = 0;
        int numRHS = 0;

        public double[][] x_v;
        public double theta_v;

        public double delta = 0;

        public List<double[]> A_Matrix;
        public List<double[]> T_Matrix;
        public List<double[]> W_Matrix;
        public List<double> c_Vector;
        public List<double> q_Vector;
        public List<double> b_Vector;
        public List<double> h_Vector;
        public List<double> R_Vector;

        public List<IAddable> OptCuts;

        public List<string> FirstStageVars;
        public List<string> SecondStageVars;
        public List<double> ScenarioProbability;

        public FisrtStageSolver _FisrtStageSolver;
        public SecondStageSolver _SecondStageSolver;

        public ILShapedMethod()
        {
            base.BinaryVarID = new List<int>();
            base.BinaryVarSln = new List<double>();
            base.myTree = new Dictionary<int, Dictionary<string, TreeNode>>();
            base.BCuts = new List<IAddable>();
            base._LBCandidate = new Dictionary<string, LBCandidate>();

            delimiterChars = new char[]{' ', '\t'};

            v = "0";
            numFirstStageVars = 0;
            numSecondStageVars = 0;
            numSecondStageCons = 0;

            x_v = new double[1][];
            theta_v = 0;

            A_Matrix = new List<double[]>();
            T_Matrix = new List<double[]>();
            W_Matrix = new List<double[]>();
            c_Vector = new List<double>();
            q_Vector = new List<double>();
            b_Vector = new List<double>();
            h_Vector = new List<double>();
            R_Vector = new List<double>();

            OptCuts = new List<IAddable>();

            FirstStageVars = new List<string>();
            SecondStageVars = new List<string>();
            ScenarioProbability = new List<double>();

            _FisrtStageSolver = new FisrtStageSolver(numFirstStageVars);                         // empty object
            _SecondStageSolver = new SecondStageSolver(numSecondStageVars, numSecondStageCons);  // empty object
        }

        public class CutGeneration
        {
            public static double qs = 0;
            public static double a = 0;

            public static void Compute_qs(ILShapedMethod _ILShapedMethod, int K)
            {
                _ILShapedMethod._SecondStageSolver.Optimality.SetOut(null);
                if (_ILShapedMethod._SecondStageSolver.Optimality.Solve())
                    qs = qs + _ILShapedMethod.ScenarioProbability[K] * _ILShapedMethod._SecondStageSolver.Optimality.ObjValue;
            }

            public static void Generate_OptimalityCut(ILShapedMethod _ILShapedMethod, string SelectedNodeID)
            {
                //_ILShapedMethod._FisrtStageSolver.MasterProb.Remove(Program._TreeNode.OptCuts.ToArray());
                //Program._TreeNode.OptCuts = new List<IRange>();

                double in_S = 0;
                double not_in_S = 0;
                ILinearNumExpr OptimalityCut = _ILShapedMethod._SecondStageSolver.Optimality.LinearNumExpr();
                for (int i = 0; i < _ILShapedMethod._FisrtStageSolver.myFirstStageVars.Length; i++)
                {
                    if (_ILShapedMethod._FisrtStageSolver.MasterProb.GetValue(_ILShapedMethod._FisrtStageSolver.myFirstStageVars[i]) == 1)
                    {
                        in_S++;
                        OptimalityCut.AddTerm((qs - _ILShapedMethod.BestLB), _ILShapedMethod._FisrtStageSolver.myFirstStageVars[i]);
                    }
                    else if (_ILShapedMethod._FisrtStageSolver.MasterProb.GetValue(_ILShapedMethod._FisrtStageSolver.myFirstStageVars[i]) == 0)
                    {
                        not_in_S++;
                        OptimalityCut.AddTerm(-(qs - _ILShapedMethod.BestLB), _ILShapedMethod._FisrtStageSolver.myFirstStageVars[i]);
                    }
                }
                _ILShapedMethod.delta = in_S - not_in_S;

                //double cut_RHS = (qs - _ILShapedMethod.BestLB) * _ILShapedMethod.delta - (qs - _ILShapedMethod.BestLB) * (in_S - 1) + _ILShapedMethod.BestLB;
                //IRange cut = _ILShapedMethod._FisrtStageSolver.MasterProb.AddGe(_ILShapedMethod._FisrtStageSolver.theta, cut_RHS, "OptCut");

                OptimalityCut.AddTerm(-1, _ILShapedMethod._FisrtStageSolver.theta);
                double cut_RHS = (qs-_ILShapedMethod.BestLB)*(in_S-1)-_ILShapedMethod.BestLB;
                IRange cut = _ILShapedMethod._FisrtStageSolver.MasterProb.AddGe(cut_RHS, OptimalityCut, "OptCut");
                _ILShapedMethod._FisrtStageSolver.MasterProb.AddCut(cut);
                Program._TreeNode.OptCuts.Add(cut);

                //string ParentNodeID = SelectedNodeID.Substring(0, SelectedNodeID.Length - 1);
                //if (ParentNodeID != "")
                //{
                //    TreeNode oldTreeNode = Program._ILShapedMethod.myTree[3][ParentNodeID];
                //    for (int numCuts = 0; numCuts < oldTreeNode.OptCuts.Count; numCuts++)
                //    {
                //        Program._TreeNode.OptCuts.Add(oldTreeNode.OptCuts[numCuts]);
                //        _ILShapedMethod._FisrtStageSolver.MasterProb.AddCut(oldTreeNode.OptCuts[numCuts]);
                //    }
                //}

                _ILShapedMethod._FisrtStageSolver.MasterProb.ExportModel("afterCut.lp");
            }

            public static double Lamda_s()
            { 
                double lamda = 0;

                return lamda;
            }

            public static void Generate_s_Neighbor(ILShapedMethod _ILShapedMethod)
            {
                double lamda = Lamda_s();
                // Proposition 6. Let xi = 1, i in S, xi = 0, i not in S be some solution with qa = E[Q(x)]
                // Let 1 <= t <= |S| be some integer, Then (2.3) holds with
                // a = max{max(qs - Lamda(s,S))/s; (qs-L)/(t+1)}
                a = Math.Max(0, (qs - _ILShapedMethod.BestLB)/2);
            }
        }
        
        public void BuildMasterProb()
        {
            numFirstStageVars = this.FirstStageVars.Count();
            this._FisrtStageSolver = new FisrtStageSolver(numFirstStageVars);
            numRHS = this.b_Vector.Count();
            
            //*************************************************************
            // Create Objective Function
            // The objective coefficient for teh first-stage variables
            // The ovjective coefficient for Theta
            //*************************************************************
            ILinearNumExpr Master_obj = this._FisrtStageSolver.MasterProb.LinearNumExpr();
            Master_obj.AddTerms(this.c_Vector.ToArray(), this._FisrtStageSolver.myFirstStageVars);
            Master_obj.AddTerm(1, this._FisrtStageSolver.theta);
            this._FisrtStageSolver.MasterProb.AddMinimize(Master_obj);

            // Create Constraints
            ILinearNumExpr Master_con = this._FisrtStageSolver.MasterProb.LinearNumExpr();
            for (int i = 0; i < numRHS; ++i)
            {
                Master_con.AddTerms(this.A_Matrix[i].ToArray(), this._FisrtStageSolver.myFirstStageVars);
                this._FisrtStageSolver.MasterProb.AddLe(Master_con, this.b_Vector[i]);
                Master_con.Clear();
            }

            this._FisrtStageSolver.MasterProb.ExportModel("Master_0.lp");
        }

        public void BuildFeasibilityProb(double[][] x)
        {
            int numVariables = this.numSecondStageVars + 2 * this.h_Vector.Count;

            //**************************
            // Create Objective Function
            // Coefficients : v+ and v-
            // int[] one = {1};
            //**************************
            int[] one_vector = new int[2 * this.numSecondStageCons];
            for (int i = 0; i < one_vector.Length; i++)
                one_vector[i] = 1;

            double[] q = new double[numVariables];
            Array.Copy(one_vector, 0, q, this.numSecondStageVars, 2 * this.numSecondStageCons);
            this._SecondStageSolver.Feasibility.AddMinimize(this._SecondStageSolver.Feasibility.ScalProd(q, this._SecondStageSolver.mySecondStageVars));

            // Create Constraints
            for (int i = 0; i < this.numSecondStageCons; i++)
            { 
                // RHS
                double myRHS = this.h_Vector[i];
                for (int j = 0; j < this.numFirstStageVars; j++)
                    myRHS = myRHS - this.T_Matrix[i][j] * x[0][j];
                
                // y coefficients
                double[] w = new double[numVariables];
                Array.Copy(this.W_Matrix[i], 0, w, 0, this.numSecondStageVars);

                // v+ and v- coefficients
                int[] one_n_negone = { 1, -1 };
                Array.Copy(one_n_negone, 0, w, this.numSecondStageVars + 2 * i, 2);
                this._SecondStageSolver.mySecondStageCons[i] = this._SecondStageSolver.Feasibility.AddGe(myRHS, this._SecondStageSolver.Feasibility.ScalProd(w, this._SecondStageSolver.mySecondStageVars));
            }
        }

        public void AddFeasibilityCut()
        {
            this.BuildFeasibilityProb(this.x_v);
        }

        public void BuildOptimalityProb()
        {
            // Create Objective Function
            ILinearNumExpr Optimality_obj = this._SecondStageSolver.Optimality.LinearNumExpr();
            Optimality_obj.AddTerm(this.q_Vector[0], this._SecondStageSolver.y_0);
            this._SecondStageSolver.Optimality.AddMinimize(Optimality_obj);

            this.x_v[0] = this._FisrtStageSolver.MasterProb.GetValues(this._FisrtStageSolver.myFirstStageVars);
            this.theta_v = this._FisrtStageSolver.MasterProb.GetValue(this._FisrtStageSolver.theta);

            // Create Constraints
            ILinearNumExpr Optimality_con2 = this._SecondStageSolver.Optimality.LinearNumExpr();
            for (int i = 0; i < this.numSecondStageCons; i++)
            {
                for (int j = 0; j < this.W_Matrix[i].Length; j++)
                {
                    if (j == 0) Optimality_con2.AddTerm(this.W_Matrix[i][j], this._SecondStageSolver.y_0);
                    else Optimality_con2.AddTerm(this.W_Matrix[i][j], this._SecondStageSolver.mySecondStageVars[j]);
                }
                

                double myRHS = this.h_Vector[i];
                for (int j = 0; j < this.numFirstStageVars; j++)
                    myRHS = myRHS + this.T_Matrix[i][j] * x_v[0][j];

                if (i == 0)
                {
                    Optimality_con2.AddTerm(R_Vector[i], _SecondStageSolver.R);
                    this._SecondStageSolver.mySecondStageCons[i] = this._SecondStageSolver.Optimality.AddEq(Optimality_con2, myRHS);
                }
                else
                {
                    Optimality_con2.AddTerm(R_Vector[i], _SecondStageSolver.R);
                    this._SecondStageSolver.mySecondStageCons[i] = this._SecondStageSolver.Optimality.AddLe(Optimality_con2, myRHS);
                }
                Optimality_con2.Clear();
            }
        }

        public void AddOptimalityCut(OPTIMALITYCUT _OPTIMALITYCUT)
        {
            switch (_OPTIMALITYCUT)
            { 
                case OPTIMALITYCUT.NORMALCUT:
                    CutGeneration.Compute_qs(Program._ILShapedMethod, 0);
                    break;
                case OPTIMALITYCUT.IMPROVEDCUT:
                    CutGeneration.Generate_s_Neighbor(Program._ILShapedMethod);
                    break;
            }  
        }

        public override void CreateRootNode(TreeNode _TreeNode)
        {
            BuildMasterProb();

            this._FisrtStageSolver.MasterProb.ImportModel(@"Master_0.lp");
            base.mEnum = this._FisrtStageSolver.MasterProb.GetLPMatrixEnumerator();
            base.mEnum.MoveNext();
            base.lpmatrix = (ILPMatrix)base.mEnum.Current;
            base.numRows = base.lpmatrix.Nrows;
            base.numCols = base.lpmatrix.Ncols;

            for (int j = 0; j < base.lpmatrix.NumVars.Length; ++j)
            {
                //String temp = ObjElement[2*j+2]; // Warning Customization Setting Param of Ceof
                //this.Variables.Add(Tuple.Create(j, Convert.ToDouble(temp), this.lpmatrix.GetNumVar(j).Type));
                if (base.lpmatrix.GetNumVar(j).Type.ToString() == "Bool")
                {
                    base.BinaryVarID.Add(j);
                    IConversion LPrelax = this._FisrtStageSolver.MasterProb.Conversion(base.lpmatrix.GetNumVar(j), NumVarType.Float);
                    this._FisrtStageSolver.MasterProb.Add(LPrelax);
                }
            }

            base.mEnum = this._FisrtStageSolver.MasterProb.GetLPMatrixEnumerator();
            base.mEnum.MoveNext();
            base.lpmatrix = (ILPMatrix)base.mEnum.Current;

            this._FisrtStageSolver.myFirstStageVars = new INumVar[numFirstStageVars];
            for (int i = 0; i < this.numCols; i++)
            {
                if (i == 2)
                    this._FisrtStageSolver.theta = base.lpmatrix.GetNumVar(i);
                else
                    this._FisrtStageSolver.myFirstStageVars[i] = base.lpmatrix.GetNumVar(i);
            }

            this._FisrtStageSolver.MasterProb.ExportModel("Root_Node.lp");
            base.myTree[0].Add("0", _TreeNode);
            
            _TreeNode.NodeStack.Push("0");
            _TreeNode.NodeQueue.Enqueue("0");
        }

        public override void CreateChildNode(int BranchingVar, string SelectedNodeID, BRANCHING _Branching){base.CreateChildNode(BranchingVar, SelectedNodeID, _Branching);}
       
        public override void CreateNode(string SelectedNodeID)
        { 
            this._FisrtStageSolver.MasterProb.Remove(base.BCuts.ToArray());
            base.BCuts = new List<IAddable>();

            for (int i = 0; i < base.myTree[0][SelectedNodeID].ParentID.Count; ++i)
            {
                if (base.myTree[0][SelectedNodeID].BCutsRHS[i] == 0)
                {
                    IRange Cut = this._FisrtStageSolver.MasterProb.AddLe(this._FisrtStageSolver.MasterProb.Prod(1.0, getLPNumVar(base.myTree[0][SelectedNodeID].SelectedVarID[i])), 0);
                    base.BCuts.Add(this._FisrtStageSolver.MasterProb.AddCut(Cut));
                }
                else if (this.myTree[0][SelectedNodeID].BCutsRHS[i] == 1)
                {
                    IRange Cut = this._FisrtStageSolver.MasterProb.AddGe(this._FisrtStageSolver.MasterProb.Prod(1.0, getLPNumVar(base.myTree[0][SelectedNodeID].SelectedVarID[i])), 1);
                    base.BCuts.Add(this._FisrtStageSolver.MasterProb.AddCut(Cut));
                }
            }

            this._FisrtStageSolver.MasterProb.ExportModel(string.Format("Node{0}.lp",SelectedNodeID));
        }
        
        public override string NodeSelection(TreeNode _TreeNode, NODESELECTIONTYPE _SELECTIONTYPE){return base.NodeSelection(_TreeNode, _SELECTIONTYPE);}

        public override int VariableSelection(VARIABLESELECTION _VARIABLESELECTION) { return base.VariableSelection(VARIABLESELECTION.MOSTFRAC); }

        public override void SolveNode(string SelectedNodeID, int algorithm)
        {
            this._FisrtStageSolver.MasterProb.SetOut(null);
            this._FisrtStageSolver.MasterProb.SetParam(Cplex.Param.RandomSeed, 0);
            this._FisrtStageSolver.MasterProb.SetParam(Cplex.Param.RootAlgorithm, algorithm);
            SetAlgorithm.Preprocessing_off(this._FisrtStageSolver.MasterProb);
            SetAlgorithm.Cuts_off(this._FisrtStageSolver.MasterProb);

            _FisrtStageSolver.MasterProb.Solve();
            double myObjValue = this._FisrtStageSolver.MasterProb.ObjValue;

            //System.Console.WriteLine(string.Format("Iteration {0}, Solution Status = {1}, Objective Value = {2}"
            //    , this.v, this._FisrtStageSolver.MasterProb.GetStatus(), myObjValue));
            //System.Console.WriteLine();
            
            //this._FisrtStageSolver.MasterProb.ExportModel("Master_" + this.v + ".lp");

            this.BinaryVarSln = new List<double>();
            foreach (int idx in this.BinaryVarID)
                this.BinaryVarSln.Add(this._FisrtStageSolver.MasterProb.GetValue(this._FisrtStageSolver.myFirstStageVars[idx]));

            TreeNode _TreeNode = this.myTree[0][SelectedNodeID];
            _TreeNode.LPR = myObjValue;
            _TreeNode.Theta = this._FisrtStageSolver.MasterProb.GetValue(this._FisrtStageSolver.theta);
            _TreeNode.BinSln = this.BinaryVarSln;

            this.myTree[0][SelectedNodeID] = _TreeNode;
            System.Console.WriteLine(string.Format("{0} objective value: {1}", SelectedNodeID, _TreeNode.LPR));
        }

        public override bool CheckIPFeasible()
        {
            double[] solutions = this._FisrtStageSolver.MasterProb.GetValues(this._FisrtStageSolver.myFirstStageVars);
            double solution = this._FisrtStageSolver.MasterProb.GetValue(this._FisrtStageSolver.theta);

            //Check if all binary variables are integer
            this.isIPFeasible = Array.TrueForAll(solutions, value => { return Math.Abs(value - Math.Round(value)) < 0.0000000001; });
            //if ((solution - Math.Round(solution) < 0.0000000001) && (this.isIPFeasible))
            //    this.isIPFeasible = true;
            //else this.isIPFeasible = false;
            return this.isIPFeasible;
        }

        public override void IdentifyNodeStatus(string SelectedNodeID)
        {
            CheckIPFeasible();

            //NodeStatus (1:Best IP, 2:Pruned IP, 3:Best LP / Regular LP, 4: Pruned LP, 5:LP Infeasible)
            TreeNode _TreeNode = base.myTree[0][SelectedNodeID];
            UpdateLB(_TreeNode, SelectedNodeID);

            if (base.isIPFeasible)
            {
                if (_TreeNode.LPR <= base.BestUB)
                {
                    base.NodeStatus = 1;

                    //Update Upper Bound & BestNodeID
                    base.BestUB = _TreeNode.LPR;
                    base.BestNodeID = SelectedNodeID;
                }
                else base.NodeStatus = 2;
            }
            else if (_TreeNode.LPR <= base.BestUB)
            {
                if (_TreeNode.LPR >= base.BestLB)
                    base.NodeStatus = 3;

                //Variable Selection and Create two Child Node
                BranchingID = VariableSelection(VARIABLESELECTION.MOSTFRAC);

                CreateChildNode(BranchingID, SelectedNodeID, BRANCHING.RIGHT);
                CreateChildNode(BranchingID, SelectedNodeID, BRANCHING.LEFT);
            }
            else base.NodeStatus = 4;

            this.BestLB = _LBCandidate.First().Value.LPR;
            base.myTree[NodeStatus].Add(SelectedNodeID, _TreeNode);
            base.myTree[0].Remove(SelectedNodeID);          // Remove incumbent node
        }
        
        public override void UpdateLB(TreeNode _TreeNode, string SelectedNodeID){ base.UpdateLB(_TreeNode, SelectedNodeID);}
       
        public override void UpdateLB(string SelectedNodeID){ base.UpdateLB(SelectedNodeID);}
        
        public override INumVar getLPNumVar(int idx){ return base.getLPNumVar(idx);}
    }
}
