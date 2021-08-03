using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;
using ILOG.Concert;
using ILOG.CPLEX;

namespace IntegerLShapedMethod
{
    #region Option for optimality cuts & addtional cuts
    enum OPTIMALITYCUTS { LAPROT, SNEIGHBOR}

    enum ADDITIONALCUTTYPE { LAPROTCUT, PROPOSITION5, PROPOSITION6}
    #endregion

    struct SecStageProb
    {
        #region Declare variables
        public Cplex OptimalityProb;
        public INumVar y_0;
        public INumVar R;
        public INumVar[] mySecStageVars;
        public IRange [] mySecStageCons;

        // SIPLIB instance
        public INumVar[][] y;
        public INumVar[] y_j_0;
        #endregion

        #region Constructor (int, int)
        public SecStageProb(int numSecStageVars, int numSecStageCons)
        {
            OptimalityProb = new Cplex();

            y_0 = OptimalityProb.NumVar(-ILShaped.infinity, ILShaped.infinity, NumVarType.Int, "y_0");
            R = OptimalityProb.NumVar(0, ILShaped.infinity, NumVarType.Int, "R");
            mySecStageVars = new INumVar[numSecStageVars];
            for (int j = 0; j < numSecStageVars; j++) mySecStageVars[j] = OptimalityProb.NumVar(0, ILShaped.infinity, NumVarType.Int, string.Format("y_{0}", j));
            mySecStageCons = new IRange[numSecStageCons];

            y = new INumVar[numSecStageVars][];
            y_j_0 = new INumVar[numSecStageVars];
        }
        #endregion

        // for SIPLIB instance (SSLP)
        #region Constructor (int, int, int) 
        public SecStageProb(int Var_i, int Var_j, int numSecStageCons)
        {
            OptimalityProb = new Cplex();

            y_0 = OptimalityProb.NumVar(-ILShaped.infinity, ILShaped.infinity, NumVarType.Int, "y_0");
            R = OptimalityProb.NumVar(0, ILShaped.infinity, NumVarType.Int, "R");
            mySecStageVars = new INumVar[Var_j];
            mySecStageCons = new IRange[numSecStageCons];

            y = new INumVar[Var_i][];
            for (int i = 0; i < Var_i; i++)
            {
                y[i] = new INumVar[Var_j];
                for (int j = 0; j < Var_j; j++) y[i][j] = OptimalityProb.NumVar(0, 1, NumVarType.Bool, string.Format("y_{0}_{1}", i, j));
            }

            y_j_0 = new INumVar[Var_j];
            for (int j = 0; j < Var_j; j++) y_j_0[j] = OptimalityProb.NumVar(0, double.MaxValue, NumVarType.Int, string.Format("y_{0}_0", j));
        }
        #endregion
    }

    class ILShaped:BranchAndBound
    {
        #region 
        public const double error = 0.0000000001;
        public const double infinity = 99;
        public char[] delimiterChars = { ' ', '\t' };

        public INumVar theta;

        public List<double[]> A_Matrix;
        public List<double[]> T_Matrix;
        public List<double[]> W_Matrix;
        public List<double> c_Vector;
        public List<double> q_Vector;
        public List<double> b_Vector;
        public List<double> h_Vector;
        public List<double> R_Vector;

        public List<double[]> d_Matrix;
        public List<double[]> q_Matrix;

        public List<IRange> OptCuts;

        public List<string> FirStageVars;
        public List<string> SecStageVars;
        public List<double> ScenarioProbability;

        public int w = 0;
        public string node_id = "0";
        public int numFirStageVars = 0;
        public int numSecStageVars = 0;
        public int numSecStageCons = 0;
        public int numFirStageRHS = 0;
        public int numSecStageRHS = 0;

        public double[][] x_v_Sln;
        public double theta_v_Sln;
        #endregion

        public SecStageProb _SecStageProb;

        #region Constructor ()
        public ILShaped()
        {
            // 1.
            #region 1. Initialize Partial Base Class
            base.RootCplex = new Cplex();
            theta = base.RootCplex.NumVar(-ILShaped.infinity, ILShaped.infinity, NumVarType.Float, "theta");

            base.BCuts = new List<IAddable>();

            base.Variables = new List<Tuple<int, double, NumVarType>>();
            base.BinaryVarID = new List<int>();
            base.BinaryVarSln = new List<double>();
            base.myTree = new Dictionary<int, Dictionary<string, TreeNode>>();
            base._LBCandidate = new Dictionary<string, LBCandidate>();
            #endregion

            // 2.
            #region Initialize ILShaped Class Variables
            A_Matrix = new List<double[]>();
            T_Matrix = new List<double[]>();
            W_Matrix = new List<double[]>();
            c_Vector = new List<double>();
            q_Vector = new List<double>();
            b_Vector = new List<double>();
            h_Vector = new List<double>();
            R_Vector = new List<double>();

            q_Matrix = new List<double[]>();
            d_Matrix = new List<double[]>();

            OptCuts = new List<IRange>();

            FirStageVars = new List<string>();
            SecStageVars = new List<string>();
            ScenarioProbability = new List<double>();

            x_v_Sln = new double[1][];

            _SecStageProb = new SecStageProb(numSecStageVars, numSecStageCons);     // empty object
            #endregion
        }
        #endregion

        #region Constructor (ILShaped, string)
        public ILShaped(ILShaped _ILShaped, string filename)
        { 
            // 1.
            #region 1. Initialize Base Class
            base.RootCplex = new Cplex();

            base.RootCplex.ImportModel(@filename);
            base.mEnum = base.RootCplex.GetLPMatrixEnumerator();
            base.mEnum.MoveNext();
            base.lpmatrix = (ILPMatrix)mEnum.Current;

            base.numRows = base.lpmatrix.Nrows;
            base.numCols = base.lpmatrix.Ncols;
            base.A = new double[numRows, numCols];
            base.RHS = new double[numRows];

            base.BCuts = new List<IAddable>();
            base._INumVar = new INumVar[base.RootCplex.Ncols-1];
            for (int i = 0; i < base.RootCplex.Ncols; i++)
            {
                if (i == base.RootCplex.Ncols - 1) this.theta = base.lpmatrix.GetNumVar(i);
                else base._INumVar[i] = base.lpmatrix.GetNumVar(i);
            }

            base.Variables = new List<Tuple<int, double, NumVarType>>();
            base.BinaryVarID = new List<int>();
            base.BinaryVarSln = new List<double>();
            base.myTree = new Dictionary<int, Dictionary<string, TreeNode>>();
            base._LBCandidate = new Dictionary<string, LBCandidate>();
            #endregion

            // 2.
            #region Initialize ILShaped Class Variables
            A_Matrix = _ILShaped.A_Matrix;
            T_Matrix = _ILShaped.T_Matrix;
            W_Matrix = _ILShaped.W_Matrix;
            c_Vector = _ILShaped.c_Vector;
            q_Vector = _ILShaped.q_Vector;
            b_Vector = _ILShaped.b_Vector;
            h_Vector = _ILShaped.h_Vector;
            R_Vector = _ILShaped.R_Vector;

            q_Matrix = _ILShaped.q_Matrix;
            d_Matrix = _ILShaped.d_Matrix;

            OptCuts = new List<IRange>();

            FirStageVars = _ILShaped.FirStageVars;
            SecStageVars = _ILShaped.SecStageVars;
            ScenarioProbability = _ILShaped.ScenarioProbability;

            x_v_Sln = new double[1][];

            _SecStageProb = new SecStageProb(numSecStageVars, numSecStageCons);     // empty object
            #endregion
        }
        #endregion

        public class OptCutGeneration
        {
            #region Declare variables for OptCutGeneration
            public static double qs = 0;
            //public static List<double> lambda_s_S = new List<double>();
            public static Dictionary<double, double> lambda_s_S = new Dictionary<double, double>();
            public static double a = 0;
            #endregion

            public static void Compute_qs(ILShaped _ILShaped, int k)
            {
                _ILShaped._SecStageProb.OptimalityProb.SetOut(null);
                if (_ILShaped._SecStageProb.OptimalityProb.Solve())
                    qs = qs + _ILShaped.ScenarioProbability[k] * _ILShaped._SecStageProb.OptimalityProb.ObjValue;
            }
            
            public static void Generate_OptimalityCut(ILShaped _ILShaped, string SelectedNodeID, ADDITIONALCUTTYPE _ADDITIONALCUTTYPE)
            {
                //**************************************
                // 1. Generate Optimality Cut
                // ( Laport Cut )
                // ( Proposition 5 )
                // ( Proposition 6 ) 
                //**************************************
                double in_S = 0;
                double not_in_S = 0;

                ILinearNumExpr OptimalityCut = _ILShaped._SecStageProb.OptimalityProb.LinearNumExpr();

                double cut_LHS = 0;
                IRange cut = null;

                switch (_ADDITIONALCUTTYPE)
                { 
                    case ADDITIONALCUTTYPE.LAPROTCUT:
                        _ILShaped.OptCuts = new List<IRange>();

                        for (int i = 0; i < _ILShaped._INumVar.Length; i++)
                        { 
                            if(_ILShaped.x_v_Sln[0][i] == 1)
                            {
                                in_S++;
                                OptimalityCut.AddTerm((qs - _ILShaped.BestLB), _ILShaped._INumVar[i]);
                            }
                            else if (_ILShaped.x_v_Sln[0][i] == 0)
                            {
                                not_in_S++;
                                OptimalityCut.AddTerm(-(qs - _ILShaped.BestLB), _ILShaped._INumVar[i]);
                            }
                        }

                        OptimalityCut.AddTerm(-1, _ILShaped.theta);
                        cut_LHS = (qs - _ILShaped.BestLB) * (in_S - 1) - _ILShaped.BestLB;
                        cut = _ILShaped.RootCplex.AddGe(cut_LHS, OptimalityCut, "OptCut");
                        break;
                    case ADDITIONALCUTTYPE.PROPOSITION5:
                        _ILShaped.RootCplex.Remove(_ILShaped.OptCuts.ToArray());
                        _ILShaped.OptCuts = new List<IRange>();

                        for (int i = 0; i < _ILShaped._INumVar.Length; i++)
                        { 
                            if(_ILShaped.x_v_Sln[0][i] == 1)
                            {
                                in_S++;
                                OptimalityCut.AddTerm(a, _ILShaped._INumVar[i]);
                            }
                            else if (_ILShaped.x_v_Sln[0][i] == 0)
                            {
                                not_in_S++;
                                OptimalityCut.AddTerm(-a, _ILShaped._INumVar[i]);
                            }
                        }

                        OptimalityCut.AddTerm(-1, _ILShaped.theta);
                        cut_LHS = -qs + a * in_S ;
                        cut = _ILShaped.RootCplex.AddGe(cut_LHS, OptimalityCut, "OptCut");
                        break;
                    case ADDITIONALCUTTYPE.PROPOSITION6:
                       
                        break;
                }

                //************************************************
                // 1. Mapping Node ID to Node Status
                // 2. Update number of optimality cut into TreeNode
                //************************************************
                // 1.
                _ILShaped.OptCuts.Add(cut);
                string ParentNodeID = "0";
                if (SelectedNodeID.Length != 1) ParentNodeID = SelectedNodeID.Substring(0,SelectedNodeID.Length-1);

                TreeNode _oldTreeNode = _ILShaped.myTree[3][ParentNodeID];

                TreeNode _TreeNode = _ILShaped.myTree[0][SelectedNodeID];
                for (int i = 0; i < _oldTreeNode.OptCuts.Count; i++) _ILShaped.OptCuts.Add(_oldTreeNode.OptCuts[i]);

                // 2.
                _TreeNode.OptCuts = new List<IRange>(_ILShaped.OptCuts);
                _ILShaped.myTree[0][SelectedNodeID] = _TreeNode;

               // _ILShaped.RootCplex.ExportModel("AfterCut.lp");
            }

            public static double Lambda_s(ILShaped _ILShaped, string SelectedNodeID, string DataFolder, ADDITIONALCUTTYPE  _ADDITIONALCUTTYPE)
            {
                //*******************************************************
                // Step 1. Solve Master Prob in Program.Main()
                // Step 2. Calculate S set
                //  << For Loop Section (3. - 5.) >>
                // Step 3. Add Additional Cut to Master Prob. 
                // Step 4. Solve new Master Prob.               
                // Step 5. for each new master sln take into subProb. 
                //         Calculate qs
                // Step 6. Calculate a
                // Step 7. Add Optimaliyt Cut into original Master Prob
                //*******************************************************
                double in_S = 0;
                double not_in_S = 0;
                double Additional_RHS = 0;
                lambda_s_S = new Dictionary<double, double>();

                Cplex N_s_s = new Cplex();
                N_s_s.SetOut(null);
                //N_s_s.ImportModel(string.Format("Node{0}.lp",SelectedNodeID));
                
                // Using CloneManager Copy Master Prob. instead of importmodel
                CloneManager _CloneManager = new SimpleCloneManager(_ILShaped.RootCplex);
                IObjective _IObjective = (IObjective)_ILShaped.RootCplex.GetObjective().MakeClone(_CloneManager);
                N_s_s.Add(_IObjective);
                IAddable _IAddable = (IAddable)(ICopyable)_ILShaped.lpmatrix.MakeClone(_CloneManager);
                for (int j = 0; j < _ILShaped.RootCplex.Ncols; j++)
                {
                    if (j < _ILShaped.RootCplex.Ncols - 1) _ILShaped._INumVar[j] = _ILShaped.lpmatrix.GetNumVar(j);
                    else _ILShaped.theta = _ILShaped.lpmatrix.GetNumVar(j);
                }    
                N_s_s.Add(_IAddable);

                for (int i = 0; i < _ILShaped.myTree[0][SelectedNodeID].ParentID.Count; i++)
                {
                    INumVar LPNumVar = _ILShaped.getLPNumVar(_ILShaped.myTree[0][SelectedNodeID].SelectedVarID[i]);
                    if (_ILShaped.myTree[0][SelectedNodeID].BCutsRHS[i] == 0)
                    { 
                        IRange cut = N_s_s.AddLe(N_s_s.Prod(1.0, LPNumVar),0);
                        N_s_s.AddCut(cut);
                    }
                    else if (_ILShaped.myTree[0][SelectedNodeID].BCutsRHS[i] == 1)
                    {
                        IRange cut = N_s_s.AddGe(N_s_s.Prod(1.0, LPNumVar),1);
                        N_s_s.AddCut(cut);
                    }
                }


                // 2.
                ILinearNumExpr AdditionalCut = N_s_s.LinearNumExpr();
                for (int i = 0; i < _ILShaped._INumVar.Length; i++)
                {
                    if (_ILShaped.myTree[0][SelectedNodeID].SelectedVarID.Contains(i)) continue;
                    else if (_ILShaped.x_v_Sln[0][i] == 1)
                    {
                        in_S++;
                        AdditionalCut.AddTerm(1, _ILShaped._INumVar[i]);
                    }
                    else if (_ILShaped.x_v_Sln[0][i] == 0)
                    {
                        not_in_S++;
                        AdditionalCut.AddTerm(-1, _ILShaped._INumVar[i]);
                    }
                }

                // Loop Section for Step 3 - Step 5
                double[][] s_x_Sln = new double[1][];
                double s_theta_Sln = 0;
                for (int s = 0; s <= in_S; s++)
                {
                    // 3. & 4.
                    Additional_RHS = in_S - s;
                    IRange cut = N_s_s.AddEq(AdditionalCut, Additional_RHS);
                    //N_s_s.ExportModel("ImprovedCut_" + s + ".lp");
                    N_s_s.Solve();
                    IEnumerator mEnum = N_s_s.GetLPMatrixEnumerator();
                    mEnum.MoveNext();
                    ILPMatrix lpmatrix = (ILPMatrix)mEnum.Current;

                    // 5.
                    s_x_Sln[0] = new double[lpmatrix.Ncols - 1];
                    for (int i = 0; i < lpmatrix.Ncols; i++)
                    {
                        if (i == lpmatrix.Ncols - 1) s_theta_Sln = N_s_s.GetValue(lpmatrix.GetNumVar(i));
                        else s_x_Sln[0][i] = N_s_s.GetValue(lpmatrix.GetNumVar(i));
                    }

                    qs = 0;
                    for (int k = 1; k <= _ILShaped.ScenarioProbability.Count; k++)
                    {
                        Data.getInputData_forScenario(DataFolder, k.ToString());
                        _ILShaped.numSecStageVars = _ILShaped.SecStageVars.Count;
                        _ILShaped.numSecStageCons = _ILShaped.W_Matrix.Count;
                        _ILShaped._SecStageProb = new SecStageProb(_ILShaped.numSecStageVars, _ILShaped.numSecStageCons);
                        Model.BuildSNeighborProb(_ILShaped, s_x_Sln);
                        Compute_qs(_ILShaped, k-1);
                    }
                    lambda_s_S.Add(s, qs);

                    N_s_s.Remove(cut);
                }

                // 6.1 Check Proposition 7. 
                double gap = 99999999;
                bool isNonIncreasing = true;
                for (int i = 0; i < lambda_s_S.Count-1; i++)
                {
                    double temp = Math.Abs(lambda_s_S[i] - lambda_s_S[i + 1]);
                    if (temp < gap) gap = temp;
                    else
                    {
                        isNonIncreasing = false;
                        break;
                    }
                }

                // 6.2 Add Improved Optimality Cuts
                double lambda = 0;
                switch (_ADDITIONALCUTTYPE)
                { 
                    case ADDITIONALCUTTYPE.PROPOSITION5:
                        if (isNonIncreasing) a = qs - lambda_s_S[1];
                        else a = Math.Max(qs - lambda_s_S[1], (qs - _ILShaped.BestLB) / 2);
                        break;
                    case ADDITIONALCUTTYPE.PROPOSITION6:
                        //lambda_s_S.Sort((x,y) => (x.CompareTo(y)));
                        //lambda = lambda_s_S[0];
                        break;
                }
                return a;
            }
        }

        public class Model
        {
            public static void BuildMasterProb()
            {
                Program._ILShaped.numFirStageVars = Program._ILShaped.FirStageVars.Count;
                Program._ILShaped.numFirStageRHS = Program._ILShaped.b_Vector.Count;
                Program._ILShaped._INumVar = new INumVar[Program._ILShaped.numFirStageVars];
                for (int j = 0; j < Program._ILShaped.numFirStageVars; j++) 
                    Program._ILShaped._INumVar[j] = Program._ILShaped.RootCplex.NumVar(0, 1, NumVarType.Bool, string.Format("x_{0}", j));

                //*************************************************************
                // Create Objective Function
                // The objective coefficient for teh first-stage variables
                // The ovjective coefficient for Theta
                //*************************************************************
                ILinearNumExpr Master_obj = Program._ILShaped.RootCplex.LinearNumExpr();
                Master_obj.AddTerms(Program._ILShaped.c_Vector.ToArray(), Program._ILShaped._INumVar);
                Master_obj.AddTerm(1, Program._ILShaped.theta);
                Program._ILShaped.RootCplex.AddMinimize(Master_obj);

                // First Stage Constraint (Knapsack)
                ILinearNumExpr Master_con = Program._ILShaped.RootCplex.LinearNumExpr();
                for (int i = 0; i < Program._ILShaped.numFirStageRHS; i++)
                {
                    Master_con.AddTerms(Program._ILShaped.A_Matrix[i].ToArray(), Program._ILShaped._INumVar);
                    Program._ILShaped.RootCplex.AddLe(Master_con, Program._ILShaped.b_Vector[i]);
                    Master_con.Clear();
                }

                Program._ILShaped.RootCplex.ExportModel("Master_0.lp");
            }

            //public static void BuildFeasibilityProb()
            //{
            
            //}

            public static void BuildOptimalityProb(ILShaped _ILShaped)
            {
                _ILShaped.numSecStageVars = _ILShaped.W_Matrix[0].Length;
                _ILShaped.numSecStageCons = _ILShaped.h_Vector.Count;

                // Create Objective Function
                ILinearNumExpr Optimality_obj = _ILShaped._SecStageProb.OptimalityProb.LinearNumExpr();
                Optimality_obj.AddTerm(_ILShaped.q_Vector[0], _ILShaped._SecStageProb.y_0);
                _ILShaped._SecStageProb.OptimalityProb.AddMinimize(Optimality_obj);

                _ILShaped.x_v_Sln[0] = _ILShaped.RootCplex.GetValues(_ILShaped._INumVar);
                _ILShaped.theta_v_Sln = _ILShaped.RootCplex.GetValue(_ILShaped.theta);

                // Create Constraints
                ILinearNumExpr Optimality_con = _ILShaped._SecStageProb.OptimalityProb.LinearNumExpr();
                for (int i = 0; i < _ILShaped.numSecStageCons; i++)
                {
                    for (int j = 0; j < _ILShaped.numSecStageVars; j++)
                    {
                        if (j == 0) Optimality_con.AddTerm(_ILShaped.W_Matrix[i][j], _ILShaped._SecStageProb.y_0);
                        else Optimality_con.AddTerm(_ILShaped.W_Matrix[i][j], _ILShaped._SecStageProb.mySecStageVars[j]);
                    }

                    double myRHS = _ILShaped.h_Vector[i];
                    for (int j = 0; j < _ILShaped.numFirStageVars; j++) myRHS = myRHS + _ILShaped.T_Matrix[i][j] * _ILShaped.x_v_Sln[0][j];

                    Optimality_con.AddTerm(_ILShaped.R_Vector[i], _ILShaped._SecStageProb.R);
                    if (i == 0) _ILShaped._SecStageProb.mySecStageCons[i] = _ILShaped._SecStageProb.OptimalityProb.AddEq(Optimality_con, myRHS);
                    else _ILShaped._SecStageProb.mySecStageCons[i] = _ILShaped._SecStageProb.OptimalityProb.AddLe(Optimality_con, myRHS);
                    Optimality_con.Clear();
                }
                //_ILShaped._SecStageProb.OptimalityProb.ExportModel(string.Format("Optimality_{0}.lp",k));
            }

            public static void BuildSIPLIBMaster(ILShaped _ILShaped)
            {
                _ILShaped.numFirStageVars = _ILShaped.FirStageVars.Count;
                _ILShaped.numFirStageRHS = _ILShaped.b_Vector.Count;
                _ILShaped._INumVar = new INumVar[_ILShaped.numFirStageVars];
                for (int j = 0; j < _ILShaped.numFirStageVars; j++) _ILShaped._INumVar[j] = _ILShaped.RootCplex.NumVar(0, 1, NumVarType.Bool, string.Format("x_{0}", j));

                //*************************************************************
                // Create Objective Function
                // The objective coefficient for teh first-stage variables
                // The ovjective coefficient for Theta
                //*************************************************************
                ILinearNumExpr Master_obj = _ILShaped.RootCplex.LinearNumExpr();
                Master_obj.AddTerms(_ILShaped.c_Vector.ToArray(), _ILShaped._INumVar);
                Master_obj.AddTerm(1, _ILShaped.theta);
                _ILShaped.RootCplex.AddMinimize(Master_obj);

                // First Stage Constraint
                ILinearNumExpr Master_con1 = _ILShaped.RootCplex.LinearNumExpr();
                for (int i = 0; i < _ILShaped.numFirStageRHS; i++)
                {
                    for (int j = 0; j < _ILShaped.numFirStageVars; j++) Master_con1.AddTerm(_ILShaped.A_Matrix[i][j], _ILShaped._INumVar[j]);
                    _ILShaped.RootCplex.AddGe(Master_con1, _ILShaped.b_Vector[i]);
                }

                ILinearNumExpr Master_con2 = _ILShaped.RootCplex.LinearNumExpr();
                int w = 3;
                for (int j = 0; j < _ILShaped.numFirStageVars; j++) Master_con2.AddTerm(1, _ILShaped._INumVar[j]);
                _ILShaped.RootCplex.AddGe(Master_con2, w);


                _ILShaped.RootCplex.ExportModel("Master_0.lp");
            }

            public static void BuildSIPLIBSubProb(ILShaped _ILShaped)
            {
                _ILShaped.numSecStageVars = _ILShaped.SecStageVars.Count;
                _ILShaped.numSecStageCons = _ILShaped.h_Vector.Count;

                // Create Objective Function
                ILinearNumExpr subProb_obj = _ILShaped._SecStageProb.OptimalityProb.LinearNumExpr();
                for (int k = 0; k < _ILShaped.ScenarioProbability.Count; k++)
                {
                    for (int i = 0; i < _ILShaped.numSecStageCons; i++)
                    {
                        for (int j = 0; j < _ILShaped.numFirStageVars; j++) subProb_obj.AddTerm(_ILShaped.ScenarioProbability[k] * _ILShaped.q_Matrix[i][j], _ILShaped._SecStageProb.y[i][j]);
                    }

                    for (int j = 0; j < _ILShaped.numFirStageVars; j++) subProb_obj.AddTerm(-1 * _ILShaped.q_Vector[j], _ILShaped._SecStageProb.y_j_0[j]);
                }
                _ILShaped._SecStageProb.OptimalityProb.AddMinimize(subProb_obj);

                // Create Constraint 1
                ILinearNumExpr subProb_con1 = _ILShaped._SecStageProb.OptimalityProb.LinearNumExpr();
                for (int j = 0; j < _ILShaped.numFirStageVars; j++)
                {
                    for (int i = 0; i < _ILShaped.numSecStageCons; i++) subProb_con1.AddTerm(_ILShaped.d_Matrix[i][j], _ILShaped._SecStageProb.y[i][j]);
                    subProb_con1.AddTerm(-1, _ILShaped._SecStageProb.y_j_0[j]);
                    double RHS = _ILShaped.T_Matrix[j][0] * _ILShaped.x_v_Sln[0][j];
                    _ILShaped._SecStageProb.OptimalityProb.AddGe(subProb_con1, RHS);

                    subProb_con1.Clear();
                }

                // Create Constraint 2
                ILinearNumExpr subProb_con2 = _ILShaped._SecStageProb.OptimalityProb.LinearNumExpr();
                for (int i = 0; i < _ILShaped.numSecStageCons; i++)
                {
                    for (int j = 0; j < _ILShaped.numFirStageVars; j++) subProb_con2.AddTerm(_ILShaped.W_Matrix[i][j], _ILShaped._SecStageProb.y[i][j]);
                    _ILShaped._SecStageProb.OptimalityProb.AddEq(subProb_con2, _ILShaped.h_Vector[i]);
                    subProb_con2.Clear();
                }

                //_ILShaped._SecStageProb.OptimalityProb.ExportModel("SIPLIBsub.lp");
            }

            // Not yet
            public static void BuildSNeighborSIPLIB(ILShaped _ILShaped, double[][] s_x_Sln)
            {
                _ILShaped.numFirStageVars = _ILShaped.FirStageVars.Count;
                _ILShaped.numFirStageRHS = _ILShaped.b_Vector.Count;
                _ILShaped._INumVar = new INumVar[_ILShaped.numFirStageVars];
                for (int j = 0; j < _ILShaped.numFirStageVars; j++) _ILShaped._INumVar[j] = _ILShaped.RootCplex.NumVar(0, 1, NumVarType.Bool, string.Format("x_{0}", j));

                //*************************************************************
                // Create Objective Function
                // The objective coefficient for teh first-stage variables
                // The ovjective coefficient for Theta
                //*************************************************************
                ILinearNumExpr Master_obj = _ILShaped.RootCplex.LinearNumExpr();
                Master_obj.AddTerms(_ILShaped.c_Vector.ToArray(), _ILShaped._INumVar);
                Master_obj.AddTerm(1, _ILShaped.theta);
                _ILShaped.RootCplex.AddMinimize(Master_obj);

                // First Stage Constraint
                ILinearNumExpr Master_con1 = _ILShaped.RootCplex.LinearNumExpr();
                for (int i = 0; i < _ILShaped.numFirStageRHS; i++)
                {
                    for (int j = 0; j < _ILShaped.numFirStageVars; j++) Master_con1.AddTerm(_ILShaped.A_Matrix[i][j], _ILShaped._INumVar[j]);
                    _ILShaped.RootCplex.AddGe(Master_con1, _ILShaped.b_Vector[i]);
                }

                ILinearNumExpr Master_con2 = _ILShaped.RootCplex.LinearNumExpr();
                int w = 3;
                for (int j = 0; j < _ILShaped.numFirStageVars; j++) Master_con2.AddTerm(1, _ILShaped._INumVar[j]);
                _ILShaped.RootCplex.AddGe(Master_con2, w);
            }

            public static void BuildSNeighborProb(ILShaped _ILShaped, double[][] s_x_Sln)
            {
                _ILShaped.numSecStageVars = _ILShaped.W_Matrix[0].Length;
                _ILShaped.numSecStageCons = _ILShaped.h_Vector.Count;

                // Create Objective Function
                ILinearNumExpr Optimality_obj = _ILShaped._SecStageProb.OptimalityProb.LinearNumExpr();
                Optimality_obj.AddTerm(_ILShaped.q_Vector[0], _ILShaped._SecStageProb.y_0);
                _ILShaped._SecStageProb.OptimalityProb.AddMinimize(Optimality_obj);

                // Create Constraints
                ILinearNumExpr Optimality_con = _ILShaped._SecStageProb.OptimalityProb.LinearNumExpr();
                for (int i = 0; i < _ILShaped.numSecStageCons; i++)
                {
                    for (int j = 0; j < _ILShaped.numSecStageVars; j++)
                    {
                        if (j == 0) Optimality_con.AddTerm(_ILShaped.W_Matrix[i][j], _ILShaped._SecStageProb.y_0);
                        else Optimality_con.AddTerm(_ILShaped.W_Matrix[i][j], _ILShaped._SecStageProb.mySecStageVars[j]);
                    }

                    double myRHS = _ILShaped.h_Vector[i];
                    for (int j = 0; j < _ILShaped.numFirStageVars; j++) myRHS = myRHS + _ILShaped.T_Matrix[i][j] * s_x_Sln[0][j];

                    Optimality_con.AddTerm(_ILShaped.R_Vector[i], _ILShaped._SecStageProb.R);
                    if (i == 0) _ILShaped._SecStageProb.mySecStageCons[i] = _ILShaped._SecStageProb.OptimalityProb.AddEq(Optimality_con, myRHS);
                    else _ILShaped._SecStageProb.mySecStageCons[i] = _ILShaped._SecStageProb.OptimalityProb.AddLe(Optimality_con, myRHS);
                    Optimality_con.Clear(); 
                }
                //_ILShaped._SecStageProb.OptimalityProb.ExportModel("AdditionalOptimality.lp");
            }

            //public static void AddFeasibilityCut()
            //{

            //}

            public static void Compute_param(int k, ILShaped _ILShaped, string SelectedNodeID, string DataFolder, OPTIMALITYCUTS _OPTIMALITYCUTS)
            {
                switch (_OPTIMALITYCUTS)
                { 
                    case OPTIMALITYCUTS.LAPROT:
                        OptCutGeneration.Compute_qs(_ILShaped, k);
                        break;
                    case OPTIMALITYCUTS.SNEIGHBOR:
                        //CutGeneration.Generate_s_Neighbor(_ILShaped);
                        OptCutGeneration.Lambda_s(_ILShaped, SelectedNodeID, DataFolder, ADDITIONALCUTTYPE.PROPOSITION5);
                        break;
                }
            }
        }

        public override void CreateRootNode(TreeNode _TreeNode)
        {
            for (int j = 0; j < base.lpmatrix.NumVars.Length; j++)
            {
                if (this.lpmatrix.GetNumVar(j).Type.ToString() == "Bool")
                {
                    base.BinaryVarID.Add(j);
                    IConversion LPrelax = base.RootCplex.Conversion(base.lpmatrix.GetNumVar(j), NumVarType.Float);
                    base.RootCplex.Add(LPrelax);
                }
            }

            base.mEnum = base.RootCplex.GetLPMatrixEnumerator();
            base.mEnum.MoveNext();
            base.lpmatrix = (ILPMatrix)base.mEnum.Current;

            //base.RootCplex.ExportModel("Root_Node.lp");
            base.myTree[0].Add("0", _TreeNode);
            _TreeNode.NodeStack.Push("0");
            _TreeNode.NodeQueue.Enqueue("0");
        }

        public override void CreateChildNode(int BranchingVar, string SelectedNodeID, BRANCHING _BRANCHING){base.CreateChildNode(BranchingVar, SelectedNodeID, _BRANCHING);}

        public override void CreateNode(string SelectedNodeID)
        {
            base.RootCplex.Remove(base.BCuts.ToArray());
            this.BCuts = new List<IAddable>();

            base.RootCplex.Remove(this.OptCuts.ToArray());

            // Add BCuts
            for (int i = 0; i < base.myTree[0][SelectedNodeID].ParentID.Count; ++i)
            {
                INumVar LPNumVar = getLPNumVar(base.myTree[0][SelectedNodeID].SelectedVarID[i]);
                if (base.myTree[0][SelectedNodeID].BCutsRHS[i] == 0)
                {
                    IRange Cut = base.RootCplex.AddLe(base.RootCplex.Prod(1.0, LPNumVar), 0);
                    base.BCuts.Add(base.RootCplex.AddCut(Cut));
                }
                else if (base.myTree[0][SelectedNodeID].BCutsRHS[i] == 1)
                {
                    IRange Cut = base.RootCplex.AddGe(base.RootCplex.Prod(1.0, LPNumVar), 1);
                    base.BCuts.Add(base.RootCplex.AddCut(Cut));
                }
            }

            // Add Optimaliyt Cuts
            for (int i = 0; i < base.myTree[0][SelectedNodeID].OptCuts.Count; i++) base.RootCplex.AddCut(base.myTree[0][SelectedNodeID].OptCuts[i]);

            //this.RootCplex.ExportModel(string.Format("Node{0}.lp", SelectedNodeID));
        }

        public override string NodeSelection(TreeNode _TreeNode, NODESELECTIONTYPE _NODESELECTIONTYPE){return base.NodeSelection(_TreeNode, _NODESELECTIONTYPE);}

        public override int VariableSelection(VARIABLESELECTION _VARIABLESELECTION){return base.VariableSelection(_VARIABLESELECTION);}

        public override void SolveNode(string SelectedNodeID, int algorithm)
        {
            base.RootCplex.SetOut(null);
            base.RootCplex.SetParam(Cplex.Param.RandomSeed, 0);
            base.RootCplex.SetParam(Cplex.Param.RootAlgorithm, algorithm);
            SetAlgorithm.Preprocessing_off(base.RootCplex);
            SetAlgorithm.Cuts_off(base.RootCplex);

            base.RootCplex.Solve();

            this.BinaryVarSln = new List<double>();

            foreach (int id in base.BinaryVarID) this.BinaryVarSln.Add(base.RootCplex.GetValue(base._INumVar[id]));

            TreeNode _TreeNode = base.myTree[0][SelectedNodeID];
            _TreeNode.LPR = base.RootCplex.ObjValue;
            _TreeNode.AllSln = base.RootCplex.GetValues(base._INumVar);
            _TreeNode.BinSln = base.BinaryVarSln;

            base.myTree[0][SelectedNodeID] = _TreeNode;
            System.Console.WriteLine(string.Format("{0} objective value: {1}", SelectedNodeID, _TreeNode.LPR));
        }

        public override bool CheckIPFeasible()
        {
            double[] solutions = base.RootCplex.GetValues(base._INumVar);
            double solution = base.RootCplex.GetValue(this.theta);

            // Check if all Binary Variables are Integer
            base.isIPFeasible = Array.TrueForAll(solutions, value => { return Math.Abs(value - Math.Round(value)) < 0.00000000001; });
            return base.isIPFeasible;
        }

        public override void IdentifyNodeStatus(string SelectedNodeID, VARIABLESELECTION _VARIABLESELECTION)
        {
            CheckIPFeasible();

            // NodeStatus (1:Best IP, 2:Pruned IP, 3:BestLP / Regular LP, 4: Pruned LP, 5: LP Infeasible)
            TreeNode _TreeNode = base.myTree[0][SelectedNodeID];
            UpdateLB(_TreeNode, SelectedNodeID);

            if (base.isIPFeasible)
            {
                if (_TreeNode.LPR <= base.BestUB) base.NodeStatus = 1;
                else base.NodeStatus = 2;
            }
            else if (_TreeNode.LPR <= base.BestUB)
            {
                if (_TreeNode.LPR >= base.BestLB) base.NodeStatus = 3;
                   
                //Variable Selection and Create two Child Node
                BranchingID = VariableSelection(VARIABLESELECTION.MOSTFRAC);
                
                CreateChildNode(BranchingID, SelectedNodeID, BRANCHING.LEFT);
                CreateChildNode(BranchingID, SelectedNodeID, BRANCHING.RIGHT);
            }
            else base.NodeStatus = 4;

            if (base.NodeStatus != 1)
            {
                base.BestLB = _LBCandidate.First().Value.LPR;
                base.myTree[NodeStatus].Add(SelectedNodeID, _TreeNode);
                base.myTree[0].Remove(SelectedNodeID);          // Remove incumbent node
            }
        }

        public override void UpdateLB(TreeNode _TreeNode, string SelectedNodeID)
        {
            //****************************************************
            // 1. Insert inCumbentNode into SolvedList
            // 2. Sort by Linq
            // 3. Delete parentNode if its children are explored
            //****************************************************
            LBCandidate inCumbentNode = new LBCandidate(_TreeNode.LPR);
            int numLBCandidate = _LBCandidate.Count;

            // 1.
            if (_LBCandidate.ContainsKey(SelectedNodeID)) _LBCandidate[SelectedNodeID] = inCumbentNode;
            else _LBCandidate.Add(SelectedNodeID, inCumbentNode);

            // 2.
            _LBCandidate = _LBCandidate.OrderBy(Date => Date.Value.LPR).ToDictionary(KeyValue => KeyValue.Key, KeyValue => KeyValue.Value);

            // 3.
            string ParentNodeID = SelectedNodeID.Substring(0, SelectedNodeID.Length - 1);
            if (_LBCandidate.ContainsKey(ParentNodeID))
            {
                LBCandidate UpdateValue;
                _LBCandidate.TryGetValue(ParentNodeID, out UpdateValue);
                UpdateValue.ChildNodes.Add(SelectedNodeID);

                if (UpdateValue.ChildNodes.Count == 2) _LBCandidate.Remove(ParentNodeID);
                else
                    _LBCandidate[ParentNodeID] = UpdateValue;
            }
        }

        public override void UpdateLB(string SelectedNodeID)
        {
            // 3.
            string ParentNodeID = SelectedNodeID.Substring(0, SelectedNodeID.Length - 1);
            if (_LBCandidate.ContainsKey(ParentNodeID))
            {
                LBCandidate UpdateValue;
                _LBCandidate.TryGetValue(ParentNodeID, out UpdateValue);
                UpdateValue.ChildNodes.Add(SelectedNodeID);


                if (UpdateValue.ChildNodes.Count == 2)
                {
                    _LBCandidate.Remove(ParentNodeID);
                    this.BestLB = _LBCandidate.First().Value.LPR;
                }
                else
                    _LBCandidate[ParentNodeID] = UpdateValue;
            }
        }

        public override INumVar getLPNumVar(int idx) { return base.getLPNumVar(idx); }
    }
}
