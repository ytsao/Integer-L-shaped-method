using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ILOG.CPLEX;
using ILOG.Concert;
using System.Collections;

namespace IntegerLShapedMethod
{
    public enum BRANCHING{ LEFT,RIGHT}

    public enum NODESELECTIONTYPE { DFS, BFS, BESTBOUND}

    public enum VARIABLESELECTION { MOSTFRAC, PSEUDOCOST}

    public struct TreeNode
    {   
        public double LPR;
        public double Theta;
        public double[] AllSln;
        public List<double> BinSln;
        public List<string> ParentID;
        public List<int> SelectedVarID;
        public List<int> BCutsRHS;

        public Stack<string> NodeStack;
        public Queue<string> NodeQueue;

        public List<double[]> DCutCoeff;
        public List<double> DCutRHS;

        public List<IRange> OptCuts;

        public TreeNode(TreeNode _TreeNode)
        {
            LPR = _TreeNode.LPR;
            Theta = _TreeNode.Theta;
            AllSln = _TreeNode.AllSln;
            BinSln = new List<double>(_TreeNode.BinSln);
            ParentID = new List<string>(_TreeNode.ParentID);
            SelectedVarID = new List<int>(_TreeNode.SelectedVarID);
            BCutsRHS = new List<int>(_TreeNode.BCutsRHS);

            NodeStack = new Stack<string>();
            NodeQueue = new Queue<string>();

            DCutCoeff = new List<double[]>(_TreeNode.DCutCoeff);
            DCutRHS = new List<double>(_TreeNode.DCutRHS);

            OptCuts = new List<IRange>();
        }

        public TreeNode(TreeNode _TreeNode, List<IRange> newOptCuts)
        {
            LPR = _TreeNode.LPR;
            Theta = _TreeNode.Theta;
            AllSln = _TreeNode.AllSln;
            BinSln = new List<double>(_TreeNode.BinSln);
            ParentID = new List<string>(_TreeNode.ParentID);
            SelectedVarID = new List<int>(_TreeNode.SelectedVarID);
            BCutsRHS = new List<int>(_TreeNode.BCutsRHS);

            NodeStack = new Stack<string>();
            NodeQueue = new Queue<string>();

            DCutCoeff = new List<double[]>(_TreeNode.DCutCoeff);
            DCutRHS = new List<double>(_TreeNode.DCutRHS);

            OptCuts = new List<IRange>(newOptCuts);
        }

        public TreeNode(BranchAndBound BB)
        { 
            //Create tree nodes in five different status
            for (int i = 0; i <= 5; i++) BB.myTree.Add(i, new Dictionary<string, TreeNode>());
            
            //Update incumbent node statistics
            LPR = 0;
            Theta = 0;
            AllSln = new double[BB.numCols];
            BinSln = new List<double>();
            ParentID = new List<string>();
            SelectedVarID = new List<int>();
            BCutsRHS = new List<int>();

            NodeStack = new Stack<string>();
            NodeQueue = new Queue<string>();
            
            DCutCoeff = new List<double[]>();
            DCutRHS = new List<double>();

            OptCuts = new List<IRange>();
        }

        public TreeNode(BranchAndBound BB, int numStatus)
        { 
            //Create tree nodes in five different status
            for (int i = 0; i <= numStatus; i++) BB.myTree.Add(i, new Dictionary<string, TreeNode>());
            
            //Update incumbent node statistics
            LPR = 0;
            Theta = 0;
            AllSln = new double[BB.numCols];
            BinSln = new List<double>();
            ParentID = new List<string>();
            SelectedVarID = new List<int>();
            BCutsRHS = new List<int>();

            NodeStack = new Stack<string>();
            NodeQueue = new Queue<string>();

            DCutCoeff = new List<double[]>();
            DCutRHS = new List<double>();

            OptCuts = new List<IRange>();
        }
    }

    public struct LBCandidate
    {
        public double LPR;
        public List<string> ChildNodes;

        public LBCandidate(double NodeLPR)
        {
            LPR = NodeLPR;
            ChildNodes = new List<string>();
        }
    }

    public class BranchAndBound
    {
        public Cplex RootCplex;
        public INumVar[] _INumVar;
        public IEnumerator mEnum;
        public ILPMatrix lpmatrix;

        public int numRows;
        public int numCols;
        public double[,] A;
        public double[] RHS;

        //*****************************************************************************************************************
        // 1. Variable = new <index, coefficient, var type>()
        // 2. Binary Variables Indices
        // 3. Binary Variables' Solution
        // 4. myTree = new Dictionary<status, Dictionary<index, TreeNode>>
        // (Status 0 = Unexplored Node; 1 = was Best UB; 2 = IP Feasible; 3 = LP was Best LB; 4 = Pruned Node(LB > UB)
        // ; 5 = Regular LP; 6 = LP Infeasible)
        // 5. Optimal Node <optimal sln, index, var sln>
        //*****************************************************************************************************************
        public List<Tuple<int, double, NumVarType>> Variables;
        public List<int> BinaryVarID;
        public List<double> BinaryVarSln;
        public Dictionary<int, Dictionary<string, TreeNode>> myTree;

        public List<IAddable> BCuts;
        public Dictionary<string, LBCandidate> _LBCandidate;

        public int BranchingID;
        public string BestNodeID = "-1";
        public int NodeStatus = -1;
        //public double BestUB = 8.5832619653e+02;    //Best IP sln
        public double BestUB = 999999;
        public double BestLB = -999999;             //Best LP sln
        public bool isIPFeasible;

        public BranchAndBound() { }

        public BranchAndBound(string filename)
        {
            RootCplex = new Cplex();

            RootCplex.ImportModel(@filename);
            mEnum = RootCplex.GetLPMatrixEnumerator();
            mEnum.MoveNext();
            lpmatrix = (ILPMatrix)mEnum.Current;

            numRows = lpmatrix.Nrows;
            numCols = lpmatrix.Ncols;
            A = new double[numRows, numCols];
            RHS = new double[numRows];

            Variables = new List<Tuple<int, double, NumVarType>>();
            BinaryVarID = new List<int>();
            BinaryVarSln = new List<double>();
            myTree = new Dictionary<int, Dictionary<string, TreeNode>>();
            BCuts = new List<IAddable>();
            _LBCandidate = new Dictionary<string, LBCandidate>();

            _INumVar = new INumVar[RootCplex.Ncols];
            for (int i = 0; i < RootCplex.Ncols; i++) _INumVar[i] = lpmatrix.GetNumVar(i);
        }

        public class SearchTree
        {
            public static string DepthFirstSearch(TreeNode _TreeNode)
            {
                if(_TreeNode.NodeStack.Count > 0) return _TreeNode.NodeStack.Pop();
                return "-1";
            }
            
            public static string BreathFirstSearch(TreeNode _TreeNode)
            {
               if(_TreeNode.NodeQueue.Count > 0) return _TreeNode.NodeQueue.Dequeue();
               return "-1";
            }

            public static string BestBound(TreeNode _TreeNode, Dictionary<string, LBCandidate> _LBCandidate, string BestNodeID)
            {
                if (_TreeNode.NodeQueue.Count == 0 && _LBCandidate.Keys.First() != BestNodeID)
                {
                    string CurBestLB_NodeID = _LBCandidate.Keys.First();
                    _TreeNode.NodeQueue.Enqueue(CurBestLB_NodeID + "0");
                    _TreeNode.NodeQueue.Enqueue(CurBestLB_NodeID + "1");
                    return _TreeNode.NodeQueue.Dequeue();
                }
                else if (_TreeNode.NodeQueue.Count > 0) return _TreeNode.NodeQueue.Dequeue();
                return "-1";
            }
        }

        public class BranchingStrategy
        {
            public static int MostFractional(List<double> BinaryVarSln)
            {
                int FracVarID = -1;
                double mostFractional = 0.5;    //Since 0 and 1

                // Get the fractional variable that is closet to 0.5
                foreach (var item in BinaryVarSln.Select((value, index) => new { value, index }))
                {
                    double FractionalValue = Math.Abs(0.5 - item.value);
                    if ( FractionalValue < mostFractional)
                    {
                        mostFractional = FractionalValue;
                        FracVarID = item.index;
                    }
                }
                return FracVarID;
            }

            public static int PseudoCost()
            {
                return 1;
            }
        }

        public virtual void CreateRootNode(TreeNode _TreeNode)
        {
            for (int j = 0; j < this.lpmatrix.NumVars.Length; ++j)
            { 
                //String temp = ObjElement[2*j+2]; // Warning Customization Setting Param of Ceof
                //this.Variables.Add(Tuple.Create(j, Convert.ToDouble(temp), this.lpmatrix.GetNumVar(j).Type));
                if (this.lpmatrix.GetNumVar(j).Type.ToString() == "Bool")
                {
                    this.BinaryVarID.Add(j);
                    IConversion LPrelax = this.RootCplex.Conversion(this.lpmatrix.GetNumVar(j), NumVarType.Float);
                    this.RootCplex.Add(LPrelax);
                }
            }
            this.mEnum = this.RootCplex.GetLPMatrixEnumerator();
            this.mEnum.MoveNext();
            this.lpmatrix = (ILPMatrix)this.mEnum.Current;

            this.RootCplex.ExportModel("Root_Node.lp");
            this.myTree[0].Add("0", _TreeNode);
            _TreeNode.NodeStack.Push("0");
            _TreeNode.NodeQueue.Enqueue("0");

            Program._DisjunctiveCuts.C3LPInitialization();
        }

        public virtual void CreateChildNode(int BranchingVar, string SelectedNodeID, BRANCHING _BRANCHING)
        {
            //******************************************************
            // Initialize left & right node
            // 1. inherit from parent node
            // 2. Add Parent node id (left : -1 , right : -2)
            // 3. Add Selected Variable index (Branching variable)
            // 4. Add Cut RHS
            //*******************************************************
            TreeNode _TreeNode = new TreeNode(this.myTree[0][SelectedNodeID]);
            switch (_BRANCHING)
            {
                case BRANCHING.LEFT:
                    _TreeNode.ParentID.Add(SelectedNodeID);
                    _TreeNode.SelectedVarID.Add(BranchingVar);
                    _TreeNode.BCutsRHS.Add(0);
                    this.myTree[0].Add(SelectedNodeID + "0", _TreeNode);
                    Program._TreeNode.NodeStack.Push(SelectedNodeID + "0");
                    Program._TreeNode.NodeQueue.Enqueue(SelectedNodeID + "0");
                    break;
                case BRANCHING.RIGHT:
                    _TreeNode.ParentID.Add(SelectedNodeID);
                    _TreeNode.SelectedVarID.Add(BranchingVar);
                    _TreeNode.BCutsRHS.Add(1);
                    this.myTree[0].Add(SelectedNodeID + "1", _TreeNode);
                    Program._TreeNode.NodeStack.Push(SelectedNodeID + "1");
                    Program._TreeNode.NodeQueue.Enqueue(SelectedNodeID + "1");
                    break;
            }
        }

        public virtual void CreateNode(string SelectedNodeID)
        {
            this.RootCplex.Remove(this.BCuts.ToArray());
            this.BCuts = new List<IAddable>();

            for (int i = 0; i < this.myTree[0][SelectedNodeID].ParentID.Count; ++i)
            {
                if (this.myTree[0][SelectedNodeID].BCutsRHS[i] == 0)
                {
                    IRange Cut = this.RootCplex.AddLe(this.RootCplex.Prod(1.0, getLPNumVar(this.myTree[0][SelectedNodeID].SelectedVarID[i])), 0);
                    this.BCuts.Add(this.RootCplex.AddCut(Cut));
                }
                else if (this.myTree[0][SelectedNodeID].BCutsRHS[i] == 1)
                {
                    IRange Cut = this.RootCplex.AddGe(this.RootCplex.Prod(1.0, getLPNumVar(this.myTree[0][SelectedNodeID].SelectedVarID[i])), 1);
                    this.BCuts.Add(this.RootCplex.AddCut(Cut));
                }
            }
            //this.RootCplex.ExportModel(string.Format("Node{0}.lp",SelectedNodeID));
        }

        public virtual string NodeSelection(TreeNode _TreeNode, NODESELECTIONTYPE _NODESELECTIONTYPE)
        {
            //foreach (var _myTree in this.myTree[0]) return _myTree.Key;
            //return -1;

            switch (_NODESELECTIONTYPE)
            {
                // 1. DFS (LIFO)
                case NODESELECTIONTYPE.DFS:
                    return SearchTree.DepthFirstSearch(_TreeNode);

                // 2. BFS
                case NODESELECTIONTYPE.BFS:
                    return SearchTree.BreathFirstSearch(_TreeNode);

                // 3. Best BOUND (FIFO)
                case NODESELECTIONTYPE.BESTBOUND:
                    return SearchTree.BestBound(_TreeNode, _LBCandidate, BestNodeID);
            }
            return "-1";
        }

        //Most fractional as example
        public virtual int VariableSelection(VARIABLESELECTION _VARIABLESELECTION)
        {
            switch (_VARIABLESELECTION)
            { 
                case VARIABLESELECTION.MOSTFRAC:
                    return BranchingStrategy.MostFractional(BinaryVarSln);
   
                case VARIABLESELECTION.PSEUDOCOST:
                    return BranchingStrategy.PseudoCost();
            }
            return 1;
        }

        public virtual void SolveNode(string SelectedNodeID, int algorithm)
        {
            this.BinaryVarSln = new List<double>();

            this.RootCplex.SetOut(null);
            this.RootCplex.SetParam(Cplex.Param.RandomSeed, 0);
            SetAlgorithm.Preprocessing_off(this.RootCplex);
            //SetAlgorithm.Cuts_off(this.RootCplex);
            this.RootCplex.SetParam(Cplex.Param.RootAlgorithm, algorithm);
            //this.RootCplex.SetParam(Cplex.Param.Simplex.DGradient, method);
            this.RootCplex.Solve();

            foreach(int idx in this.BinaryVarID)
                this.BinaryVarSln.Add(this.RootCplex.GetValue(this._INumVar[idx]));
            
            TreeNode _TreeNode = this.myTree[0][SelectedNodeID];
            _TreeNode.LPR = this.RootCplex.ObjValue;
            _TreeNode.AllSln = this.RootCplex.GetValues(this._INumVar);
            _TreeNode.BinSln = this.BinaryVarSln;

            this.myTree[0][SelectedNodeID] = _TreeNode;
            System.Console.WriteLine(string.Format("{0} objective value: {1}", SelectedNodeID, _TreeNode.LPR));
        }

        public virtual bool CheckIPFeasible()
        {
            double[] solutions = this.BinaryVarSln.ToArray();

            //Check if all binary variables are integer
            this.isIPFeasible = Array.TrueForAll(solutions, value => { return Math.Abs(value - Math.Round(value)) < 0.0000000001;});
            return this.isIPFeasible;
        }

        public virtual void IdentifyNodeStatus(string SelectedNodeID, VARIABLESELECTION _VARIABLESELECTION)
        {
            //NodeStatus (1:Best IP, 2:Pruned IP, 3:Best LP / Regular LP, 4: Pruned LP, 5:LP Infeasible)
            TreeNode _TreeNode = this.myTree[0][SelectedNodeID];
            UpdateLB(_TreeNode, SelectedNodeID);

            if (this.isIPFeasible)
            {
                if (_TreeNode.LPR <= this.BestUB)
                {
                    this.NodeStatus = 1;

                    //Update Upper Bound & BestNodeID
                    this.BestUB = _TreeNode.LPR;
                    this.BestNodeID = SelectedNodeID;
                }
                else this.NodeStatus = 2;
            }
            else if (_TreeNode.LPR <= this.BestUB)
            {
                if (_TreeNode.LPR >= this.BestLB)
                    this.NodeStatus = 3;

                //Variable Selection and Create two Child Node
                BranchingID = VariableSelection(_VARIABLESELECTION);

                CreateChildNode(BranchingID, SelectedNodeID, BRANCHING.LEFT);
                CreateChildNode(BranchingID, SelectedNodeID, BRANCHING.RIGHT);
            }
            else this.NodeStatus = 4;

            this.BestLB = _LBCandidate.First().Value.LPR;
            this.myTree[NodeStatus].Add(SelectedNodeID, _TreeNode);
            this.myTree[0].Remove(SelectedNodeID);          // Remove incumbent node
        }

        public virtual void UpdateLB(TreeNode _TreeNode, string SelectedNodeID)
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
            string ParentNodeID = SelectedNodeID.Substring(0,SelectedNodeID.Length-1);
            if (_LBCandidate.ContainsKey(ParentNodeID))
            {
                LBCandidate UpdateValue;
                _LBCandidate.TryGetValue(ParentNodeID, out UpdateValue);
                UpdateValue.ChildNodes.Add(SelectedNodeID);

                if (UpdateValue.ChildNodes.Count == 2)
                {
                    _LBCandidate.Remove(ParentNodeID);
                    this.myTree[3].Remove(ParentNodeID);
                }
                else
                    _LBCandidate[ParentNodeID] = UpdateValue;
            }
        }

        public virtual void UpdateLB(string SelectedNodeID)
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
                    this.myTree[3].Remove(ParentNodeID);
                    this.BestLB = _LBCandidate.First().Value.LPR;
                }
                else
                    _LBCandidate[ParentNodeID] = UpdateValue;
            }
        }

        public virtual INumVar getLPNumVar(int idx) { return this.lpmatrix.GetNumVar(idx); }
    }
}
