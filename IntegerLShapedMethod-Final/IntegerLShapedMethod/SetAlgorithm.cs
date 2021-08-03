using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace IntegerLShapedMethod
{
    class SetAlgorithm
    {
        public static void Cuts_off(ILOG.CPLEX.Cplex i)
        {
            //cuts off
            i.SetParam(ILOG.CPLEX.Cplex.Param.MIP.Cuts.Covers, -1);
            i.SetParam(ILOG.CPLEX.Cplex.Param.MIP.Cuts.Cliques, -1);
            i.SetParam(ILOG.CPLEX.Cplex.Param.MIP.Cuts.Disjunctive, -1);
            i.SetParam(ILOG.CPLEX.Cplex.Param.MIP.Cuts.FlowCovers, -1);
            i.SetParam(ILOG.CPLEX.Cplex.Param.MIP.Cuts.Gomory, -1);
            i.SetParam(ILOG.CPLEX.Cplex.Param.MIP.Cuts.GUBCovers, -1);
            i.SetParam(ILOG.CPLEX.Cplex.Param.MIP.Cuts.Implied, -1);
            i.SetParam(ILOG.CPLEX.Cplex.Param.MIP.Cuts.LiftProj, -1);
            i.SetParam(ILOG.CPLEX.Cplex.Param.MIP.Cuts.MCFCut, -1);
            i.SetParam(ILOG.CPLEX.Cplex.Param.MIP.Cuts.MIRCut, -1);
            i.SetParam(ILOG.CPLEX.Cplex.Param.MIP.Cuts.PathCut, -1);
            i.SetParam(ILOG.CPLEX.Cplex.Param.MIP.Cuts.ZeroHalfCut, -1);
        }

        public static void Preprocessing_off(ILOG.CPLEX.Cplex i)
        {
            //preprocessing off
            i.SetParam(ILOG.CPLEX.Cplex.Param.Preprocessing.BoundStrength, 0);
            i.SetParam(ILOG.CPLEX.Cplex.Param.Preprocessing.CoeffReduce, 0);
            i.SetParam(ILOG.CPLEX.Cplex.Param.Preprocessing.Dependency, 0);
            i.SetParam(ILOG.CPLEX.Cplex.Param.Preprocessing.Dual, -1);
            i.SetParam(ILOG.CPLEX.Cplex.Param.Preprocessing.Fill, 0);
            i.SetParam(ILOG.CPLEX.Cplex.Param.Preprocessing.Linear, 0);
            i.SetParam(ILOG.CPLEX.Cplex.Param.Preprocessing.NumPass, 0);
            i.SetParam(ILOG.CPLEX.Cplex.Param.Preprocessing.Presolve, false);
            i.SetParam(ILOG.CPLEX.Cplex.Param.Preprocessing.QPMakePSD, false);
            i.SetParam(ILOG.CPLEX.Cplex.Param.Preprocessing.Reduce, 0);
            i.SetParam(ILOG.CPLEX.Cplex.Param.Preprocessing.Relax, 0);
            i.SetParam(ILOG.CPLEX.Cplex.Param.Preprocessing.RepeatPresolve, 0);
            i.SetParam(ILOG.CPLEX.Cplex.Param.Preprocessing.Symmetry, 0);
        }

        public static void Turn_off(ILOG.CPLEX.Cplex i, int k)
        {
            if (k == 0)
            {
                Cuts_off(i);
            }
            else if (k == 1)
            {
                Preprocessing_off(i);
            }
            else if (k == 2)
            {
                Cuts_off(i);
                Preprocessing_off(i);
            }
        }
    }
}
