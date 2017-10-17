/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */



#if defined(SOPLEX_WITH_QPOASES) || defined(SOPLEX_WITH_QORE)
#ifndef SOPLEX_LEGACY

#include <iostream>
#include <assert.h>
#ifdef SOPLEX_WITH_QPOASES
#include <qpOASES.hpp>
#endif
#ifdef SOPLEX_WITH_QORE
extern "C" {
#include "qpsolver.h"
}
#endif
#include <vector>

#include "soplex.h"
#include "statistics.h"
#include "slufactor_rational.h"

#ifdef SOPLEX_WITH_QPOASES
USING_NAMESPACE_QPOASES
#endif

namespace soplex
{

#ifdef SOPLEX_WITH_QPOASES
   /// copies QP from soplex into qpOASES arrays
   void SoPlex::_copyQPintoqpOASES(
      const int var_num,
      const int constr_num)
   {
      assert( var_num == _realLP->nCols() );
      assert( constr_num == _realLP->nRows() );

      MSG_INFO1( spxout, spxout << "Preparing QP data . . .\n" );

      // calculate the sparcity
      const double sparsity_A = static_cast<double>(_realLP->nNzos())/var_num/constr_num;
      const double sparsity_H = static_cast<double>(_realLP->nNzos_qp())/var_num/var_num;

      MSG_DEBUG( std::cout << "Sparsity of constraints and objective: " << sparsity_A << ", " << sparsity_H << "\n" );

      // decide for sparce qpOASES call
      if( sparsity_A + sparsity_H < 2 * realParam(SoPlex::QP_SPARSITY_THRESHOLD) )
#ifdef SOPLEX_WITH_QPOASES_SPARSE
         _useSparseQP = 1;
      else
#else
         MSG_WARNING( spxout, spxout << "Warning: compile with QPOASES=sparse to enable sparse solving with qpOASES, use dense mode instead.\n" );
#endif
      {
         _useSparseQP = 0;
      }

      // allocate qpOASES variables
#ifdef SOPLEX_WITH_QPOASES_SPARSE
      if( _useSparseQP == 1 )
      {
         spx_alloc(H_jc, var_num+1);
         spx_alloc(H_ir, _realLP->nNzos_qp());
         spx_alloc(H_val, _realLP->nNzos_qp());
         spx_alloc(A_jc, var_num+1);
         spx_alloc(A_ir, _realLP->nNzos());
         spx_alloc(A_val, _realLP->nNzos());
      }
      else
#endif
      {
         spx_alloc(H, var_num*var_num);
         spx_alloc(A, constr_num*var_num);
      }
      spx_alloc(g, var_num);
      spx_alloc(lb, var_num);
      spx_alloc(ub, var_num);
      spx_alloc(lbA, constr_num);
      spx_alloc(ubA, constr_num);

      // copy constraint matrix, sparse  or dense
#ifdef SOPLEX_WITH_QPOASES_SPARSE
      if( _useSparseQP == 1 )
      {
	 int k = 0;
         for( int i = 0; i < var_num; i++ )
         {
            A_jc[i] = k;
	    const SVectorBase<Real>& col = _realLP->colVector(i);
	    for( int j = 0; j < col.size(); j++ )
	    {
	       A_ir[k] = col.index(j);
	       A_val[k] = col.value(j);
               k++;
            }
	 }
         A_jc[var_num] = k;
         A_spar = new SparseMatrix(constr_num, var_num, A_ir, A_jc, A_val);
      }
      else
#endif
      {
         for( int i = 0; i < constr_num; i++ )
         {
	    const SVectorBase<Real>& row = _realLP->rowVector(i);
	    for( int j = 0; j < var_num; j++ )
	       A[i*var_num + j] = row[j];
	 }
      }


      // copy quadratic objective matrix, sparse or dense
#ifdef SOPLEX_WITH_QPOASES_SPARSE
      if( _useSparseQP == 1 )
      {
	 int n;
	 int k = 0;
         for( int i = 0; i < var_num; i++ )
         {
            H_jc[i] = k;
            for( int j = 0; j < i; j++ )
            {
               const SVectorBase<Real>& rowcol = _realLP->rowcolVector(j);
               n = rowcol.number(i);
               if( n >= 0 )
	       {
                  H_ir[k] = j;
	          H_val[k] = rowcol.value(n);
                  k++;
               }
	    }
	    const SVectorBase<Real>& rowcol = _realLP->rowcolVector(i);
            if( rowcol.size() > 0 )
	    {
	       if( rowcol.index(0) != i )
	       {
                  H_ir[k] = i;
                  H_val[k] = 0;
                  k++;
               }
               for( int j = 0; j < rowcol.size(); j++ )
	       {
	          H_ir[k] = rowcol.index(j);
	          H_val[k] = rowcol.value(j);
                  k++;
               }
            }
            else
	    {
               H_ir[k] = i;
	       H_val[k] = 0;
               k++;
            }
	 }
         H_jc[var_num] = k;
         H_spar = new SymSparseMat(var_num, var_num, H_ir, H_jc, H_val);
         H_spar->createDiagInfo();
      }
      else
#endif
      {
         for( int i = 0; i < var_num; i++ )
         {
            const SVectorBase<Real>& rowcol = _realLP->rowcolVector(i);
            H[i*var_num + i] = rowcol[i]; // diagonal entry
            for( int j = i+1; j < var_num; j++ ) // symmetric non diagonal entries
            {
               H[i*var_num + j] = rowcol[j];
               H[j*var_num + i] = rowcol[j];
            }
         }
      }

      // copy linear objective and bounds
      for( int i = 0; i < var_num; i++ )
      {
         g[i] = _realLP->obj(i);
         if( _realLP->lower(i) > -qpOASES::INFTY )
            lb[i] = _realLP->lower(i);
         else
            lb[i] = -qpOASES::INFTY;
         if( _realLP->upper(i) < qpOASES::INFTY )
            ub[i] = _realLP->upper(i);
         else
            ub[i] = qpOASES::INFTY;
      }

      // copy sides
      for( int i = 0; i < constr_num; i++ )
      {
         if( _realLP->lhs(i) > -qpOASES::INFTY )
	    lbA[i] = _realLP->lhs(i);
	 else
	    lbA[i] = -qpOASES::INFTY;
	 if( _realLP->rhs(i) < qpOASES::INFTY )
	    ubA[i] = _realLP->rhs(i);
	 else
	    ubA[i] = qpOASES::INFTY;
      }
   }
#endif



#ifdef SOPLEX_WITH_QORE
   /// copies QP from soplex into QORE arrays
   void SoPlex::_copyQPintoQORE(
      const int var_num,
      const int constr_num)
   {
      assert( var_num == _realLP->nCols() );
      assert( constr_num == _realLP->nRows() );

      MSG_INFO1( spxout, spxout << "Preparing QP data . . .\n" );

      // allocate QORE QP variables
      spx_alloc(Hcp, var_num+1);
      spx_alloc(Hri, _realLP->nNzos_qp());
      spx_alloc(Hnz, _realLP->nNzos_qp());
      spx_alloc(Atcp, constr_num+1);
      spx_alloc(Atri, _realLP->nNzos());
      spx_alloc(Atnz, _realLP->nNzos());
      spx_alloc(g, var_num);
      spx_alloc(lbt, var_num+constr_num);
      spx_alloc(ubt, var_num+constr_num);

      // copy constraint matrix, sparse (QORE)
      int k = 0;
      for( int i = 0; i < constr_num; i++ )
      {
	 Atcp[i] = k;
	 const SVectorBase<Real>& row = _realLP->rowVector(i);
	 for( int j = 0; j < row.size(); j++ )
	 {
	    Atri[k] = row.index(j);
	    Atnz[k] = row.value(j);
	    k++;
	 }
      }
      Atcp[constr_num] = k;

      // copy quadratic objective matrix, sparse (QORE)
      int n;
      k = 0;
      for( int i = 0; i < var_num; i++ )
      {
	 Hcp[i] = k;
	 for( int j = 0; j < i; j++ )
	 {
	    const SVectorBase<Real>& rowcol = _realLP->rowcolVector(j);
	    n = rowcol.number(i);
	    if( n >= 0 )
	    {
	       Hri[k] = j;
	       Hnz[k] = rowcol.value(n);
	       k++;
	    }
	 }
	 const SVectorBase<Real>& rowcol = _realLP->rowcolVector(i);
	 if( rowcol.size() > 0 )
	 {
	    if( rowcol.index(0) != i )
	    {
	       Hri[k] = i;
	       Hnz[k] = 0;
	       k++;
	    }
	    for( int j = 0; j < rowcol.size(); j++ )
	    {
	       Hri[k] = rowcol.index(j);
	       Hnz[k] = rowcol.value(j);
	       k++;
	    }
	 }
	 else
	 {
	    Hri[k] = i;
	    Hnz[k] = 0;
	    k++;
	 }
      }
      Hcp[var_num] = k;

      // copy linear objective and bounds
      for( int i = 0; i < var_num; i++ )
      {
	 g[i] = _realLP->obj(i);
	 if( _realLP->lower(i) > -INFINITY )
	    lbt[i] = _realLP->lower(i);
	 else
	    lbt[i] = -INFINITY;
	 if( _realLP->upper(i) < INFINITY )
	    ubt[i] = _realLP->upper(i);
	 else
	    ubt[i] = INFINITY;
      }

      // copy sides
      for( int i = var_num; i < constr_num+var_num; i++ )
      {
	 if( _realLP->lhs(i-var_num) > -INFINITY )
	    lbt[i] = _realLP->lhs(i-var_num);
	 else
	    lbt[i] = -INFINITY;
	 if( _realLP->rhs(i-var_num) < INFINITY )
	    ubt[i] = _realLP->rhs(i-var_num);
	 else
	    ubt[i] = INFINITY;
      }
   }
#endif



#ifdef SOPLEX_WITH_QPOASES
   /// Copies the solution of qpOASES back into soplex
   void SoPlex::_copySolintoSPX(
      const QProblem& qp,
      VectorReal& primal,
      VectorReal& dual,
      const int var_num,
      const int constr_num)
   {

      assert( var_num == _realLP->nCols() );
      assert( constr_num == _realLP->nRows() );

      // get solution
      real_t* xOpt = 0;
      real_t* yOpt = 0;

      spx_alloc(xOpt, var_num);
      spx_alloc(yOpt, var_num + constr_num);

      qp.getPrimalSolution( xOpt );
      qp.getDualSolution( yOpt );

      // store primal solution
      _solReal._hasPrimal = true;
      _solReal._primal.reDim(var_num);
      for( int i = 0; i < var_num; i++ )
      {
         primal[i] = xOpt[i];
         _solReal._primal[i] = xOpt[i];
      }

      // store dual solution
      _solReal._hasDual = true;
      _solReal._dual.reDim(constr_num);
      _solReal._redCost.reDim(var_num);
      for( int i = 0; i < var_num; i++ )
         _solReal._redCost[i] = yOpt[i];
      for( int i = 0; i < constr_num; i++ )
      {
         dual[i] = yOpt[var_num + i];
         _solReal._dual[i] = yOpt[var_num + i];
      }

      // get objective value
      real_t ObjVal = qp.getObjVal() + _realLP->objOffset();
      _solReal._primalObjVal = ObjVal;
      _solReal._dualObjVal = ObjVal;

      // get active bounds and constraints
      _hasBasis = true;
      Bounds Bounds;
      Constraints Constraints;
      qp.getBounds( Bounds );
      qp.getConstraints( Constraints );
      SubjectToStatus Status;
      _basisStatusCols.reSize( var_num );
      for( int i = 0; i < var_num; i++ )
      {
         Status = Bounds.getStatus( i );
         if( Status == ST_LOWER )
         {
            _basisStatusCols[i] = SPxSolver::ON_LOWER;
         }
         else if( Status == ST_INACTIVE )
         {
            _basisStatusCols[i] = SPxSolver::BASIC;
         }
	 else if( Status == ST_UPPER )
	 {
	    _basisStatusCols[i] = SPxSolver::ON_UPPER;
	 }
	 else
	 {
	    _basisStatusCols[i] = SPxSolver::UNDEFINED;
	 }
      }
      _basisStatusRows.reSize( constr_num );
      for( int i = 0; i < constr_num; i++ )
      {
	 Status = Constraints.getStatus( i );
	 if( Status == ST_LOWER )
	 {
	    _basisStatusRows[i] = SPxSolver::ON_LOWER;
	 }
	 else if( Status == ST_INACTIVE )
	 {
	    _basisStatusRows[i] = SPxSolver::BASIC;
	 }
	 else if( Status == ST_UPPER )
	 {
	    _basisStatusRows[i] = SPxSolver::ON_UPPER;
	 }
	 else
	 {
	    _basisStatusRows[i] = SPxSolver::UNDEFINED;
	 }
      }

      spx_free(xOpt);
      spx_free(yOpt);
   }
#endif



#ifdef SOPLEX_WITH_QORE
   /// Copies the solution of QORE back into soplex
   void SoPlex::_copySolintoSPX(
      QoreProblem* pproblem,
      VectorReal& primal,
      VectorReal& dual,
      const int var_num,
      const int constr_num)
   {

      assert( var_num == _realLP->nCols() );
      assert( constr_num == _realLP->nRows() );

      // get solution
      double* xOpt = 0;
      double* yOpt = 0;
      int* workingset = 0;

      // internally QORE adds slack variables as well and writes them to xOpt, we don't use them (yet?)
      spx_alloc(xOpt, var_num + constr_num);
      spx_alloc(yOpt, var_num + constr_num);
      spx_alloc(workingset, var_num + constr_num);

#ifndef NDEBUG
      qp_int rv = QPSOLVER_OK;
      rv = QPGetDblVector( pproblem, "primalsol", xOpt );
      assert( rv == QPSOLVER_OK );
      rv = QPGetDblVector( pproblem, "dualsol", yOpt );
      assert( rv == QPSOLVER_OK );
#else
      (void) QPGetDblVector( pproblem, "primalsol", xOpt );
      (void) QPGetDblVector( pproblem, "dualsol", yOpt );
#endif

      // store primal solution
      _solReal._hasPrimal = true;
      _solReal._primal.reDim(var_num);
      for( int i = 0; i < var_num; i++ )
      {
         primal[i] = xOpt[i];
         _solReal._primal[i] = xOpt[i];
      }

      // store dual solution
      _solReal._hasDual = true;
      _solReal._dual.reDim(constr_num);
      _solReal._redCost.reDim(var_num);
      for( int i = 0; i < var_num; i++ )
         _solReal._redCost[i] = yOpt[i];
      for( int i = 0; i < constr_num; i++ )
      {
         dual[i] = yOpt[var_num + i];
         _solReal._dual[i] = yOpt[var_num + i];
      }

      // get objective value
      //real_t ObjVal = qp.getObjVal() + _realLP->objOffset();
      //_solReal._primalObjVal = ObjVal;
      //_solReal._dualObjVal = ObjVal;

      // get active bounds and constraints
      _hasBasis = true;
#ifndef NDEBUG
      rv = QPGetIntVector(pproblem, "workingset", workingset);
      assert( rv == QPSOLVER_OK );
#else
      (void) QPGetIntVector(pproblem, "workingset", workingset);
#endif
      _basisStatusCols.reSize( var_num );
      for( int i = 0; i < var_num; i++ )
      {
	if( workingset[i] == 1 )
         {
            _basisStatusCols[i] = SPxSolver::ON_LOWER;
         }
         else if( workingset[i] == 0 )
         {
            _basisStatusCols[i] = SPxSolver::BASIC;
         }
	 else if( workingset[i] == -1 )
	 {
	    _basisStatusCols[i] = SPxSolver::ON_UPPER;
	 }
	 else
	 {
	    _basisStatusCols[i] = SPxSolver::UNDEFINED;
	 }
      }
      _basisStatusRows.reSize( constr_num );
      for( int i = var_num; i < constr_num+var_num; i++ )
      {
	 if( workingset[i] == 1 )
	 {
	    _basisStatusRows[i-var_num] = SPxSolver::ON_LOWER;
	 }
	 else if( workingset[i] == 0 )
	 {
	    _basisStatusRows[i-var_num] = SPxSolver::BASIC;
	 }
	 else if( workingset[i] == -1 )
	 {
	    _basisStatusRows[i-var_num] = SPxSolver::ON_UPPER;
	 }
	 else
	 {
	    _basisStatusRows[i-var_num] = SPxSolver::UNDEFINED;
	 }
      }

      spx_free(xOpt);
      spx_free(yOpt);
      spx_free(workingset);
   }
#endif



   /// solves QP with iterative refinement
   SPxSolver::Status SoPlex::solveQP()
   {
      bool infeasibilityNotCertified = false;
      bool unboundednessNotCertified = false;

      // clear solving statistics
      _statistics->clearSolvingData();

      // ensure rational LP is available
      if( intParam(SoPlex::SYNCMODE) == SYNCMODE_ONLYREAL )
         _syncLPRational(true);

      // start timing
      _statistics->solvingTime->start();
      _statistics->preprocessingTime->start();

      // remember that last solve was rational
      _lastSolveMode = SOLVEMODE_RATIONAL;

      bool primalFeasible = false;
      bool dualFeasible = false;
      bool infeasible = false;
      bool unbounded = false;
      bool stoppedTime = false;
      bool stoppedIter = false;
      bool error = false;

      _realLP = 0;
      spx_alloc(_realLP);
      _realLP = new (_realLP) SPxLPReal(_solver);
      _isRealLPLoaded = false;

      // store objective, bounds, and sides of real LP in case they will be modified during iterative refinement
      _storeLPReal();

      // stop timing of prerpocessing
      _statistics->preprocessingTime->stop();

      // introduce slack variables to transform inequality constraints into equations
      if( boolParam(SoPlex::EQTRANS) == true )
      {
         _transformEquality();
      }
      else
         assert( intParam(SoPlex::SOLVEMODE) == SOLVEMODE_REAL || (intParam(SoPlex::SOLVEMODE) == SOLVEMODE_AUTO && GE(realParam(SoPlex::FEASTOL), 1e-15) && GE(realParam(SoPlex::OPTTOL), 1e-15)) );

      // clear rational LU factorization of basis system
      // _rationalLUSolver.clear();

      // the solution is no longer valid
      _invalidateSolution();

      // solve problem with iterative refinement (no testing for infeas. or unbound.)
      _performOptIRQP(_solRational, !unboundednessNotCertified, !infeasibilityNotCertified, 0,
		      primalFeasible, dualFeasible, infeasible, unbounded, stoppedTime, stoppedIter, error);

      // case: an unrecoverable error occured
      if( error )
      {
         _status = SPxSolver::ERROR;
      }
      // case: stopped due to some limit
      else if( stoppedTime )
      {
         _status = SPxSolver::ABORT_TIME;
      }
      else if( stoppedIter )
      {
         _status = SPxSolver::ABORT_ITER;
      }


      // restore original problem
      if( boolParam(SoPlex::EQTRANS) == true )
      {
         _untransformEquality(_solRational);
      }
      else
         assert( intParam(SoPlex::SOLVEMODE) == SOLVEMODE_REAL || (intParam(SoPlex::SOLVEMODE) == SOLVEMODE_AUTO && GE(realParam(SoPlex::FEASTOL), 1e-12) && GE(realParam(SoPlex::OPTTOL), 1e-12)) );

      // restore objective, bounds, and sides of real LP in case they have been modified during iterative refinement
      _restoreLPReal();

      // stop timing
      _statistics->solvingTime->stop();

      MSG_INFO1( spxout, spxout << "\n";
         printShortStatistics(spxout.getStream(SPxOut::INFO1));
         spxout << "\n" );

      // free memory
      if( intParam(SoPlex::QP_SOLVER) == SOLVER_QPOASES )
      {
#ifdef SOPLEX_WITH_QPOASES
#ifdef SOPLEX_WITH_QPOASES_SPARSE
	 if( _useSparseQP == 1 )
	 {
	    spx_free(H_jc);
	    spx_free(H_ir);
	    spx_free(H_val);
	    delete H_spar;
	    spx_free(A_jc);
	    spx_free(A_ir);
	    spx_free(A_val);
	    delete A_spar;
	 }
	 else
#endif
	 {
	    spx_free(H);
	    spx_free(A);
	 }
	 spx_free(g);
	 spx_free(lb);
	 spx_free(ub);
	 spx_free(lbA);
	 spx_free(ubA);
#else
         MSG_ERROR( std::cerr << "Error: Cannot free QPOASES interface variables when not linked to QPOASES.\n" );
         return status();
#endif
      }
      else if( intParam(SoPlex::QP_SOLVER) == SOLVER_QORE )
      {
#ifdef SOPLEX_WITH_QORE
	 spx_free(g);
	 spx_free(lbt);
	 spx_free(ubt);
#else
         MSG_ERROR( std::cerr << "Error: Cannot free QORE interface variables when not linked to QORE.\n" );
         return status();
#endif
      }
      return status();
   }



   /// solves current QP problem with iterative refinement
   void SoPlex::_performOptIRQP(SolRational& sol,
      bool acceptUnbounded,
      bool acceptInfeasible,
      int minRounds,
      bool& primalFeasible,
      bool& dualFeasible,
      bool& infeasible,
      bool& unbounded,
      bool& stoppedTime,
      bool& stoppedIter,
      bool& error)
   {

      // setting up QProblem object/struct
#ifdef SOPLEX_WITH_QORE
      QoreProblem * pproblem = 0;
#endif
#ifdef SOPLEX_WITH_QPOASES
      QProblem soplex_ir_qp;
#endif
      if( intParam(SoPlex::QP_SOLVER) == SOLVER_QPOASES )
      {
#ifdef SOPLEX_WITH_QPOASES
         soplex_ir_qp = QProblem(_realLP->nCols(), _realLP->nRows());
#else
         MSG_ERROR( std::cerr << "Error: Cannot use QPOASES when not linked to it.\n" );
         error = true;
         return;
#endif
      }
      else if( intParam(SoPlex::QP_SOLVER) == SOLVER_QORE )
      {
#ifdef SOPLEX_WITH_QORE
#ifndef NDEBUG
         qp_int rv = QPSOLVER_OK;
         rv = QPNew( &pproblem, _realLP->nCols(), _realLP->nRows(), _realLP->nNzos(), _realLP->nNzos_qp() );
         assert( rv == QPSOLVER_OK );
         assert( pproblem != 0 );
#else
         (void) QPNew( &pproblem, _realLP->nCols(), _realLP->nRows(), _realLP->nNzos(), _realLP->nNzos_qp() );
#endif
#else
         MSG_ERROR( std::cerr << "Error: Cannot use QORE when not linked to it.\n" );
         error = true;
         return;
#endif
      }

      // start rational solving timing
      _statistics->rationalTime->start();

      primalFeasible = false;
      dualFeasible = false;
      infeasible = false;
      unbounded = false;
      stoppedTime = false;
      stoppedIter = false;
      error = false;

      // declare vectors and variables
      SPxSolver::Status result = SPxSolver::UNKNOWN;

      _modLower.reDim(numColsRational(), false);
      _modUpper.reDim(numColsRational(), false);
      _modLhs.reDim(numRowsRational(), false);
      _modRhs.reDim(numRowsRational(), false);
      _modObj.reDim(numColsRational(), false);

      DVectorReal primalReal(numColsRational());
      DVectorReal dualReal(numRowsRational());

      Rational temporary_rational;
      Rational boundsViolation;
      Rational sideViolation;
      Rational redCostViolation;
      Rational dualViolation;
      Rational pdScale;
      Rational maxScale;
      Rational minScale;

      // solve original QP
      MSG_INFO1( spxout, spxout << "Initial floating-point solve . . .\n" );

      _statistics->rationalTime->stop();
      if( intParam(SoPlex::QP_SOLVER) == SOLVER_QPOASES )
      {
#ifdef SOPLEX_WITH_QPOASES
         result = _solveRealQP(soplex_ir_qp, acceptUnbounded, acceptInfeasible, primalReal, dualReal, _basisStatusRows, _basisStatusCols, _hasBasis);
#endif
      }
      else if( intParam(SoPlex::QP_SOLVER) == SOLVER_QORE )
      {
#ifdef SOPLEX_WITH_QORE
         result = _solveRealQP(pproblem, acceptUnbounded, acceptInfeasible, primalReal, dualReal, _basisStatusRows, _basisStatusCols, _hasBasis);
#endif
      }

      // evaluate result
      switch( result )
      {
      case SPxSolver::OPTIMAL:
         MSG_INFO1( spxout, spxout << "Floating-point optimal.\n" );
         break;
      case SPxSolver::INFEASIBLE:
         MSG_INFO1( spxout, spxout << "Floating-point infeasible.\n" );
         // the floating-point solve returns a Farkas ray if and only if the simplifier was not used, which is exactly
         // the case when a basis could be returned
         if( _hasBasis )
         {
            sol._dualFarkas = dualReal;
            sol._hasDualFarkas = true;
         }
         else
            sol._hasDualFarkas = false;
         infeasible = true;
         return;
      case SPxSolver::UNBOUNDED:
         MSG_INFO1( spxout, spxout << "Floating-point unbounded.\n" );
         unbounded = true;
         return;
      case SPxSolver::ABORT_TIME:
         stoppedTime = true;
         return;
      case SPxSolver::ABORT_ITER:
         stoppedIter = true;
         return;
      default:
         error = true;
         return;
      }

      // start rational timing again
      _statistics->rationalTime->start();

      // store floating-point solution of original LP as current rational solution and ensure that solution vectors have
      // right dimension; ensure that solution is aligned with basis
      _realSolToRational( sol, primalReal, dualReal );

      // initial scaling factor is one
      pdScale = Rational::POSONE;

      // control progress
      Rational maxViolation;
      Rational bestViolation = _rationalPosInfty;
      const Rational violationImprovementFactor = 0.9;
      const Rational errorCorrectionFactor = 1.1;
      Rational errorCorrection = 2;

      // store basis status in case solving modified problem failed
      DataArray< SPxSolver::VarStatus > basisStatusRowsFirst;
      DataArray< SPxSolver::VarStatus > basisStatusColsFirst;

      // refinement loop
      const bool maximizing = (intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MAXIMIZE);
      const int maxDimRational = numColsRational() > numRowsRational() ? numColsRational() : numRowsRational();
      SolRational factorSol;
      bool factorSolNewBasis = true;
      int lastStallRefinements = 0;
      do
      {
         // decrement minRounds counter
         minRounds--;

         MSG_DEBUG( std::cout << "Computing primal violations.\n" );

	 // compute primal violation
	 _primViolation( sol, boundsViolation, sideViolation );

         MSG_DEBUG( std::cout << "Computing dual violations.\n" );

	 // compute dual violation
	 _dualViolation( sol, redCostViolation, dualViolation, maximizing );

         // output violations; the reduced cost violations for artificially introduced slack columns are actually violations of the dual multipliers
         MSG_INFO1( spxout, spxout
            << "Max. bound violation = " << rationalToString(boundsViolation) << "\n"
            << "Max. row violation = " << rationalToString(sideViolation) << "\n"
            << "Max. reduced cost violation = " << rationalToString(redCostViolation) << "\n"
            << "Max. dual violation = " << rationalToString(dualViolation) << "\n" );

         MSG_DEBUG( spxout
            << std::fixed << std::setprecision(2) << std::setw(10)
            << "Progress table: "
            << std::setw(10) << _statistics->refinements << " & "
            << std::setw(10) << _statistics->iterations << " & "
            << std::setw(10) << _statistics->solvingTime->time() << " & "
            << std::setw(10) << _statistics->rationalTime->time() << " & "
            << std::setw(10) << rationalToString(boundsViolation > sideViolation ? boundsViolation : sideViolation, 2) << " & "
            << std::setw(10) << rationalToString(redCostViolation > dualViolation ? redCostViolation : dualViolation, 2) << "\n");

         // terminate if tolerances are satisfied
         primalFeasible = (boundsViolation <= _rationalFeastol && sideViolation <= _rationalFeastol);
         dualFeasible = (redCostViolation <= _rationalOpttol && dualViolation <= _rationalOpttol);
         if( primalFeasible && dualFeasible )
         {
            if( minRounds < 0 )
            {
               MSG_INFO1( spxout, spxout << "Tolerances reached.\n" );
               _hasSolRational = true;
               break;
            }
            else
            {
               MSG_INFO1( spxout, spxout << "Tolerances reached but minRounds forcing additional refinement rounds.\n" );
            }
         }

         // terminate if some limit is reached
         if( _isSolveStopped(stoppedTime, stoppedIter) )
         {
            break;
         }

         // check progress
         maxViolation = boundsViolation;
         if( sideViolation > maxViolation )
            maxViolation = sideViolation;
         if( redCostViolation > maxViolation )
            maxViolation = redCostViolation;
         if( dualViolation > maxViolation )
            maxViolation = dualViolation;
         bestViolation *= violationImprovementFactor;
         if( maxViolation > bestViolation )
         {
            MSG_INFO2( spxout, spxout << "Failed to reduce violation significantly.\n" );
         }
         else
            bestViolation = maxViolation;

         // decide whether to perform rational reconstruction and/or factorization
         bool performRatfac = boolParam(SoPlex::RATFAC)
            && lastStallRefinements >= intParam(SoPlex::RATFAC_MINSTALLS) && _hasBasis && factorSolNewBasis;

	 // ensure that artificial slack columns are basic and inequality constraints are nonbasic; otherwise we may end
	 // up with dual violation on inequality constraints after removing the slack columns
	 if( _slackCols.num() > 0 && _hasBasis )
         {
            int numOrigCols = numColsRational() - _slackCols.num();
            for( int i = 0; i < _slackCols.num(); i++ )
            {
	       int row = _slackCols.colVector(i).index(0);
	       int col = numOrigCols + i;

	       assert(row >= 0);
	       assert(row < numRowsRational());

	       if( _basisStatusRows[row] == SPxSolver::BASIC && _basisStatusCols[col] != SPxSolver::BASIC )
	       {
	          _basisStatusRows[row] = _basisStatusCols[col];
	          _basisStatusCols[col] = SPxSolver::BASIC;
	       }
	    }
	 }

         // solve basis systems exactly
         if( performRatfac && maxViolation > 0 )
         {
            MSG_INFO1( spxout, spxout << "Performing rational factorization . . .\n" );

            bool optimal;
            _factorizeColumnRowRational_qp(sol, _basisStatusRows, _basisStatusCols, stoppedTime, stoppedIter, error, optimal);
            factorSolNewBasis = false;

            if( stoppedTime )
            {
               MSG_INFO1( spxout, spxout << "Stopped rational factorization.\n" );
            }
            else if( error )
            {
               // message was already printed; reset error flag and continue without factorization
               error = false;
            }
            else if( optimal )
            {
               MSG_INFO1( spxout, spxout << "Tolerances reached.\n" );
               _hasSolRational = true;
               primalFeasible = true;
               dualFeasible = true;
               break;
            }
            else if( boolParam(SoPlex::RATFACJUMP) )
            {
               MSG_INFO1( spxout, spxout << "Jumping to exact basic solution.\n" );
               minRounds++;
               continue;
            }
         }

         // start refinement

         // compute primal/dual scaling factor; limit increase in scaling by tolerance used in floating point solve
         minScale = pdScale;
         maxScale = pdScale;
         maxScale *= _rationalMaxscaleincr;

         pdScale = boundsViolation > sideViolation ? boundsViolation : sideViolation;
         if( pdScale < redCostViolation )
            pdScale = redCostViolation;
         assert(pdScale >= 0);

         if( pdScale > 0 )
         {
            pdScale.invert();
            if( pdScale > maxScale )
               pdScale = maxScale;
         }
         else
            pdScale = maxScale;

         if( boolParam(SoPlex::POWERSCALING) )
            pdScale.powRound();

         // variable (decreasing) scaling in refinement if floating-point solve fails
         for( int i = 0; i <= intParam(SoPlex::QP_REF_MAX_NUM_BACKSTEP); i++ )
         {
	    // print scaling factor
            MSG_INFO1( spxout, spxout << "Primal/dual scaling factor = " << rationalToString(pdScale) << "\n");

	    // apply primal scaleing factor
	    _appPrimScale( pdScale );

            // apply dual scaleing factor
            _appDualScale( sol, pdScale );

            MSG_INFO1( spxout, spxout << "Refined floating-point solve . . .\n" );

	    // solve modified problem
	    int prevIterations = _statistics->iterations;
	    _statistics->rationalTime->stop();
            if( intParam(SoPlex::QP_SOLVER) == SOLVER_QPOASES )
	    {
#ifdef SOPLEX_WITH_QPOASES
               result = _resolveRealQP(soplex_ir_qp, acceptUnbounded, acceptInfeasible, primalReal, dualReal, _basisStatusRows, _basisStatusCols, _hasBasis, pdScale > 1e20);
#endif
	    }
	    else if( intParam(SoPlex::QP_SOLVER) == SOLVER_QORE )
	    {
#ifdef SOPLEX_WITH_QORE
               result = _resolveRealQP(pproblem, acceptUnbounded, acceptInfeasible, primalReal, dualReal, _basisStatusRows, _basisStatusCols, _hasBasis, pdScale > 1e20);
#endif
	    }

	    // count refinements and remember whether we moved to a new basis
	    _statistics->refinements++;
	    if( _statistics->iterations <= prevIterations )
	    {
	       lastStallRefinements++;
	       _statistics->stallRefinements++;
	    }
	    else
	    {
	       factorSolNewBasis = true;
	       lastStallRefinements = 0;
	       _statistics->pivotRefinements = _statistics->refinements;
	    }

	    // If scaled problem was not solved we reduce the scaling factor, otherwise we procede with iterative refinement
            temporary_rational = pdScale*1e-2;
            if( result != SPxSolver::OPTIMAL && temporary_rational > minScale  )
	    {
               pdScale = temporary_rational;
	       _statistics->backstepping++;
	    }
            else
               break;
         }

         // evaluate result; if modified problem was not solved to optimality, stop refinement
         switch( result )
         {
         case SPxSolver::OPTIMAL:
            MSG_INFO1( spxout, spxout << "Floating-point optimal.\n" );
            break;
         case SPxSolver::INFEASIBLE:
            MSG_INFO1( spxout, spxout << "Floating-point infeasible.\n" );
            sol._dualFarkas = dualReal;
            sol._hasDualFarkas = true;
            infeasible = true;
            _solver.clearRowObjs();
            return;
         case SPxSolver::UNBOUNDED:
            MSG_INFO1( spxout, spxout << "Floating-point unbounded.\n" );
            unbounded = true;
            _solver.clearRowObjs();
            return;
         case SPxSolver::ABORT_TIME:
            stoppedTime = true;
            return;
         case SPxSolver::ABORT_ITER:
            stoppedIter = true;
            _solver.clearRowObjs();
            return;
         default:
            error = true;
            _solver.clearRowObjs();
            return;
         }

         _statistics->rationalTime->start();

         // correct primal solution and align with basis
         MSG_DEBUG( std::cout << "Correcting primal solution.\n" );

	 // do correct primal solution
	 (void) _corrPrimSol( sol, pdScale, maxDimRational, primalReal );

         // update or recompute slacks depending on which looks faster
         _rationalLP->computePrimalActivity(sol._primal, sol._slacks);
         const int numCorrectedPrimals = _primalDualDiff.size();

         // correct dual solution and align with basis
         MSG_DEBUG( std::cout << "Correcting dual solution.\n" );

#ifndef NDEBUG
	 _dbgRedAndDualViolation( sol, dualReal, primalReal, maximizing, pdScale);
#endif

	 // do correct dual solution
	 (void) _corrDualSol( sol, pdScale, maxDimRational, dualReal );

         // update or recompute reduced cost values depending on which looks faster; adding one to the length of the
         // dual vector accounts for the objective function vector
         {
            // we assume that the objective function vector has less nonzeros than the reduced cost vector, and so multiplying
            // with -1 first and subtracting the dual activity should be faster than adding the dual activity and negating
            // afterwards
            _rationalLP->getObj(sol._redCost);
            _rationalLP->subDualActivity(sol._dual, sol._redCost);
            _rationalLP->addPrimalQPlin(sol._primal, sol._redCost);
         }
         const int numCorrectedDuals = _primalDualDiff.size();

         if( numCorrectedPrimals + numCorrectedDuals > 0 )
         {
            MSG_INFO2( spxout, spxout << "Corrected " << numCorrectedPrimals << " primal variables and " << numCorrectedDuals << " dual values.\n" );
         }
      }
      while( true );

      // destruct solver instance
#ifdef SOPLEX_WITH_QORE
      if( intParam(SoPlex::QP_SOLVER) == SOLVER_QORE )
         QPFree( &pproblem );
#endif

      // correct basis status for restricted inequalities
      if( _hasBasis )
      {
         for( int r = numRowsRational() - 1; r >= 0; r-- )
         {
            assert((lhsRational(r) == rhsRational(r)) == (_rowTypes[r] == RANGETYPE_FIXED));
            if( _rowTypes[r] != RANGETYPE_FIXED && _basisStatusRows[r] == SPxSolver::FIXED )
               _basisStatusRows[r] = (maximizing == (sol._dual[r] < 0))
                  ? SPxSolver::ON_LOWER
                  : SPxSolver::ON_UPPER;
         }
      }

      // compute objective function values
      assert(sol._hasPrimal == sol._hasDual);
      if( sol._hasPrimal )
      {
	 DVectorRational linearized_objective(_rationalLP->maxObj());
         if( intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE )
            linearized_objective *= -1;
         linearized_objective *= 2;
	 _rationalLP->addPrimalQPlin(sol._primal, linearized_objective);
         linearized_objective *= 0.5;
         sol._primalObjVal = sol._primal * linearized_objective + _rationalLP->objOffset();
         sol._dualObjVal = sol._primalObjVal;
      }

      // set objective coefficients for all rows to zero
      _solver.clearRowObjs();

      // stop rational solving time
      _statistics->rationalTime->stop();
   }



   /// correct dual solution
   int SoPlex::_corrDualSol(SolRational& sol,
      Rational& dualScale,
      const int maxDimRational,
      DVectorReal& dualReal)
   {
      Rational dualScaleInverse = dualScale;
      dualScaleInverse.invert();
      _primalDualDiff.clear();
      int dualSize = 0;
      for( int r = numRowsRational() - 1; r >= 0; r-- )
      {
         SPxSolver::VarStatus& basisStatusRow = _basisStatusRows[r];

         // it may happen that left-hand and right-hand side are different in the rational, but equal in the real LP,
         // leading to a fixed basis status; this is critical because rows with fixed basis status are ignored in the
         // computation of the dual violation; to avoid rational comparisons we do not check this but simply switch
         // to the left-hand side status
         if( basisStatusRow == SPxSolver::FIXED )
            basisStatusRow = SPxSolver::ON_LOWER;

         {
            if( dualReal[r] != 0 )
            {
               sol._dual[r] += dualReal[r]*dualScaleInverse;

               dualSize++;
            }
            else
            {
               // we do not check whether the dual value is nonzero, because it probably is; this gives us an
               // overestimation of the number of nonzeros in the dual solution
               dualSize++;
            }
         }
      }
      return dualSize;
   }



   /// debug: calculate reduced cost and dual violation
   void SoPlex::_dbgRedAndDualViolation(SolRational& sol,
      DVectorReal& dualReal,
      DVectorReal& primalReal,
      const bool maximizing,
      Rational& dualScale)
   {
      {
         // compute reduced cost violation; _sol._dual is not yet corrected by primal solution but _sol._primal is already corrected!
         DVectorRational debugRedCost(numColsRational());
         Rational dualScaleinvers=dualScale;
         dualScaleinvers.invert();
         _rationalLP->getObj(debugRedCost);
         _rationalLP->subDualActivity(sol._dual + (dualScaleinvers * DVectorRational(dualReal)), debugRedCost);
         _rationalLP->addPrimalQPlin(sol._primal, debugRedCost);

         Rational debugRedCostViolation = 0;
         for( int c = numColsRational() - 1; c >= 0; c-- )
         {
            if( _colTypes[c] == RANGETYPE_FIXED )
               continue;

            const SPxSolver::VarStatus& basisStatusCol = _basisStatusCols[c];
            assert(basisStatusCol != SPxSolver::FIXED);

            if( ((maximizing && basisStatusCol != SPxSolver::ON_LOWER) || (!maximizing && basisStatusCol != SPxSolver::ON_UPPER))
               && debugRedCost[c] < -debugRedCostViolation )
            {
               MSG_DEBUG( std::cout << "basisStatusCol = " << basisStatusCol
                  << ", lower tight = " << bool(sol._primal[c] <= lowerRational(c))
                  << ", upper tight = " << bool(sol._primal[c] >= upperRational(c))
                  << ", obj[c] = " << _realLP->obj(c)
                  << ", debugRedCost[c] = " << rationalToString(debugRedCost[c])
                  << "\n" );
               debugRedCostViolation = -debugRedCost[c];
            }

            if( ((maximizing && basisStatusCol != SPxSolver::ON_UPPER) || (!maximizing && basisStatusCol != SPxSolver::ON_LOWER))
               && debugRedCost[c] > debugRedCostViolation )
            {
               MSG_DEBUG( std::cout << "basisStatusCol = " << basisStatusCol
                  << ", lower tight = " << bool(sol._primal[c] <= lowerRational(c))
                  << ", upper tight = " << bool(sol._primal[c] >= upperRational(c))
                  << ", obj[c] = " << _realLP->obj(c)
                  << ", debugRedCost[c] = " << rationalToString(debugRedCost[c])
                  << "\n" );
               debugRedCostViolation = debugRedCost[c];
            }
         }

         // compute dual violation
         Rational debugDualViolation = 0;
         Rational debugBasicDualViolation = 0;
         for( int r = numRowsRational() - 1; r >= 0; r-- )
         {
            if( _rowTypes[r] == RANGETYPE_FIXED )
               continue;

            const SPxSolver::VarStatus& basisStatusRow = _basisStatusRows[r];
            assert(basisStatusRow != SPxSolver::FIXED);

            Rational val = (-dualScale * sol._dual[r]) - Rational(dualReal[r]);

            if( ((maximizing && basisStatusRow != SPxSolver::ON_LOWER) || (!maximizing && basisStatusRow != SPxSolver::ON_UPPER))
               && val > debugDualViolation )
            {
               MSG_DEBUG( std::cout << "basisStatusRow = " << basisStatusRow
                  << ", lower tight = " << bool(sol._slacks[r] <= lhsRational(r))
                  << ", upper tight = " << bool(sol._slacks[r] >= rhsRational(r))
                  << ", dualReal[r] = " << rationalToString(val)
                  << ", dualReal[r] = " << dualReal[r]
                  << "\n" );
               debugDualViolation = val;
            }

            if( ((maximizing && basisStatusRow != SPxSolver::ON_UPPER) || (!maximizing && basisStatusRow != SPxSolver::ON_LOWER))
               && val < -debugDualViolation )
            {
               MSG_DEBUG( std::cout << "basisStatusRow = " << basisStatusRow
                  << ", lower tight = " << bool(sol._slacks[r] <= lhsRational(r))
                  << ", upper tight = " << bool(sol._slacks[r] >= rhsRational(r))
                  << ", dualReal[r] = " << rationalToString(val)
                  << ", dualReal[r] = " << dualReal[r]
                  << "\n" );
               debugDualViolation = -val;
            }

            if( basisStatusRow == SPxSolver::BASIC && spxAbs(val) > debugBasicDualViolation )
            {
               MSG_DEBUG( std::cout << "basisStatusRow = " << basisStatusRow
                  << ", lower tight = " << bool(sol._slacks[r] <= lhsRational(r))
                  << ", upper tight = " << bool(sol._slacks[r] >= rhsRational(r))
                  << ", dualReal[r] = " << rationalToString(val)
                  << ", dualReal[r] = " << dualReal[r]
                  << "\n" );
               debugBasicDualViolation = spxAbs(val);
            }
         }

         if( debugRedCostViolation > _solver.opttol() || debugDualViolation > _solver.opttol() || debugBasicDualViolation > 1e-9 )
         {
            MSG_WARNING( spxout, spxout << "Warning: floating-point dual solution with violation "
               << rationalToString(debugRedCostViolation) << " / "
               << rationalToString(debugDualViolation) << " / "
               << rationalToString(debugBasicDualViolation)
               << " (red. cost, dual, basic).\n" );
         }
      }
   }



   /// correct primal solution
   int SoPlex::_corrPrimSol(SolRational& sol,
      Rational& primalScale,
      const int maxDimRational,
      DVectorReal& primalReal)
   {
      int primalSize = 0;
      Rational primalScaleInverse = primalScale;
      primalScaleInverse.invert();
      _primalDualDiff.clear();
      for( int c = numColsRational() - 1; c >= 0; c-- )
      {
         // force values of nonbasic variables to bounds
         SPxSolver::VarStatus& basisStatusCol = _basisStatusCols[c];

         if( basisStatusCol == SPxSolver::ON_LOWER )
         {
            if( sol._primal[c] != lowerRational(c) )
            {
               sol._primal[c] = lowerRational(c);
            }
         }
         else if( basisStatusCol == SPxSolver::ON_UPPER )
         {
            if( sol._primal[c] != upperRational(c) )
            {
               sol._primal[c] = upperRational(c);
            }
         }
         else if( basisStatusCol == SPxSolver::FIXED )
         {
            // it may happen that lower and upper are only equal in the real LP but different in the rational LP; we
            // do not check this to avoid rational comparisons, but simply switch the basis status to the lower
            // bound; this is necessary, because for fixed variables any reduced cost is feasible
            basisStatusCol = SPxSolver::ON_LOWER;
            if( sol._primal[c] != lowerRational(c) )
            {
               sol._primal[c] = lowerRational(c);
            }
         }
         else if( basisStatusCol == SPxSolver::ZERO )
         {
            if( sol._primal[c] != 0 )
            {
               sol._primal[c] = 0;
            }
         }
         else
         {
            if( primalReal[c] == 1.0 )
            {
               sol._primal[c] += primalScaleInverse;
            }
            else if( primalReal[c] == -1.0 )
            {
               sol._primal[c] -= primalScaleInverse;
            }
            else if( primalReal[c] != 0.0 )
            {
               sol._primal[c] +=  primalReal[c]*primalScaleInverse;
            }
         }

         if( sol._primal[c] != 0 )
            primalSize++;
      }
      return primalSize;
   }



   /// applies dual scaling factor to objective function of QP
   void SoPlex::_appDualScale( SolRational& sol,
      Rational& dualScale )
   {
      Rational temporary_rational;

      if( dualScale < 1 )
         dualScale = 1;
      else
         MSG_INFO2( spxout, spxout << "Scaling dual by " << rationalToString(dualScale) << ".\n" );

      if( intParam(SoPlex::QP_SOLVER) == SOLVER_QPOASES )
      {
#ifdef SOPLEX_WITH_QPOASES
	 // apply scaled objective function
	 for( int c = numColsRational() - 1; c >= 0; c-- )
	 {
	    // perform dual scaling
	    temporary_rational = _modObj[c] * dualScale;
	    if( temporary_rational >= qpOASES::INFTY )
	       g[c] = qpOASES::INFTY;
	    else if( temporary_rational <= -qpOASES::INFTY )
	       g[c] = -qpOASES::INFTY;
	    else
	       g[c] = Real(temporary_rational);
	 }
#endif
      }
      else if( intParam(SoPlex::QP_SOLVER) == SOLVER_QORE )
      {
#ifdef SOPLEX_WITH_QORE
	 // apply scaled objective function
	 for( int c = numColsRational() - 1; c >= 0; c-- )
	 {
	    // perform dual scaling
	    temporary_rational = _modObj[c] * dualScale;
	    if( Real(temporary_rational) >= INFINITY )
	       g[c] = INFINITY;
	    else if( Real(temporary_rational) <= -INFINITY )
	       g[c] = -INFINITY;
	    else
	       g[c] = Real(temporary_rational);
	 }
#endif
      }
   }



   /// applies primal scaling factor to bounds and sides of QP
  void SoPlex::_appPrimScale( Rational& primalScale )
   {
      Rational temporary_rational;

      if( intParam(SoPlex::QP_SOLVER) == SOLVER_QPOASES )
      {
#ifdef SOPLEX_WITH_QPOASES
         // apply scaled bounds
	 if( primalScale <= 1 )
	 {
	    if( primalScale < 1 )
	       primalScale = 1;
	    for( int c = numColsRational() - 1; c >= 0; c-- )
	    {
	       if( _lowerFinite(_colTypes[c]) )
	       {
		  if( _modLower[c] <= -qpOASES::INFTY )
		     lb[c] = -qpOASES::INFTY;
		  else
		     lb[c] = Real(_modLower[c]);
	       }
	       if( _upperFinite(_colTypes[c]) )
	       {
		  if( _modUpper[c] >= qpOASES::INFTY )
		     ub[c] = qpOASES::INFTY;
		  else
		     ub[c] = Real(_modUpper[c]);
	       }
	    }
	 }
	 else
	 {
	    MSG_INFO2( spxout, spxout << "Scaling primal by " << rationalToString(primalScale) << ".\n" );

	    for( int c = numColsRational() - 1; c >= 0; c-- )
	    {
	       if( _lowerFinite(_colTypes[c]) )
	       {
		  temporary_rational = _modLower[c] * primalScale;
		  if( temporary_rational <= -qpOASES::INFTY )
		     lb[c] = -qpOASES::INFTY;
		  else
		     lb[c] = Real(temporary_rational);
	       }
	       if( _upperFinite(_colTypes[c]) )
	       {
		  temporary_rational = _modUpper[c] * primalScale;
		  if( temporary_rational >= qpOASES::INFTY )
		     ub[c] = qpOASES::INFTY;
		  else
		     ub[c] = Real(temporary_rational);
	       }
	    }
	 }

	 // apply scaled sides (constraints)
	 assert(primalScale >= 1);
	 if( primalScale == 1 )
	 {
	    for( int r = numRowsRational() - 1; r >= 0; r-- )
	    {
	       if( _lowerFinite(_rowTypes[r]) )
	       {
		  if( _modLhs[r] <= -qpOASES::INFTY )
		     lbA[r] = -qpOASES::INFTY;
		  else
		     lbA[r] = Real(_modLhs[r]);
	       }
	       if( _upperFinite(_rowTypes[r]) )
	       {
		  if( _modRhs[r] >= qpOASES::INFTY )
		     ubA[r] = qpOASES::INFTY;
		  else
		     ubA[r] = Real(_modRhs[r]);
	       }
	    }
	 }
	 else
	 {
	    for( int r = numRowsRational() - 1; r >= 0; r-- )
	    {
	       if( _lowerFinite(_rowTypes[r]) )
	       {
		  temporary_rational = _modLhs[r] * primalScale;
		  if( temporary_rational <= -qpOASES::INFTY )
		     lbA[r] = -qpOASES::INFTY;
		  else
		     lbA[r] = Real(temporary_rational);
	       }
	       if( _upperFinite(_rowTypes[r]) )
	       {
		  temporary_rational = _modRhs[r] * primalScale;
		  if( temporary_rational >= qpOASES::INFTY )
		     ubA[r] = qpOASES::INFTY;
		  else
		     ubA[r] = Real(temporary_rational);
	       }
	    }
	 }
#endif
      }
      else if( intParam(SoPlex::QP_SOLVER) == SOLVER_QORE )
      {
#ifdef SOPLEX_WITH_QORE
	 // apply scaled bounds
         if( primalScale <= 1 )
	 {
	    if( primalScale < 1 )
	       primalScale = 1;
	    for( int c = numColsRational() - 1; c >= 0; c-- )
	    {
	       if( _lowerFinite(_colTypes[c]) )
	       {
		  if( Real(_modLower[c]) <= -INFINITY )
		     lbt[c] = -INFINITY;
		  else
		     lbt[c] = Real(_modLower[c]);
	       }
	       if( _upperFinite(_colTypes[c]) )
	       {
		  if( Real(_modUpper[c]) >= INFINITY )
		     ubt[c] = INFINITY;
		  else
		     ubt[c] = Real(_modUpper[c]);
	       }
	    }
	 }
	 else
	 {
	    MSG_INFO2( spxout, spxout << "Scaling primal by " << rationalToString(primalScale) << ".\n" );

	    for( int c = numColsRational() - 1; c >= 0; c-- )
	    {
	       if( _lowerFinite(_colTypes[c]) )
	       {
		  temporary_rational = _modLower[c] * primalScale;
		  if( Real(temporary_rational) <= -INFINITY )
		     lbt[c] = -INFINITY;
		  else
		     lbt[c] = Real(temporary_rational);
	       }
	       if( _upperFinite(_colTypes[c]) )
	       {
		  temporary_rational = _modUpper[c] * primalScale;
		  if( Real(temporary_rational) >= INFINITY )
		     ubt[c] = INFINITY;
		  else
		     ubt[c] = Real(temporary_rational);
	       }
	    }
	 }

	 // apply scaled sides (constraints)
         int col_num = numColsRational();
	 assert(primalScale >= 1);
	 if( primalScale == 1 )
	 {
	    for( int r = numRowsRational() - 1; r >= 0; r-- )
	    {
	       if( _lowerFinite(_rowTypes[r]) )
	       {
		  if( Real(_modLhs[r]) <= -INFINITY )
		     lbt[r+col_num] = -INFINITY;
		  else
		     lbt[r+col_num] = Real(_modLhs[r]);
	       }
	       if( _upperFinite(_rowTypes[r]) )
	       {
		  if( Real(_modRhs[r]) >= INFINITY )
		     ubt[r+col_num] = INFINITY;
		  else
		     ubt[r+col_num] = Real(_modRhs[r]);
	       }
	    }
	 }
	 else
	 {
	    for( int r = numRowsRational() - 1; r >= 0; r-- )
	    {
	       if( _lowerFinite(_rowTypes[r]) )
	       {
		  temporary_rational = _modLhs[r] * primalScale;
		  if( Real(temporary_rational) <= -INFINITY )
		     lbt[r+col_num] = -INFINITY;
		  else
		     lbt[r+col_num] = Real(temporary_rational);
	       }
	       if( _upperFinite(_rowTypes[r]) )
	       {
		  temporary_rational = _modRhs[r] * primalScale;
		  if( Real(temporary_rational) >= INFINITY )
		     ubt[r+col_num] = INFINITY;
		  else
		     ubt[r+col_num] = Real(temporary_rational);
	       }
	    }
	 }
#endif
      }
   }



   /// computes dual violation of current solution to the refined QP
   void SoPlex::_dualViolation( SolRational& sol,
      Rational& redCostViolation,
      Rational& dualViolation,
      const bool maximizing )
   {
      // compute reduced cost violation
      redCostViolation = 0;
      for( int c = numColsRational() - 1; c >= 0; c-- )
      {
         if( _colTypes[c] == RANGETYPE_FIXED )
            continue;

         const SPxSolver::VarStatus& basisStatusCol = _basisStatusCols[c];
         assert(basisStatusCol != SPxSolver::FIXED);

         if( ((maximizing && basisStatusCol != SPxSolver::ON_LOWER) || (!maximizing && basisStatusCol != SPxSolver::ON_UPPER))
            && sol._redCost[c] < -redCostViolation )
         {
            MSG_DEBUG( std::cout << "basisStatusCol = " << basisStatusCol
               << ", lower tight = " << bool(sol._primal[c] <= lowerRational(c))
               << ", upper tight = " << bool(sol._primal[c] >= upperRational(c))
               << ", sol._redCost[c] = " << rationalToString(sol._redCost[c])
               << ", _solReal._redCost[c] = " << _solReal._redCost[c]
               << "\n" );
            redCostViolation = -sol._redCost[c];
         }

         if( ((maximizing && basisStatusCol != SPxSolver::ON_UPPER) || (!maximizing && basisStatusCol != SPxSolver::ON_LOWER))
            && sol._redCost[c] > redCostViolation )
         {
            MSG_DEBUG( std::cout << "basisStatusCol = " << basisStatusCol
               << ", lower tight = " << bool(sol._primal[c] <= lowerRational(c))
               << ", upper tight = " << bool(sol._primal[c] >= upperRational(c))
               << ", sol._redCost[c] = " << rationalToString(sol._redCost[c])
               << ", _solReal._redCost[c] = " << _solReal._redCost[c]
               << "\n" );
            redCostViolation = sol._redCost[c];
         }
      }

      // compute dual violation
      dualViolation = 0;
      for( int r = numRowsRational() - 1; r >= 0; r-- )
      {
         if( _rowTypes[r] == RANGETYPE_FIXED )
            continue;

         const SPxSolver::VarStatus& basisStatusRow = _basisStatusRows[r];
         assert(basisStatusRow != SPxSolver::FIXED);

         if( ((maximizing && basisStatusRow != SPxSolver::ON_LOWER) || (!maximizing && basisStatusRow != SPxSolver::ON_UPPER))
            && sol._dual[r] < -dualViolation )
         {
            MSG_DEBUG( std::cout << "basisStatusRow = " << basisStatusRow
               << ", lower tight = " << bool(sol._slacks[r] <= lhsRational(r))
               << ", upper tight = " << bool(sol._slacks[r] >= rhsRational(r))
               << ", sol._dual[r] = " << rationalToString(sol._dual[r])
               << "\n" );
            dualViolation = -sol._dual[r];
         }

         if( ((maximizing && basisStatusRow != SPxSolver::ON_UPPER) || (!maximizing && basisStatusRow != SPxSolver::ON_LOWER))
            && sol._dual[r] > dualViolation )
         {
            MSG_DEBUG( std::cout << "basisStatusRow = " << basisStatusRow
               << ", lower tight = " << bool(sol._slacks[r] <= lhsRational(r))
               << ", upper tight = " << bool(sol._slacks[r] >= rhsRational(r))
               << ", sol._dual[r] = " << rationalToString(sol._dual[r])
               << "\n" );
            dualViolation = sol._dual[r];
         }
      }

      _modObj = sol._redCost;
   }



   /// computes primal violation of current soulution to the refined QP
   void SoPlex::_primViolation( SolRational& sol,
      Rational& boundsViolation,
      Rational& sideViolation )
   {
      // compute violation of bounds
      boundsViolation = 0;
      for( int c = numColsRational() - 1; c >= 0; c-- )
      {
         // lower bound
         assert((lowerRational(c) > _rationalNegInfty) == _lowerFinite(_colTypes[c]));
         if( _lowerFinite(_colTypes[c]) )
         {
            if( lowerRational(c) == 0 )
            {
               _modLower[c] = sol._primal[c];
               _modLower[c] *= -1;
               if( _modLower[c] > boundsViolation )
                  boundsViolation = _modLower[c];
            }
            else
            {
               _modLower[c] = lowerRational(c);
               _modLower[c] -= sol._primal[c];
               if( _modLower[c] > boundsViolation )
                  boundsViolation = _modLower[c];
            }
         }

         // upper bound
         assert((upperRational(c) < _rationalPosInfty) == _upperFinite(_colTypes[c]));
         if( _upperFinite(_colTypes[c]) )
         {
            if( upperRational(c) == 0 )
            {
               _modUpper[c] = sol._primal[c];
               _modUpper[c] *= -1;
               if( _modUpper[c] < -boundsViolation )
                  boundsViolation = -_modUpper[c];
            }
            else
            {
               _modUpper[c] = upperRational(c);
               _modUpper[c] -= sol._primal[c];
               if( _modUpper[c] < -boundsViolation )
                  boundsViolation = -_modUpper[c];
            }
         }
      }

      // compute violation of sides
      sideViolation = 0;
      for( int r = numRowsRational() - 1; r >= 0; r-- )
      {
         const SPxSolver::VarStatus& basisStatusRow = _basisStatusRows[r];

         // left-hand side
         assert((lhsRational(r) > _rationalNegInfty) == _lowerFinite(_rowTypes[r]));
         if( _lowerFinite(_rowTypes[r]) )
         {
            if( lhsRational(r) == 0 )
            {
               _modLhs[r] = sol._slacks[r];
               _modLhs[r] *= -1;
            }
            else
            {
               _modLhs[r] = lhsRational(r);
               _modLhs[r] -= sol._slacks[r];
            }

            if( _modLhs[r] > sideViolation )
               sideViolation = _modLhs[r];
            // if the activity is feasible, but too far from the bound, this violates complementary slackness; we
            // count it as side violation here
            else if( basisStatusRow == SPxSolver::ON_LOWER && _modLhs[r] < -sideViolation )
               sideViolation = -_modLhs[r];
         }

         // right-hand side
         assert((rhsRational(r) < _rationalPosInfty) == _upperFinite(_rowTypes[r]));
         if( _upperFinite(_rowTypes[r]) )
         {
            if( rhsRational(r) == 0 )
            {
               _modRhs[r] = sol._slacks[r];
               _modRhs[r] *= -1;
            }
            else
            {
               _modRhs[r] = rhsRational(r);
               _modRhs[r] -= sol._slacks[r];
            }

            if( _modRhs[r] < -sideViolation )
               sideViolation = -_modRhs[r];
            // if the activity is feasible, but too far from the bound, this violates complementary slackness; we
            // count it as side violation here
            else if( basisStatusRow == SPxSolver::ON_UPPER && _modRhs[r] > sideViolation )
               sideViolation = _modRhs[r];
         }
      }
   }



#ifdef SOPLEX_WITH_QPOASES
   /// call qpOASES to solve real QP
   SPxSolver::Status SoPlex::_solveRealQP(QProblem& soplex_ir_qp,
      const bool acceptUnbounded,
      const bool acceptInfeasible,
      VectorReal& primal,
      VectorReal& dual,
      DataArray< SPxSolver::VarStatus >& basisStatusRows,
      DataArray< SPxSolver::VarStatus >& basisStatusCols,
      bool& returnedBasis,
      const bool forceNoSimplifier)
   {
      // qp dimensions
      int var_num = _realLP->nCols();
      int constr_num = _realLP->nRows();

      // flag if we want reliable options
      bool reliable;

      // make copy of qp problem to recover failed solving
      QProblem soplex_ir_qp_copy;
      if( SoPlex::boolParam(SoPlex::RELIABLERESOLVE) )
      {
         soplex_ir_qp_copy = QProblem(var_num, constr_num);
         soplex_ir_qp_copy = soplex_ir_qp;
      }

      // fill qpOASES arrays with QP from soplex
      _copyQPintoqpOASES(var_num, constr_num);

      // set options for qpOASES
      Options options;
      // default
      reliable = false;
      _setQPoasesOptions( options, reliable);
      if( SoPlex::boolParam(SoPlex::RELIABLERESOLVE) )
         soplex_ir_qp_copy.setOptions( options );
      else
	 soplex_ir_qp.setOptions( options );

      // prepare first QP solve
      MSG_INFO1( spxout, spxout << "Calling qpOASES initially . . .\n" );
      int nWSR;
      if( intParam(SoPlex::ITERLIMIT) >= 0 )
         nWSR = intParam(SoPlex::ITERLIMIT);
      else
	 nWSR = 10000;
      real_t* cputime = 0;
      real_t cputime_var;
      if( realParam(SoPlex::TIMELIMIT) < realParam(SoPlex::INFTY) )
      {
	 cputime = &cputime_var;
         *cputime = realParam(SoPlex::TIMELIMIT) - _statistics->solvingTime->time();
         if( *cputime <= 0 )
         {
            MSG_INFO2( spxout, spxout << "Timelimit reached, abort calling qpOASES.\n" );
            _status = SPxSolver::ABORT_TIME;
            _hasSolReal = false;
	    return SPxSolver::ABORT_TIME;
         }
      }

      // solve first QP
      returnValue QP_STATUS;
      if( SoPlex::boolParam(SoPlex::RELIABLERESOLVE) )
      {
         if( _useSparseQP == 0 )
            QP_STATUS = soplex_ir_qp_copy.init(H, g, A, lb, ub, lbA, ubA, nWSR, cputime);
#ifdef SOPLEX_WITH_QPOASES_SPARSE
         else
	    QP_STATUS = soplex_ir_qp_copy.init(H_spar, g, A_spar, lb, ub, lbA, ubA, nWSR, cputime);
#endif
         assert( nWSR >= 0 );
         if( QP_STATUS == SUCCESSFUL_RETURN || nWSR < 10000 )
            _statistics->iterations += nWSR;
      }
      else
      {
         if( _useSparseQP == 0 )
            QP_STATUS = soplex_ir_qp.init(H, g, A, lb, ub, lbA, ubA, nWSR, cputime);
#ifdef SOPLEX_WITH_QPOASES_SPARSE
         else
	    QP_STATUS = soplex_ir_qp.init(H_spar, g, A_spar, lb, ub, lbA, ubA, nWSR, cputime);
#endif
         assert( nWSR >= 0 );
         if( QP_STATUS == SUCCESSFUL_RETURN || nWSR < 10000 )
            _statistics->iterations += nWSR;
      }

      // try again with reliable settings
      if( QP_STATUS != SUCCESSFUL_RETURN && SoPlex::boolParam(SoPlex::RELIABLERESOLVE) )
      {
	 // count resolves
	 _statistics->resolve += 1;

	 // get original QP
         soplex_ir_qp_copy = soplex_ir_qp;

         // set reliable options for qpOASES
         reliable = true;
         _setQPoasesOptions( options, reliable);
         soplex_ir_qp_copy.setOptions( options );

         // prepare first QP solve again
         MSG_INFO1( spxout, spxout << "Unsuccessful, calling qpOASES initially with reliable settings . . .\n" );
         if( intParam(SoPlex::ITERLIMIT) >= 0 )
            nWSR = intParam(SoPlex::ITERLIMIT);
         else
            nWSR = 10000;
         cputime = 0;
         if( realParam(SoPlex::TIMELIMIT) < realParam(SoPlex::INFTY) )
         {
            cputime = &cputime_var;
            *cputime = realParam(SoPlex::TIMELIMIT) - _statistics->solvingTime->time();
            if( *cputime <= 0 )
            {
               MSG_INFO2( spxout, spxout << "Timelimit reached, abort calling qpOASES.\n" );
               _status = SPxSolver::ABORT_TIME;
               _hasSolReal = false;
               return SPxSolver::ABORT_TIME;
            }
         }

         // solve first QP again
         if( _useSparseQP == 0 )
            QP_STATUS = soplex_ir_qp_copy.init(H, g, A, lb, ub, lbA, ubA, nWSR, cputime);
#ifdef SOPLEX_WITH_QPOASES_SPARSE
         else
            QP_STATUS = soplex_ir_qp_copy.init(H_spar, g, A_spar, lb, ub, lbA, ubA, nWSR, cputime);
#endif
         assert( nWSR >= 0 );
         if( QP_STATUS == SUCCESSFUL_RETURN || nWSR < 10000 )
            _statistics->iterations += nWSR;
      }

      // copy solution back to original qp problem
      if( SoPlex::boolParam(SoPlex::RELIABLERESOLVE) )
         soplex_ir_qp = soplex_ir_qp_copy;

      // get qpOASES tolerances with SolutionAnalysis class
      real_t maxStat;
      real_t maxFeas;
      real_t maxCmpl;
      SolutionAnalysis solutionAnalysis;
      real_t maxKktViolation = solutionAnalysis.getKktViolation( &soplex_ir_qp, &maxStat, &maxFeas, &maxCmpl);
      MSG_INFO1( spxout, spxout << "Maximal violation of KKT system (qpOASES) : " << double(maxKktViolation) << "\n"  );

      if( QP_STATUS == SUCCESSFUL_RETURN )
      {
         _status = SPxSolver::OPTIMAL;
         _hasSolReal = true;
	 returnedBasis = true;
	 _copySolintoSPX(soplex_ir_qp, primal, dual, var_num, constr_num);
	 return SPxSolver::OPTIMAL;
      }
      else if( QP_STATUS == RET_MAX_NWSR_REACHED )
      {
         _status = SPxSolver::ABORT_ITER;
         _hasSolReal = false;
	 return SPxSolver::ABORT_ITER;
      }
      else // QP_STATUS == RET_INIT_FAILED
      {
         _status = SPxSolver::ERROR;
         _hasSolReal = false;
	 return SPxSolver::ERROR;
      }
   }
#endif



#ifdef SOPLEX_WITH_QORE
   /// call QORE to solve real QP
   SPxSolver::Status SoPlex::_solveRealQP(QoreProblem* pproblem,
      const bool acceptUnbounded,
      const bool acceptInfeasible,
      VectorReal& primal,
      VectorReal& dual,
      DataArray< SPxSolver::VarStatus >& basisStatusRows,
      DataArray< SPxSolver::VarStatus >& basisStatusCols,
      bool& returnedBasis,
      const bool forceNoSimplifier)
   {
      // qp dimensions
      int var_num = _realLP->nCols();
      int constr_num = _realLP->nRows();

      // fill QORE arrays with QP from soplex
      _copyQPintoQORE(var_num, constr_num);

      // set options for QORE
      QPSetDefault(pproblem);

      // prepare first QP solve
      MSG_INFO1( spxout, spxout << "Calling QORE initially . . .\n" );
      double* cputime = 0;
      double cputime_var;
      if( realParam(SoPlex::TIMELIMIT) < realParam(SoPlex::INFTY) )
      {
	 cputime = &cputime_var;
         *cputime = realParam(SoPlex::TIMELIMIT) - _statistics->solvingTime->time();
         if( *cputime <= 0 )
         {
            MSG_INFO2( spxout, spxout << "Timelimit reached, abort calling QORE.\n" );
            _status = SPxSolver::ABORT_TIME;
            _hasSolReal = false;
	    return SPxSolver::ABORT_TIME;
         }
      }
      qp_int rv = QPSOLVER_OK;
      rv = QPSetData( pproblem, var_num, constr_num, Atcp, Atri, Atnz, Hcp, Hri, Hnz );
      assert( rv == QPSOLVER_OK );

      // free matrix arrays, not used anymore (QPSetData() made deep copies)
      spx_free(Hcp);
      spx_free(Hri);
      spx_free(Hnz);
      spx_free(Atcp);
      spx_free(Atri);
      spx_free(Atnz);

      // solve first QP
      rv = QPOptimize( pproblem, lbt, ubt, g, 0, 0 );

      if( rv == QPSOLVER_OK )
      {
         _status = SPxSolver::OPTIMAL;
         _hasSolReal = true;
	 returnedBasis = true;
	 _copySolintoSPX(pproblem, primal, dual, var_num, constr_num);
	 return SPxSolver::OPTIMAL;
      }
      else
      {
         _status = SPxSolver::ERROR;
         _hasSolReal = false;
	 return SPxSolver::ERROR;
      }
   }
#endif



#ifdef SOPLEX_WITH_QPOASES
   /// call qpOASES to resolve real QP
   SPxSolver::Status SoPlex::_resolveRealQP(QProblem& soplex_ir_qp,
      const bool acceptUnbounded,
      const bool acceptInfeasible,
      VectorReal& primal,
      VectorReal& dual,
      DataArray< SPxSolver::VarStatus >& basisStatusRows,
      DataArray< SPxSolver::VarStatus >& basisStatusCols,
      bool& returnedBasis,
      const bool forceNoSimplifier)
   {
      // qp dimensions
      int var_num = numColsRational();
      int constr_num = numRowsRational();

      // flag if we want reliable options
      bool reliable;

      // make copies of qp problem to recover failed solving
      QProblem soplex_ir_qp_copy_one;
      if( SoPlex::boolParam(SoPlex::RELIABLERESOLVE) )
      {
         soplex_ir_qp_copy_one = QProblem(var_num, constr_num);
         soplex_ir_qp_copy_one = soplex_ir_qp;
      }

      // declare Bounds and Constraints object for qpOASES
      Bounds qpBounds;
      Constraints qpConstraints;
      if( SoPlex::boolParam(SoPlex::RELIABLERESOLVE) )
      {
         soplex_ir_qp_copy_one.getBounds( qpBounds );
         soplex_ir_qp_copy_one.getConstraints( qpConstraints );
      }
      else
      {
         soplex_ir_qp.getBounds( qpBounds );
         soplex_ir_qp.getConstraints( qpConstraints );
      }

      // load basis into qpOASES objects
      if( returnedBasis )
      {
         SPxSolver::VarStatus varStatus;
         for( int i = 0; i < var_num; i++ )
            {
            varStatus = _basisStatusCols[i];
            if( varStatus == SPxSolver::ON_LOWER )
            {
               qpBounds.setStatus( i, ST_LOWER );
            }
            else if( varStatus == SPxSolver::BASIC )
            {
               qpBounds.setStatus( i, ST_INACTIVE );
            }
            else if( varStatus == SPxSolver::ON_UPPER )
            {
               qpBounds.setStatus( i, ST_UPPER );
            }
            else
            {
               qpBounds.setStatus( i, ST_UNDEFINED);
            }
         }
         for( int i = 0; i < constr_num; i++ )
         {
            varStatus = _basisStatusRows[i];
            if( varStatus == SPxSolver::ON_LOWER )
            {
               qpConstraints.setStatus( i, ST_LOWER );
            }
            else if( varStatus == SPxSolver::BASIC )
            {
               qpConstraints.setStatus( i, ST_INACTIVE );
            }
            else if( varStatus == SPxSolver::ON_UPPER )
            {
               qpConstraints.setStatus( i, ST_UPPER );
            }
            else
            {
               qpConstraints.setStatus( i, ST_UNDEFINED );
            }
         }
      }

      // set options for qpOASES
      Options options;
      reliable = false;
      _setQPoasesOptions( options, reliable );
      if( SoPlex::boolParam(SoPlex::RELIABLERESOLVE) )
         soplex_ir_qp_copy_one.setOptions( options );
      else
         soplex_ir_qp.setOptions( options );

      // prepare next QP solve
      MSG_INFO1( spxout, spxout << "Calling qpOASES again . . .\n" );
      int nWSR;
      if( intParam(SoPlex::ITERLIMIT) >= _statistics->iterations  && intParam(SoPlex::ITERLIMIT) >= 0 )
         nWSR = intParam(SoPlex::ITERLIMIT) - _statistics->iterations;
      else
	 nWSR = 10000;
      real_t* cputime = 0;
      real_t cputime_var;
      if( realParam(SoPlex::TIMELIMIT) < realParam(SoPlex::INFTY) )
      {
         cputime = &cputime_var;
         *cputime = realParam(SoPlex::TIMELIMIT) - _statistics->solvingTime->time();
         if( *cputime <= 0 )
         {
            MSG_INFO2( spxout, spxout << "Timelimit reached, abort calling qpOASES.\n" );
            _status = SPxSolver::ABORT_TIME;
            _hasSolReal = false;
	    return SPxSolver::ABORT_TIME;
         }
      }

      // solve next QP
      returnValue QP_STATUS;
      if( SoPlex::boolParam(SoPlex::RELIABLERESOLVE) )
         QP_STATUS = soplex_ir_qp_copy_one.hotstart( g, lb, ub, lbA, ubA, nWSR, cputime, &qpBounds, &qpConstraints );
      else
         QP_STATUS = soplex_ir_qp.hotstart( g, lb, ub, lbA, ubA, nWSR, cputime, &qpBounds, &qpConstraints );
      assert( nWSR >= 0 );
      if( QP_STATUS == SUCCESSFUL_RETURN || nWSR < 10000 )
         _statistics->iterations += nWSR;

      // try with reliable settings again if solver failed
      if( QP_STATUS != SUCCESSFUL_RETURN && SoPlex::boolParam(SoPlex::RELIABLERESOLVE) )
      {
	 // count resolves
	 _statistics->resolve += 1;

         // copy QP for repeated next solve
	 soplex_ir_qp_copy_one = soplex_ir_qp;

         // set qpOASES options
         reliable = true;
         _setQPoasesOptions( options, reliable );
         soplex_ir_qp_copy_one.setOptions( options );

         // prepare next qp solve again
         MSG_INFO1( spxout, spxout << "Unsuccessful, calling qpOASES again with reliable settings . . .\n" );
         if( intParam(SoPlex::ITERLIMIT) >= _statistics->iterations  && intParam(SoPlex::ITERLIMIT) >= 0 )
            nWSR = intParam(SoPlex::ITERLIMIT) - _statistics->iterations;
         else
	    nWSR = 10000;
         cputime = 0;
         if( realParam(SoPlex::TIMELIMIT) < realParam(SoPlex::INFTY) )
         {
            cputime = &cputime_var;
            *cputime = realParam(SoPlex::TIMELIMIT) - _statistics->solvingTime->time();
            if( *cputime <= 0 )
            {
               MSG_INFO2( spxout, spxout << "Timelimit reached, abort calling qpOASES.\n" );
               _status = SPxSolver::ABORT_TIME;
               _hasSolReal = false;
               return SPxSolver::ABORT_TIME;
            }
         }

         // perform next QP solve again
         QP_STATUS = soplex_ir_qp_copy_one.hotstart( g, lb, ub, lbA, ubA, nWSR, cputime);//, &qpBounds, &qpConstraints );
         assert( nWSR >= 0 );
         if( QP_STATUS == SUCCESSFUL_RETURN || nWSR < 10000 )
            _statistics->iterations += nWSR;
      }

      // copy solution back to original qp problem
      if( SoPlex::boolParam(SoPlex::RELIABLERESOLVE) )
         soplex_ir_qp = soplex_ir_qp_copy_one;

      // get qpOASES tolerances with SolutionAnalysis class
      real_t maxStat;
      real_t maxFeas;
      real_t maxCmpl;
      SolutionAnalysis solutionAnalysis;
      real_t maxKktViolation = solutionAnalysis.getKktViolation( &soplex_ir_qp, &maxStat, &maxFeas, &maxCmpl);
      MSG_INFO1( spxout, spxout << "Maximal violation of KKT system (qpOASES) : " << double(maxKktViolation) << "\n"  );

      if( QP_STATUS == SUCCESSFUL_RETURN )
      {
         _status = SPxSolver::OPTIMAL;
         _hasSolReal = false;
	 returnedBasis = true;
	 _copySolintoSPX(soplex_ir_qp, primal, dual, var_num, constr_num);
	 return SPxSolver::OPTIMAL;
      }
      else if( QP_STATUS == RET_MAX_NWSR_REACHED )
      {
         _status = SPxSolver::ABORT_ITER;
         _hasSolReal = false;
	 return SPxSolver::ABORT_ITER;
      }
      else // QP_STATUS == RET_INIT_FAILED
      {
         _status = SPxSolver::ERROR;
         _hasSolReal = false;
	 return SPxSolver::ERROR;
      }
   }
#endif



#ifdef SOPLEX_WITH_QORE
   /// call QORE to resolve real QP
   SPxSolver::Status SoPlex::_resolveRealQP(QoreProblem* pproblem,
      const bool acceptUnbounded,
      const bool acceptInfeasible,
      VectorReal& primal,
      VectorReal& dual,
      DataArray< SPxSolver::VarStatus >& basisStatusRows,
      DataArray< SPxSolver::VarStatus >& basisStatusCols,
      bool& returnedBasis,
      const bool forceNoSimplifier)
   {
      // qp dimensions
      int var_num = numColsRational();
      int constr_num = numRowsRational();

      // set options for QORE
      QPSetDefault(pproblem);

      // prepare next QP solve
      MSG_INFO1( spxout, spxout << "Calling QORE again . . .\n" );
      double* cputime = 0;
      double cputime_var;
      if( realParam(SoPlex::TIMELIMIT) < realParam(SoPlex::INFTY) )
      {
         cputime = &cputime_var;
         *cputime = realParam(SoPlex::TIMELIMIT) - _statistics->solvingTime->time();
         if( *cputime <= 0 )
         {
            MSG_INFO2( spxout, spxout << "Timelimit reached, abort calling qpOASES.\n" );
            _status = SPxSolver::ABORT_TIME;
            _hasSolReal = false;
	    return SPxSolver::ABORT_TIME;
         }
      }

      // solve next QP
      qp_int rv = QPSOLVER_OK;
      rv = QPOptimize(pproblem, lbt, ubt, g, 0, 0 );

      if( rv == QPSOLVER_OK )
      {
         _status = SPxSolver::OPTIMAL;
         _hasSolReal = false;
	 returnedBasis = true;
	 _copySolintoSPX(pproblem, primal, dual, var_num, constr_num);
	 return SPxSolver::OPTIMAL;
      }
      else
      {
         _status = SPxSolver::ERROR;
         _hasSolReal = false;
	 return SPxSolver::ERROR;
      }
   }
#endif


   /// stores floating-point solution of original LP as current rational solution and ensures that solution vectors have
   /// right dimension; ensure that solution is aligned with basis
   void SoPlex::_realSolToRational( SolRational& sol,
      DVectorReal& primalReal,
      DVectorReal& dualReal)
   {
      sol._primal.reDim(numColsRational(), false);
      sol._slacks.reDim(numRowsRational(), false);
      sol._dual.reDim(numRowsRational(), false);
      sol._redCost.reDim(numColsRational(), false);
      sol._hasPrimal = true;
      sol._hasDual = true;

      for( int c = numColsRational() - 1; c >= 0; c-- )
      {
         SPxSolver::VarStatus& basisStatusCol = _basisStatusCols[c];

         if( basisStatusCol == SPxSolver::ON_LOWER )
            sol._primal[c] = lowerRational(c);
	 else if( basisStatusCol == SPxSolver::ON_UPPER )
            sol._primal[c] = upperRational(c);
         else if( basisStatusCol == SPxSolver::FIXED )
         {
            // it may happen that lower and upper are only equal in the real LP but different in the rational LP; we do
            // not check this to avoid rational comparisons, but simply switch the basis status to the lower bound; this
            // is necessary, because for fixed variables any reduced cost is feasible
            sol._primal[c] = lowerRational(c);
            basisStatusCol = SPxSolver::ON_LOWER;
         }
         else if( basisStatusCol == SPxSolver::ZERO )
            sol._primal[c] = 0;
         else
            sol._primal[c] = primalReal[c];
      }
      _rationalLP->computePrimalActivity(sol._primal, sol._slacks);

      int dualSize = 0;
      for( int r = numRowsRational() - 1; r >= 0; r-- )
      {
         SPxSolver::VarStatus& basisStatusRow = _basisStatusRows[r];

         // it may happen that left-hand and right-hand side are different in the rational, but equal in the real LP,
         // leading to a fixed basis status; this is critical because rows with fixed basis status are ignored in the
         // computation of the dual violation; to avoid rational comparisons we do not check this but simply switch to
         // the left-hand side status
         if( basisStatusRow == SPxSolver::FIXED )
            basisStatusRow = SPxSolver::ON_LOWER;

         {
            sol._dual[r] = dualReal[r];
            if( dualReal[r] != 0.0 )
               dualSize++;
         }
      }
      // we assume that the objective function vector has less nonzeros than the reduced cost vector, and so multiplying
      // with -1 first and subtracting the dual activity should be faster than adding the dual activity and negating
      // afterwards
      _rationalLP->getObj(sol._redCost);
      _rationalLP->subDualActivity(sol._dual, sol._redCost);
      _rationalLP->addPrimalQPlin(sol._primal, sol._redCost);
   }



#ifdef SOPLEX_WITH_QPOASES
   /// set appropriate options for all qpOASES calls
   void SoPlex::_setQPoasesOptions( Options& options,
      const bool reliable)
   {
      if( reliable == false )
      {
         if( QPOASES_SETTINGS_DEFAULT == intParam(SoPlex::QPOASES_SETTINGS) )
	 {
            options.setToMPC();
            options.enableNZCTests = BT_TRUE;
            options.enableDriftCorrection = BT_TRUE;
            options.enableRamping = BT_TRUE;
            options.terminationTolerance = 1e-3;
            if( boolParam(SoPlex::LITESTFULL) )
               options.enableFullLITests = BT_TRUE;
	 }
	 else if( QPOASES_SETTINGS_FAST == intParam(SoPlex::QPOASES_SETTINGS) )
	 {
            options.setToMPC();
	 }
	 else if( QPOASES_SETTINGS_MEDIUM == intParam(SoPlex::QPOASES_SETTINGS) )
	 {
            options.setToDefault();
	 }
	 else if( QPOASES_SETTINGS_RELIABLE == intParam(SoPlex::QPOASES_SETTINGS) )
	 {
            options.setToReliable();
	 }
      }
      else
      {
         options.setToReliable();
         options.numRefinementSteps = 10;
         if( boolParam(SoPlex::LITESTFULL) )
            options.enableFullLITests = BT_TRUE;
      }
   }
#endif



   /// factorizes rational basis matrix in column representation
   void SoPlex::_factorizeColumnRowRational_qp(SolRational& sol,
      DataArray< SPxSolver::VarStatus >& basisStatusRows,
      DataArray< SPxSolver::VarStatus >& basisStatusCols,
      bool& stoppedTime,
      bool& stoppedIter,
      bool& error,
      bool& optimal)
   {
      // start rational solving time
      _statistics->rationalTime->start();

      stoppedTime = false;
      stoppedIter = false;
      error = false;
      optimal = false;

      assert(basisStatusRows.size() == numRowsRational());
      assert(basisStatusCols.size() == numColsRational());

      const bool maximizing = (intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MAXIMIZE);
      int matrixdim = basisStatusRows.size() + basisStatusCols.size();
      int numBasicRows, numBasicCols;

      _workSol._primal.reDim(basisStatusCols.size());
      _workSol._slacks.reDim(basisStatusRows.size());
      _workSol._dual.reDim(basisStatusRows.size());
      _workSol._redCost.reDim(basisStatusCols.size());

      DVectorRational& Primal = _workSol._primal;
      DVectorRational& Dual = _workSol._dual;

      DVectorRational righthandside(basisStatusCols.size()+basisStatusRows.size());
      DVectorRational solution;

      Rational violation;
      Rational primalViolation;
      Rational dualViolation;
      bool primalFeasible = false;
      bool dualFeasible = false;

      int j;

      _computeAndLoadRationalQPBasis( error, basisStatusRows, basisStatusCols, righthandside, matrixdim, numBasicCols, numBasicRows);
      if( error == true )
         goto TERMINATE;

      if( _rationalLUSolver.status() == SLinSolverRational::TIME )
      {
         MSG_INFO2( spxout, spxout << "Rational factorization hit time limit.\n" );
         stoppedTime = true;
         return;
      }
      else if( _rationalLUSolver.status() != SLinSolverRational::OK )
      {
         MSG_ERROR( std::cerr << "Error performing rational LU factorization.\n" );
         error = true;
         return;
      }

      // solve basis system (KKT: primal/dual)
      solution.reDim(matrixdim);
      if( realParam(SoPlex::TIMELIMIT) < realParam(SoPlex::INFTY) )
         _rationalLUSolver.setTimeLimit(realParam(SoPlex::TIMELIMIT) - _statistics->solvingTime->time());
      else
         _rationalLUSolver.setTimeLimit(-1.0);
      _rationalLUSolver.solveRight(solution, righthandside);

      // record statistics
      _statistics->luSolveTimeRational += _rationalLUSolver.getSolveTime();
      _rationalLUSolver.resetCounters();

      if( _isSolveStopped(stoppedTime, stoppedIter) )
      {
         MSG_INFO2( spxout, spxout << "Rational factorization hit time limit while solving rational basis system.\n" );
         return;
      }

      // reconstruct full primal solution
      j = 0;
      for( int i = 0; i < basisStatusCols.size(); i++ )
      {
         if( basisStatusCols[i] == SPxSolver::BASIC )
         {
            assert(j < matrixdim);
            Primal[i] = solution[j];
            j++;
	 }
         else if( basisStatusCols[i] == SPxSolver::ON_LOWER )
         {
            Primal[i] = lowerRational(i);
         }
         else if( basisStatusCols[i] == SPxSolver::ON_UPPER )
         {
            Primal[i] = upperRational(i);
         }
         else if( basisStatusCols[i] == SPxSolver::ZERO )
         {
	    Primal[i] = 0;
         }
         else if( basisStatusCols[i] == SPxSolver::FIXED )
         {
            assert(lowerRational(i) == upperRational(i));
            Primal[i] = lowerRational(i);
         }
         else if( basisStatusCols[i] == SPxSolver::UNDEFINED )
         {
            MSG_ERROR( std::cerr << "Undefined basis status of column in rational factorization.\n" );
            error = true;
            goto TERMINATE;
         }
      }
      assert( j == numBasicCols );

      for( int i = 0; i < basisStatusRows.size(); i++ )
      {
         assert(j < matrixdim);
         Dual[i] = solution[j];
         j++;
      }
      assert( j == matrixdim );

      // check bound violation on basic columns (primal and slack)
      j = 0;
      primalViolation = 0;
      primalFeasible = true;
      for( int i = 0; i < basisStatusCols.size(); i++ )
      {
         if( basisStatusCols[i] == SPxSolver::BASIC )
         {
            assert(j < matrixdim);
            assert(_rationalLUSolverBind[j] == i);
            if( solution[j] < lowerRational(i) )
            {
               violation = lowerRational(i);
               violation -= solution[j];
               if( violation > primalViolation )
               {
                  primalFeasible = false;
                  primalViolation = violation;
               }
            }
            if( solution[j] > upperRational(i) )
            {
               violation = solution[j];
               violation -= upperRational(i);
               if( violation > primalViolation )
               {
                  primalFeasible = false;
                  primalViolation = violation;
               }
            }
            j++;
         }
      }
      assert( j == numBasicCols );

      if( !primalFeasible )
      {
         MSG_INFO1( spxout, spxout << "Rational solution primal infeasible.\n" );
      }

      // check reduced cost violation on nonbasic columns (primal and slack)
      dualViolation = 0;
      dualFeasible = true;
      for( int i = 0; i < basisStatusCols.size(); i++ )
      {
         if( _colTypes[i] == RANGETYPE_FIXED
            && (basisStatusCols[i] == SPxSolver::ON_LOWER || basisStatusCols[i] == SPxSolver::ON_UPPER) )
         {
            assert(lowerRational(i) == upperRational(i));
            basisStatusCols[i] = SPxSolver::FIXED;
         }

         assert(basisStatusCols[i] != SPxSolver::BASIC || Dual * colVectorRational(i) - _rationalLP->PrimalQPlinentry(Primal, i) == objRational(i));
         if( basisStatusCols[i] == SPxSolver::BASIC || basisStatusCols[i] == SPxSolver::FIXED )
            continue;
         else
         {
            _workSol._redCost[i] = Dual * colVectorRational(i);
            _workSol._redCost[i] -= _rationalLP->PrimalQPlinentry(Primal, i);
            _workSol._redCost[i] -= objRational(i);
            if( _workSol._redCost[i] > 0 )
            {
               if( ((maximizing && basisStatusCols[i] != SPxSolver::ON_LOWER) || (!maximizing && basisStatusCols[i] != SPxSolver::ON_UPPER))
                  && (basisStatusCols[i] != SPxSolver::ZERO || upperRational(i) != 0) )
               {
                  dualFeasible = false;
                  if( _workSol._redCost[i] > dualViolation )
                     dualViolation = _workSol._redCost[i];
               }
               _workSol._redCost[i] *= -1;
            }
            else if( _workSol._redCost[i] < 0 )
            {
               _workSol._redCost[i] *= -1;
               if( ((maximizing && basisStatusCols[i] != SPxSolver::ON_UPPER) || (!maximizing && basisStatusCols[i] != SPxSolver::ON_LOWER))
                  && (basisStatusCols[i] != SPxSolver::ZERO || lowerRational(i) != 0) )
               {
                  dualFeasible = false;
                  if( _workSol._redCost[i] > dualViolation )
                     dualViolation = _workSol._redCost[i];
               }
            }
            else
               _workSol._redCost[i] *= -1;
         }
      }

      if( !dualFeasible )
      {
         MSG_INFO1( spxout, spxout << "Rational solution dual infeasible.\n" );
      }

      // store solution
      optimal = primalFeasible && dualFeasible;
      if( optimal || boolParam(SoPlex::RATFACJUMP) )
      {
         _hasBasis = true;
         if( &basisStatusRows != &_basisStatusRows )
            _basisStatusRows = basisStatusRows;
         if( &basisStatusCols != &_basisStatusCols )
            _basisStatusCols = basisStatusCols;

         sol._primal.reDim(numColsRational());
         for( int i = 0; i < basisStatusCols.size(); i++ )
         {
            sol._primal[i] = Primal[i];
         }
         sol._slacks.reDim(numRowsRational());
         _rationalLP->computePrimalActivity(sol._primal, sol._slacks);
         sol._hasPrimal = true;

         sol._dual = Dual;
         for( int i = 0; i < numColsRational(); i++ )
         {
            if( basisStatusCols[i] == SPxSolver::BASIC )
               sol._redCost[i] = 0;
            else if( basisStatusCols[i] == SPxSolver::FIXED )
            {
               sol._redCost[i] = Dual * colVectorRational(i);
               sol._redCost[i] -= _rationalLP->PrimalQPlinentry(Primal, i);
               sol._redCost[i] -= objRational(i);
               sol._redCost[i] *= -1;
            }
            else
               sol._redCost[i] = _workSol._redCost[i];
         }
         sol._hasDual = true;
      }

   TERMINATE:
      // stop rational solving time
      _statistics->rationalTime->stop();
      return;
   }



   /// computes and loads rational QP-basis into rationalLUSolver
   void SoPlex::_computeAndLoadRationalQPBasis(bool& error, DataArray< SPxSolver::VarStatus >& basisStatusRows, DataArray< SPxSolver::VarStatus >& basisStatusCols, DVectorRational& righthandside, int& matrixdim, int& numBasicCols, int& numBasicRows)
   {
      DataArray< DSVectorRational > matrix(matrixdim);
      DataArray< const SVectorRational* > pntr_matrix(matrixdim);

      _rationalLUSolverBind.reSize(basisStatusCols.size()+basisStatusRows.size());

      // compute upper left part of basis system Q_BB basic rows and columns of Hessian and corresponding rhs (c_B - Q_BN l_N)
      int j = 0;
      int l;
      for( int i = 0; i < basisStatusCols.size(); i++ )
      {
         if( basisStatusCols[i] == SPxSolver::BASIC && j < matrixdim )
         {
	    _rationalLUSolverBind[j] = i;
            righthandside[j] = objRational(i);
            l = 0;
            for( int k = 0; k < basisStatusCols.size(); k++ )
	    {
	       if( basisStatusCols[k] == SPxSolver::BASIC && l < matrixdim )
	       {
                  if( j == 0 )
                     new(&matrix[l]) DSVectorRational();
	          if( rowcolVectorRational(i).number(k) >= 0 )
                  {
		     matrix[j].add(l,-rowcolVectorRational(i)[k]);
                     if( i != k )
		     {
		        assert(k > i);
                        matrix[l].add(j,-rowcolVectorRational(i)[k]);
                     }
                  }
                  l++;
	       }
	       else if( basisStatusCols[k] == SPxSolver::ON_LOWER )
	       {
	          if( rowcolVectorRational(i).number(k) >= 0 )
                     righthandside[j] += lowerRational(k) * rowcolVectorRational(i)[k];
                  else if( rowcolVectorRational(k).number(i) >= 0 )
                     righthandside[j] += lowerRational(k) * rowcolVectorRational(k)[i];
	       }
	       else if( basisStatusCols[k] == SPxSolver::ON_UPPER )
	       {
	          if( rowcolVectorRational(i).number(k) >= 0 )
                     righthandside[j] += upperRational(k) * rowcolVectorRational(i)[k];
                  else if( rowcolVectorRational(k).number(i) >= 0 )
                     righthandside[j] += upperRational(k) * rowcolVectorRational(k)[i];
	       }
	       else if( basisStatusCols[k] == SPxSolver::ZERO )
               {}
               else if( basisStatusCols[k] == SPxSolver::FIXED )
               {
                  assert(lowerRational(k) == upperRational(k));
	          if( rowcolVectorRational(i).number(k) >= 0 )
                     righthandside[j] += lowerRational(k) * rowcolVectorRational(i)[k];
                  else if( rowcolVectorRational(k).number(i) >= 0 )
                     righthandside[j] += lowerRational(k) * rowcolVectorRational(k)[i];
               }
               else if( basisStatusCols[k] == SPxSolver::UNDEFINED )
               {
                  MSG_ERROR( std::cerr << "Undefined basis status of column in rational factorization.\n" );
                  error = true;
                  return;
               }
               else
               {
                  assert(basisStatusCols[k] == SPxSolver::BASIC);
                  MSG_ERROR( std::cerr << "Too many basic columns in rational factorization (QP objective, row).\n" );
                  error = true;
                  return;
               }
	    }
            pntr_matrix[j] = &matrix[j];
            j++;
         }
         else if( basisStatusCols[i] == SPxSolver::ON_LOWER || basisStatusCols[i] == SPxSolver::ON_UPPER || basisStatusCols[i] == SPxSolver::ZERO || basisStatusCols[i] == SPxSolver::FIXED )
         {}
         else if( basisStatusCols[i] == SPxSolver::UNDEFINED )
         {
            MSG_ERROR( std::cerr << "Undefined basis status of column in rational factorization.\n" );
            error = true;
            return;
         }
         else
         {
            assert(basisStatusCols[i] == SPxSolver::BASIC);
            MSG_ERROR( std::cerr << "Too many basic columns in rational factorization (QP objective, column).\n" );
            error = true;
            return;
         }
      }
      numBasicCols = j;

      // compute symmetric of diagonal part of rows (without basic cols) and corresponding rhs (b - A_.N * l_N)
      for( int i = 0; i < basisStatusRows.size(); i++ )
      {
         new(&matrix[j]) DSVectorRational();
	 _rationalLUSolverBind[j] = -1 - i;
         if( basisStatusRows[i] == SPxSolver::ON_LOWER )
            righthandside[j] = lhsRational(i);
         else if( basisStatusRows[i] == SPxSolver::ON_UPPER )
            righthandside[j] = rhsRational(i);
         else if( basisStatusRows[i] == SPxSolver::ZERO )
            righthandside[j] = 0;
         else if( basisStatusRows[i] == SPxSolver::FIXED )
         {
            assert(lhsRational(i) == rhsRational(i));
            righthandside[j] = lhsRational(i);
         }
         else if( basisStatusRows[i] == SPxSolver::UNDEFINED )
         {
            MSG_ERROR( std::cerr << "Undefined basis status of row in rational factorization.\n" );
            error = true;
            return;
         }
         else
         {
            assert(basisStatusRows[i] == SPxSolver::BASIC);
            MSG_ERROR( std::cerr << "Basic row in rational factorization.\n" );
            error = true;
            return;
         }

         l = 0;
	 for( int k = 0; k < basisStatusCols.size(); k++ )
	 {
	    if( basisStatusCols[k] == SPxSolver::BASIC && j < matrixdim)
	    {
	       if( rowVectorRational(i)[k] != 0 )
               {
                  matrix[j].add(l,rowVectorRational(i)[k]);
                  matrix[l].add(j,matrix[j][l]);
               }
	       l++;
	    }
	    else if( basisStatusCols[k] == SPxSolver::ON_LOWER )
	    {
	       righthandside[j] -= lowerRational(k) * rowVectorRational(i)[k];
	    }
	    else if( basisStatusCols[k] == SPxSolver::ON_UPPER )
	    {
	       righthandside[j] -= upperRational(k) * rowVectorRational(i)[k];
	    }
	    else if( basisStatusCols[k] == SPxSolver::ZERO )
	    {}
	    else if( basisStatusCols[k] == SPxSolver::FIXED )
	    {
	       assert(lowerRational(i) == upperRational(i));
	       righthandside[j] -= lowerRational(k) * rowVectorRational(i)[k];
	    }
	    else if( basisStatusCols[k] == SPxSolver::UNDEFINED )
	    {
	       MSG_ERROR( std::cerr << "Undefined basis status of column in rational factorization.\n" );
	       error = true;
               return;
	    }
            else
            {
               assert(basisStatusCols[k] == SPxSolver::BASIC);
               MSG_ERROR( std::cerr << "Too many basic columns in rational factorization (constraint matrix).\n" );
               error = true;
               return;
            }
	 }
         pntr_matrix[j] = &matrix[j];
         j++;
      }
      _rationalLUSolverBind.reSize(j);
      assert( j - numBasicCols == basisStatusRows.size());

      if( numBasicCols < numRowsRational() )
      {
         MSG_ERROR( std::cerr << "Too few basic columns in rational factorization.\n" );
         error = true;
         return;
      }

      assert(matrixdim >= j);
      matrixdim = j;
      matrix.reSize(matrixdim);
      pntr_matrix.reSize(matrixdim);
      righthandside.reDim(matrixdim);
      righthandside.reSize(matrixdim+1);

      // load rational basis matrix
      if( realParam(SoPlex::TIMELIMIT) < realParam(SoPlex::INFTY) )
         _rationalLUSolver.setTimeLimit(realParam(SoPlex::TIMELIMIT) - _statistics->solvingTime->time());
      else
         _rationalLUSolver.setTimeLimit(-1.0);
      _rationalLUSolver.load(pntr_matrix.get_ptr(), matrixdim);

      // record statistics
      _statistics->luFactorizationTimeRational += _rationalLUSolver.getFactorTime();
      _statistics->luFactorizationsRational += _rationalLUSolver.getFactorCount();
      _rationalLUSolver.resetCounters();
   }
} // namespace soplex
#endif
#else
#include "soplex.h"

namespace soplex
{
   /// solves bQP with iterative refinement
   SPxSolver::Status SoPlex::solveQP()
   {
      MSG_ERROR( std::cerr << "Error: Cannot solve QPs without being linked to a QP solcer. Compile with QPOASES=true or QORE=true.\n" );
      return status();
   }
}
#endif
