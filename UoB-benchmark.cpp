#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include "ParmParse.H"
#include "LoadBalance.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "BRMeshRefine.H"
#include "MultilevelLinearOp.H"
#include "VCAMRPoissonOp2.H"
#include "RelaxSolver.H"
#include "AMRIO.H"

#include "UsingNamespace.H"

struct VCPoissonParameters
{

  int nLevels, nCoarsestCells;
  Real coarsestDx;
  ProblemDomain coarsestDomain;
  Vector<int> refRatio;
  int verbosity;
  Real alpha, beta;
  int coefficient_average_type;
  int maxBoxSize, blockFactor, nestingRadius;
  Real fillRatio;


  VCPoissonParameters()
  {
    ParmParse pp;
    
    nLevels = 2;
    pp.query("n_levels",nLevels);

    int n_cells = 128;
    pp.query("coarsest_n_cells", n_cells);
    nCoarsestCells = n_cells;
    n_cells /= 2;
    Box domBox(IntVect::Unit*(-n_cells), IntVect::Unit*(n_cells-1));
    bool periodic[CH_SPACEDIM];
    for (int dir = 0; dir < SpaceDim; dir++) 
      { 
	periodic[dir]= false; 
      }
    coarsestDomain = ProblemDomain(domBox, periodic );
    coarsestDx = 1.0 / Real(n_cells);
    
    coefficient_average_type = 1;
    alpha = 1.0;
    beta = 1.0;

    for (int lev = 0; lev <= nLevels; lev++)
      {
	refRatio.push_back(2);
      }

    maxBoxSize = 32;
    pp.query("max_box_size", maxBoxSize);

    blockFactor = 8;
    pp.query("block_factor", blockFactor);

    fillRatio = 0.95;
    pp.query("fill_ratio", fillRatio);

    nestingRadius = 1;
    pp.query("nesting_radius", nestingRadius);

    verbosity = 4;
  }

};

ostream& operator<<(ostream& os, VCPoissonParameters p)
{
  return os << "n_levels = " << p.nLevels << endl 
	    << "coarsest_n_cells = " << p.nCoarsestCells << endl;
    
}

extern
AMRLevelOpFactory<LevelData<FArrayBox> >*
defineOperatorFactory(
                      const Vector<DisjointBoxLayout>&             a_grids,
                      const Vector<ProblemDomain>&                 a_vectDomain,
                      Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_aCoef,
                      Vector<RefCountedPtr<LevelData<FluxBox> > >& a_bCoef,
                      const VCPoissonParameters&                     a_params)
{
  ParmParse pp2;

  VCAMRPoissonOp2Factory* opFactory = new VCAMRPoissonOp2Factory;
  //Homogeneous Dirichlett BCs
  BCHolder bc(ConstDiriNeumBC(IntVect::Unit, RealVect::Zero,
  			      IntVect::Unit, RealVect::Zero));
  opFactory->define(a_params.coarsestDomain,
                    a_grids,
                    a_params.refRatio,
                    a_params.coarsestDx,
                    bc,
                    a_params.alpha,
                    a_aCoef,
                    a_params.beta,
                    a_bCoef);

  if (a_params.coefficient_average_type >= 0)
    {
      opFactory->m_coefficient_average_type
        = a_params.coefficient_average_type;
    }

  return (AMRLevelOpFactory<LevelData<FArrayBox> >*) opFactory;
}


/******/
void poissonSolve(Vector<LevelData<FArrayBox>* >& a_phi,
                 Vector<LevelData<FArrayBox>* >& a_rhs,
                 const Vector< DisjointBoxLayout >&   a_grids,
                 const VCPoissonParameters&  a_params)
{
  CH_TIME("poissonSolve");
  ParmParse pp;

  int nlevels = a_params.nLevels;
  a_phi.resize(nlevels);
  a_rhs.resize(nlevels);
  Vector<RefCountedPtr<LevelData<FArrayBox> > > aCoef(nlevels);
  Vector<RefCountedPtr<LevelData<FluxBox> > > bCoef(nlevels);
  Vector<ProblemDomain> vectDomains(nlevels);
  Vector<RealVect> vectDx(nlevels);

  RealVect dxLev = RealVect::Unit;
  dxLev *=a_params.coarsestDx;
  ProblemDomain domLev(a_params.coarsestDomain);
  for (int ilev = 0; ilev < nlevels; ilev++)
    {
      a_rhs[ilev] = new LevelData<FArrayBox>(a_grids[ilev], 1,IntVect::Zero);
      a_phi[ilev] = new LevelData<FArrayBox>(a_grids[ilev], 1,IntVect::Unit);
      aCoef[ilev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(a_grids[ilev], 1, IntVect::Zero));
      bCoef[ilev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(a_grids[ilev], 1, IntVect::Zero));
      vectDomains[ilev] = domLev;
      vectDx[ilev] = dxLev;

      for (DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          (*a_phi[ilev])[dit()].setVal(0.);
	  (*a_rhs[ilev])[dit()].setVal(0.);
	  const Box& b = a_grids[ilev][dit];
	  for (BoxIterator bit(b); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      Real r = 0.0;
	      for (int dir = 0; dir < SpaceDim; dir++)
		{
		  Real x = ( Real(iv[dir])+0.5 ) * dxLev[dir];
		  r += x*x;
		}
	      r = std::sqrt(r) - 0.5;
	      (*a_rhs[ilev])[dit()](iv) = exp(-4.0*r*r);
	    }
	  (*aCoef[ilev])[dit()].setVal(1.0e-10);
	  for (int dir = 0; dir < SpaceDim; dir++)
	    {
	     
	      (*bCoef[ilev])[dit()][dir].setVal(1.0);
	    }
	}
      //prepare dx, domain for next level
      dxLev /=      a_params.refRatio[ilev];
      domLev.refine(a_params.refRatio[ilev]);
    }

  // set up solver
  RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > > opFactory
    = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >
    (defineOperatorFactory(a_grids, vectDomains, aCoef, bCoef, a_params));

  int lBase = 0;
  MultilevelLinearOp<FArrayBox> mlOp;
  int nMGIter = 1;
  pp.query("nMGIterations", nMGIter);
  mlOp.m_num_mg_iterations = nMGIter;
  int nMGSmooth = 4;
  pp.query("nMGsmooth", nMGSmooth);
  mlOp.m_num_mg_smooth = nMGSmooth;
  int preCondSolverDepth = -1;
  pp.query("preCondSolverDepth", preCondSolverDepth);
  mlOp.m_preCondSolverDepth = preCondSolverDepth;

  Real tolerance = 0.0;
  pp.query("tolerance", tolerance);

  int max_iter = 2;
  pp.query("max_iterations", max_iter);

  mlOp.define(a_grids, a_params.refRatio, vectDomains,
              vectDx, opFactory, lBase);

  RelaxSolver<Vector<LevelData<FArrayBox>* > > solver;

  bool homogeneousBC = false;
  solver.define(&mlOp, homogeneousBC);
  solver.m_verbosity = a_params.verbosity;
  solver.m_normType = 0;
  solver.m_eps = tolerance;
  solver.m_imax = max_iter;

  solver.solve(a_phi, a_rhs);

}

string ioName(const VCPoissonParameters& a_params, 
	      std::string a_suffix, int a_rank, int a_lev, int a_box)
{

  ostringstream os;
  os << "uob-benchmark" 
     << "-coarse-" << a_params.nCoarsestCells
     << "-nlev-" << a_params.nLevels
     << "-nrank-" << numProc();
  if (a_rank >= 0)
    os << "-rank-" << a_rank;
  if (a_lev >= 0)
    os << "-lev-" << a_lev;
  if (a_box >= 0)
    os << "-box-" << a_box;
  os << "." <<  CH_SPACEDIM << "d." <<  a_suffix;
  return os.str();

}
void ioDataManyFile(const Vector<LevelData<FArrayBox>* >&   a_phi,
			const Vector< DisjointBoxLayout >&  a_grids,
		    const VCPoissonParameters&  a_params, bool a_out)
{
  CH_TIME("ioDataManyFile");
  
  Interval interval(0,0);

  int n_buf = sizeof(Real);
  for (int dir = 0; dir < SpaceDim; dir++)  
    {
      n_buf *= (4+a_params.maxBoxSize);
    }
  char* buf = new char[n_buf];

  for (int lev = 0; lev < a_grids.size(); lev++)
    {
      int boxid = 0;
      for (DataIterator dit(a_grids[lev]); dit.ok(); ++dit, ++boxid)
	{
	  FArrayBox& fab = (*a_phi[lev])[dit];
	  int n = fab.size(fab.box(),interval);
	  if (n > n_buf)
	    {
	      //shouldn't really happen, but...
	      n_buf = n;
	      if (buf) delete[] buf;
	      buf = new char[n_buf];
	    }
	  if (buf)
	    {
	      if (a_out)
		{
		  ofstream os(ioName( a_params, "bin", procID(),  lev, boxid).c_str(),ios::binary);
		  fab.linearOut(buf, fab.box(),interval);
		  os.write(buf, n);
		  os.close();
		}
	      else
		{
		  ifstream is(ioName( a_params, "bin",  procID(), lev, boxid).c_str(),ios::binary);
		  is.read(buf, n);
		  fab.linearIn(buf, fab.box(),interval);
		  is.close();
		}
	    }
	  
	}
    }
  if (buf) delete[] buf;
}



void outputDataHDF5(const Vector<LevelData<FArrayBox>* >&   a_phi,
		    const Vector< DisjointBoxLayout >&       a_grids,
		    const VCPoissonParameters&                 a_params)
{
#ifdef CH_USE_HDF5
  CH_TIME("outputDataHDF5");
  string fileName = ioName( a_params, "hdf5",-1,-1,-1);
  Vector<string> phiNames(1,"phi");  
  Real fakeTime = 1.0;
  Real fakeDt = 1.0;
  WriteAMRHierarchyHDF5(fileName, a_grids,
			a_phi, phiNames,
			a_params.coarsestDomain.domainBox(),
			a_params.coarsestDx,
			fakeDt, fakeTime,
			a_params.refRatio,
			a_params.nLevels);

#endif
}


void refineGrids(Vector<DisjointBoxLayout>& a_grids, 
		 const VCPoissonParameters& a_param, 
		 int a_finest_level)
{
 
  a_grids.resize(a_finest_level + 1);

 
  // tags parts of the mesh for refinement; 
  // convert current Vector<DisjointBoxLayout> into Vector<Vector<Box> >
  Vector<IntVectSet> tagSetVect(a_finest_level+1);
  Vector<Vector<Box> > oldGrids(a_finest_level+1);
  Vector<Vector<Box> > newGrids;
  RealVect dxLev = RealVect::Unit;
  dxLev *=a_param.coarsestDx;
  Real refineRadius = 0.125;

  for (int lev = 0; lev < a_finest_level; lev++)
    {
      IntVectSet& tagSet = tagSetVect[lev];
      const DisjointBoxLayout& dbl = a_grids[lev];
      oldGrids[lev].resize(dbl.size());
      LayoutIterator lit = dbl.layoutIterator();
      int boxIndex = 0;
     
      for (lit.begin(); lit.ok(); ++lit, ++boxIndex) 
	{
	  const Box& b = dbl[lit()];
	  oldGrids[lev][boxIndex] = b;
	  for (BoxIterator bit(b); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      Real r = 0.0;
	      for (int dir = 0; dir < SpaceDim; dir++)
		{
		  Real x = ( Real(iv[dir])+0.5 ) * dxLev[dir];
		  r += x*x;
		}
	      r = std::sqrt(r);
	      if ( Abs( r - 0.5) < refineRadius ) tagSet |= iv;
	    }
	}
      dxLev /=      a_param.refRatio[lev];
      refineRadius /=  a_param.refRatio[lev];
    }
 
  //Mesh generator
  BRMeshRefine meshrefine(a_param.coarsestDomain, a_param.refRatio,
			  a_param.fillRatio, a_param.blockFactor,  
			  a_param.nestingRadius, a_param.maxBoxSize);
  
  int fin = meshrefine.regrid(newGrids, tagSetVect, 0, a_finest_level, oldGrids);

  ProblemDomain domain = a_param.coarsestDomain;
  for (int lev = 1; lev <= a_finest_level; lev++)
    {
      domain.refine( a_param.refRatio[lev] );
      Vector<int> procID(newGrids[lev].size());
      LoadBalance(procID, newGrids[lev]);
      a_grids[lev] = DisjointBoxLayout(newGrids[lev], procID, domain);
    }


}

void setGrids(Vector<DisjointBoxLayout>& a_grids,  
	      const VCPoissonParameters& a_param)
{
  CH_TIME("setGrids");
  //coarsest level
  Vector<Box>  baseBoxes;
  domainSplit(a_param.coarsestDomain, baseBoxes, 
	      a_param.maxBoxSize, a_param.blockFactor);
  Vector<int> procAssign(baseBoxes.size());
  LoadBalance(procAssign,baseBoxes);
  a_grids.push_back(DisjointBoxLayout(baseBoxes, procAssign, a_param.coarsestDomain));
  
  for (int finest_level = 1; finest_level < a_param.nLevels; finest_level++)
    {
      refineGrids(a_grids, a_param,finest_level );
    }

}
time_t initial_time;
void logTime(string a_item)
{
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  pout() << "time: " << std::difftime(std::time(0) , initial_time) <<  " clock: " << float(std::clock())/CLOCKS_PER_SEC << " " << a_item << endl;
}

int main(int argc, char* argv[])
{
  int status = 0;
  setenv("CH_TIMER","1",1);
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  // Scoping trick
  {
#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif
    int rank, nrank;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
    MPI_Comm_size(Chombo_MPI::comm, &nrank);
#else
    rank=0;
    nrank=1;
#endif
    initial_time = std::time(0);
    // command line params only - no input file for now
    ParmParse pp(argc-1,argv+1,NULL,NULL);
    
    VCPoissonParameters param;
    setPoutBaseName(ioName(param, "pout", procID(),  -1, -1));
    pout() << param;

    pout() << "n_rank = " <<  nrank << endl;

    Vector<DisjointBoxLayout> grids;
    Vector<LevelData<FArrayBox>* > phi(param.nLevels, NULL);
    Vector<LevelData<FArrayBox>* > rhs(param.nLevels, NULL);

    setGrids(grids,  param);
    logTime("end initialize");

    poissonSolve(phi, rhs, grids,  param);
    logTime("end poissonSolve");

    bool out(true);
    ioDataManyFile(phi, grids, param, out);
    logTime("end outputDataManyFile");

    ioDataManyFile(phi, grids, param, !out);
    logTime("end inputDataManyFile");

    outputDataHDF5(phi, grids, param);
    logTime("end outputDataHDF5");
    

    // clear memory
    for (int level = 0; level<phi.size(); level++)
      {
        if (phi[level] != NULL)
          {
            delete phi[level];
            phi[level] = NULL;
          }
        if (rhs[level] != NULL)
          {
            delete rhs[level];
            rhs[level] = NULL;
          }
      }

    CH_TIMER_REPORTNAME(pout(), "poissonSolve"); pout() << endl;
    CH_TIMER_REPORTNAME(pout(), "ioDataManyFile"); pout() << endl;
    CH_TIMER_REPORTNAME(pout(), "outputDataHDF5"); pout() << endl;
  }// End scoping trick

CH_TIMER_REPORT();

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return status;
}
