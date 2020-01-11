/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * kmetis.c
 *
 * This is the entry point of Moc_PARMETIS_PartGraphKway
 *
 * Started 10/19/96
 * George
 *
 * $Id: kmetis_withtimers.c,v 1.1 2007/02/12 18:46:50 kddevin Exp $
 *
 */

#include <parmetislib.h>

#define KDDKDD
#ifdef KDDKDD
#include "zoltan_timer.h"
struct Zoltan_Timer *zt;
int zt0, zt1, zt2, zt3, zt4, zt5;
#endif /* KDDKDD */

/***********************************************************************************
* This function is the entry point of the parallel k-way multilevel partitionioner. 
* This function assumes nothing about the graph distribution.
* It is the general case.
************************************************************************************/
void ParMETIS_V3_PartKway(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
              idxtype *adjwgt, int *wgtflag, int *numflag, int *ncon, int *nparts, 
	      float *tpwgts, float *ubvec, int *options, int *edgecut, idxtype *part, 
	      MPI_Comm *comm)
{
  int h, i;
  int nvtxs = -1, npes, mype;
  CtrlType ctrl;
  WorkSpaceType wspace;
  GraphType *graph;
  float avg, maximb, *mytpwgts;
  int moptions[10];
  int seed, dbglvl = 0;
  int iwgtflag, inumflag, incon, inparts, ioptions[10];
  float *itpwgts, iubvec[MAXNCON];

#ifdef KDDKDD
static int kdd_firsttime=1;
#endif /* KDDKDD */
  MPI_Comm_size(*comm, &npes);
  MPI_Comm_rank(*comm, &mype);

#ifdef KDDKDD
  if (kdd_firsttime) {
    zt = Zoltan_Timer_Create(ZOLTAN_TIME_WALL);
    zt0 = Zoltan_Timer_Init(zt, 1, "SETUP");
    zt1 = Zoltan_Timer_Init(zt, 1, "MATCH_AND_COARSEN");
    zt2 = Zoltan_Timer_Init(zt, 1, "COARSE_PARTITION");
    zt3 = Zoltan_Timer_Init(zt, 1, "REFINE");
    zt4 = Zoltan_Timer_Init(zt, 1, "REMAP");
    zt5 = Zoltan_Timer_Init(zt, 1, "CLEANUP");
    kdd_firsttime=0;
  }
  ZOLTAN_TIMER_START(zt, zt0, MPI_COMM_WORLD);
#endif /* KDDKDD */

  /********************************/
  /* Try and take care bad inputs */
  /********************************/
  if (options != NULL && options[0] == 1)
    dbglvl = options[PMV3_OPTION_DBGLVL];

  CheckInputs(STATIC_PARTITION, npes, dbglvl, wgtflag, &iwgtflag, numflag, &inumflag, ncon, 
              &incon, nparts, &inparts, tpwgts, &itpwgts, ubvec, iubvec, NULL, NULL, 
	      options, ioptions, part, comm);


  /*********************************/
  /* Take care the nparts = 1 case */
  /*********************************/
  if (inparts <= 1) {
    idxset(vtxdist[mype+1]-vtxdist[mype], 0, part); 
    *edgecut = 0;
    return;
  }

  /******************************/
  /* Take care of npes = 1 case */
  /******************************/
  if (npes == 1 && inparts > 1) {
    moptions[0] = 0;
    nvtxs = vtxdist[1];

    if (incon == 1) {
      METIS_WPartGraphKway(&nvtxs, xadj, adjncy, vwgt, adjwgt, &iwgtflag, &inumflag, 
            &inparts, itpwgts, moptions, edgecut, part);
    }
    else {
      /* ADD: this is because METIS does not support tpwgts for all constraints */
      mytpwgts = fmalloc(inparts, "mytpwgts");
      for (i=0; i<inparts; i++)
        mytpwgts[i] = itpwgts[i*incon];

      moptions[7] = -1;
      METIS_mCPartGraphRecursive2(&nvtxs, &incon, xadj, adjncy, vwgt, adjwgt, &iwgtflag, 
            &inumflag, &inparts, mytpwgts, moptions, edgecut, part);

      free(mytpwgts);
    }
 
    return;
  }


  if (inumflag == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 1);

  /*****************************/
  /* Set up control structures */
  /*****************************/
  if (ioptions[0] == 1) {
    dbglvl = ioptions[PMV3_OPTION_DBGLVL];
    seed = ioptions[PMV3_OPTION_SEED];
  }
  else {
    dbglvl = GLOBAL_DBGLVL;
    seed = GLOBAL_SEED;
  }
  SetUpCtrl(&ctrl, inparts, dbglvl, *comm);
  ctrl.CoarsenTo = amin(vtxdist[npes]+1, 25*incon*amax(npes, inparts));
  ctrl.seed = (seed == 0) ? mype : seed*mype;
  ctrl.sync = GlobalSEMax(&ctrl, seed);
  ctrl.partType = STATIC_PARTITION;
  ctrl.ps_relation = -1;
  ctrl.tpwgts = itpwgts;
  scopy(incon, iubvec, ctrl.ubvec);

  graph = Moc_SetUpGraph(&ctrl, incon, vtxdist, xadj, vwgt, adjncy, adjwgt, &iwgtflag);

  PreAllocateMemory(&ctrl, graph, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  /*******************************************/
  /* Check for funny cases                   */
  /* 	-graph with no edges                 */
  /* 	-graph with self edges               */
  /* 	-graph with poor vertex distribution */
  /* 	-graph with less than 2*npe nodes    */
  /*******************************************/
#ifdef KDDKDD
  ZOLTAN_TIMER_STOP(zt, zt0, MPI_COMM_WORLD);
#endif /* KDDKDD */
  if (vtxdist[npes] < SMALLGRAPH || vtxdist[npes] < npes*20 || GlobalSESum(&ctrl, graph->nedges) == 0) {
    IFSET(ctrl.dbglvl, DBG_INFO, rprintf(&ctrl, "Partitioning a graph of size %d serially\n", vtxdist[npes]));
    PartitionSmallGraph(&ctrl, graph, &wspace);
  }
  else {
    /***********************/
    /* Partition the graph */
    /***********************/
    Moc_Global_Partition(&ctrl, graph, &wspace);
#ifdef KDDKDD
    ZOLTAN_TIMER_START(zt, zt4, MPI_COMM_WORLD);
#endif /* KDDKDD */
    ParallelReMapGraph(&ctrl, graph, &wspace);
#ifdef KDDKDD
    ZOLTAN_TIMER_STOP(zt, zt4, MPI_COMM_WORLD);
#endif /* KDDKDD */
  }

#ifdef KDDKDD
  ZOLTAN_TIMER_START(zt, zt5, MPI_COMM_WORLD);
#endif /* KDDKDD */
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));

  idxcopy(graph->nvtxs, graph->where, part);
  *edgecut = graph->mincut;

  /*******************/
  /* Print out stats */
  /*******************/
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));

  if (ctrl.dbglvl&DBG_INFO) {
    rprintf(&ctrl, "Final %d-way CUT: %6d \tBalance: ", inparts, graph->mincut);
    avg = 0.0;
    for (h=0; h<incon; h++) {
      maximb = 0.0;
      for (i=0; i<inparts; i++)
        maximb = amax(maximb, graph->gnpwgts[i*incon+h]/itpwgts[i*incon+h]);
      avg += maximb;
      rprintf(&ctrl, "%.3f ", maximb);
    }
    rprintf(&ctrl, "  avg: %.3f\n", avg/(float)incon);
  }

  GKfree((void **)&itpwgts, (void **)&graph->lnpwgts, (void **)&graph->gnpwgts, (void **)&graph->nvwgt, LTERM);
  FreeInitialGraphAndRemap(graph, iwgtflag);
  FreeWSpace(&wspace);
  FreeCtrl(&ctrl);

  if (inumflag == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 0);
#ifdef KDDKDD
  ZOLTAN_TIMER_STOP(zt, zt5, MPI_COMM_WORLD);
  Zoltan_Timer_PrintAll(zt, 0, MPI_COMM_WORLD, stdout);
#endif /* KDDKDD */

}



/*************************************************************************
* This function is the driver to the multi-constraint partitioning algorithm.
**************************************************************************/
void Moc_Global_Partition(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace)
{
  int i, ncon, nparts;
  float ftmp, ubavg, lbavg, lbvec[MAXNCON];
 
#ifdef KDDKDD
ZOLTAN_TIMER_START(zt, zt1, MPI_COMM_WORLD);
#endif /* KDDKDD */
  ncon = graph->ncon;
  nparts = ctrl->nparts;
  ubavg = savg(graph->ncon, ctrl->ubvec);

  SetUp(ctrl, graph, wspace);

  if (ctrl->dbglvl&DBG_PROGRESS) {
    rprintf(ctrl, "[%6d %8d %5d %5d] [%d] [", graph->gnvtxs, GlobalSESum(ctrl, graph->nedges),
	    GlobalSEMin(ctrl, graph->nvtxs), GlobalSEMax(ctrl, graph->nvtxs), ctrl->CoarsenTo);
    for (i=0; i<ncon; i++)
      rprintf(ctrl, " %.3f", GlobalSEMinFloat(ctrl,graph->nvwgt[samin_strd(graph->nvtxs, graph->nvwgt+i, ncon)*ncon+i]));  
    rprintf(ctrl, "] [");
    for (i=0; i<ncon; i++)
      rprintf(ctrl, " %.3f", GlobalSEMaxFloat(ctrl, graph->nvwgt[samax_strd(graph->nvtxs, graph->nvwgt+i, ncon)*ncon+i]));  
    rprintf(ctrl, "]\n");
  }

  if (graph->gnvtxs < 1.3*ctrl->CoarsenTo ||
	(graph->finer != NULL &&
	graph->gnvtxs > graph->finer->gnvtxs*COARSEN_FRACTION)) {

    /* Done with coarsening. Find a partition */
#ifdef KDDKDD
ZOLTAN_TIMER_STOP(zt, zt1, MPI_COMM_WORLD);
ZOLTAN_TIMER_START(zt, zt2, MPI_COMM_WORLD);
#endif /* KDDKDD */
    graph->where = idxmalloc(graph->nvtxs+graph->nrecv, "graph->where");
    Moc_InitPartition_RB(ctrl, graph, wspace);

    if (ctrl->dbglvl&DBG_PROGRESS) {
      Moc_ComputeParallelBalance(ctrl, graph, graph->where, lbvec);
      rprintf(ctrl, "nvtxs: %10d, balance: ", graph->gnvtxs);
      for (i=0; i<graph->ncon; i++) 
        rprintf(ctrl, "%.3f ", lbvec[i]);
      rprintf(ctrl, "\n");
    }
#ifdef KDDKDD
ZOLTAN_TIMER_STOP(zt, zt2, MPI_COMM_WORLD);
#endif /* KDDKDD */

    /* In case no coarsening took place */
#ifdef KDDKDD
ZOLTAN_TIMER_START(zt, zt3, MPI_COMM_WORLD);
#endif /* KDDKDD */
    if (graph->finer == NULL) {
      Moc_ComputePartitionParams(ctrl, graph, wspace);
      Moc_KWayFM(ctrl, graph, wspace, NGR_PASSES);
    }
#ifdef KDDKDD
ZOLTAN_TIMER_STOP(zt, zt3, MPI_COMM_WORLD);
#endif /* KDDKDD */
  }
  else {
    Moc_GlobalMatch_Balance(ctrl, graph, wspace);

#ifdef KDDKDD
ZOLTAN_TIMER_STOP(zt, zt1, MPI_COMM_WORLD);
#endif /* KDDKDD */
    Moc_Global_Partition(ctrl, graph->coarser, wspace);
#ifdef KDDKDD
ZOLTAN_TIMER_START(zt, zt1, MPI_COMM_WORLD);
#endif /* KDDKDD */

    Moc_ProjectPartition(ctrl, graph, wspace);
    Moc_ComputePartitionParams(ctrl, graph, wspace);

    if (graph->ncon > 1 && graph->level < 3) {
      for (i=0; i<ncon; i++) {
        ftmp = ssum_strd(nparts, graph->gnpwgts+i, ncon);
        if (ftmp != 0.0)
          lbvec[i] = (float)(nparts) *
          graph->gnpwgts[samax_strd(nparts, graph->gnpwgts+i, ncon)*ncon+i]/ftmp;
        else
          lbvec[i] = 1.0;
      }
      lbavg = savg(graph->ncon, lbvec);

      if (lbavg > ubavg + 0.035) {
        if (ctrl->dbglvl&DBG_PROGRESS) {
          Moc_ComputeParallelBalance(ctrl, graph, graph->where, lbvec);
          rprintf(ctrl, "nvtxs: %10d, cut: %8d, balance: ", graph->gnvtxs, graph->mincut);
          for (i=0; i<graph->ncon; i++) 
            rprintf(ctrl, "%.3f ", lbvec[i]);
          rprintf(ctrl, "\n");
	}

        Moc_KWayBalance(ctrl, graph, wspace, graph->ncon);
      }
    }

#ifdef KDDKDD
ZOLTAN_TIMER_STOP(zt, zt1, MPI_COMM_WORLD);
ZOLTAN_TIMER_START(zt, zt3, MPI_COMM_WORLD);
#endif /* KDDKDD */
    Moc_KWayFM(ctrl, graph, wspace, NGR_PASSES);
#ifdef KDDKDD
ZOLTAN_TIMER_STOP(zt, zt3, MPI_COMM_WORLD);
ZOLTAN_TIMER_START(zt, zt1, MPI_COMM_WORLD);
#endif /* KDDKDD */

    if (ctrl->dbglvl&DBG_PROGRESS) {
      Moc_ComputeParallelBalance(ctrl, graph, graph->where, lbvec);
      rprintf(ctrl, "nvtxs: %10d, cut: %8d, balance: ", graph->gnvtxs, graph->mincut);
      for (i=0; i<graph->ncon; i++) 
        rprintf(ctrl, "%.3f ", lbvec[i]);
      rprintf(ctrl, "\n");
    }

    if (graph->level != 0)
      GKfree((void **)&graph->lnpwgts, (void **)&graph->gnpwgts, LTERM);
#ifdef KDDKDD
ZOLTAN_TIMER_STOP(zt, zt1, MPI_COMM_WORLD);
#endif /* KDDKDD */
  }

  return;
}


