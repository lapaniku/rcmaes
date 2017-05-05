//
//  cmaes_R.cpp
//  cma-es
//
//  Created by Андрей Лапаник on 3/8/17.
//  Copyright © 2017 RocketScience.ai. All rights reserved.
//

#include <Rcpp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "boundary_transformation.h"
#include "cmaes_interface.h"

using namespace Rcpp;

//#define DEBUG 1
#ifdef DEBUG
#define D if(1) 
#else
#define D if(0) 
#endif

// Utility
double* vectorToArrayP(std::vector<double> vec) {
  double* result;
  result = (double *)malloc(vec.size()*sizeof(double));
  std::copy(vec.begin(), vec.end(), result);
  return result;
}

long int* vectorToLongIntArrayP(std::vector<double> vec) {
  long int* result;
  std::vector<long int> vecLong(vec.begin(), vec.end());
  result = (long int *)malloc(vec.size()*sizeof(long int));
  std::copy(vecLong.begin(), vecLong.end(), result);
  return result;
}

double** matrixToArrayP(NumericMatrix M) {
  double** result;
  result = (double **)malloc(M.nrow()*sizeof(double*));
  for (unsigned i=0; i < M.nrow(); ++i) {
    double *row = (double *)malloc(M.ncol()*sizeof(double));
    result[i] = row;
    for (unsigned j = 0; j < M.ncol(); ++j) {
      result[i][j] = M(i, j); 
    }
  }
  return result;
}

NumericMatrix arrayPToMatrix(double** m, int rows, int cols) {
  NumericMatrix result(rows, cols);
  for (unsigned i=0; i<rows; ++i) {
    for (unsigned j=0; j<cols; ++j) {
      result(i, j) = m[i][j];
      D printf("Matrix[%d, %d]=%f\n", i, j, m[i][j]);
    }
  }
  return result;
}



List cmaesReadParaToList(cmaes_readpara_t sp) {
    D printf("The value of sp is: %p\n", (void *) &sp);
    D printf("The value of filename is: %p\n", &(sp.filename));
    D printf("The value of flgsupplemented is: %p\n", &(sp.flgsupplemented));
  
    List result = List::create();

    if(sp.filename != NULL) {
      std::string filename;
      filename = filename + std::string(sp.filename);
      result["filename"] = filename;
      D printf("Dimension: %s\n", filename.c_str());
    }
     
    D printf("flgsupplemented: %d\n", sp.flgsupplemented);
    result["flgsupplemented"] = sp.flgsupplemented;
    result["N"] = sp.N;
    result["seed"] = sp.seed;
    
    D printf("Dimension: %d\n", sp.N);

    std::vector<double> xstart(sp.xstart, sp.xstart + sp.N);
    result["xstart"] = xstart;
    
    if(sp.typicalX != NULL) {
      std::vector<double> typicalX(sp.typicalX, sp.typicalX + sp.N);
      result["typicalX"] = typicalX;
    }
    
    result["typicalXcase"] = sp.typicalXcase;
    
    std::vector<double> rgInitialStds(sp.rgInitialStds, sp.rgInitialStds + sp.N);
    result["rgInitialStds"] = rgInitialStds;
     
    // ???
    if(sp.rgDiffMinChange != NULL) {
      std::vector<double> rgDiffMinChange(sp.rgDiffMinChange, sp.rgDiffMinChange + sp.N);
      result["rgDiffMinChange"] = rgDiffMinChange;
    }
     
    result["stopMaxFunEvals"] = sp.stopMaxFunEvals;
    result["facmaxeval"] = sp.facmaxeval;
    result["stopMaxIter"] = sp.stopMaxIter;
    result["stStopFitness"] = List::create(Rcpp::Named("flg") = sp.stStopFitness.flg,
                                    Rcpp::Named("val") = sp.stStopFitness.val); 
    result["stopTolFun"] = sp.stopTolFun;
    result["stopTolFunHist"] = sp.stopTolFunHist;
    result["stopTolX"] = sp.stopTolX;
    result["stopTolUpXFactor"] = sp.stopTolUpXFactor;
    result["lambda"] = sp.lambda;
    result["mu"] = sp.mu;
    result["mucov"] = sp.mucov;
    result["mueff"] = sp.mueff;
    
    std::vector<double> weights(sp.weights, sp.weights + sp.mu);
    result["weights"] = weights;

    result["damps"] = sp.damps;
    result["cs"] = sp.cs;
    result["ccumcov"] = sp.ccumcov;
    result["ccov"] = sp.ccov;
    result["diagonalCov"] = sp.diagonalCov;
    result["updateCmode"] = List::create(Rcpp::Named("flgalways") = sp.updateCmode.flgalways,
                                      Rcpp::Named("modulo") = sp.updateCmode.modulo,
                                      Rcpp::Named("maxtime") = sp.updateCmode.maxtime);
    result["facupdateCmode"] = sp.facupdateCmode;

    std::string weigkey;
    weigkey = weigkey + std::string(sp.weigkey);
    result["weigkey"] = weigkey;

//    std::string resumefile;
//    resumefile = resumefile + std::string(sp.resumefile);
//    char resumefile[99];
//    const char **rgsformat;
//    void **rgpadr;
//    const char **rgskeyar;
//    double ***rgp2adr;
    
    result["n1para"] = sp.n1para;
    result["n1outpara"] = sp.n1outpara;
    result["n2para"] = sp.n2para;
    return result;
}

cmaes_readpara_t listToCmaesReadPara(List sp) {
    cmaes_readpara_t result;
    
    if(sp.containsElementNamed("filename")) {
      std::string filenameStr = sp["filename"];
      D printf("filename length: %lu\n", filenameStr.length());
      char *filename = new char[filenameStr.length() + 1];
      strcpy(filename, filenameStr.c_str());
      result.filename = filename;
    } else {
      D printf("filename is null\n");
      result.filename = NULL;
    }
    
    result.flgsupplemented = sp["flgsupplemented"];
    result.N = sp["N"];
    result.seed = sp["seed"];
    
    std::vector<double> xstartVec = sp["xstart"];
    D printf("xstart vec length: %lu\n", xstartVec.size());
    result.xstart = (double *)malloc(xstartVec.size()*sizeof(double));
    std::copy(xstartVec.begin(), xstartVec.end(), result.xstart);

    if(sp.containsElementNamed("typicalX")) {
      D printf("typicalX\n");
      std::vector<double> typicalXVec = sp["typicalX"];
      result.typicalX = (double *)malloc(typicalXVec.size()*sizeof(double));
      std::copy(typicalXVec.begin(), typicalXVec.end(), result.typicalX);
    } else {
      result.typicalX = NULL;
    }

    result.typicalXcase = sp["typicalXcase"];
    D printf("typicalXcase: %d\n", result.typicalXcase);
    
    std::vector<double> rgInitialStdsVec = sp["rgInitialStds"];
    D printf("rgInitialStds vec length: %lu\n", rgInitialStdsVec.size());
    result.rgInitialStds = (double *)malloc(rgInitialStdsVec.size()*sizeof(double));
    std::copy(rgInitialStdsVec.begin(), rgInitialStdsVec.end(), result.rgInitialStds);

    if(sp.containsElementNamed("rgDiffMinChange")) {
      std::vector<double> rgDiffMinChangeVec = sp["rgDiffMinChange"];
      result.rgDiffMinChange = (double *)malloc(rgDiffMinChangeVec.size()*sizeof(double));
      std::copy(rgDiffMinChangeVec.begin(), rgDiffMinChangeVec.end(), result.rgInitialStds);
    } else {
      result.rgDiffMinChange = NULL;
    }

    result.stopMaxFunEvals = sp["stopMaxFunEvals"];
    D printf("stopMaxFunEvals: %f\n", result.stopMaxFunEvals);
    result.facmaxeval = sp["facmaxeval"];
    D printf("facmaxeval: %f\n", result.facmaxeval);
    result.stopMaxIter = sp["stopMaxIter"];
    D printf("stopMaxIter: %f\n", result.stopMaxIter);
    if(sp.containsElementNamed("stStopFitness")) {
      List stStopFitness = sp["stStopFitness"];
      result.stStopFitness.flg = stStopFitness["flg"];
      D printf("stStopFitness.flg: %d\n", result.stStopFitness.flg);
      result.stStopFitness.val = stStopFitness["val"];
      D printf("stStopFitness.val: %f\n", result.stStopFitness.val);
    }
    result.stopTolFun = sp["stopTolFun"];
    D printf("stopTolFun: %f\n", result.stopTolFun);
    result.stopTolFunHist = sp["stopTolFunHist"];
    D printf("stopTolFunHist: %f\n", result.stopTolFunHist);
    result.stopTolUpXFactor = sp["stopTolUpXFactor"];
    D printf("stopTolUpXFactor: %f\n", result.stopTolUpXFactor);
    result.lambda = sp["lambda"];
    D printf("lambda: %d\n", result.lambda);
    result.mu = sp["mu"];
    D printf("mu: %d\n", result.mu);
    result.mucov = sp["mucov"];
    D printf("mucov: %f\n", result.mucov);
    result.mueff = sp["mueff"];
    D printf("mueff: %f\n", result.mueff);
    
    std::vector<double> weightsVec = sp["weights"];
    D printf("weights vec length: %lu\n", weightsVec.size());
    result.weights = (double *)malloc(weightsVec.size()*sizeof(double));
    std::copy(weightsVec.begin(), weightsVec.end(), result.weights);

    result.damps = sp["damps"];
    D printf("damps: %f\n", result.damps);
    result.cs = sp["cs"];
    D printf("cs: %f\n", result.cs);
    result.ccumcov = sp["ccumcov"];
    D printf("ccumcov: %f\n", result.ccumcov);
    result.ccov = sp["ccov"];
    D printf("ccov: %f\n", result.ccov);
    result.diagonalCov = as<double>(sp["diagonalCov"]);
    D printf("diagonalCov: %f\n", result.diagonalCov);
    
    if(sp.containsElementNamed("updateCmode")) {
      List updateCmode = sp["updateCmode"];
      result.updateCmode.flgalways = updateCmode["flgalways"];
      D printf("updateCmode.flgalways: %d\n", result.updateCmode.flgalways);
      result.updateCmode.modulo = updateCmode["modulo"];
      D printf("updateCmode.modulo: %f\n", result.updateCmode.modulo);
      result.updateCmode.maxtime = updateCmode["maxtime"];
      D printf("updateCmode.maxtime: %f\n", result.updateCmode.maxtime);
    }
    
    result.facupdateCmode = sp["facupdateCmode"];
    D printf("facupdateCmode: %f\n", result.facupdateCmode);
    
    std::string weigkeyStr = sp["weigkey"];
    char *weigkey = new char[weigkeyStr.length() + 1];
    strcpy(weigkey, weigkeyStr.c_str());
    result.weigkey = weigkey;
    
    result.n1para = sp["n1para"];
    result.n1outpara = sp["n1outpara"];
    result.n2para = sp["n2para"];
    
    return result;
}

cmaes_random_t listToCmaesRandom(List rand) {
  cmaes_random_t result;
  
  result.startseed = as<long int>(rand["startseed"]);
  result.aktseed = as<long int>(rand["aktseed"]);
  result.aktrand = as<long int>(rand["aktrand"]);
  
  std::vector<double> rgrandDouble = rand["rgrand"];

  std::vector<long int> rgrandLong(rgrandDouble.begin(), rgrandDouble.end());
  result.rgrand = (long int *)malloc(rgrandLong.size()*sizeof(long int));
  std::copy(rgrandLong.begin(), rgrandLong.end(), result.rgrand);

  result.flgstored = as<short>(rand["flgstored"]);
  result.hold = as<double>(rand["hold"]);

  return result;
}

List cmaesRandomToList(cmaes_random_t rand) {
  std::vector<long int> rgrand(rand.rgrand, rand.rgrand + 32);

  List result = List::create();
  result["startseed"] = rand.startseed; 
  result["aktseed"] = rand.aktseed;
  result["aktrand"] = rand.aktrand;
  result["rgrand"] = rgrand;
  result["flgstored"] = rand.flgstored;
  result["hold"] = rand.hold;

  return result;
}

cmaes_timings_t listToCmaesTimings(List timing) {
  cmaes_timings_t result;
  
  result.totaltime = as<double>(timing["totaltime"]);
  result.totaltotaltime = as<double>(timing["totaltotaltime"]);
  result.tictoctime = as<double>(timing["tictoctime"]);
  result.lasttictoctime = as<double>(timing["lasttictoctime"]);
  
  result.lastclock = as<double>(timing["lastclock"]);
  result.lasttime = as<double>(timing["lasttime"]);
  result.ticclock = as<double>(timing["ticclock"]);
  result.tictime = as<double>(timing["tictime"]);
    
  result.istic = as<short>(timing["istic"]);
  result.isstarted = as<short>(timing["isstarted"]);
  result.lastdiff = as<double>(timing["lastdiff"]);
  result.tictoczwischensumme = as<double>(timing["tictoczwischensumme"]);
    
  return result;
}

List cmaesTimingsToList(cmaes_timings_t timing) {
  List result = List::create();
  result["totaltime"] = timing.totaltime; 
  result["totaltotaltime"] = timing.totaltotaltime;
  result["tictoctime"] = timing.tictoctime;
  result["lasttictoctime"] = timing.lasttictoctime;
  
  result["lastclock"] = timing.lastclock;
  result["lasttime"] = timing.lasttime;
  result["ticclock"] = timing.ticclock;
  result["tictime"] = timing.tictime;

  result["istic"] = timing.istic;
  result["isstarted"] = timing.isstarted;
  result["lastdiff"] = timing.lastdiff;
  result["tictoczwischensumme"] = timing.tictoczwischensumme;
  
  return result;
}

cmaes_t listToCmaes(List cmaes) {

    cmaes_t result;
    
    std::string versionStr = as<std::string>(cmaes["version"]);
    char *version = (char *) malloc(versionStr.length()+1); 
    strncpy(version, versionStr.c_str(), versionStr.length()+1);
    result.version = version;
    D printf("Version: %s\n", result.version);

    D printf("Converting: sp\n");
    result.sp = listToCmaesReadPara(cmaes["sp"]);
    D printf("Done: sp\n");
    result.rand = listToCmaesRandom(cmaes["rand"]);
    D printf("Done: rand\n");
    
    result.sigma = as<double>(cmaes["sigma"]);
    
    result.rgxmean = vectorToArrayP(as<std::vector<double> >(cmaes["rgxmean"]));
    D printf("Done: rgxmean\n");
    result.rgxbestever = vectorToArrayP(as<std::vector<double> >(cmaes["rgxbestever"]));
    D printf("Done: rgxbestever\n");
    
    result.rgrgx = matrixToArrayP(as<NumericMatrix>(cmaes["rgrgx"]));
    D printf("Done: rgrgx\n");
    
    std::vector<double> indexDouble = cmaes["index"];
    
    std::vector<int> indexInt(indexDouble.begin(), indexDouble.end());
    result.index = (int *)malloc(indexInt.size()*sizeof(int));
    std::copy(indexInt.begin(), indexInt.end(), result.index);
    D printf("Done: index\n");
    
    // special trick for size in zero element
    std::vector<double> arFuncValueHist = as<std::vector<double> >(cmaes["arFuncValueHist"]);
    result.arFuncValueHist = (double *)malloc((arFuncValueHist.size()+1)*sizeof(double));
    result.arFuncValueHist[0] = arFuncValueHist.size();
    result.arFuncValueHist++;
    std::copy(arFuncValueHist.begin(), arFuncValueHist.end(), result.arFuncValueHist);

    D printf("Done: arFuncValueHist\n");
    
    result.flgIniphase = as<short>(cmaes["flgIniphase"]);
    result.flgStop = as<short>(cmaes["flgStop"]);

    result.chiN = as<double>(cmaes["chiN"]);

    result.C = matrixToArrayP(as<NumericMatrix>(cmaes["C"]));
    D printf("Done: C\n");
    result.B = matrixToArrayP(as<NumericMatrix>(cmaes["B"]));
    D printf("Done: B\n");
    result.rgD = vectorToArrayP(as<std::vector<double> >(cmaes["rgD"]));
    D printf("Done: rgD\n");
    
    result.rgpc = vectorToArrayP(as<std::vector<double> >(cmaes["rgpc"]));
    D printf("Done: rgpc\n");
    result.rgps = vectorToArrayP(as<std::vector<double> >(cmaes["rgps"]));
    D printf("Done: rgps\n");
    result.rgxold = vectorToArrayP(as<std::vector<double> >(cmaes["rgxold"]));
    D printf("Done: rgxold\n");
    result.rgout = vectorToArrayP(as<std::vector<double> >(cmaes["rgout"]));
    D printf("Done: rgout\n");
    result.rgBDz = vectorToArrayP(as<std::vector<double> >(cmaes["rgBDz"]));
    D printf("Done: rgBDz\n");
    result.rgdTmp = vectorToArrayP(as<std::vector<double> >(cmaes["rgdTmp"]));
    D printf("Done: rgdTmp\n");
    result.rgFuncValue = vectorToArrayP(as<std::vector<double> >(cmaes["rgFuncValue"]));
    D printf("Done: rgFuncValue\n");
    result.publicFitness = vectorToArrayP(as<std::vector<double> >(cmaes["publicFitness"]));
    D printf("Done: publicFitness\n");
    
    result.gen = as<double>(cmaes["gen"]);
    result.countevals = as<double>(cmaes["countevals"]);
    result.state = as<double>(cmaes["state"]);
    
    result.maxdiagC = as<double>(cmaes["maxdiagC"]);
    result.mindiagC = as<double>(cmaes["mindiagC"]);
    result.maxEW = as<double>(cmaes["maxEW"]);
    result.minEW = as<double>(cmaes["minEW"]);

    result.flgEigensysIsUptodate = as<short>(cmaes["flgEigensysIsUptodate"]);
    result.flgCheckEigen = as<short>(cmaes["flgCheckEigen"]);
    result.genOfEigensysUpdate = as<double>(cmaes["genOfEigensysUpdate"]);
    
    if(cmaes.containsElementNamed("eigenTimings")) {
      result.eigenTimings = listToCmaesTimings(cmaes["eigenTimings"]);
      D printf("Done: eigenTimings\n");
    }

    result.dMaxSignifKond = as<double>(cmaes["dMaxSignifKond"]);
    result.dLastMinEWgroesserNull = as<double>(cmaes["dLastMinEWgroesserNull"]);

    result.flgresumedone = as<short>(cmaes["flgresumedone"]);

    result.printtime = as<short>(cmaes["printtime"]);
    result.writetime = as<short>(cmaes["writetime"]);
    result.firstwritetime = as<short>(cmaes["firstwritetime"]);
    result.firstprinttime = as<short>(cmaes["firstprinttime"]);

    
    D printf("Done: listToCmaes\n");
    
    return result;
    
    //char sOutString[330]; /* 4x80 */
}

List cmaesToList(cmaes_t cmaes) {
  D printf("Start: cmaesToList\n");
  
    std::vector<double> rgpc(cmaes.rgpc, cmaes.rgpc + cmaes.sp.N);
    D printf("Done: rgpc\n");
    std::vector<double> rgps(cmaes.rgps, cmaes.rgps + cmaes.sp.N);
    D printf("Done: rgps\n");
    std::vector<double> rgdTmp(cmaes.rgdTmp, cmaes.rgdTmp + cmaes.sp.N + 1);
    D printf("Done: rgdTmp\n");
    std::vector<double> rgBDz(cmaes.rgBDz, cmaes.rgBDz + cmaes.sp.N);
    D printf("Done: rgBDz\n");
    std::vector<double> rgxmean(cmaes.rgxmean, cmaes.rgxmean + cmaes.sp.N + 1);
    D printf("Done: rgxmean\n");
    std::vector<double> rgxold(cmaes.rgxold, cmaes.rgxold + cmaes.sp.N + 1);
    D printf("Done: rgxold\n");
    std::vector<double> rgxbestever(cmaes.rgxbestever, cmaes.rgxbestever + cmaes.sp.N + 2);
    D printf("Done: rgxbestever\n");
    std::vector<double> rgout(cmaes.rgout, cmaes.rgout + cmaes.sp.N + 1);
    D printf("Done: rgout\n");
    std::vector<double> rgD(cmaes.rgD, cmaes.rgD + cmaes.sp.N);
    D printf("Done: rgD\n");
    std::vector<double> publicFitness(cmaes.publicFitness, cmaes.publicFitness + cmaes.sp.lambda);
    D printf("Done: publicFitness\n");
    std::vector<double> rgFuncValue(cmaes.rgFuncValue, cmaes.rgFuncValue + cmaes.sp.lambda);
    D printf("Done: rgFuncValue\n");
    std::vector<int> index(cmaes.index, cmaes.index + cmaes.sp.lambda);
    D printf("Done: index\n");
    std::vector<double> arFuncValueHist(cmaes.arFuncValueHist, cmaes.arFuncValueHist + (10 + (int)ceil(3.*10.*cmaes.sp.N/cmaes.sp.lambda)));
    D printf("Done: arFuncValueHist\n");
    
    D printf("Start: result\n");
    List result = List::create();
    
    std::string version;
    version = version + std::string(cmaes.version);
    result["version"] = version;
    D printf("Version: %s\n", cmaes.version);
    
    D printf("Converting: sp\n");
    result["sp"] = cmaesReadParaToList(cmaes.sp);
    D printf("Done: sp\n");
    
    result["rand"] = cmaesRandomToList(cmaes.rand);
    D printf("Done: rand\n");

    result["sigma"] = cmaes.sigma;
    result["rgxmean"] = rgxmean;
    result["rgxbestever"] = rgxbestever;
    
    D printf("Converting rgrgx\n");
    result["rgrgx"] = arrayPToMatrix(cmaes.rgrgx, cmaes.sp.lambda, cmaes.sp.N+1);
    /*
    for (int i = 0; i < cmaes.sp.lambda; ++i) {
      for (int j = 0; j < cmaes.sp.N+2; ++j) {
        D printf("rgrgx[%d][%d]=%f", i, j, cmaes.rgrgx[i][j]);
      }
    }
     */
    
    
    result["index"] = index;
    
    
    result["arFuncValueHist"] = arFuncValueHist;
    
    result["flgIniphase"] = cmaes.flgIniphase;
    result["flgStop"] = cmaes.flgStop;
    
    result["chiN"] = cmaes.chiN;
    D printf("Converting C\n");
    
    // cmaes_WriteToFile(&cmaes, "C", "C_conv.dat");
    
    result["C"] = arrayPToMatrix(cmaes.C, cmaes.sp.N, cmaes.sp.N);
    D printf("Converting B\n");
    result["B"] = arrayPToMatrix(cmaes.B, cmaes.sp.N, cmaes.sp.N);
    result["rgD"] = rgD;
    
    result["rgpc"] = rgpc;
    result["rgps"] = rgps;
    result["rgxold"] = rgxold;
    result["rgout"] = rgout;
    result["rgBDz"] = rgBDz;
    result["rgdTmp"] = rgdTmp;
    result["rgFuncValue"] = rgFuncValue;
    result["publicFitness"] = rgFuncValue;
    
    result["gen"] = cmaes.gen;
    result["countevals"] = cmaes.countevals;
    result["state"] = cmaes.state;
    
    result["maxdiagC"] = cmaes.maxdiagC;
    result["mindiagC"] = cmaes.mindiagC;
    result["maxEW"] = cmaes.maxEW;
    result["minEW"] = cmaes.minEW;
    
    std::string sOutString;
    sOutString = sOutString + std::string(cmaes.sOutString);
    result["sOutString"] = sOutString;
    
    result["flgEigensysIsUptodate"] = cmaes.flgEigensysIsUptodate;
    result["flgCheckEigen"] = cmaes.flgCheckEigen;
    result["genOfEigensysUpdate"] = cmaes.genOfEigensysUpdate;
    result["eigenTimings"] = cmaesTimingsToList(cmaes.eigenTimings);
    
    result["dMaxSignifKond"] = cmaes.dMaxSignifKond;
    result["dLastMinEWgroesserNull"] = cmaes.dLastMinEWgroesserNull;
    
    result["flgresumedone"] = cmaes.flgresumedone;
    
    result["printtime"] = cmaes.printtime;
    result["writetime"] = cmaes.writetime;
    result["firstwritetime"] = cmaes.firstwritetime;
    result["firstprinttime"] = cmaes.firstprinttime;
    
    return result;
}

List boundaryTranfrormToList(cmaes_boundary_transformation_t boundaries) {
  std::vector<double> lowerBounds(boundaries.lower_bounds, boundaries.lower_bounds + boundaries.len_of_bounds);
  std::vector<double> upperBounds(boundaries.upper_bounds, boundaries.upper_bounds + boundaries.len_of_bounds);
  std::vector<double> al(boundaries.al, boundaries.al + boundaries.len_of_bounds);
  std::vector<double> au(boundaries.au, boundaries.au + boundaries.len_of_bounds);
  
  
  List result = List::create();
  
  result["lower_bounds"] = lowerBounds;
  result["upper_bounds"] = upperBounds;
  result["len_of_bounds"] = boundaries.len_of_bounds;
  result["al"] = al;
  result["au"] = au;
  
  return result;
}

cmaes_boundary_transformation_t boundaryTransformationToStruct(List boundaries) {
  cmaes_boundary_transformation_t result;
  
  result.al = vectorToArrayP(as<std::vector<double> >(boundaries["al"]));
  result.au = vectorToArrayP(as<std::vector<double> >(boundaries["au"]));
  result.lower_bounds = vectorToArrayP(as<std::vector<double> >(boundaries["lower_bounds"]));
  result.upper_bounds = vectorToArrayP(as<std::vector<double> >(boundaries["upper_bounds"]));
  result.len_of_bounds = as<unsigned long>(boundaries["len_of_bounds"]);
  
  return result;
}

/*
 Testing section
 */

void printCmaesReadPara(cmaes_readpara_t t, long dimension) {
    printf("\nfilename: %s are:\n", t.filename);
    printf("N: %d\n",t.N);
    for (int i = 0; i < dimension; ++i)
        printf("%f\t", t.xstart[i]);
    printf("\n");
}

// [[Rcpp::export]]
void testCmaesReadPara(List sp, long dim) {
    
    cmaes_readpara_t t = listToCmaesReadPara(sp);
    printCmaesReadPara(t, dim);
}

// [[Rcpp::export]]
List testCmaesRand(List rand) {
  
  cmaes_random_t t = listToCmaesRandom(rand);
  return cmaesRandomToList(t);
}

// [[Rcpp::export]]
NumericMatrix testMatrix(NumericMatrix M, long rows, long cols) {
  double** m = matrixToArrayP(M);
  return arrayPToMatrix(m, rows, cols);
}

// [[Rcpp::export]]
List testCmaes(List cmaes) {
  cmaes_t cmaes_struct = listToCmaes(cmaes);
  return cmaesToList(cmaes_struct);
}

/* 
 * CMAES related
 */

//' Boundary transformation rules init
//' @param List of $max and 4min values
//' @return List of boundary transformation parameters
// [[Rcpp::export]]
List boundaryTransformationInit(List parameters) {
  
  std::vector<double> max_valuesV = as<std::vector<double> >(parameters["max"]);
  std::vector<double> min_valuesV = as<std::vector<double> >(parameters["min"]);
  
  int dimension = max_valuesV.size();
  D printf("Dimension: %d\n", dimension);
  
  double *max_values = (double*) malloc(dimension * sizeof(double));
  double *min_values = (double*) malloc(dimension * sizeof(double));
  
  D printf("Max values size: %lu\n", max_valuesV.size());
  D printf("Min values size: %lu\n", min_valuesV.size());
  
  max_values = vectorToArrayP(max_valuesV);
  min_values = vectorToArrayP(min_valuesV);
  
  cmaes_boundary_transformation_t boundaries;
  
  D printf("Start transformation\n");
  cmaes_boundary_transformation_init(&boundaries, min_values, max_values, dimension);
  D printf("Done\n");
  
  return boundaryTranfrormToList(boundaries);
}

//' Linear transformation of parameter vector to fit values to boundaries
//' @param boundaries List of boundary transformation parameters
//' @param pop Population for boundary transformation
//' @return Vector of population values after boundary transformation
// [[Rcpp::export]]
NumericVector boundaryTransformation(List boundaries, NumericVector pop) {
  cmaes_boundary_transformation_t boundary_struct = boundaryTransformationToStruct(boundaries);
  
  unsigned long pop_size = pop.size();
  double *params_in_bounds = (double *)malloc(pop_size * sizeof(double));
  double *values = vectorToArrayP(as<std::vector<double> >(pop)); 
  cmaes_boundary_transformation(&boundary_struct, values, params_in_bounds, boundary_struct.len_of_bounds);
  
  //Round the values according to their step values
  // ???
  //enforce_steps(params_in_bounds, steps, dimension);
  
  std::vector<double> result(params_in_bounds, params_in_bounds + pop_size);
  
  return Rcpp::wrap(result);
}

// [[Rcpp::export]]
List cmaesInit(NumericVector values, NumericVector stdDevs, long int inseed) {
  
  D printf("Parameters conversion\n");
  double *parameters = vectorToArrayP(as<std::vector<double> >(values));
  double *std_devs = vectorToArrayP(as<std::vector<double> >(stdDevs));
  D printf("Done P:%lu, SD:%lu\n", values.size(), stdDevs.size());
  
  
  double *arFunvals;
  cmaes_t evo;
  D printf("CMA-ES init\n");
  arFunvals = cmaes_init(&evo, values.size(), parameters, std_devs, inseed, 0, "no");
  D printf("Done\n");
  
  return cmaesToList(evo);
  //return List::create();
}

// [[Rcpp::export]]
List cmaesSamplePopulation(List cmaes) {
  D printf("-- Call: cmaesSamplePopulation\n");
  cmaes_t cmaes_struct = listToCmaes(cmaes);
  
  cmaes_WriteToFile(&cmaes_struct, "clock", "clock.dat");
  
  
  D printf("Done: listToCmaes\n");
  D printf("sp.N: %d\n", cmaes_struct.sp.N);
  D printf("sp.diagonalCov: %f\n", cmaes_struct.sp.diagonalCov);
  D printf("gen: %f\n", cmaes_struct.gen);
  
  /*
   
  double dpopsize = cmaes_Get(&cmaes_struct, "lambda");
  int popsize = (int) (dpopsize + 1e-5); //in case of a floating point rounding issue
  D printf("Popsize: %d\n", popsize);
  
  double *const *pop = cmaes_SamplePopulation(&cmaes_struct);
  
  
  double *popCopy = (double *)malloc(popsize * sizeof(double));
  memcpy(popCopy, pop, popsize * sizeof(double));
  //free(pop);
  
  std::vector<double> popVec(popCopy, popCopy + popsize);
  
  List result = List::create();
  
  result["cmaes"] = cmaesToList(cmaes_struct);
  result["pop"] = popVec;
  return result;
  */

  cmaes_SamplePopulation(&cmaes_struct);
  return cmaesToList(cmaes_struct);
}

// [[Rcpp::export]]
List cmaesUpdateDistribution(List cmaes, NumericVector rgFunValVec) {
  D printf("-- Call: cmaesUpdateDistribution\n");
  
  cmaes_t cmaes_struct = listToCmaes(cmaes);
  
  double *rgFunVal = vectorToArrayP(as<std::vector<double> >(rgFunValVec));
  
  D printf("Begin cmaes_UpdateDistribution\n");
  cmaes_UpdateDistribution(&cmaes_struct, rgFunVal);
  D printf("end cmaes_UpdateDistribution\n");
  
  return cmaesToList(cmaes_struct);
}

// [[Rcpp::export]]
List cmaesInitWithSamplePopulation(NumericVector values, NumericVector stdDevs) {
  
  D printf("Parameters conversion\n");
  double *parameters = vectorToArrayP(as<std::vector<double> >(values));
  double *std_devs = vectorToArrayP(as<std::vector<double> >(stdDevs));
  D printf("Done P:%lu, SD:%lu\n", values.size(), stdDevs.size());
  
  
  double *arFunvals;
  cmaes_t evo;
  D printf("CMA-ES init\n");
  arFunvals = cmaes_init(&evo, values.size(), parameters, std_devs, 0, 0, "no");
  D printf("Done\n");
  cmaes_SamplePopulation(&evo);
  
  return cmaesToList(evo);
  //return List::create();
}

double testFunc(double x) {
  return 2 * x * x + 8 * x + 3;
}

double testCostFunc(double a, double b, double c) {
  double sum = 0;
  for(double x = -10; x <= 10; x += 0.1) {
    double d = testFunc(x) - a*x*x - b*x - c;
    sum += d*d;
  }
  return sqrt(sum);
}

// [[Rcpp::export]]
NumericMatrix testConvergence(int times) {
  cmaes_t evo;
  double *arFunvals, *const*pop, *xfinal;
  cmaes_boundary_transformation_t boundaries;
  
  static int dim = 3;
  std::vector<double> max;
  max.push_back(10);
  max.push_back(10);
  max.push_back(10);
  
  std::vector<double> min;
  min.push_back(1);
  min.push_back(1);
  min.push_back(1);

  std::vector<double> std;
  for(int i = 0; i < dim; i++) {
    std.push_back((max[i] - min[i])/3);
  }
  std::vector<double> params; // some initial values in bounds
  params.push_back(3);
  params.push_back(3);
  params.push_back(3);
  
  double *min_values = vectorToArrayP(min);  
  double *max_values = vectorToArrayP(max);
  double *std_devs = vectorToArrayP(std);
  double *parameters = vectorToArrayP(params);
  
  D printf("Boundary transformation\n");
  cmaes_boundary_transformation_init(&boundaries, min_values, max_values, dim);
  
  //initialize cmaes
  D printf("Cmaes init\n");
  arFunvals = cmaes_init(&evo, dim, parameters, std_devs, 123, 0, "no");
  
  int runs = evo.sp.lambda;
  NumericMatrix result(times, runs); 
  for(int tm = 0; tm < times; tm++) {
    D printf("Iteration: %d\n", tm);
    pop = cmaes_SamplePopulation(&evo);

    std::vector<double> cost(runs);
    for (int rn = 0; rn < runs; rn++) {
      D printf("Run: %d\n", rn);

      //Adjust the parameters to be within in the bounds
      D printf("Boundary transformation\n");
      double* params_in_bounds = (double *)malloc(3 * sizeof(double));
      cmaes_boundary_transformation(&boundaries, pop[rn], params_in_bounds, dim);
      
      // TODO check enforce steps 
      //Round the values according to their step values
      //enforce_steps(params_in_bounds, steps, dimension);
      
      D printf("Params: %f %f %f \n", 
               round(params_in_bounds[0]), 
               round(params_in_bounds[1]), 
               round(params_in_bounds[2]));
      D printf("Cost function\n");
      cost[rn] = testCostFunc(round(params_in_bounds[0]), 
                                  round(params_in_bounds[1]), 
                                  round(params_in_bounds[2]));
      
      D printf("Cost: %f \n", cost[rn]);
      D printf("Saving result\n");
      result(tm, rn) = cost[rn];
      free(params_in_bounds);
    }
    arFunvals = vectorToArrayP(cost);
    cmaes_UpdateDistribution(&evo, arFunvals);
  }
  
  D printf("--- BEGIN Params ---\n");
  D printf("cmaes$sp$ccumcov: %f\n", (float)evo.sp.ccumcov);
  D printf("cmaes$sp$cs: %f\n", (float)evo.sp.cs);
  D printf("cmaes$gen: %f\n", (float)evo.gen);
  D printf("cmaes$chiN: %f\n", (float)evo.chiN);
  D printf("cmaes$sp$mu: %f\n", (float)evo.sp.mu);
  D printf("cmaes$sigma: %f\n", (float)evo.sigma);
  D printf("cmaes$sp$mucov: %f\n", (float)evo.sp.mucov);
  D printf("cmaes$sp$damps: %f\n", (float)evo.sp.damps);
  D printf("cmaes$sp$ccov: %f\n", (float)evo.sp.ccov);
  D printf("cmaes$mindiagC: %f\n", (float)evo.mindiagC);
  D printf("cmaes$maxdiagC: %f\n", (float)evo.maxdiagC);
  D printf("--- END Params ---\n");
  
  /*
  cmaes_WriteToFile(&evo, "all", "all.dat");
  cmaes_WriteToFile(&evo, "xmean", "xmean.dat");
  cmaes_WriteToFile(&evo, "C", "C.dat");
  cmaes_WriteToFile(&evo, "B", "B.dat");
  cmaes_WriteToFile(&evo, "diag(D)", "D.dat");
  */
  
  return result;
}

