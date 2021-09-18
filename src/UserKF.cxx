#include <math.h>
#include <stdio.h>
#include <fstream>
#include <string>
#include "TMath.h"
#include "TMinuit.h"
#include "TRandom3.h"
#include "TMatrixD.h"
#include "TGraph.h"
#include "TCanvas.h"
#include <iostream>
#include <random>
#include <chrono>

using namespace std;

class UserKF{

public:
  static double GetNpar(){ return fNpar; }

  static void SetPar(const double par[], const int npar=1){
    fLambda = par[0];

    if(npar==fNpar){
      fX   = par[1];
      fY   = par[2];
    }
  }

  static void SetPar(const TMinuit * mnt){
    double dummy;
    mnt->GetParameter(0, fLambda, dummy);
    mnt->GetParameter(1, fX, dummy);
    mnt->GetParameter(2, fY, dummy);
  }
  
  static double Constraint(){
    return fX*fX + fY*fY -1;
  }

  static double FullLikelihood(){
    return CoreLikelihood() + fLambda * Constraint();
  }
  
  static void IniCoreMIN(TMinuit * mnt){
    mnt->DefineParameter(0, "lambda", fLambda, 1e-2, -1e6, 1e6);
    mnt->FixParameter(0);  // Fix Lambda
    
    //user previous fit result as initial values
    mnt->DefineParameter(1, "x", fX, 1e-2, -1e6, 1e6);
    mnt->DefineParameter(2, "y", fY, 1e-2, -1e6, 1e6);
  }

  static bool IsConstraintGood(){
    const double eps = 1e-4;
    return (TMath::Abs(Constraint())<eps);
  }
  

  static void Print(const TString tag){
    printf("%20s lambda %10.6e X %10.6e Y %10.6e, core %20.6e constraint %20.6e full %20.6e\n", tag.Data(), fLambda, fX, fY, CoreLikelihood(), Constraint(), FullLikelihood());
  }

private:
  static double fLambda;
  static double fX;
  static double fY;
  static int fNpar;
  
  static double CoreLikelihood(){
    return fX + fY;
  }
};

double UserKF::fLambda = -999;
double UserKF::fX = -999;
double UserKF::fY = -999;
int UserKF::fNpar = 3;

void CoreFCN(int &npars, double *grad, double &value, double *par, int flag)
{
  UserKF::SetPar(par, UserKF::GetNpar());
  
  value = UserKF::FullLikelihood();
}

void LambdaFCN(int &npars, double *grad, double &value, double *par, int flag)
{
  //parameters and results are passed in this order
  //par -> UserKF -> CoreMIN -> UserKF -> value
  
  UserKF::SetPar(par, 1);//only setting lambda

  UserKF::Print("\nLambdaFCN before fit");

  // Second Minimization
  TMinuit * CoreMIN = new TMinuit(3);
  CoreMIN->SetPrintLevel(-1);
  
  CoreMIN->SetFCN(CoreFCN);

  UserKF::IniCoreMIN(CoreMIN);

  int flagL = CoreMIN->Command("MIGRAD");
  UserKF::Print("LambdaFCN after fit");
  
  int irun = 1;
  const int maxnrun = 2;//no need to try many times, fail alwasy if lambda is bad
  while(flagL!=0){
    printf("CoreMIN bad fit! %d ---- run once more [%d]\n", flagL, irun++); 

    flagL = CoreMIN->Command("MIGRAD");
    UserKF::Print("LambdaFCN after fit");
    
    if(irun>=maxnrun){
      break;
    }
  }
  
  if(flagL==0){
    value = TMath::Abs(UserKF::Constraint());
  }
  else{
    printf("CoreMIN giving up now... %d\n", flagL);
    value = 1E50;//must be large enough wrt possible FullLikelihood when fit fail
  }

  //already set during CoreMIN
  //UserKF::SetPar(CoreMIN);//only save fit values to MIN when converge

  delete CoreMIN;
}

bool DoubleMin(const double iniLambda, const double lmin, const double lmax)
{
  UserKF::Print("DoubleMin before fit");

  TMinuit * LambdaMIN = new TMinuit(1);
  LambdaMIN->SetPrintLevel(-1);

  LambdaMIN->SetFCN(LambdaFCN);

  LambdaMIN->DefineParameter(0, "lambda", iniLambda, 1e-2, lmin, lmax);
  
  int flag = LambdaMIN->Command("MIGRAD");
  UserKF::Print("DoubleMin after fit");
  
  int irun = 1;
  const int maxnrun=20;//if time permits, the larger the better  
  while(flag!=0 || ! UserKF::IsConstraintGood()){
    printf("LambdaMIN bad fit! %d %e ------- run once more! [%d]\n", flag, UserKF::Constraint(), irun++);
    
    flag = LambdaMIN->Command("MIGRAD");
    UserKF::Print("DoubleMin after fit");
    
    if(irun>=maxnrun){
      break;
    }
  }

  delete LambdaMIN;

  if(flag==0 && UserKF::IsConstraintGood()){
    printf("DoubleMin finishes: it works for %f %f %f\n", iniLambda, lmin, lmax);
    return true;
  }
  else{
    printf("DoubleMin finishes: giving up now... %d %f for  %f %f %f\n", flag, UserKF::Constraint(), iniLambda, lmin, lmax);
    return false;
  }
}


int main()
{
  DoubleMin(0.5, 0, 1000);
  DoubleMin(  5, 0, 100);
  DoubleMin( 10, 0, 100);
  DoubleMin( 10, 0, 1000);//fail -> works now after changing the range and step of X and Y to general ones
  DoubleMin(  5, -1, 100);
  DoubleMin(5, -10, 100);//fail at maxnrun =4, but works with maxnrun >= 9
  DoubleMin(5, 0, 100);
  DoubleMin(5, -100, 100);//works 3 tries
  DoubleMin(10, -1e6, 1e6);//don't start with 0, no sensitivity, fail eventually: lambda can't be too large
  DoubleMin(10, -1e3, 1e3);
  DoubleMin(0, -1e3, 1e3);//fail: lambda can't start as 0

  return 0;
}

