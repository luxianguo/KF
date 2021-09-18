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
    mnt->DefineParameter(0, "lambda", fLambda, 0.01, -1000, 1000);
    mnt->FixParameter(0);  // Fix Lambda
    
    //user previous fit result as initial values
    mnt->DefineParameter(1, "x", fX, 0.1, -2000, 3000);
    mnt->DefineParameter(2, "y", fY, 0.1, -2000, 3000);
  }

  static bool IsConstraintGood(){
    const double eps = 1e-4;
    return (TMath::Abs(Constraint())<eps);
  }
  

  static void Print(const TString tag){
    printf("%20s lambda %10.6f X %10.6f Y %10.6f, core %20.6f constraint %20.6f full %20.6f\n", tag.Data(), fLambda, fX, fY, CoreLikelihood(), Constraint(), FullLikelihood());
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
  if(flagL!=0){
    printf("CoreMIN bad fit! %d ---- run once mroe\n", flagL); //exit(1);
    flagL = CoreMIN->Command("MIGRAD");
  }
  
  UserKF::SetPar(CoreMIN);

  delete CoreMIN;

  UserKF::Print("LambdaFCN after fit");

  if(flagL!=0){
    printf("CoreMIN bad fit! %d\n", flagL); //exit(1);
    value = 1E20;
  }
  else{
    value = TMath::Abs(UserKF::Constraint());
  }
}

bool DoubleMin()
{
  UserKF::Print("DoubleMin before fit");

  TMinuit * LambdaMIN = new TMinuit(1);
  LambdaMIN->SetPrintLevel(-1);

  LambdaMIN->SetFCN(LambdaFCN);

  //LambdaMIN->DefineParameter(0, "lambda", 0.5, 0.01, 0, 1000);
  //LambdaMIN->DefineParameter(0, "lambda",   5, 0.01, 0, 100);//works
  //LambdaMIN->DefineParameter(0, "lambda", 0.5, 0.01, 0, 1000);//works
  //LambdaMIN->DefineParameter(0, "lambda",   10, 0.01, 0, 100);//works
  //LambdaMIN->DefineParameter(0, "lambda",   10, 0.01, 0, 1000);//fail
  //LambdaMIN->DefineParameter(0, "lambda",   5, 0.01, -1, 100);//works
  //LambdaMIN->DefineParameter(0, "lambda",   5, 0.01, -10, 100);//fail at maxnrun =4, but works with maxnrun >= 9
  //LambdaMIN->DefineParameter(0, "lambda",   5, 0.01, 0, 100);//works
  LambdaMIN->DefineParameter(0, "lambda",   5, 0.01, -100, 100);//
  
  int flag = LambdaMIN->Command("MIGRAD");
  int irun = 1;
  const int maxnrun=10;//if time permits, the larger the better
  
  while(flag!=0 || ! UserKF::IsConstraintGood()){
    printf("LambdaMIN bad fit! %d %e ------- run once more! [%d]\n", flag, UserKF::Constraint(), irun++);
    
    flag = LambdaMIN->Command("MIGRAD");

    if(irun>=maxnrun){
      printf("giving up now... %d %f\n", flag, UserKF::Constraint());
      break;
    }
  }

  delete LambdaMIN;

  UserKF::Print("DoubleMin after fit");

  return (flag==0 && UserKF::IsConstraintGood());
}


int main()
{
  DoubleMin();
}

