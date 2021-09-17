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
  static double GetLambda(){ return fLambda; }
  static double GetX(){ return fX; }
  static double GetY(){ return fY; }
  static double GetNpar(){ return fNpar; }

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
  
  static void SetPar(const double par[], const int npar=1){
    fLambda = par[0];

    if(npar==fNpar){
      fX   = par[1];
      fY   = par[2];
    }
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

  //printf("CoreFCN %f %f %f %f\n", gXfit, gYfit, lambda, value);
}


void LambdaFCN(int &npars, double *grad, double &value, double *par, int flag)
{
  UserKF::SetPar(par, 1);//only setting lambda

  UserKF::Print("\nLambdaFCN before fit");

  // Second Minimization
  TMinuit * CoreMIN = new TMinuit(3);
  CoreMIN->SetPrintLevel(-1);
  
  CoreMIN->SetFCN(CoreFCN);

  CoreMIN->DefineParameter(0, "lambda", UserKF::GetLambda(), 0.01, -1000, 1000);
  //user previous fit result as initial values
  CoreMIN->DefineParameter(1, "gXfit", UserKF::GetX(), 0.1, -2000, 3000);
  CoreMIN->DefineParameter(2, "gYfit", UserKF::GetY(), 0.1, -2000, 3000);
  CoreMIN->FixParameter(0);  // Fix Lambda

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

void DoubleMin()
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
  //LambdaMIN->DefineParameter(0, "lambda",   5, 0.01, -10, 100);//fail
  LambdaMIN->DefineParameter(0, "lambda",   5, 0.01, 0, 100);//works
  
  int flag = LambdaMIN->Command("MIGRAD");
  double finalConstraint = UserKF::Constraint();
  
  const double EPS = 1e-4;
  if(flag!=0 || TMath::Abs(finalConstraint) > EPS){
    printf("LambdaMIN bad fit! %d %e ------- run once more!\n", flag, finalConstraint);
    
    flag = LambdaMIN->Command("MIGRAD");
    finalConstraint = UserKF::Constraint();
    
    if(flag!=0 || TMath::Abs(finalConstraint) > EPS){
      printf("LambdaMIN bad fit! %d %e ------- run once more!\n", flag, finalConstraint);
      
      flag = LambdaMIN->Command("MIGRAD");
      finalConstraint = UserKF::Constraint();

      if(flag!=0 || TMath::Abs(finalConstraint) > EPS){
        printf("final fit not converge! %d %e exit!\n", flag, finalConstraint); exit(1);
      }
    }
  }

  delete LambdaMIN;

  UserKF::Print("DoubleMin after fit");
}


int main()
{
  DoubleMin();
}

