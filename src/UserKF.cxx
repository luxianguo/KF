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
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TLegend.h"
#include "TLegendEntry.h"

#include <iostream>
#include <random>
#include <chrono>

using namespace std;


class UserKF{

public:
  static void SetPar(const double par[], const int npar){
    if(npar!=fNpar-2){//npar inside TMinuit is only for free par
      printf("npar mismatch! %d %d\n", npar, fNpar); exit(1);
    }
    
    fLambda = par[0];
    fX      = par[1];
    fY      = par[2];
    fZ      = E3;//par[3];

  }

  static double Constraint(){
    //return fX + fY + fZ - 29;
    return 2*fX*fY*(1-cos(fZ)) - 0.134977*0.134977;
  }

  static double FullLikelihood(){
    return CoreLikelihood(E1,E2,E3) + fLambda * Constraint();
  }
  
  static void IniCoreMIN(TMinuit * mnt, const double inputl){
    mnt->DefineParameter(0, "lambda", inputl, 1e-2, -1e6, 1e6);//overwrite lambda from input
    mnt->FixParameter(0);  // Fix Lambda
    
    //user previous fit result as initial values
    mnt->DefineParameter(1, "x", fX, 1e-2, -1e6, 1e6);
    mnt->DefineParameter(2, "y", fY, 1e-2, -1e6, 1e6);
    mnt->DefineParameter(3, "z", fZ, 1e-2, -1e6, 1e6);
    mnt->FixParameter(3);
  }

  static bool IsConstraintGood(){
    const double eps = 1e-4;
    return (TMath::Abs(Constraint())<eps);
  }

  static TMatrixD CovarianceMatrix(){
    // Defines the covariance matrix and the variables
    double V[9];
    // Set the covariance matrix elements 
    for(int i=0;i<3;i++) {
      for(int j=0;j<3;j++) {
        V[3*i+j]= 10E-5; //off diagonal elements have no co-correlation here but a small number here does no effect the algorithm
        if(i==j && (i==0)) {V[3*i+j]=TMath::Power(0.5,2.);}
        if(i==j && (i==1)) {V[3*i+j]=TMath::Power(1,2.);}
        if(i==j && (i==2)) {V[3*i+j]=TMath::Power(1,2.);}
      }
    }
    // Create the matrix
    int nparameters = fNpar-1;
    TMatrixD CovMatrix(nparameters,nparameters);
    for(int i=0;i<nparameters*nparameters;i++) {
      int x=floor(i/nparameters);
      int y=i%nparameters;
      CovMatrix[x][y]=V[i];
    }
    return CovMatrix;
  }

  static void Print(const TString tag){
    printf("%20s lambda %10.6e X %10.6e Y %10.6e, core %20.6e constraint %20.6e full %20.6e\n", tag.Data(), fLambda, fX, fY, CoreLikelihood(E1,E2,E3), Constraint(), FullLikelihood());
  }
  static double GetfX(){
    return fX;
  }
  static double GetfY(){
    return fY;
  }
  static double GetfZ(){
    return fZ;
  }
  static double GetE1(){
    return E1;
  }
  static double GetE2(){
    return E2;
  }
  static double GetE3(){
    return E3;
  }
  static double GetfLambda(){
    return fLambda;
  }
  static void Set(double e1, double e2, double e3){
    E1 = e1;
    E2 = e2;
    E3 = e3;
  }
  static void SetCVM(vector<double> V){
    int size = V.size();
    int dim = sqrt(size);
    //CovMatrix.ResizeTo(dim,dim);
    for(int i=0;i<size;i++) {
      int x=floor(i/dim);
      int y=i%dim;
      CovMatrix[x][y]=V[i];
    }
    CovMatrix[0][2] = 0;
    CovMatrix[2][0] = 0;
    CovMatrix[2][2] = 1;
  }
private:
  static double fLambda;
  static double fX;
  static double fY;
  static double fZ;
  static const int fNpar;

  static double E1;
  static double E2;
  static double E3;

  static TMatrixD CovMatrix;

  
  static double CoreLikelihood(double E1, double E2, double E3){

    TMatrixD CovMatrix = CovarianceMatrix();
    TMatrixD CovMatrixInverse = CovMatrix.Invert();

    int nparameters = fNpar-1;

    TMatrixD Diff(nparameters,1);
    Diff[0][0]=fX-E1;
    Diff[1][0]=fY-E2;
    Diff[2][0]=fZ-E3;

    TMatrixD DiffT(1,nparameters);
    DiffT.Transpose(Diff);

    TMatrixD Chi2 = (DiffT*CovMatrixInverse*Diff);
    
    return Chi2[0][0];

  }
};

double UserKF::fLambda = -999;
double UserKF::fX = -999;
double UserKF::fY = -999;
double UserKF::fZ = -999;
const int UserKF::fNpar = 4;

double UserKF::E1 = -999;
double UserKF::E2 = -999;
double UserKF::E3 = -999;
TMatrixD tmpMatrix(3,3);
TMatrixD UserKF::CovMatrix = tmpMatrix;


vector<double> IniE(){

  std::default_random_engine engine; 
  engine.seed(std::chrono::system_clock::now().time_since_epoch().count());

  // Gaussian mean followed by stdiv
  std::normal_distribution<double> nd1(4*0.8, 0.5); 
  std::normal_distribution<double> nd2(9*0.8, 1); 
  std::normal_distribution<double> nd3(16*0.8, 1); 
  // Generate the intial E's value
  double E1 = nd1(engine);
  double E2 = nd2(engine);
  double E3 = nd3(engine);
  vector<double> vet;
  vet.push_back(E1);
  vet.push_back(E2);
  vet.push_back(E3);
  UserKF::Set(E1,E2,E3);
  return vet;

}

void CoreFCN(int &npars, double *grad, double &value, double *par, int flag)
{
  UserKF::SetPar(par, npars);
  
  value = UserKF::FullLikelihood();
}

void LambdaFCN(int &npars, double *grad, double &value, double *par, int flag)
{
  //parameters and results are passed in this order
  //par/UserKF -> CoreMIN -> UserKF -> value
  
  //UserKF::Print("\nLambdaFCN before fit");

  // Second Minimization
  TMinuit * CoreMIN = new TMinuit(3);
  CoreMIN->SetPrintLevel(-1);
  
  CoreMIN->SetFCN(CoreFCN);

  UserKF::IniCoreMIN(CoreMIN, par[0]);

  int flagL = CoreMIN->Command("MIGRAD");
  //UserKF::Print("LambdaFCN after fit");
  
  int irun = 1;
  const int maxnrun = 2;//no need to try many times, fail alwasy if lambda is bad
  while(flagL!=0){
    printf("CoreMIN bad fit! %d ---- run once more [%d]\n", flagL, irun++); 

    flagL = CoreMIN->Command("MIGRAD");
    //UserKF::Print("LambdaFCN after fit");
    
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

void SetInitials(double iniE1, double iniE2, double iniOA, vector<double> CVM){
  UserKF::Set(iniE1,iniE2,iniOA);
  UserKF::SetCVM(CVM);
}

bool DoubleMin(const double iniLambda, const double lmin, const double lmax, double iniE1, double iniE2, double iniOA, vector<double> CVM, double &optE1, double &optE2, int &fg)
{
  //UserKF::Print("DoubleMin before fit\n");
  SetInitials(iniE1, iniE2, iniOA, CVM);

  TMinuit * LambdaMIN = new TMinuit(1);
  LambdaMIN->SetPrintLevel(-1);

  LambdaMIN->SetFCN(LambdaFCN);

  LambdaMIN->DefineParameter(0, "lambda", iniLambda, 1e-2, lmin, lmax);
  
  int flag = LambdaMIN->Command("MIGRAD");
  //UserKF::Print("DoubleMin after fit\n");
  
  int irun = 1;
  const int maxnrun=20;//if time permits, the larger the better  
  while(flag!=0 || ! UserKF::IsConstraintGood()){
    printf("LambdaMIN bad fit! %d %e ------- run once more! [%d]\n", flag, UserKF::Constraint(), irun++);
    
    flag = LambdaMIN->Command("MIGRAD");
    //UserKF::Print("DoubleMin after fit\n");
    
    if(irun>=maxnrun){
      break;
    }
  }

  delete LambdaMIN;

  if(flag==0 && UserKF::IsConstraintGood()){
    printf("DoubleMin finishes: it works for %f %f %f\n", iniLambda, lmin, lmax);
    cout << "E1: " <<  UserKF::GetE1() << " E2: " <<  UserKF::GetE2() << " E3: " <<  UserKF::GetE3()<< endl;
    cout << "X: " <<  UserKF::GetfX() << " Y: " <<  UserKF::GetfY() << " Z: " <<  UserKF::GetfZ()<< " lambda: " <<   UserKF::GetfLambda() << endl;
    optE1 = UserKF::GetfX();
    optE2 = UserKF::GetfY();
    double massPost = sqrt(2*optE1*optE2*(1-cos(UserKF::GetfZ())));
    cout << "massPost: " << massPost << endl;
    if(massPost < 0.1) cout << "small mass!!!" << endl;
    fg = 0;
    return true;
  }
  else{
    printf("DoubleMin finishes: giving up now... %d %f for  %f %f %f\n", flag, UserKF::Constraint(), iniLambda, lmin, lmax);
    fg = 4;
    return false;
  }
}



/*
int main()
{
  
  TH1F *hBefore = new TH1F("hBefore","Energy - Before Fitting",80,0,20);
  TH1F *hAfter = new TH1F("hAfter","Energy - After Fitting",80,0,20);

  for(int i = 0; i < 1000; i++){
    vector<double> vet = IniE();

    for(unsigned int i = 0; i < vet.size(); i++){
      cout << "vet: " << vet[i] << endl;
    }
    hBefore->Fill(vet[0]);
    hBefore->Fill(vet[1]);
    hBefore->Fill(vet[2]);

    DoubleMin(5, -10, 100);

    cout << "xxX: " <<  UserKF::GetfX() << " Y: " <<  UserKF::GetfY() << " Z: " <<  UserKF::GetfZ()<< " lambda: " <<   UserKF::GetfLambda() << endl;
    double E1 = UserKF::GetfX();
    double E2 = UserKF::GetfY();
    double E3 = UserKF::GetfZ();
    hAfter->Fill(E1);
    hAfter->Fill(E2);
    hAfter->Fill(E3);
  }

  TCanvas * c1 = new TCanvas("c1", "", 1200, 800);
  auto legend = new TLegend(0.5,0.7,0.68,0.88);
  hBefore->SetMaximum(330);
  //hBefore->SetStats(0);
  hBefore->SetFillStyle(4050);
  hBefore->SetFillColor(24);
  hBefore->SetLineColor(24);
  hBefore->Draw("hist");
  hAfter->SetFillStyle(3001);
  hAfter->SetFillColor(46);
  hAfter->SetLineColor(46);
  hAfter->Draw("SAMES hist");
  legend->AddEntry(hBefore,"Before Fitting","f");
  legend->AddEntry(hAfter,"After Fitting","f");
  legend->Draw("same");
  c1->Print("hEnergyFitting.png");

  Other tests
  DoubleMin(0.5, 0, 1000);
  DoubleMin(  5, 0, 100);
  DoubleMin( 10, 0, 100);
  DoubleMin( 10, 0, 1000);//fail -> works now after changing the range and step of X and Y to general ones
  DoubleMin(  5, -1, 100);
  DoubleMin(5, -10, 100);//fail at maxnrun =4, but works with maxnrun >= 9
  DoubleMin(5, 0, 100);
  DoubleMin(5, -100, 100);//works 3 tries
  DoubleMin(10, -1e6, 1e6);//fail eventually: lambda can't be too large
  DoubleMin(10, -1e3, 1e3);
  DoubleMin(0, -1e3, 1e3);//fail: lambda can't start as 0
  DoubleMin(1e-2, -1e3, 1e3);//
  */

  
  //return 0;
//}

