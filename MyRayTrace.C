#include "TSystem.h"
#include <iostream>
#include <fstream>
#include <string>
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TString.h"
#include "TNtuple.h"
#include "TAxis.h"
#include "TStopwatch.h"
#include "TStyle.h"

const Double_t pi=4.0*atan(1.0); /**< Gives back value of Pi */
const Double_t spedc=299792458.0; /**< Speed of Light in m/s */

const Double_t A=1.78;/**< Value of Parameter A for SP model of refractive index */
//AraSim model
const Double_t B=-0.43; /**< Value of Parameter B for SP model of refractive index */
const Double_t C=0.0132;/**< Value of Parameter A for SP model of refractive index. There will be a negative sign for postive depth depth. */ 

//ARIANNA model
// const Double_t B=-0.427; /**< Value of Parameter B for SP model of refractive index */
// const Double_t C=1.0/71.0;/**< Value of Parameter A for SP model of refractive index. There will be a negative sign for postive depth depth. */ 

//More temperature and Lattenuation models in TF1s fFr1 and fFr2:
//Temperature model:Takes in value of depth z and returns the value of temperature.
//The temperature model has been taken from AraSim which also took it from here http://icecube.wisc.edu/~mnewcomb/radio/atten/ . This is basically Matt Newcomb's icecube directory which has alot of information, plots and codes about South Pole Ice activities. Please read it if you find it interesting.
// http://icecube.wisc.edu/~mnewcomb/radio/#iceabsorbtion

//Lattenuation model: Takes in value of frequency in Ghz and depth z and returns you the value of attenuation length.

Double_t *MyRayTrace(){
//Double_t *MyRayTrace(Double_t x0, Double_t z0, Double_t x1, Double_t z1){
  
  TStopwatch watch;
  
  Double_t *output=new Double_t[9];

  //Plot the ray solutions
  Bool_t Plot=true;
  //calculate the attenuation
  Bool_t attcal=false;

  //My solution can not launch rays below 90 degrees for that I have to switch/flip Tx and Rx depths, trace the rays and then put the Tx and Rx back in their original positions
  Bool_t Flip=false;

  Double_t x0=0;
  Double_t z0=-100;
  Double_t x1=100;
  Double_t z1=-5;

  // Double_t x0=0;
  // Double_t z0=-67.5489;
  // Double_t x1=1338.3;
  // Double_t z1=-800;

  Double_t D56cor[3]={x0,0,z0};//Tx positions
  Double_t antcor[3]={x1,0,z1};//Rx Positions

  Double_t lowerz=0;
  Double_t upz=0;
  if(antcor[2]<D56cor[2]){
    lowerz=antcor[2];
    upz=D56cor[2];
  }
  if(antcor[2]>D56cor[2]){
    lowerz=D56cor[2];
    upz=antcor[2];
  }
  if(antcor[2]==D56cor[2]){
    lowerz=D56cor[2];
    upz=antcor[2];
  }
  lowerz=lowerz;

  Double_t dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }
  
  //Ray path function x(z)
  //for reflected ray
  TF1 * fR=new TF1("fR","([3]/[2])*(1.0/sqrt([0]*[0]-[3]*[3]))*([2]*x-log([0]*([0]+[1]*exp([2]*x))-[3]*[3]+sqrt([0]*[0]-[3]*[3])*sqrt(pow([0]+[1]*exp([2]*x),2)-[3]*[3])))",1*pow(10,-9),20000);
  fR->FixParameter(0,A);
  fR->FixParameter(1,B);
  fR->FixParameter(2,-C);
  fR->SetParameter(3,1.5);

  //For direct ray
  TF1 * fD=new TF1("fD","([3]/[2])*(1.0/sqrt([0]*[0]-[3]*[3]))*([2]*x-log([0]*([0]+[1]*exp([2]*x))-[3]*[3]+sqrt([0]*[0]-[3]*[3])*sqrt(pow([0]+[1]*exp([2]*x),2)-[3]*[3])))",-200000, -1*pow(10,-9));
  fD->FixParameter(0,A);
  fD->FixParameter(1,B);
  fD->FixParameter(2,C);
  fD->SetParameter(3,1.5);

  ////Calculate launch angle for Direct ray
  TF1 * fDa=new TF1("fDa","(x/[2])*(1.0/sqrt([0]*[0]-x*x))*([2]*[3]-log([0]*([0]+[1]*exp([2]*[3]))-x*x+sqrt([0]*[0]-x*x)*sqrt(pow([0]+[1]*exp([2]*[3]),2)-x*x))) - (x/[2])*(1.0/sqrt([0]*[0]-x*x))*([2]*[4]-log([0]*([0]+[1]*exp([2]*[4]))-x*x+sqrt([0]*[0]-x*x)*sqrt(pow([0]+[1]*exp([2]*[4]),2)-x*x))) - [5]",1*pow(10,-9),(A+B*exp(-(C)*(-z0)))*sin(90*(pi/180.0)));
  fDa->FixParameter(0,A);
  fDa->FixParameter(1,B);
  fDa->FixParameter(2,C);
  fDa->FixParameter(3,z1);
  fDa->FixParameter(4,z0);
  fDa->FixParameter(5,x1);
  Double_t langD=asin((fDa->GetX(-0.0))/(A+B*exp(-(C)*(-z0))))*(180.0/pi);
  Double_t checkzeroD=fDa->Eval(fDa->GetX(-0.0));  
  Double_t lvalueD=fDa->GetX(-0.0);

  ////Calculate travel time for Direct ray
  TF1 * ftimeD=new TF1("ftimeD","(1.0/([3]*[2]*sqrt(pow([0]+[1]*exp(x*[2]),2)-[4]*[4])))*(pow([0]+[1]*exp(x*[2]),2)-[4]*[4]+([2]*x-log([0]*([0]+[1]*exp([2]*x))-[4]*[4]+sqrt([0]*[0]-[4]*[4])*sqrt(pow([0]+[1]*exp([2]*x),2)-[4]*[4])))*([0]*[0]*sqrt(pow([0]+[1]*exp(x*[2]),2)-[4]*[4]))/sqrt([0]*[0]-[4]*[4]) +[0]*sqrt(pow([0]+[1]*exp(x*[2]),2)-[4]*[4])*log([0]+[1]*exp(x*[2])+sqrt(pow([0]+[1]*exp(x*[2]),2)-[4]*[4])) )",1*pow(10,-9),200000);
  ftimeD->FixParameter(0,A);
  ftimeD->FixParameter(1,B);
  ftimeD->FixParameter(2,-C);
  ftimeD->FixParameter(3,spedc);
  ftimeD->SetParameter(4,lvalueD);

  //(1.0/(c*C*sqrt(pow(A+B*exp(x*C),2)-[4]*[4])))*(pow(A+B*exp(x*C),2)-[4]*[4]+(C*x-log(A*(A+B*exp(C*x))-[4]*[4]+sqrt(A*A-[4]*[4])*sqrt(pow(A+B*exp(C*x),2)-[4]*[4])))*(A*A*sqrt(pow(A+B*exp(x*C),2)-[4]*[4]))/sqrt(A*A-[4]*[4]) +A*sqrt(pow(A+B*exp(x*C),2)-[4]*[4])*log(A+B*exp(x*C)+sqrt(pow(A+B*exp(x*C),2)-[4]*[4])) )
  
  Double_t timeD=ftimeD->Eval(-z0)-ftimeD->Eval(-z1);
  //cout<<"The D time is: "<<timeD<<endl;

  ////Calculate launch angle for Reflected ray
  TF1 * fRa=new TF1("fRa","(x/[2])*(1.0/sqrt([0]*[0]-x*x))*([2]*[3]-log([0]*([0]+[1]*exp([2]*[3]))-x*x+sqrt([0]*[0]-x*x)*sqrt(pow([0]+[1]*exp([2]*[3]),2)-x*x))) - (x/[2])*(1.0/sqrt([0]*[0]-x*x))*([2]*[4]-log([0]*([0]+[1]*exp([2]*[4]))-x*x+sqrt([0]*[0]-x*x)*sqrt(pow([0]+[1]*exp([2]*[4]),2)-x*x))) - 2*( (x/[2])*(1.0/sqrt([0]*[0]-x*x))*([2]*0.0 -log([0]*([0]+[1]*exp([2]*0.0 ))-x*x+sqrt([0]*[0]-x*x)*sqrt(pow([0]+[1]*exp([2]*0.0 ),2)-x*x)))  - (x/[2])*(1.0/sqrt([0]*[0]-x*x))*([2]*[4]-log([0]*([0]+[1]*exp([2]*[4]))-x*x+sqrt([0]*[0]-x*x)*sqrt(pow([0]+[1]*exp([2]*[4]),2)-x*x))) ) - [5]",1*pow(10,-9),(A+B*exp(-(C)*(-z0)))*sin(90*(pi/180.0)));
  fRa->FixParameter(0,A);
  fRa->FixParameter(1,B);
  fRa->FixParameter(2,-(C));
  fRa->FixParameter(3,-z1);
  fRa->FixParameter(4,-z0);
  fRa->FixParameter(5,x1);
  Double_t langR=asin((fRa->GetX(0.0))/(A+B*exp(-(C)*(-z0))))*(180.0/pi);
  Double_t checkzeroR=fRa->Eval(fRa->GetX(0.0));
  Double_t lvalueR=fRa->GetX(0.0);

  ////Calculate travel time for Reflected ray
  TF1 * ftimeR=new TF1("ftimeR","(1.0/([3]*[2]*sqrt(pow([0]+[1]*exp(x*[2]),2)-[4]*[4])))*(pow([0]+[1]*exp(x*[2]),2)-[4]*[4]+([2]*x-log([0]*([0]+[1]*exp([2]*x))-[4]*[4]+sqrt([0]*[0]-[4]*[4])*sqrt(pow([0]+[1]*exp([2]*x),2)-[4]*[4])))*([0]*[0]*sqrt(pow([0]+[1]*exp(x*[2]),2)-[4]*[4]))/sqrt([0]*[0]-[4]*[4]) +[0]*sqrt(pow([0]+[1]*exp(x*[2]),2)-[4]*[4])*log([0]+[1]*exp(x*[2])+sqrt(pow([0]+[1]*exp(x*[2]),2)-[4]*[4])) )",-200000, -1*pow(10,-9));
  ftimeR->FixParameter(0,A);
  ftimeR->FixParameter(1,B);
  ftimeR->FixParameter(2,C);
  ftimeR->FixParameter(3,spedc);
  ftimeR->SetParameter(4,lvalueR);
  
  Double_t timeR=ftimeR->Eval(-0.000001)-ftimeR->Eval(z0)+ftimeR->Eval(-0.000001)-ftimeR->Eval(z1);
  //cout<<"The R time is: "<<timeR<<endl;

  ////Calculate launch angle for Refracted ray
  TF1 * fRaa=new TF1("fRaa","(x/[2])*(1.0/sqrt([0]*[0]-x*x))*([2]*[3]-log([0]*([0]+[1]*exp([2]*[3]))-x*x+sqrt([0]*[0]-x*x)*sqrt(pow([0]+[1]*exp([2]*[3]),2)-x*x))) - (x/[2])*(1.0/sqrt([0]*[0]-x*x))*([2]*[4]-log([0]*([0]+[1]*exp([2]*[4]))-x*x+sqrt([0]*[0]-x*x)*sqrt(pow([0]+[1]*exp([2]*[4]),2)-x*x))) - 2*( (x/[2])*(1.0/sqrt([0]*[0]-x*x))*([2]*((1.0/[2])*log(fabs((x-[0])/[1]))) -log([0]*([0]+[1]*exp([2]*((1.0/[2])*log(fabs((x-[0])/[1]))) ))-x*x+sqrt([0]*[0]-x*x)*sqrt(pow([0]+[1]*exp([2]*((1.0/[2])*log(fabs((x-[0])/[1]))) ),2)-x*x)))  - (x/[2])*(1.0/sqrt([0]*[0]-x*x))*([2]*[4]-log([0]*([0]+[1]*exp([2]*[4]))-x*x+sqrt([0]*[0]-x*x)*sqrt(pow([0]+[1]*exp([2]*[4]),2)-x*x))) ) - [5]",1*pow(10,-9),(A+B*exp(-(C)*(-z0)))*sin(90*(pi/180.0)));
  fRaa->FixParameter(0,A);
  fRaa->FixParameter(1,B);
  fRaa->FixParameter(2,-(C));
  fRaa->FixParameter(3,-z1);
  fRaa->FixParameter(4,-z0);
  fRaa->FixParameter(5,x1);

  Double_t langRa=asin((fRaa->GetX(0.0))/(A+B*exp(-(C)*(-z0))))*(180.0/pi);
  Double_t checkzeroRa=fRaa->Eval(fRaa->GetX(0.0));
  Double_t lvalueRa=fRaa->GetX(0.0);
  //Double_t angrange=fabs(langR-((pi/2)-atan(fabs((z1-z0)/x1)))*(180.0/pi));
  Double_t angrange=fabs(langR-90);
  Double_t angstep=angrange/20;

  Double_t tempval[2];
  Int_t itemp=0;
  for(Int_t iang=0;iang<21;iang++){
    fRaa->SetRange((A+B*exp(-(C)*(-z0)))*sin(((langR+iang*angstep)*(pi/180.0))),(A+B*exp(-(C)*(-z0)))*sin((90)*(pi/180.0)));
    //cout<<((langR+iang*angstep)*(pi/180.0))*(180.0/pi)<<endl;
    lvalueRa=fRaa->GetX(0.0);
    langRa=langR+iang*angstep;
    checkzeroRa=fRaa->Eval(lvalueRa);;
    //cout<<"the derivative is "<<fR->Derivative(fabs((1.0/C)*log((lvalueRa-A)/B))+0.001,0,1*pow(10,-9))<<" "<<fD->Derivative(-fabs((1.0/C)*log((lvalueRa-A)/B))-0.001,0,1*pow(10,-9))<<" "<<checkzeroRa<<" "<<(1.0/(-C))*log((lvalueRa-A)/B)<<" "<<asin((fRaa->GetX(0.0))/(A+B*exp(-(C)*(-z0))))*(180.0/pi)<<" "<<lvalueRa<<endl;
    tempval[itemp]=lvalueRa;
    langRa=asin((fRaa->GetX(0.0))/(A+B*exp(-(C)*(-z0))))*(180.0/pi);
  
    itemp++;
    if(itemp>1 && itemp<4){
      itemp=0;
    }
    if(iang>0 && fabs(tempval[0]-tempval[1])<1*pow(10,-7) && (1.0/(-C))*log((lvalueRa-A)/B)>0){
      iang=20;
      itemp=20;
    }
 
  }

  if(TMath::IsNaN(checkzeroRa)==true){
    checkzeroRa=1000;
  }
  //cout<<"lvalueRa is "<<lvalueRa<<" "<<itemp<<endl; 

  ////Calculate travel time for Refracted ray
  TF1 * ftimeRa=new TF1("ftimeRa","(1.0/([3]*[2]*sqrt(pow([0]+[1]*exp(x*[2]),2)-[4]*[4])))*(pow([0]+[1]*exp(x*[2]),2)-[4]*[4]+([2]*x-log([0]*([0]+[1]*exp([2]*x))-[4]*[4]+sqrt([0]*[0]-[4]*[4])*sqrt(pow([0]+[1]*exp([2]*x),2)-[4]*[4])))*([0]*[0]*sqrt(pow([0]+[1]*exp(x*[2]),2)-[4]*[4]))/sqrt([0]*[0]-[4]*[4]) +[0]*sqrt(pow([0]+[1]*exp(x*[2]),2)-[4]*[4])*log([0]+[1]*exp(x*[2])+sqrt(pow([0]+[1]*exp(x*[2]),2)-[4]*[4])) )",-200000, -1*pow(10,-9));
  ftimeRa->FixParameter(0,A);
  ftimeRa->FixParameter(1,B);
  ftimeRa->FixParameter(2,C);
  ftimeRa->FixParameter(3,spedc);
  ftimeRa->SetParameter(4,lvalueRa);
  
  Double_t timeRa=0;
  Double_t raytime=0;
  if((z0<-fabs((-1.0/C)*log((lvalueRa-A)/B)) || fabs((-1.0/C)*log((lvalueRa-A)/B))<-z1) && itemp==20){
    raytime=ftimeRa->Eval(-fabs((-1.0/C)*log((lvalueRa-A)/B))-0.00001)-ftimeRa->Eval(z0)+ftimeRa->Eval(-fabs((-1.0/C)*log((lvalueRa-A)/B))-0.00001)-ftimeRa->Eval(z1);
    //cout<<"The Ra time is: "<<raytime<<endl;
  }
  timeRa=raytime;
   
  ////Callculate the receive angles
  fD->SetParameter(3,lvalueD);
  Double_t RangD=atan(fD->Derivative(z1,0,1*pow(10,-9)))*(180.0/pi);
  fR->SetParameter(3,lvalueR);
  Double_t RangR=180-atan(fR->Derivative(-z1,0,1*pow(10,-9)))*(180.0/pi);
  fR->SetParameter(3,lvalueRa);
  Double_t RangRa=180-atan(fR->Derivative(-z1,0,1*pow(10,-9)))*(180.0/pi);

  if(z1==z0 && TMath::IsNaN(RangRa)==true){
    RangRa=180-langRa;
  }
  if(z1==z0 && TMath::IsNaN(RangR)==true){
    RangR=180-langR;
  }
  if(z1==z0 && TMath::IsNaN(RangD)==true){
    RangD=180-langD;
  }

  if(z1!=z0 && TMath::IsNaN(RangRa)==true){
    RangRa=90;
  }
  if(z1!=z0 && TMath::IsNaN(RangR)==true){
    RangR=90;
  }
  if(z1!=z0 && TMath::IsNaN(RangD)==true){
    RangD=90;
  }

  if(attcal==true){
    const Double_t f0=0.0001, f2=3.16;
    //Calculate the attenuation
    // for frequency<1.Ghz. This is the integrand
    TF1 *fFr1=new TF1("fFr1","sqrt(1+pow(tan(asin([6]/([3]+[4]*exp([5]*x)))),2))*exp((((-6.22121-(-51.5 + x*(-4.5319e-3 + 5.822e-6*x))*(0.070927+(-51.5 + x*(-4.5319e-3 + 5.822e-6*x))*0.001773))*[0]-(-6.74890+(-51.5 + x*(-4.5319e-3 + 5.822e-6*x))*(0.026709-(-51.5 + x*(-4.5319e-3 + 5.822e-6*x))*0.000884))*[1])/([0]-[1]))+(((-6.22121-(-51.5 + x*(-4.5319e-3 + 5.822e-6*x))*(0.070927+(-51.5 + x*(-4.5319e-3 + 5.822e-6*x))*0.001773))-(-6.74890+(-51.5 + x*(-4.5319e-3 + 5.822e-6*x))*(0.026709-(-51.5 + x*(-4.5319e-3 + 5.822e-6*x))*0.000884)))/([1]-[0]))*[2])",-20000,0);
    
    fFr1->SetParameter(0,log10(f0));
    fFr1->SetParameter(1,0.0);
    fFr1->FixParameter(3,A);
    fFr1->FixParameter(4,B);
    fFr1->FixParameter(5,C);

    //sqrt(1+pow(tan(asin(lval/(A+B*exp(C*x)))),2))*exp((((-6.22121-(-51.5 + x*(-4.5319e-3 + 5.822e-6*x))*(0.070927+(-51.5 + x*(-4.5319e-3 + 5.822e-6*x))*0.001773))*log10(f0))/log10(f0))+(((-6.22121-(-51.5 + x*(-4.5319e-3 + 5.822e-6*x))*(0.070927+(-51.5 + x*(-4.5319e-3 + 5.822e-6*x))*0.001773))-(-6.74890+(-51.5 + x*(-4.5319e-3 + 5.822e-6*x))*(0.026709-(-51.5 + x*(-4.5319e-3 + 5.822e-6*x))*0.000884)))/(-log10(f0)))*freq)
    
    //loop over the different frequencies
    for(Int_t iat=0;iat<10;iat++){
      //set the frequency in Ghz
      fFr1->SetParameter(2,log10(0.1+iat*0.1));

      //Calculate attenuation for direct
      fFr1->SetParameter(6,lvalueD);
      Double_t totattD=1-fFr1->Integral(z0,z1);

      //Calculate attenuation for reflected
      fFr1->SetParameter(6,lvalueR);
      Double_t totattR=1-(fFr1->Integral(z0,-0.00001)+fFr1->Integral(z1,-0.00001));

      //Calculate attenuation for refracted
      fFr1->SetParameter(6,lvalueRa);
      Double_t totattRa=1-(fFr1->Integral(z0,-fabs((-1.0/C)*log((lvalueRa-A)/B))-0.00001)+fFr1->Integral(z1,-fabs((-1.0/C)*log((lvalueRa-A)/B))-0.00001));

      cout<<"the attenuation is "<<totattD<<" "<<totattR<<" "<<totattRa<<endl;
      delete fFr1;
    }
  
  // //for frequency>1.0Ghz
  // TF1 *fFr2=new TF1("fFr2","exp((((-4.09468-(-51.5 + x*(-4.5319e-3 + 5.822e-6*x))*(0.002213+(-51.5 + x*(-4.5319e-3 + 5.822e-6*x))*0.000332))*[0]-(-6.22121-(-51.5 + x*(-4.5319e-3 + 5.822e-6*x))*(0.070927+(-51.5 + x*(-4.5319e-3 + 5.822e-6*x))*0.001773))*[1])/([0]-[1]))+(((-4.09468-(-51.5 + x*(-4.5319e-3 + 5.822e-6*x))*(0.002213+(-51.5 + x*(-4.5319e-3 + 5.822e-6*x))*0.000332))-(-6.22121-(-51.5 + x*(-4.5319e-3 + 5.822e-6*x))*(0.070927+(-51.5 + x*(-4.5319e-3 + 5.822e-6*x))*0.001773)))/([1]-[0]))*[2])",1,3);
  // fFr2->SetParameter(0,0.0);
  // fFr2->SetParameter(1,log10(f2));
  // fFr2->SetParameter(2,log10(1.5));
  }
  
  output[0]=langD;
  output[1]=langR;
  output[2]=langRa;
  output[3]=timeD;
  output[4]=timeR;
  output[5]=timeRa;
  output[6]=RangD;
  output[7]=RangR;
  output[8]=RangRa;
  
  if(Flip==true){
    output[0]=180-RangD;
    output[1]=180-RangR;
    output[2]=180-RangRa;
    output[3]=timeD;
    output[4]=timeR;
    output[5]=timeRa;
    output[6]=180-langD;
    output[7]=180-langR;
    output[8]=180-langRa;
  }

  ////Plot the possible ray paths between 2 points
  if(Plot==true){

    Double_t h=1;
    Int_t dmax=-round(lowerz/h);
    dmax=100000;
    Double_t zn=z1;

    fD->SetParameter(3,lvalueD);
    TNtuple *nt1=new TNtuple("nt1","nt1","x:z");
    TNtuple *nt1ch=new TNtuple("nt1ch","nt1ch","x:z");
    for(Int_t i=0;i<dmax;i++){
      if(TMath::IsNaN(fD->Eval(zn))==false && Flip==false){
	nt1->Fill(fD->Eval(zn)-fD->Eval(z0),zn);
      }

      if(TMath::IsNaN(fD->Eval(zn))==false && Flip==true){
      	nt1->Fill(x1-(fD->Eval(zn)-fD->Eval(z0)),zn);
      }

      zn=zn-h;
      if(zn<z0){
	zn=z0;
	i=dmax+2;      
      }  
    }

    zn=z1;
    fR->SetParameter(3,lvalueR);
    TNtuple *nt2=new TNtuple("nt2","nt2","x:z");
    for(Int_t i=0;i<dmax;i++){
      if(TMath::IsNaN(fR->Eval(-zn))==false && zn<=0 && Flip==false){
	nt2->Fill((fR->Eval(-zn)-fR->Eval(-z0)+2*fabs(fR->Eval(1*pow(10,-10))-fR->Eval(-z0))),zn);
      }

      if(TMath::IsNaN(fR->Eval(-zn))==false && zn<=0 && Flip==true){
	nt2->Fill(x1-(fR->Eval(-zn)-fR->Eval(-z0)+2*fabs(fR->Eval(1*pow(10,-10))-fR->Eval(-z0))),zn);
      }
      
      zn=zn+h;
      if(zn>0){
	i=dmax+2;      
      }
    }

    fD->SetParameter(3,lvalueR);
    zn=-1*pow(10,-9);
    for(Int_t i=0;i<dmax;i++){  
      if(TMath::IsNaN(fD->Eval(zn))==false && Flip==false){
	nt2->Fill(fD->Eval(zn)-fD->Eval(z0),zn);
      }
      
      if(TMath::IsNaN(fD->Eval(zn))==false && Flip==true){
	nt2->Fill(x1-(fD->Eval(zn)-fD->Eval(z0)),zn);
      }

      zn=zn-h;
      if(zn<z0){
	zn=z0;
	i=dmax+2;
      }
    }

    zn=z1;
    fR->SetParameter(3,lvalueRa);    
    TNtuple *nt3=new TNtuple("nt3","nt3","x:z");
    for(Int_t i=0;i<dmax;i++){
      if(TMath::IsNaN(fR->Eval(-zn))==false && zn<=0 && Flip==false){
	nt3->Fill((fR->Eval(-zn)-fR->Eval(-z0)+2*fabs(fR->Eval(fabs((-1.0/C)*log((lvalueRa-A)/B)))-fR->Eval(-z0))),zn);
      }

      if(TMath::IsNaN(fR->Eval(-zn))==false && zn<=0 && Flip==true){
	nt3->Fill(x1-(fR->Eval(-zn)-fR->Eval(-z0)+2*fabs(fR->Eval(fabs((-1.0/C)*log((lvalueRa-A)/B)))-fR->Eval(-z0))),zn);
      }
    
      zn=zn+h;
      if(zn>-fabs((-1.0/C)*log((lvalueRa-A)/B))){
	i=dmax+2;      
      }
    }

    fD->SetParameter(3,lvalueRa);
    zn=-fabs((-1.0/C)*log((lvalueRa-A)/B));
    for(Int_t i=0;i<dmax;i++){  
      if(TMath::IsNaN(fD->Eval(zn))==false && Flip==false){
	nt3->Fill(fD->Eval(zn)-fD->Eval(z0),zn);
      }

      if(TMath::IsNaN(fD->Eval(zn))==false && Flip==true){
	nt3->Fill(x1-(fD->Eval(zn)-fD->Eval(z0)),zn);
      }
    
      zn=zn-h;
      if(zn<z0){
	zn=z0;
	i=dmax+2;      
      }
    }

    if(Flip==true){
      dsw=z0;
      z0=z1;
      z1=dsw;
    }

    TNtuple *nt4=new TNtuple("nt4","nt4","x:z");
    //nt4->Fill(0,z0-100);
    nt4->Fill(x1,z1);
    nt4->SetMarkerStyle(20);
    nt4->SetMarkerColor(2);

    nt1->SetMarkerStyle(20);
    nt1->SetMarkerColor(2);

    nt2->SetMarkerStyle(20);
    nt2->SetMarkerColor(2);

    TNtuple *nt5=new TNtuple("nt5","nt5","x:z");
    nt5->Fill(0,lowerz-50);
    nt5->Fill(x1+50,0);
    
    gStyle->SetTitleX(0.2);
    gStyle->SetTitleAlign(23);
    TCanvas *c2=new TCanvas("c2","c2");
    c2->cd();
    nt5->Draw("z:x","","P");
    TH2D *htemp = (TH2D*)gPad->GetPrimitive("htemp");
    // htemp->GetXaxis()->SetTitle("my new X title");
    // htemp->GetYaxis()->SetTitle("my new Y title"); 
    htemp->GetXaxis()->SetNdivisions(20);
    c2->SetGridx();
    TString title="Depth vs Distance, Tx at x=";
    title+=x0;
    title+=" m,z=";
    title+=z0;
    title+=" m, Rx at x=";
    title+=x1;
    title+=" m,z=";
    title+=z1;
    title+=" m; Distance (m);Depth (m)";
    htemp->SetTitle(title);
    
    nt4->Draw("z:x","","SAME P");
    if(fabs(checkzeroR)<0.5){
      nt2->Draw("z:x","","SAME L");
    }
    if(fabs(checkzeroRa)<0.5 && timeRa!=0){
      //cout<<"checkzeroRa is "<<checkzeroRa<<endl;
      nt3->Draw("z:x","","SAME L");
     }
    if(fabs(checkzeroD)<0.5){
      nt1->Draw("z:x","","SAME L");
    }
    
  }
  
  cout<<0<<" ,x0= "<<x0<<" ,z0= "<<z0<<" ,x1= "<<x1<<" ,z1= "<<z1<<" ,langRa= "<<output[2]<<" ,langR= "<<output[1]<<" ,langD= "<<output[0]<<" ,langD-langR= "<<output[0]-output[1]<<" ,langD-langRa= "<<output[0]-output[2]<<" ,RangRa= "<<output[8]<<" ,RangR= "<<output[7]<<" ,RangD= "<<output[6]<<" ,RangR-RangD= "<<output[7]-output[6]<<" ,RangRa-RangD= "<<output[8]-output[6]<<" ,timeRa= "<<output[5]<<" ,timeR= "<<output[4]<<" ,timeD= "<<output[3]<<" ,timeR-timeD= "<<output[4]-output[3]<<" ,timeRa-timeD= "<<output[5]-output[3]<<" ,lvalueRa "<<lvalueRa<<" ,lvalueR "<<lvalueR<<" "<<" ,lvalueD "<<lvalueD<<" ,checkzeroRa "<<checkzeroRa<<" ,checkzeroR "<<checkzeroR<<" ,checkzeroD "<<checkzeroD<<endl;
  
  if(fabs(checkzeroD)>0.5){
    output[6]=0;
  }
  if(fabs(checkzeroR)>0.5){
    output[7]=0;
  }
  if(fabs(checkzeroRa)>0.5){
    output[8]=0;
  }
  
  watch.Stop();
  //watch.Print("u");

  delete fR;
  delete fD;
  delete fDa;
  delete ftimeD;
  delete fRa;
  delete ftimeR;
  delete fRaa;
  delete ftimeRa;
  
  return output;
}
