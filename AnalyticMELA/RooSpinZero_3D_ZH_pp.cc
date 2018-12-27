/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "RooSpinZero_3D_ZH_pp.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 

#include "TLorentzRotation.h"
#include "TLorentzVector.h"

#include <math.h>
#include "TMath.h"

using namespace TMath;

enum parameterizationList {kMagPhase_As=0,kRealImag_Gs=1,kFracPhase_Gs=2,kNUMparameterizations=3};

ClassImp(RooSpinZero_3D_ZH_pp) 

RooSpinZero_3D_ZH_pp::RooSpinZero_3D_ZH_pp(const char *name, const char *title, 
                                           RooAbsReal& _h1,
                                           RooAbsReal& _h2,
                                           RooAbsReal& _Phi,
                                           RooAbsReal& _m,
                                           RooAbsReal& _Y,                                                                          
                                           RooAbsReal& _sqrts,
                                           RooAbsReal& _mX,
                                           RooAbsReal& _mZ,
                                           RooAbsReal& _R1Val,
                                           RooAbsReal& _R2Val,
                                           int _parameterization,
                                           RooAbsReal& _a1Val,
                                           RooAbsReal& _phi1Val,
                                           RooAbsReal& _a2Val,
                                           RooAbsReal& _phi2Val,
                                           RooAbsReal& _a3Val,
                                           RooAbsReal& _phi3Val,
                                           RooAbsReal& _g1Val,
                                           RooAbsReal& _g2Val,
                                           RooAbsReal& _g3Val,
                                           RooAbsReal& _g4Val,
                                           RooAbsReal& _g1ValIm,
                                           RooAbsReal& _g2ValIm,
                                           RooAbsReal& _g3ValIm,
                                           RooAbsReal& _g4ValIm,
                                           RooAbsReal& _fa2,
                                           RooAbsReal& _fa3,
                                           RooAbsReal& _phia2,
                                           RooAbsReal& _phia3,
                                           bool      _withAcc):
RooAbsPdf(name,title), 
h1("h1","h1",this,_h1),
h2("h2","h2",this,_h2),
Phi("Phi","Phi",this,_Phi),
m("m","m",this,_m),
Y("Y","Y",this,_Y),
sqrts("sqrts","sqrts",this,_sqrts),
mX("mX","mX",this,_mX),
mZ("mZ","mZ",this,_mZ),
R1Val("R1Val","R1Val",this,_R1Val),
R2Val("R2Val","R2Val",this,_R2Val),
parameterization(_parameterization),
a1Val("a1Val","a1Val",this,_a1Val),
phi1Val("phi1Val","phi1Val",this,_phi1Val),
a2Val("a2Val","a2Val",this,_a2Val),
phi2Val("phi2Val","phi2Val",this,_phi2Val),
a3Val("a3Val","a3Val",this,_a3Val),
phi3Val("phi3Val","phi3Val",this,_phi3Val),
g1Val("g1Val","g1Val",this,_g1Val),
g2Val("g2Val","g2Val",this,_g2Val),
g3Val("g3Val","g3Val",this,_g3Val),
g4Val("g4Val","g4Val",this,_g4Val),
g1ValIm("g1ValIm","g1ValIm",this,_g1ValIm),
g2ValIm("g2ValIm","g2ValIm",this,_g2ValIm),
g3ValIm("g3ValIm","g3ValIm",this,_g3ValIm),
g4ValIm("g4ValIm","g4ValIm",this,_g4ValIm),
fa2("fa2","fa2",this,_fa2),
fa3("fa3","fa3",this,_fa3),
phia2("phia2","phia2",this,_phia2),
phia3("phia3","phia3",this,_phia3),
withAcc(_withAcc)
{ 
} 


RooSpinZero_3D_ZH_pp::RooSpinZero_3D_ZH_pp(const RooSpinZero_3D_ZH_pp& other, const char* name) :  
RooAbsPdf(other,name), 
h1("h1",this,other.h1),
h2("h2",this,other.h2),
Phi("Phi",this,other.Phi),
m("m",this,other.m),
Y("Y",this,other.Y),
sqrts("sqrts",this,other.sqrts),
mX("mX",this,other.mX),
mZ("mZ",this,other.mZ),
R1Val("R1Val",this,other.R1Val),
R2Val("R2Val",this,other.R2Val),
parameterization(other.parameterization),
a1Val("a1Val",this,other.a1Val),
phi1Val("phi1Val",this,other.phi1Val),
a2Val("a2Val",this,other.a2Val),
phi2Val("phi2Val",this,other.phi2Val),
a3Val("a3Val",this,other.a3Val),
phi3Val("phi3Val",this,other.phi3Val),
g1Val("g1Val",this,other.g1Val),
g2Val("a2Val",this,other.g2Val),
g3Val("g3Val",this,other.g3Val),
g4Val("g4Val",this,other.g4Val),
g1ValIm("g1ValIm",this,other.g1ValIm),
g2ValIm("a2ValIm",this,other.g2ValIm),
g3ValIm("g3ValIm",this,other.g3ValIm),
g4ValIm("g4ValIm",this,other.g4ValIm),
fa2("fa2",this,other.fa2),
fa3("fa3",this,other.fa3),
phia2("phia2",this,other.phia2),
phia3("phia3",this,other.phia3),
withAcc(other.withAcc)
{ 
} 



Double_t RooSpinZero_3D_ZH_pp::evaluate() const 
{ 
    
    // these amplitudes are calculated based on comparing equations to the PRD paper
    // the beta and gamma are velocities of the Z in the C.O.M
    // Double_t betaValSquared =(pow(sqrts,2)-(pow(mX+mZ,2)))*(pow(sqrts,2)-(pow(mX-mZ,2)))/pow(pow(sqrts,2)-mX*mX+mZ*mZ,2);
    // Double_t gamma = 1./sqrt(1-betaValSquared);
    // http://prd.aps.org/pdf/PRD/v49/i1/p79_1
    // Double_t f00 = gamma*gamma/(2+gamma*gamma);
    
    
    // check whether event is in acceptance or not
    
    
    if ( withAcc ) {
        
        vector<TLorentzVector> lep_4vecs = Calculate4Momentum(m,mZ,mX,acos(h1),acos(h2),acos(0),Phi,0);
	// now boost the 4leptons to the original frame
	TLorentzVector pZstar_new;
	// calculate pz and E based on m and Y
	double E_Zstar_new =  m*(exp(2*Y)+1)/(2*exp(Y));
	double pz_Zstar_new = (exp(2*Y)-1)/(exp(2*Y)+1)*E_Zstar_new;
	pZstar_new.SetPxPyPzE(0, 0, pz_Zstar_new, E_Zstar_new);
	TVector3 boost_pZstar = pZstar_new.BoostVector();
		
	for (int i = 0 ; i < 4 ; i++ )
	  lep_4vecs[i].Boost(boost_pZstar); 
        
        double pt_plus_sq = pow(lep_4vecs[1].Px(),2) + pow(lep_4vecs[1].Py(),2);
        double pt_minus_sq = pow(lep_4vecs[0].Px(),2) + pow(lep_4vecs[0].Py(),2);
        
        double sinh2_eta_plus = pow(lep_4vecs[1].Pz(),2)/pt_plus_sq;
        double sinh2_eta_minus = pow(lep_4vecs[0].Pz(),2)/pt_minus_sq;
        
	double higgsPtSq = pow((lep_4vecs[2]+lep_4vecs[3]).Pt(),2);
	
	// std::cout << __LINE__ << "\t " << higgsPtSq << "\n";
	if( pt_minus_sq<400.0 
           || pt_plus_sq<400.0 
           || sinh2_eta_minus>pow(sinh(2.4),2) 
           || sinh2_eta_plus>pow(sinh(2.4),2) 
           || higgsPtSq < (200.*200)
           ) return 1e-30;
    }
    //-------------------------------------------------
    
    // below calcualtions are based on the H->ZZ amplitudes 
    double s=-(mX*mX-m*m-mZ*mZ)/2.;
    double kappa=s/(1000*1000);
    
    double a1=0,a2=0,a3=0,a1Im=0,a2Im=0,a3Im=0;
    double g1(1.0), g1Im(0.), g2(0.), g2Im(0.), g3(0.), g3Im(0.), g4(0.), g4Im(0.);
    
    if(parameterization==kMagPhase_As){
        a1=a1Val;
        a1Im=phi1Val;
        a2=a2Val;
        a2Im=phi2Val;
        a3=a3Val;
        a3Im=phi3Val;
    } else {
        
        if(parameterization==kFracPhase_Gs){
            
            double nanval = sqrt(1 - fa2 - fa3);
            if (nanval != nanval) return 0.0;
            
	    // LHC numbers at mH= 125 GeV (JHUGen) based on Decay
	    Double_t sigma1_p = 1.860351; // was 2.03971 at 126 GeV
	    Double_t sigma2_p = 0.672859; // was 0.77517 at 126 GeV
	    Double_t sigma4_p = 0.284353; // was 0.32689 at 126 GeV
            
            double g1Mag = 1.;
            double g2Mag = 0.;
            double g4Mag = 0.;
            
            if ( fa2 == 1. ) {
                g1Mag = 0.;
                g2Mag = 1.;
                g4Mag = 0.; 
            } else if ( fa3 == 1. ) {
                g1Mag = 0.;
                g2Mag = 0.;
                g4Mag = 1.; 
            } else {
                g2Mag = sqrt(fa2/(1.-fa2-fa3))*sqrt(sigma1_p/sigma2_p); 
                g4Mag = sqrt(fa3/(1.-fa2-fa3))*sqrt(sigma1_p/sigma4_p); 
            }
            
            g1   = g1Mag;
            g1Im = 0.0;
            g2   = - g2Mag*cos(phia2);
            g2Im = - g2Mag*sin(phia2);
            g3   = 0.0;
            g3Im = 0.0;
            g4   = - g4Mag*cos(phia3);
            g4Im = - g4Mag*sin(phia3);
            
            
        }else if(parameterization==kRealImag_Gs){
            
            g1   =  g1Val;
            g1Im =  g1ValIm;
            g2   =  - g2Val;
            g2Im =  - g2ValIm;
            g3   = g3Val;
            g3Im =  g3ValIm;
            g4   = - g4Val;
            g4Im = - g4ValIm;
        }
        
        a1 = g1*mZ*mZ/(mX*mX) + g2*2.*s/(mX*mX) + g3*kappa*s/(mX*mX);
        a1Im = g1Im*mZ*mZ/(mX*mX) + g2Im*2.*s/(mX*mX) + g3Im*kappa*s/(mX*mX);
        a2 = -2.*g2 - g3*kappa;
        a2Im = -2.*g2Im - g3Im*kappa;
        a3 = -2.*g4;
        a3Im = -2.*g4Im;
        
    }
    
    
    Double_t x = pow((mX*mX-m*m-mZ*mZ)/(2.*m*mZ),2)-1;
    
    Double_t A00Real = - (a1*sqrt(1+x) + a2*(mZ*m)/(mX*mX)*x);
    Double_t A00Imag = - (a1Im*sqrt(1+x) + a2Im*(mZ*m)/(mX*mX)*x);
    
    Double_t Ap0Real = a1 - a3Im*(mZ*m)/(mX*mX)*sqrt(x);
    Double_t Ap0Imag = a1Im + a3*(mZ*m)/(mX*mX)*sqrt(x);
    
    Double_t Am0Real = a1 + a3Im*(mZ*m)/(mX*mX)*sqrt(x);
    Double_t Am0Imag = a1Im - a3*(mZ*m)/(mX*mX)*sqrt(x);
    
    Double_t f00 = A00Real*A00Real + A00Imag*A00Imag;
    Double_t fp0 = Ap0Real*Ap0Real + Ap0Imag*Ap0Imag;
    Double_t fm0 = Am0Real*Am0Real + Am0Imag*Am0Imag;
    
    Double_t phi00=atan2(A00Imag,A00Real);
    Double_t phip0=atan2(Ap0Imag,Ap0Real)-phi00;
    Double_t phim0=atan2(Am0Imag,Am0Real)-phi00;
    
    Double_t value = 0;
    
    
    value += (f00*(-1 + Power(h1,2))*(-1 + Power(h2,2)))/4.;
    
    value += (fp0*(1 + Power(h1,2) - 2*h1*R1Val)*(1 + Power(h2,2) + 2*h2*R2Val))/16.;
    
    value += (fm0*(1 + Power(h1,2) + 2*h1*R1Val)*(1 + Power(h2,2) - 2*h2*R2Val))/16.;
    
    value += -(Sqrt(f00)*Sqrt(fp0)*Sqrt(1 - Power(h1,2))*Sqrt(1 - Power(h2,2))*(h1 - R1Val)*(h2 + R2Val)*Cos(Phi + phip0))/4.;
    
    value += -(Sqrt(f00)*Sqrt(fm0)*Sqrt(1 - Power(h1,2))*Sqrt(1 - Power(h2,2))*(h1 + R1Val)*(h2 - R2Val)*Cos(Phi - phim0))/4.;
    
    value += (Sqrt(fm0)*Sqrt(fp0)*(-1 + Power(h1,2))*(-1 + Power(h2,2))*Cos(2*Phi - phim0 + phip0))/8.;
    
    // Get partonic lumi
    value *= partonicLuminosity( m, Y, sqrts );
    // Get phasespace factor
    if ((m*m-(mZ+mX)*(mZ+mX))*(m*m-(mZ-mX)*(mZ-mX)) < 0 ){
        return 1.e-9;
    }
    double beta_psf = sqrt( (m*m-(mZ+mX)*(mZ+mX))*(m*m-(mZ-mX)*(mZ-mX)))/(m*m);
    double sigma_SM_psf = beta_psf*(beta_psf*beta_psf+12.*mZ*mZ/(m*m))/(m*m*(1-mZ*mZ/(m*m))*(1-mZ*mZ/(m*m)));

    double angularPartIntegral_SM = 0.;
    // calculate the SM integrals
    // use g1 = 1., and all others 0.
    double a1_SM = mZ*mZ/(mX*mX);
    double a1Im_SM =  0.;
    double a2_SM  = 0.;
    double a2Im_SM = 0.;
    double a3_SM = 0.;
    double a3Im_SM = 0.;

    Double_t A00Real_SM = - (a1_SM*sqrt(1+x) + a2_SM*(mZ*m)/(mX*mX)*x);
    Double_t A00Imag_SM = - (a1Im_SM*sqrt(1+x) + a2Im_SM*(mZ*m)/(mX*mX)*x);
    
    Double_t Ap0Real_SM = a1_SM - a3Im_SM*(mZ*m)/(mX*mX)*sqrt(x);
    Double_t Ap0Imag_SM = a1Im_SM + a3_SM*(mZ*m)/(mX*mX)*sqrt(x);
    
    Double_t Am0Real_SM = a1_SM + a3Im_SM*(mZ*m)/(mX*mX)*sqrt(x);
    Double_t Am0Imag_SM = a1Im_SM - a3_SM*(mZ*m)/(mX*mX)*sqrt(x);
    
    Double_t f00_SM = A00Real_SM*A00Real_SM + A00Imag_SM*A00Imag_SM;
    Double_t fp0_SM = Ap0Real_SM*Ap0Real_SM + Ap0Imag_SM*Ap0Imag_SM;
    Double_t fm0_SM = Am0Real_SM*Am0Real_SM + Am0Imag_SM*Am0Imag_SM;
    
    angularPartIntegral_SM += (32*f00_SM*Power(Pi(),2))/9.;
    angularPartIntegral_SM += (32*fp0_SM*Power(Pi(),2))/9.;
    angularPartIntegral_SM += (32*fm0_SM*Power(Pi(),2))/9.;
    double phaseSpaceFactor = sigma_SM_psf/angularPartIntegral_SM;

    value *= phaseSpaceFactor;
    
    //std::cout << "phaseSpaceFactor = " << phaseSpaceFactor << ", sigma_SM_psf = " << sigma_SM_psf << ", angularPartIntegral = " << angularPartIntegral << std::endl;
    
    return value; 
} 

Int_t RooSpinZero_3D_ZH_pp::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
    
    if ( !withAcc ) {
        if (matchArgs(allVars,analVars,RooArgSet(*h1.absArg(),*h2.absArg(),*Phi.absArg()))) return 4 ;
        if (matchArgs(allVars,analVars,h1,h2)) return 1 ;
        if (matchArgs(allVars,analVars,h1,Phi)) return 2 ;
        if (matchArgs(allVars,analVars,h2,Phi)) return 3 ;
    }
    return 0 ;
}

Double_t RooSpinZero_3D_ZH_pp::analyticalIntegral(Int_t code, const char* rangeName) const
{
    
    // check whether event is in acceptance or not
    
    if ( withAcc ) {
        
        vector<TLorentzVector> lep_4vecs = Calculate4Momentum(m,mZ,mX,acos(h1),acos(h2),acos(0),Phi,0);
        
	// now boost the 4leptons to the original frame
	TLorentzVector pZstar_new;
	// calculate pz and E based on m and Y
	double E_Zstar_new =  m*(exp(2*Y)+1)/(2*exp(Y));
	double pz_Zstar_new = (exp(2*Y)-1)/(exp(2*Y)+1)*E_Zstar_new;
	pZstar_new.SetPxPyPzE(0, 0, pz_Zstar_new, E_Zstar_new);
	TVector3 boost_pZstar = pZstar_new.BoostVector();
	
	for (int i = 0 ; i < 4 ; i++ )
	  lep_4vecs[i].Boost(boost_pZstar); 
	
        double pt_plus_sq = pow(lep_4vecs[1].Px(),2) + pow(lep_4vecs[1].Py(),2);
        double pt_minus_sq = pow(lep_4vecs[0].Px(),2) + pow(lep_4vecs[0].Py(),2);
        
        double sinh2_eta_plus = pow(lep_4vecs[1].Pz(),2)/pt_plus_sq;
        double sinh2_eta_minus = pow(lep_4vecs[0].Pz(),2)/pt_minus_sq;

	double higgsPtSq = pow((lep_4vecs[2]+lep_4vecs[3]).Pt(),2);
	
	if( pt_minus_sq<400.0 
           || pt_plus_sq<400.0 
           || sinh2_eta_minus>pow(sinh(2.4),2) 
           || sinh2_eta_plus>pow(sinh(2.4),2) 
           || higgsPtSq < (200.*200)
           ) return 1e-30;

    }
    
    // the beta and gamma are velocities of the Z in the C.O.M
    Double_t betaValSquared =(pow(m,2)-(pow(mX+mZ,2)))*(pow(m,2)-(pow(mX-mZ,2)))/pow(pow(m,2)-mX*mX+mZ*mZ,2);
    Double_t gamma = 1./sqrt(1-betaValSquared);
    
    // these amplitudes are calculated based on comparing equations to the PRD paper
    // http://prd.aps.org/pdf/PRD/v49/i1/p79_1
    // Double_t f00 = gamma*gamma/(2+gamma*gamma);
    
    // below calcualtions are based on the H->ZZ amplitudes 
    double s=-(mX*mX-m*m-mZ*mZ)/2.;
    double kappa=s/(1000*1000);
    
    double a1=0,a2=0,a3=0,a1Im=0,a2Im=0,a3Im=0;
    double g1(1.0), g1Im(0.), g2(0.), g2Im(0.), g3(0.), g3Im(0.), g4(0.), g4Im(0.);
    
    if(parameterization==kMagPhase_As){
        a1=a1Val;
        a1Im=phi1Val;
        a2=a2Val;
        a2Im=phi2Val;
        a3=a3Val;
        a3Im=phi3Val;
    } else {
        
        if(parameterization==kFracPhase_Gs){
            // 
            // place holder, not correct
            // 
            double nanval = sqrt(1 - fa2 - fa3);
            if (nanval != nanval) return 0.0;
            
	    // LHC numbers at mH= 125 GeV (JHUGen) based on Decay
	    Double_t sigma1_p = 1.860351; // was 2.03971 at 126 GeV
	    Double_t sigma2_p = 0.672859; // was 0.77517 at 126 GeV
	    Double_t sigma4_p = 0.284353; // was 0.32689 at 126 GeV
            
            double g1Mag = 1.;
            double g2Mag = 0.;
            double g4Mag = 0.;
            
            if ( fa2 == 1. ) {
                g1Mag = 0.;
                g2Mag = 1.;
                g4Mag = 0.; 
            } else if ( fa3 == 1. ) {
                g1Mag = 0.;
                g2Mag = 0.;
                g4Mag = 1.; 
            } else {
                g2Mag = sqrt(fa2/(1.-fa2-fa3))*sqrt(sigma1_p/sigma2_p); 
                g4Mag = sqrt(fa3/(1.-fa2-fa3))*sqrt(sigma1_p/sigma4_p); 
            }
            
            g1   = g1Mag;
            g1Im = 0.0;
            g2   = - g2Mag*cos(phia2);
            g2Im = - g2Mag*sin(phia2);
            g3   = 0.0;
            g3Im = 0.0;
            g4   = - g4Mag*cos(phia3);
            g4Im = - g4Mag*sin(phia3);
            
            
        }else if(parameterization==kRealImag_Gs){
            
            g1   =  g1Val;
            g1Im =  g1ValIm;
            g2   =  -g2Val;
            g2Im =  -g2ValIm;
            g3   = g3Val;
            g3Im = g3ValIm;
            g4   = - g4Val;
            g4Im = - g4ValIm;
        }
        
        a1 = g1*mZ*mZ/(mX*mX) + g2*2.*s/(mX*mX) + g3*kappa*s/(mX*mX);
        a1Im = g1Im*mZ*mZ/(mX*mX) + g2Im*2.*s/(mX*mX) + g3Im*kappa*s/(mX*mX);
        a2 = -2.*g2 - g3*kappa;
        a2Im = -2.*g2Im - g3Im*kappa;
        a3 =  - 2.*g4;
        a3Im = -2.*g4Im;
        
    }
    
    Double_t x = pow((mX*mX-m*m-mZ*mZ)/(2.*m*mZ),2)-1;
    
    Double_t A00Real = - (a1*sqrt(1+x) + a2*(mZ*m)/(mX*mX)*x);
    Double_t A00Imag = - (a1Im*sqrt(1+x) + a2Im*(mZ*m)/(mX*mX)*x);
    
    Double_t Ap0Real = a1 - a3Im*(mZ*m)/(mX*mX)*sqrt(x);
    Double_t Ap0Imag = a1Im + a3*(mZ*m)/(mX*mX)*sqrt(x);
    
    Double_t Am0Real = a1 + a3Im*(mZ*m)/(mX*mX)*sqrt(x);
    Double_t Am0Imag = a1Im - a3*(mZ*m)/(mX*mX)*sqrt(x);
    
    Double_t f00 = A00Real*A00Real + A00Imag*A00Imag;
    Double_t fp0 = Ap0Real*Ap0Real + Ap0Imag*Ap0Imag;
    Double_t fm0 = Am0Real*Am0Real + Am0Imag*Am0Imag;
    
    Double_t phi00=atan2(A00Imag,A00Real);
    Double_t phip0=atan2(Ap0Imag,Ap0Real)-phi00;
    Double_t phim0=atan2(Am0Imag,Am0Real)-phi00;
    
    // Get phasespace factor
    if ((m*m-(mZ+mX)*(mZ+mX))*(m*m-(mZ-mX)*(mZ-mX)) < 0 ){
        return 1.e-9;
    }
    double beta_psf = sqrt( (m*m-(mZ+mX)*(mZ+mX))*(m*m-(mZ-mX)*(mZ-mX)))/(m*m);
    double sigma_SM_psf = beta_psf*(beta_psf*beta_psf+12.*mZ*mZ/(m*m))/(m*m*(1-mZ*mZ/(m*m))*(1-mZ*mZ/(m*m)));
    double angularPartIntegral_SM = 0.;
    // calculate the SM integrals
    // use g1 = 1., and all others 0.
    double a1_SM = mZ*mZ/(mX*mX);
    double a1Im_SM =  0.;
    double a2_SM  = 0.;
    double a2Im_SM = 0.;
    double a3_SM = 0.;
    double a3Im_SM = 0.;

    Double_t A00Real_SM = - (a1_SM*sqrt(1+x) + a2_SM*(mZ*m)/(mX*mX)*x);
    Double_t A00Imag_SM = - (a1Im_SM*sqrt(1+x) + a2Im_SM*(mZ*m)/(mX*mX)*x);
    
    Double_t Ap0Real_SM = a1_SM - a3Im_SM*(mZ*m)/(mX*mX)*sqrt(x);
    Double_t Ap0Imag_SM = a1Im_SM + a3_SM*(mZ*m)/(mX*mX)*sqrt(x);
    
    Double_t Am0Real_SM = a1_SM + a3Im_SM*(mZ*m)/(mX*mX)*sqrt(x);
    Double_t Am0Imag_SM = a1Im_SM - a3_SM*(mZ*m)/(mX*mX)*sqrt(x);
    
    Double_t f00_SM = A00Real_SM*A00Real_SM + A00Imag_SM*A00Imag_SM;
    Double_t fp0_SM = Ap0Real_SM*Ap0Real_SM + Ap0Imag_SM*Ap0Imag_SM;
    Double_t fm0_SM = Am0Real_SM*Am0Real_SM + Am0Imag_SM*Am0Imag_SM;
    
    angularPartIntegral_SM += (32*f00_SM*Power(Pi(),2))/9.;
    angularPartIntegral_SM += (32*fp0_SM*Power(Pi(),2))/9.;
    angularPartIntegral_SM += (32*fm0_SM*Power(Pi(),2))/9.;
    double phaseSpaceFactor = sigma_SM_psf/angularPartIntegral_SM;
    
    //std::cout << "phaseSpaceFactor = " << phaseSpaceFactor << ", sigma_SM_psf = " << sigma_SM_psf << ", angularPartIntegral = " << angularPartIntegral << std::endl;

    
    switch(code)
    {
        // projections to h2
        case 2:
        {
        
        
        double value = 0.;
        
        value += (-8*f00*(-1 + Power(h2,2))*Power(Pi(),2))/3.;
        
        value += (4*fp0*Power(Pi(),2)*(1 + Power(h2,2) + 2*h2*R2Val))/3.;
        
        value += (4*fm0*Power(Pi(),2)*(1 + Power(h2,2) - 2*h2*R2Val))/3.;
        
        value *= partonicLuminosity( m, Y, sqrts );
        value *= phaseSpaceFactor;
                
        return value/(4*Pi());
        
        }
        // projections to h1
        case 3:
        {
        
        double value = 0;
        
        value += (-8*f00*(-1 + Power(h1,2))*Power(Pi(),2))/3.;
        
        value += (4*fp0*Power(Pi(),2)*(1 + Power(h1,2) - 2*h1*R1Val))/3.;
        
        value += (4*fm0*Power(Pi(),2)*(1 + Power(h1,2) + 2*h1*R1Val))/3.;
        
        value *= partonicLuminosity( m, Y, sqrts );
        value *= phaseSpaceFactor;
                
        return value/(4*Pi());
        
        }
        // projections to Phi
        case 1:
        {
        
        
        double value = 0.;
        
        value += (16*f00*Pi())/9.;
        
        value += (16*fp0*Pi())/9.;
        
        value += (16*fm0*Pi())/9.;
        
        value += (Sqrt(f00)*Sqrt(fp0)*Power(Pi(),3)*R1Val*R2Val*Cos(Phi + phip0))/4.;
        
        value += (Sqrt(f00)*Sqrt(fm0)*Power(Pi(),3)*R1Val*R2Val*Cos(Phi - phim0))/4.;
        
        value += (8*sqrt(fm0)*sqrt(fp0)*Pi()*Cos(2*Phi - phim0 + phip0))/9.;

        value *= partonicLuminosity( m, Y, sqrts );
        value *= phaseSpaceFactor;
                
        return value/(4*Pi());
        
        }
        // projected everything
        case 4:
        {
        double value = 0.;
        
        value += (32*f00*Power(Pi(),2))/9.;
        
        value += (32*fp0*Power(Pi(),2))/9.;
        
        value += (32*fm0*Power(Pi(),2))/9.;
        
        value *= partonicLuminosity( m, Y, sqrts ); 
        value *= phaseSpaceFactor;
                
        return value/(4*Pi());
        }
    }
    assert(0) ;
    return 0 ;
}

vector<TLorentzVector> RooSpinZero_3D_ZH_pp::Calculate4Momentum(float Mx,float M1,float M2,float theta,float theta1,float theta2,float _Phi1_,float _Phi_) const 
{
    
    float phi1,phi2;
    phi1=TMath::Pi()-_Phi1_;
    phi2=_Phi1_+_Phi_;
    
    
    float gamma1,gamma2,beta1,beta2;
    
    gamma1=(Mx*Mx+M1*M1-M2*M2)/(2*Mx*M1);
    gamma2=(Mx*Mx-M1*M1+M2*M2)/(2*Mx*M2);
    beta1=sqrt(1-1/(gamma1*gamma1));
    beta2=sqrt(1-1/(gamma2*gamma2));
    
    //gluon 4 vectors
    TLorentzVector p1CM(0,0,Mx/2,Mx/2);
    TLorentzVector p2CM(0,0,-Mx/2,Mx/2);
    
    //vector boson 4 vectors
    TLorentzVector kZ1(gamma1*M1*sin(theta)*beta1,0, gamma1*M1*cos(theta)*beta1,gamma1*M1*1);   
    TLorentzVector kZ2(-gamma2*M2*sin(theta)*beta2,0, -gamma2*M2*cos(theta)*beta2,gamma2*M2*1);
    
    //Rotation and Boost matrices. Note gamma1*beta1*M1=gamma2*beta2*M2.
    
    TLorentzRotation Z1ToZ,Z2ToZ;
    
    Z1ToZ.Boost(0,0,beta1);
    Z2ToZ.Boost(0,0,beta2);
    Z1ToZ.RotateY(theta);
    Z2ToZ.RotateY(TMath::Pi()+theta);
    
    //fermons 4 vectors in vector boson rest frame
    
    TLorentzVector p3Z1((M1/2)*sin(theta1)*cos(phi1),(M1/2)*sin(theta1)*sin(phi1),(M1/2)*cos(theta1),(M1/2)*1);
    
    TLorentzVector p4Z1(-(M1/2)*sin(theta1)*cos(phi1),-(M1/2)*sin(theta1)*sin(phi1),-(M1/2)*cos(theta1),(M1/2)*1);
    
    TLorentzVector p5Z2((M2/2)*sin(theta2)*cos(phi2),(M2/2)*sin(theta2)*sin(phi2),(M2/2)*cos(theta2),(M2/2)*1);
    
    TLorentzVector p6Z2(-(M2/2)*sin(theta2)*cos(phi2),-(M2/2)*sin(theta2)*sin(phi2),-(M2/2)*cos(theta2),(M2/2)*1);
    
    
    // fermions 4 vectors in CM frame
    
    TLorentzVector p3CM,p4CM,p5CM,p6CM;
    
    p3CM=Z1ToZ*p3Z1;
    p4CM=Z1ToZ*p4Z1;
    p5CM=Z2ToZ*p5Z2;
    p6CM=Z2ToZ*p6Z2;
    
    vector<TLorentzVector> p;
    
    p.push_back(p3CM);
    p.push_back(p4CM);
    p.push_back(p5CM);
    p.push_back(p6CM);
    
    return p;
}

float RooSpinZero_3D_ZH_pp::partonicLuminosity(float mVal, float YVal, float sqrtsVal) const
{
    
    Double_t s0 = sqrtsVal*sqrtsVal;
	Double_t s = mVal*mVal;
	Double_t Q = mVal;
	Double_t xa = exp(YVal)*sqrt(s/s0);
	Double_t xb = exp(-YVal)*sqrt(s/s0);
    
	
	Double_t weightu = 0.5;
	Double_t weightd = 0.5;
	Double_t weightc = 1.0;
	Double_t weights = 1.0;
	Double_t weightb = 1.0;
    

	// PDF parameters
	// up params
	Double_t u0par0 = 0.03134; Double_t u0par1 =-2.068e-05; Double_t u0par2 = 1.283e-08; 
	Double_t u1par0 = 0.9; Double_t u1par1 =-0.0004307; Double_t u1par2 = 2.458e-07;
	Double_t u2par0 =-0.1369; Double_t u2par1 = 0.003423; Double_t u2par2 =-2.155e-06;
	Double_t u3par0 =-0.4013; Double_t u3par1 =-0.0002574; Double_t u3par2 = 1.561e-07;
	Double_t u4par0 = 0.5782; Double_t u4par1 =-0.004728; Double_t u4par2 = 2.906e-06;
	Double_t ubar0par0 = 0.02856; Double_t ubar0par1 =-2.112e-05; Double_t ubar0par2 = 1.272e-08;
	Double_t ubar1par0 =-0.06822; Double_t ubar1par1 = 3.172e-05; Double_t ubar1par2 =-2.008e-08;
	Double_t ubar2par0 = 0.1967; Double_t ubar2par1 =-0.000118; Double_t ubar2par2 = 6.871e-08;
	Double_t ubar3par0 =-0.2251; Double_t ubar3par1 = 0.0001295; Double_t ubar3par2 =-7.181e-08;
	Double_t ubar4par0 =-0.4068; Double_t ubar4par1 =-0.0002956; Double_t ubar4par2 = 1.783e-07;
	Double_t ubar5par0 =-2.251; Double_t ubar5par1 =-0.0001699; Double_t ubar5par2 = 1.492e-07;
	// down params
	Double_t d0par0 = 0.03278; Double_t d0par1 =-2.915e-05; Double_t d0par2 = 1.809e-08;
	Double_t d1par0 = 0.479; Double_t d1par1 =-0.0002559; Double_t d1par2 = 1.557e-07;
	Double_t d2par0 =-0.5972; Double_t d2par1 = 0.0003118; Double_t d2par2 =-1.905e-07;
	Double_t d3par0 =-0.3892; Double_t d3par1 =-0.000317; Double_t d3par2 = 1.944e-07;
	Double_t d4par0 = 0.5007; Double_t d4par1 =-0.001665; Double_t d4par2 = 9.895e-07;
	Double_t dbar0par0 = 0.02328; Double_t dbar0par1 =-1.367e-05; Double_t dbar0par2 = 8.246e-09;
	Double_t dbar1par0 = 0.09422; Double_t dbar1par1 =-0.0001019; Double_t dbar1par2 = 6.375e-08;
	Double_t dbar2par0 =-0.5296; Double_t dbar2par1 = 0.000466; Double_t dbar2par2 =-2.896e-07;
	Double_t dbar3par0 = 0.5354; Double_t dbar3par1 =-0.0004404; Double_t dbar3par2 = 2.728e-07;
	Double_t dbar4par0 =-0.4386; Double_t dbar4par1 =-0.0002605; Double_t dbar4par2 = 1.582e-07;
	Double_t dbar5par0 =-1.289; Double_t dbar5par1 =-0.001618; Double_t dbar5par2 = 9.601e-07;
	// charm, strange, bottom params
	Double_t c0par0 = 0.01829; Double_t c0par1 =-6.93e-06; Double_t c0par2 = 3.796e-09;
	Double_t c1par0 = 0.03081; Double_t c1par1 = 4.325e-05; Double_t c1par2 =-3.95e-08;
	Double_t c2par0 = 0.5398; Double_t c2par1 =-4.284e-05; Double_t c2par2 =-1.362e-08;
	Double_t c3par0 =-0.5986; Double_t c3par1 = 0.002565; Double_t c3par2 =-1.937e-06;
	Double_t c4par0 =-0.4534; Double_t c4par1 =-0.0002329; Double_t c4par2 = 1.343e-07;
	Double_t c5par0 =-8.657; Double_t c5par1 =-0.005157; Double_t c5par2 = 3.68e-06;
	Double_t s0par0 = 0.01312; Double_t s0par1 =-3.743e-06; Double_t s0par2 = 2.076e-09;
	Double_t s1par0 =-0.001416; Double_t s1par1 =-7.649e-06; Double_t s1par2 = 4.757e-09;
	Double_t s2par0 = 0.2864; Double_t s2par1 =-6.693e-05; Double_t s2par2 = 3.566e-08;
	Double_t s3par0 =-0.4857; Double_t s3par1 =-0.000253; Double_t s3par2 = 1.541e-07;
	Double_t s4par0 =-10.33; Double_t s4par1 =-0.001601; Double_t s4par2 = 9.718e-07;
	Double_t b0par0 = 0.005934; Double_t b0par1 = 2.516e-06; Double_t b0par2 =-1.828e-09;
	Double_t b1par0 =-0.003063; Double_t b1par1 =-6.761e-06; Double_t b1par2 = 4.298e-09;
	Double_t b2par0 = 0.1174; Double_t b2par1 = 3.752e-05; Double_t b2par2 =-2.863e-08;
	Double_t b3par0 =-0.5549; Double_t b3par1 =-0.0002205; Double_t b3par2 = 1.334e-07;
	Double_t b4par0 =-10.18; Double_t b4par1 =-0.001136; Double_t b4par2 = 6.931e-07;
    
	
	// PDF definition
	Double_t up0 = u0par0 + u0par1*Q + u0par2*Q*Q;
	Double_t up1 = u1par0 + u1par1*Q + u1par2*Q*Q;
	Double_t up2 = u2par0 + u2par1*Q + u2par2*Q*Q;
	Double_t up3 = u3par0 + u3par1*Q + u3par2*Q*Q;
	Double_t up4 = u4par0 + u4par1*Q + u4par2*Q*Q;
	Double_t antiup0 = ubar0par0 + ubar0par1*Q + ubar0par2*Q*Q;
	Double_t antiup1 = ubar1par0 + ubar1par1*Q + ubar1par2*Q*Q;
	Double_t antiup2 = ubar2par0 + ubar2par1*Q + ubar2par2*Q*Q;
	Double_t antiup3 = ubar3par0 + ubar3par1*Q + ubar3par2*Q*Q;
	Double_t antiup4 = ubar4par0 + ubar4par1*Q + ubar4par2*Q*Q;
	Double_t antiup5 = ubar5par0 + ubar5par1*Q + ubar5par2*Q*Q;
	Double_t down0 = d0par0 + d0par1*Q + d0par2*Q*Q;
	Double_t down1 = d1par0 + d1par1*Q + d1par2*Q*Q;
	Double_t down2 = d2par0 + d2par1*Q + d2par2*Q*Q;
	Double_t down3 = d3par0 + d3par1*Q + d3par2*Q*Q;
	Double_t down4 = d4par0 + d4par1*Q + d4par2*Q*Q;
	Double_t antidown0 = dbar0par0 + dbar0par1*Q + dbar0par2*Q*Q;
	Double_t antidown1 = dbar1par0 + dbar1par1*Q + dbar1par2*Q*Q;
	Double_t antidown2 = dbar2par0 + dbar2par1*Q + dbar2par2*Q*Q;
	Double_t antidown3 = dbar3par0 + dbar3par1*Q + dbar3par2*Q*Q;
	Double_t antidown4 = dbar4par0 + dbar4par1*Q + dbar4par2*Q*Q;
	Double_t antidown5 = dbar5par0 + dbar5par1*Q + dbar5par2*Q*Q;
	Double_t charm0 = c0par0 + c0par1*Q + c0par2*Q*Q;
	Double_t charm1 = c1par0 + c1par1*Q + c1par2*Q*Q;
	Double_t charm2 = c2par0 + c2par1*Q + c2par2*Q*Q;
	Double_t charm3 = c3par0 + c3par1*Q + c3par2*Q*Q;
	Double_t charm4 = c4par0 + c4par1*Q + c4par2*Q*Q;
	Double_t charm5 = c5par0 + c5par1*Q + c5par2*Q*Q;
	Double_t strange0 = s0par0 + s0par1*Q + s0par2*Q*Q;
	Double_t strange1 = s1par0 + s1par1*Q + s1par2*Q*Q;
	Double_t strange2 = s2par0 + s2par1*Q + s2par2*Q*Q;
	Double_t strange3 = s3par0 + s3par1*Q + s3par2*Q*Q;
	Double_t strange4 = s4par0 + s4par1*Q + s4par2*Q*Q;
	Double_t bottom0 = b0par0 + b0par1*Q + b0par2*Q*Q;
	Double_t bottom1 = b1par0 + b1par1*Q + b1par2*Q*Q;
	Double_t bottom2 = b2par0 + b2par1*Q + b2par2*Q*Q;
	Double_t bottom3 = b3par0 + b3par1*Q + b3par2*Q*Q;
	Double_t bottom4 = b4par0 + b4par1*Q + b4par2*Q*Q;
    
	Double_t FuncAu1 = (up0+up1*xa+up2*pow(xa,2))*pow((1-xa),4)*pow(xa,up3)*exp(1.0+up4*xa);
	Double_t FuncBu1 = (antiup0+antiup1*xb+antiup2*pow(xb,2)+antiup3*pow(xb,3))*pow((1-xb),4)*pow(xb,antiup4)*exp(1.0+antiup5*xb);
	Double_t FuncAu2 = (up0+up1*xb+up2*pow(xb,2))*pow((1-xb),4)*pow(xb,up3)*exp(1.0+up4*xb);
	Double_t FuncBu2 = (antiup0+antiup1*xa+antiup2*pow(xa,2)+antiup3*pow(xa,3))*pow((1-xa),4)*pow(xa,antiup4)*exp(1.0+antiup5*xa);
	Double_t FuncABu = FuncAu1/xa*FuncBu1/xb+FuncAu2/xa*FuncBu2/xb;
    
	
	Double_t FuncAd1 = (down0+down1*xa+down2*pow(xa,2))*pow((1-xa),4)*pow(xa,down3)*exp(1.0+down4*xa);
	Double_t FuncBd1 = (antidown0+antidown1*xb+antidown2*pow(xb,2)+antidown3*pow(xb,3))*pow((1-xb),4)*pow(xb,antidown4)*exp(1.0+antidown5*xb);
	Double_t FuncAd2 = (down0+down1*xb+down2*pow(xb,2))*pow((1-xb),4)*pow(xb,down3)*exp(1.0+down4*xb);
	Double_t FuncBd2 = (antidown0+antidown1*xa+antidown2*pow(xa,2)+antidown3*pow(xa,3))*pow((1-xa),4)*pow(xa,antidown4)*exp(1.0+antidown5*xa);
	Double_t FuncABd = FuncAd1/xa*FuncBd1/xb+FuncAd2/xa*FuncBd2/xb;
    
	Double_t Funcca = (charm0+charm1*xa+charm2*pow(xa,2)+charm3*pow(xa,3))*pow((1-xa),4)*pow(xa,charm4)*exp(1.0+charm5*xa);
	Double_t Funccb = (charm0+charm1*xb+charm2*pow(xb,2)+charm3*pow(xb,3))*pow((1-xb),4)*pow(xb,charm4)*exp(1.0+charm5*xb);
	Double_t Funcsa = (strange0+strange1*xa+strange2*pow(xa,2))*pow((1-xa),4)*pow(xa,strange3)*exp(1.0+strange4*xa);
	Double_t Funcsb = (strange0+strange1*xb+strange2*pow(xb,2))*pow((1-xb),4)*pow(xb,strange3)*exp(1.0+strange4*xb);
	Double_t Funcba = (bottom0+bottom1*xa+bottom2*pow(xa,2))*pow((1-xa),4)*pow(xa,bottom3)*exp(1.0+bottom4*xa);
	Double_t Funcbb = (bottom0+bottom1*xb+bottom2*pow(xb,2))*pow((1-xb),4)*pow(xb,bottom3)*exp(1.0+bottom4*xb);
	Double_t FuncABc = Funcsa*Funcsb/xa/xb;
	Double_t FuncABs = Funcca*Funccb/xa/xb;
	Double_t FuncABb = Funcba*Funcbb/xa/xb;
    
    
	Double_t totSec = 2*mVal*(
                           (FuncABu)*weightu
                           +(FuncABd)*weightd
                           +(FuncABc)*weightc
                           +(FuncABs)*weights
                           +(FuncABb)*weightb
                           );
    
	if(( mVal <= 600. && TMath::Abs(YVal) > 20*pow(float(mVal),float(-0.32))) || ( mVal > 600. && TMath::Abs(YVal) > 21*pow(float(mVal),float(-0.34))))
        {
	    //Find totSec when mZZ, YVal=0
	    Double_t xa0 = sqrt(s/s0); //at YVal=0 xa=xb
        
	    //up
	    //if xa=xb then FuncAu1=FuncAu2 and FuncBu1=FuncBu2
	    FuncAu1 = (up0+up1*xa0+up2*pow(xa0,2))*pow((1-xa0),4)*pow(xa0,up3)*exp(1.0+up4*xa0);
	    FuncBu1 = (antiup0+antiup1*xa0+antiup2*pow(xa0,2)+antiup3*pow(xa0,3))*pow((1-xa0),4)*pow(xa0,antiup4)*exp(1.0+antiup5*xa0);
	    FuncABu = 2*(FuncAu1/xa0*FuncBu1/xa0);
        
	    //down
	    //if xa=xb then FuncAd1=FuncAd2 and FuncBd1=FuncBd2
	    FuncAd1 = (down0+down1*xa0+down2*pow(xa0,2))*pow((1-xa0),4)*pow(xa0,down3)*exp(1.0+down4*xa0);
	    FuncBd1 = (antidown0+antidown1*xa0+antidown2*pow(xa0,2)+antidown3*pow(xa0,3))*pow((1-xa0),4)*pow(xa0,antidown4)*exp(1.0+antidown5*xa0);
	    FuncABd = 2*(FuncAd1/xa0*FuncBd1/xa0);
        
	    //sea
	    Funcca = (charm0+charm1*xa0+charm2*pow(xa0,2)+charm3*pow(xa0,3))*pow((1-xa0),4)*pow(xa0,charm4)*exp(1.0+charm5*xa0); //Funcca=Funccb
	    Funcsa = (strange0+strange1*xa0+strange2*pow(xa0,2))*pow((1-xa0),4)*pow(xa0,strange3)*exp(1.0+strange4*xa0); //Funcsa=Funcsb
	    Funcba = (bottom0+bottom1*xa0+bottom2*pow(xa0,2))*pow((1-xa0),4)*pow(xa0,bottom3)*exp(1.0+bottom4*xa0); //Funcba=Funcbb
	    FuncABc = Funcsa*Funcsa/xa0/xa0;
	    FuncABs = Funcca*Funcca/xa0/xa0;
	    FuncABb = Funcba*Funcba/xa0/xa0;
	    Double_t totSec0 = 2*mVal*(
                                (FuncABu)*weightu
                                +(FuncABd)*weightd
                                +(FuncABc)*weightc
                                +(FuncABs)*weights
                                +(FuncABb)*weightb
                                );
	    totSec = 1.e-5*totSec0;
	    }
    
    if (totSec <= 0.) totSec = 0.00001;
	return totSec;
    
}




