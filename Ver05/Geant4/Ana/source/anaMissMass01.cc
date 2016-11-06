/*
  Event Selection 01
  D*- -> D0 pi- -> K+ pi- pi- analysis
*/

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <signal.h>

#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TNtuple.h>
#include <TRandom.h>
#include <TTree.h>
#include <TBranch.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TGraph.h>

#include "ThreeVector.hh"
#include "LorentzVector.hh"
#include "RootHelper.hh"

#include "Random/Randomize.h"

#include "Particle.hh"
#include "ParticleID.hh"

#define FILE 0

#define AllDetect 0

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);
const double PI = acos(-1.);

const double PionMass        = 0.1395701;
const double KaonMass        = 0.493677;
const double ProtonMass      = 0.93827200;
const double DMesonZMass     = 1.86480;
const double DMesonStarMMass = 2.010270;
const double LambdacPMass    = 2.286460;

#if !AllDetect
//RICH PID
const double KMom1   = 0.0;
const double PiMom1  = 0.0;
const double PrMom1  = 0.0;

//PID1 PID
const double KMom2   = 2.5;
const double PiMom2  = 3.5;
const double PrMom2  = 2.5;

//TOF PID
const double KMom3   = 2.0;
const double PiMom3  = 2.0;
const double PrMom3  = 2.5;

//PID2 PID
const double KMom4   = 2.5;
const double PiMom4  = 2.5;
const double PrMom4  = 2.5;

//ITOF PID
const double KMom5   = 2.0;
const double PiMom5  = 2.0;
const double PrMom5  = 2.5;

//PAD PID
const double KMom6   = 1.5;
const double PiMom6  = 1.5;
const double PrMom6  = 2.3;

//MF PID
const double KMom7   = 1.5;
const double PiMom7  = 1.5;
const double PrMom7  = 0.0;
#endif

#if AllDetect
//RICH PID
const double KMom1   = 0.0;
const double PiMom1  = 0.0;
const double PrMom1  = 0.0;

//PID1 PID
const double KMom2   = 20.0;
const double PiMom2  = 20.0;
const double PrMom2  = 20.0;

//TOF PID
const double KMom3   = 20.0;
const double PiMom3  = 20.0;
const double PrMom3  = 20.0;

//PID2 PID
const double KMom4   = 20.0;
const double PiMom4  = 20.0;
const double PrMom4  = 20.0;

//ITOF PID
const double KMom5   = 20.0;
const double PiMom5  = 20.0;
const double PrMom5  = 20.0;

//PAD PID
const double KMom6   = 20.0;
const double PiMom6  = 20.0;
const double PrMom6  = 20.0;

//MF PID
const double KMom7   = 0.0;
const double PiMom7  = 0.0;
const double PrMom7  = 0.0;
#endif

#if FILE
const double MinDMMass = 1.852;
const double MaxDMMass = 1.878;

const double MinQV = 0.0043;
const double MaxQV = 0.0075;
#endif

#if !FILE
const double MinDMMass =  0.000;
const double MaxDMMass = 10.000;

const double MinQV =-100.000;
const double MaxQV = 100.000;
#endif

const double ResolDM   = 5.00E-03;
const double ResolDSM  = 0.60E-03;
const double ResolDSM0 = 6.00E-03;
const double ResolQ    = 0.60E-03;
const double ResolMM1  = 1.00E-02;

#ifndef MaxHits 
#define MaxHits 30
#endif

////For Reading data
struct Event{
  //Primary
  double d_priposx, d_priposy, d_priposz;
  double d_prim1, d_prim2;
  double d_prip1, d_prip2, d_prip1cm, d_prip2cm;
  double d_pritheta1, d_pritheta2, d_pritheta1cm, d_pritheta2cm;
  double d_priphi1, d_priphi2, d_priphi1cm, d_priphi2cm;
  double d_pricos1, d_pricos2, d_pricos1cm, d_pricos2cm;
  double d_priu1, d_priu2, d_priu1cm, d_priu2cm;
  double d_priv1, d_priv2, d_priv1cm, d_priv2cm;
  double d_pribeta1, d_pribeta2;
  double d_pbeam, d_ubeam, d_vbeam;
  double d_prip3, d_priu3, d_priv3;

  //Beam
  int d_ntB;
  double d_pBx[MaxHits], d_pBy[MaxHits], d_pBz[MaxHits];

  //PID
  int d_ntKp, d_ntKm;
  int d_ntPip, d_ntPim;
  int d_ntP, d_ntPb;
  int d_ntMup, d_ntMum;
  int d_ntEp, d_ntEm;

  double d_pkpx[MaxHits],  d_pkpy[MaxHits],  d_pkpz[MaxHits];
  double d_pkmx[MaxHits],  d_pkmy[MaxHits],  d_pkmz[MaxHits];
  double d_ppipx[MaxHits], d_ppipy[MaxHits], d_ppipz[MaxHits];
  double d_ppimx[MaxHits], d_ppimy[MaxHits], d_ppimz[MaxHits];
  double d_ppx[MaxHits],   d_ppy[MaxHits],   d_ppz[MaxHits];
  double d_ppbx[MaxHits],  d_ppby[MaxHits],  d_ppbz[MaxHits];
  double d_pmupx[MaxHits], d_pmupy[MaxHits], d_pmupz[MaxHits];
  double d_pmumx[MaxHits], d_pmumy[MaxHits], d_pmumz[MaxHits];
  double d_pepx[MaxHits],  d_pepy[MaxHits],  d_pepz[MaxHits];
  double d_pemx[MaxHits],  d_pemy[MaxHits],  d_pemz[MaxHits];

  int d_flkp[MaxHits], d_flkm[MaxHits];
  int d_flpip[MaxHits], d_flpim[MaxHits];
  int d_flp[MaxHits], d_flpb[MaxHits];
  int d_flmup[MaxHits], d_flmum[MaxHits];
  int d_flep[MaxHits], d_flem[MaxHits];
};
static Event event;

//For making tree
struct Ana{  
  //Primary
  double priposx, priposy, priposz;
  double prim1, prim2;
  double prip1, prip2, prip1cm, prip2cm;
  double pritheta1, pritheta2, pritheta1cm, pritheta2cm;
  double priphi1, priphi2, priphi1cm, priphi2cm;
  double pricos1, pricos2, pricos1cm, pricos2cm;
  double priu1, priu2, priu1cm, priu2cm;
  double priv1, priv2, priv1cm, priv2cm;
  double pribeta1, pribeta2;
  double pbeam, ubeam, vbeam;

  //PID
  int ntKp,  ntKm;
  int ntPip, ntPim;
  int ntP,   ntPb;
  int ntMup, ntMum;
  int ntEp,  ntEm;
  std::vector<double> pkp,  ukp,  vkp,  thetakp,  phikp;
  std::vector<double> pkm,  ukm,  vkm,  thetakm,  phikm;
  std::vector<double> ppip, upip, vpip, thetapip, phipip;
  std::vector<double> ppim, upim, vpim, thetapim, phipim;
  std::vector<double> pp,   up,   vp,   thetap,   phip;
  std::vector<double> ppb,  upb,  vpb,  thetapb,  phipb;
  std::vector<double> pmup, umup, vmup, thetamup, phimup;
  std::vector<double> pmum, umum, vmum, thetamum, phimum;
  std::vector<double> pep,  uep,  vep,  thetaep,  phiep;
  std::vector<double> pem,  uem,  vem,  thetaem,  phiem;

  //PID flag
  std::vector<int> fkp, fkm, fpip, fpim, fp, fpb, fmup, fmum, fep, fem;

  //Reaction
  int ntB;
  std::vector<double> pB,  uB,  vB;

  //(K+, pi-, pis-)
  std::vector<double> pkpD, thetakpD,  phikpD;
  std::vector<double> ppimD, thetapimD, phipimD;
  std::vector<double> ppimDs, thetapimDs, phipimDs;
  std::vector<double> thetakpi,  thetadpi;
  std::vector<double> thetakpis, thetapipis;
  std::vector<double> pDz, uDz, vDz;
  std::vector<double> thetaDz, phiDz;
  std::vector<double> DMass, DSMass, DSSMass;
  std::vector<double> DSMass0, M12, M23, M31;
  std::vector<double> Q_Value;
  std::vector<double> rDMass, rDSMass, rDSSMass;
  std::vector<double> rDSMass0, rM12, rM23, rM31;
  std::vector<double> rQ_Value;
  std::vector<double> costcmkp, costcmpim;
  std::vector<double> phicmkp,  phicmpim;
  std::vector<double> costcmD, costcmpims;
  std::vector<double> phicmD,  phicmpims;

  //Pion+ check
  std::vector<double> M_PipPim, M_PipPims, M_PipPimPims;
  std::vector<double> p_PipPim, p_PipPims, p_PipPimPims;
  std::vector<double> M_KpPipPim, M_KpPipPims, M_KpPipPimPims;
  std::vector<double> p_KpPipPim, p_KpPipPims, p_KpPipPimPims;
  std::vector<double> costcmpippim, costcmpippims;

  //Proton check
  std::vector<double> M_PrPim, M_PrPims, M_PrPimPims;
  std::vector<double> p_PrPim, p_PrPims, p_PrPimPims;
  std::vector<double> M_KpPrPim, M_KpPrPims, M_KpPrPimPims;
  std::vector<double> p_KpPrPim, p_KpPrPims, p_KpPrPimPims;
  std::vector<double> costcmprpim, costcmprpims;

  int ntD;
  std::vector<double> pD,  uD,  vD;
  std::vector<double> theta, phi;
  std::vector<double> MissMass, rMissMass;
  std::vector<double> RecoilMom;
  std::vector<double> costDcm, costYccm;
  std::vector<double> phiDcm, phiYccm;

  //Decay
  int    ntYcdp;
  std::vector<double> pYcdp,  uYcdp,  vYcdp;
  std::vector<double> thetaYcdp;
  std::vector<double> ppipYcd, thetapipYcd, phipipYcd;
  std::vector<double> MissMassdp, rMissMassdp;
  std::vector<double> costYcdpcm, phiYcdpcm;
  std::vector<double> ppipd, thetapipd, phipipd;

  int    ntYcdm;
  std::vector<double> pYcdm,  uYcdm,  vYcdm;
  std::vector<double> thetaYcdm;
  std::vector<double> MissMassdm, rMissMassdm;
  std::vector<double> costYcdmcm, phiYcdmcm;
  std::vector<double> ppimd, thetapimd, phipimd;

  int    ntYcdpm;
  std::vector<double> pYcdpm,  uYcdpm,  vYcdpm;
  std::vector<double> thetaYcdpm;
  std::vector<double> MissMassdpm, rMissMassdpm;
  std::vector<double> costYcdpmcm, phiYcdpmcm;
  std::vector<double> ppipd2, thetapipd2, phipipd2;
  std::vector<double> ppimd2, thetapimd2, phipimd2;

  int    ntYcdpd;
  std::vector<double> pYcdpd,  uYcdpd,  vYcdpd;
  std::vector<double> thetaYcdpd;
  std::vector<double> MissMassdpd, rMissMassdpd;
  std::vector<double> costYcdpdcm, phiYcdpdcm;
  std::vector<double> ppd, thetapd, phipd;
};
static Ana ana;

void SetBranchAddress(TChain *chain);
void InitializeAna();
void DefineHistograms();

void closeFile( int sig ){
  if(gFile){
    gFile->Write();
    gFile->Close();
  }
}

////Main
int main(Int_t argc, char **argv)
{
  if( argc<3 ){
    std::cerr << "Usage: " << argv[0] << " in_rootfile out_tootfile" << std::endl;
    exit(-1);
  }

  TFile *fout = new TFile( argv[2], "recreate");
  DefineHistograms();

  TChain *chain = new TChain("tree");
  chain->Add( argv[1] );

  SetBranchAddress(chain);

  signal(SIGINT,closeFile);

  int n_entries = chain->GetEntries();
  std::cout << n_entries << std::endl;

  for( int i=0; i<n_entries; i++ ){
    TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
    InitializeAna();

    if(i%10000==0)
      std::cout<<i<< " : " << n_entries <<std::endl;
    chain->GetEntry(i);

    /////////////////////Reading data
    //Beam  
    std::vector <ThreeVector> BeamPCont, BeamXCont; 
    ThreeVector primPos;
    int ntB=event.d_ntB;
    for( int ib=0; ib<ntB; ib++ ){
      double pbx=event.d_pBx[ib], pby=event.d_pBy[ib], pbz=event.d_pBz[ib];
      double x=event.d_priposx, y=event.d_priposy, z=event.d_priposz; 
      ThreeVector Pos( x, y, z );
      ThreeVector bMom( pbx, pby, pbz );
      primPos=Pos;
      BeamPCont.push_back(bMom); BeamXCont.push_back(primPos);   
    }

    ////PID
    Particle *part = new Particle();
    int id=0;

    //Number of Particles
    int ntKp = event.d_ntKp, ntKm = event.d_ntKm;
    int ntPip= event.d_ntPip,ntPim= event.d_ntPim;
    int ntP  = event.d_ntP,  ntPb = event.d_ntPb;
    int ntMup= event.d_ntMup,ntMum= event.d_ntMum;
    int ntEp = event.d_ntEp, ntEm = event.d_ntEm;

    //PID

    //Kaon+
    for( int it=0; it<ntKp; it++ ){
      double px=event.d_pkpx[it], py=event.d_pkpy[it], pz=event.d_pkpz[it];
      int fdet=event.d_flkp[it];
      double p = sqrt(px*px+py*py+pz*pz);

      if( fdet==1 && p>KMom1 ){
	id++;
	part->SetId(id);part->SetPid(IdKaonP);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==2 && p<KMom2 ){
	id++;
	part->SetId(id);part->SetPid(IdKaonP);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==3 && p<KMom3 ){
	id++;
	part->SetId(id);part->SetPid(IdKaonP);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==4 && p<KMom4 ){
	id++;
	part->SetId(id);part->SetPid(IdKaonP);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==5 && p<KMom5 ){
	id++;
	part->SetId(id);part->SetPid(IdKaonP);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==6 && p<KMom6 ){
      	id++;
      	part->SetId(id);part->SetPid(IdKaonP);
      	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
      	part->SetDetFlag(fdet);
      }
      // if( fdet==7 && p>KMom7 ){
      // 	id++;
      // 	part->SetId(id);part->SetPid(IdKaonP);
      // 	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
      // 	part->SetDetFlag(fdet);
      // }
    }

    //Kaon-
    for( int it=0; it<ntKm; it++ ){
      double px=event.d_pkmx[it], py=event.d_pkmy[it], pz=event.d_pkmz[it];
      int fdet=event.d_flkm[it];
      double p = sqrt(px*px+py*py+pz*pz);

      if( fdet==1 && p>KMom1 ){
	id++;
	part->SetId(id);part->SetPid(IdKaonM);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==2 && p<KMom2 ){
	id++;
	part->SetId(id);part->SetPid(IdKaonM);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==3 && p<KMom3 ){
	id++;
	part->SetId(id);part->SetPid(IdKaonM);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==4 && p<KMom4 ){
	id++;
	part->SetId(id);part->SetPid(IdKaonM);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==5 && p<KMom5 ){
	id++;
	part->SetId(id);part->SetPid(IdKaonM);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==6 && p<KMom6 ){
      	id++;
      	part->SetId(id);part->SetPid(IdKaonM);
      	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
      	part->SetDetFlag(fdet);
      }
      // if( fdet==7 && p>KMom7 ){
      // 	id++;
      // 	part->SetId(id);part->SetPid(IdKaonM);
      // 	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
      // 	part->SetDetFlag(fdet);
      // }
    }

    //Pion+
    for( int it=0; it<ntPip; it++ ){
      double px=event.d_ppipx[it], py=event.d_ppipy[it], pz=event.d_ppipz[it];
      int fdet=event.d_flpip[it];
      double p = sqrt(px*px+py*py+pz*pz);

      if( fdet==1 && p>PiMom1 ){
	id++;
	part->SetId(id);part->SetPid(IdPionP);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==2 && p<PiMom2 ){
	id++;
	part->SetId(id);part->SetPid(IdPionP);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==3 && p<PiMom3 ){
	id++;
	part->SetId(id);part->SetPid(IdPionP);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==4 && p<PiMom4 ){
	id++;
	part->SetId(id);part->SetPid(IdPionP);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==5 && p<PiMom5 ){
	id++;
	part->SetId(id);part->SetPid(IdPionP);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==6 && p<PiMom6 ){
      	id++;
      	part->SetId(id);part->SetPid(IdPionP);
      	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
      	part->SetDetFlag(fdet);
      }
      // if( fdet==7 && p>PiMom7 ){
      // 	id++;
      // 	part->SetId(id);part->SetPid(IdPionP);
      // 	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
      // 	part->SetDetFlag(fdet);
      // }
    }

    //Pion-
    for( int it=0; it<ntPim; it++ ){
      double px=event.d_ppimx[it], py=event.d_ppimy[it], pz=event.d_ppimz[it];
      int fdet=event.d_flpim[it];
      double p = sqrt(px*px+py*py+pz*pz);

      if( fdet==1 && p>PiMom1 ){
	id++;
	part->SetId(id);part->SetPid(IdPionM);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==2 && p<PiMom2 ){
	id++;
	part->SetId(id);part->SetPid(IdPionM);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==3 && p<PiMom3 ){
	id++;
	part->SetId(id);part->SetPid(IdPionM);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==4 && p<PiMom4 ){
	id++;
	part->SetId(id);part->SetPid(IdPionM);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==5 && p<PiMom5 ){
	id++;
	part->SetId(id);part->SetPid(IdPionM);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==6 && p<PiMom6 ){
      	id++;
      	part->SetId(id);part->SetPid(IdPionM);
      	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
      	part->SetDetFlag(fdet);
      }
      // if( fdet==7 && p>PiMom7 ){
      // 	id++;
      // 	part->SetId(id);part->SetPid(IdPionM);
      // 	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
      // 	part->SetDetFlag(fdet);
      // }
    }

    //Proton
    for( int it=0; it<ntP; it++ ){
      double px=event.d_ppx[it], py=event.d_ppy[it], pz=event.d_ppz[it];
      int fdet=event.d_flp[it];
      double p = sqrt(px*px+py*py+pz*pz);

      if( fdet==1 && p>PrMom1 ){
	id++;
	part->SetId(id);part->SetPid(IdProton);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==2 && p<PrMom2 ){
	id++;
	part->SetId(id);part->SetPid(IdProton);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==3 && p<PrMom3 ){
	id++;
	part->SetId(id);part->SetPid(IdProton);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==4 && p<PrMom4 ){
	id++;
	part->SetId(id);part->SetPid(IdProton);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==5 && p<PrMom5 ){
	id++;
	part->SetId(id);part->SetPid(IdProton);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==6 && p<PrMom6 ){
      	id++;
      	part->SetId(id);part->SetPid(IdProton);
      	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
      	part->SetDetFlag(fdet);
      }
      // if( fdet==7 && p>PrMom7 ){
      // 	id++;
      // 	part->SetId(id);part->SetPid(IdProton);
      // 	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
      // 	part->SetDetFlag(fdet);
      // }
    }

    //Protonb
    for( int it=0; it<ntPb; it++ ){
      double px=event.d_ppbx[it], py=event.d_ppby[it], pz=event.d_ppbz[it];
      int fdet=event.d_flpb[it];
      double p = sqrt(px*px+py*py+pz*pz);

      if( fdet==1 && p>PrMom1 ){
	id++;
	part->SetId(id);part->SetPid(IdProtonb);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==2 && p<PrMom2 ){
	id++;
	part->SetId(id);part->SetPid(IdProtonb);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==3 && p<PrMom3 ){
	id++;
	part->SetId(id);part->SetPid(IdProtonb);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==4 && p<PrMom4 ){
	id++;
	part->SetId(id);part->SetPid(IdProtonb);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==5 && p<PrMom5 ){
	id++;
	part->SetId(id);part->SetPid(IdProtonb);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==6 && p<PrMom6 ){
      	id++;
      	part->SetId(id);part->SetPid(IdProtonb);
      	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
      	part->SetDetFlag(fdet);
      }
      // if( fdet==7 && p>PrMom7 ){
      // 	id++;
      // 	part->SetId(id);part->SetPid(IdProtonb);
      // 	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
      // 	part->SetDetFlag(fdet);
      // }
    }

    //Muon+
    for( int it=0; it<ntMup; it++ ){
      double px=event.d_pmupx[it], py=event.d_pmupy[it], pz=event.d_pmupz[it];
      int fdet=event.d_flmup[it];
      double p = sqrt(px*px+py*py+pz*pz);
      
      if( fdet==1 && p>PiMom1 ){
	id++;
	part->SetId(id);part->SetPid(IdMuonP);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==2 && p<PiMom2 ){
	id++;
	part->SetId(id);part->SetPid(IdMuonP);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==3 && p<PiMom3 ){
	id++;
	part->SetId(id);part->SetPid(IdMuonP);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==4 && p<PiMom4 ){
	id++;
	part->SetId(id);part->SetPid(IdMuonP);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==5 && p<PiMom5 ){
	id++;
	part->SetId(id);part->SetPid(IdMuonP);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==6 && p<PiMom6 ){
      	id++;
      	part->SetId(id);part->SetPid(IdMuonP);
      	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
      	part->SetDetFlag(fdet);
      }
      // if( fdet==7 && p>PiMom7 ){
      // 	id++;
      // 	part->SetId(id);part->SetPid(IdMuonP);
      // 	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
      // 	part->SetDetFlag(fdet);
      // }
    }

    //Muon-
    for( int it=0; it<ntMum; it++ ){
      double px=event.d_pmumx[it], py=event.d_pmumy[it], pz=event.d_pmumz[it];
      int fdet=event.d_flmum[it]; 
      double p = sqrt(px*px+py*py+pz*pz);
      
      if( fdet==1 && p>PiMom1 ){
	id++;
	part->SetId(id);part->SetPid(IdMuonM);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==2 && p<PiMom2 ){
	id++;
	part->SetId(id);part->SetPid(IdMuonM);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==3 && p<PiMom3 ){
	id++;
	part->SetId(id);part->SetPid(IdMuonM);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==4 && p<PiMom4 ){
	id++;
	part->SetId(id);part->SetPid(IdMuonM);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==5 && p<PiMom5 ){
	id++;
	part->SetId(id);part->SetPid(IdMuonM);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==6 && p<PiMom6 ){
      	id++;
      	part->SetId(id);part->SetPid(IdMuonM);
      	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
      	part->SetDetFlag(fdet);
      }
      // if( fdet==7 && p>PiMom7 ){
      // 	id++;
      // 	part->SetId(id);part->SetPid(IdMuonM);
      // 	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
      // 	part->SetDetFlag(fdet);
      // }
    }

    //Electron
    for( int it=0; it<ntEp; it++ ){
      double px=event.d_pepx[it], py=event.d_pepy[it], pz=event.d_pepz[it];
      int fdet=event.d_flep[it];
      double p = sqrt(px*px+py*py+pz*pz);
      
      if( fdet==1 && p>PiMom1 ){
	id++;
	part->SetId(id);part->SetPid(IdPosi);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==2 && p<PiMom2 ){
	id++;
	part->SetId(id);part->SetPid(IdPosi);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==3 && p<PiMom3 ){
	id++;
	part->SetId(id);part->SetPid(IdPosi);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==4 && p<PiMom4 ){
	id++;
	part->SetId(id);part->SetPid(IdPosi);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==5 && p<PiMom5 ){
	id++;
	part->SetId(id);part->SetPid(IdPosi);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==6 && p<PiMom6 ){
      	id++;
      	part->SetId(id);part->SetPid(IdPosi);
      	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
      	part->SetDetFlag(fdet);
      }
      // if( fdet==7 && p>PiMom7 ){
      // 	id++;
      // 	part->SetId(id);part->SetPid(IdPosi);
      // 	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
      // 	part->SetDetFlag(fdet);
      // }
    }

    //Positron
    for( int it=0; it<ntEm; it++ ){
      double px=event.d_pemx[it], py=event.d_pemy[it], pz=event.d_pemz[it];
      int fdet=event.d_flem[it];
      double p = sqrt(px*px+py*py+pz*pz);
      
      if( fdet==1 && p>PiMom1 ){
	id++;
	part->SetId(id);part->SetPid(IdElec);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==2 && p<PiMom2 ){
	id++;
	part->SetId(id);part->SetPid(IdElec);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==3 && p<PiMom3 ){
	id++;
	part->SetId(id);part->SetPid(IdElec);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==4 && p<PiMom4 ){
	id++;
	part->SetId(id);part->SetPid(IdElec);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==5 && p<PiMom5 ){
	id++;
	part->SetId(id);part->SetPid(IdElec);
	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
	part->SetDetFlag(fdet);
      }
      if( fdet==6 && p<PiMom6 ){
      	id++;
      	part->SetId(id);part->SetPid(IdElec);
      	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
      	part->SetDetFlag(fdet);
      }
      // if( fdet==7 && p>PiMom7 ){
      // 	id++;
      // 	part->SetId(id);part->SetPid(IdElec);
      // 	part->SetMomX(px);part->SetMomY(py);part->SetMomZ(pz);
      // 	part->SetDetFlag(fdet);
      // }
    }

    //////////////////////////////////////
    /////////////////////Setting data

    //Primary Information
    ana.priposx = event.d_priposx, ana.priposy = event.d_priposy, ana.priposx = event.d_priposz; 
    ana.prim1 = event.d_prim1, ana.prim2 = event.d_prim2; 
    ana.prip1 = event.d_prip1, ana.prip2 = event.d_prip2; 
    ana.prip1cm = event.d_prip1cm, ana.prip2cm = event.d_prip2cm; 
    ana.pritheta1 = event.d_pritheta1, ana.pritheta2 = event.d_pritheta2; 
    ana.pritheta1cm = event.d_pritheta1cm, ana.pritheta2cm = event.d_pritheta2cm; 
    ana.priphi1 = event.d_priphi1, ana.priphi2 = event.d_priphi2; 
    ana.priphi1cm = event.d_priphi1cm, ana.priphi2cm = event.d_priphi2cm; 
    ana.pricos1 = event.d_pricos1, ana.pricos2 = event.d_pricos2; 
    ana.pricos1cm = event.d_pricos1cm, ana.pricos2cm = event.d_pricos2cm; 
    ana.priu1 = event.d_priu1, ana.priu2 = event.d_priu2; 
    ana.priu1cm = event.d_priu1cm, ana.priu2cm = event.d_priu2cm; 
    ana.priv1 = event.d_priv1, ana.priv2 = event.d_priv2; 
    ana.priv1cm = event.d_priv1cm, ana.priv2cm = event.d_priv2cm; 
    ana.pribeta1 = event.d_pribeta1, ana.pribeta2 = event.d_pribeta2; 
    ana.pbeam = event.d_pbeam, ana.ubeam = event.d_ubeam, ana.vbeam = event.d_vbeam; 

    //Beam
    int nB=BeamPCont.size();
    ana.ntB = nB;
    for( int ib=0; ib<nB; ib++ ){
      ThreeVector pb=BeamPCont[ib];
      double ub=pb.x()/pb.z(), vb=pb.y()/pb.z();
      ana.pB.push_back(pb.mag());
      ana.uB.push_back(ub);ana.vB.push_back(vb);
    }

    //Number of Particles
    int nKp=0, nKm=0, nPip=0, nPim=0, nP=0, nPb=0;
    int nMup=0, nMum=0, nEp=0, nEm=0;

    int pnum =part->GetSize();

    for( int i=0; i<pnum; i++ ){
      int pid=part->GetPid(i);
      double px=part->GetMomX(i), py=part->GetMomY(i), pz=part->GetMomZ(i);
      ThreeVector mom(px, py, pz);
      double p=mom.mag(), u=mom.x()/p, v=mom.y()/p;
      double theta = mom.theta()*Rad2Deg, phi=mom.phi()*Rad2Deg;
      int fdet=part->GetDetFlag(i);

      //Kaon+
      if( pid == IdKaonP ){
      	ana.pkp.push_back(p);ana.ukp.push_back(u);ana.vkp.push_back(v);
      	ana.thetakp.push_back(theta);ana.phikp.push_back(phi);
	ana.fkp.push_back(fdet);
	nKp++;
      }
      //Kaon-
      if( pid == IdKaonM ){
	ana.pkm.push_back(p);ana.ukm.push_back(u);ana.vkm.push_back(v);
	ana.thetakm.push_back(theta);ana.phikm.push_back(phi);
	ana.fkm.push_back(fdet);
	nKm++;
      }
      //Pion+
      if( pid == IdPionP ){
	ana.ppip.push_back(p);ana.upip.push_back(u);ana.vpip.push_back(v);
	ana.thetapip.push_back(theta);ana.phipip.push_back(phi);
	ana.fpip.push_back(fdet);
	nPip++;
      }
      //Pion-
      if( pid == IdPionM ){
	ana.ppim.push_back(p);ana.upim.push_back(u);ana.vpim.push_back(v);
	ana.thetapim.push_back(theta);ana.phipim.push_back(phi); 
	ana.fpim.push_back(fdet);
	nPim++;
      }
      //Proton
      if( pid == IdProton ){
	ana.pp.push_back(p);ana.up.push_back(u);ana.vp.push_back(v);
	ana.thetap.push_back(theta);ana.phip.push_back(phi);
	ana.fp.push_back(fdet);
	nP++;
      }
      //Protonb
      if( pid == IdProtonb ){
	ana.ppb.push_back(p);ana.upb.push_back(u);ana.vpb.push_back(v);
	ana.thetapb.push_back(theta);ana.phipb.push_back(phi);
	ana.fpb.push_back(fdet);
	nPb++;
      }
      //Muon+
      if( pid == IdMuonP ){
	ana.pmup.push_back(p);ana.umup.push_back(u);ana.vmup.push_back(v);
	ana.thetamup.push_back(theta);ana.phimup.push_back(phi);
	ana.fmup.push_back(fdet);
	nMup++;
      }
      //Muon-
      if( pid == IdMuonM ){
	ana.pmum.push_back(p);ana.umum.push_back(u);ana.vmum.push_back(v);
	ana.thetamum.push_back(theta);ana.phimum.push_back(phi);
	ana.fmum.push_back(fdet);
	nMum++;
      }
      //Positron
      if( pid == IdPosi ){
	ana.pep.push_back(p);ana.uep.push_back(u);ana.vep.push_back(v);
	ana.thetaep.push_back(theta);ana.phiep.push_back(phi);
	ana.fep.push_back(fdet);
	nEp++;
      }
      //Electron
      if( pid == IdElec ){
	ana.pem.push_back(p);ana.uem.push_back(u);ana.vem.push_back(v);
	ana.thetaem.push_back(theta);ana.phiem.push_back(phi);
	ana.fem.push_back(fdet);
	nEm++;
      }
    }
  
    ana.ntKp  = nKp;  ana.ntKm  = nKm; 
    ana.ntPip = nPip; ana.ntPim = nPim; 
    ana.ntP   = nP;   ana.ntPb  = nPb; 
    ana.ntMup = nMup; ana.ntMum = nMum; 
    ana.ntEp  = nEp;  ana.ntEm  = nEm; 


    /////////////////////////////////////////////////////////////////////////////
    //std::cout<< "**********" << std::endl; 
    //D* Meson Mass
    int ntDM=0;
    bool fileflag = false; 
    for( int ik1=0; ik1<pnum; ik1++ ){
      int id1=part->GetId(ik1);
      int pid1=part->GetPid(ik1);
      if( pid1 != IdKaonP ) continue;
      double px1=part->GetMomX(ik1), py1=part->GetMomY(ik1), pz1=part->GetMomZ(ik1);
      ThreeVector pkp1(px1, py1, pz1);
      LorentzVector LvKP1( pkp1, sqrt(KaonMass*KaonMass+pkp1.mag2()));

      for( int ipim1=0; ipim1<pnum; ipim1++ ){
	int id2=part->GetId(ipim1);
	int pid2=part->GetPid(ipim1);
	if( pid2 != IdPionM ) continue;
	double px2=part->GetMomX(ipim1), py2=part->GetMomY(ipim1), pz2=part->GetMomZ(ipim1);
	ThreeVector ppim1(px2, py2, pz2);
	LorentzVector LvPiM1( ppim1, sqrt(PionMass*PionMass+ppim1.mag2()));

	double costkpi=ppim1*pkp1/(ppim1.mag()*pkp1.mag());
	ThreeVector DMpos = primPos;

	LorentzVector LvDM( pkp1+ppim1, LvKP1.e()+LvPiM1.e() );
	double DMMass = LvDM.mag();
	ThreeVector DMmom = LvDM.vect();
	
	double pDz=DMmom.mag();
	double uDz=DMmom.x()/DMmom.z(), vDz=DMmom.y()/DMmom.z();
	
	double costDz = DMmom.cosTheta();
	double thetaDz = acos(costDz)*Rad2Deg;
	double phiDz = DMmom.phi()*Rad2Deg;
	
	//D invariant angle cut
	ThreeVector betaD((LvKP1+LvPiM1).vect()/(LvKP1+LvPiM1).e());
	//std::cout<<"betaD= "<< beta <<std::endl;
	LvKP1.boost(-1.*betaD);
	LvPiM1.boost(-1.*betaD);
	
	ThreeVector pKPcm = LvKP1.vect();
	double costKPcm   = pKPcm.cosTheta();
	double uKPcm =pKPcm.x()/pKPcm.z(), vKPcm=pKPcm.y()/pKPcm.z();
	double phiKPcm    =  pKPcm.phi()*Rad2Deg;
	
	ThreeVector pPiMcm = LvPiM1.vect();
	double costPiMcm   = pPiMcm.cosTheta();
	double uPiMcm =pPiMcm.x()/pPiMcm.z(), vPiMcm=pPiMcm.y()/pPiMcm.z();
	double phiPiMcm   = pPiMcm.phi()*Rad2Deg;
	
	// std::cout<<"costKPcm= "<< costKPcm <<std::endl;
	// std::cout<<"costPimcm= "<< costPiMcm <<std::endl;
	
	LvKP1.boost(betaD);
	LvPiM1.boost(betaD);
	
	for( int ipim2=0; ipim2<pnum; ipim2++ ){
	  int id3=part->GetId(ipim2);
	  int pid3=part->GetPid(ipim2);
	  if( pid3 != IdPionM ) continue;
	  if( id3 == id2 ) continue;
	  double px3=part->GetMomX(ipim2), py3=part->GetMomY(ipim2), pz3=part->GetMomZ(ipim2);
	  ThreeVector ppim2(px3, py3, pz3);
	  LorentzVector LvPiM2( ppim2, sqrt(PionMass*PionMass+ppim2.mag2()));
	  
	  double costkpis=ppim2*pkp1/(ppim2.mag()*pkp1.mag());
	  double costpipis=ppim1*ppim2/(ppim1.mag()*ppim2.mag());
	  double costdpi=DMmom*ppim2/(DMmom.mag()*ppim2.mag());
	  ThreeVector DSMpos = primPos;
	
	  LorentzVector LvDSM( DMmom+ppim2, sqrt(DMesonZMass*DMesonZMass+DMmom.mag2())+LvPiM2.e());
	  double DSMMass = LvDSM.mag();
	  ThreeVector DSMmom = LvDSM.vect();
	  
	  LorentzVector LvDSM0( pkp1+ppim1+ppim2, LvKP1.e()+LvPiM1.e()+LvPiM2.e());
	  double DSMMass0 = LvDSM0.mag();
	  ThreeVector DSMmom0 = LvDSM0.vect();	  
	  
	  double q_value = DSMMass0-DMMass-PionMass;
	  
	  LorentzVector LvM12( pkp1+ppim1,   LvKP1.e()+LvPiM1.e() );
	  LorentzVector LvM23( ppim1+ppim2, LvPiM1.e()+LvPiM2.e() );
	  LorentzVector LvM31( ppim2+pkp1, LvPiM2.e()+LvKP1.e() );
	  double M12 = LvM12.mag();
	  double M23 = LvM23.mag();
	  double M31 = LvM31.mag();

	  if( MinDMMass<DMMass && DMMass<MaxDMMass && MinQV<q_value && q_value<MaxQV ){
	    ntDM++;
	    // std::cout<< LvKP1.mag() << " " << LvKP1.vect() << std::endl;
	    // std::cout<< LvPiM1.mag() << " " << LvPiM1.vect() << std::endl;
	    // std::cout<< LvDM.mag()  << " " << LvDM.vect() << std::endl;
	    // std::cout<< LvPiM2.mag() << " " << LvPiM2.vect() << std::endl;
	    // std::cout<< LvDSM.mag() << " " << LvDSM.vect() << std::endl;
	    
	    // std::cout<<"******* D0 Meson Mass= "<< DMMass <<std::endl;
	    // std::cout<<"******* D* Meson Mass= "<< DSMMass <<std::endl;
	    	  
	    ana.pkpD.push_back(pkp1.mag());
	    ana.thetakpD.push_back(pkp1.theta()*Rad2Deg);
	    ana.phikpD.push_back(pkp1.phi()*Rad2Deg);
	    ana.ppimD.push_back(ppim1.mag());
	    ana.thetapimD.push_back(ppim1.theta()*Rad2Deg);
	    ana.phipimD.push_back(ppim1.phi()*Rad2Deg);
	    ana.ppimDs.push_back(ppim2.mag());
	    ana.thetapimDs.push_back(ppim2.theta()*Rad2Deg);
	    ana.phipimDs.push_back(ppim2.phi()*Rad2Deg);

	    ana.pDz.push_back(pDz);
	    ana.uDz.push_back(uDz);
	    ana.vDz.push_back(vDz);
	    ana.thetaDz.push_back(thetaDz);
	    ana.phiDz.push_back(phiDz);
	    ana.DMass.push_back(DMMass);
	    ana.rDMass.push_back( DMMass+CLHEP::RandGauss::shoot( 0.0, ResolDM ) );

	    ana.thetakpi.push_back(acos(costkpi)*Rad2Deg);
	    ana.thetadpi.push_back(acos(costdpi)*Rad2Deg);
	    ana.thetakpis.push_back(acos(costkpis)*Rad2Deg);
	    ana.thetapipis.push_back(acos(costpipis)*Rad2Deg);
	    
	    ana.DSMass.push_back(DSMMass);
	    ana.DSMass0.push_back(DSMMass0);
	    ana.M12.push_back(M12);
	    ana.M23.push_back(M23);
	    ana.M31.push_back(M31);
	    ana.Q_Value.push_back(q_value);
	    
	    ana.rDSMass.push_back( DSMMass+CLHEP::RandGauss::shoot( 0.0, ResolDSM ) );
	    ana.rDSMass0.push_back( DSMMass0+CLHEP::RandGauss::shoot( 0.0, ResolDSM0 ) );
	    ana.rM12.push_back( M12+CLHEP::RandGauss::shoot( 0.0, ResolDM ));
	    ana.rM23.push_back( M23+CLHEP::RandGauss::shoot( 0.0, ResolDM ));
	    ana.rM31.push_back( M31+CLHEP::RandGauss::shoot( 0.0, ResolDM ));
	    ana.rQ_Value.push_back( q_value+CLHEP::RandGauss::shoot( 0.0, ResolQ ) );
	    
	    ana.costcmkp.push_back(costKPcm);
	    ana.costcmpim.push_back(costPiMcm);
	    ana.phicmkp.push_back(phiKPcm);
	    ana.phicmpim.push_back(phiPiMcm);
	    
	    //D* invariant angle cut
	    ThreeVector betaDS( (LvDM+LvPiM2).vect()/(LvDM+LvPiM2).e());
	    //std::cout<<"betaDS= "<< beta <<std::endl;
	    LvDM.boost(-1.*betaDS);
	    LvPiM2.boost(-1.*betaDS);
	    
	    ThreeVector pDcm = LvDM.vect();
	    double costDcm   = pDcm.cosTheta();
	    double uDcm =pDcm.x()/pDcm.z(),   vDcm=pDcm.y()/pDcm.z();
	    double phiDcm    = pDcm.phi()*Rad2Deg;
	    
	    ThreeVector pPiMscm =  LvPiM2.vect();
	    double costPiMscm   = pPiMscm.cosTheta();
	    double uPiMscm=pPiMscm.x()/pPiMscm.z(), vPiMscm=pPiMscm.y()/pPiMscm.z();
	    double phiPiMscm   = pPiMscm.phi()*Rad2Deg;
	    
	    LvDM.boost(betaDS);
	    LvPiM2.boost(betaDS);
	    
	    ana.costcmD.push_back(costDcm);
	    ana.costcmpims.push_back(costPiMscm);
	    ana.phicmD.push_back(phiDcm);
	    ana.phicmpims.push_back(phiPiMscm);
	    
	    //Background check
	    for( int ipip=0; ipip<pnum; ipip++ ){
	      int id4=part->GetId(ipip);
	      int pid4=part->GetPid(ipip);
	      if( pid4 != IdPionP ) continue;
	      double px4=part->GetMomX(ipip), py4=part->GetMomY(ipip), pz4=part->GetMomZ(ipip);
	      ThreeVector ppip(px4, py4, pz4);
	      LorentzVector LvPip( ppip, sqrt(PionMass*PionMass+ppip.mag2()) );
	      
	      LorentzVector LvM_PipPim( ppip+ppim1, LvPip.e()+LvPiM1.e() );
	      LorentzVector LvM_PipPims( ppip+ppim2, LvPip.e()+LvPiM2.e() );
	      LorentzVector LvM_PipPimPims( ppip+ppim1+ppim2, LvPip.e()+LvPiM1.e()+LvPiM2.e() );
	      
	      ana.M_PipPim.push_back(LvM_PipPim.mag());
	      ana.M_PipPims.push_back(LvM_PipPims.mag());
	      ana.M_PipPimPims.push_back(LvM_PipPimPims.mag());
	      ana.p_PipPim.push_back(LvM_PipPim.vect().mag());
	      ana.p_PipPims.push_back(LvM_PipPims.vect().mag());
	      ana.p_PipPimPims.push_back(LvM_PipPimPims.vect().mag());
	      
	      LorentzVector LvM_KpPipPim( pkp1+ppip+ppim1, LvKP1.e()+LvPip.e()+LvPiM1.e() );
	      LorentzVector LvM_KpPipPims( pkp1+ppip+ppim2, LvKP1.e()+LvPip.e()+LvPiM2.e() );
	      LorentzVector LvM_KpPipPimPims( pkp1+ppip+ppim1+ppim2, LvKP1.e()+LvPip.e()+LvPiM1.e()+LvPiM2.e() );
	      
	      ana.M_KpPipPim.push_back(LvM_KpPipPim.mag());
	      ana.M_KpPipPims.push_back(LvM_KpPipPims.mag());
	      ana.M_KpPipPimPims.push_back(LvM_KpPipPimPims.mag());
	      ana.p_KpPipPim.push_back(LvM_KpPipPim.vect().mag());
	      ana.p_KpPipPims.push_back(LvM_KpPipPims.vect().mag());
	      ana.p_KpPipPimPims.push_back(LvM_KpPipPimPims.vect().mag());
	      
	      ThreeVector beta1((LvPip+LvPiM1).vect()/(LvPip+LvPiM1).e());
	      LvPip.boost(-1.*beta1);
	      ThreeVector pPipcm1 = LvPip.vect();
	      double costPipcm1   = pPipcm1.cosTheta();
	      ana.costcmpippim.push_back(costPipcm1);	  
	      LvPip.boost(beta1);
	      
	      ThreeVector beta2((LvPip+LvPiM2).vect()/(LvPip+LvPiM2).e());
	      LvPip.boost(-1.*beta2);
	      ThreeVector pPipcm2 = LvPip.vect();
	      double costPipcm2   = pPipcm2.cosTheta();
	      ana.costcmpippims.push_back(costPipcm2);	  
	      LvPip.boost(beta2);
	      
	      LorentzVector LvDSS( pkp1+ppim1+ppim2+ppip,
				   LvKP1.e()+LvPiM1.e()+LvPiM2.e()+LvPip.e() );
	      double DSSMMass = LvDSS.mag();
	      ana.DSSMass.push_back(DSSMMass);
	    }
	    for( int ipr=0; ipr<pnum; ipr++ ){
	      int id5=part->GetId(ipr);
	      int pid5=part->GetPid(ipr);
	      if( pid5 != IdProton ) continue;
	      double px5=part->GetMomX(ipr), py5=part->GetMomY(ipr), pz5=part->GetMomZ(ipr);
	      ThreeVector ppr(px5, py5, pz5);
	      LorentzVector LvPr( ppr, sqrt(ProtonMass*ProtonMass+ppr.mag2()) );
	      
	      LorentzVector LvM_PrPim( ppr+ppim1, LvPr.e()+LvPiM1.e() );
	      LorentzVector LvM_PrPims( ppr+ppim2, LvPr.e()+LvPiM2.e() );
	      LorentzVector LvM_PrPimPims( ppr+ppim1+ppim2, LvPr.e()+LvPiM1.e()+LvPiM2.e() );
	      
	      ana.M_PrPim.push_back(LvM_PrPim.mag());
	      ana.M_PrPims.push_back(LvM_PrPims.mag());
	      ana.M_PrPimPims.push_back(LvM_PrPimPims.mag());
	      ana.p_PrPim.push_back(LvM_PrPim.vect().mag());
	      ana.p_PrPims.push_back(LvM_PrPims.vect().mag());
	      ana.p_PrPimPims.push_back(LvM_PrPimPims.vect().mag());
	      
	      LorentzVector LvM_KpPrPim( pkp1+ppr+ppim1, LvKP1.e()+LvPr.e()+LvPiM1.e() );
	      LorentzVector LvM_KpPrPims( pkp1+ppr+ppim2, LvKP1.e()+LvPr.e()+LvPiM2.e() );
	      LorentzVector LvM_KpPrPimPims( pkp1+ppr+ppim1+ppim2, LvKP1.e()+LvPr.e()+LvPiM1.e()+LvPiM2.e() );
	      
	      ana.M_KpPrPim.push_back(LvM_KpPrPim.mag());
	      ana.M_KpPrPims.push_back(LvM_KpPrPims.mag());
	      ana.M_KpPrPimPims.push_back(LvM_KpPrPimPims.mag());
	      ana.p_KpPrPim.push_back(LvM_KpPrPim.vect().mag());
	      ana.p_KpPrPims.push_back(LvM_KpPrPims.vect().mag());
	      ana.p_KpPrPimPims.push_back(LvM_KpPrPimPims.vect().mag());
	      
	      ThreeVector beta1((LvPr+LvPiM1).vect()/(LvPr+LvPiM1).e());
	      LvPr.boost(-1.*beta1);
	      ThreeVector pPrcm1 = LvPr.vect();
	      double costPrcm1   = pPrcm1.cosTheta();
	      ana.costcmprpim.push_back(costPrcm1);	  
	      LvPr.boost(beta1);
	      
	      ThreeVector beta2((LvPr+LvPiM2).vect()/(LvPr+LvPiM2).e());
	      LvPr.boost(-1.*beta2);
	      ThreeVector pPrcm2 = LvPr.vect();
	      double costPrcm2   = pPrcm2.cosTheta();
	      ana.costcmprpims.push_back(costPrcm2);	  
	      LvPr.boost(beta2);
	    }
	    
	    /////////////////////////////////////////////////////////////////////////////
	    //Yc* Mass
	    for( int ib=0; ib<nB; ib++ ){
	      ThreeVector pb=BeamPCont[ib], xb=BeamXCont[ib];
	      ThreeVector pd=DSMmom, xd=primPos;
	      
	      double ud=pd.x()/pd.z(), vd=pd.y()/pd.z();
	      double cost=pb*pd/(pb.mag()*pd.mag());
	      
	      LorentzVector LvT( 0., 0., 0., ProtonMass );
	      LorentzVector LvD( pd, sqrt(DMesonStarMMass*DMesonStarMMass+pd.mag2()) );     
	      LorentzVector LvB( pb, sqrt(PionMass*PionMass+pb.mag2()) );
	      LorentzVector LvRc = LvB+LvT-LvD;
	      double MissMass = LvRc.mag();
	      double RecoilP = LvRc.vect().mag();
	      
	      ThreeVector betaI( (LvB+LvT).vect()/(LvB+LvT).e() );
	      LvD.boost(-1.*betaI);
	      LvRc.boost(-1.*betaI);
	      
	      ThreeVector pDcm  = LvD.vect();
	      ThreeVector pYccm = LvRc.vect();
	      double costDcm = pDcm.cosTheta();
	      double costYccm = pYccm.cosTheta();
	      double phiDcm = pDcm.phi()*Rad2Deg;
	      double phiYccm = pYccm.phi()*Rad2Deg;
	      
	      //std::cout<<"******* Missing Mass= " << MissMass <<std::endl;
	      
	      ana.pD.push_back(pd.mag());
	      ana.uD.push_back(ud);
	      ana.vD.push_back(vd);
	      ana.theta.push_back(acos(cost)*Rad2Deg);
	      ana.MissMass.push_back(MissMass);
	      ana.rMissMass.push_back( MissMass+CLHEP::RandGauss::shoot( 0.0, ResolMM1 ) );
	      ana.RecoilMom.push_back(RecoilP);
	      ana.costDcm.push_back(costDcm);
	      ana.costYccm.push_back(costYccm);
	      ana.phiDcm.push_back(phiDcm);
	      ana.phiYccm.push_back(phiYccm);
	      
	      LvD.boost(betaI);
	      LvRc.boost(betaI);
	      
#if FILE
	      fileflag = true; 
#endif	
	      
	      ////Decay Analysis
	      //Yc*+ -> Yc0 pi+
	      for( int ipipd=0; ipipd<pnum; ipipd++ ){
		int id6=part->GetId(ipipd);
		int pid6=part->GetPid(ipipd);
		if( pid6 != IdPionP ) continue;
		double px6=part->GetMomX(ipipd), py6=part->GetMomY(ipipd), pz6=part->GetMomZ(ipipd);
		ThreeVector ppipd(px6, py6, pz6);
		LorentzVector LvPiPd( ppipd, sqrt(PionMass*PionMass+ppipd.mag2()) );
		
		LorentzVector LvYcdp = LvRc-LvPiPd;
		double MissMassdp = LvYcdp.mag();
		double costYcdp = LvYcdp.vect().cosTheta();
		
		ThreeVector betaYcdp( LvRc.vect()/LvRc.e() );
		LvYcdp.boost(-1.*betaYcdp);
		ThreeVector pYcdpcm = LvYcdp.vect();
		double costYcdpcm   = pYcdpcm.cosTheta();
		double phiYcdpcm   = pYcdpcm.phi()*Rad2Deg;
		LvYcdp.boost(betaYcdp);
		
		// std::cout<<"******* Missing Mass Decay (Pi+) = "
		// 	 << MissMassdp <<std::endl;
		
		ana.pYcdp.push_back(LvYcdp.vect().mag());
		ana.uYcdp.push_back(LvYcdp.vect().x()/LvYcdp.vect().z());
		ana.vYcdp.push_back(LvYcdp.vect().y()/LvYcdp.vect().z());
		ana.thetaYcdp.push_back(acos(costYcdp)*Rad2Deg);

		ana.MissMassdp.push_back(MissMassdp); 
		ana.rMissMassdp.push_back( MissMassdp+CLHEP::RandGauss::shoot( 0.0, ResolMM1 ) );
		ana.costYcdpcm.push_back(costYcdpcm);
		ana.phiYcdpcm.push_back(phiYcdpcm);
		ana.ppipd.push_back(ppipd.mag());
		ana.thetapipd.push_back(ppipd.theta()*Rad2Deg);
		ana.phipipd.push_back(ppipd.phi()*Rad2Deg);
	      }
	      //Yc*+ -> Yc++ pi-
	      for( int ipimd=0; ipimd<pnum; ipimd++ ){
		int id7=part->GetId(ipimd);
		int pid7=part->GetPid(ipimd);
		if( pid7 != IdPionM ) continue;
		if( id7 == id2 || id7 == id3 ) continue; 
		double px7=part->GetMomX(ipimd), py7=part->GetMomY(ipimd), pz7=part->GetMomZ(ipimd);
		ThreeVector ppimd(px7, py7, pz7);
		LorentzVector LvPiMd( ppimd, sqrt(PionMass*PionMass+ppimd.mag2()) );
		
		LorentzVector LvYcdm = LvRc-LvPiMd;
		double MissMassdm = LvYcdm.mag();
		double costYcdm = LvYcdm.vect().cosTheta();
		
		ThreeVector betaYcdm( LvRc.vect()/LvRc.e() );
		LvYcdm.boost(-1.*betaYcdm);
		ThreeVector pYcdmcm = LvYcdm.vect();
		double costYcdmcm   = pYcdmcm.cosTheta();
		double phiYcdmcm   = pYcdmcm.phi()*Rad2Deg;
		LvYcdm.boost(betaYcdm);
		
		// std::cout<<"******* Missing Mass Decay (Pi-) = "
		// 	   << MissMassdm <<std::endl;
		
		ana.pYcdm.push_back(LvYcdm.vect().mag());
		ana.uYcdm.push_back(LvYcdm.vect().x()/LvYcdm.vect().z());
		ana.vYcdm.push_back(LvYcdm.vect().y()/LvYcdm.vect().z());
		ana.thetaYcdm.push_back(acos(costYcdm)*Rad2Deg);

		ana.MissMassdm.push_back(MissMassdm); 
		ana.rMissMassdm.push_back( MissMassdm+CLHEP::RandGauss::shoot( 0.0, ResolMM1 ) );
		ana.costYcdmcm.push_back(costYcdmcm);
		ana.phiYcdmcm.push_back(phiYcdmcm);
		ana.ppimd.push_back(ppimd.mag());
		ana.thetapimd.push_back(ppimd.theta()*Rad2Deg);
		ana.phipimd.push_back(ppimd.phi()*Rad2Deg);
	      }
	      //Yc*+ -> Lc+ pi+ pi-
	      for( int ipipd=0; ipipd<pnum; ipipd++ ){
		int id8=part->GetId(ipipd);
		int pid8=part->GetPid(ipipd);
		if( pid8 != IdPionP ) continue;
		double px8=part->GetMomX(ipipd), py8=part->GetMomY(ipipd), pz8=part->GetMomZ(ipipd);
		ThreeVector ppipd(px8, py8, pz8);
		LorentzVector LvPiPd( ppipd, sqrt(PionMass*PionMass+ppipd.mag2()) );
		
		for( int ipimd=0; ipimd<pnum; ipimd++ ){
		  int id9=part->GetId(ipimd);
		  int pid9=part->GetPid(ipimd);
		  if( pid9 != IdPionM ) continue; 
		  if( id9 == id2 || id9 == id3 ) continue;
		  double px9=part->GetMomX(ipimd), py9=part->GetMomY(ipimd), pz9=part->GetMomZ(ipimd);
		  ThreeVector ppimd(px9, py9, pz9);
		  LorentzVector LvPiMd( ppimd, sqrt(PionMass*PionMass+ppimd.mag2()) );
		  
		  LorentzVector LvYcdpm = LvRc-(LvPiPd+LvPiMd);
		  double MissMassdpm = LvYcdpm.mag();
		  double costYcdpm = LvYcdpm.vect().cosTheta();
		  
		  ThreeVector betaYcdpm( LvRc.vect()/LvRc.e() );
		  LvYcdpm.boost(-1.*betaYcdpm);
		  ThreeVector pYcdpmcm = LvYcdpm.vect();
		  double costYcdpmcm  = pYcdpmcm.cosTheta();
		  double phiYcdpmcm  = pYcdpmcm.phi()*Rad2Deg;
		  LvYcdpm.boost(betaYcdpm);
		  
		  // std::cout<<"******* Missing Mass Decay (Pi+-)= "
		  // 	       << MissMassdpm <<std::endl;
		  
		  ana.pYcdpm.push_back(LvYcdpm.vect().mag());
		  ana.uYcdpm.push_back(LvYcdpm.vect().x()/LvYcdpm.vect().z());
		  ana.vYcdpm.push_back(LvYcdpm.vect().y()/LvYcdpm.vect().z());
		  ana.thetaYcdpm.push_back(acos(costYcdpm)*Rad2Deg);

		  ana.MissMassdpm.push_back(MissMassdpm);
		  ana.rMissMassdpm.push_back( MissMassdpm+CLHEP::RandGauss::shoot( 0.0, ResolMM1 ) );
		  ana.costYcdpmcm.push_back(costYcdpmcm);
		  ana.phiYcdpmcm.push_back(phiYcdpmcm);
		  ana.ppipd2.push_back(ppipd.mag());
		  ana.thetapipd2.push_back(ppipd.theta()*Rad2Deg);
		  ana.phipipd2.push_back(ppipd.phi()*Rad2Deg);
		  ana.ppimd2.push_back(ppimd.mag());
		  ana.thetapimd2.push_back(ppimd.theta()*Rad2Deg);
		  ana.phipimd2.push_back(ppimd.phi()*Rad2Deg);
		}
	      }
	      //Yc*+ ->  p + D0
	      for( int ipr=0; ipr<pnum; ipr++ ){
		int id10=part->GetId(ipr);
		int pid10=part->GetPid(ipr);
		if( pid10 != IdProton ) continue;
		double px10=part->GetMomX(ipr), py10=part->GetMomY(ipr), pz10=part->GetMomZ(ipr);
		ThreeVector ppd(px10, py10, pz10);
		LorentzVector LvPd( ppd, sqrt(ProtonMass*ProtonMass+ppd.mag2()) );
		
		LorentzVector LvYcdpd = LvRc-LvPd;
		double MissMassdpd = LvYcdpd.mag();
		double costYcdpd = LvYcdpd.vect().cosTheta();
		
		ThreeVector betaYcdpd( LvRc.vect()/LvRc.e() );
		LvYcdpd.boost(-1.*betaYcdpd);
		ThreeVector pYcdpdcm = LvYcdpd.vect();
		double costYcdpdcm  = pYcdpdcm.cosTheta();
		double phiYcdpdcm  = pYcdpdcm.phi()*Rad2Deg;
		LvYcdpd.boost(betaYcdpd);
		
		// std::cout<<"******* Missing Mass Decay (P) = "
		// 	     << MissMassdpd <<std::endl;
		
		ana.pYcdpd.push_back(LvYcdpd.vect().mag());
		ana.uYcdpd.push_back(LvYcdpd.vect().x()/LvYcdpd.vect().z());
		ana.vYcdpd.push_back(LvYcdpd.vect().y()/LvYcdpd.vect().z());
		ana.thetaYcdpd.push_back(acos(costYcdpd)*Rad2Deg);

		ana.MissMassdpd.push_back(MissMassdpd);
		ana.rMissMassdpd.push_back( MissMassdpd+CLHEP::RandGauss::shoot( 0.0, ResolMM1 ) );
		ana.costYcdpdcm.push_back(costYcdpdcm);
		ana.phiYcdpdcm.push_back(phiYcdpdcm);
		ana.ppd.push_back(ppd.mag());
		ana.thetapd.push_back(ppd.theta()*Rad2Deg);
		ana.phipd.push_back(ppd.phi()*Rad2Deg);
	      }
	    }
	  }
	}
      }
    }
    ana.ntD = ntDM; 

#if FILE    
    if( fileflag ) tree->Fill();
#endif

#if !FILE
    tree->Fill();
#endif

    delete part;
  }
  ////End analysis

  fout->Write();
  fout->Close();
}

void SetBranchAddress( TChain *chain )
{
  //Primary
  chain->SetBranchAddress("priposx", &event.d_priposx); 
  chain->SetBranchAddress("priposy", &event.d_priposy); 
  chain->SetBranchAddress("priposz", &event.d_priposz); 
  chain->SetBranchAddress("prim1", &event.d_prim1); 
  chain->SetBranchAddress("prim2", &event.d_prim2); 
  chain->SetBranchAddress("prip1", &event.d_prip1); 
  chain->SetBranchAddress("prip2", &event.d_prip2); 
  chain->SetBranchAddress("prip1cm", &event.d_prip1cm); 
  chain->SetBranchAddress("prip2cm", &event.d_prip2cm); 
  chain->SetBranchAddress("pritheta1", &event.d_pritheta1); 
  chain->SetBranchAddress("pritheta2", &event.d_pritheta2); 
  chain->SetBranchAddress("pritheta1cm", &event.d_pritheta1cm); 
  chain->SetBranchAddress("pritheta2cm", &event.d_pritheta2cm); 
  chain->SetBranchAddress("priphi1", &event.d_priphi1); 
  chain->SetBranchAddress("priphi2", &event.d_priphi2); 
  chain->SetBranchAddress("priphi1cm", &event.d_priphi1cm); 
  chain->SetBranchAddress("priphi2cm", &event.d_priphi2cm); 
  chain->SetBranchAddress("pricos1", &event.d_pricos1); 
  chain->SetBranchAddress("pricos2", &event.d_pricos2); 
  chain->SetBranchAddress("pricos1cm", &event.d_pricos1cm); 
  chain->SetBranchAddress("pricos2cm", &event.d_pricos2cm); 
  chain->SetBranchAddress("priu1", &event.d_priu1); 
  chain->SetBranchAddress("priu2", &event.d_priu2); 
  chain->SetBranchAddress("priu1cm", &event.d_priu1cm); 
  chain->SetBranchAddress("priu2cm", &event.d_priu2cm); 
  chain->SetBranchAddress("priv1", &event.d_priv1); 
  chain->SetBranchAddress("priv2", &event.d_priv2); 
  chain->SetBranchAddress("priv1cm", &event.d_priv1cm); 
  chain->SetBranchAddress("priv2cm", &event.d_priv2cm); 
  chain->SetBranchAddress("pribeta1", &event.d_pribeta1); 
  chain->SetBranchAddress("pribeta2", &event.d_pribeta2); 
  chain->SetBranchAddress("pbeam", &event.d_pbeam); 
  chain->SetBranchAddress("ubeam", &event.d_ubeam); 
  chain->SetBranchAddress("vbeam", &event.d_vbeam); 
		 
  //PID		 
  chain->SetBranchAddress("ntKp", &event.d_ntKp);
  chain->SetBranchAddress("ntKm", &event.d_ntKm);
  chain->SetBranchAddress("ntPip", &event.d_ntPip);
  chain->SetBranchAddress("ntPim", &event.d_ntPim);
  chain->SetBranchAddress("ntP", &event.d_ntP);
  chain->SetBranchAddress("ntPb", &event.d_ntPb);
  chain->SetBranchAddress("ntMup", &event.d_ntMup);
  chain->SetBranchAddress("ntMum", &event.d_ntMum);
  chain->SetBranchAddress("ntEp", &event.d_ntEp);
  chain->SetBranchAddress("ntEm", &event.d_ntEm);

  chain->SetBranchAddress("pkpx", event.d_pkpx);
  chain->SetBranchAddress("pkpy", event.d_pkpy);
  chain->SetBranchAddress("pkpz", event.d_pkpz);
  chain->SetBranchAddress("pkmx", event.d_pkmx);
  chain->SetBranchAddress("pkmy", event.d_pkmy);
  chain->SetBranchAddress("pkmz", event.d_pkmz);
  chain->SetBranchAddress("ppipx", event.d_ppipx);
  chain->SetBranchAddress("ppipy", event.d_ppipy);
  chain->SetBranchAddress("ppipz", event.d_ppipz);
  chain->SetBranchAddress("ppimx", event.d_ppimx);
  chain->SetBranchAddress("ppimy", event.d_ppimy);
  chain->SetBranchAddress("ppimz", event.d_ppimz);
  chain->SetBranchAddress("ppx", event.d_ppx);
  chain->SetBranchAddress("ppy", event.d_ppy);
  chain->SetBranchAddress("ppz", event.d_ppz);
  chain->SetBranchAddress("ppbx", event.d_ppbx);
  chain->SetBranchAddress("ppby", event.d_ppby);
  chain->SetBranchAddress("ppbz", event.d_ppbz);
  chain->SetBranchAddress("pmupx", event.d_pmupx);
  chain->SetBranchAddress("pmupy", event.d_pmupy);
  chain->SetBranchAddress("pmupz", event.d_pmupz);
  chain->SetBranchAddress("pmumx", event.d_pmumx);
  chain->SetBranchAddress("pmumy", event.d_pmumy);
  chain->SetBranchAddress("pmumz", event.d_pmumz);
  chain->SetBranchAddress("pepx", event.d_pepx);
  chain->SetBranchAddress("pepy", event.d_pepy);
  chain->SetBranchAddress("pepz", event.d_pepz);
  chain->SetBranchAddress("pemx", event.d_pemx);
  chain->SetBranchAddress("pemy", event.d_pemy);
  chain->SetBranchAddress("pemz", event.d_pemz);

  //Beam
  chain->SetBranchAddress("ntB",&event.d_ntB);
  chain->SetBranchAddress("pBx",  event.d_pBx);
  chain->SetBranchAddress("pBy",  event.d_pBy);
  chain->SetBranchAddress("pBz",  event.d_pBz);

  //PID flag
  chain->SetBranchAddress("flkp", event.d_flkp);
  chain->SetBranchAddress("flkm", event.d_flkm);
  chain->SetBranchAddress("flpip", event.d_flpip);
  chain->SetBranchAddress("flpim", event.d_flpim);
  chain->SetBranchAddress("flp", event.d_flp);
  chain->SetBranchAddress("flpb", event.d_flpb);
  chain->SetBranchAddress("flmup", event.d_flmup);
  chain->SetBranchAddress("flmum", event.d_flmum);
  chain->SetBranchAddress("flep", event.d_flep);
  chain->SetBranchAddress("flem", event.d_flem);
}

//////Histograms
void InitializeAna( void )
{
  //Primary
  ana.priposx = -9999.0; 
  ana.priposy = -9999.0; 
  ana.priposz = -9999.0;
  ana.prim1 = -9999.0; 
  ana.prim2 = -9999.0;
  ana.prip1 = -9999.0; 
  ana.prip2 = -9999.0;
  ana.prip1cm = -9999.0; 
  ana.prip2cm = -9999.0;
  ana.pritheta1 = -9999.0; 
  ana.pritheta2 = -9999.0;
  ana.pritheta1cm = -9999.0; 
  ana.pritheta2cm = -9999.0;
  ana.priphi1 = -9999.0; 
  ana.priphi2 = -9999.0;
  ana.priphi1cm = -9999.0; 
  ana.priphi2cm = -9999.0;
  ana.pricos1 = -9999.0; 
  ana.pricos2 = -9999.0;
  ana.pricos1cm = -9999.0; 
  ana.pricos2cm = -9999.0;
  ana.priu1 = -9999.0; 
  ana.priu2 = -9999.0;
  ana.priu1cm = -9999.0; 
  ana.priu2cm = -9999.0;
  ana.priv1 = -9999.0; 
  ana.priv2 = -9999.0;
  ana.priv1cm = -9999.0; 
  ana.priv2cm = -9999.0;
  ana.pribeta1 = -9999.0; 
  ana.pribeta2 = -9999.0;
  ana.pbeam = -9999.0; 
  ana.ubeam = -9999.0;
  ana.vbeam = -9999.0;

  //PID
  ana.ntKp = -1; 
  ana.ntKm = -1;
  ana.ntPip = -1; 
  ana.ntPim = -1;
  ana.ntP = -1; 
  ana.ntPb = -1;
  ana.ntMup = -1; 
  ana.ntMum = -1;
  ana.ntEp = -1; 
  ana.ntEm = -1;

  ana.pkp.clear();
  ana.ukp.clear();
  ana.vkp.clear();
  ana.thetakp.clear();
  ana.phikp.clear();

  ana.pkm.clear();
  ana.ukm.clear();
  ana.vkm.clear();
  ana.thetakm.clear();
  ana.phikm.clear();

  ana.ppip.clear();
  ana.upip.clear();
  ana.vpip.clear();
  ana.thetapip.clear();
  ana.phipip.clear();

  ana.ppim.clear();
  ana.upim.clear();
  ana.vpim.clear();
  ana.thetapim.clear();
  ana.phipim.clear();

  ana.pp.clear();
  ana.up.clear();
  ana.vp.clear();
  ana.thetap.clear();
  ana.phip.clear();

  ana.ppb.clear();
  ana.upb.clear();
  ana.vpb.clear();
  ana.thetapb.clear();
  ana.phipb.clear();

  ana.pmup.clear();
  ana.umup.clear();
  ana.vmup.clear();
  ana.thetamup.clear();
  ana.phimup.clear();

  ana.pmum.clear();
  ana.umum.clear();
  ana.vmum.clear();
  ana.thetamum.clear();
  ana.phimum.clear();

  ana.pep.clear();
  ana.uep.clear();
  ana.vep.clear();
  ana.thetaep.clear();
  ana.phiep.clear();

  ana.pem.clear();
  ana.uem.clear();
  ana.vem.clear();
  ana.thetaem.clear();
  ana.phiem.clear();

  ana.fkp.clear();
  ana.fkm.clear();
  ana.fpip.clear();
  ana.fpim.clear();
  ana.fp.clear();
  ana.fpb.clear();
  ana.fmup.clear();
  ana.fmum.clear();
  ana.fep.clear();
  ana.fem.clear();

  //Reacrion
  ana.ntB = -1;
  ana.pB.clear();
  ana.uB.clear();
  ana.vB.clear();

  //(K+, pi-, pis-)
  ana.pkpD.clear();
  ana.thetakpD.clear();
  ana.phikpD.clear();
  ana.ppimD.clear();
  ana.thetapimD.clear();
  ana.phipimD.clear();
  ana.ppimDs.clear();
  ana.thetapimDs.clear();
  ana.phipimDs.clear();
  ana.thetakpi.clear();
  ana.thetadpi.clear();
  ana.thetakpis.clear();
  ana.thetapipis.clear();
  ana.pDz.clear();
  ana.uDz.clear();
  ana.vDz.clear();
  ana.thetaDz.clear();
  ana.phiDz.clear();
  ana.DMass.clear();
  ana.DSMass.clear();
  ana.DSSMass.clear();
  ana.DSMass0.clear();
  ana.M12.clear();
  ana.M23.clear();
  ana.M31.clear();
  ana.Q_Value.clear();
  ana.rDMass.clear();
  ana.rDSMass.clear();
  ana.rDSSMass.clear();
  ana.rDSMass0.clear();
  ana.rM12.clear();
  ana.rM23.clear();
  ana.rM31.clear();
  ana.rQ_Value.clear();
  ana.costcmkp.clear();
  ana.costcmpim.clear();
  ana.phicmkp.clear();
  ana.phicmpim.clear();
  ana.costcmD.clear();
  ana.costcmpims.clear();
  ana.phicmD.clear();
  ana.phicmpims.clear();

  ana.M_PipPim.clear();
  ana.M_PipPims.clear();
  ana.M_PipPimPims.clear();
  ana.p_PipPim.clear();
  ana.p_PipPims.clear();
  ana.p_PipPimPims.clear();
  ana.M_KpPipPim.clear();
  ana.M_KpPipPims.clear();
  ana.M_KpPipPimPims.clear();
  ana.p_KpPipPim.clear();
  ana.p_KpPipPims.clear();
  ana.p_KpPipPimPims.clear();
  ana.costcmpippim.clear();
  ana.costcmpippims.clear();

  ana.M_PrPim.clear();
  ana.M_PrPims.clear();
  ana.M_PrPimPims.clear();
  ana.p_PrPim.clear();
  ana.p_PrPims.clear();
  ana.p_PrPimPims.clear();
  ana.M_KpPrPim.clear();
  ana.M_KpPrPims.clear();
  ana.M_KpPrPimPims.clear();
  ana.p_KpPrPim.clear();
  ana.p_KpPrPims.clear();
  ana.p_KpPrPimPims.clear();
  ana.costcmprpim.clear();
  ana.costcmprpims.clear();

  ana.ntD = -1;
  ana.pD.clear();
  ana.uD.clear();
  ana.vD.clear();
  ana.theta.clear();
  ana.MissMass.clear();
  ana.rMissMass.clear();
  ana.RecoilMom.clear();
  ana.costDcm.clear();
  ana.costYccm.clear();
  ana.phiDcm.clear();
  ana.phiYccm.clear();

  //Decay
  ana.ntYcdp = -1;
  ana.ntYcdm = -1;
  ana.ntYcdpm = -1;
  ana.ntYcdpd = -1;
  ana.pYcdp.clear();
  ana.uYcdp.clear();
  ana.vYcdp.clear();
  ana.thetaYcdp.clear();
  ana.MissMassdp.clear();
  ana.rMissMassdp.clear();
  ana.costYcdpcm.clear();
  ana.phiYcdpcm.clear();
  ana.ppipd.clear();
  ana.thetapipd.clear();
  ana.phipipd.clear();

  ana.pYcdm.clear();
  ana.uYcdm.clear();
  ana.vYcdm.clear();
  ana.thetaYcdm.clear();
  ana.MissMassdm.clear();
  ana.rMissMassdm.clear();
  ana.costYcdmcm.clear();
  ana.phiYcdmcm.clear();
  ana.ppimd.clear();
  ana.thetapimd.clear();
  ana.phipimd.clear();

  ana.pYcdpm.clear();
  ana.uYcdpm.clear();
  ana.vYcdpm.clear();
  ana.thetaYcdpm.clear();
  ana.MissMassdpm.clear();
  ana.rMissMassdpm.clear();
  ana.costYcdpmcm.clear();
  ana.phiYcdpmcm.clear();
  ana.ppipd2.clear();
  ana.thetapipd2.clear();
  ana.phipipd2.clear();
  ana.ppimd2.clear();
  ana.thetapimd2.clear();
  ana.phipimd2.clear();

  ana.pYcdpd.clear();
  ana.uYcdpd.clear();
  ana.vYcdpd.clear();
  ana.thetaYcdpd.clear();
  ana.MissMassdpd.clear();
  ana.rMissMassdpd.clear();
  ana.costYcdpdcm.clear();
  ana.phiYcdpdcm.clear();
  ana.ppd.clear();
  ana.thetapd.clear();
  ana.phipd.clear();

  return;
}

void DefineHistograms()
{ 
  //Tree
  HBTree("tree","tree of Ana");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
   
  //Primary
  tree->Branch("priposx", &ana.priposx); 
  tree->Branch("priposy", &ana.priposy); 
  tree->Branch("priposz", &ana.priposz); 
  tree->Branch("prim1", &ana.prim1); 
  tree->Branch("prim2", &ana.prim2); 
  tree->Branch("prip1", &ana.prip1); 
  tree->Branch("prip2", &ana.prip2); 
  tree->Branch("prip1cm", &ana.prip1cm); 
  tree->Branch("prip2cm", &ana.prip2cm); 
  tree->Branch("pritheta1", &ana.pritheta1); 
  tree->Branch("pritheta2", &ana.pritheta2); 
  tree->Branch("pritheta1cm", &ana.pritheta1cm); 
  tree->Branch("pritheta2cm", &ana.pritheta2cm); 
  tree->Branch("priphi1", &ana.priphi1); 
  tree->Branch("priphi2", &ana.priphi2); 
  tree->Branch("priphi1cm", &ana.priphi1cm); 
  tree->Branch("priphi2cm", &ana.priphi2cm); 
  tree->Branch("pricos1", &ana.pricos1); 
  tree->Branch("pricos2", &ana.pricos2); 
  tree->Branch("pricos1cm", &ana.pricos1cm); 
  tree->Branch("pricos2cm", &ana.pricos2cm); 
  tree->Branch("priu1", &ana.priu1); 
  tree->Branch("priu2", &ana.priu2); 
  tree->Branch("priu1cm", &ana.priu1cm); 
  tree->Branch("priu2cm", &ana.priu2cm); 
  tree->Branch("priv1", &ana.priv1); 
  tree->Branch("priv2", &ana.priv2); 
  tree->Branch("priv1cm", &ana.priv1cm); 
  tree->Branch("priv2cm", &ana.priv2cm); 
  tree->Branch("pribeta1", &ana.pribeta1); 
  tree->Branch("pribeta2", &ana.pribeta2); 
  tree->Branch("pbeam", &ana.pbeam); 
  tree->Branch("ubeam", &ana.ubeam); 
  tree->Branch("vbeam", &ana.vbeam); 

  //PID
  tree->Branch("ntKp", &ana.ntKp);
  tree->Branch("ntKm", &ana.ntKm);
  tree->Branch("ntPip", &ana.ntPip);
  tree->Branch("ntPim", &ana.ntPim);
  tree->Branch("ntP", &ana.ntP);
  tree->Branch("ntPb", &ana.ntPb);
  tree->Branch("ntMup", &ana.ntMup);
  tree->Branch("ntMum", &ana.ntMum);
  tree->Branch("ntEp", &ana.ntEp);
  tree->Branch("ntEm", &ana.ntEm);

  tree->Branch("pkp", &ana.pkp);
  tree->Branch("ukp", &ana.ukp);
  tree->Branch("vkp", &ana.vkp);
  tree->Branch("thetakp", &ana.thetakp);
  tree->Branch("phikp", &ana.phikp);

  tree->Branch("pkm", &ana.pkm);
  tree->Branch("ukm", &ana.ukm);
  tree->Branch("vkm", &ana.vkm);
  tree->Branch("thetakm", &ana.thetakm);
  tree->Branch("phikm", &ana.phikm);

  tree->Branch("ppip", &ana.ppip);
  tree->Branch("upip", &ana.upip);
  tree->Branch("vpip", &ana.vpip);
  tree->Branch("thetapip", &ana.thetapip);
  tree->Branch("phipip", &ana.phipip);

  tree->Branch("ppim", &ana.ppim);
  tree->Branch("upim", &ana.upim);
  tree->Branch("vpim", &ana.vpim);
  tree->Branch("thetapim", &ana.thetapim);
  tree->Branch("phipim", &ana.phipim);

  tree->Branch("pp", &ana.pp);
  tree->Branch("up", &ana.up);
  tree->Branch("vp", &ana.vp);
  tree->Branch("thetap", &ana.thetap);
  tree->Branch("phip", &ana.phip);

  tree->Branch("ppb", &ana.ppb);
  tree->Branch("upb", &ana.upb);
  tree->Branch("vpb", &ana.vpb);
  tree->Branch("thetapb", &ana.thetapb);
  tree->Branch("phipb", &ana.phipb);

  tree->Branch("pmup", &ana.pmup);
  tree->Branch("umup", &ana.umup);
  tree->Branch("vmup", &ana.vmup);
  tree->Branch("thetamup", &ana.thetamup);
  tree->Branch("phimup", &ana.phimup);

  tree->Branch("pmum", &ana.pmum);
  tree->Branch("umum", &ana.umum);
  tree->Branch("vmum", &ana.vmum);
  tree->Branch("thetamum", &ana.thetamum);
  tree->Branch("phimum", &ana.phimum);

  tree->Branch("pep", &ana.pep);
  tree->Branch("uep", &ana.uep);
  tree->Branch("vep", &ana.vep);
  tree->Branch("thetaep", &ana.thetaep);
  tree->Branch("phiep", &ana.phiep);

  tree->Branch("pem", &ana.pem);
  tree->Branch("uem", &ana.uem);
  tree->Branch("vem", &ana.vem);
  tree->Branch("thetaem", &ana.thetaem);
  tree->Branch("phiem", &ana.phiem);

  //PID flag
  tree->Branch("fkp", &ana.fkp);
  tree->Branch("fkm", &ana.fkm);
  tree->Branch("fpip", &ana.fpip);
  tree->Branch("fpim", &ana.fpim);
  tree->Branch("fp", &ana.fp);
  tree->Branch("fpb", &ana.fpb);
  tree->Branch("fmup", &ana.fmup);
  tree->Branch("fmum", &ana.fmum);
  tree->Branch("fep", &ana.fep);
  tree->Branch("fem", &ana.fem);

  //Reaction
  tree->Branch("ntB", &ana.ntB);
  tree->Branch("pB", &ana.pB);
  tree->Branch("uB", &ana.uB);
  tree->Branch("vB", &ana.vB);

  //(K+, pi-, pis-)
  tree->Branch("pkpD", &ana.pkpD);
  tree->Branch("thetakpD", &ana.thetakpD);
  tree->Branch("phikpD", &ana.phikpD);
  tree->Branch("ppimD", &ana.ppimD);
  tree->Branch("thetapimD", &ana.thetapimD);
  tree->Branch("phipimD", &ana.phipimD);
  tree->Branch("ppimDs", &ana.ppimDs);
  tree->Branch("thetapimDs", &ana.thetapimDs);
  tree->Branch("phipimDs", &ana.phipimDs);
  tree->Branch("thetakpi", &ana.thetakpi);
  tree->Branch("thetadpi", &ana.thetadpi);
  tree->Branch("thetakpis", &ana.thetakpis);
  tree->Branch("thetapipis", &ana.thetapipis);
  tree->Branch("pDz", &ana.pDz);
  tree->Branch("uDz", &ana.uDz);
  tree->Branch("vDz", &ana.vDz);
  tree->Branch("thetaDz", &ana.thetaDz);
  tree->Branch("phiDz", &ana.phiDz);
  tree->Branch("DMass", &ana.DMass);
  tree->Branch("DSMass", &ana.DSMass);
  tree->Branch("DSSMass", &ana.DSSMass);
  tree->Branch("DSMass0", &ana.DSMass0);
  tree->Branch("M12", &ana.M12);
  tree->Branch("M23", &ana.M23);
  tree->Branch("M31", &ana.M31);
  tree->Branch("Q_Value", &ana.Q_Value);
  tree->Branch("rDMass", &ana.rDMass);
  tree->Branch("rDSMass", &ana.rDSMass);
  tree->Branch("rDSSMass", &ana.rDSSMass);
  tree->Branch("rDSMass0", &ana.rDSMass0);
  tree->Branch("rM12", &ana.rM12);
  tree->Branch("rM23", &ana.rM23);
  tree->Branch("rM31", &ana.rM31);
  tree->Branch("rQ_Value", &ana.rQ_Value);
  tree->Branch("costcmkp", &ana.costcmkp);
  tree->Branch("costcmpim", &ana.costcmpim);
  tree->Branch("phicmkp", &ana.phicmkp);
  tree->Branch("phicmpim", &ana.phicmpim);
  tree->Branch("costcmD", &ana.costcmD);
  tree->Branch("costcmpims", &ana.costcmpims);
  tree->Branch("phicmD", &ana.phicmD);
  tree->Branch("phicmpims", &ana.phicmpims);

  tree->Branch("M_PipPim", &ana.M_PipPim);
  tree->Branch("M_PipPims", &ana.M_PipPims);
  tree->Branch("M_PipPimPims", &ana.M_PipPimPims);
  tree->Branch("p_PipPim", &ana.p_PipPim);
  tree->Branch("p_PipPims", &ana.p_PipPims);
  tree->Branch("p_PipPimPims", &ana.p_PipPimPims);
  tree->Branch("M_KpPipPim", &ana.M_KpPipPim);
  tree->Branch("M_KpPipPims", &ana.M_KpPipPims);
  tree->Branch("M_KpPipPimPims", &ana.M_KpPipPimPims);
  tree->Branch("p_KpPipPim", &ana.p_KpPipPim);
  tree->Branch("p_KpPipPims", &ana.p_KpPipPims);
  tree->Branch("p_KpPipPimPims", &ana.p_KpPipPimPims);
  tree->Branch("costcmpippim", &ana.costcmpippim);
  tree->Branch("costcmpippims", &ana.costcmpippims);

  tree->Branch("M_PrPim", &ana.M_PrPim);
  tree->Branch("M_PrPims", &ana.M_PrPims);
  tree->Branch("M_PrPimPims", &ana.M_PrPimPims);
  tree->Branch("p_PrPim", &ana.p_PrPim);
  tree->Branch("p_PrPims", &ana.p_PrPims);
  tree->Branch("p_PrPimPims", &ana.p_PrPimPims);
  tree->Branch("M_KpPrPim", &ana.M_KpPrPim);
  tree->Branch("M_KpPrPims", &ana.M_KpPrPims);
  tree->Branch("M_KpPrPimPims", &ana.M_KpPrPimPims);
  tree->Branch("p_KpPrPim", &ana.p_KpPrPim);
  tree->Branch("p_KpPrPims", &ana.p_KpPrPims);
  tree->Branch("p_KpPrPimPims", &ana.p_KpPrPimPims);
  tree->Branch("costcmprpim", &ana.costcmprpim);
  tree->Branch("costcmprpims", &ana.costcmprpims);

  tree->Branch("ntD", &ana.ntD);
  tree->Branch("pD", &ana.pD);
  tree->Branch("uD", &ana.uD);
  tree->Branch("vD", &ana.vD);
  tree->Branch("theta", &ana.theta);
  tree->Branch("MissMass", &ana.MissMass);
  tree->Branch("rMissMass", &ana.rMissMass);
  tree->Branch("RecoilMom", &ana.RecoilMom);
  tree->Branch("costDcm", &ana.costDcm);
  tree->Branch("costYccm", &ana.costYccm);
  tree->Branch("phiDcm", &ana.phiDcm);
  tree->Branch("phiYccm", &ana.phiYccm);

  //Decay
  tree->Branch("ntYcdp", &ana.ntYcdp);
  tree->Branch("ntYcdm", &ana.ntYcdm);
  tree->Branch("ntYcdpm", &ana.ntYcdpm);
  tree->Branch("ntYcdpd", &ana.ntYcdpd);
  tree->Branch("pYcdp", &ana.pYcdp);
  tree->Branch("uYcdp", &ana.uYcdp);
  tree->Branch("vYcdp", &ana.vYcdp);
  tree->Branch("thetaYcdp", &ana.thetaYcdp);
  tree->Branch("MissMassdp", &ana.MissMassdp);
  tree->Branch("rMissMassdp", &ana.rMissMassdp);
  tree->Branch("costYcdpcm", &ana.costYcdpcm);
  tree->Branch("phiYcdpcm", &ana.phiYcdpcm);
  tree->Branch("ppipd", &ana.ppipd);
  tree->Branch("thetapipd", &ana.thetapipd);
  tree->Branch("phipipd", &ana.phipipd);

  tree->Branch("pYcdm", &ana.pYcdm);
  tree->Branch("uYcdm", &ana.uYcdm);
  tree->Branch("vYcdm", &ana.vYcdm);
  tree->Branch("thetaYcdm", &ana.thetaYcdm);
  tree->Branch("MissMassdm", &ana.MissMassdm);
  tree->Branch("rMissMassdm", &ana.rMissMassdm);
  tree->Branch("costYcdmcm", &ana.costYcdmcm);
  tree->Branch("phiYcdmcm", &ana.phiYcdmcm);
  tree->Branch("ppimd", &ana.ppimd);
  tree->Branch("thetapimd", &ana.thetapimd);
  tree->Branch("phipimd", &ana.phipimd);

  tree->Branch("pYcdpm", &ana.pYcdpm);
  tree->Branch("uYcdpm", &ana.uYcdpm);
  tree->Branch("vYcdpm", &ana.vYcdpm);
  tree->Branch("thetaYcdpm", &ana.thetaYcdpm);
  tree->Branch("MissMassdpm", &ana.MissMassdpm);
  tree->Branch("rMissMassdpm", &ana.rMissMassdpm);
  tree->Branch("costYcdpmcm", &ana.costYcdpmcm);
  tree->Branch("phiYcdpmcm", &ana.phiYcdpmcm);
  tree->Branch("ppipd2", &ana.ppipd2);
  tree->Branch("thetapipd2", &ana.thetapipd2);
  tree->Branch("phipipd2", &ana.phipipd2);
  tree->Branch("ppimd2", &ana.ppimd2);
  tree->Branch("thetapimd2", &ana.thetapimd2);
  tree->Branch("phipimd2", &ana.phipimd2);

  tree->Branch("pYcdpd", &ana.pYcdpd);
  tree->Branch("uYcdpd", &ana.uYcdpd);
  tree->Branch("vYcdpd", &ana.vYcdpd);
  tree->Branch("thetaYcdpd", &ana.thetaYcdpd);
  tree->Branch("MissMassdpd", &ana.MissMassdpd);
  tree->Branch("rMissMassdpd", &ana.rMissMassdpd);
  tree->Branch("costYcdpdcm", &ana.costYcdpdcm);
  tree->Branch("phiYcdpdcm", &ana.phiYcdpdcm);
  tree->Branch("ppd", &ana.ppd);
  tree->Branch("thetapd", &ana.thetapd);
  tree->Branch("phipd", &ana.phipd);
}
