using namespace std;
#include "stdio.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <TChain.h>
#include "TLeaf.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TObjArray.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TROOT.h"
#include <thread>
#include "TStopwatch.h"
#include <ROOT/TThreadExecutor.hxx>
#include "TLorentzVector.h"
#include <algorithm>

const float PI = TMath::Pi();
const float pt_trig_up = 2;
const float pt_trig_lo = .15;
const float EtaCut = 1.0;
const float Eta_ESE_Cut = 0.5;
const float cenDef[9] = {12.86, 12.02, 11.13, 10.17, 9.1, 7.89, 6.44, 4.56, 3.23}; //{13.14,12.29,11.38,10.39,9.29,8.05,6.57,4.65,3.29};   // 200 GeV Au+Au
const float v2_averaged[9] = {0.0711383, 0.0549932, 0.0539245, 0.0526822, 0.0517047, 0.0462665, 0.0357129, 0.0276133, 0.0227029};

const Int_t total_jobs = 10;
const Int_t total_min = 60;

void Calculate_Gamma_lamp(int cen = 0, int job = 9, int min = 92)
{ // main_function
        // Each centrality has 10 jobs, each job has 60 mins
        TStopwatch timer; // Times the current min
        timer.Start();
        TStopwatch timer2; // Times every 10000 events
        timer2.Start();
        cout << job << "\t" << min << endl;

        TChain *chain = new TChain("tree");
        int nfile = 0;
        char fname_in1[200];
        // Wherever your data is stored
        sprintf(fname_in1, "/media/Disk_YIN/AMPT-27GeV/v27-9977/*%d%02d1.root", job, min);

        // nfile counts number of files ending in "[job][min].root"
        nfile += chain->Add(fname_in1);

        cout << "Added " << nfile << " files" << endl;
        cout << "# entries in chain: " << chain->GetEntries() << endl;

        char fname_out[200];
        sprintf(fname_out, "cen%d.AMPT_Correlation_job%d%02d.root", cen, job, min);
        TFile fout(fname_out, "NEW");
        if (fout.IsZombie() || !fout.IsOpen())
        {
                std::cout << "File already exists: " << fname_out << std::endl;
                delete chain;
                return;
        }

        std::vector<float> *px_vec = nullptr;
        std::vector<float> *py_vec = nullptr;
        std::vector<float> *pz_vec = nullptr;
        std::vector<int> *pid_vec = nullptr;

        chain->SetBranchAddress("px", &px_vec);
        chain->SetBranchAddress("py", &py_vec);
        chain->SetBranchAddress("pz", &pz_vec);
        chain->SetBranchAddress("pid", &pid_vec);

        // defining variables
        Int_t Centrality, NPTracks;                //
        Float_t psi, b, Eta, Theta, Phi, Pt;       // track info
        Float_t psi2, b2, Eta2, Theta2, Phi2, Pt2; // track info
        Float_t Px, Py, Pz;
        Float_t Px2, Py2, Pz2;

        // defining histograms
        //  Anything you need add here
        TH1D *hCentrality = new TH1D("hCentrality", "hCentrality", 10, 0, 10);
        TH1D *hRefMult = new TH1D("hRefMult", "hRefMult", 1000, -0.5, 999.5);
        TH1D *hMult = new TH1D("hMult", "hMult", 1000, -0.5, 999.5);
        TH2D *hQaQb2 = new TH2D("hQaQb2", "hQaQb2", 250, 0, 25, 250, 0, 25);
        TH2D *hQ2Qc2 = new TH2D("hQ2Qc2", "hQ2Qc2", 250, 0, 25, 250, 0, 25);
        TProfile *p_Qa2_Qb2 = new TProfile("p_Qa2_Qb2", "p_Qa2_Qb2", 250, 0, 25, 0, 25, "");
        TProfile *p_Q2_Qc2 = new TProfile("p_Q2_Qc2", "p_Q2_Qc2", 250, 0, 25, 0, 25, "");
        TProfile *p_RefMult_Q2 = new TProfile("p_RefMult_Q2", "p_RefMult_Q2", 250, 0, 25, 0, 1000, "");
        TProfile *p_RefMult_Qa2 = new TProfile("p_RefMult_Qa2", "p_RefMult_Qa2", 250, 0, 25, 0, 1000, "");
        TProfile *p_RefMult_Qb2 = new TProfile("p_RefMult_Qb2", "p_RefMult_Qb2", 250, 0, 25, 0, 1000, "");
        TProfile *p_RefMult_Qc2 = new TProfile("p_RefMult_Qc2", "p_RefMult_Qc2", 250, 0, 25, 0, 1000, "");
        TProfile *Hist_cos = new TProfile("Hist_cos", "Hist_cos", 3, 0.5, 3.5, -1, 1, "");
        TProfile *p_cos_Q2 = new TProfile("p_cos_Q2", "p_cos_Q2", 250, 0, 25, -1, 1, "");
        TProfile *p_cos_Qa2 = new TProfile("p_cos_Qa2", "p_cos_Qa2", 250, 0, 25, -1, 1, "");
        TProfile *p_cos_Qc2 = new TProfile("p_cos_Qc2", "p_cos_Qc2", 250, 0, 25, -1, 1, "");
        TProfile *p_cos_Q = new TProfile("p_cos_Q", "p_cos_Q", 250, 0, 25, -1, 1, "");
        TProfile *p_cos_Qa = new TProfile("p_cos_Qa", "p_cos_Qa", 250, 0, 25, -1, 1, "");
        TProfile *p_cos_Qc = new TProfile("p_cos_Qc", "p_cos_Qc", 250, 0, 25, -1, 1, "");
        TProfile *pTemp_v2 = new TProfile("pTemp_v2", "pTemp_v2", 2, 0.5, 2.5, -100, 100, "");

        TProfile *p_v2_Q = new TProfile("p_v2_Q", "p_v2_Q", 250, 0, 5, -100, 100, "");
        TProfile *p_v2_Q_obs = new TProfile("p_v2_Q_obs", "p_v2_Q_obs", 250, 0, 5, -100, 100, "");
        TProfile *p_v2_Q2 = new TProfile("p_v2_Q2", "p_v2_Q2", 250, 0, 5, -100, 100, "");
        TProfile *p_v2_Q2_obs = new TProfile("p_v2_Q2_obs", "p_v2_Q2_obs", 250, 0, 5, -100, 100, "");
        TProfile *p_v2_Qc = new TProfile("p_v2_Qc", "p_v2_Qc", 250, 0, 5, -100, 100, "");
        TProfile *p_v2_Qc_obs = new TProfile("p_v2_Qc_obs", "p_v2_Qc_obs", 250, 0, 5, -100, 100, "");
        TProfile *p_v2_Qc2 = new TProfile("p_v2_Qc2", "p_v2_Qc2", 250, 0, 5, -100, 100, "");
        TProfile *p_v2_Qc2_obs = new TProfile("p_v2_Qc2_obs", "p_v2_Qc2_obs", 250, 0, 5, -100, 100, "");
        TProfile *p_v2_Qa = new TProfile("p_v2_Qa", "p_v2_Qa", 250, 0, 5, -100, 100, "");
        TProfile *p_v2_Qa_obs = new TProfile("p_v2_Qa_obs", "p_v2_Qa_obs", 250, 0, 5, -100, 100, "");
        TProfile *p_v2_Qa2 = new TProfile("p_v2_Qa2", "p_v2_Qa2", 250, 0, 5, -100, 100, "");
        TProfile *p_v2_Qa2_obs = new TProfile("p_v2_Qa2_obs", "p_v2_Qa2_obs", 250, 0, 5, -100, 100, "");

        TProfile *Gamma112 = new TProfile("Gamma112", "Gamma112", 8, 0.5, 8.5, -100, 100, "");
        TProfile *Gamma132 = new TProfile("Gamma132", "Gamma132", 8, 0.5, 8.5, -100, 100, "");
        TProfile *Delta = new TProfile("Delta", "Delta", 4, 0.5, 4.5, -100, 100, "");
        TProfile *v2_lam_p = new TProfile("v2_lam_p", "v2_lam_p", 4, 0.5, 4.5, -100, 100, "");

        Int_t nentries = chain->GetEntries();
        // loop through events
        for (int i = 0; i < nentries; i++)
        { // loop through all events

                if ((i + 1) % 10000 == 0)
                {
                        cout << "Processing entry == " << i + 1 << " == out of " << nentries << " in job " << job << ", min " << min << ".";
                        timer2.Stop();
                        timer2.Print();
                        timer2.Start();
                }

                chain->GetEntry(i);

                TLeaf *leaf_b = chain->GetLeaf("imp");
                TLeaf *leaf_RefMult = chain->GetLeaf("refmult");
                TLeaf *leaf_Np_p = chain->GetLeaf("npp");
                TLeaf *leaf_Np_t = chain->GetLeaf("npt");

                hRefMult->Fill(leaf_RefMult->GetValue(0));
                int Np = leaf_Np_p->GetValue(0) + leaf_Np_t->GetValue(0);
                if (Np < 3)
                        continue;

                psi = 0;
                b = leaf_b->GetValue(0);
                NPTracks = px_vec->size();

                Centrality = -1;
                // for(int j=0;j<9;j++) if(b<cenDef[j]) Centrality = j+1;
                for (int j = 0; j < 9; j++)
                        if (b < cenDef[j])
                                Centrality = j;

                hCentrality->Fill(Centrality);
                //cout << "Centrality = " << Centrality << endl;
		//cout << "cen = " << cen << endl;
		if (Centrality != cen)
                        continue;

		//cout << "Centrality = " << Centrality << endl;
                //cout << "cen = " << cen << endl;

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // construct TPC EP
                // calculate 4 Q: Q0 [0,0.5], Qa-1/3 Qb-2/3[0,0.5] Qc[0.5,1]
                TVector2 mQ1, mQ2, mQ, mQ3, mQ4;                                                                                     //> 1.5 range used for EP
                Double_t mQx1 = 0., mQy1 = 0., mQx2 = 0., mQy2 = 0., mQx = 0., mQy = 0., mQx3 = 0., mQy3 = 0., mQx4 = 0., mQy4 = 0.; // EP
                int Ecount = 0, Wcount = 0, Fcount = 0, Scount = 0;
                float mQQx = 0., mQQy = 0., mQQx1 = 0., mQQy1 = 0., mQQx2 = 0., mQQy2 = 0., mQQx3 = 0., mQQy3 = 0.; // 4 Q regions
                int Qcount = 0, Qcount1 = 0, Qcount2 = 0, Qcount3 = 0;

                for (int trki = 0; trki < NPTracks; trki++){
                        Px = px_vec->at(trki);
                        Py = py_vec->at(trki);
                        Pz = pz_vec->at(trki);
                        Pt = sqrt(Px * Px + Py * Py);
                        float PID = pid_vec->at(trki);

                        //if ((abs(PID) == 3122) || (abs(PID) == 2212))
                         if ((abs(PID) == 3122) || (abs(PID) == 211))
			        continue;
                        if (Pt > pt_trig_up || Pt < pt_trig_lo)
                                continue;

                        Fcount++;
                }
                std::vector<int> iTrack;
                iTrack.clear();
                Scount = Fcount/2;
                // cout << "First Fcount = " << Fcount << endl;
                for (int q = 0; q < Fcount; q++)
                        iTrack.push_back(q);
                random_shuffle(iTrack.begin(), iTrack.end());
                // cout << "size of iTrack = " << iTrack.size() << endl;

                // loop through matched primary tracks (the POI)
                //  First loop for Q and EP
                Fcount = 0;
                for (int trki = 0; trki < NPTracks; trki++)
                {
                        Px = px_vec->at(trki);
                        Py = py_vec->at(trki);
                        Pz = pz_vec->at(trki);
                        float PID = pid_vec->at(trki);

                        Pt = sqrt(Px * Px + Py * Py);
                        Phi = atan2(Py, Px);
                        Theta = atan2(Pt, Pz);
                        Eta = -log(tan(Theta / 2.));

                        //if ((abs(PID) == 3122) || (abs(PID) == 2212))
                        if ((abs(PID) == 3122) || (abs(PID) == 211))
			        continue;
                        if (Pt > pt_trig_up || Pt < pt_trig_lo)
                                continue;

                        // calculate 4 Q
                        if (fabs(Eta) < Eta_ESE_Cut)
                        {
                                mQQx += cos(Phi * 2.);
                                mQQy += sin(Phi * 2.);
                                Qcount++;
                        }
                        else
                        {
                                mQQx3 += cos(Phi * 2.);
                                mQQy3 += sin(Phi * 2.);
                                Qcount3++;
                        }

                        if (fabs(Eta) < Eta_ESE_Cut)
                        {
                                if (trki < int(NPTracks / 3.))
                                {
                                        mQQx1 += cos(Phi * 2.);
                                        mQQy1 += sin(Phi * 2.);
                                        Qcount1++;
                                }
                                else
                                {
                                        mQQx2 += cos(Phi * 2.);
                                        mQQy2 += sin(Phi * 2.);
                                        Qcount2++;
                                }
                        }

                        // if(Eta > 1.5)
                        if (Eta > 0.1)
                        {
                                mQx1 += Pt * cos(Phi * 2.);
                                mQy1 += Pt * sin(Phi * 2.);
                                Ecount++;
                        }
                        else if (Eta < -0.1)
                        {
                                mQx2 += Pt * cos(Phi * 2.);
                                mQy2 += Pt * sin(Phi * 2.);
                                Wcount++;
                        }

                        mQx += Pt * cos(Phi * 2.);
                        mQy += Pt * sin(Phi * 2.);

                        if (iTrack[Fcount] < (Scount))
                        {
                                mQx3 += Pt * cos(Phi * 2.);
                                mQy3 += Pt * sin(Phi * 2.);
                        }
                        else
                        {
                                mQx4 += Pt * cos(Phi * 2.);
                                mQy4 += Pt * sin(Phi * 2.);
                        }

                        // cout << "Fcount = " << Fcount << endl;
                        Fcount++;

                } // Track

                if ((mQx1 == 0 || mQy1 == 0 || mQx2 == 0 || mQy2 == 0))
                        continue;

                mQ1.Set(mQx1, mQy1);
                mQ2.Set(mQx2, mQy2);
                mQ.Set(mQx, mQy);
                float TPC_EP_east_new = 0.5 * mQ1.Phi(); // Sub event plane psi_n
                float TPC_EP_west_new = 0.5 * mQ2.Phi();
                float TPC_EP_full_new = 0.5 * mQ.Phi();
                float cos_fb = cos(2. * TPC_EP_east_new - 2. * TPC_EP_west_new);
                Hist_cos->Fill(1, cos(2. * TPC_EP_east_new - 2. * TPC_EP_west_new)); // resolution value for this event. will do for all events then sqrt

                mQ3.Set(mQx3, mQy3);
                mQ4.Set(mQx4, mQy4);
                float TPC_EP_for_new = 0.5 * mQ3.Phi();
                float TPC_EP_bac_new = 0.5 * mQ4.Phi();
                float cos_ew = cos(2. * TPC_EP_for_new - 2. * TPC_EP_bac_new);
                Hist_cos->Fill(2, cos_ew);

                // Second loop

                for (int trki = 0; trki < NPTracks; trki++)
                {
                        Px = px_vec->at(trki);
                        Py = py_vec->at(trki);
                        Pz = pz_vec->at(trki);
                        float PID = pid_vec->at(trki);

                        Pt = sqrt(Px * Px + Py * Py);
                        Phi = atan2(Py, Px);
                        Theta = atan2(Pt, Pz);
                        Eta = -log(tan(Theta / 2.));

                        if (abs(PID) != 3122)
                                continue;
                        if (Pt > 2.0 || Pt < 0.5)
                                continue;
			
			// cout << "PID = " << PID << endl;

                        float v2a = cos(2. * Phi - 2. * psi) * 100;
                        float v2b1 = cos(2. * Phi - 2. * TPC_EP_east_new) * 100;
                        float v2b2 = cos(2. * Phi - 2. * TPC_EP_west_new) * 100;

                        v2_lam_p->Fill(1, v2a);
                        if (Eta > 0.1)
                                v2_lam_p->Fill(2, v2b2);
                        else if (Eta < -0.1)
                                v2_lam_p->Fill(2, v2b1);

                        for (int trkj = 0; trkj < NPTracks; trkj++)
                        {
                                Px2 = px_vec->at(trkj);
                                Py2 = py_vec->at(trkj);
                                Pz2 = pz_vec->at(trkj);
                                float PID2 = pid_vec->at(trkj);

                                Pt2 = sqrt(Px2 * Px2 + Py2 * Py2);
                                Phi2 = atan2(Py2, Px2);
                                Theta2 = atan2(Pt2, Pz2);
                                Eta2 = -log(tan(Theta2 / 2.));

                                //if (abs(PID2) != 2212)
                                if (abs(PID2) != 211)
				        continue;
                                if (Pt2 > 2.0 || Pt2 < 0.4)
                                        continue;

                                float v2a_p = cos(2. * Phi2 - 2. * psi) * 100;
                                float v2b1_p = cos(2. * Phi2 - 2. * TPC_EP_east_new) * 100;
                                float v2b2_p = cos(2. * Phi2 - 2. * TPC_EP_west_new) * 100;

                                v2_lam_p->Fill(3, v2a_p);
                                if (Eta2 > 0.1)
                                        v2_lam_p->Fill(4, v2b2_p);
                                else if (Eta2 < -0.1)
                                        v2_lam_p->Fill(4, v2b1_p);

                                // correlations
                                float delta = cos(Phi - Phi2) * 100;

                                float gamma112e = cos(Phi + Phi2 - 2 * TPC_EP_east_new) * 100;
                                float gamma112w = cos(Phi + Phi2 - 2 * TPC_EP_west_new) * 100;
                                float gamma112 = cos(Phi + Phi2 - 2 * TPC_EP_full_new) * 100;
                                float gamma132a = cos(Phi - 3 * Phi2 + 2 * TPC_EP_full_new) * 100;
                                float gamma132b = cos(Phi2 - 3 * Phi + 2 * TPC_EP_full_new) * 100;
                                float gamma132 = (gamma132a + gamma132b) / 2.0;

                                float gamma112_rp = cos(Phi + Phi2) * 100;
                                float gamma132a_rp = cos(Phi - 3 * Phi2) * 100;
                                float gamma132b_rp = cos(Phi2 - 3 * Phi) * 100;
                                float gamma132_rp = (gamma132a_rp + gamma132b_rp) / 2.0;

                                //if((PID == 3122) && (PID2 == 2212)){
                                if((PID == 3122) && (PID2 == 211)){
				        Gamma112->Fill(1, gamma112);
                                        Gamma132->Fill(1, gamma132);
                                        Delta->Fill(1, delta);
                                        Gamma112->Fill(3, gamma112);
                                        Gamma132->Fill(3, gamma132);
                                        Delta->Fill(3, delta);

                                        Gamma112->Fill(5, gamma112_rp);
                                        Gamma132->Fill(5, gamma132_rp);
                                        Gamma112->Fill(7, gamma112_rp);
                                        Gamma132->Fill(7, gamma132_rp);
                                }
				else if((PID == -3122) && (PID2 == -211)){
                                //else if((PID == -3122) && (PID2 == -2212)){
                                        Gamma112->Fill(2, gamma112);
                                        Gamma132->Fill(2, gamma132);
                                        Delta->Fill(2, delta);
                                        Gamma112->Fill(3, gamma112);
                                        Gamma132->Fill(3, gamma132);
                                        Delta->Fill(3, delta);

                                        Gamma112->Fill(6, gamma112_rp);
                                        Gamma132->Fill(6, gamma132_rp);
                                        Gamma112->Fill(7, gamma112_rp);
                                        Gamma132->Fill(7, gamma132_rp);
                                }
                                else{
                                        Gamma112->Fill(4, gamma112);
                                        Gamma132->Fill(4, gamma132);
                                        Delta->Fill(4, delta);

                                        Gamma112->Fill(8, gamma112_rp);
                                        Gamma132->Fill(8, gamma132_rp);
                                }

                        } // track 2
                }         // track 1

                pTemp_v2->Reset();

        } // Event

        timer.Stop();
        std::cout << "The time of min " << min << ": ";
        timer.Print();
        fout.Write();
        delete chain;
        return;
}

void file_looper(Int_t centrality_select = 0, Int_t job_select = 0)
{
        static Int_t min_select = 0;
        int centrality_select_tmp = centrality_select;

        // Loops through all the centralities
        // If you want a specfic centrality, change the parameter value and value in while loop condition
        // e.g. for centrality 5: make parameter = 5, change while loop to <= 5
        // for centralities 4, 5, and 6: make parameter = 4, change while loop to <= 6
        while (centrality_select <= centrality_select_tmp)
        {
                // Directory where the data will be stored
                char folder_name[200];
                sprintf(folder_name, "/home/gang/Brian/New_Cuts/new_cen%dresults_newCuts", centrality_select);

                // Command to create the directory if it doesn't exist
                char create_folder[200];
                sprintf(create_folder, "mkdir /home/gang/Brian/New_Cuts/new_cen%dresults_newCuts", centrality_select);

                gSystem->Exec(create_folder);
                gSystem->cd(folder_name);

                // Loop through all jobs
                while (job_select < 10) // should be total_jobs
                {
                        for (int min_select = 0; min_select < 100; min_select++)
                        {
                                Calculate_Gamma_lamp(centrality_select, job_select, min_select);
                        };
                        job_select++;
                }
                // Start over for next centrality
                job_select = 0;
                min_select = 0;
                gSystem->cd("/media/Disk_YIN/AMPT-27GeV/");
                centrality_select++;
        }
        return;
}
