#if !defined(__CINT__) || defined(__CLING__)

#include <iostream>
#include <vector>

#include "yaml-cpp/yaml.h"

#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TPythia8Decayer.h>
#include <TPaveStats.h>

#endif

using std::cout;
using std::endl;
using std::vector;

vector<TParticle*> SearchForDaughters(TClonesArray* array, vector<int> pdgCodeD);

void SimulateDecays(TString cfgfilename = "config_Bhadron_decay.yml")
{
    //upload config
    YAML::Node config = YAML::LoadFile(cfgfilename.Data());
    vector<int> pdgCodeB = config["pdgCodeB"].as<vector<int> >();
    vector<int> pdgCodeD = config["pdgCodeD"].as<vector<int> >();
    int nDecaysPerSpecies = config["nDecaysPerSpecies"].as<int>();

    //define decayer
    TPythia8Decayer *pdec = new TPythia8Decayer();
    pdec->Init();

    //define tree
    TTree* fTreeDecays = new TTree("fTreeDecays", "fTreeDecays");

    int pdgB = -9999;
    float ptB = -1.;
    float pB = -1.;
    float yB = -1.;
    vector<float> ptD;
    vector<float> pD;
    vector<float> yD;
    vector<int> pdgD;

    fTreeDecays->Branch("pdgB", &pdgB);
    fTreeDecays->Branch("ptB", &ptB);
    fTreeDecays->Branch("pB", &pB);
    fTreeDecays->Branch("yB", &yB);
    fTreeDecays->Branch("ptD", &ptD);
    fTreeDecays->Branch("pD", &pD);
    fTreeDecays->Branch("yD", &yD);
    fTreeDecays->Branch("pdgD", &pdgD);

    //decay simulation
    TRandom3 *gener = new TRandom3(0);
    TLorentzVector *vec = new TLorentzVector();
    TClonesArray *array = new TClonesArray("TParticle", 100);
    TDatabasePDG* db = TDatabasePDG::Instance();

    for (auto pdgCode: pdgCodeB)
    {
        for (auto iGen = 0; iGen < nDecaysPerSpecies; iGen++)
        {
            if (iGen % 100000 == 0)
                cout << "Generation number " << iGen << endl;

            float massB0 = db->GetParticle(pdgCode)->Mass();
            vec->SetPxPyPzE(0, 0, 0, massB0); //decay at rest
            pdec->Decay(pdgCode, vec);
            array->Clear();

            int nentries = pdec->ImportParticles(array);
            TParticle *Bmes = dynamic_cast<TParticle *>(array->At(0));
            pdgB = Bmes->GetPdgCode();
            ptB = Bmes->Pt();
            pB = Bmes->P();
            yB = Bmes->Y();
            vector<TParticle*> dautokeep = SearchForDaughters(array, pdgCodeD);
            for(unsigned int idau=0; idau<dautokeep.size(); idau++)
            {
                int pdgCodeDau = dautokeep[idau]->GetPdgCode();
                ptD.push_back(dautokeep[idau]->Pt());
                pD.push_back(dautokeep[idau]->P());
                yD.push_back(dautokeep[idau]->Y());
                pdgD.push_back(pdgCodeDau);
            }
            if(ptD.size() == 0)
            {
                ptD.push_back(-1);
                pD.push_back(-1);
                yD.push_back(-1);
                pdgD.push_back(-1);
            }

            fTreeDecays->Fill();
            pdgB = -1;
            ptD.clear();
            pD.clear();
            yD.clear();
            pdgD.clear();
        }
    }

    TFile outfile("Bhadrons_diffBRtoD.root", "recreate");
    fTreeDecays->Write();
    outfile.Close();
}

//________________________________________________________________________
vector<TParticle*> SearchForDaughters(TClonesArray* array, vector<int> pdgCodeD)
{
    vector<TParticle*> parttokeep;
    for (int idau = 1; idau < array->GetEntriesFast(); idau++) //the first one is the B
    {
        TParticle* dau = dynamic_cast<TParticle *>(array->At(idau));
        if(dau == nullptr)
            continue;
        int pdgDau = TMath::Abs(dau->GetPdgCode());
        if(std::find(pdgCodeD.begin(), pdgCodeD.end(), pdgDau) != pdgCodeD.end())
        {
            parttokeep.push_back(dau);
        }
    }

    return parttokeep;
}