#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <strstream>

#include <TH1F.h>
#include <TTree.h>
#include <TROOT.h>
#include <TFile.h>
#include <TSystem.h>
#include <TLorentzVector.h>
#include <math.h>
#include <TMatrixD.h>
#include <TMatrixDEigen.h>
#include <TTreeReader.h>

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/FWLite/interface/FWLiteEnabler.h"

#include "DataFormats/FWLite/interface/Run.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/FWLite/interface/LuminosityBlock.h"

#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoProduct.h"

#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"

#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "GeneratorInterface/LHEInterface/interface/LHERunInfo.h"

#include "DataFormats/Common/interface/MergeableCounter.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

const double eMass = 0.0005109989; //GeV
const double muMass = 0.1056583715; //GeV
bool debug = false;
bool mcDebug = false;

std::pair<int, int> getJets(edm::Handle<std::vector<pat::Jet>> selectedjets,
	std::string mode="CSV", std::string bTagAlgo="pfCombinedInclusiveSecondaryVertexV2BJetTags") {
	/*Selects pair of jets ordered by pT*/
	std::pair<int, int> pair;
	if (mode == "mass") { //Invariant mass closest to 125 GeV, 
		double deltaMin = -1;
		double delta;
		int iMin = -1;
		int jMin = -1;
		for(size_t i = 0; i < selectedjets->size(); ++i) {
			for(size_t j = 0; j < selectedjets->size(); ++j) {
				if (i == j) continue;
				delta = std::abs(125-(selectedjets->at(i).p4() + selectedjets->at(j).p4()).M());
				if (deltaMin > delta || deltaMin < 0) {
					deltaMin = delta;
					iMin = i;
					jMin = j;
				}
			}
		}
		pair.first = iMin;
		pair.second = jMin;
		if (selectedjets->at(pair.first).p4().Pt() < selectedjets->at(pair.second).p4().Pt()) {
			pair.second = iMin;
			pair.first = jMin;
		}
	} else if (mode == "pT") { //Highest pT
		std::pair<int, double> leading = std::pair<int, double>(-1,-1.0);
		std::pair<int, double> subLeading = std::pair<int, double>(-1,-1.0);
		for(size_t i = 0; i < selectedjets->size(); ++i) {
			if (selectedjets->at(i).p4().Pt() > leading.second) {
				subLeading.first = leading.first;
				subLeading.second = leading.second;
				leading.first = i;
				leading.second = selectedjets->at(i).p4().Pt();
			} else if (selectedjets->at(i).p4().Pt() > subLeading.second) {
				subLeading.first = i;
				subLeading.second = selectedjets->at(i).p4().Pt();
			}
		}
		pair.first = leading.first;
		pair.second = subLeading.first;
	}  else if (mode == "CSV") { //Two highest CSV, ordered by pT
		std::pair<int, double> leading = std::pair<int, double>(-1, -9999);
		std::pair<int, double> subLeading = std::pair<int, double>(-1, -9999);
		for(size_t i = 0; i < selectedjets->size(); ++i) {
			if (selectedjets->at(i).bDiscriminator(bTagAlgo) > leading.second) {
				subLeading.first = leading.first;
				subLeading.second = leading.second;
				leading.first = i;
				leading.second = selectedjets->at(i).bDiscriminator(bTagAlgo);
			} else if (selectedjets->at(i).bDiscriminator(bTagAlgo) > subLeading.second) {
				subLeading.first = i;
				subLeading.second = selectedjets->at(i).bDiscriminator(bTagAlgo);
			}
		}
		pair.first = leading.first;
		pair.second = subLeading.first;
		if (selectedjets->at(pair.first).p4().Pt() < selectedjets->at(pair.second).p4().Pt()) {
			pair.second = leading.first;
			pair.first = subLeading.first;
		}
	}
	return pair;
}

double muonMatch(const reco::Candidate*& particle, pat::Muon* target) {
	/*Performs matching checks between perticles. Returns dR for positive ID*/
	if (particle->pdgId() != target->pdgId()) return -1;
	double dR = ROOT::Math::VectorUtil::DeltaR(particle->p4(), target->p4());
	if (dR > 0.5) {
		if (mcDebug) std::cout << "Muon match failed on DR: " << dR << "\n";
		return -1;
	}
	double momDiff = std::abs(particle->p4().Pt()-target->p4().Pt());
	if (momDiff > 20) {
		if (mcDebug) std::cout << "Muon match failed on pT: " << momDiff << "\n";
		return -1;
	}
	return dR;
}

double muonSearch(const reco::Candidate*& particle, pat::Muon* target) {
	/*Recursive search through particle's decays for a particle matching reco particle. Returns dR separation or -1*/
	double match = -1;
	if (particle->numberOfDaughters() >= 2) {
		const reco::Candidate* d0 = particle->daughter(0);
		const reco::Candidate* d1 = particle->daughter(1);
		if (muonMatch(d0, target)) {
			return true;
		} else {
			match = muonSearch(d0, target);
		}
		if (match == -1 && muonMatch(d1, target)) {
			return true;
		} else {
			match = muonSearch(d1, target);
		}
	}
	return match;
}

bool checkBJets(pat::Jet* bjet0, pat::Jet* bjet1,
	const reco::Candidate*& gen_bjet0, const reco::Candidate*& gen_bjet1, double R) {
	/*Checks whether the particles are within their nearest jet*/
	//Associate particles to closest found jet___
	if (ROOT::Math::VectorUtil::DeltaR(gen_bjet0->p4(), bjet0->p4()) > 
		ROOT::Math::VectorUtil::DeltaR(gen_bjet0->p4(), bjet1->p4())) { //Wrong assignemnt; swap
		const reco::Candidate* temp = gen_bjet0;
		gen_bjet0 = gen_bjet1;
		gen_bjet1 = temp;
	}
	//___________________________________________
	//Check jets_________________________________
	double dR_0 = ROOT::Math::VectorUtil::DeltaR(gen_bjet0->p4(), bjet0->p4());
	double dR_1 = ROOT::Math::VectorUtil::DeltaR(gen_bjet1->p4(), bjet1->p4());
	if (dR_0 > R || dR_1 > R) { //particle(s) outside jet
		return false;
	} else {
		return true;
	}
}

bool truthFlag(TH1D* mcPlots,
	const reco::GenParticle* gen_hBB, const reco::GenParticle* gen_hTauTau,
	const reco::Candidate*& gen_bjet0, const reco::Candidate*& gen_bjet1, const reco::Candidate*& gen_tau0, const reco::Candidate*& gen_tau1,
	pat::Jet* bjet0, pat::Jet* bjet1, pat::Tau* tau, pat::Muon* muon) {
	/*Checks whether selected final states are correct*/
	double jetRadius = 0.5;
	mcPlots->Fill("MC-truth check", 1);
	//Check b jets_______________________________
	mcPlots->Fill("h->bb check", 1);
	if (debug) std::cout << "Checking b-jets\n";
	if (!checkBJets(bjet0, bjet1, gen_bjet0, gen_bjet1, jetRadius)) {
		if (mcDebug) std::cout << "MC check fails due to di-Jet on b-jets check\n";
		return false; //b-jet selection incorrect
	}
	if (debug) std::cout << "Both b jets confirmed\n";
	mcPlots->Fill("h->bb pass", 1);
	//___________________________________________
	//Check taus_________________________________
	if (debug) std::cout << "Checking taus\n";
	// std::vector<std::string> options;
	// boost::split(options, mode, boost::is_any_of(":"));
	mcPlots->Fill("h->#tau#tau check", 1);
	//if (options[0] == "tau" && options[1] == "tau") {
		//h->tau_h tau_h_________________________
		// if (!checkDiJet(branchJet, branchParticle, l_0, l_1, hTauTau, 15, &swap, (*plots)["tauMatch"], jetRadius)) {
		// 	if (debug) std::cout << "MC check fails due to di-Jet on tau-jets check\n";
		// 	chain->Delete();
		// 	return false; //tau-jet selection incorrect
		// }
		// if (swap) {
		// 	tau_0 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D2, branchParticle));
		// 	tau_1 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D1, branchParticle));
		// } else {
		// 	tau_0 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D1, branchParticle));
		// 	tau_1 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D2, branchParticle));
		// }
		// mcPlots->Fill(("h->#tau#tau->" + typeLookup(mode) + " pass").c_str(), 1);
		//_______________________________________
	//} else if ((options[0] == "tau" && options[1] == "muon") || (options[0] == "muon" && options[1] == "tau")) {
	//h->tau_h light-lepton__________________
	double dRMuon0 = muonSearch(gen_tau0, muon);
	double dRMuon1 = muonSearch(gen_tau1, muon);
	mcPlots->Fill("#tau->#mu check", 1);
	if ((dRMuon0 == -1) && (dRMuon1 == -1)) { //Neither taus decay to matched muons
		if (mcDebug) std::cout << "MC match failed due to neither tau decaying to matched muon\n";
		return false;
	}
	mcPlots->Fill("#tau->#mu pass", 1);
	mcPlots->Fill("#tau->#tau_{h} check", 1);
	double dRJet0 = ROOT::Math::VectorUtil::DeltaR(gen_tau0->p4(), tau->p4());
	double dRJet1 = ROOT::Math::VectorUtil::DeltaR(gen_tau1->p4(), tau->p4());
	if ((dRJet0 > jetRadius) && (dRJet1 > jetRadius)) { //Neither taus within tau jet
		if (mcDebug) std::cout << "MC match failed due to neither tau being within tau jet\n";
		return false; 
	}
	mcPlots->Fill("#tau->#tau_{h} pass", 1);
	mcPlots->Fill("#tau#tau assignement check", 1);
	if ((dRMuon0 == -1) && (dRMuon1 != -1)) { //tau1 decays to muon and tau0 does not
		if (dRJet0 > jetRadius) { //tau0 not within reco jet
			if (mcDebug) std::cout << "MC match failed due to non tau_l being within tau jet\n";
			return false;
		}
	}
	if ((dRMuon0 != -1) && (dRMuon1 == -1)) { //tau0 decays to muon and tau1 does not
		if (dRJet1 <= jetRadius) { //tau1 within reco jet
			const reco::Candidate* temp = gen_tau0; //Reassociate tau0 to tau_h
			gen_tau0 = gen_tau1;
			gen_tau1 = temp;
		} else {
			if (mcDebug) std::cout << "MC match failed due to non tau_l being within tau jet\n";
			return false;
		}
	}
	if ((dRMuon0 != -1) && (dRMuon1 != -1)) { //Both taus decay matched muon
		if ((dRJet0 <= jetRadius) && (dRJet1 <= jetRadius)) { //Both taus within tau jet
			if (dRMuon0 < dRMuon1) { //Choose by smallest angle to muon
				const reco::Candidate* temp = gen_tau0; //Reassociate tau0 to tau_h
				gen_tau0 = gen_tau1;
				gen_tau1 = temp;
			}
		}  else if ((dRJet0 > jetRadius) && (dRJet1 <= jetRadius)) { //Only tau1 within tau jet
			const reco::Candidate* temp = gen_tau0; //Reassociate tau0 to tau_h
			gen_tau0 = gen_tau1;
			gen_tau1 = temp;
		}
	}
	mcPlots->Fill("#tau#tau assignement pass", 1);
	mcPlots->Fill("h->#tau#tau pass", 1);
	//_______________________________________
	// } else {
		//h->light-lepton light-lepton___________
		//Load objects___________________________
		// GenParticle* higgs = (GenParticle*)branchParticle->At(hTauTau);
		// tau_0 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D1, branchParticle));
		// tau_1 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D2, branchParticle));
		// GenParticle *lightLepton_0, *lightLepton_1;
		// if (options[0] == "muon") {
		// 	lightLepton_0 = (GenParticle*)((Muon*)branchMuon->At(l_0))->Particle.GetObject();
		// } else if (options[0] == "electron") {
		// 	lightLepton_0 = (GenParticle*)((Electron*)branchElectron->At(l_0))->Particle.GetObject();
		// }
		// if (options[1] == "muon") {
		// 	lightLepton_1 = (GenParticle*)((Muon*)branchMuon->At(l_1))->Particle.GetObject();
		// } else if (options[1] == "electron") {
		// 	lightLepton_1 = (GenParticle*)((Electron*)branchElectron->At(l_1))->Particle.GetObject();
		// }
		// //_______________________________________
		// //Check taus_____________________________
		// int leptonMother_0 = ancestrySearch(lightLepton_0, tau_0, tau_1, branchParticle);
		// if (leptonMother_0 == -1) {
		// 	if (debug) std::cout << "MC check fails due to ancestry check\n";
		// 	chain->Delete();
		// 	return false; //Light lepton 0 did not come from tau decay
		// }
		// int leptonMother_1 = ancestrySearch(lightLepton_1, tau_0, tau_1, branchParticle);
		// if (leptonMother_1 == -1) {
		// 	if (debug) std::cout << "MC check fails due to ancestry check\n";
		// 	chain->Delete();
		// 	return false; //Light lepton 1 did not come from tau decay
		// }
		// if (leptonMother_0 == leptonMother_1) {
		// 	if (debug) std::cout << "MC check fails due to both leptons coming from same tau\n";
		// 	chain->Delete();
		// 	return false; //Leptons both came from same mother (somehow)
		// }
		// mcPlots->Fill(("h->#tau#tau->" + typeLookup(mode) + " pass").c_str(), 1);
		// if ((lightLepton_0->PT > lightLepton_1->PT & leptonMother_0 == 1) |
		// 	(lightLepton_0->PT < lightLepton_1->PT & leptonMother_0 == 0)) {
		// 	tau_0 = tau_1;
		// tau_1 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D1, branchParticle));
	// }
		//_______________________________________
		//_______________________________________
	if (debug) std::cout << "Both taus confirmed\n";
	//___________________________________________
	if (debug) std::cout << "Event accepted\n";
	mcPlots->Fill("MC-truth pass", 1);
	return true;
}

bool getGenParticles(edm::Handle<reco::GenParticleCollection> genParticles,
	int* gen_hBB_key, int* gen_hTauTau_key, TH1D* mcCuts, TH1D* higgDecay) {
	/*Point hbb and htautau to the Higgs*/
	bool hBBFound = false, hTauTauFound = false;
	int nHiggs = 0;
	mcCuts->Fill("hh->bb#tau#tau check", 1);
	for(size_t i = 0; i < genParticles->size(); ++ i) {
		const reco::GenParticle& p = (*genParticles)[i];
		if (std::abs(p.pdgId()) == 25) { //Particle is Higgs
			if (p.numberOfDaughters() >= 2) { //Daughters exists
				if (p.numberOfDaughters() > 2) std::cout << "N daughters: " << p.numberOfDaughters() << "\n";
				const reco::Candidate* d0 = p.daughter(0);
				const reco::Candidate* d1 = p.daughter(1);
				if (d0->pdgId() != 25 && d1->pdgId() != 25) {
					nHiggs++;
					higgDecay->Fill(std::abs(d0->pdgId()));
					higgDecay->Fill(std::abs(d1->pdgId()));
					if (std::abs(d0->pdgId()) == 5 && std::abs(d1->pdgId()) == 5) { //Daughters are b quarks
						hBBFound = true;
						*gen_hBB_key = i; //Point to Higgs
						if (hBBFound && hTauTauFound) { //h->bb and h->tautau found, so accept event
							mcCuts->Fill("hh->bb#tau#tau pass", 1);
							return true;
						}
					}
					if (std::abs(d0->pdgId()) == 15 && std::abs(d1->pdgId()) == 15) { //Daughters are taus
						hTauTauFound = true;
						*gen_hTauTau_key = i; //Point to Higgs
						if (hBBFound && hTauTauFound) { //h->bb and h->tautau found, so accept event
							mcCuts->Fill("hh->bb#tau#tau pass", 1);
							return true;
						}
					}
				}
			}
			if (nHiggs >= 2) break; //Both Higgs found
		}
	}
	std::cout << nHiggs << " Higgs found in signal event\n";
	return false; //Both h->bb and h->tautau not found
}

TMatrixD decomposeVector(math::XYZTLorentzVector* in) {
	TMatrixD out(3, 3);
	out(0, 0) = in->Px()*in->Px();
	out(0, 1) = in->Px()*in->Py();
	out(0, 2) = in->Px()*in->Pz();
	out(1, 0) = in->Py()*in->Px();
	out(1, 1) = in->Py()*in->Py();
	out(1, 2) = in->Py()*in->Pz();
	out(2, 0) = in->Pz()*in->Px();
	out(2, 1) = in->Pz()*in->Py();
	out(2, 2) = in->Pz()*in->Pz();
	return out;
}

void appendSphericity(TMatrixD* mat, double* div, math::XYZTLorentzVector* mom) {
	TMatrixD decomp = decomposeVector(mom);
	*mat += decomp;
	*div += pow(mom->P(), 2);
}

void appendSpherocity(TMatrixD* mat, double* div, math::XYZTLorentzVector* mom) {
	TMatrixD decomp = decomposeVector(mom);
	decomp *= 1/std::abs(mom->P());
	*mat += decomp;
	*div += std::abs(mom->P());
}

std::vector<double> getEigenValues(TMatrixD in) {
	/*Return vector of sorted, nomalised eigenvalues of parssed matrix*/
	TMatrixD eigenMatrix = TMatrixDEigen(in).GetEigenValues();
	std::vector<double> eigenValues(3);
	eigenValues[0] = eigenMatrix(0, 0);
	eigenValues[1] = eigenMatrix(1, 1);
	eigenValues[2] = eigenMatrix(2, 2);
	std::sort(eigenValues.begin(), eigenValues.end(), std::greater<double>());
	double sum = 0;
	for (double n : eigenValues)
		sum += n;
	std::for_each(eigenValues.begin(), eigenValues.end(), [sum](double i) { return i/sum; });
	return eigenValues;
}

void getEventShapes(std::vector<double> sphericityV, std::vector<double> spherocityV,
	double* sphericity, double* spherocity,
	double* aplanarity, double* aplanority,
	double* upsilon, double* dShape) {
	*sphericity = (3/2)*(sphericityV[1]+sphericityV[2]);
	*spherocity = (3/2)*(spherocityV[1]+spherocityV[2]);
	*aplanarity = 3*sphericityV[2]/2;
	*aplanority = 3*spherocityV[2]/2;
	*upsilon = sqrt(3.0)*(sphericityV[1]-sphericityV[2])/2;
	*dShape = 27*spherocityV[0]*spherocityV[1]*spherocityV[2];
}

void getGlobalEventInfo(math::XYZTLorentzVector* v_tau_0, math::XYZTLorentzVector* v_tau_1,
	math::XYZTLorentzVector* v_bJet_0, math::XYZTLorentzVector* v_bJet_1, math::XYZTLorentzVector* v_met,
	double*  hT, double*  sT, double* centrality, double* eVis) {
	/*Fills referenced variables with global event information*/
	if (debug) std::cout << "Getting global event info\n";
	//Reset variables____________________________
	*hT = 0;
	*sT = 0;
	*centrality = 0;
	*eVis = 0;
	//____________________________________________
	//HT__________________________________________
	*hT += v_bJet_0->Et();
	*hT += v_bJet_1->Et();
	*hT += v_tau_0->Et();
	//____________________________________________
	//ST__________________________________________
	*sT += *hT;
	*sT += v_tau_1->Pt();
	*sT += v_met->Pt();
	//____________________________________________
	//Centrality__________________________________
	*eVis += v_tau_0->E();
	*centrality += v_tau_0->Pt();
	*eVis += v_tau_1->E();
	*centrality += v_tau_1->Pt();
	*eVis += v_bJet_0->E();
	*centrality += v_bJet_0->Pt();
	*eVis += v_bJet_1->E();
	*centrality += v_bJet_1->Pt();
	*centrality /= *eVis;
	//___________________________________________
}

void getPrimaryEventShapes(math::XYZTLorentzVector* v_tau_0, math::XYZTLorentzVector* v_tau_1,
	math::XYZTLorentzVector* v_bJet_0, math::XYZTLorentzVector* v_bJet_1,
	double* sphericity, double* spherocity,
	double* aplanarity, double* aplanority,
	double* upsilon, double* dShape,
	double* sphericityEigen0, double* sphericityEigen1, double* sphericityEigen2,
	double* spherocityEigen0, double* spherocityEigen1, double* spherocityEigen2) {
	/*Sets values of referenced event-shape variables for final-states*/
	//Reset values________________________________
	*sphericity = 0;
	*spherocity = 0;
	*aplanarity = 0;
	*aplanority = 0;
	*upsilon = 0;
	*dShape = 0;
	*sphericityEigen0 = 0;
	*sphericityEigen1 = 0;
	*sphericityEigen2 = 0;
	*spherocityEigen0 = 0;
	*spherocityEigen1 = 0;
	*spherocityEigen2 = 0;
	//____________________________________________
	if (debug) std::cout << "Getting primary event shapes\n";
	TMatrixD sphericityT(3, 3), spherocityT(3, 3);
	double sphericityD = 0, spherocityD = 0;
	//Populate tensors___________________________
	appendSphericity(&sphericityT, &sphericityD, v_tau_0);
	appendSpherocity(&spherocityT, &spherocityD, v_tau_0);
	appendSphericity(&sphericityT, &sphericityD, v_tau_1);
	appendSpherocity(&spherocityT, &spherocityD, v_tau_1);
	appendSphericity(&sphericityT, &sphericityD, v_bJet_0);
	appendSpherocity(&spherocityT, &spherocityD, v_bJet_0);
	appendSphericity(&sphericityT, &sphericityD, v_bJet_1);
	appendSpherocity(&spherocityT, &spherocityD, v_bJet_1);
	sphericityT *= 1/sphericityD;
	spherocityT *= 1/spherocityD;
	//___________________________________________
	//Calculate event shapes_____________________
	if (debug) std::cout << "Calculating primary event shapes\n";
	std::vector<double> sphericityV = getEigenValues(sphericityT);
	std::vector<double> spherocityV = getEigenValues(spherocityT);
	getEventShapes(sphericityV, spherocityV,
		sphericity, spherocity,
		aplanarity, aplanority,
		upsilon, dShape);
	*sphericityEigen0 = sphericityV[0];
	*sphericityEigen1 = sphericityV[1];
	*sphericityEigen2 = sphericityV[2];
	*spherocityEigen0 = spherocityV[0];
	*spherocityEigen1 = spherocityV[1];
	*spherocityEigen2 = spherocityV[2];
	//___________________________________________
}

int main(int argc, char* argv[])
{
	using pat::Muon;

	gSystem->Load( "libFWCoreFWLite" );
	FWLiteEnabler::enable();

	// command line parser:
	optutl::CommandLineParser parser ("Analyze FWLite Histograms");
	parser.integerValue ("maxEvents"  ) =    -1;
	parser.integerValue ("outputEvery") = 10000;
	parser.stringValue  ("outputFile" ) = "analyzeFWLiteHistograms.root";
	parser.addOption ("runOnData",     optutl::CommandLineParser::kBool, "", false);
	parser.addOption ("setdebug",      optutl::CommandLineParser::kBool, "", false);
	parser.addOption ("monitorFile",   optutl::CommandLineParser::kString, "", "");

	parser.parseArguments (argc, argv);
	std::string outputFile_ = parser.stringValue("outputFile");
	std::vector<std::string> inputFiles_ = parser.stringVector("inputFiles");
	//bool runOnData = parser.boolValue("runOnData");
	bool runOnSignal = false;
	//debug = parser.boolValue("setdebug");

	fwlite::TFileService fs = fwlite::TFileService(outputFile_.c_str());

	double weight;

	math::XYZTLorentzVector t_0_p4, t_1_p4, bjet0_p4, bjet1_p4, met_p4, svFit_p4, hbb_p4, htt_p4, hh_p4;
	math::XYZTLorentzVector gen_tau0_p4, gen_tau1_p4, gen_bjet0_p4, gen_bjet1_p4, gen_hbb_p4, gen_htt_p4, gen_hh_p4;

	//Low-level variables________________________
	double t_0_px, t_0_py, t_0_pz, t_0_P, t_0_E, t_0_mass, t_0_mT; //Tau 0 variables
	double t_1_px, t_1_py, t_1_pz, t_1_P, t_1_E, t_1_mass, t_1_mT; //Tau 1 variables
	double b_0_px, b_0_py, b_0_pz, b_0_P, b_0_E, b_0_mass, b_0_csv, b_0_rawf, b_0_mva; //b-jet 0 variables
	double b_1_px, b_1_py, b_1_pz, b_1_P, b_1_E, b_1_mass, b_1_csv, b_1_rawf, b_1_mva; //b-jet 1 variables
	double mPT_px, mPT_py; //mPT_cov_00, mPT_cov_01, mPT_cov_10, mPT_cov_11; //Missing ET variables
	//___________________________________________
	//Reconstructed variables____________________
	double h_tt_px, h_tt_py, h_tt_pz, h_tt_P, h_tt_E, h_tt_mass; //Higgs->tau tau variables
	double h_tt_svFit_px, h_tt_svFit_py, h_tt_svFit_pz, h_tt_svFit_P, h_tt_svFit_E, h_tt_svFit_mass, h_tt_svFit_mT; //SVFit vector
	double h_bb_px, h_bb_py, h_bb_pz, h_bb_P, h_bb_E, h_bb_mass; //Higgs->bb variables
	double diH_px, diH_py, diH_pz, diH_P, diH_E, diH_mass; //di-Higgs variables
	double diH_kinFit_mass, diH_kinFit_chi2, diH_kinFit_conv; //Kinfit variables
	//___________________________________________
	//Twist______________________________________
	double twist_b_0_b_1, twist_b_0_t_0, twist_b_0_t_1;
	double twist_b_1_t_0, twist_b_1_t_1;
	double twist_t_0_t_1;
	double twist_h_bb_h_tt;
	//___________________________________________
	//DeltaR_____________________________________
	double dR_b_0_b_1, dR_b_0_t_0, dR_b_0_t_1;
	double dR_b_1_t_0, dR_b_1_t_1;
	double dR_t_0_t_1;
	double dR_h_bb_h_tt;
	//___________________________________________
	//Global event variables_____________________
	double nJets;
	double hT, hT_jets, sT, centrality, eVis; //Global kinematics
	double sphericity, spherocity, aplanarity, aplanority, upsilon, dShape; //Event shapes for primary objects
	double sphericityEigen0, sphericityEigen1, sphericityEigen2; //Eigenvalues for sphericity of primary objects
	double spherocityEigen0, spherocityEigen1, spherocityEigen2; //Eigenvalues for spherocity of primary objects
	//___________________________________________
	//Generator-level variables for regression and cuts
	/*double gen_t_0_px, gen_t_0_py, gen_t_0_pz, gen_t_0_E;*/ bool gen_t_0_match;//Tau 0 variables
	/*double gen_t_1_px, gen_t_1_py, gen_t_1_pz, gen_t_1_E;*/ bool gen_t_1_match; //Tau 1 variables
	/*double gen_b_0_px, gen_b_0_py, gen_b_0_pz, gen_b_0_E; //bool gen_b_0_match; //b-jet 0 variables
	double gen_b_1_px, gen_b_1_py, gen_b_1_p, gen_b_1_E; //bool gen_b_1_match; //b-jet 1 variables
	double gen_diH_px, gen_diH_py, gen_diH_pz, gen_diH_E, gen_diH_mass; //diHiggs variables
	double gen_h_bb_px, gen_h_bb_py, gen_h_bb_pz, gen_h_bb_E; //Higgs->bb variables
	double gen_h_tt_px, gen_h_tt_py, gen_h_tt_pz, gen_h_tt_E; //Higgs->tau tau variables*/
	bool gen_mctMatch; //MC truth match
	//___________________________________________

	//Data tree__________________________________
	TFileDirectory flatData = fs.mkdir("data");
	TTree* mu_tau_b_b = flatData.make<TTree>("mu_tau_b_b", "#mu #tau_{h} b #bar{b}");
	mu_tau_b_b->Branch("t_0_px", &t_0_px);
	mu_tau_b_b->Branch("t_0_py", &t_0_py);
	mu_tau_b_b->Branch("t_0_pz", &t_0_pz);
	mu_tau_b_b->Branch("t_0_P", &t_0_P);
	mu_tau_b_b->Branch("t_0_E", &t_0_E);
	mu_tau_b_b->Branch("t_0_mass", &t_0_mass);
	mu_tau_b_b->Branch("t_0_mT", &t_0_mT);

	mu_tau_b_b->Branch("t_1_px", &t_1_px);
	mu_tau_b_b->Branch("t_1_py", &t_1_py);
	mu_tau_b_b->Branch("t_1_pz", &t_1_pz);
	mu_tau_b_b->Branch("t_1_P", &t_1_P);
	mu_tau_b_b->Branch("t_1_E", &t_1_E);
	mu_tau_b_b->Branch("t_1_mass", &t_1_mass);
	mu_tau_b_b->Branch("t_1_mT", &t_1_mT);

	mu_tau_b_b->Branch("b_0_px", &b_0_px);
	mu_tau_b_b->Branch("b_0_py", &b_0_py);
	mu_tau_b_b->Branch("b_0_pz", &b_0_pz);
	mu_tau_b_b->Branch("b_0_P", &b_0_P);
	mu_tau_b_b->Branch("b_0_E", &b_0_E);
	mu_tau_b_b->Branch("b_0_mass", &b_0_mass);
	mu_tau_b_b->Branch("b_0_csv", &b_0_csv);
	mu_tau_b_b->Branch("b_0_rawf", &b_0_rawf);
	mu_tau_b_b->Branch("b_0_mva", &b_0_mva);

	mu_tau_b_b->Branch("b_1_px", &b_1_px);
	mu_tau_b_b->Branch("b_1_py", &b_1_py);
	mu_tau_b_b->Branch("b_1_pz", &b_1_pz);
	mu_tau_b_b->Branch("b_1_P", &b_1_P);
	mu_tau_b_b->Branch("b_1_E", &b_1_E);
	mu_tau_b_b->Branch("b_1_mass", &b_1_mass);
	mu_tau_b_b->Branch("b_1_csv", &b_1_csv);
	mu_tau_b_b->Branch("b_1_rawf", &b_1_rawf);
	mu_tau_b_b->Branch("b_1_mva", &b_1_mva);

	mu_tau_b_b->Branch("mPT_px", &mPT_px);
	mu_tau_b_b->Branch("mPT_py", &mPT_py);
	/*mu_tau_b_b->Branch("mPT_cov_00", &mPT_cov_00);
	mu_tau_b_b->Branch("mPT_cov_01", &mPT_cov_01);
	mu_tau_b_b->Branch("mPT_cov_10", &mPT_cov_10);
	mu_tau_b_b->Branch("mPT_cov_11", &mPT_cov_11);*/

	mu_tau_b_b->Branch("h_tt_px", &h_tt_px);
	mu_tau_b_b->Branch("h_tt_py", &h_tt_py);
	mu_tau_b_b->Branch("h_tt_pz", &h_tt_pz);
	mu_tau_b_b->Branch("h_tt_P", &h_tt_P);
	mu_tau_b_b->Branch("h_tt_E", &h_tt_E);
	mu_tau_b_b->Branch("h_tt_mass", &h_tt_mass);

	mu_tau_b_b->Branch("h_tt_svFit_px", &h_tt_svFit_px);
	mu_tau_b_b->Branch("h_tt_svFit_py", &h_tt_svFit_py);
	mu_tau_b_b->Branch("h_tt_svFit_pz", &h_tt_svFit_pz);
	mu_tau_b_b->Branch("h_tt_svFit_P", &h_tt_svFit_P);
	mu_tau_b_b->Branch("h_tt_svFit_E", &h_tt_svFit_E);
	mu_tau_b_b->Branch("h_tt_svFit_mass", &h_tt_svFit_mass);
	mu_tau_b_b->Branch("h_tt_svFit_mT", &h_tt_svFit_mT);

	mu_tau_b_b->Branch("h_bb_px", &h_bb_px);
	mu_tau_b_b->Branch("h_bb_py", &h_bb_py);
	mu_tau_b_b->Branch("h_bb_pz", &h_bb_pz);
	mu_tau_b_b->Branch("h_bb_P", &h_bb_P);
	mu_tau_b_b->Branch("h_bb_E", &h_bb_E);
	mu_tau_b_b->Branch("h_bb_mass", &h_bb_mass);

	mu_tau_b_b->Branch("diH_px", &diH_px);
	mu_tau_b_b->Branch("diH_py", &diH_py);
	mu_tau_b_b->Branch("diH_pz", &diH_pz);
	mu_tau_b_b->Branch("diH_P", &diH_P);
	mu_tau_b_b->Branch("diH_E", &diH_E);
	mu_tau_b_b->Branch("diH_mass", &diH_mass);

	mu_tau_b_b->Branch("diH_kinFit_mass", &diH_kinFit_mass);
	mu_tau_b_b->Branch("diH_kinFit_chi2", &diH_kinFit_chi2);
	mu_tau_b_b->Branch("diH_kinFit_conv", &diH_kinFit_conv);

	mu_tau_b_b->Branch("twist_b_0_b_1", &twist_b_0_b_1);
	mu_tau_b_b->Branch("twist_b_0_t_0", &twist_b_0_t_0);
	mu_tau_b_b->Branch("twist_b_0_t_1", &twist_b_0_t_1);
	mu_tau_b_b->Branch("twist_b_1_t_0", &twist_b_1_t_0);
	mu_tau_b_b->Branch("twist_b_1_t_1", &twist_b_1_t_1);
	mu_tau_b_b->Branch("twist_t_0_t_1", &twist_t_0_t_1);
	mu_tau_b_b->Branch("twist_h_bb_h_tt", &twist_h_bb_h_tt);

	mu_tau_b_b->Branch("twist_b_0_b_1", &dR_b_0_b_1);
	mu_tau_b_b->Branch("dR_b_0_t_0", &dR_b_0_t_0);
	mu_tau_b_b->Branch("dR_b_0_t_1", &dR_b_0_t_1);
	mu_tau_b_b->Branch("dR_b_1_t_0", &dR_b_1_t_0);
	mu_tau_b_b->Branch("dR_b_1_t_1", &dR_b_1_t_1);
	mu_tau_b_b->Branch("dR_t_0_t_1", &dR_t_0_t_1);
	mu_tau_b_b->Branch("dR_h_bb_h_tt", &dR_h_bb_h_tt);

	mu_tau_b_b->Branch("nJets", &nJets);
	mu_tau_b_b->Branch("hT", &hT);
	mu_tau_b_b->Branch("hT_jets", &hT_jets);
	mu_tau_b_b->Branch("sT", &sT);
	mu_tau_b_b->Branch("centrality", &centrality);
	mu_tau_b_b->Branch("eVis", &eVis);

	mu_tau_b_b->Branch("sphericity", &sphericity);
	mu_tau_b_b->Branch("spherocity", &spherocity);
	mu_tau_b_b->Branch("aplanarity", &aplanarity);
	mu_tau_b_b->Branch("aplanority", &aplanority);
	mu_tau_b_b->Branch("upsilon", &upsilon);
	mu_tau_b_b->Branch("dShape", &dShape);

	mu_tau_b_b->Branch("sphericityEigen0", &sphericityEigen0);
	mu_tau_b_b->Branch("sphericityEigen1", &sphericityEigen1);
	mu_tau_b_b->Branch("sphericityEigen2", &sphericityEigen2);
	mu_tau_b_b->Branch("spherocityEigen0", &spherocityEigen0);
	mu_tau_b_b->Branch("spherocityEigen1", &spherocityEigen1);
	mu_tau_b_b->Branch("spherocityEigen2", &spherocityEigen2);

/*	mu_tau_b_b->Branch("gen_t_0_px", &gen_t_0_px);
	mu_tau_b_b->Branch("gen_t_0_py", &gen_t_0_py);
	mu_tau_b_b->Branch("gen_t_0_pz", &gen_t_0_pz);
	mu_tau_b_b->Branch("gen_t_0_E", &gen_t_0_E); */
	mu_tau_b_b->Branch("gen_t_0_match", &gen_t_0_match);

/*	mu_tau_b_b->Branch("gen_t_1_px", &gen_t_1_px);
	mu_tau_b_b->Branch("gen_t_1_py", &gen_t_1_py);
	mu_tau_b_b->Branch("gen_t_1_pz", &gen_t_1_pz);
	mu_tau_b_b->Branch("gen_t_1_E", &gen_t_1_E);*/
	mu_tau_b_b->Branch("gen_t_1_match", &gen_t_1_match);

/*	mu_tau_b_b->Branch("gen_b_0_px", &gen_b_0_px);
	mu_tau_b_b->Branch("gen_b_0_py", &gen_b_0_py);
	mu_tau_b_b->Branch("gen_b_0_pz", &gen_b_0_pz);
	mu_tau_b_b->Branch("gen_b_0_E", &gen_b_0_E);*/
	//mu_tau_b_b->Branch("gen_b_0_match", &gen_b_0_match);

/*	mu_tau_b_b->Branch("gen_b_1_px", &gen_b_1_px);
	mu_tau_b_b->Branch("gen_b_1_py", &gen_b_1_py);
	mu_tau_b_b->Branch("gen_b_1_pz", &gen_b_1_pz);
	mu_tau_b_b->Branch("gen_b_1_E", &gen_b_1_E);*/
	//mu_tau_b_b->Branch("gen_b_1_match", &gen_b_1_match);

/*	mu_tau_b_b->Branch("gen_diH_px", &gen_diH_px);
	mu_tau_b_b->Branch("gen_diH_py", &gen_diH_py);
	mu_tau_b_b->Branch("gen_diH_pz", &gen_diH_pz);
	mu_tau_b_b->Branch("gen_diH_E", &gen_diH_E);
	mu_tau_b_b->Branch("gen_diH_mass", &gen_diH_mass);

	mu_tau_b_b->Branch("gen_h_bb_px", &gen_h_bb_px);
	mu_tau_b_b->Branch("gen_h_bb_py", &gen_h_bb_py);
	mu_tau_b_b->Branch("gen_h_bb_pz", &gen_h_bb_pz);
	mu_tau_b_b->Branch("gen_h_bb_E", &gen_h_bb_E);

	mu_tau_b_b->Branch("gen_h_tt_px", &gen_h_tt_px);
	mu_tau_b_b->Branch("gen_h_tt_py", &gen_h_tt_py);
	mu_tau_b_b->Branch("gen_h_tt_pz", &gen_h_tt_pz);
	mu_tau_b_b->Branch("gen_h_tt_E", &gen_h_tt_E);*/

	mu_tau_b_b->Branch("gen_mctMatch", &gen_mctMatch);
	mu_tau_b_b->Branch("gen_weight", &weight);
	//___________________________________________

	// loop over all input files:
	for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile){
		std::string filename = inputFiles_[iFile];
		TFile* inFile = TFile::Open(filename.c_str());

		//if(filename.find("Data") != std::string::npos || filename.find("Run201") != std::string::npos ) runOnData = true;
		if(filename.find("ToHHTo2B2Tau") != std::string::npos) runOnSignal = true;

		if( inFile ){

			TTreeReader reader("muTau", inFile);

			//Define event features______________
			//General info_______________________
			TTreeReaderValue<double> r_weight(reader, "weight_total");
			TTreeReaderValue<unsigned int> r_njets(reader, "n_jets");
			TTreeReaderValue<float> r_jet_HT(reader, "ht_other_jets");
			//Tau_0______________________________
			TTreeReaderValue<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> > > r_t_0_p4(reader, "p4_2");
			TTreeReaderValue<int> r_gen_t_0_match(reader, "gen_match_2");
			//Tau_1______________________________
			TTreeReaderValue<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> > > r_t_1_p4(reader, "p4_1");
			TTreeReaderValue<int> r_gen_t_1_match(reader, "gen_match_1");
			//MET________________________________
			TTreeReaderValue<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> > > r_met_p4(reader, "pfMET_p4");
			//TTreeReaderValue<ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepStd<double,2,2> > > r_met_cov(reader, "pfMET_cov");
			//Jets_______________________________
			TTreeReaderValue<std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > > > r_jets_p4(reader, "jets_p4");
			TTreeReaderValue<std::vector<float> > r_jets_csv(reader, "jets_csv");
			TTreeReaderValue<std::vector<float> > r_jets_rawf(reader, "jets_rawf");
			TTreeReaderValue<std::vector<float> > r_jets_mva(reader, "jets_mva");
			//SVFit______________________________
			TTreeReaderValue<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> > > r_svfit_p4(reader, "SVfit_p4");
			TTreeReaderValue<float> r_svfit_mT(reader, "SVfit_mt");
			//KinFit______________________________
			TTreeReaderValue<std::vector<float> > r_kinFit_mH(reader, "kinFit_m");
			TTreeReaderValue<std::vector<float> > r_kinFit_chi2(reader, "kinFit_chi2");
			TTreeReaderValue<std::vector<int> > r_kinFit_conv(reader, "kinFit_convergence");
			//___________________________________

			while (reader.Next()) {
				//Getting objects____________________
				//General info_______________________
				weight = *r_weight;
				nJets = *r_njets;
				hT_jets = *r_jet_HT;
				//MET________________________________
				met_p4 = *r_met_p4;
				mPT_px = met_p4.Px();
				mPT_py = met_p4.Py();
				/*mPT_cov_00 = (*r_met_cov)(0, 0);
				mPT_cov_01 = (*r_met_cov)(0, 1);
				mPT_cov_10 = (*r_met_cov)(1, 0);
				mPT_cov_11 = (*r_met_cov)(1, 1);*/
				//Tau_0______________________________
				t_0_p4 = *r_t_0_p4;
				t_0_px = t_0_p4.Px();
				t_0_py = t_0_p4.Py();
				t_0_pz = t_0_p4.Pz();
				t_0_P = t_0_p4.P();
				t_0_E = t_0_p4.E();
				t_0_mass = t_0_p4.M();
				t_0_mT = sqrt(2 * t_0_p4.Pt() * met_p4.Pt() * (1-cos(t_0_p4.Phi()-met_p4.Phi()))); //TODO: Generalise this
				gen_t_0_match = *r_gen_t_0_match;
				//Tau_1______________________________
				t_1_p4 = *r_t_1_p4;
				t_1_px = t_1_p4.Px();
				t_1_py = t_1_p4.Py();
				t_1_pz = t_1_p4.Pz();
				t_1_P = t_1_p4.P();
				t_1_E = t_1_p4.E();
				t_1_mass = t_1_p4.M();
				t_1_mT = sqrt(2 * t_1_p4.Pt() * met_p4.Pt() * (1-cos(t_1_p4.Phi()-met_p4.Phi()))); //TODO: Generalise this
				gen_t_1_match = *r_gen_t_1_match;
				//Jets_______________________________
				bjet0_p4 = (*r_jets_p4)[0];
				bjet1_p4 = (*r_jets_p4)[1];
				if (bjet0_p4.Pt() > bjet1_p4.Pt()) { //Order jets by pT
					bjet0_p4 = (*r_jets_p4)[1];
					b_0_csv = (*r_jets_csv)[1];
					b_0_rawf = (*r_jets_rawf)[1];
					b_0_mva = (*r_jets_mva)[1];

					bjet1_p4 = (*r_jets_p4)[0];
					b_1_csv = (*r_jets_csv)[0];
					b_1_rawf = (*r_jets_rawf)[0];
					b_1_mva = (*r_jets_mva)[0];
				} else {
					b_0_csv = (*r_jets_csv)[0];
					b_0_rawf = (*r_jets_rawf)[0];
					b_0_mva = (*r_jets_mva)[0];

					b_1_csv = (*r_jets_csv)[1];
					b_1_rawf = (*r_jets_rawf)[1];
					b_1_mva = (*r_jets_mva)[1];
				}

				b_0_px = bjet0_p4.Px();
				b_0_py = bjet0_p4.Py();
				b_0_pz = bjet0_p4.Pz();
				b_0_P = bjet0_p4.P();
				b_0_E = bjet0_p4.E();
				b_0_mass = bjet0_p4.M();

				b_1_px = bjet1_p4.Px();
				b_1_py = bjet1_p4.Py();
				b_1_pz = bjet1_p4.Pz();
				b_1_P = bjet1_p4.P();
				b_1_E = bjet1_p4.E();
				b_1_mass = bjet1_p4.M();
				//SVFit_______________________________
				svFit_p4 = *r_svfit_p4;
				h_tt_svFit_px = svFit_p4.Px();
				h_tt_svFit_py = svFit_p4.Py();
				h_tt_svFit_pz = svFit_p4.Pz();
				h_tt_svFit_P = svFit_p4.P();
				h_tt_svFit_E = svFit_p4.E();
				h_tt_svFit_mass = svFit_p4.M();
				h_tt_svFit_mT = *r_svfit_mT;
				//KinFit______________________________
				diH_kinFit_mass = (*r_kinFit_mH)[0];
				diH_kinFit_chi2 = (*r_kinFit_chi2)[0];
				diH_kinFit_conv = (*r_kinFit_conv)[0];
				//h->bb_______________________________
				hbb_p4 = bjet0_p4+bjet1_p4;
				h_bb_px = hbb_p4.Px();
				h_bb_py = hbb_p4.Py();
				h_bb_pz = hbb_p4.Pz();
				h_bb_P = hbb_p4.P();
				h_bb_E = hbb_p4.E();
				h_bb_mass = hbb_p4.M();
				//h->tautau___________________________
				htt_p4 = t_0_p4+t_1_p4+met_p4;
				h_tt_px = htt_p4.Px();
				h_tt_py = htt_p4.Py();
				h_tt_pz = htt_p4.Pz();
				h_tt_P = htt_p4.P();
				h_tt_E = htt_p4.E();
				h_tt_mass = htt_p4.M();
				//Di-higgs____________________________
				hh_p4 = hbb_p4+htt_p4;
				diH_px = hh_p4.Px();
				diH_py = hh_p4.Py();
				diH_pz = hh_p4.Pz();
				diH_P = hh_p4.P();
				diH_E = hh_p4.E();
				diH_mass = hh_p4.M();
				//Shapes__________________________
				getGlobalEventInfo(&t_0_p4, &t_1_p4, &bjet0_p4, &bjet0_p4, &met_p4,
					&hT, &sT, &centrality, &eVis);
				getPrimaryEventShapes(&t_0_p4, &t_1_p4, &bjet0_p4, &bjet0_p4,
					&sphericity, &spherocity,
					&aplanarity, &aplanority,
					&upsilon, &dShape,
					&sphericityEigen0, &sphericityEigen1, &sphericityEigen2,
					&spherocityEigen0, &spherocityEigen1, &spherocityEigen2);
				//________________________________
				//Twist___________________________
				twist_b_0_b_1 = atan(std::abs(ROOT::Math::VectorUtil::DeltaPhi(bjet0_p4, bjet1_p4)/(bjet0_p4.Eta()-bjet1_p4.Eta())));
				twist_b_0_t_0 = atan(std::abs(ROOT::Math::VectorUtil::DeltaPhi(bjet0_p4, t_0_p4)/(bjet0_p4.Eta()-t_0_p4.Eta())));
				twist_b_0_t_1 = atan(std::abs(ROOT::Math::VectorUtil::DeltaPhi(bjet0_p4, t_1_p4)/(bjet0_p4.Eta()-t_1_p4.Eta())));
				twist_b_1_t_0 = atan(std::abs(ROOT::Math::VectorUtil::DeltaPhi(bjet1_p4, t_0_p4)/(bjet1_p4.Eta()-t_0_p4.Eta())));
				twist_b_1_t_1 = atan(std::abs(ROOT::Math::VectorUtil::DeltaPhi(bjet1_p4, t_1_p4)/(bjet1_p4.Eta()-t_1_p4.Eta())));
				twist_t_0_t_1 = atan(std::abs(ROOT::Math::VectorUtil::DeltaPhi(t_0_p4, t_1_p4)/(t_0_p4.Eta()-t_1_p4.Eta())));
				twist_h_bb_h_tt = atan(std::abs(ROOT::Math::VectorUtil::DeltaPhi(hbb_p4, htt_p4)/(hbb_p4.Eta()-htt_p4.Eta())));
				//____________________________________
				//dR__________________________________
				dR_b_0_b_1 = ROOT::Math::VectorUtil::DeltaR(bjet0_p4, bjet1_p4);
				dR_b_0_t_0 = ROOT::Math::VectorUtil::DeltaR(bjet0_p4, t_0_p4);
				dR_b_0_t_1 = ROOT::Math::VectorUtil::DeltaR(bjet0_p4, t_1_p4);
				dR_b_1_t_0 = ROOT::Math::VectorUtil::DeltaR(bjet1_p4, t_0_p4);
				dR_b_1_t_1 = ROOT::Math::VectorUtil::DeltaR(bjet1_p4, t_1_p4);
				dR_t_0_t_1 = ROOT::Math::VectorUtil::DeltaR(t_0_p4, t_1_p4);
				dR_h_bb_h_tt = ROOT::Math::VectorUtil::DeltaR(hbb_p4, htt_p4);
				//____________________________________
				
				//MC truth_____________________
				//Reset variables______________
				/*gen_t_0_px = 0;
				gen_t_0_py = 0;
				gen_t_0_pz = 0;
				gen_t_0_E = 0;
				gen_t_1_px = 0;
				gen_t_1_py = 0;
				gen_t_1_pz = 0;
				gen_t_1_E = 0;
				gen_b_0_px = 0;
				gen_b_0_py = 0;
				gen_b_0_pz = 0;
				gen_b_0_E = 0;
				//gen_b_0_match = 0;
				gen_b_1_px = 0;
				gen_b_1_py = 0;
				gen_b_1_pz = 0;
				gen_b_1_E = 0;
				//gen_b_1_match = 0;
				gen_diH_px = 0;
				gen_diH_py = 0;
				gen_diH_pz = 0;
				gen_diH_E = 0;
				gen_diH_mass = 0;
				gen_h_bb_px = 0;
				gen_h_bb_py = 0;
				gen_h_bb_pz = 0;
				gen_h_bb_E = 0;
				gen_h_tt_px = 0;
				gen_h_tt_py = 0;
				gen_h_tt_pz = 0;
				gen_h_tt_E = 0;*/
				//_____________________________
				//_____________________________
				if (runOnSignal) { //TODO:Add this
					continue;
					/*//Get gen info______________
					int gen_hBB_key, gen_hTauTau_key;
					if (getGenParticles(genParticles, &gen_hBB_key, &gen_hTauTau_key, mcCuts, higgsDecay)) { //If both Higgs found
						const reco::GenParticle& gen_hBB = (*genParticles)[gen_hBB_key];
						const reco::GenParticle& gen_hTauTau = (*genParticles)[gen_hTauTau_key];
						const reco::Candidate* gen_bjet0 = gen_hBB.daughter(0);
						const reco::Candidate* gen_bjet1 = gen_hBB.daughter(1);
						const reco::Candidate* gen_tau0 = gen_hTauTau.daughter(0);
						const reco::Candidate* gen_tau1 = gen_hTauTau.daughter(1);
						//__________________________
						//Check FSs_________________
						gen_mctMatch = truthFlag(mcCuts,//Checks final-state selection was correct
							&gen_hBB, &gen_hTauTau, gen_bjet0, gen_bjet1, gen_tau0, gen_tau1,
							&bjet1, &bjet2, &tau, &muon);
						//__________________________
						//Get 4-momenta_____________
						gen_hbb_p4 = gen_hBB.p4();
						gen_htt_p4 = gen_hTauTau.p4();
						gen_hh_p4 = gen_hbb_p4 + gen_htt_p4;
						gen_tau0_p4 = gen_tau0->p4();
						gen_tau1_p4 = gen_tau1->p4();
						gen_bjet0_p4 = gen_bjet0->p4();
						gen_bjet1_p4 = gen_bjet1->p4();
						//__________________________
						//Decompose info____________
						gen_t_0_px = gen_tau0_p4.Px();
						gen_t_0_py = gen_tau0_p4.Py();
						gen_t_0_pz = gen_tau0_p4.Pz();
						gen_t_0_E = gen_tau0_p4.E();
						gen_t_1_px = gen_tau1_p4.Px();
						gen_t_1_py = gen_tau1_p4.Py();
						gen_t_1_pz = gen_tau1_p4.Pz();
						gen_t_1_E = gen_tau1_p4.E();
						gen_b_0_px = gen_bjet0_p4.Px();
						gen_b_0_py = gen_bjet0_p4.Py();
						gen_b_0_pz = gen_bjet0_p4.Pz();
						gen_b_0_E = gen_bjet0_p4.E();
						gen_b_1_px = gen_bjet1_p4.Px();
						gen_b_1_py = gen_bjet1_p4.Py();
						gen_b_1_pz = gen_bjet1_p4.Pz();
						gen_b_1_E = gen_bjet1_p4.E();
						gen_diH_px = gen_hh_p4.Px();
						gen_diH_py = gen_hh_p4.Py();
						gen_diH_pz = gen_hh_p4.Pz();
						gen_diH_E = gen_hh_p4.E();
						gen_diH_mass = gen_hh_p4.M();
						gen_h_bb_px = gen_hbb_p4.Px();
						gen_h_bb_py = gen_hbb_p4.Py();
						gen_h_bb_pz = gen_hbb_p4.Pz();
						gen_h_bb_E = gen_hbb_p4.E();
						gen_h_tt_px = gen_htt_p4.Px();
						gen_h_tt_py = gen_htt_p4.Py();
						gen_h_tt_pz = gen_htt_p4.Pz();
						gen_h_tt_E = gen_htt_p4.E();
						//__________________________
					}*/
				}
				//_____________________________
				//_____________________________
				mu_tau_b_b->Fill();
				//________________________________

			}
		}
		return 0;
	}
}