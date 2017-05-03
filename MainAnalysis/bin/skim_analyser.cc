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

std::pair<int, int> getJets(edm::Handle<std::vector<pat::Jet>> selectedjets, std::string mode="mass") {
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
	}  else if (mode == "CSV") { //Highest CSV
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

bool truthFlag(edm::Handle<reco::GenParticleCollection> genParticles, TH1D* mcPlots,
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
	int maxEvents_ = parser.integerValue("maxEvents");
	unsigned int outputEvery_ = parser.integerValue("outputEvery");
	std::string outputFile_ = parser.stringValue("outputFile");
	std::string monitorFile_ = parser.stringValue("monitorFile");
	std::vector<std::string> inputFiles_ = parser.stringVector("inputFiles");
	bool runOnData = parser.boolValue("runOnData");
	//debug = parser.boolValue("setdebug");

	// book a set of histograms:
	fwlite::TFileService fs = fwlite::TFileService(outputFile_.c_str());
	TFileDirectory dir_noPU = fs.mkdir("noPU");
	TH1F* hnVert_noPU  = dir_noPU.make<TH1F>("nVert" ,   "nVert" ,  50,    0,   50);

	TFileDirectory dir_weighted = fs.mkdir("weighted");
	TH1F* hmuonPt      = dir_weighted.make<TH1F>("muonPt"  , "pt"  ,   30,   0., 300.);
	TH1F* hmuonEta     = dir_weighted.make<TH1F>("muonEta" , "eta" ,   100,  -3.,   3.);
	TH1F* hmuonPhi     = dir_weighted.make<TH1F>("muonPhi" , "phi" ,   100,  -5.,   5.); 
	TH1F* hnVert       = dir_weighted.make<TH1F>("nVert" ,   "nVert" ,  50,    0,   50);   
	TH1F* hSVFit       = dir_weighted.make<TH1F>("svFitMass"  , "svFitMass"  ,   50,   0., 500.);
	TH1F* h_m_kinfit   = dir_weighted.make<TH1F>("m_kinfit", "KinFit Mass" ,  50,    0,   500);
	TH1F* h_P_fit      = dir_weighted.make<TH1F>("P_chi2" ,  "FitProb" ,  22,    0,   1.1);
	TH1F* h_m_bbar     = dir_weighted.make<TH1F>("m_bbar" ,  "invariant di-bjet mass" ,  40,    0,   400);
	TH1F* h_m_tauvis   = dir_weighted.make<TH1F>("m_tauvis" ,"invariant visible di-tau mass" ,  40,    0,   400);
	TH1F* h_m_t_mu     = dir_weighted.make<TH1F>("m_t_mu" ,  "transverse mass" ,  30,    0,   300);
	TH1F* h_met        = dir_weighted.make<TH1F>("m_met" ,  "missing transverse energy" ,  30,    0,   300);

	TH1F* h_cov_xx     = dir_weighted.make<TH1F>("h_cov_xx" ,  "h_cov_xx" ,  20,    0,   100);
	TH1F* h_cov_xy     = dir_weighted.make<TH1F>("h_cov_xy" ,  "h_cov_xy" ,  40,    -200,   200);
	TH1F* h_cov_yx     = dir_weighted.make<TH1F>("h_cov_yx" ,  "h_cov_yx" ,  40,    -200,   200);
	TH1F* h_cov_yy     = dir_weighted.make<TH1F>("h_cov_yy" ,  "h_cov_yy" ,  20,    0,   100);

	TH1F* h_monitor       = dir_weighted.make<TH1F>("h_monitor" ,  "h_monitor" ,  10,    0,  10);

	TFileDirectory dir_weights = fs.mkdir("weights");
	TH1F* hPUWeight = dir_weights.make<TH1F>("hPUWeight" ,   "hPUWeight" ,  400,    -2,   2);   
	TH1F* hBTaggingSF   = dir_weights.make<TH1F>("hBTaggingSF" ,   "hBTaggingSF" ,  400,    -2,   2);   
	TH1F* hMuonIdScaleFactor = dir_weights.make<TH1F>("hMuonIdScaleFactor" , "hMuonIdScaleFactor" ,  400,    -2,   2);
	TH1F* hMuonTriggerScaleFactor = dir_weights.make<TH1F>("hMuonTriggerScaleFactor" , "hMuonTriggerScaleFactor" ,  400,    -2,   2);
	TH1F* hLumiWeight = dir_weights.make<TH1F>("hLumiWeight" ,   "hLumiWeight" ,  10000,    0,   0.1);

	TFileDirectory dir_mcTruth = fs.mkdir("mcChecks");   
	TH1D* higgsDecay = dir_mcTruth.make<TH1D>("mcTruth_higgsDecay", "Higgs product |PID|", 50, 0, 50);
	TH1D* mcCuts = dir_mcTruth.make<TH1D>("mcTruth_cutFlow", "MC Truth Cuts", 14, -7.0, 7.0);
	mcCuts->GetXaxis()->SetBinLabel(1, "hh->bb#tau#tau check");
	mcCuts->GetXaxis()->SetBinLabel(2, "hh->bb#tau#tau pass");
	mcCuts->GetXaxis()->SetBinLabel(3, "MC-truth check");
	mcCuts->GetXaxis()->SetBinLabel(4, "MC-truth pass");
	mcCuts->GetXaxis()->SetBinLabel(5, "h->bb check");
	mcCuts->GetXaxis()->SetBinLabel(6, "h->bb pass");
	mcCuts->GetXaxis()->SetBinLabel(7, "h->#tau#tau check");
	mcCuts->GetXaxis()->SetBinLabel(8, "h->#tau#tau pass");
	mcCuts->GetXaxis()->SetBinLabel(9, "#tau->#mu check");
	mcCuts->GetXaxis()->SetBinLabel(10, "#tau->#mu pass");
	mcCuts->GetXaxis()->SetBinLabel(11, "#tau->#tau_{h} check");
	mcCuts->GetXaxis()->SetBinLabel(12, "#tau->#tau_{h} pass");
	mcCuts->GetXaxis()->SetBinLabel(13, "#tau#tau assignement check");
	mcCuts->GetXaxis()->SetBinLabel(14, "#tau#tau assignement pass");
	mcCuts->GetXaxis()->SetTitle("Cuts");
	mcCuts->GetYaxis()->SetTitle("Events");

	double xsec = -1;
	int nev_total = 0;
	int nev_selected = 0;
	double weight = 1;
	double weight_noPU = 1;
	double weight_lumi = 1;
	double weight_PU = 1;
	double weight_btag = 1;
	double svfitmass = 1;
	double weight_muon_IdScaleFactor = 1;
	double weight_muon_triggerScaleFactor = 1;

	double m_kinfit = -1;
	double FitProb = -1;
	double m_t_mu = -1;

	bool invIso = false;
	bool LS = false;
	bool runOnSignal = false;

	math::XYZTLorentzVector muon_p4, tau_p4, bjet0_p4, bjet1_p4, met_p4, hbb_p4, htt_p4, hh_p4;
	math::XYZTLorentzVector gen_tau0_p4, gen_tau1_p4, gen_bjet0_p4, gen_bjet1_p4, gen_hbb_p4, gen_htt_p4, gen_hh_p4;

	//Low-level variables________________________
	double t_0_pT, t_0_eta, t_0_phi, t_0_mass; //Tau 0 variables
	double t_1_pT, t_1_eta, t_1_phi, t_1_mass; //Tau 1 variables
	double b_0_pT, b_0_eta, b_0_phi, b_0_mass; //b-jet 0 variables
	double b_1_pT, b_1_eta, b_1_phi, b_1_mass; //b-jet 1 variables
	double mPT_pT, mPT_phi; //Missing ET variables
	//___________________________________________
	//Reconstructed variables____________________
	double h_tt_pT, h_tt_eta, h_tt_phi, h_tt_mass, h_tt_svFit_mass; //Higgs 0 variables
	double h_bb_pT, h_bb_eta, h_bb_phi, h_bb_mass; //Higgs 1 variables
	double diH_pT, diH_eta, diH_phi, diH_mass, diH_kinFit_mass, diH_kinFit_prob; //di-Higgs variables
	double mT; //Transverse mass
	//___________________________________________
	//Global event variables_____________________
	double hT, sT, centrality, eVis; //Global kinematics
	double sphericity, spherocity, aplanarity, aplanority, upsilon, dShape; //Event shapes for primary objects
	double sphericityEigen0, sphericityEigen1, sphericityEigen2; //Eigenvalues for sphericity of primary objects
	double spherocityEigen0, spherocityEigen1, spherocityEigen2; //Eigenvalues for spherocity of primary objects
	//___________________________________________
	//Generator-level variables for regression and cuts
	double gen_t_0_pT, gen_t_0_eta, gen_t_0_phi, gen_t_0_E; //Tau 0 variables
	double gen_t_1_pT, gen_t_1_eta, gen_t_1_phi, gen_t_1_E; //Tau 1 variables
	double gen_b_0_pT, gen_b_0_eta, gen_b_0_phi, gen_b_0_E; //b-jet 0 variables
	double gen_b_1_pT, gen_b_1_eta, gen_b_1_phi, gen_b_1_E; //b-jet 1 variables
	double gen_diH_pT, gen_diH_eta, gen_diH_phi, gen_diH_E, gen_diH_mass; //diHiggs variables
	double gen_h_bb_pT, gen_h_bb_eta, gen_h_bb_phi, gen_h_bb_E; //Higgs->bb variables
	double gen_h_tt_pT, gen_h_tt_eta, gen_h_tt_phi, gen_h_tt_E; //Higgs->tau tau variables
	bool gen_mctMatch; //MC truth match
	//___________________________________________

	//Data tree__________________________________
	TFileDirectory flatData = fs.mkdir("data");
	TTree* mu_tau_b_b = flatData.make<TTree>("mu_tau_b_b", "#mu #tau_{h} b #bar{b}");
	mu_tau_b_b->Branch("t_0_pT", &t_0_pT);
	mu_tau_b_b->Branch("t_0_eta", &t_0_eta);
	mu_tau_b_b->Branch("t_0_phi", &t_0_phi);
	mu_tau_b_b->Branch("t_0_mass", &t_0_mass);
	mu_tau_b_b->Branch("t_1_pT", &t_1_pT);
	mu_tau_b_b->Branch("t_1_eta", &t_1_eta);
	mu_tau_b_b->Branch("t_1_phi", &t_1_phi);
	mu_tau_b_b->Branch("t_1_mass", &t_1_mass);
	mu_tau_b_b->Branch("b_0_pT", &b_0_pT);
	mu_tau_b_b->Branch("b_0_eta", &b_0_eta);
	mu_tau_b_b->Branch("b_0_phi", &b_0_phi);
	mu_tau_b_b->Branch("b_0_mass", &b_0_mass);
	mu_tau_b_b->Branch("b_1_pT", &b_1_pT);
	mu_tau_b_b->Branch("b_1_eta", &b_1_eta);
	mu_tau_b_b->Branch("b_1_phi", &b_1_phi);
	mu_tau_b_b->Branch("b_1_mass", &b_1_mass);
	mu_tau_b_b->Branch("mPT_pT", &mPT_pT);
	mu_tau_b_b->Branch("mPT_phi", &mPT_phi);
	mu_tau_b_b->Branch("h_tt_pT", &h_tt_pT);
	mu_tau_b_b->Branch("h_tt_eta", &h_tt_eta);
	mu_tau_b_b->Branch("h_tt_phi", &h_tt_phi);
	mu_tau_b_b->Branch("h_tt_mass", &h_tt_mass);
	mu_tau_b_b->Branch("h_tt_svFit_mass", &h_tt_svFit_mass);
	mu_tau_b_b->Branch("h_bb_pT", &h_bb_pT);
	mu_tau_b_b->Branch("h_bb_eta", &h_bb_eta);
	mu_tau_b_b->Branch("h_bb_phi", &h_bb_phi);
	mu_tau_b_b->Branch("h_bb_mass", &h_bb_mass);
	mu_tau_b_b->Branch("diH_pT", &diH_pT);
	mu_tau_b_b->Branch("diH_eta", &diH_eta);
	mu_tau_b_b->Branch("diH_phi", &diH_phi);
	mu_tau_b_b->Branch("diH_mass", &diH_mass);
	mu_tau_b_b->Branch("diH_kinFit_mass", &diH_kinFit_mass);
	mu_tau_b_b->Branch("diH_kinFit_prob", &diH_kinFit_prob);
	mu_tau_b_b->Branch("mT", &mT);
	mu_tau_b_b->Branch("hT", &hT);
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
	mu_tau_b_b->Branch("gen_t_0_pT", &gen_t_0_pT);
	mu_tau_b_b->Branch("gen_t_0_eta", &gen_t_0_eta);
	mu_tau_b_b->Branch("gen_t_0_phi", &gen_t_0_phi);
	mu_tau_b_b->Branch("gen_t_0_E", &gen_t_0_E);
	mu_tau_b_b->Branch("gen_t_1_pT", &gen_t_1_pT);
	mu_tau_b_b->Branch("gen_t_1_eta", &gen_t_1_eta);
	mu_tau_b_b->Branch("gen_t_1_phi", &gen_t_1_phi);
	mu_tau_b_b->Branch("gen_t_1_E", &gen_t_1_E);
	mu_tau_b_b->Branch("gen_b_0_pT", &gen_b_0_pT);
	mu_tau_b_b->Branch("gen_b_0_eta", &gen_b_0_eta);
	mu_tau_b_b->Branch("gen_b_0_phi", &gen_b_0_phi);
	mu_tau_b_b->Branch("gen_b_0_E", &gen_b_0_E);
	mu_tau_b_b->Branch("gen_b_1_pT", &gen_b_1_pT);
	mu_tau_b_b->Branch("gen_b_1_eta", &gen_b_1_eta);
	mu_tau_b_b->Branch("gen_b_1_phi", &gen_b_1_phi);
	mu_tau_b_b->Branch("gen_b_1_E", &gen_b_1_E);
	mu_tau_b_b->Branch("gen_diH_pT", &gen_diH_pT);
	mu_tau_b_b->Branch("gen_diH_eta", &gen_diH_eta);
	mu_tau_b_b->Branch("gen_diH_phi", &gen_diH_phi);
	mu_tau_b_b->Branch("gen_diH_E", &gen_diH_E);
	mu_tau_b_b->Branch("gen_diH_mass", &gen_diH_mass);
	mu_tau_b_b->Branch("gen_h_bb_pT", &gen_h_bb_pT);
	mu_tau_b_b->Branch("gen_h_bb_eta", &gen_h_bb_eta);
	mu_tau_b_b->Branch("gen_h_bb_phi", &gen_h_bb_phi);
	mu_tau_b_b->Branch("gen_h_bb_E", &gen_h_bb_E);
	mu_tau_b_b->Branch("gen_h_tt_pT", &gen_h_tt_pT);
	mu_tau_b_b->Branch("gen_h_tt_eta", &gen_h_tt_eta);
	mu_tau_b_b->Branch("gen_h_tt_phi", &gen_h_tt_phi);
	mu_tau_b_b->Branch("gen_h_tt_E", &gen_h_tt_E);
	mu_tau_b_b->Branch("gen_mctMatch", &gen_mctMatch);
	mu_tau_b_b->Branch("gen_weight", &weight);
	//___________________________________________

	// loop over all input files:
	int ievt=0;
	for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile){
		std::string filename = inputFiles_[iFile];
		TFile* inFile = TFile::Open(filename.c_str());

		if(monitorFile_ != "none") {
			TFile *monitorfile = new TFile(monitorFile_.c_str());
			TH1F *tmphisto = new TH1F("tmphisto","" , 10, 0, 10);
			// fill monitor histogram:	
			tmphisto = (TH1F*)monitorfile->Get("mon1/yield");
			double normpar = 2402;
			h_monitor->SetBinContent(1, tmphisto->GetEntries()/normpar);

			tmphisto = (TH1F*)monitorfile->Get("mon2/yield");
			h_monitor->SetBinContent(2, tmphisto->GetEntries()/normpar);

			tmphisto = (TH1F*)monitorfile->Get("mon3/yield");
			h_monitor->SetBinContent(3, tmphisto->GetEntries()/normpar);

			tmphisto = (TH1F*)monitorfile->Get("mon4/yield");
			h_monitor->SetBinContent(4, tmphisto->GetEntries()/normpar);

			tmphisto = (TH1F*)monitorfile->Get("mon5/yield");
			h_monitor->SetBinContent(5, tmphisto->GetEntries()/normpar);

			tmphisto = (TH1F*)monitorfile->Get("mon6/yield");
			h_monitor->SetBinContent(6, tmphisto->GetEntries()/normpar);

			monitorfile->Close();
		}

		if(filename.find("Data") != std::string::npos || filename.find("Run201") != std::string::npos ) runOnData = true;
		if(filename.find("InvIso") != std::string::npos) invIso = true;
		if(filename.find("ToHHTo2B2Tau") != std::string::npos) runOnSignal = true;
		if(filename.find("LS_") != std::string::npos) LS = true;

		if( inFile ){

			// get luminosity weights:
			if(!runOnData){
				// get cross section from GenInfo:
				TTree *runs = (TTree*)inFile->Get("Runs");
				runs->GetEntry();
				GenRunInfoProduct myRunInfoProduct;
				TBranch* genInfoBranch = runs->GetBranch("GenRunInfoProduct_generator__SIM.obj");

				if(genInfoBranch){
					genInfoBranch->SetAddress(&myRunInfoProduct);
					genInfoBranch->GetEntry(0);
					xsec = (double)myRunInfoProduct.crossSection();
				}

				// loop over all luminosity blocks and get number of events:
				fwlite::LuminosityBlock lumiBlock(inFile);
				int ilumiBlock = 0;
				for(lumiBlock.toBegin(); !lumiBlock.atEnd(); ++lumiBlock, ++ilumiBlock){
					edm::Handle<edm::MergeableCounter> hCnt;
					lumiBlock.getByLabel(std::string("nEventsTotal"), hCnt);
					nev_total = nev_total + hCnt->value;

					edm::Handle<edm::MergeableCounter> hCntSel;
					lumiBlock.getByLabel(std::string("nEventsSelected"), hCntSel);
					nev_selected = nev_selected + hCntSel->value;
				}

				if(xsec>0) weight_lumi = xsec / nev_total;

				// for now, hardcoded xsections for all MC:

				if(filename.find("GluGluToRadionToHHTo2B2Tau") != std::string::npos){
						weight_lumi = 100*0.58*0.063*2 / nev_total;   // HARD CODED WEIGTHS FOR MCSIGNAL //*0.58*0.063
					} 

					if(filename.find("TT_TuneCUETP8M1") != std::string::npos) {
						weight_lumi = 831.76 / nev_total;
					}
					if(filename.find("WJetsToLNu") != std::string::npos) {
						weight_lumi = 61526.7 / nev_total;
					}

					if(filename.find("GluGluToRadionToHHTo2B2Tau") == std::string::npos){
						weight_lumi *= 0.55;
					} 

					weight_lumi *= 12.876*1000;

					if(debug){
						std::cout << " xsec = " << xsec << "\n";
						std::cout << " nev_total = " << nev_total << "\n";
						std::cout << " nev_selected = " << nev_selected << "\n";
						std::cout << " weight_lumi = " << weight_lumi << "\n";
					}
				}

			// loop over all events:
				fwlite::Event ev(inFile);
				for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievt){

					try{

						edm::EventBase const & event = ev;
						if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
						if(outputEvery_!=0 ? (ievt>0 && ievt%outputEvery_==0) : false)
							std::cout << "  processing event: " << ievt << std::endl;

						// get all handles present in both data and MC:
						edm::Handle<std::vector<Muon>> selectedmuons;
						try {
							event.getByLabel(edm::InputTag("selectedMuons"), selectedmuons);
						} catch (cms::Exception& iException) {
							std::cout << "  critical: no selectedMuons in event! " << std::endl;
							continue;
						}

						edm::Handle<std::vector<pat::Tau>> selectedtaus;
						if(!invIso) event.getByLabel(edm::InputTag("selectedTaus"), selectedtaus);
						else event.getByLabel(edm::InputTag("selectedTausinvertedIso"), selectedtaus);

						edm::Handle<std::vector<pat::Jet>> selectedjets;
						if(!invIso) event.getByLabel(edm::InputTag("selectedJets"), selectedjets);
						else event.getByLabel(edm::InputTag("selectedJetsInvTauIso"), selectedjets);

						edm::Handle<std::vector<pat::Jet>> updatedPatJetsUpdatedJEC;
						event.getByLabel(edm::InputTag("updatedPatJetsUpdatedJEC", ""), updatedPatJetsUpdatedJEC);

						edm::Handle<std::vector<pat::MET>> slimmedMET;
						event.getByLabel(edm::InputTag("slimmedMETs","","SKIM"), slimmedMET);

						edm::Handle<edm::ValueMap<float>> offlineSlimmedPrimaryVertices;
						event.getByLabel(edm::InputTag("offlineSlimmedPrimaryVertices"), offlineSlimmedPrimaryVertices);

						edm::Handle<double> bTaggingSF;
						if(!invIso){
							if(!LS) event.getByLabel(edm::InputTag("bTaggingSF", "BTaggingSF"), bTaggingSF);
							else event.getByLabel(edm::InputTag("bTaggingSFLSIso", "BTaggingSF"), bTaggingSF);
						} else {
							if(!LS) event.getByLabel(edm::InputTag("bTaggingSFOSInvIso", "BTaggingSF"), bTaggingSF);
							else event.getByLabel(edm::InputTag("bTaggingSFLSInvIso", "BTaggingSF"), bTaggingSF);
						}

						edm::Handle<double> METSignificance;
						event.getByLabel(edm::InputTag("METSignificance"), METSignificance);

						edm::Handle<double> kinfit_p;
						edm::Handle<double> kinfit_chi2;
						edm::Handle<double> kinfit_mH;
						edm::Handle<double> kinfit_convergence;

						if(!invIso){
							event.getByLabel(edm::InputTag("kinfit", "P"), kinfit_p);
							event.getByLabel(edm::InputTag("kinfit", "chi2"), kinfit_chi2);
							event.getByLabel(edm::InputTag("kinfit", "mH"), kinfit_mH);
							event.getByLabel(edm::InputTag("kinfit", "convergence"), kinfit_convergence);
						} else {
							event.getByLabel(edm::InputTag("kinfitInvTauIso", "P"), kinfit_p);
							event.getByLabel(edm::InputTag("kinfitInvTauIso", "chi2"), kinfit_chi2);
							event.getByLabel(edm::InputTag("kinfitInvTauIso", "mH"), kinfit_mH);
							event.getByLabel(edm::InputTag("kinfitInvTauIso", "convergence"), kinfit_convergence);
						}

						edm::Handle<ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2>>> covMatrixHandle;
						event.getByLabel(edm::InputTag("METSignificance", "METCovariance"), covMatrixHandle);

						FitProb = *kinfit_p.product();
						edm::Handle<double> svfitmassHandle;
						if (!invIso) event.getByLabel(edm::InputTag("SVFit", ""), svfitmassHandle);
						else event.getByLabel(edm::InputTag("SVFitInvTauIso", ""), svfitmassHandle);

					// get muon weights:
						edm::Handle<double> muon_triggerSF;
						edm::Handle<double> muon_IDSF;
						event.getByLabel(edm::InputTag("muonTriggerScaleFactor", "ScaleFactorValue"), muon_triggerSF);
						weight_muon_triggerScaleFactor = *muon_triggerSF.product();
						event.getByLabel(edm::InputTag("muonIDScaleFactor", "ScaleFactorValue"), muon_IDSF);
						weight_muon_IdScaleFactor = *muon_IDSF.product();

					// get specific handles for MC:
						if(!runOnData){
							edm::Handle<double> PUweight;
							event.getByLabel(edm::InputTag("PUWeightProducer", "PUWeight"), PUweight);
							weight_PU = *PUweight.product();
							if(debug) std::cout << " weight_PU = " << weight_PU << std::endl;
						}

					// read product with handles:
						try {
							svfitmass = *svfitmassHandle.product();
							if(debug) std::cout << " svfitmass = " << svfitmass << std::endl;
						} catch (cms::Exception& iException) {
							svfitmass = -1;
						}

					// get btag weights:
						try {
							weight_btag = *bTaggingSF.product();
							if(debug) std::cout << " bTaggingSF = " << weight_btag << std::endl;
						} catch (cms::Exception& iException) {
							weight_btag = -1;
						}

						if(debug) std::cout << " b-tagging weight = " << weight_btag << std::endl;
						if(debug) std::cout << " weight_muon_IdScaleFactor * weight_btag = " << weight_muon_IdScaleFactor * weight_btag << std::endl;

					// construct weighting factor, don't actually apply negative weights:
						weight_noPU = weight_lumi * weight_muon_IdScaleFactor * weight_muon_triggerScaleFactor;
						if(weight_btag>0) weight_noPU *= weight_btag;
						weight = weight_noPU * weight_PU;

						if(debug) std::cout << " weight = " << weight << std::endl;

					//Getting objects
						pat::Muon muon = selectedmuons->at(0);
						pat::Tau  tau = selectedtaus->at(0);
						std::pair<int, int> jets = getJets(selectedjets);
						pat::Jet  bjet1 = selectedjets->at(jets.first);
						pat::Jet  bjet2 = selectedjets->at(jets.second);
						pat::MET  met = slimmedMET->at(0);

						muon_p4 = muon.p4();
						tau_p4 = tau.p4();
						bjet0_p4 = bjet1.p4();
						bjet1_p4 = bjet2.p4();
						met_p4 = met.p4();
						hbb_p4 = bjet0_p4+bjet1_p4;
						htt_p4 = tau_p4+muon_p4+met_p4;
						hh_p4 = hbb_p4+htt_p4;

					//Observables
						m_t_mu = sqrt(2 * muon_p4.Pt() * met_p4.Pt() * (1-cos(muon_p4.Phi()-met_p4.Phi())) );


					// fill histograms:
						for(std::vector<Muon>::const_iterator mu1=selectedmuons->begin(); mu1!=selectedmuons->end(); ++mu1){
							hmuonPt->Fill( mu1->pt (), weight );
							hmuonEta->Fill( mu1->eta(), weight );
							hmuonPhi->Fill( mu1->phi(), weight );
						}

						h_m_bbar->Fill(hbb_p4.M(), weight);  
						h_m_tauvis->Fill(htt_p4.M(), weight);
						h_m_t_mu->Fill(m_t_mu, weight);    

					// fill number of vertices:
						double nVert = offlineSlimmedPrimaryVertices->size();
						hnVert->Fill( nVert, weight );
						hnVert_noPU->Fill( nVert, weight_noPU );

						if (svfitmass != -1) hSVFit->Fill(svfitmass, weight);

					// plot Fitprob
						m_kinfit = *kinfit_mH.product();
						h_m_kinfit->Fill(m_kinfit, weight);
						h_P_fit->Fill(FitProb, weight);

					// fill weight histograms:
						hLumiWeight->Fill(weight_lumi);
						hPUWeight->Fill(weight_PU);
						hMuonIdScaleFactor->Fill(weight_muon_IdScaleFactor);
						hMuonTriggerScaleFactor->Fill(weight_muon_triggerScaleFactor);
						hBTaggingSF->Fill(weight_btag);

					// fill covmatrix histos:
						ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2>> mymatrix;
						mymatrix = *covMatrixHandle.product();
						h_cov_xx->Fill( sqrt(mymatrix[0][0]) );
						h_cov_xy->Fill( mymatrix[0][1] );
						h_cov_yx->Fill( mymatrix[1][0] );
						h_cov_yy->Fill( sqrt(mymatrix[1][1]) );

						h_met->Fill(met_p4.Pt(), weight);

						//MC truth_____________________
						//Reset variables______________
						gen_t_0_pT = 0;
						gen_t_0_eta = 0;
						gen_t_0_phi = 0;
						gen_t_0_E = 0;
						gen_t_1_pT = 0;
						gen_t_1_eta = 0;
						gen_t_1_phi = 0;
						gen_t_1_E = 0;
						gen_b_0_pT = 0;
						gen_b_0_eta = 0;
						gen_b_0_phi = 0;
						gen_b_0_E = 0;
						gen_b_1_pT = 0;
						gen_b_1_eta = 0;
						gen_b_1_phi = 0;
						gen_b_1_E = 0;
						gen_diH_pT = 0;
						gen_diH_eta = 0;
						gen_diH_phi = 0;
						gen_diH_E = 0;
						gen_diH_mass = 0;
						gen_h_bb_pT = 0;
						gen_h_bb_eta = 0;
						gen_h_bb_phi = 0;
						gen_h_bb_E = 0;
						gen_h_tt_pT = 0;
						gen_h_tt_eta = 0;
						gen_h_tt_phi = 0;
						gen_h_tt_E = 0;
						//_____________________________
						if (runOnSignal) {
							//Get gen info______________
							edm::Handle<reco::GenParticleCollection> genParticles;
							event.getByLabel(edm::InputTag("prunedGenParticles"), genParticles);
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
								gen_mctMatch = truthFlag(genParticles, mcCuts,//Checks final-state selection was correct
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
								gen_t_0_pT = gen_tau0_p4.Pt();
								gen_t_0_eta = gen_tau0_p4.Eta();
								gen_t_0_phi = gen_tau0_p4.Phi();
								gen_t_0_E = gen_tau0_p4.E();
								gen_t_1_pT = gen_tau1_p4.Pt();
								gen_t_1_eta = gen_tau1_p4.Eta();
								gen_t_1_phi = gen_tau1_p4.Phi();
								gen_t_1_E = gen_tau1_p4.E();
								gen_b_0_pT = gen_bjet0_p4.Pt();
								gen_b_0_eta = gen_bjet0_p4.Eta();
								gen_b_0_phi = gen_bjet0_p4.Phi();
								gen_b_0_E = gen_bjet0_p4.E();
								gen_b_1_pT = gen_bjet1_p4.Pt();
								gen_b_1_eta = gen_bjet1_p4.Eta();
								gen_b_1_phi = gen_bjet1_p4.Phi();
								gen_b_1_E = gen_bjet1_p4.E();
								gen_diH_pT = gen_hh_p4.Pt();
								gen_diH_eta = gen_hh_p4.Eta();
								gen_diH_phi = gen_hh_p4.Phi();
								gen_diH_E = gen_hh_p4.E();
								gen_diH_mass = gen_hh_p4.M();
								gen_h_bb_pT = gen_hbb_p4.Pt();
								gen_h_bb_eta = gen_hbb_p4.Eta();
								gen_h_bb_phi = gen_hbb_p4.Phi();
								gen_h_bb_E = gen_hbb_p4.E();
								gen_h_tt_pT = gen_htt_p4.Pt();
								gen_h_tt_eta = gen_htt_p4.Eta();
								gen_h_tt_phi = gen_htt_p4.Phi();
								gen_h_tt_E = gen_htt_p4.E();
								//__________________________
							}
						}
						//_____________________________
						//_____________________________


					//Fill flat data__________________
					//Reco FS info____________________
						t_0_pT = tau_p4.Pt();
						t_0_eta = tau_p4.Eta();
						t_0_phi = tau_p4.Phi();
						t_0_mass = tau_p4.M();
						t_1_pT = muon_p4.Pt();
						t_1_eta = muon_p4.Eta();
						t_1_phi = muon_p4.Phi();
						t_1_mass = muMass;
						b_0_pT = bjet0_p4.Pt();
						b_0_eta = bjet0_p4.Eta();
						b_0_phi = bjet0_p4.Phi();
						b_0_mass = bjet0_p4.M();
						b_1_pT = bjet1_p4.Pt();
						b_1_eta = bjet1_p4.Eta();
						b_1_phi = bjet1_p4.Phi();
						b_1_mass = bjet1_p4.M();
						mPT_pT = met_p4.Pt();
						mPT_phi = met_p4.Phi();
						h_tt_pT = htt_p4.Pt();
						h_tt_eta = htt_p4.Eta();
						h_tt_phi = htt_p4.Phi();
						h_tt_mass = htt_p4.M();
						h_tt_svFit_mass = svfitmass;
						h_bb_pT = hbb_p4.Pt();
						h_bb_eta = hbb_p4.Eta();
						h_bb_phi = hbb_p4.Phi();
						h_bb_mass = hbb_p4.M();
						diH_pT = hh_p4.Pt();
						diH_eta = hh_p4.Eta();
						diH_phi = hh_p4.Phi();
						diH_mass = hh_p4.M();
						diH_kinFit_mass = m_kinfit;
						diH_kinFit_prob = FitProb;
						mT = m_t_mu;
					//________________________________
					//Shapes__________________________
						getGlobalEventInfo(&tau_p4, &muon_p4, &bjet0_p4, &bjet0_p4, &met_p4,
							&hT, &sT, &centrality, &eVis);
						getPrimaryEventShapes(&tau_p4, &muon_p4, &bjet0_p4, &bjet0_p4,
							&sphericity, &spherocity,
							&aplanarity, &aplanority,
							&upsilon, &dShape,
							&sphericityEigen0, &sphericityEigen1, &sphericityEigen2,
							&spherocityEigen0, &spherocityEigen1, &spherocityEigen2);
					//________________________________
					//________________________________
						mu_tau_b_b->Fill();
					//________________________________

					} catch (cms::Exception& iException) {
						std::cout << iException << "\n";
						std::cout << " skipping single event\n\n";
						continue;
					}

				}

				inFile->Close();
			}
			if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
		}
		return 0;
	}
