#ifndef INCLUDED_EXTRAS
#define INCLUDED_EXTRAS 1
#include "graph/Molecule.hpp"
namespace specific
{
	class ExtendableForceField : public ForceField
	{
	public:
		ExtendableForceField(const ForceField& ff) : ForceField(ff) {}

		inline void addExtraAtomType(Atom* a) {atom_types.push_back(a);}
		inline void addLengthParameters(const std::string& a1,const std::string& a2,double l,double k) {length_parameters.push_back(BondData(a1,a2,l,k));}
		inline void addLengthParameters(const Atom* a1,const Atom* a2,double l,double k) {addLengthParameters(a1->toString(),a2->toString(),l,k);}
		inline void addAngleParameters(const std::string& a1,const std::string& a2,const std::string& a3,double eq,double kt,double ub1=0,double ub2=0)
		{
			std::vector<std::string> atl(1,a1); atl.push_back(a2); atl.push_back(a3);
			std::vector<double> vall; vall.push_back(eq); vall.push_back(kt); vall.push_back(ub1); vall.push_back(ub2);
			angle_parameters.push_back(make_pair(atl,vall));
		}
		inline void addAngleParameters(const Atom* a1,const Atom* a2,const Atom* a3,double eq,double kt,double ub1=0,double ub2=0) {addAngleParameters(a1->toString(),a2->toString(),a3->toString(),eq,kt,ub1,ub2);}

		inline void addDihedralParameters(const std::string& a1,const std::string& a2,const std::string& a3,const std::string& a4,double eq,double k,double mult)
		{
			std::vector<std::string> atl(1,a1); atl.push_back(a2); atl.push_back(a3); atl.push_back(a4);
			std::vector<double> vall; vall.push_back(eq); vall.push_back(k); vall.push_back(mult);
			dihedral_parameters.push_back(make_pair(atl,vall));
		}
		inline void addDihedralParameters(const Atom* a1,const Atom* a2,const Atom* a3,const Atom* a4,double eq,double k,double mult) {addDihedralParameters(a1->toString(),a2->toString(),a3->toString(),a4->toString(),eq,k,mult);}

		inline void addAllowedBonds(const std::string& ref,const std::vector<BondData>& list) {atombondmap.push_back(make_pair(ref,list));}
	};
	class ConstrictedForceField : public ExtendableForceField
	{
		Molecule* refer=nullptr;
		const ForceField* referff;
		const Category* aromcat=nullptr;
		const Category* halocat=nullptr;
		const Category* ringcat=nullptr;
		std::map<std::string,std::string> nameref;

	public:
		std::vector<Atom*> seedAtoms;
	public:
		ConstrictedForceField(Molecule* src,const ForceField& srcff,bool silent=false) : ExtendableForceField(ForceField())
		{
			refer=new Molecule(src);
			referff=&srcff;
			//cout <<"Loading category: " << srcff.getCategory("ring").toString()<<"\n";
			//cout << "Has category: " << srcff.hasCategory("ring") <<"\n";
			if(srcff.hasCategory("aromatic")) aromcat=&(srcff.getCategory("aromatic"));
			if(srcff.hasCategory("halogen")) halocat=&(srcff.getCategory("halogen"));
			if(srcff.hasCategory("ring")) ringcat=&(srcff.getCategory("ring"));
			cout << "Ring category: "<< ringcat->toString()<< "\n";
			seedAtoms=buildForceField(silent);
		}

		const std::vector<Atom*>& getSeedAtoms() const {return seedAtoms;}
		std::vector<Atom*> getSeedAtoms() {return seedAtoms;}

		std::vector<Atom*> buildForceField(bool silent=false)
		{
			int refs[256]; for(int i=0;i<256;i++) refs[i]=0;
			std::map<long int,std::string> oldnames;
			addAngleParameters("C","C","C",109.47,410.0,0,0);
			addDihedralParameters("C","C","C","C",0,5.384,6);
			for(Atom* a : refer->getAtoms())
			{
				char el=a->toString()[0];
				oldnames[(long)a]=a->toString();
				if(halocat && halocat->contains(a))
				{
					if(el=='B') {refs['A']++; a->rename("BR"+to_string(refs['A'])); }
					else if(el=='C') {refs['D']++; a->rename("CL"+to_string(refs['D']));}
				}
				else {refs[el]++; a->rename(std::string(1,el)+to_string(refs[el]));}
				nameref[a->toString()]=oldnames[(long)a];
			}
			refer->dumpMol("templateDock.pdb");
			std::vector<Atom*> seeds;
			rules=std::vector<std::pair<std::string,std::vector<chemtools::Rule>>>();
			std::string aromlist="";
			for(Atom* rA : refer->getAtoms())
			{
				Atom* eA=new Atom(*rA);
				eA->setCyclized(false); eA->setValency(eA->getStandardValency()); eA->setEnergyContribution(0); eA->setHotspotEnergyContribution(0);
				if(aromcat && aromcat->contains(oldnames[(long)rA]))
				{
					if(aromlist!="") aromlist+="|";
					aromlist+=eA->toString();
				}
				eA->setPosition(Eigen::Vector3d(0,0,0));
				if(ringcat && ringcat->contains(oldnames[(long)rA])) eA->mustCycl=true;
				addExtraAtomType(eA);
				if(!eA->isHydrogen()) seeds.push_back(eA);
				std::vector<Atom*> bats=refer->getBondedAtoms(rA);
				std::string rulestring="";
				std::vector<BondData> allbonds;
				for(Atom* fA : bats)
				{
					if(!fA->isHydrogen())
					{
						if(rulestring!="") rulestring+=";";
						rulestring+=fA->toString()+" 1 1";
					}
					#ifndef QUIET
					if(!silent) cout << rA->toString()<<"-"<<fA->toString()<<"\t"<<rA->distanceFrom(fA)<<"\n";
					#endif
					double kb=get<1>(referff->getLengthParameters(oldnames[(long)rA],oldnames[(long)fA]));
					addLengthParameters(rA,fA,rA->distanceFrom(fA),kb);
					allbonds.push_back(BondData(rA->toString(),fA->toString(),rA->distanceFrom(fA),kb));
					std::vector<Atom*> aats=refer->getBondedAtoms(fA);
					for(Atom* sA : aats)
					{
						if(sA==rA) continue;
						std::vector<double> apars=referff->getAngleParamtersFailsafe(oldnames[(long)rA],oldnames[(long)fA],oldnames[(long)sA]);
						#ifndef QUIET
						if(!silent) cout << rA->toString()<<"-"<<fA->toString()<<"-"<<sA->toString()<<"\t"<<apars[0]<<"\n";
						#endif
						addAngleParameters(rA,fA,sA,apars[0],apars[1],apars[2],apars[3]);
						std::vector<Atom*> dats=refer->getBondedAtoms(sA);
						for(Atom* tA : dats)
						{
							if(tA==rA || tA==fA) continue;
							std::vector<double> dpars=referff->getDihedralParametersFailsafe(oldnames[(long)rA],oldnames[(long)fA],oldnames[(long)sA],oldnames[(long)tA]);
							#ifndef QUIET
							if(!silent) cout << rA->toString()<<"-"<<fA->toString()<<"-"<<sA->toString()<<"-"<<tA->toString()<<"\t"<<dpars[0]<<"\n";
							#endif
							addDihedralParameters(rA,fA,sA,tA,dpars[0],dpars[1],dpars[2]);
						}
					}
				}
				for(Atom* fA : bats)
				{
					if(fA->isHydrogen())
					{
						if(rulestring!="") rulestring+=";";
						rulestring+=fA->toString()+" 1 1";
					}
				}
				cout << rA->toString()<<" "<<rulestring<<"\t"<<allbonds.size()<<"\n";
				addAllowedBonds(rA->toString(),allbonds);
				rules.push_back(make_pair(rA->toString(),chemtools::buildRules(rA->toString()+" "+rulestring,referff)));
			}
			if(aromlist!="")
			{
				Category narom("aromatic",aromlist,this);
				categories.push_back(narom);
				#ifndef QUIET
				if(!silent) cout << "Aromatic atom-types detected: " << aromlist << "\n";
				#endif
			}
			#ifndef QUIET
			if(!silent)
			{
				cout << "Added extra parameters\n";
				cout << "Will attempt to add definitions\n";
			}
			#endif

			Topology temptop(refer,*this);
			cout << "Topology built\n";
			#ifndef QUIET
				cout << "Created seeds:\n";
				for(Atom* sa : seeds) cout <<" "<< sa->toString()<<" "<<sa->mustCycl<<"\n";
			#endif
			return seeds;
		}

		inline void renameMolecule(Molecule* mol) {for(Atom* na : mol->getAtoms()) na->rename(nameref[na->toString()]);}
	};

	static std::vector<Molecule*> getSeedsInRegion(quickgeom::Container* reg,Protein& p, std::vector<Atom*> seeds, int nseed,double threshold=1)
	{
		//if(ff.hasCategory("ring")) cyccat=&(ff.getCategory("ring"));
		std::vector<Molecule*> ret;
		Molecule* fullmol=new Molecule(0);
		long int K=0;
		while(ret.size()<nseed)
		{
			Eigen::Vector3d vec=getRandomPositionInSpace(reg);
			Atom* selat=randomSelect(seeds);
			Atom* dum=new Atom(*selat);
			dum->setPosition(vec);
			bool OK=false;
			double eva=get<0>(p.calculateNonBondingEnergy(dum));
			if(eva>-threshold) continue;
			for(Atom* pa : p.getAtoms())
			{
				if(pa->isHydrogen()) continue;
				if(dum->distanceFrom(pa)<RESTRAINDISTANCE/2) {OK=true;}
			}
			K++;
			if(K%100==0) cout << K << "\t"<<ret.size()<<"\n";
			if(OK)
			{
				cout << "Energy: " << eva << "\t";
				Atom* nat=new Atom(*dum); nat->setPosition(dum->getPosition());
				fullmol->addAtom(nat); nat->setEnergyContribution(eva);
				if(selat->mustCycl) nat->mustCycl=true;
				ret.push_back(new Molecule(nat));
				cout << get<0>(p.calculateNonBondingEnergy(ret[ret.size()-1])) << "\t";
				cout << "Accepted\n";
			}
			delete dum;
		}
		fullmol->dumpMol("seeds.pdb");
		cin >> K;
		return ret;
	}

	static std::vector<Molecule*> dockMolecule(Protein& p,Molecule* drug,const ForceField& ff,quickgeom::Container* reg,int poses=80,const std::string& progname="dnv_dock",int seednum=96,int osc=5)
	{
		specific::ConstrictedForceField cff(drug,ff);
		std::vector<Atom*> seeds=cff.getSeedAtoms();
		std::vector<Molecule*> rootseeds=getSeedsInRegion(reg,p,seeds,seednum);
		//generateLigands(ff,n,sizemin,logfilename,rootmols,double rmax=2,double ecut=1,int osccount=5,int tottries=250,int tries_per_lig=45,double temp=298,bool restr=false,bool superseed=false,bool reseed=false,bool deseed=false,const quickgeom::Container* contbox=nullptr,Molecule* initmol=nullptr,bool prune=false,int prunedepth=3)
	  std::vector<System> mols=p.generateLigands(cff,poses,drug->getEffectiveSize()-1,progname+".log",rootseeds,1.75,1,osc,300,50,298,true,false,false,false);
		std::vector<Molecule*> resl;
	  for(System& s : mols)
	  {
	    Molecule* m=s.getDrugMolecule();
	    for(Atom* a : m->getAtoms()) cout <<a->toString()<<"\t"<<cff.isSatisfied(a,m->getBondedAtoms(a))<<","<<a->isCyclized()<<"\n";
			cff.renameMolecule(m);
			resl.push_back(m);
	  }
		return resl;
	}

	class SolvatedProtein : public Protein
	{
		std::vector<Atom*> solvent;
		bool hsopt=false;
	protected:
		SolvatedProtein() {}
	public:
		SolvatedProtein(const Molecule* m) : Protein(m) {cleanSolvent();}
		SolvatedProtein(const std::string& str) : Protein(str.c_str()) {cleanSolvent();}
		SolvatedProtein(const std::string& protfile,const std::string& hostspotfile) : Protein(protfile) {hsopt=true; cleanSolvent();}
	  SolvatedProtein(const std::string& infile,const ForceField& ff,std::vector<int> hotspotreses,const std::string& format="pdb") : Protein(infile,ff,hotspotreses,format) {if(hotspotreses.size()) hsopt=true; cleanSolvent();}

		void cleanSolvent()
		{
			cout << "Solvated Protein loading\n";
			const double THRESH=RESTRAINDISTANCE+0.35;
			Molecule* solvmol=new Molecule(0);
			for(int i=0;i<atoms.size();i++)
			{
				//cout << "'"<<atoms[i]->getResidue()<<"'\n";
				if(atoms[i]->getResidue()=="SOL")
				{
					if(atoms[i]->isHydrogen()) delete atoms[i];
					else
					{
						bool close=false;
						//cout << atoms[i]->toString()<<" "<<atoms[i]->getPosition().transpose()<<"\n";
						for(Atom* hsa : this->getHotspotAtoms())
						{
							if(hsa->distanceFrom(atoms[i])<THRESH) {close=true; break;}
						}
						if(close) {solvent.push_back(atoms[i]); solvmol->addAtom(atoms[i]);}
						else delete atoms[i];
					}
					atoms.erase(atoms.begin()+i);
					i--;
				}
			}
			cout << this->getHotspotAtoms().size() << " hotspot atoms and have used threshold: "<<THRESH<<"\n";
			int binds=0;
			for(Atom* sola : solvent)
			{
				sola->bindpoints=0;
				for(Atom* pa : atoms)
				{
					if(pa->distanceFrom(sola)<pa->seek_sigma()+0.1) {sola->bindpoints++; if(sola->bindpoints>=3) break;}
				}
				binds+=sola->bindpoints;
			}
			cout << "Average bind-points per solvent: "<< ((double)binds)/solvent.size() << "\n";
			cout << "Loaded solvated molecule: "<<solvent.size()<<" solvent atoms kept\n";
			solvmol->dumpMol("solvent.pdb");
		}
		bool hasSolvent() const {return solvent.size();}

		std::pair<float,float> getExposedSolvation(const Molecule* m)
		{
			float ret=0;
			int overlap=0,inter=0;
			float ent=0;
			std::vector<int> bps;
			for(Atom* sol : solvent)
			{
				sol->genpurpflag=true;
				int obp=sol->bindpoints;
				bps.push_back(obp);
				for(Atom* ma : m->getAtoms())
				{
					if(ma->distanceFrom(sol)<ma->seek_radius())
					{
						sol->genpurpflag=false; overlap++;
						ent+=(float)(-WATENT*((float)obp/3));
						break;
					}
				}
			}
			for(Atom* tA : m->getAtoms())
			{
				int nb=0;
				double rad=tA->seek_sigma();
				for(int si=0;si<solvent.size();si++)
				{
					Atom* sol=solvent[si];
					if(!(sol->genpurpflag)) continue;
					if(sol->distanceFrom(tA)<rad) {nb++; if(bps[si]<3) bps[si]++;}
				}
				if(!nb) continue;
				inter+=nb;
				int lim=0;
				switch(tA->getStandardValency())
				{
					case 1:
						lim=5;
						break;
					case 2:
						lim=4;
						break;
					case 3:
						lim=2;
						break;
					case 4:
						lim=1;
						break;
				}
				float frac=(float)nb/lim;
				if(frac>1) frac=1;
				ret+=tA->getSolvationCorrection()*frac;
			}
			cout << overlap<<" overlapping solvent molecules omitted\n";
			cout << "They produce a total entropy difference of: "<<ent<<" kJ/mol\n";
			for(int si=0;si<solvent.size();si++)
			{
				if(!(solvent[si]->genpurpflag)) continue;
				if(bps[si]!=solvent[si]->bindpoints) ent+=WATENT*(((float)(bps[si]-solvent[si]->bindpoints))/3);
			}
			cout << "Finally: Total entropy difference of: "<<ent<<" kJ/mol\n";
			cout << inter << " interaction sites involved finally\n";
			return make_pair(ret,ent);
		}
	};

	std::pair<double,double> stats_mean(const std::vector<double>& v,const std::vector<double>& w=std::vector<double>())
	{
	  if(!v.size()) return make_pair(0,0);
	  double res=0,res2;
	  double den=0;
	  int K=0;
	  for(const double& d : v)
	  {
	    if(w.size()) {res+=d*w[K]; res2+=d*d*w[K]; den+=w[K++];}
	    else {res+=d; res2+=d*d*w[K]; den++;}
	  }
	  double mean=res/den;
		double stdev=::sqrt(res2/den-mean*mean);
		return make_pair(mean,stdev);
	}
	std::vector<double> boltzmannWeights(const std::vector<double>& v)
	{
	  std::vector<double> ret;
	  for(double n : v) ret.push_back(::pow(2.71828,(-::beta*n)/100.00));
	  return ret;
	}

	class ProteinEnsemble : public SolvatedProtein
	{
		std::vector<SolvatedProtein*> ensemble;
		std::pair<double,double> (*met)(const std::vector<double>&,const std::vector<double>&);
		std::vector<double> (*wts)(const std::vector<double>&);
	public:
		ProteinEnsemble(const std::string& files,const ForceField& ff,std::vector<int> hotspotreses,const std::string& format="pdb",std::pair<double,double> (*metric)(const std::vector<double>&,const std::vector<double>&)=stats_mean,std::vector<double> (*wtsf)(const std::vector<double>&)=boltzmannWeights)
		{
			wts=wtsf;
			met=metric;
			std::ifstream filelist; filelist.open(files);
			std::string fname;
			while(true)
			{
				filelist >> fname;
				if(filelist.eof()) break;
				SolvatedProtein* prot=new SolvatedProtein(fname,ff,hotspotreses,format);
				ensemble.push_back(prot);
			}
			if(!ensemble.size()) {cout << "ERR: No proteins in ensemble\n";}
			HN=ensemble[0]->getHotspotAtomCount();
			for(int i=0;i<HN;i++)
			{
				Atom* hsat=new Atom(*(ensemble[0]->getHotspotAtoms()[i]));
				Eigen::Vector3d pos(0,0,0);
				for(SolvatedProtein* prot : ensemble) pos+=prot->getHotspotAtoms()[i]->getPosition();
				pos/=ensemble.size();
				hotspotatoms.push_back(hsat);
			}
			int N=ensemble[0]->getAtoms().size();
			for(int i=0;i<N;i++)
			{
				Atom* nat=new Atom(*(ensemble[0]->getAtoms()[i]));
				Eigen::Vector3d pos(0,0,0);
				for(SolvatedProtein* prot : ensemble) pos+=prot->getAtoms()[i]->getPosition();
				pos/=ensemble.size();
				atoms.push_back(nat);
				bonds.push_back(std::vector<int>());
			}
			cout << ensemble.size()<<" proteins in ensemble\n";
		}

		virtual std::pair<double,double> calculateNonBondingEnergy(Atom* dummy,bool write=false) override
		{
			std::vector<double> l1,l2;
			for(SolvatedProtein* prot : ensemble)
			{
				l1.push_back(get<0>(Molecule::calculateNonBondingEnergy(dummy,write)));
				double thsbe=0;
				for(int i=0;i<prot->getHotspotAtomCount();i++) thsbe+=chemtools::getNonbondingPotential(dummy,prot->getHotspotAtoms()[i]);
				l2.push_back(thsbe);
			}
			std::pair<double,double> meanEn=met(l1,wts(l1)),meanHEn=met(l2,wts(l2));
			if(get<0>(meanEn)!=get<0>(meanEn)) return make_pair(1e10,1e10);
			dummy->setEnergyContribution(get<0>(meanEn)); dummy->setHotspotEnergyContribution(get<1>(meanHEn));
			return make_pair(dummy->getEnergyContribution(),dummy->getHotspotEnergyContribution());
		}
		virtual std::pair<double,double> calculateNonBondingEnergy(Molecule* drug,bool write=false) override
		{
			cout << "Ensemble-based energy calculation\n";
			std::vector<double> l1,l2;
			for(SolvatedProtein* prot : ensemble)
			{
				std::pair<double,double> eners=prot->calculateNonBondingEnergy(drug,write);
				l1.push_back(get<0>(eners)); l2.push_back(get<1>(eners));
			}
			std::pair<double,double> meanEn=met(l1,wts(l1)),meanHEn=met(l2,wts(l2));
			drug->myBE=get<0>(meanEn); drug->errBE=get<1>(meanEn);
			drug->myHBE=get<0>(meanHEn); drug->errHBE=get<1>(meanHEn);
			return make_pair(drug->myBE,drug->myHBE);
		}

		virtual inline int getHotspotAtomCount() const override {return ensemble[0]->getHotspotAtomCount();}
		virtual inline void loadDevations(const std::vector<int>& resnums,const std::vector<double>& devs) override {for(SolvatedProtein* prot : ensemble) prot->loadDevations(resnums,devs);}
	};
}
#endif
