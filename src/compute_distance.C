//#include "TrieArray.h"
//#include "Sequences.h"
//#include "Matrix.h"
#include <map>
#include <cstdio>
#include <iterator>
#include <time.h>
#include "TrieFunctions.h"

int main(int argc, char** argv)
{	
	size_t nof_peptides = 10;
	if (argc < 6)
	{
		std::cout << "Usage: compute_distance <matrix> <trie> <input> <output> <affinity_file> [<nof_peptides>]" << std::endl;
		return -1;
	}

	if (argc > 6) // nof_peptides is optional
	{
		nof_peptides = atoi(argv[6]);
	}

	string matrix(argv[1]);
	string trie(argv[2]);
	string filename(argv[3]);
	string outname(argv[4]);
	string affinity_file(argv[5]);
	std::cout << matrix << ", " << trie << ", " << filename << ", " << outname << ", " << nof_peptides << std::endl;

	std::cout << "Reading trie..." << std::endl;
	TrieArray ta;
	{
	std::ifstream ifs(trie.c_str()); //trie is a string containing the path and filename of the trie file.
	boost::archive::text_iarchive ia(ifs);
	ta.load(ia,1);
	}

	
	Matrix m(matrix);	 
//# read in affinities for all human peptides for combined output of distance and affinties.
	
	map<string,string> affinities;
	{
//		string affinity_file = "/share/usr/feldhahn/Projects/immunopeptidomics/proteome_distance/raw_data/human_9mer_binders_A0201_with_affinities.txt";
//		string affinity_file = "/share/usr/feldhahn/Projects/immunopeptidomics/proteome_distance/raw_data/human_9mer_binders_B3501_NETMHC_with_affinities.txt";
		std::cout << "Reading affinites for peptides in trie..." << std::endl;
		ifstream is_aff(affinity_file.c_str());
		if (not is_aff)
			throw "Cannot open affinity File!";
		string line;
		while (getline(is_aff,line))
		{
			std::stringstream   linestream(line);
    			std::string         data;
    			string	peptide;
			string affinity;

                     std::getline(linestream, data, '\t');  // read up-to the first tab (discard tab).
    
                     // Read the integers using the operator >>
                         linestream >> peptide  >> affinity;
			//cout << peptide << ' ' << affinity << endl;
    
			//string::size_type tab_loc = line.find("\t");
			//size_t l = line.length();
			//string peptide = line.substr(0,tab_loc);
			//string affinity = line.substr(tab_loc+1,l);
			//cout << peptide << affinity << endl;
			//string peptide = 
			//string affinity = 
			affinities[peptide] = affinity;
		}
		is_aff.close();
	}


//# read in info on peptides:
//#peptide\treactivity\tNetMHC-score\tvirality\nYVFKRYLSI	0	0.667000	1	

	map<string,string> query_reactivity;
	map<string,string> query_affinity;
	map<string,string> query_virality;
	set<string> peptides;
	
	
	{
//		string peptide_info_file = "/share/usr/feldhahn/Projects/immunopeptidomics/proteome_distance/raw_data/peptides_with_reactivities_affinities_virality.txt";
		std::cout << "Reading search peptides and additional information..." << std::endl;
		
		ifstream is(filename.c_str());
		if (not is)
			throw "Cannot open info File!";
		string line;
		while (getline(is,line))
		{	
			string::size_type comment = line.find("#");
			if (comment == string::npos)
			{	
//				std::cout << line << std::endl;
				string::size_type tab_loc = line.find("\t"); //split after peptide
				if (tab_loc != string::npos)
				{
					size_t l = line.length();	
					string peptide = line.substr(0,tab_loc);
					string rest = line.substr(tab_loc+1,l);
						
					tab_loc = rest.find("\t"); //split after reactivity
					l = rest.length();	
					string reactivity = rest.substr(0,tab_loc);
					rest = rest.substr(tab_loc+1,l);
		
						
					tab_loc = rest.find("\t"); //split after NetMHC score
					l = rest.length();	
					string affinity = rest.substr(0,tab_loc);
					string virality = rest.substr(tab_loc+1,l);
					
					peptides.insert(peptide);
					query_reactivity[peptide] = reactivity;
					query_affinity[peptide] = affinity;
					query_virality[peptide] = virality;
				}
				else
				{
					peptides.insert(line);
				}
			}
		}
		is.close();
	}
	
	
	std::cout << "Computing distances..." << std::endl;
	ofstream os( outname.c_str() );	
	for( set<string>::const_iterator iter = peptides.begin(); iter != peptides.end(); ++iter ) 
	{
    string s = *iter;
    std::cout << s << std::endl;
    //std::cout << query_reactivity[s] << "," << query_affinity[s] << "," << query_virality[s] << std::endl;
    cout << affinities[s] << endl; 
  	
		std::cout << "." ;
		flush(cout);		
		Node n (0,0);
		Peptide p;
		Peptide seq;
		m.translate(s, seq);
		
		multiset<pair<double,string> > dist_min;
		multiset<pair<double,string> > dt;
		double dist = 0.0;
		dist_min = DFS_BnB_x_pair(ta,n,dist,m,p,seq,dt,nof_peptides);	
		//os << s << "," << query_reactivity[s] << "," << query_affinity[s] << "," << query_virality[s] <<":";
		for (multiset<pair<double,string> >::iterator it=dist_min.begin() ; it != dist_min.end(); it++ )
//							{os << (*it).second <<"," << (*it).first << ";";}
			{os << (*it).second <<"," << (*it).first << "," << affinities[(*it).second] << ";";}
		os << endl;	

		

	}
	std::cout << endl;
	os.close();
	return 0;
	
	
	
}

