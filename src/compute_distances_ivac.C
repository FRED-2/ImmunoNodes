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
	if (argc < 4)
	{
		std::cout << "Usage: compute_distance <matrix> <trie> <input> [<nof_peptides>]" << std::endl;
		return -1;
	}
	//cout << argc << endl;
	if (argc > 4) // nof_peptides is optional
	{
		nof_peptides = atoi(argv[4]);
	}

	string matrix(argv[1]);
	string trie(argv[2]);
	string filename(argv[3]);
	//std::cout << matrix << ", " << trie << ", " << filename << ", " << outname << ", " << nof_peptides << std::endl;

	//std::cout << "Reading trie..." << std::endl;
	TrieArray ta;
	{
	std::ifstream ifs(trie.c_str()); //trie is a string containing the path and filename of the trie file.
	boost::archive::text_iarchive ia(ifs);
	ta.load(ia,1);
	}

	
	Matrix m(matrix);	 

	set<string> peptides;
	
	
	// Read petides! One peptide sequence per line
	
	{ 

		//std::cout << "Reading search peptides and additional information from file " <<   std::endl;
		
		ifstream is(filename.c_str());
		if (not is)
			throw "Cannot open info File!";
		string line;
		while (getline(is,line))
		{	
			string::size_type comment = line.find("#");
			if (comment == string::npos)
			{	
						peptides.insert(line);
					//	std::cout << line << std::endl;
			}
		}
		is.close();
	}
	
	
	//std::cout << "Computing distances..." << std::endl;
	
	//ofstream os( outname.c_str() );	
	for( set<string>::const_iterator iter = peptides.begin(); iter != peptides.end(); ++iter ) 
	{
    string s = *iter;
    //std::cout << s << std::endl; 
  	
		//std::cout << "." ;
		flush(cout);		
		Node n (0,0); //start at top of the trie
		Peptide p; //
		Peptide seq;
		//std::cout << s << std::endl;
		m.translate(s, seq); //translate peptide sequence to matrix indices. seq contains the translated peptide sequence. 
		
		multiset<pair<double,string> > dist_min;
		multiset<pair<double,string> > dt;
		double dist = 0.0;
		dist_min = DFS_BnB_x_pair(ta,n,dist,m,p,seq,dt,nof_peptides);	
//		os << s << "," << query_reactivity[s] << "," << query_affinity[s] << "," << query_virality[s] <<":";
		//os << s << ":";
		cout << s << ":";
		for (multiset<pair<double,string> >::iterator it=dist_min.begin() ; it != dist_min.end(); it++ )
							//{os << (*it).second <<"," << (*it).first << ";";}
							{cout << (*it).second <<"," << (*it).first << ";";}
					//cout << (*it).second << (*it).first << endl;}
			//{os << (*it).second <<"," << (*it).first << "," << affinities[(*it).second] << ";";}
		//os << std::endl;	
		cout << std::endl;
		

	}
	//std::cout << std::endl;
//	os.close();
	return 0;
	
	
	
}

