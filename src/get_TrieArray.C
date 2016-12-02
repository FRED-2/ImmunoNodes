#include "TrieArray.h"
#include "Sequences.h"
#include "Matrix.h"
#include <iterator>

// Takes a fasta file and computes a trie of all 9mer subsequences in the file. 
//
// The matrix (line 63) is needed because the letters are converted into integers. 
// All our matrices have the same ordering of letters. If you use a new matrix, pleas make sure to use the same ordering!! Otherwise the tries have to be recomputed!

void generateAllSubstrings(Sequences& substrings, const Sequences& s, size_t len)
{

  	
	for (size_t i = 0; i < s.size(); ++i)
	{
		for (int j = 0; j <= ((int)s[i].size() - (int)len); ++j)
		{
			string sub(string(&(s[i][j]), len));
			substrings.push_back(sub);
			if (sub.size() != len)
			  cout << "generateAllSubstrings size: " << sub.size() << endl;
		}
	}
}


//---------------------------------------------------------------------------------
//
//typedef pair<size_t, size_t> Node;
//typedef vector<Node> Children;
//typedef vector<int> Peptide;




//---------------------------------------------------------------------------------



int main(int argc, char** argv)
{
	if (argc < 4)
	{
		std::cout << "Usage: get_TrieArray <fastafile> <matrix> <peptideLength> <outfile> " << std::endl;
		return -1;
	}

	string fastafile(argv[1]);
	string matrix(argv[2]);
	cout << "test" << endl;
	int peptideLength(atoi(argv[3]));
	string outname(argv[4]);

	
	//std::cout << fastafile << "\t" << outname << std::endl;
		
//----------------------------------------------------------------------------------------
	cout << "Reading FASTA file..." << endl;

	Sequences s(fastafile);
	cout << "Read " << s.size() << " sequences." << endl;

	cout << "Generating peptides..." << endl;
	Sequences ninemers;
	generateAllSubstrings(ninemers, s, peptideLength);
	cout << "Generated " << ninemers.size() << " peptides." << endl;

	s.clear();


  	//Matrix m("/abi-projects/dist2self/matrices/BLOSUM45_distance_normal.dat"); cout << "Initializing trie. " << endl; Trie t; 
  	Matrix m(matrix); cout << "Initializing trie. " << endl; Trie t;
	Matrix::IndexSequence indices; 
  	for (size_t i = 0; i < ninemers.size(); ++i) {
 
  	 m.translate(ninemers[i], indices); t.add(indices);
  	} 
    t.dump();


	cout << "Converting to trie array." << endl;
  	TrieArray ta(t, peptideLength);

	cout << "Done." << endl;

//	std::ofstream ofs("test.trie");
	std::ofstream ofs(outname.c_str());
	boost::archive::text_oarchive oa(ofs);
	ta.save(oa,1);


}


