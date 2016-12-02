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
			if (sub.size() != 9)
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
	if (argc < 2)
	{
		std::cout << "Usage: get_TrieArray <fastafile> <outfile> " << std::endl;
		return -1;
	}

	ofstream myfile;
	string fastafile(argv[1]);
	const char* outname(argv[2]);
	
	myfile.open(outname);

	std::cout << fastafile << "\t" << outname << std::endl;
		
//----------------------------------------------------------------------------------------
	cout << "Reading FASTA file..." << endl;

	Sequences s(fastafile);
	cout << "Read " << s.size() << " sequences." << endl;

	cout << "Generating ninemers..." << endl;
	Sequences ninemers;
	generateAllSubstrings(ninemers, s, 9);
	cout << "Generated " << ninemers.size() << " ninemers." << endl;

	s.clear();


  	for (size_t i = 0; i < ninemers.size(); ++i) {
  		myfile << "> " << i << endl;
  		myfile << ninemers[i] << endl; 
  	} 

  	myfile.close();

	cout << "Done." << endl;


}


