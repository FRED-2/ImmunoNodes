#include "TrieArray.h"
#include "Sequences.h"
#include "Matrix.h"
#include <iterator>

int main(int argc, char** argv)
{
	if (argc < 2)
	{
		std::cout << "Usage: debug <triefile> " << std::endl;
		return -1;
	}

	string triefile(argv[1]);

	TrieArray ta;

  	std::ofstream ofs(triefile.c_str());

	boost::archive::text_oarchive oa(ofs);

	ta.load(oa, 1);

	cout << "Done." << endl;
}