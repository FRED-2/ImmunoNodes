#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <tr1/cstdint>
#include "Trie.h"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

using namespace std;



class TrieArray
{
 	public:
  typedef size_t PackedIndex;
  typedef vector<PackedIndex> LevelArray;
  typedef vector<LevelArray> PackedTrie;

	PackedTrie data;
	TrieArray();
	TrieArray(const Trie& t, size_t size);
 	void add(string s);
	
	template<class Archive>
		void save(Archive & ar, const unsigned int version) const;
		
	template<class Archive>
		void load(Archive & ar, const unsigned int version);
	BOOST_SERIALIZATION_SPLIT_MEMBER()

	private:
	friend class boost::serialization::access;

};
BOOST_CLASS_VERSION(TrieArray, 1)


TrieArray::TrieArray()
	: data ()
	{}


TrieArray::TrieArray(const Trie& t, size_t size)
  {
		vector<Trie::TrieGraph::vertex_descriptor> nodes;
		nodes.push_back(t.root);

		data.push_back(LevelArray(1, 0));
		while (!nodes.empty())
		{
			size_t number_of_nodes = 0;
			for (size_t i = 0; i < nodes.size(); ++i)
			{
				data.back()[i] |= number_of_nodes;
				number_of_nodes += out_degree(nodes[i], t.g);
			}
			data.push_back(LevelArray());
			data.back().resize(number_of_nodes);
			cout << "Layer " << data.size() << " has " << number_of_nodes << " nodes." << endl;
			size_t j = 0;
			vector<Trie::TrieGraph::vertex_descriptor> new_nodes(number_of_nodes);
			for (size_t i = 0; i < nodes.size(); ++i)
			{
				Trie::TrieGraph::adjacency_iterator ai, a_end;
        for (tie(ai, a_end) = adjacent_vertices(nodes[i], t.g); ai != a_end; ++ai, ++j)
				{
					data.back()[j] = get(vertex_index_t(), t.g, *ai) << 56;
					new_nodes[j] = *ai;
				}
			}
			nodes = new_nodes;
		}
  }
  


void TrieArray::add(string s)
  {
    cout << "Adding: " << s << endl;
  }


template<class Archive>
	void TrieArray::save(Archive & ar, const unsigned int version) const
	{		

	    for (size_t i=0; i<data.size(); ++i)
  	{
  						cout << data.size() << endl;
						ar & data[i];
  	} 	
	}

template<class Archive>
	void TrieArray::load(Archive & ar, const unsigned int version)
	{	
		
		for (size_t i=0; i<11; ++i)
		{
			LevelArray v;
			ar & v;
			data.push_back(v);
		}
		
	}
