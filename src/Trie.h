#include <boost/graph/adjacency_list.hpp>
using namespace std;
using namespace boost;

class Trie
{
	public:
	typedef property<vertex_index_t, size_t> VertexIndex;
	typedef adjacency_list<listS, listS, directedS, VertexIndex> TrieGraph;
	typedef property_map<TrieGraph, vertex_index_t>::type PropertyMap;
	typedef TrieGraph::adjacency_iterator VertexIterator;
	typedef graph_traits<TrieGraph>::vertex_descriptor Vertex;
	TrieGraph g;
	ofstream myfile;
	PropertyMap map;
	Vertex root;

	Trie() 
	{
		myfile.open("indices.txt");
		root = add_vertex(0, g);
	}

	void dump()
	{
		cout << "Graph has " << num_edges(g) << " edges and " << num_vertices(g) << " vertices." << endl;
		myfile.close();
	}

	void add(const vector<size_t>& sequence)
	{

		 //cout << "Adding sequence ";
		//td::copy(sequence.begin(), sequence.end(), ostream_iterator<int>(myfile, ", "));
		//myfile << sequence.size() << endl;
			 //copy(sequence.begin(), sequence.end(), ostream_iterator<int>(cout, ", ")) << endl;
		 //cout << endl;
		Vertex u = root;
		PropertyMap idx = get(vertex_index_t(), g);
		for (size_t i = 0; i < sequence.size(); ++i)
		{
			graph_traits<TrieGraph>::adjacency_iterator vi, vi_end;
			tie(vi, vi_end) = adjacent_vertices(u, g);
			if (vi != vi_end)
			{  // There are vertices on this level, now we need to figure out
         // whether there is a node labeled with sequence[i] already.
				for (; vi != vi_end; ++vi)
				{
					if (idx[*vi] == sequence[i]) 
					{
						myfile << sequence[i] << ',' ;
						u = *vi;
						break;
					}
				}
				if (vi == vi_end)
				{
					// We didn't find a matching node, add it to the trie.
					Vertex v = add_vertex(sequence[i], g);
					myfile << sequence[i] << ',';
					add_edge(u, v, g);
					u = v;
				}
			}
			else
			{ // No node on this level exists, so we can safely add a node with
        // sequence[i] on this level.
				Vertex v = add_vertex(sequence[i], g);
				myfile << sequence[i] << ',';
				add_edge(u, v, g);
				u = v;
			}
		}

			myfile  << endl;
	}

};
