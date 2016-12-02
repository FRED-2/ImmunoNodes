#include "TrieArray.h"
#include "Sequences.h"
#include "Matrix.h"
#include <algorithm>
#include <set>
#include <utility>



void generateAllSubstrings(Sequences& substrings, const Sequences& s, size_t len)
{
	for (size_t i = 0; i < s.size(); ++i)
	{
		for (int j = 0; j <= ((int)s[i].size() - (int)len); ++j)
		{
			string sub(string(&(s[i][j]), len));
			substrings.push_back(sub);
		}
	}
}


typedef pair<size_t, size_t> Node;
typedef vector<Node> Children;
typedef vector<size_t> Peptide;

//------------------------------------------------------------------------------------------------------------------------
Children getChildren(const TrieArray& ta, Node n, Matrix& m)
{
	Children children;
	if (n.first == ta.data.size()-1) //leaf
	{
		return children;
	}
	else
	{
		size_t v_e;

		size_t v = ta.data[n.first][n.second]; //Inhalt aktueller Knoten
		size_t v_s = (v & 0x00ffffffffffffff); //Zeiger auf erstes Kind bzw. Index des ersten Kindes in nächster Ebene.

		if (n.second >= ta.data[n.first].size()-1) //last sibling
		{
			v_e = ta.data[n.first + 1].size();
		}
		else
		{
	        size_t v_next = ta.data[n.first][n.second+1]; //Inhalt naechster Geschwisterknoten, braucht man für den Zeiger auf Ende der Kinder vom aktuellen Knoten
			v_e = (v_next & 0x00ffffffffffffff); //Zeiger auf letztes Kind bzw. Index des letzten Kindes in nächster Ebene.
		}

		if (v_s > v_e)
		{
			cout << v << " " << v_s << " " << v_e << " n1: " << n.first << " n2: " << n.second << " " << ta.data[n.first].size() << endl;
		}

		for( size_t i = v_s; i < v_e; i++ )
		{
			size_t level = n.first + 1;
			size_t sibling = i;
			Node c (level,sibling);
			children.push_back(c);
		}
		return children;
	}
}


//------------------------------------------------------------------------------------------------------------------------
//double print_peptides(const TrieArray& ta, Node n, double& dist, Matrix& m, Peptide& p, Peptide& s, double& dist_min)
//{	
//		Children c = getChildren(ta, n, m);
//		if (c.size() == 0) //at leaf
//		{
//			if (dist < dist_min)//peptide with new minimal distance found. 
//			{	
//				dist_min = dist;
//			}	
//			return dist_min;
//		}
//		for (size_t i=0; i<c.size(); ++i)
//		{
//		int a = (ta.data[c[i].first][c[i].second] >> 24);
//		int b = s[(c[i].first-1)];
//		p.push_back(a);
//		double d  = m.data_[a*20 + b];
//		if (abs(d)<0.000001)
//			{d = 0.0;}
//		dist += d;
//		cout << i << " " <<"node: " << c[i].first << " " << c[i].second << endl;
//		DFS_BnB(ta,c[i],dist,m,p,s, dist_min);
//		p.pop_back();
//		dist -= d;
//		if (abs(dist)<0.0000001)
//				{dist = 0.0;}
//		if (p.size() == 0)
//			{dist = 0.0;}
//		}
//		return dist_min;
//}



//

double DFS_BnB(const TrieArray& ta, Node n, double& dist, Matrix& m, Peptide& p, Peptide& s, double& dist_min)
/*
Retruns the minimal distance 
Parameters:
TrieArray: The search trie
dist: distance of search peptide s to 
m: Distance matrix 
p: path to the current node in the trie. keeps track of the words stored in the nodes. 
s: search peptide
dist_min: minimal distance foud so far. Required for break in branch and bound. 
*/
{	
	if (dist >= dist_min)//check if the bound criteria is fullfilled (distance is already greater or equal to minimal distancs foud so far)
			{	
				return dist_min;
			}
	
	else//Branch - recursion
	{	
		Children c = getChildren(ta, n, m);
		if (c.size() == 0) //at leaf
		{
			if (dist < dist_min)//peptide with new minimal distance found. 
			{	
				dist_min = dist;
			}	
			return dist_min;
		}
		for (size_t i=0; i<c.size(); ++i)
		{
		int a = (ta.data[c[i].first][c[i].second] >> 24);
		int b = s[(c[i].first-1)];
		p.push_back(a);
//		int d1 = m.data_[a*20 + b];
//		int d2 = m.data_[b*20 + a];
//		double d = (d1 + d2)/2;
		double d  = m.data_[a*20 + b];
		if (abs(d)<0.000001)
			{d = 0.0;}
		dist += d;
//		cout << i << " " <<"node: " << c[i].first << " " << c[i].second << endl;
		DFS_BnB(ta,c[i],dist,m,p,s, dist_min);
		p.pop_back();
		dist -= d;
		if (abs(dist)<0.0000001)
				{dist = 0.0;}
		if (p.size() == 0)
			{dist = 0.0;}
		}
  }

return dist_min;
}

//------------------------------------------------------------------------------------------------------------------------


multiset<pair<double, string> > DFS_BnB_x_pair(const TrieArray& ta, Node n, double& dist, Matrix& m, Peptide& p, Peptide& s, multiset<pair<double, string> >& dist_min, size_t x)
/*
Returns the x clostest peptides
Returns pairs of distance and peptide
Parameters:
TrieArray: The search trie
dist: distance of search peptide s to 
m: Distance matrix //we have to see what we will do with similarity matrices... 
p: path to the current node in the trie. keeps track of the words stored in the nodes. 
s: search peptide
dist_min: minimal distance foud so far. Required for break in branch and bound. 
*/
{	
		if (dist_min.size() >= x) //If not enough distances found so far
		{	
			if (dist > (*--dist_min.end()).first)//check if the bound criteria is fullfilled (distance is already greater or equal to maximal distance selected so far)		
			{	
				//cout << "break: " << *--dist_min.end() << " (" << dist << ") " <<  "-----";
				return dist_min;
			}
		}

	
		Children c = getChildren(ta, n, m);
		if (c.size() == 0) //at leaf
		{	
			string seq;
			for (size_t pos=0; pos<p.size();++pos)
			{
				seq += m.indexToChar(p[pos]); 
			}
			if (dist_min.size() < x)
			{
				pair<double, string> temp_pair(dist,seq);
				dist_min.insert(temp_pair);
				//cout << dist << "\t" << seq << endl;
			}
			else if (dist < (*--dist_min.end()).first)//peptide with new minimal distance found. 
			{	
				pair<double, string> temp_pair(dist,seq);
				//cout << seq << endl;
				dist_min.insert(temp_pair);
				//cout << dist << "\t" << seq << endl;
				dist_min.erase(--dist_min.end()); //remove largest distance from list
				//cout << --dist_min.end() << endl;
			}	
		
			return dist_min;
	}
		
		for (size_t i=0; i<c.size(); ++i)
		{
			size_t a = (ta.data[c[i].first][c[i].second] >> 56);
			size_t b = s[(c[i].first-1)];
			p.push_back(a);
			double d  = m.data_[a*20 + b];
			if (abs(d)<0.0000001)
				{d = 0.0;}
			dist += d;
			multiset<pair<double,string> > dist_min_tmp;
			dist_min_tmp = DFS_BnB_x_pair(ta,c[i],dist,m,p,s, dist_min, x);
			dist_min = dist_min_tmp;
			p.pop_back();
			dist -= d;
			if (abs(dist)<0.0000001)
				{dist = 0.0;}
		}
return dist_min;
}





multiset<double> DFS_BnB_x(const TrieArray& ta, Node n, double& dist, Matrix& m, Peptide& p, Peptide& s, multiset<double>& dist_min, size_t x)
/*
Returns the distances to the x clostest peptides

Parameters:
TrieArray: The search trie
dist: distance of search peptide s to 
m: Distance matrix //we have to see what we will do with similarity matrices... 
p: path to the current node in the trie. keeos track of the words stored in the nodes. 
s: search peptide
dist_min: minimal distance foud so far. Required for break in branch and bound. 
*/
{	
		if (dist_min.size() >= x) //If not enough distances found so far
		{	
			if (dist > *--dist_min.end())//check if the bound criteria is fullfilled (distance is already greater or equal to maximal distance selected so far)		
			{	
	//			cout << "break: " << *--dist_min.end() << " (" << dist << ") " <<  "-----";
				return dist_min;
			}
		}

	
		Children c = getChildren(ta, n, m);
		if (c.size() == 0) //at leaf
		{	
//			string seq;
//			for (size_t pos=0; pos<p.size();++pos)
//			{
//				seq += m.indexToChar(p[pos]); 
//			}
			if (dist_min.size() < x)
			{
				dist_min.insert(dist);
//				cout << dist << "\t" << seq << endl;
			}
			else if (dist < *--dist_min.end())//peptide with new minimal distance found. 
			{	
				dist_min.insert(dist); //add new distance
//				cout << dist << "\t" << seq << endl;
				dist_min.erase(--dist_min.end()); //remove largest distance from list
			}	
		
			return dist_min;
	}
		
		for (size_t i=0; i<c.size(); ++i)
		{
			int a = (ta.data[c[i].first][c[i].second] >> 24);
			int b = s[(c[i].first-1)];
			p.push_back(a);
			double d  = m.data_[a*20 + b];
			if (abs(d)<0.0000001)
				{d = 0.0;}
			dist += d;
			multiset<double> dist_min_tmp;
			dist_min_tmp = DFS_BnB_x(ta,c[i],dist,m,p,s, dist_min, x);
			dist_min = dist_min_tmp;
			p.pop_back();
			dist -= d;
			if (abs(dist)<0.0000001)
				{dist = 0.0;}
		}
return dist_min;
}

multiset<double> DFS_BnB_x_weighted(const TrieArray& ta, Node n, double& dist, Matrix& m, Peptide& p, Peptide& s, multiset<double>& dist_min, size_t x, int positions[])
/*Parameters:
TrieArray: The search trie
dist: distance of search peptide s to 
m: Distance matrix 
p: path to the current node in the trie. keeps track of the words stored in the nodes. 
s: search peptide
dist_min: minimal distance foud so far. Required for break in branch and bound. 
positions: a vector of length 9 with 0/1 indicationg the positions that shall be considered for distance
*/
{	
		if (dist_min.size() >= x) //If not enough distances found so far
		{	
			if (dist > *--dist_min.end())//check if the bound criteria is fullfilled (distance is already greater or equal to maximal distance selected so far)		
			{	
				return dist_min;
			}
		}
		Children c = getChildren(ta, n, m);
		if (c.size() == 0) //at leaf
		{	
//			string seq;
//			for (size_t pos=0; pos<p.size();++pos)
//			{
//				seq += m.indexToChar(p[pos]); 
//			}
			if (dist_min.size() < x)
			{
				dist_min.insert(dist);
//				cout << dist << "\t" << seq << endl;
			}
			else if (dist < *--dist_min.end())//peptide with new minimal distance found. 
			{	
				dist_min.insert(dist); //add new distance
//				cout << dist << "\t" << seq << endl;
				dist_min.erase(--dist_min.end()); //remove largest distance from list
			}	
			return dist_min;
	}
		
		for (size_t i=0; i<c.size(); ++i)
		{
			int a = (ta.data[c[i].first][c[i].second] >> 24);
			int b = s[(c[i].first-1)];
			p.push_back(a);
			double d;
			if (positions[p.size()-1] == 1) //only consider some positions
			{
				d  = m.data_[a*20 + b];
				if (abs(d)<0.0000001)
					{d = 0.0;}
			}
			else //do not consider this position for distance 
			{
				d = 0.0;
			}
			dist += d;
			multiset<double> dist_min_tmp;
			dist_min_tmp = DFS_BnB_x_weighted(ta,c[i],dist,m,p,s, dist_min, x, positions);
			dist_min = dist_min_tmp;
			dist -= d;		
			p.pop_back();
			if (abs(dist)<0.0000001)
				{dist = 0.0;}			
		}
return dist_min;
}













//--------------------------------------------------------------------------------------------------------------------

double find_closest(vector<Peptide>& haystack, Peptide& needle, Matrix& m)
	{
		double dist_min = 100;
		
		for (size_t i=0; i<haystack.size(); ++i )
		{
			double dist = 0.0;
			Peptide pep = haystack[i];
			for (size_t pos=0; pos<pep.size();++pos)
			{
				double d = m.data_[pep[pos]*20 + needle[pos]];
				dist += d; 		
			}	
			if (dist < dist_min){dist_min = dist;}	
		}
		return dist_min;
	}

multiset<double> find_closest_x(vector<Peptide>& haystack, Peptide& needle, Matrix& m, size_t x)
{
		multiset<double> dist_min;
		
		for (size_t i=0; i<haystack.size();++i)
		{			
			double dist = 0.0;
			Peptide pep = haystack[i];
			
			string seq;
			for (size_t pos=0; pos<pep.size();++pos)
			{
				double d = m.data_[pep[pos]*20 + needle[pos]];
				dist += d;
				seq += m.indexToChar(pep[pos]); 			
			}	
			if (dist_min.size() < x)
				{
					dist_min.insert(dist);
				}
			else if (dist < *--dist_min.end() )
			{
				dist_min.insert(dist);
				dist_min.erase(--dist_min.end());
				
			}
		}
		return dist_min;
}

//------------------------------------------------------------------------------------------------------------------------
