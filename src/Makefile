CFLAGS= -I/usr/include/ -O3 -Wall
CXXFLAGS=$(CFLAGS)
LDFLAGS=-L/usr/lib/x86_64-linux-gnu/boost/

all:	compute_distance get_TrieArray compute_distances_ivac debug generate_Peptides

compute_distance.o: compute_distance.C Matrix.h Sequences.h TrieArray.h Trie.h 

compute_distance:	compute_distance.o
	$(CXX) $(LDFLAGS) -o compute_distance $< -lboost_serialization

compute_distances_ivac.o: compute_distances_ivac.C Matrix.h Sequences.h TrieArray.h Trie.h 

compute_distances_ivac:	compute_distances_ivac.o
	$(CXX) $(LDFLAGS) -o compute_distances_ivac  $< -lboost_serialization

get_TrieArray.o: get_TrieArray.C Matrix.h Sequences.h TrieArray.h Trie.h 	

get_TrieArray: get_TrieArray.o
	$(CXX)  $(LDFLAGS) -o get_TrieArray $< -lboost_serialization

debug.o: debug.C Matrix.h Sequences.h TrieArray.h Trie.h

debug: debug.o
	$(CXX)  $(LDFLAGS) -o debug $< -lboost_serialization
	
generate_Peptides.o: generate_Peptides.C Matrix.h Sequences.h TrieArray.h Trie.h

generate_Peptides: generate_Peptides.o
	$(CXX)  $(LDFLAGS) -o generate_Peptides $< -lboost_serialization

