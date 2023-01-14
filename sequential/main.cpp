#include <iostream>
#include <fstream>
#include <set>
#include <iterator>
#include <map>
#include <queue>
#include <vector>
#include <string>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <cstdio>
#include <map>
#include <algorithm>
#include <cmath>
#include <omp.h>
#include <chrono>
#include <mpi.h>
#include <metis.h>

using namespace std;

/*
rth row and jth column (r and j starting from 0)
a. increase row_begin(r+1) by 1
b. sum from 0 to r-1. values[sum+1] is the new node's position
c. col_indices[sum+1] is its equivalent column reference

3 r 4 c
in-degree --> row sum
oudegree --> column sum

If there is an edge from k to l, the kth column lth row is 1

*/

unordered_map<int, vector<int>> all_edges_new;
#define FILE "graph.txt"
#define DEBUG_FILE_READ 0
#define CSR_PRINT 0

class CSR{
    public:
        vector<int> row_begin;
        CSR(int nodeNum);
        double* values;
        vector<int> col_indices;
};

CSR::CSR(int nodeNum){
    this->row_begin = vector<int>(nodeNum+1);
    this->col_indices = vector<int>();
}

CSR* csr;


int main(int arg, char* argv[])
{
    // For MPI
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    // Inıtialize the variables
    string first,second;
    int index = 0;
    // This holdes the all edges as source and destination
    vector<pair<int, int> > all_edges; 
    //This holds the all nodes as name and index
    unordered_map<string, int> dictionary;
    unordered_map<int, vector<int> >  all_edges_new;
    unordered_map<int, int> out_arrows;

    freopen(FILE, "r", stdin);  // read from file
    //READING FILE
    while (cin >> first >> second)
    {   // We first used vector<vector<int>> all_edges, but it was too slow.
        // So we used vector<pair<string, string>> all_edges, it is faster to read.
        // Since it is hard to compare string than int, we used dictionary to map string to int.
        if (dictionary.find(first) == dictionary.end()) {  // Add the first node to dictionary if it is not in the dictionary
            dictionary.emplace(first, index);
            index++; // Increase the index of the dictionary items
        }
        if (dictionary.find(second) == dictionary.end()) {  // Add the second node to dictionary if it is not in the dictionary
            dictionary.emplace(second, index);
            index++;
        }
    
        if (out_arrows.find(dictionary[first]) == out_arrows.end()) {  // Add the first node to dictionary if it is not in the dictionary
            out_arrows.emplace(dictionary[first], 1);
        }
        else {
            out_arrows[dictionary[first]]++;
        }
    
        //all_edges.emplace_back(dictionary[first], dictionary[second]);
        if (find(all_edges_new[dictionary[second]].begin(), all_edges_new[dictionary[second]].end(), dictionary[first])==all_edges_new[dictionary[second]].end()) {
            // it is to make the incoming edges unique, i.e it is as if set.
          all_edges_new[dictionary[second]].push_back(dictionary[first]);
        }
        // Since all of them added to dictionary on the first if,
        //  we do not need to check for first again.
        /*
        if (dictionary.find(first) == dictionary.end()) { 
            dictionary.emplace(first, dictionary.size()); 
        }
        */
    }

    fclose(stdin); //Close the file
    //READIN FILE ENDS
    printf("--> File is read\n");
    int n_number_of_nodes = dictionary.size(); // Number of nodes

    printf("--> Out_arrows are calculated\n");
    csr = new CSR(n_number_of_nodes);
    csr->row_begin[0] = 0;
    int priorLen = 0;
    int priorNonzero = 0;

    if (rank == 0) {
        for(int i = 0; i < n_number_of_nodes; i++){
            if(all_edges_new.find(i) == all_edges_new.end()){
                csr->row_begin[i] = 0; 
            }else{
                vector<int>* adjList = &all_edges_new[i];
                csr->row_begin[i] = priorNonzero + priorLen;
                priorNonzero = priorNonzero + priorLen;
                priorLen = adjList->size();
            
                sort(adjList->begin(), adjList->end());
                csr->col_indices.insert(csr->col_indices.end() , adjList->begin(), adjList->end());
            }
        }

        csr->row_begin[n_number_of_nodes] = priorNonzero + priorLen;

        int filledSize = csr->col_indices.size();
        csr->values = new double[filledSize];
        for(int i = 0; i < filledSize; i++){
            if(out_arrows[csr->col_indices[i]]==0){
                printf("Error: out_arrows[csr->col_indices[i]]==0");
                continue;
            }
            csr->values[i] = 1. / out_arrows[csr->col_indices[i]];
        }
   }

    vector<double> r0(n_number_of_nodes, 1.0); // This holds the initial rank of the nodes
    vector<double> r_next(n_number_of_nodes, 0.0); // This holds the next iteration of the pagerank
    int iterations = 0; // This holds the number of iterations
    double alpha = 0.20; // Alpha value
    double diff = 0.0; // Difference between the r0 and r_next
    double epsilon = 0.000001; // Epsilon value
    //Calculate the pagerank
    double x = (1 - alpha) * (1. / n_number_of_nodes); //!!
    double beginTime = omp_get_wtime();
    for (int i = 0; i < n_number_of_nodes; i++)
        r0[i] = 1.0 / n_number_of_nodes;
    
    while(true)
    {
        iterations++;
        diff = 0.0;
        y = 1;
        for (int i = 0; i < n_number_of_nodes; i++)
        {
            r_next[i] = 0.0;
            int y=1;
            for (int j = csr->row_begin[i]; j < csr->row_begin[i + 1]; j++)
            {
                r_next[i] += csr->values[j] * r0[csr->col_indices[j]];
            }
        
            r_next[i] = (x + (alpha)*r_next[i])/n_number_of_nodes;
            diff += abs(r_next[i] - r0[i]);
        }
        r0 = r_next;
        printf("Iteration: %d, Diff: %.*e\n", iterations,10, diff);
        if (diff<epsilon)
            break;
    }
    double endTime = omp_get_wtime();

    printf("--> Pagerank is calculated\n");
    printf("Number of iterations: %d\n", iterations);
    printf("Top 5 nodes: \n");

    vector<pair<double, int>> top5;
    for (int i = 0; i < n_number_of_nodes; i++) 
        top5.emplace_back(r_next[i], i);

    sort(top5.begin(), top5.end(), greater<pair<double, int>>());
    
    for(auto it = dictionary.begin(); it != dictionary.end(); ++it){
        for (int i = 0; i < 5; i++) {
            if(it->second == top5[i].second){
                printf("%.*e ", 10,top5[i].first);
                cout << it->first<<endl;
            }
        }
    }
	return 0;
}
