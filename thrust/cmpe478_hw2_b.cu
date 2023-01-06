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
#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/inner_product.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/reduce.h>
#include <thrust/transform.h>

template<typename T>
struct abs_ : public unary_function<T,T>
{
  __host__ __device__ T operator()(const T &x,const T &x) const
  {
    return abs(x - y);
  }
};

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
    #if DEBUG_FILE_READ
     printf("--> Dictionary: \n");
    {    
        for (auto it = dictionary.begin(); it != dictionary.end(); ++it)
    {
        cout << it->first << " => " << it->second << '\n';
    }
        printf("---------------\n");
    }
    #endif
    
    #if DEBUG_FILE_READ
    {
        printf("--> N_number_of_nodes: %d\n", n_number_of_nodes);
        printf("--> All_edges: \n");

            for (auto k : all_edges_new){
                sort(k.second.begin(), k.second.end());
                cout << k.first<< " => ";
                for(int j = 0; j < k.second.size(); j++){
                    cout << k.second[j] << " ";
                }
             cout << endl;
        }
        printf("---------------\n");
    }
    #endif

    printf("--> Out_arrows are calculated\n");
    csr = new CSR(n_number_of_nodes);
    csr->row_begin[0] = 0;
    int priorLen = 0;
    int priorNonzero = 0;
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

    #if CSR_PRINT
    printf("row_begins: \n");
    for(int i = 0; i < n_number_of_nodes+1; i++){
        printf("%d ",csr->row_begin[i]);
    }
    printf("\n");

    printf("col_indices: \n");
    for(int i = 0; i < csr->col_indices.size(); i++){
        printf("%d ",csr->col_indices[i]);
    }
    printf("\n");

    printf("values: \n");
    for(int i = 0; i < filledSize; i++){
        printf("%f ",csr->values[i]);
    }
    printf("\n");
    #endif
    
    //Transfer the CSR to GPU
    thrust::device_vector<double> values=csr->values; // This holds the initial rank of the nodes
    thrust::device_vector<int> col_indices=csr->col_indices; // This holds the initial rank of the nodes
    thrust::device_vector<int> row_begin=csr->row_begin; // This holds the initial rank of the nodes

    // We need to r0 and r_next to GPU.
    //vector<double> r0_(n_number_of_nodes, 1.0); // This holds the initial rank of the nodes
    thrust::device_vector<double> r0(n_number_of_nodes,1/n_number_of_nodes); // This holds the initial rank of the nodes
    //vector<double> r_next(n_number_of_nodes, 0.0); // This holds the next iteration of the pagerank
    thrust::device_vector<double> r_next(n_number_of_nodes, 0.0); // This holds the next iteration of the pagerank
    
    int iterations = 0; // This holds the number of iterations
    double alpha = 0.20; // Alpha value
    //double diff = 0.0; // Difference between the r0 and r_next
    thrust::device_vector<double> diff(n_number_of_nodes, 0.0); // Difference between the r0 and r_next
    double epsilon = 0.000001; // Epsilon value
    //Calculate the pagerank
    double x = (1 - alpha) * (1. / n_number_of_nodes); //!!

    double bt = omp_get_wtime();

    
    while(true)
    {
        iterations++;
        for (int i = 0; i < n_number_of_nodes; i++)
        {
            //make_permutation_iterator(
            //    InputIterator first,
            //    Iterator map
            //)
            // Returns an iterator that dereferences to first[map[i]] for each i.
            auto p_iter = thrust::make_permutation_iterator(r_next.begin(), col_indices.begin()+row_begin[i]);
            //inner_product(
            //    InputIterator1 first1,
            //    InputIterator1 last1,
            //    InputIterator2 first2,
            //    T init)
            // Returns the inner product of the two ranges . The inner product is defined as:
            // init + (first1[0] * first2[0]) + (first1[1] * first2[1]) + ... + (first1[n-1] * first2[n-1])
            //r_next[i] += csr->values[j] * r0[csr->col_indices[j]];
            auto i_prod = thrust::inner_product(row_begin[i] + values.begin(), values.begin() + rowBegin[i + 1], p_iter, 0.0);
            r_next[i] = ((i_prod*alpha)+x)/n_number_of_nodes;
        }

       
        // transform_reduce()
        // first	The beginning of the sequence.
        // last	The end of the sequence.
        // unary_op	The function to apply to each element of the input sequence.
        // init	The result is initialized to this value.
        // binary_op	The reduction operation.
         //diff += abs(r_next[i] - r0[i]);
        double result = sthrust::transform_reduce(r_next.begin(), r_next.end(), r0.begin(), diff.begin(), abs_<double>(), 0.0, thrust::plus<double>());
        if (result <= epsilon) 
            break;
         else {
            //r0 = r_next; device to variable
            thrust::copy(r_next.begin(), r_next.end(), r0.begin());
        }
        printf("Iteration: %d, Diff: %.*e\n", iterations,10, diff);
    }

    double et = omp_get_wtime();
    cout<<"Elapsed duration: " <<et-bt<<endl;


    
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