// Compile with c++ ece650-a2cpp -std=c++11 -o ece650-a2
#include <iostream>
#include <sstream>
#include <vector>
#include <regex>
#include <algorithm>
// defined std::unique_ptr
#include <memory>
// defines Var and Lit
#include "minisat/minisat/core/SolverTypes.h"
// defines Solver
#include "minisat/minisat/core/Solver.h"

using namespace std;
                  
class Graph{
public:
    size_t Vert;
    std::vector<unsigned> Edges;
    std::vector<int> *graph;
    // constructor for class
    Graph(int Vert){
        this -> Vert = Vert;
        graph = new std::vector<int> [Vert];
    }
    //other functions
    void StoreEdges(std::vector<unsigned> edges);
    void Build(); 
    std::vector<int> vertexCover(size_t V_cover, std::vector<unsigned> Edges);
};

void Graph::StoreEdges(std::vector<unsigned> edges){
    Edges = edges;
}


void Graph::Build(){
    //loop through vertices and figure out what vert is connected to what
    for(size_t i = 0; i < Edges.size(); i = i + 1){
        if( (i % 2) != 0){
            graph[Edges[i]].push_back(Edges[i-1]);
        }
        else{
            graph[Edges[i]].push_back(Edges[i+1]);
        }
    }
}
std::vector<int> Graph::vertexCover(size_t V_cover, std::vector<unsigned> Edges){
    // vertex cover code
    // -- allocate on the heap so that we can reset later if needed
    Minisat::Solver solver;

    std::vector<std::vector<Minisat::Lit>> Vertices(Vert);

    // create positive literals in the new vertice vector
    for(size_t i = 0; i < Vert; ++i){
        for(size_t j = 0; j < V_cover; ++j){
            Minisat::Lit literal = Minisat::mkLit(solver.newVar());
            Vertices[i].push_back(literal);
        }
    }

    // Rule 1
    for(size_t i = 0; i < V_cover; ++i){
        Minisat::vec<Minisat::Lit> Literals_n;
        for(size_t j = 0; j < Vert; ++j){
            Literals_n.push(Vertices[j][i]);
        }
        solver.addClause(Literals_n);
        Literals_n.clear();
    }
    // Rule 2
    for(size_t m = 0; m < Vert; ++m){
        for(size_t p = 0; p < (V_cover-1); ++p){
            for(size_t q = p+1; q < V_cover; ++q){
                solver.addClause(~Vertices[m][p],~Vertices[m][q]);
            }
        }
    }
    // Rule 3
    for(size_t m = 0; m < V_cover; ++m){
        for(size_t p = 0; p < (Vert-1); ++p){
            for(size_t q = p+1; q < Vert; ++q){
                solver.addClause(~Vertices[p][m],~Vertices[q][m]);
            }
        }
    }
    // Rule 4
    for(size_t i = 0; i < Edges.size(); i = i + 2){
        Minisat::vec<Minisat::Lit> Literals_n;
        for(size_t k = 0; k < V_cover; ++k){
            Literals_n.push(Vertices[Edges[i]][k]);
            Literals_n.push(Vertices[Edges[i+1]][k]);
        }
        solver.addClause(Literals_n);
        Literals_n.clear();
    }
    auto res = solver.solve();
    std::vector<int> cover;

    if(res){
        for(size_t i = 0; i < Vert; ++i){
            for(size_t j = 0; j < V_cover; ++j){
                int result = Minisat::toInt(solver.modelValue(Vertices[i][j]));
                if(result == 0){
                    cover.push_back(i);
                }
            }
        }
        return cover;
    }
    else{
        return {-1};
    }
}

std::vector<unsigned> ExtractNumber(std::string line,std::vector<unsigned> nums){
    int i = 1;
    std::smatch ematch;
    regex E("([0-9]+)");
    while (regex_search(line, ematch, E)) {   
        nums.push_back(std::stoi(ematch.str(0)));
        i++; 
        // suffix to find the rest of the string. 
        line = ematch.suffix().str();    
    }  
    return(nums); 
}


int main() {

    // Variables
    std::string line;
    std::smatch ematch;
    Graph street(0);

    // read from stdin until EOF
    while (!std::cin.eof()) {
        std::vector<unsigned> nums;
        std::vector<int> SP;
        // read a line of input until EOL and store in a string 
        std::getline(std::cin, line);

        // Use regex to detect command characters and numbers
        regex C("^\\w");
        regex E("([0-9]+)");
        regex V("[0-9]+");
        std::smatch cmatch;

        if (regex_search(line,cmatch,C)){

            // Depending on command character, do different regex extractions
            if (cmatch[0] == 'E'){
                nums = ExtractNumber(line,nums);
                if(nums.size()){
                    size_t Vertcomp = street.Vert;
                    if ((street.Vert > 1) && (*std::max_element(std::begin(nums) ,std::end(nums))<Vertcomp)){
                        street.StoreEdges(nums);
                        // Build connected vertice graph (adjacent nodes)
                        street.Build();
                        std::vector<int> Cover;
                        size_t i = 1;
                        while(i < street.Vert){
                            Cover = street.vertexCover(i,street.Edges);
                            if(Cover[0] == -1){
                                i = i + 1;
                                continue;
                            }
                            else{
                                break;
                            }
                        }
                        sort(Cover.begin(),Cover.end());
                        for(size_t i = 0; i < Cover.size(); ++i){
                            std::cout << Cover[i] << " ";
                        }
                        std::cout << std::endl;
                    }
                    else if(street.Vert == 1){
                        continue;
                    }
                    else{
                        std::cerr << "Error: Not enough vertices to compute edges" << std::endl;
                    }
                }
                else{
                    std::vector<unsigned> zero;
                    street.StoreEdges(zero);
                    std::cout << std::endl;
                }
            }
            else if (cmatch[0] == 'V'){
                nums = ExtractNumber(line,nums);
                if(nums[0]>0){
                    street = Graph(nums[0]);
                }
                else{
                    std::cerr << "Error: Need to enter at least 1 Vertice" << std::endl;
                }
            }
            else{
                std::cerr << "Error: No command match" << std::endl;
            }
        }
        else{
            continue;
        }
    }
    return 0;
}