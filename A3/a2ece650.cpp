// Compile with c++ ece650-a2cpp -std=c++11 -o ece650-a2
#include <iostream>
#include <sstream>
#include <vector>
#include <regex>
#include <algorithm>

using namespace std;

// Algorithm references: 1. https://en.wikipedia.org/wiki/Breadth-first_search
//                       2. https://www.geeksforgeeks.org/breadth-first-search-or-bfs-for-a-graph/
class Graph{
public:
    int Vert;
    std::vector<unsigned> Edges;
    std::vector<int> *graph;
    // constructor for class
    Graph(int Vert){
        this -> Vert = Vert;
        graph = new std::vector<int> [Vert];
    }
    //other functions
    void StoreEdges(std::vector<unsigned> edges);
    std::vector<int> Shortestpath(int a, int b);
    void Build(); 
};

void Graph::StoreEdges(std::vector<unsigned> edges){
    Edges = edges;
}


std::vector<int> Graph::Shortestpath(int a, int b){   

    // If the entered values aren't just equal, find the path using modified BFS!
    if(a != b){
        // Use implementation of breadth first search to find shortest path between a and b (pseudo code from Wikipedia)
        // Set all vertices to be unvisited or "false"
        bool *Vert_visited = new bool[Vert];

        for(int i = 0; i < Vert; i++)
            Vert_visited[i] = false;
        // use queue implementation to hold data sequences
        std::vector<std::vector<int>> Q;
        std::vector<int> Current_Node;
        std::vector<int> Nearby;
        // Mark the starting point "a" as true or visited
        Current_Node.push_back(a);
        Q.push_back(Current_Node); 
        int Node;

        // While the queue vector is not empty
        while(!Q.empty())
        {
            // dequeue the vertex and examine current Node
            Current_Node = Q.front();
            Q.erase(Q.begin()); 
            Node = Current_Node.back();
            // retrieve adjacent vertices (from build graph) and examine each adjacent one
            if (!Vert_visited[Node])
            {
                // Pull vector of adjacent vertices at this node
                Nearby = graph[Node];
                // Loop through nodes adjacent vertices
                int vertex;
                for (size_t i = 0; i < Nearby.size(); i = i + 1)
                {
                    vertex = Nearby[i];
                    // Define a placeholder for the current Node/vertex (so you don't add extra vertices)
                    std::vector<int> New_Node;
                    // enqueue the node
                    New_Node = Current_Node;
                    New_Node.push_back(vertex);
                    // Check if we still need to travel to find destination Vertex
                    if(vertex != b){
                        Q.push_back(New_Node);
                    }
                    else{
                        return New_Node;
                    }
                }
                Vert_visited[Node] = true;
            } 
        }
        // if nothing found return empty vector = false
        std::vector<int> None;
        return None;
    }
    // If equal return 
    else{
        std::vector<int> Same;
        return Same;
    }
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

bool SameCheck(std::vector<unsigned> edges){
    if(edges.size() > 2){
        for(size_t i = 0; i < edges.size(); i = i + 2){
            if(edges[i] == edges[i+1]){
                return true;
            }
            else{
                continue;
            }
        }
    }
    else{
        return false;
    }
    return false;    
}

// checks edges to find duplicate
bool DupeCheck(std::vector<unsigned> edges){
    // initialize boolian flags
    bool Flag_1 = false;
    bool Flag_2 = false;
    // double forloop to check every vertex pair with each other; throw flags if match found
    if(edges.size() > 2){
        for(size_t i = 0; i < edges.size(); i = i + 2){
            for(size_t j = i + 2; j < edges.size(); j = j + 2){
                if(edges[i] == edges[j] || edges[i] == edges[j+1]){
                    Flag_1 = true;
                }
                if(edges[i+1] == edges[j] || edges[i+1] == edges[j+1]){
                    Flag_2 = true;
                }
                // check flag status
                if(Flag_1 && Flag_2){
                    return true;
                }
                else{
                    Flag_1 = false;
                    Flag_2 = false;
                }
            }
        }
    }
    else{
        return false;
    }
    return false;
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
        regex s("\\s[0-9]+\\s[0-9]+");
        std::smatch cmatch;

        if (regex_search(line,cmatch,C)){

            // Depending on command character, do different regex extractions
            if (cmatch[0] == 'E'){
                std::cout << line << std::endl;
                nums = ExtractNumber(line,nums);
                if(nums.size()){
                    size_t Vertcomp = street.Vert;
                    if ((street.Vert > 1) && (*std::max_element(std::begin(nums) ,std::end(nums))<Vertcomp)){
                        if(!DupeCheck(nums) && !SameCheck(nums)){  
                            street.StoreEdges(nums);
                            // Build connected vertice graph (adjacent nodes)
                            street.Build();
                            } 
                        else{
                            std::cerr << "Error: Duplicate edge entry" << std::endl;
                            return(1);
                            }   
                    }
                    else if(street.Vert == 1){
                        continue;
                    }
                    else{
                        std::cerr << "Error: Not enough vertices to compute edges" << std::endl;
                        return(1);
                    }
                }
                else{
                    std::vector<unsigned> zero;
                    street.StoreEdges(zero);
                }
            }
            else if (cmatch[0] == 'V'){
                std::cout << line << std::endl;
                nums = ExtractNumber(line,nums);
                if(nums[0]>0){
                    street = Graph(nums[0]);
                }
                else{
                    std::cerr << "Error: Need to enter at least 1 Vertice" << std::endl;
                    return(1);
                }
            }
            else if (cmatch[0] == 's'){
                nums = ExtractNumber(line,nums);
                if (street.Vert > 1 && street.Edges.size()>1){
                    size_t Vertcomp = street.Vert;
                    if((*std::max_element(std::begin(nums) ,std::end(nums))<Vertcomp)){
                        SP = street.Shortestpath(nums[0],nums[1]);
                        if(!SP.empty()){
                            std::string OUT;
                            int X;
                            for (size_t i = 0; i < SP.size(); i = i + 1)
                            {
                                X = SP[i];
                                if(!OUT.empty()){
                                    OUT = OUT + "-"; 
                                }
                                OUT = OUT + to_string(X);
                            }   
                            std::cout << OUT << std::endl;
                        }
                        else{
                            std::cerr << "Error: No path between vertices" << std::endl;
                            return(1);
                        } 
                    }
                    else{
                        std::cerr << "Error: Requested a path to a vertex that does not exist" << std::endl;
                        return(1);
                    }    
                }
                else{
                    std::cerr << "Error: Not enough vertices to compute edges" << std::endl;
                    return(1);
                }
            }
            else{
                std::cerr << "Error: No command match" << std::endl;
                return(1);
            }
        }
        else{
            continue;
        }
    }
    return 0;
}