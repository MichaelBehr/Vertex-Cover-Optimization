// Compile with c++ ece650-a2cpp -std=c++11 -o ece650-a2
#include <algorithm>
#include <iostream>
#include <iterator>
#include <pthread.h>
#include <regex>
#include <sstream>
#include <vector>
#include <time.h>
// defined std::unique_ptr
#include <memory>
// defines Var and Lit
#include "minisat/minisat/core/SolverTypes.h"
// defines Solver
#include "minisat/minisat/core/Solver.h"
#define NUM_THREADS 3


// Taken from Linux example page on pthread_getcpuclockid //

#define handle_error(msg) \
       do { perror(msg); exit(EXIT_FAILURE); } while (0)

#define handle_error_en(en, msg) \
       do { errno = en; perror(msg); exit(EXIT_FAILURE); } while (0)

timespec pclock(clockid_t cid){
  struct timespec ts;
  if (clock_gettime(cid, &ts) == -1){
    handle_error("clock_gettime");
  }
  //std::cout << "Microseconds taken for thread: " << (ts.tv_sec * 1000000 + ts.tv_nsec / 1000) << std::endl;
  //printf("%4ld.%03ld\n", ts.tv_sec, ts.tv_nsec / 1000000);
  return ts;
}


class Graph {
public:
  size_t Vert;
  std::vector<unsigned> Edges;
  std::vector<int> *graph;

  std::vector<int> cnf_sat_cover;
  size_t V_cover;
  std::vector<int> approx_1_cover;
  std::vector<int> approx_2_cover;

  // Clockid variables
  clockid_t start_CNF;
  clockid_t start_AP_1;
  clockid_t start_AP_2;
  clockid_t end_CNF;
  clockid_t end_AP_1;
  clockid_t end_AP_2;

  // Timing information variables

  // constructor for class
  Graph(int Vert) {
    this->Vert = Vert;
    graph = new std::vector<int>[Vert];
  }
  // other functions
  void StoreEdges(std::vector<unsigned> edges);
  void Build();
  std::vector<int> Highest_Degree_Vert();

  // Vertex cover algorithms
  void CNF_SAT_VC(void);
  void APPROX_VC_1(void);
  void APPROX_VC_2(void);
};

void Graph::StoreEdges(std::vector<unsigned> edges) { Edges = edges; }

void Graph::Build() {
  // loop through vertices and figure out what vert is connected to what
  for (size_t i = 0; i < Edges.size(); i = i + 1) {
    if ((i % 2) != 0) {
      graph[Edges[i]].push_back(Edges[i - 1]);
    } else {
      graph[Edges[i]].push_back(Edges[i + 1]);
    }
  }
}

std::vector<int> Graph::Highest_Degree_Vert() {
  // Loop through the edge list and create vertex incident size vector
  std::vector<int> incident_size_vector;
  for (size_t i = 0; i < Vert; i = i + 1) {
    // std::cout << "Vertex " << i << " has" << graph[i].size() << "incident
    // vertices" << std::endl;
    incident_size_vector.push_back(graph[i].size());
  }
  return (incident_size_vector);
}

void Graph::CNF_SAT_VC() {
  // vertex cover code
  // -- allocate on the heap so that we can reset later if needed
  std::vector<unsigned> Edges = this->Edges;
  size_t V_cover = this->V_cover;
  Minisat::Solver solver;

  std::vector<std::vector<Minisat::Lit>> Vertices(Vert);

  // create positive literals in the new vertice vector
  for (size_t i = 0; i < Vert; ++i) {
    for (size_t j = 0; j < V_cover; ++j) {
      Minisat::Lit literal = Minisat::mkLit(solver.newVar());
      Vertices[i].push_back(literal);
    }
  }

  // Rule 1
  for (size_t i = 0; i < V_cover; ++i) {
    Minisat::vec<Minisat::Lit> Literals_n;
    for (size_t j = 0; j < Vert; ++j) {
      Literals_n.push(Vertices[j][i]);
    }
    solver.addClause(Literals_n);
    Literals_n.clear();
  }
  // Rule 2
  for (size_t m = 0; m < Vert; ++m) {
    for (size_t p = 0; p < (V_cover - 1); ++p) {
      for (size_t q = p + 1; q < V_cover; ++q) {
        solver.addClause(~Vertices[m][p], ~Vertices[m][q]);
      }
    }
  }
  // Rule 3
  for (size_t m = 0; m < V_cover; ++m) {
    for (size_t p = 0; p < (Vert - 1); ++p) {
      for (size_t q = p + 1; q < Vert; ++q) {
        solver.addClause(~Vertices[p][m], ~Vertices[q][m]);
      }
    }
  }
  // Rule 4
  for (size_t i = 0; i < Edges.size(); i = i + 2) {
    Minisat::vec<Minisat::Lit> Literals_n;
    for (size_t k = 0; k < V_cover; ++k) {
      Literals_n.push(Vertices[Edges[i]][k]);
      Literals_n.push(Vertices[Edges[i + 1]][k]);
    }
    solver.addClause(Literals_n);
    Literals_n.clear();
  }
  auto res = solver.solve();
  std::vector<int> cover;

  if (!res) {
    this->cnf_sat_cover = {-1};
  } else {
    for (size_t i = 0; i < Vert; ++i) {
      for (size_t j = 0; j < V_cover; ++j) {
        int result = Minisat::toInt(solver.modelValue(Vertices[i][j]));
        if (result == 0) {
          cover.push_back(i);
        }
      }
    }
    this->cnf_sat_cover = cover;
  }
}

void Graph::APPROX_VC_1() {
  // NOTE: nums is a vector representing the edge list
  std::vector<unsigned> nums = Edges;
  // Initialize vertex cover vector
  std::vector<int> Vertex_Cover;

  // create static adjacency size vector
  std::vector<int> incident_size_vector = Highest_Degree_Vert();

  // create chancing adjacency size vector
  std::vector<int> changing_vector = Highest_Degree_Vert();

  size_t idx;
  int max;
  // While this code executes, we will remove edges until all gone -> nums ==
  // empty
  while (!nums.empty()) {

    // determine the location and size of the vertex with most incident vertices
    std::vector<int>::iterator result = std::max_element(
        std::begin(changing_vector), std::end(changing_vector));
    idx = std::distance(std::begin(changing_vector), result);
    max = incident_size_vector[idx];

    Vertex_Cover.push_back(idx);
    // std::cout<< "Adding " << idx << " to Vertex Cover" << std::endl;
    // Now add the edges belonging to this index to our vertex cover list

    // std::cout << "Starting approx loop" << std::endl;
    for (int i = 0; i < max; i = i + 1) {

      size_t Vertex = graph[idx][i];

      // Decrement incident size at both connected vertex locations
      changing_vector[Vertex] = changing_vector[Vertex] - 1;
      changing_vector[idx] = changing_vector[idx] - 1;

      // Now create a temporary vector to hold the new edge list after removing
      // edge
      std::vector<unsigned> temp;
      for (size_t j = 0; (j + 1) < (nums.size()); j = j + 2) {
        // std::cout << "j: " << j << std::endl;
        // std::cout << "Checking edge [" << nums[j] << "," << nums[j+1] << "]"
        // << "against [" << idx << "," << Vertex << "]" << std::endl;

        if (((nums[j] == Vertex) && (nums[j + 1] == idx)) ||
            ((nums[j] == idx) && (nums[j + 1] == Vertex))) {
          // std::cout<< "Found edge [" << nums[j] << "," << nums[j+1] << "]" <<
          // std::endl;
          continue;
        } else {
          temp.push_back(nums[j]);
          temp.push_back(nums[j + 1]);
        }
      }
      nums = temp;
      temp.clear();
    }
  }
  this->approx_1_cover = Vertex_Cover;
}

void Graph::APPROX_VC_2() {
  // NOTE: nums is a vector representing the edge list
  std::vector<unsigned> nums = Edges;
  std::vector<int> incident_size_vector = Highest_Degree_Vert();
  // Initialize vertex cover vector
  std::vector<int> Vertex_Cover;

  size_t idx1;
  size_t idx2;
  // While this code executes, we will remove edges until all gone -> nums ==
  // empty
  while (!nums.empty()) {

    // determine the location and size of the vertex with most incident vertices
    idx1 = nums[0];
    idx2 = nums[1];

    Vertex_Cover.push_back(idx1);
    Vertex_Cover.push_back(idx2);

    // Now create a temporary vector to hold the new edge list after removing
    // edge
    std::vector<unsigned> temp;
    for (size_t j = 0; (j + 1) < (nums.size()); j = j + 2) {
      // std::cout << "Checking edge [" << nums[j] << "," << nums[j+1] << "]" <<
      // "against [" << idx << "," << Vertex << "]" << std::endl;

      if ((nums[j] == idx1) || (nums[j + 1] == idx2) || (nums[j] == idx2) ||
          (nums[j + 1] == idx1)) {
        // std::cout<< "Found edge [" << nums[j] << "," << nums[j+1] << "]" <<
        // std::endl;
        continue;
      } else {
        temp.push_back(nums[j]);
        temp.push_back(nums[j + 1]);
      }
    }
    nums = temp;
    temp.clear();
  }
  this->approx_2_cover = Vertex_Cover;
}

std::vector<unsigned> ExtractNumber(std::string line,
                                    std::vector<unsigned> nums) {
  int i = 1;
  std::smatch ematch;
  std::regex E("([0-9]+)");
  while (regex_search(line, ematch, E)) {
    nums.push_back(std::stoi(ematch.str(0)));
    i++;
    // suffix to find the rest of the string.
    line = ematch.suffix().str();
  }
  return (nums);
}

void *vnc_sat_function(void *initialized_graph) {
  // Thread function
  Graph *graph_pointer = (Graph *)initialized_graph;
  size_t i = 1;
  int s;
  std::vector<int> Cover_CNFSAT;
  while (i < graph_pointer->Vert) {
    graph_pointer->V_cover = i;
    graph_pointer->CNF_SAT_VC();
    Cover_CNFSAT = graph_pointer->cnf_sat_cover;
    if (Cover_CNFSAT[0] == -1) {
      i = i + 1;
      continue;
    } else {
      break;
    }
  }
  graph_pointer->CNF_SAT_VC();

  // Clock timing using getcpuclockid just like in the linux example page
  int state = pthread_getcpuclockid(pthread_self(), &(graph_pointer->start_CNF));
  if (state != 0){
    handle_error_en(s, "pthread_getcpuclockid");
  }

  timespec CNF_tspec =  pclock(graph_pointer->start_CNF);


  pthread_exit(NULL);
}

void *approx_1_function(void *initialized_graph) {
  Graph *graph_pointer = (Graph *)initialized_graph;
  graph_pointer->APPROX_VC_1();
  int s;

  // Clock timing using getcpuclockid just like in the linux example page
  int state = pthread_getcpuclockid(pthread_self(), &(graph_pointer->start_AP_1));
  if (state != 0){
    handle_error_en(s, "pthread_getcpuclockid");
  }

  timespec AP_1_tspec =  pclock(graph_pointer->start_AP_1);


  pthread_exit(NULL);
}

void *approx_2_function(void *initialized_graph) {
  Graph *graph_pointer = (Graph *)initialized_graph;
  graph_pointer->APPROX_VC_2();
  int s;
  // Clock timing using getcpuclockid just like in the linux example page
  int state = pthread_getcpuclockid(pthread_self(), &(graph_pointer->start_AP_2));
  if (state != 0){
    handle_error_en(s, "pthread_getcpuclockid");
  }

  timespec AP_2_tspec =  pclock(graph_pointer->start_AP_2);

  pthread_exit(NULL);
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
    std::regex C("^\\w");
    std::regex E("([0-9]+)");
    std::regex V("[0-9]+");
    std::smatch cmatch;

    if (regex_search(line, cmatch, C)) {

      // Depending on command character, do different regex extractions
      if (cmatch[0] == 'E') {
        nums = ExtractNumber(line, nums);
        if (nums.size()) {
          size_t Vertcomp = street.Vert;
          if ((street.Vert > 1) &&
              (*std::max_element(std::begin(nums), std::end(nums)) <
               Vertcomp)) {
            street.StoreEdges(nums);
            // Build connected vertice graph (adjacent nodes)
            street.Build();
            std::vector<int> Cover_CNFSAT;
            std::vector<int> Cover_VC1;
            std::vector<int> Cover_VC2;

            // Create threads calling vertex cover functions.
            pthread_t cnf_thread;
            pthread_t approx_1_thread;
            pthread_t approx_2_thread;

            clockid_t CNF_cid;
            clockid_t APPROX_1_cid;
            clockid_t APPROX_2_cid;

            int thread_action_return[NUM_THREADS];
            thread_action_return[0] = pthread_create(
                &cnf_thread, NULL, vnc_sat_function, (void *)&street);
            thread_action_return[1] = pthread_create(
                &approx_1_thread, NULL, approx_1_function, (void *)&street);
            thread_action_return[2] = pthread_create(
                &approx_2_thread, NULL, approx_2_function, (void *)&street);

            for (int i = 0; i < NUM_THREADS; i++) {
              if (thread_action_return[i]) {
                std::cerr << "Error: Unable to create thread: " << i
                          << std::endl;
                exit(-1);
              }
            }

            // Join all threads, no zombie threads.
            void *status;
            bool flag = false;
            //std::cout << "System clock time BEFORE addition: " << ts_CNF.tv_sec << std::endl;
            struct timespec ts_CNF;
            if (clock_gettime(CLOCK_REALTIME, &ts_CNF) == -1) {
              std::cerr << "Clock error" << std::endl;
            }
            ts_CNF.tv_sec += 120;
            //std::cout << "System clock time AFTER addition: " << ts_CNF.tv_sec << std::endl;
            thread_action_return[0] =
                pthread_timedjoin_np(cnf_thread, NULL, &ts_CNF);
            if (thread_action_return[0]) {
              flag = true;
            }
            thread_action_return[1] = pthread_join(approx_1_thread, &status);
            if (thread_action_return[1]) {
              std::cerr << "Error: Unable to join," << thread_action_return[1]
                        << std::endl;
              exit(-1);
            }
            thread_action_return[2] = pthread_join(approx_2_thread, &status);
            if (thread_action_return[2]) {
              std::cerr << "Error: Unable to join," << thread_action_return[2]
                        << std::endl;
              exit(-1);
            }

            // All threads executed and ended successfully here.

            // Get covers.
            Cover_CNFSAT = street.cnf_sat_cover;
            Cover_VC1 = street.approx_1_cover;
            Cover_VC2 = street.approx_2_cover;
            // SORT all 3 covers
            sort(Cover_CNFSAT.begin(), Cover_CNFSAT.end());
            sort(Cover_VC1.begin(), Cover_VC1.end());
            sort(Cover_VC2.begin(), Cover_VC2.end());

            // Print CNF-SAT-VC
            std::cout << "CNF-SAT-VC: ";
            if(flag==false){
              for (size_t i = 0; i < (Cover_CNFSAT.size()); i = i + 1) {
                if (i == (Cover_CNFSAT.size() - 1)) {
                  std::cout << Cover_CNFSAT[i];
                } else {
                  std::cout << Cover_CNFSAT[i] << ",";
                }
              }
            }
            else{
              std::cerr << "timeout";
            }
            std::cout << std::endl;

            // Print APPROX-VC-1
            std::cout << "APPROX-VC-1: ";
            for (size_t i = 0; i < (Cover_VC1.size()); i = i + 1) {
              if (i == (Cover_VC1.size() - 1)) {
                std::cout << Cover_VC1[i];
              } else {
                std::cout << Cover_VC1[i] << ",";
              }
            }
            std::cout << std::endl;

            // Print APPROX-VC-2
            std::cout << "APPROX-VC-2: ";
            for (size_t i = 0; i < (Cover_VC2.size()); i = i + 1) {
              if (i == (Cover_VC2.size() - 1)) {
                std::cout << Cover_VC2[i];
              } else {
                std::cout << Cover_VC2[i] << ",";
              }
            }
            std::cout << std::endl;
          } else if (street.Vert == 1) {
            continue;
          } else {
            std::cerr << "Error: Not enough vertices to compute edges"
                      << std::endl;
          }
        } else {
          std::vector<unsigned> zero;
          street.StoreEdges(zero);
          std::cout << "CNF-SAT-VC: " << std::endl;
          std::cout << "APPROX-VC-1: " << std::endl;
          std::cout << "APPROX-VC-2: " << std::endl;
        }
      } else if (cmatch[0] == 'V') {
        nums = ExtractNumber(line, nums);
        if (nums[0] > 0) {
          street = Graph(nums[0]);
        } else {
          std::cerr << "Error: Need to enter at least 1 Vertice" << std::endl;
        }
      } else {
        std::cerr << "Error: No command match" << std::endl;
      }
    } else {
      continue;
    }
  }
  return 0;
}
