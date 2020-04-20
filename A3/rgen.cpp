// Compile with c++ ece650-a2cpp -std=c++11 -o ece650-a2
#include <iostream>
#include <sstream>
#include <vector>
#include <regex>
#include <algorithm>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
using namespace std;

unsigned int Random_init(){
    unsigned int randomnumber = 0;
    std::ifstream urandom("/dev/urandom", std::ios::in | std::ios::binary);

    if (urandom.fail()) {
        std::cerr << "Error: cannot open /dev/urandom\n";
        return 1;
    }
    urandom.read((char *) &randomnumber, 1);
    urandom.close();
    return randomnumber;
}
int Rand_Num(unsigned int rand_init, int a,int k){
    int randomnumber = round((((float)rand_init / 255)) * (float)(k-a)) + a;
    return randomnumber;
}

bool Onsegment(std::vector<int> p3,std::vector<int> p4, std::vector<int> p2){
    //check if point is within bounds of line segment. If yes return true
    if ((p2[0] <= std::max(p3[0],p4[0])) && (p2[0] >= std::min(p3[0],p4[0])) && (p2[1] <= std::max(p3[1],p4[1])) && (p2[1] >= std::min(p3[1],p4[1]))){
        return true;
    }
    return false;
}

bool overlap(std::vector<int> p1,std::vector<int> p2, std::vector<int> p3, std::vector<int> p4){
    // first determine if points are in the bounding box of one another
    if((std::max(p1[0],p2[0]) < std::min(p3[0],p4[0])) || (std::min(p1[0],p2[0]) > std::max(p3[0],p4[0])) || (std::max(p1[1],p2[1]) < std::min(p3[1],p4[1])) || (std::min(p1[1],p2[1]) > std::max(p3[1],p4[1]))){
        return false;
    }
    else{
        //use cramers rule and logic from assignment 1
        int a1 = p1[1] - p2[1];
        int b1 = p2[0] - p1[0];
        int a2 = p3[1] - p4[1];
        int b2 = p4[0] - p3[0];

        int Det = a1*b2 - b1*a2;

        int Testp1 = (p3[0]-p1[0])*(p2[1]-p4[1]) - (p2[0]-p4[0])*(p3[1]-p1[1]);
        int Testp2 = (p4[0]-p1[0])*(p3[1]-p2[1]) - (p3[0]-p2[0])*(p4[1]-p1[1]);
        if(Det == 0 && Testp1 == 0 && Testp2 == 0){
            int Testp3 = (p2[0]-p1[0])*(p3[1]-p4[1]) - (p3[0]-p4[0])*(p2[1]-p1[1]);
            int Testp4 = (p4[0]-p1[0])*(p2[1]-p3[1]) - (p2[0]-p3[0])*(p4[1]-p1[1]);
            if (Testp3 == 0 && Testp4 == 0){
                if(((((p2[0] == p3[0]) && (p2[1] == p3[1])) && (!Onsegment(p1,p2,p4))) || ((!Onsegment(p1,p2,p3)) && ((p2[0] == p4[0]) && (p2[1] == p4[1]))))){
                    return false;
                }
                else{
                    return true;
                }
            }
            else{
                return false;
            }
        }
        if(Det != 0){
            return false;
        }
        else{
            return false;
        }

    }
    return false;
}

bool checkoverlap(int x2,int y2,std::vector<std::vector<int>> streets,std::vector<int> segments,int index){
    if(index > 0){
        // compare point x2,y2 against all other points and return true if matched
        std::vector<int> p2;
        p2.push_back(x2);
        p2.push_back(y2);
        int x1 = segments[segments.size()-2];
        int y1 = segments[segments.size()-1];
        std::vector<int> p1;
        p1.push_back(x1);
        p1.push_back(y1);
        std::vector<int> p4;
        std::vector<int> p3;
        std::vector<int> oldstreet;
        //define other point variables
        int x3;
        int y3;
        int x4;
        int y4;

        for(int s = 0; s < index; ++s){
            oldstreet.clear();
            oldstreet = streets[s];

            // now look for overlap of segment with the other street segments
            for(unsigned int i = 0; i < (oldstreet.size()-2);i = i + 2){
                x3 = oldstreet[i];
                y3 = oldstreet[i+1];
                x4 = oldstreet[i+2];
                y4 = oldstreet[i+3];
                p3.clear();
                p4.clear();
                p3.push_back(x3);
                p3.push_back(y3);
                p4.push_back(x4);
                p4.push_back(y4);
                if(overlap(p1,p2,p3,p4)){
                    return true;
                }
                else{
                    continue;
                }
            }
        }
    }
    else{
        return false;
    }
    return false;
}

//checks if segments intersect
bool intersect(std::vector<int> p1,std::vector<int> p2, std::vector<int> p3, std::vector<int> p4){
    // first determine if points are in the bounding box of one another
    if((std::max(p1[0],p2[0]) < std::min(p3[0],p4[0])) || (std::min(p1[0],p2[0]) > std::max(p3[0],p4[0])) || (std::max(p1[1],p2[1]) < std::min(p3[1],p4[1])) || (std::min(p1[1],p2[1]) > std::max(p3[1],p4[1]))){
        return false;
    }
    else{
        //use cramers rule and logic from assignment 1
        int a1 = p1[1] - p2[1];
        int b1 = p2[0] - p1[0];
        int c1 = p2[0]*p1[1] - p1[0]*p2[1];
        int a2 = p3[1] - p4[1];
        int b2 = p4[0] - p3[0];
        int c2 = p4[0]*p3[1] - p3[0]*p4[1];

        int Det = a1*b2 - b1*a2;
        int Detx = c1*b2 - b1*c2;
        int Dety = a1*c2 - c1*a2;

        int Testp1 = (p3[0]-p1[0])*(p2[1]-p4[1]) - (p2[0]-p4[0])*(p3[1]-p1[1]);
        int Testp2 = (p4[0]-p1[0])*(p3[1]-p2[1]) - (p3[0]-p2[0])*(p4[1]-p1[1]);
        if(Det == 0 && Testp1 == 0 && Testp2 == 0){
            int Testp3 = (p2[0]-p1[0])*(p3[1]-p4[1]) - (p3[0]-p4[0])*(p2[1]-p1[1]);
            int Testp4 = (p4[0]-p1[0])*(p2[1]-p3[1]) - (p2[0]-p3[0])*(p4[1]-p1[1]);
            if (Testp3 == 0 && Testp4 == 0){
                if((p4[0] == p1[0] && p4[1] == p1[1]) && !Onsegment(p3,p4,p2)){
                    return false;
                }
                else{
                    return true;
                }
            }
            else{
                return false;
            }
        }
        if(Det != 0){
            float X = (float)Detx / (float)Det;
            float Y = (float)Dety / (float)Det;
            if((X < (std::max(std::min(p1[0],p2[0]),std::min(p3[0],p4[0])))) || (X > (std::min(std::max(p1[0],p2[0]),std::max(p3[0],p4[0])))) || (Y < (std::max(std::min(p1[1],p2[1]),std::min(p3[1],p4[1])))) || (Y > (std::min(std::max(p1[1],p2[1]),std::max(p3[1],p4[1]))))){
                return false;
            }
            else{
                if((p4[0] == p1[0] && p4[1] == p1[1])){
                    return false;
                }
                else{
                    return true;
                }
            }    
        }
        else{
            return false;
        }

    }
    return false;
}

bool checkintersect(int x2,int y2,std::vector<int> segments){
    // compare point x2,y2 against all other points and return true if matched
    std::vector<int> p2;
    p2.push_back(x2);
    p2.push_back(y2);
    for(unsigned int i = 0; i < segments.size();i = i + 2){
        if(x2 == segments[i] && y2 == segments[i+1]){
            return true;
        }
    }
    int x1 = segments[segments.size()-2];
    int y1 = segments[segments.size()-1];
    std::vector<int> p1;
    p1.push_back(x1);
    p1.push_back(y1);
    std::vector<int> p3;
    std::vector<int> p4;

    //define other point variables
    int x3;
    int y3;
    int x4;
    int y4;

    // now look for intersections of segments with it's own street segments + overlap
    for(unsigned int i = 0; i < (segments.size()-2);i = i + 2){
        x3 = segments[i];
        y3 = segments[i+1];
        x4 = segments[i+2];
        y4 = segments[i+3];
        p3.clear();
        p4.clear();
        p3.push_back(x3);
        p3.push_back(y3);
        p4.push_back(x4);
        p4.push_back(y4);

        if(intersect(p1,p2,p3,p4)){
            return true;
        }
        else{
            continue;
        }
    }
    return false;
}

std::vector<int> generateSegment(std::vector<std::vector<int>> streets, std::vector<int> segments,int index,int max_coords){
    std::vector<int> points;
    unsigned int randx2 = Random_init();
    unsigned int randy2 = Random_init();
    int x2 = Rand_Num(randx2,-max_coords,max_coords);
    int y2 = Rand_Num(randy2,-max_coords,max_coords);
    if(segments.size() == 0){
        points.push_back(x2);
        points.push_back(y2);
        return(points);
    }
    else if(segments.size() == 2){
        if(segments[0] != x2 || segments[1] != y2){
            if(!checkoverlap(x2,y2,streets,segments,index)){
                points.push_back(x2);
                points.push_back(y2);
                return(points);
            }
        }
    }
    else{
        if(!checkoverlap(x2,y2,streets,segments,index) && !checkintersect(x2,y2,segments)){
            points.push_back(x2);
            points.push_back(y2);
            return(points);
        }
    }
    return(points);
}

char intToAlphabet( int i ){
   return static_cast<char>('A'+ i);
}

bool ifintersectExists(std::vector<std::vector<int>> streets){
    std::vector<int> p1;
    std::vector<int> p2;
    std::vector<int> p3;
    std::vector<int> p4;
    for(size_t i = 0; i < (streets.size()-1); ++i){
        for(size_t j = i+1; j < streets.size(); ++j){
            for(size_t m = 0; m < (streets[i].size()-2); m = m + 2){
                p1.clear();
                p2.clear();
                p1.push_back(streets[i][m]);
                p1.push_back(streets[i][m+1]);
                p2.push_back(streets[i][m+2]);
                p2.push_back(streets[i][m+3]);
                for(size_t n = 0; n < (streets[j].size()-2); n = n +2){
                    p3.clear();
                    p4.clear();
                    p3.push_back(streets[j][n]);
                    p3.push_back(streets[j][n+1]);
                    p4.push_back(streets[j][n+2]);
                    p4.push_back(streets[j][n+3]);
                    if(p4[0] == p1[0] && p4[1] == p1[1]){
                        return true;
                    }
                    else{
                        bool flag = (intersect(p1,p2,p3,p4));
                        if(flag == true){
                            return flag;
                        }
                        else{
                            continue;
                        }
                    }
                }
            }
        }
    }
    return false;
}

int main(int argc,char **argv){

	// Define default values for program
	int numstreets_k = 10;
	int linesegments_k = 5;
	int wait_time_k = 5;
	int max_coords = 20;
    int c;
    int count = 0;

    int numstreets;
    int wait_time;
    opterr = 0;

    // use switch cases to grab command line arguments from the function call
	while ((c = getopt (argc, argv, "s:n:l:c:")) != -1)
    	switch (c){
    		case 's':
    			numstreets_k = atoi(optarg);
                if((numstreets_k >= 2)){
                    break;
                }
                else{
                    std::cerr << "Error: parameter entered for option s invalid" << std::endl;
                    return 1;
                }
      		case 'n':
        		linesegments_k = atoi(optarg);
                if((linesegments_k >= 1)){
                    break;
                }
                else{
                    std::cerr << "Error: parameter entered for option n invalid" << std::endl;
                    return 1;
                }
      		case 'l':
        		wait_time_k = atoi(optarg);
                if((wait_time_k >= 5)){
                    break;
                }
                else{
                    std::cerr << "Error: parameter entered for option l invalid" << std::endl;
                    return 1;
                }
      		case 'c':
      			max_coords = atoi(optarg);
                if((max_coords >= 1)){
                    break;
                }
                else{
                    std::cerr << "Error: parameter entered for option c invalid" << std::endl;
                    return 1;
                }
            case '?':
                if(optopt == 's'){
                    std::cerr << "Error: option -s requires an argument" << std::endl;
                    return 1;
                }

                else if(optopt == 'n'){
                    std::cerr << "Error: option -n requires an argument" << std::endl;
                    return 1;
                }

                else if(optopt == 'l'){
                    std::cerr << "Error: option -l requires an argument" << std::endl;
                    return 1;
                }
 
                else if(optopt == 'c'){
                    std::cerr << "Error: option -c requires an argument" << std::endl;
                    return 1;
                }
                else{
                    std::cerr << "Error: unknown option " << std::endl;
                    return 1;
                }
      		default:
        		return 0;
    }
    while(!std::cin.eof()){
        //// generate random street amount
        unsigned int rand = Random_init();
        numstreets = Rand_Num(rand,2,numstreets_k);
        std::vector<std::vector<int>> streets;
        std::vector<string> Names;
        for(int i = 0; i < numstreets; ++i){
            //std::cout << "\nInitial street loop.";
            // generate random # of segments this street should have
            unsigned int rand = Random_init();
            int Num_linesegments = Rand_Num(rand,1,linesegments_k);
            //std::cout << "\n" << Num_linesegments << " linesegments for this street.";
            std::vector<int> segments;
            for(int j = 0; j < (Num_linesegments+1); ++j){
                count = 0;
                while(count < 25){
                    // generate a random point and add it to the street segment list
                    std::vector<int> points = generateSegment(streets,segments,i,max_coords);
                    if(points.size()>0){
                        if((i == (numstreets-1))&&(j==(Num_linesegments))){
                            segments.push_back(points[0]);
                            segments.push_back(points[1]);
                            streets.push_back(segments);
                            if(ifintersectExists(streets)){
                                streets.pop_back();
                                break;
                            }
                            else{
                                streets.pop_back();
                                segments.pop_back();
                                segments.pop_back();
                            }
                        }
                        else{
                            segments.push_back(points[0]);
                            segments.push_back(points[1]);
                            break;
                        }
                    }
                    count = count + 1;
                }
                if(count >= 25){
                    std::cerr << "Error: Unable to generate specification after 25 attempts" << std::endl;
                    return(1);                   
                }
                else{
                    continue;
                }
            }
            std::string output;
            for(size_t m = 0; m < segments.size(); m = m + 2){
                output.append("(");
                output.append(std::to_string(segments[m]));
                output.append(",");
                output.append(std::to_string(segments[m+1]));
                if ((m+1) != segments.size() - 1) {
                    output.append(") ");
                }
                else {
                    output.append(")");
                }
            }
            std::string sletter(1,intToAlphabet(i));
            std::string streetName = "\"Bubba" + sletter + "\" ";
            Names.push_back(streetName);
            std::cout << "a " << streetName << output << std::endl;

            output.clear();   
            //std::cout << "\n# of points on street " << i << ": " << segments.size()/2;
            streets.push_back(segments);
            segments.clear();
        }
        std::cout << "g" << std::endl;

        unsigned int wait = Random_init();
        wait_time = Rand_Num(wait,5,wait_time_k);
        sleep(wait_time);

        for (size_t n = 0; n < Names.size(); ++n) {
            std::cout << "r " << Names[n] << std::endl;
        }
        Names.clear();

        // std::cout << "\n# of streets: " << streets.size();
        // std::cout << "\n";

    }
	return 0;
}