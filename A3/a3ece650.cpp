#include <iostream>
#include <vector>
#include <unistd.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/wait.h>


int main(int argc, char **argv)
{
    // create pid vector
    std::vector<pid_t> kids;
    // create a pipe
    int toA1[2];
    pipe(toA1);
    int toA2[2];
    pipe(toA2);

    pid_t child_pid;

    child_pid = fork();
    kids.push_back(child_pid);

    // lets do A1 first
    if (child_pid == 0)
    {
        // redirect stdout to the pipe
        dup2(toA1[0],STDIN_FILENO);
        dup2(toA2[1], STDOUT_FILENO);
        char *args[]={(char*)"python3",(char*)"a1ece650.py",(char*)NULL}; 

        close(toA1[1]);
        close(toA1[0]);
        close(toA2[1]);
        close(toA2[0]);  
        int state = execvp("python3",args); 
        if(state == -1){
            std::cerr << "Error: Problem with executing a1ece650" << std::endl;
        }
    }
    else if (child_pid < 0) {
        std::cerr << "Error: could not fork\n";
        return 1;
    }

    // Now we setup A2 pipe connections
    child_pid = fork();
    if (child_pid == 0)
    {
        // redirect stdin from the 2nd pipe
        dup2(toA2[0], STDIN_FILENO);
        close(toA2[1]);
        close(toA2[0]);
        char *args[]={(char*)"./a2ece650",NULL}; 
        // Execute the 2nd assignment
        int state = execvp(args[0],args);
        if(state == -1){
            std::cerr << "Error: Problem with executing a2ece650" << std::endl;
        }
    }
    else if (child_pid < 0) {
        std::cerr << "Error: could not fork\n";
        return 1;
    }

    kids.push_back(child_pid);
    // Now we setup rgen pipe connections
    child_pid = fork();
    if (child_pid == 0)
    {
        // redirect std from the 2nd pipe
        dup2(toA1[1], STDOUT_FILENO);
        close(toA1[0]);
        close(toA1[1]);
        argv[0] = (char*)"rgen";
        // Execute the 2nd assignment
        int state = execvp("./rgen",argv);
        if(state == -1){
            std::cerr << "Error: Problem with executing rgen" << std::endl;
        }
    }
    else if (child_pid < 0) {
        std::cerr << "Error: could not fork\n";
        exit(1);
    }

    kids.push_back(child_pid);

    // Now we setup EXTRA process to collect input through standard input line
    child_pid = fork();
    if (child_pid == 0){
        dup2(toA2[1],STDOUT_FILENO);
        close(toA2[1]);
        close(toA2[0]);
        while (!std::cin.eof()) {
            std::string line;
            getline(std::cin, line);
            if (line.size() > 0) {
                std::cout << line << std::endl;
            }
        }
        return 0;
    }
    else if (child_pid < 0) {
        std::cerr << "Error: could not fork\n";
        exit(1);
    }

    kids.push_back(child_pid);

    // MONITOR PROCESSES
    int monitor;
    wait(&monitor);

    // send kill signal to all children
    for (pid_t k : kids) {
        int status;
        kill (k, SIGTERM);
        waitpid(k, &status, 0);
    }

    // exit with return code of process B
    return 0;

}