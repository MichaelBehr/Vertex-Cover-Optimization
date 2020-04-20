// Compile with:
// g++ -pthread -o threadtest -std=c++11 thread_prototype.cpp

#include <iostream>
#include <pthread.h>
#define NUM_THREADS 3

class DummyClass {
public:
  DummyClass() { this->a = 0; };
  void SomeMethod(void) {
    std::cout << "SomeMethod ran!" << std::endl;
    this->a = 21;
  };
  int a;
};

void *DummyThreadFunction(void *dummy_object) {
  DummyClass *function_dummy_object = (DummyClass *)dummy_object;
  std::cout << "Address of dummy 2:" << function_dummy_object << std::endl;
  function_dummy_object->SomeMethod();
  pthread_exit(NULL);
}

void *function1(void *threadarg) {
  long tid;
  tid = (long)threadarg;
  std::cout << "Hello World! Thread ID, " << tid << std::endl;
  pthread_exit(NULL);
}

void *function2(void *threadarg) {
  long tid;
  tid = (long)threadarg;
  std::cout << "Hello World! Thread ID, " << tid << std::endl;
  pthread_exit(NULL);
}

void *function3(void *threadarg) {
  long tid;
  tid = (long)threadarg;
  std::cout << "Hello World! Thread ID, " << tid << std::endl;
  pthread_exit(NULL);
}

int main(int argc, char *argv[]) {
  DummyClass dummy;
  std::cout << "Address of dummy 1:" << &dummy << std::endl;

  // Create 3 threads
  pthread_t threads[NUM_THREADS];
  int rc[NUM_THREADS];

  std::cout << "Main: creating thread, " << std::endl;
  rc[0] =
      pthread_create(&threads[0], NULL, DummyThreadFunction, (void *)&dummy);

  if (rc[0]) {
    std::cout << "Error:unable to create thread," << rc << std::endl;
  }

  void *status;
  rc[0] = pthread_join(threads[0], &status);
  if (rc[0]) {
    std::cout << "Error:unable to join," << rc << std::endl;
    exit(-1);
  }
  std::cout << "Main: completed thread id :" << threads[0];
  std::cout << "  exiting with status :" << status << std::endl;

  std::cout << dummy.a << std::endl;

  // for (i = 0; i < NUM_THREADS; i++) {
  //   std::cout << "Main: creating thread, " << i << std::endl;
  //   rc[i] = pthread_create(&threads[i], NULL, function1, &threads[i]);

  //   if (rc[i]) {
  //     std::cout << "Error:unable to create thread," << rc << std::endl;
  //   }
  // }
  // free attribute and wait for the other threads
  // void *status;
  // for (i = 0; i < NUM_THREADS; i++) {
  //   rc[i] = pthread_join(threads[i], &status);
  //   if (rc[i]) {
  //     std::cout << "Error:unable to join," << rc << std::endl;
  //     exit(-1);
  //   }
  //   std::cout << "Main: completed thread id :" << i;
  //   std::cout << "  exiting with status :" << status << std::endl;
  // }

  std::cout << "Main: program exiting." << std::endl;
  pthread_exit(NULL);
}
