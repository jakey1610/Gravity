#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <limits>
#include <iomanip>
#include <cmath>
#include <string.h>

int main(int argc, char** argv)
{
  //./assignment 0.01 100.0 1e-8 0 0 0 0 0 0 4 3 0 0 0 0 0 5 3 4 0 0 0 0 3
  char* times[] = {"1e-6","2e-6","3e-6","4e-6","5e-6","6e-6","7e-6","8e-6","9e-6","1e-7","2e-7","3e-7","4e-7","5e-7","6e-7","7e-7","8e-7","9e-7","1e-8","2e-8","3e-8","4e-8"};//,"5e-8","6e-8","7e-8","8e-8","9e-8"};
  char name[] = "./assignment ";
  char child[];
  char* args[] = { "0.01", "1", "0.01", "0", "0", "0", "0", "0", "0", "4" , "3", "0", "0", "0", "0", "0", "5" , "3", "4", "0", "0", "0", "0", "3" };
  for (int i = 0; i < 22; i++){
    args[3] = times[i];
    strcpy(child,name);
    char* argsSt
    strcat(child,args)
    printf(child);
    int pr = system(child);
    if(pr == -1) {
        printf("\nfailed connection\n");
    }
    // return 1;
  }

}
