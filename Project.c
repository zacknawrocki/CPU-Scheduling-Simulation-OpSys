#include <stdio.h>
#include <stdlib.h>

struct process{
	int pid; //process ID
	int AT; //arrival time
	int WT; //Waittime
	int BT; //burst time
	int preemptions;
	int TAT;
};

int main(int argc, char ** argv)
{
	int rand_num = atoi(argv[1]); //random number
	srand48(RandNum); 
	double lambda = atof(argv[2]);
	double tail = atoi(argv[3]);
	int num_of_proc = atoi(argv[4]); //Process IDs are assigned in alphabetical order A through Z. therefore max 26
	int context_switch = atoi(argv[5]);
	double alpha = atoi(argv[6]);
	int time_slice = atoi(argv[7]);
	char * begin_or_end = "END"; //DEFAULT
	if (argv[8]!=NULL) begin_or_end = argv[8];
	
	
	
	
	return 0;
}
