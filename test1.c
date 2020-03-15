#include <stdio.h>
#include <math.h>
#include<time.h>
#include<stdlib.h>
#include "dislin.h"

double get_rand(void)
{
	srand(time(NULL));
	//This fuction will be used to generat random doubles between 0 and 1.
	int num = rand() % 10001;
	return (double)num * 0.0001;
}

int main()
{
	//srand(time(NULL));
	int small = 0, big = 0;
	{
		double j = get_rand();
		{
			small ++;
		else if (j > 0.5)
		{
			big ++;
		}
	}
	printf("small= %d, big = %d\n", small, big);
	
	return 0;
}