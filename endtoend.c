#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>


// typedef struct centerOfMass
// {
// 	float x, y, z;
// } CENTER_OF_MASS;

// typedef struct coordinates
// {
// 	float x, y, z;
// } COORDINATES;

int calculateNatoms (FILE *input)
{
	rewind (input);
	int nAtoms = 0;
	char lineString[2000];

	for (int i = 0; i < 4; ++i) {
		fgets (lineString, 2000, input); }

	sscanf (lineString, "%d", &nAtoms);

	return nAtoms;
}

int calculateNtimeframes (FILE *input, int nAtoms)
{
	rewind (input);
	int nTimeframes, nLines = 0;
	char lineString[2000];

	while (fgets (lineString, 2000, input) != NULL)
	{
		nLines++;
	}

	nTimeframes = nLines / (nAtoms + 9);

	return nTimeframes;
}

float **computeEndtoEnddistance (float **endtoenddistance, FILE *input, int nAtoms, int degreeOfPolymerization, int nChains)
{
	
	rewind (input);
	char lineString[2000];
	int currentLine_1 = 0, currentTimestep = 0;
	int currentLine_2 = 0;
	float dx,dy,dz, distance = 0;
	int currentChain = 0;
	float x1,y1,z1,x2,y2,z2 =0;


	while (fgets (lineString, 2000, input) != NULL)
	{
		currentLine_1++;
		//printf("%s", lineString);
		//usleep(100000);

		if (currentLine_1 > 9 && currentLine_1 <= 9 + nAtoms)
		{
		
			if(currentLine_2 < degreeOfPolymerization)
			{

				if (currentLine_2 == 0)
				{
		  			sscanf (lineString, "%*d %*d %*f %*f %*f %f %f %f\n", &x1, &y1, &z1);
					//currentLine_2 ++;
				}

				currentLine_2 ++;
			}

			if(currentLine_2 == degreeOfPolymerization)
			{
				sscanf (lineString, "%*d %*d %*f %*f %*f %f %f %f\n", &x2, &y2, &z2); 
				dx = x2 - x1;
				dy = y2 - y1;
				dz = z2 - z1;
				//printf("%s", lineString);
				//usleep(100000);
			
				distance = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));
				//printf("%f\n", distance);
				//usleep(100000);
				endtoenddistance[currentTimestep][currentChain] = distance;

				currentChain ++;
				x1 = 0;
				y1 = 0;
				z1 = 0;
				x2 = 0;
				y2 = 0;
				z2 = 0;
				currentLine_2 = 0;
				distance = 0;
			}

			if (currentLine_1 == (nAtoms +9))
			{
				currentChain = 0;
				currentLine_1 = 0;
				currentTimestep ++;
			}
		}

	}

	return endtoenddistance;
}

void printEndtoEnd (float **endtoenddistance, int nTimeframes, int nChains, const char *inputFilename)
{
	char *outputFilename;
	outputFilename = (char *) malloc (100 * sizeof (char));
	snprintf (outputFilename, 100, "%s.ete", inputFilename);

	FILE *output;
	output = fopen (outputFilename, "w");

	for (int i = 0; i < nTimeframes; ++i)
	{
		for (int j = 0; j < nChains; ++j)
		{
			fprintf(output, "%d %d %f\t \n ", i, j, endtoenddistance[i][j]);		
		}
		
	}

	fclose (output);
	free (outputFilename);
}


int main (int argc, char const *argv[])
	{
        if (argc == 1) {
		(void)printf("\nARGUMENTS:\n~~~~~~~~~~\n\n   {~} argv[0] = program\n   {~} argv[1] = input dump filename\n   {~} argv[2] = Degree of polymerization\n   {~} argv[3] = Number of chains present in the dump file\n\nHOW THIS WORKS?\n~~~~~~~~~~~~~~~\n\n   The program will calculate the average end to end distance for each chain as a function of time. All molecules of interest must be present initially, followed by other molecules (solvent/nanoparticles). Output is saved as \"argv[1].ete\". Column 1 in the output contains the si.no. for the time (multiply this number by your dump frequency) and column 2 contains the end to end distance values. \n\nNOTE:\n~~~~~\n\n   Download a fresh copy of this code from https://github.com/raghurame/sameer before usage.\n\n");
		exit (1); }
		FILE *input;
		input = fopen (argv[1], "r");
		int degreeOfPolymerization = atoi (argv[2]);
		int nChains = atoi (argv[3]);
		int nAtoms = calculateNatoms (input);
		int nTimeframes = calculateNtimeframes (input, nAtoms);

		float **endtoenddistance;
		endtoenddistance = (float **) malloc (nTimeframes * sizeof (float *));
		for (int i = 0; i < nTimeframes; ++i) 
		{
		endtoenddistance[i] = (float *) malloc (nChains * sizeof (float)); 
	    }

	    endtoenddistance = computeEndtoEnddistance (endtoenddistance, input, nAtoms, degreeOfPolymerization, nChains);

	    // printEndtoEnd (endtoenddistance, nTimeframes, nChains, argv[1]);

	    for (int i = 0; i < nTimeframes; ++i)
	    {
	    	for (int j = 0; j < nChains; ++j)
			{
				printf("%d %d %f\n", i, j, endtoenddistance[i][j] );
	    		usleep(100000);
			}	    	
			
	    }
		return 0;
	} 