#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

// This code calcukates the radius of gyration for each polymer chain and returns the ensemble averaged value

/*
ARGUMENTS TO PASS:
~~~~~~~~~~~~~~~~~
argv[0] = program
argv[1] = input dump file name
argv[2] = degree of polymerization (number of beads per chain)
argv[3] = number of chains
*/


//Before the main program function, we define all our functions and user defined data structures to be used.

typedef struct centerOfMass
{
	float x, y, z;
} CENTER_OF_MASS;

typedef struct coordinates
{
	float x, y, z;
} COORDINATES;

int calculateNatoms (FILE *input) //Calculates the number of atoms in the dump file.
{
	rewind (input);
	int nAtoms;
	char lineString[2000];

	for (int i = 0; i < 4; ++i) //The atom information is on the 4th line in the dump file. 
     {
		fgets (lineString, 2000, input); //The fgets comes to the 4th line and stops.
     } 

	sscanf (lineString, "%d", &nAtoms); //sscanf scans the line string where fgets stopped, and returns the string there as integer.

	return nAtoms;
}

int calculateNtimeframes (FILE *input, int nAtoms) //Finds the total timeframes in the dump file. Takes the input file and nAtoms as input
{
	rewind (input); //Sets the file position to the starting line of the file.
	int nTimeframes, nLines = 0; //Initializing the variables
	char lineString[2000]; // Choose random high number. The value of '2000' must be larger than the longest line in dump file.

	while (fgets (lineString, 2000, input) != NULL) //Reads a text/line from the file 'input' and stores it in the variable lineString.
	{
		nLines++;
	}

	nTimeframes = nLines / (nAtoms + 9);

	rewind (input);
	return nTimeframes;
}

// chainCOMs = computeCenterOfMass (chainCOMs, input, nAtoms, degreeOfPolymerization, nChains);
CENTER_OF_MASS **computeCenterOfMass (CENTER_OF_MASS **chainCOMs, FILE *input, int nAtoms, int degreeOfPolymerization, int nChains)
{
	rewind (input);

	char lineString[2000];
	int currentLine = 0, currentTimestep = 0, currentArrayCounter = 0;

	COORDINATES *chainCoordinates;
	chainCoordinates = (COORDINATES *) malloc (nAtoms * sizeof (COORDINATES));

	while (fgets (lineString, 2000, input) != NULL)
	{
		currentLine++;

		if (currentLine > 9 && currentLine <= 9 + nAtoms)
		{
			sscanf (lineString, "%*d %*d %*f %*f %*f %f %f %f\n", &chainCoordinates[currentLine - 10].x, &chainCoordinates[currentLine - 10].y, &chainCoordinates[currentLine - 10].z);
		}

		if (currentLine == (nAtoms + 9))
		{
			currentArrayCounter = 0;

			for (int j = 0; j < nChains; ++j)
			{
				chainCOMs[currentTimestep][j].x = 0;
				chainCOMs[currentTimestep][j].y = 0;
				chainCOMs[currentTimestep][j].z = 0;
			}

			for (int j = 0; j < nChains; ++j)
			{
				for (int k = 0; k < degreeOfPolymerization; ++k)
				{
					chainCOMs[currentTimestep][j].x += chainCoordinates[currentArrayCounter].x; 
					chainCOMs[currentTimestep][j].y += chainCoordinates[currentArrayCounter].y;
					chainCOMs[currentTimestep][j].z += chainCoordinates[currentArrayCounter].z;

					currentArrayCounter++;
				}

				chainCOMs[currentTimestep][j].x /= degreeOfPolymerization;
				chainCOMs[currentTimestep][j].y /= degreeOfPolymerization;
				chainCOMs[currentTimestep][j].z /= degreeOfPolymerization;
			}

			currentLine = 0;
			currentTimestep++;
		}
	}

	return chainCOMs;
}


float computeDistanceSq (float x1, float y1, float z1, float x2, float y2, float z2)
{
	float distance;
	distance = pow ((x2 - x1), 2) + pow ((y2 - y1), 2) + pow ((z2 - z1), 2);
	return distance;
}

float **computeRadiusOfGyration (float **radiusOfGyration, FILE *input, int nAtoms, int nTimeframes, int degreeOfPolymerization, CENTER_OF_MASS **chainCOMs, int nChains)
{
	COORDINATES **dump;
	dump = (COORDINATES **) malloc (nChains * sizeof (COORDINATES *));

	float x, y, z, distanceSq;
	int progress = 0, progressCheckFreq = (int) ceil (nChains * nTimeframes * 0.01);

	for (int i = 0; i < nChains;  ++i) {
		dump[i] = (COORDINATES *) malloc (degreeOfPolymerization * sizeof (COORDINATES)); }

	rewind (input);
	char lineString[2000];

	printf("Computing radius of gyration ...\n");
	for (int i = 0; i < nTimeframes; ++i)
	{
		// Skipping the header lines
		for (int j = 0; j < 9; ++j)
		{
			fgets (lineString, 2000, input);
		}
		// Reading the atom entries
		for (int j = 0; j < nChains; ++j)
		{
			for (int k = 0; k < degreeOfPolymerization; ++k)
			{
				fgets (lineString, 2000, input);
				sscanf (lineString, "%*d %*d %*f %*f %*f %f %f %f", &x, &y, &z);
				distanceSq = computeDistanceSq (x, y, z, chainCOMs[i][j].x, chainCOMs[i][j].y, chainCOMs[i][j].z);
				radiusOfGyration[i][j] += distanceSq;

				progress++;
				if ((progress % progressCheckFreq) == 0) {
					printf("%.0f%% completed               \r", (float) (progress * 100)/(nChains * degreeOfPolymerization * nTimeframes));
					fflush (stdout); }
			}

			radiusOfGyration[i][j] /= degreeOfPolymerization;
		}
	}

	printf("100%% completed...                    \n");

	return radiusOfGyration;
}

float *computeTimeavgRg (float *TimeavgRg, float **radiusOfGyration, int nTimeframes, int nChains)
{
	// Initialize the values to zero
	for (int i = 0; i < nTimeframes; ++i)
	{
		TimeavgRg[i] = 0;
	}
	for (int i = 0; i < nTimeframes; ++i)
	{
		for (int j = 0; j < nChains; ++j)
		{
			TimeavgRg[i] += radiusOfGyration[i][j];
		}

		TimeavgRg[i] /= nChains;
	}
	return TimeavgRg;
}

float *computeEnsembleavgRg (float *EnsembleavgRg, float **radiusOfGyration, int nTimeframes, int nChains)
{
	// Initialize the values to zero
	for (int i = 0; i < nChains; ++i)
	{
		EnsembleavgRg[i] = 0;
	}
	for (int i = 0; i < nChains; ++i)
	{
		for (int j = 0; j < nTimeframes; ++j)
		{
			EnsembleavgRg[i] += radiusOfGyration[j][i];
		}

		EnsembleavgRg[i] /= nTimeframes;
	}
	return EnsembleavgRg;
}

void printTimeavgRg (float *TimeavgRg, int nTimeframes, const char *inputFilename, int nChains)
{
	char *outputFilename_1;
	outputFilename_1 = (char *) malloc (100 * sizeof (char));
	snprintf (outputFilename_1, 100, "%s.timeavgrg", inputFilename);

	FILE *output_1;
	output_1 = fopen (outputFilename_1, "w");

	for (int i = 0; i < nTimeframes; ++i)
	{
		fprintf(output_1, "%d %f\n", i, TimeavgRg[i]);
	}

	fclose (output_1);
	free (outputFilename_1);
}

void printEnsemblevgRg (float *EnsembleavgRg, int nTimeframes, const char *inputFilename, int nChains)
{
	char *outputFilename_2;
	outputFilename_2 = (char *) malloc (100 * sizeof (char));
	snprintf (outputFilename_2, 100, "%s.ensembleavgrg", inputFilename);

	FILE *output_2;
	output_2 = fopen (outputFilename_2, "w");

	for (int i = 0; i < nChains; ++i)
	{
		fprintf(output_2, "%d %f\n", i, EnsembleavgRg[i]);
	}

	fclose (output_2);
	free (outputFilename_2);
}

float **zeroRG (float **radiusOfGyration, int nTimeframes, int nChains)
{
	for (int i = 0; i < nTimeframes; ++i)
	{
		for (int j = 0; j < nChains; ++j)
		{
			radiusOfGyration[i][j] = 0;
		}
	}

	return radiusOfGyration;
}

int main(int argc, char const *argv[]) //By default, agrv[0] is the program name
{
	if (argc == 1)
	{
		(void)printf("ARGUMENTS:~~~~~~~~~~\n\n{~} argv[0] = program\n{~} argv[1] = input dump file\n{~} argv[2] = Degree of polymerization\n{~} argv[3] = Number of chains\n{~} argv[4] = Output \n");
		exit(1);
	}
    FILE *input; //This is declaring a pointer called input of type 'File type'. 
    input = fopen (argv[1], "r"); //argv[1] is the input dump file name. It will open the file and assign it to the pointer input.

    int degreeOfPolymerization = atoi (argv[2]); // Since argv[2] is a string, atoi converts it to an integer and stores. 
    int nChains = atoi (argv[3]); //Similar as above.
    int nAtoms = calculateNatoms(input);
    int nTimeframes = calculateNtimeframes(input, nAtoms);

    CENTER_OF_MASS **chainCOMs;
	chainCOMs = (CENTER_OF_MASS **) malloc (nTimeframes * sizeof (CENTER_OF_MASS *));

	for (int i = 0; i < nTimeframes; ++i) {
		chainCOMs[i] = (CENTER_OF_MASS *) malloc (nChains * sizeof (CENTER_OF_MASS)); }

	chainCOMs = computeCenterOfMass (chainCOMs, input, nAtoms, degreeOfPolymerization, nChains);

	float **radiusOfGyration;
	radiusOfGyration = (float **) malloc (nTimeframes * sizeof (float *));

	for (int i = 0; i < nTimeframes; ++i)
	{
		radiusOfGyration[i] = (float *) malloc (nChains * sizeof (float));
	}

	radiusOfGyration = zeroRG (radiusOfGyration, nTimeframes, nChains);

	radiusOfGyration = computeRadiusOfGyration (radiusOfGyration, input, nAtoms, nTimeframes, degreeOfPolymerization, chainCOMs, nChains);

	// for (int i = 0; i < nChains; ++i)
	// {
	// 	printf("%f\n", radiusOfGyration[nTimeframes-1][i]);
	// 	sleep (1);
	// }

	float *TimeavgRg;

	TimeavgRg = (float *) malloc (nTimeframes * sizeof (float));
    TimeavgRg = computeTimeavgRg (TimeavgRg, radiusOfGyration, nTimeframes, nChains);

    float *EnsembleavgRg;

    EnsembleavgRg = (float *) malloc (nChains * sizeof (float));
    EnsembleavgRg = computeEnsembleavgRg (EnsembleavgRg, radiusOfGyration, nTimeframes, nChains);

	printTimeavgRg (TimeavgRg, nTimeframes, argv[1], nChains);
	printEnsemblevgRg (EnsembleavgRg, nTimeframes, argv[1], nChains);

    for (int i = 0; i < nTimeframes; ++i)
    {
    	for (int j = 0; j < nChains; ++j)
    	{
    	 	printf("%d %d %f\n", 	i, j, radiusOfGyration[i][j]);
			usleep(100000);
    	}
    	
    }
	
	free (chainCOMs);
	free (radiusOfGyration);
	free (TimeavgRg);
	fclose (input);

   return (0);

}