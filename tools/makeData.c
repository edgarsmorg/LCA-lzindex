#include <stdio.h>
#include <stdlib.h>
#include "mtwister.h"

#define GENOME_LENGTH 100000
#define NUM_GENOMES 1000
#define GENOME_NOISE 0.001
#define NUM_PATTERNS 1000
#define PATTERN_LENGTH 1000
#define PATTERN_NOISE 0.1


void main(int argc, char *argv[]) {

	MTRand mtx = seedRand(2024);

	char *reference = (char *) malloc(GENOME_LENGTH);

	for (int i = 0; i < GENOME_LENGTH; i++) {
		double x = genRand(&mtx);

		if (x < 0.25) {
			reference[i] = 'A';
		} else if (x < 0.5) {
			reference[i] = 'C';
		} else if (x < 0.75) {
			reference[i] = 'G';
		} else {
			reference[i] = 'T';
		}
	}

	FILE *textFile = fopen(argv[1], "w");

	fprintf(textFile, ">T\n");
	
	for (int i = 0; i < NUM_GENOMES; i++) {
		for (int j = 0; j < GENOME_LENGTH; j++) {
			double x = genRand(&mtx);
			if (x < GENOME_NOISE / 3) {
				fputc('A', textFile);
			} else if (x < 2 * GENOME_NOISE / 3) {
				fputc('C', textFile);
			} else if (x < GENOME_NOISE) {
				fputc('G', textFile);
			} else if (x < 4 * GENOME_NOISE / 3) {
				fputc('T', textFile);
			} else {
				fputc((int) reference[j], textFile);
			}
		}
		fputc('\n', textFile);
	}

	fclose(textFile);

	FILE *patternFile = fopen(argv[2], "w");

	for (int i = 0; i < NUM_PATTERNS; i++) {
		fprintf(patternFile, ">P%i\n", i);
		
		for (int j = 0; j < PATTERN_LENGTH; j++) {
			double x = genRand(&mtx);
			if (x < PATTERN_NOISE / 3) {
				fputc('A', patternFile);
			} else if (x < 2 * PATTERN_NOISE / 3) {
				fputc('C', patternFile);
			} else if (x < PATTERN_NOISE) {
				fputc('G', patternFile);
			} else if (x < 4 * PATTERN_NOISE / 3) {
				fputc('T', patternFile);
			} else {
				fputc((int) reference[j], patternFile);
			}
		}
		fputc('\n', patternFile);
	}

	fclose(patternFile);

	free(reference);

	return;
}
