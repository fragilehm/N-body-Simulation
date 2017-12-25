#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
//#include <mpi.h>
#define PI 3.14159265358979323846

typedef struct _body {
	float x, y; // position
	float ax, ay; // acceleration
	float vx, vy; // velocity
	float mass; // mass
} body;

void initializeBodies (int nBodies, body *bodies);
float randValue();

int main(int argc, char **argv) {
	float dt = atof(argv[1]);
	int nBodies = atoi(argv[2]);
	body *bodies = (body*) malloc(nBodies * sizeof(*bodies));
	initializeBodies(nBodies, bodies);
	return 0;
}
float randValue(){
	return ((float) rand() / RAND_MAX);
}
void initializeBodies (int nBodies, body *bodies) {
	srand(time(NULL));
	const float accelerationScale = 100.0;
	for (int i = 0; i < nBodies; i++) {
		float angle = 	((float) i / nBodies) * 2.0 * PI + 
						((randValue() - 0.5) * 0.5);
		float initialMass = 2;
		body object = 	{	
							.x = randValue(), .y = randValue(),
					  		.vx = cos(angle) * accelerationScale * randValue(), 
							.vy = sin(angle) * accelerationScale * randValue(), 
							.mass = randValue() + initialMass * 0.5
					  	};
        bodies[i] = object;
    }
}




























