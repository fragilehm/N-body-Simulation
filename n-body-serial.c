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

void integrate(body *planet, float deltaTime);
void calculateNewtonGravityAcceleration(body *a, body *b, float *ax, float *ay);
void simulateWithBruteforce(int nBodies, body *bodies, float dt);
void initializeBodies (int nBodies, body *bodies);
float randValue();

int main(int argc, char **argv) {
	float dt = 0.01;
	int nBodies = 100;
	if (argc > 2) {
		dt = atof(argv[1]);
		nBodies = atoi(argv[2]);
	}
	printf("%d\n", argc);
	body *bodies = (body*) malloc(nBodies * sizeof(*bodies));
	initializeBodies(nBodies, bodies);
	simulateWithBruteforce(nBodies, bodies, dt);
	return 0;
}
void integrate(body *body, float deltaTime) {
	body->x += body->vx * deltaTime + (1 / 2) * body->ax * deltaTime * deltaTime;
	body->y += body->vy * deltaTime + (1 / 2) * body->ay * deltaTime * deltaTime;
	body->vx += body->ax * deltaTime;
	body->vy += body->ay * deltaTime;
}
void calculateNewtonGravityAcceleration(body *a, body *b, float *ax, float *ay) {
	float distanceX = fabsf(b->x - a->x);
	float distanceY = fabsf(b->y - a->y);
	float vectorDistance = sqrt(a->x * a->x + a->y * a->y);
	float vectorDistanceCubed = vectorDistance * vectorDistance * vectorDistance;
    float inverse = 1.0 / vectorDistanceCubed;
    float scale = b->mass * inverse;
	*ax = (distanceX * scale);
	*ay = (distanceY * scale);
}
void simulateWithBruteforce(int nBodies, body *bodies, float dt) {
	double totalTime = 0.0; 
	for(size_t i = 0; i < nBodies; i++) {
		double start = clock();
		float total_ax = 0, total_ay = 0;
		for (size_t j = 0; j < nBodies; j++) {
			if (i == j) {
				continue;
			}
			float ax, ay;
			calculateNewtonGravityAcceleration(&bodies[i], &bodies[j], &ax, &ay);
			total_ax += ax;
			total_ay += ay;
		}
		bodies[i].ax = total_ax;
		bodies[i].ay = total_ay;
		integrate(&bodies[i], dt);
		double timeElapsed = ((double) clock()  - start) / CLOCKS_PER_SEC;
		totalTime += timeElapsed;
	}
	printf("%f\n", totalTime);
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
		body object = {	
			.x = randValue(), .y = randValue(),
			.vx = cos(angle) * accelerationScale * randValue(), 
			.vy = sin(angle) * accelerationScale * randValue(), 
			.mass = randValue() + initialMass * 0.5
		};
        bodies[i] = object;
    }
}


