build: process_image.c
	mpicc -Wall $^ -o process_image -lm

clean:
	rm process_image

