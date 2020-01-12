		Image processing

	Image processing using MPI.
	A distributed algorithm that process a pnm/pgm image using filters.

	Available filters:
		- blur
		- smooth
		- sharpen
		- emboss
		- mean

	Usage:
		make
		mpirun -np number_of_processes ./process_image input_image out_image list_of_filters