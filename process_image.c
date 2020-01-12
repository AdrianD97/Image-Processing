#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utils.h"

/*
	Elibereaza resursele(memoria) utilizate
*/
void freeMemory(int n, Pixel** image) {
	for (int i = 0; i < n; ++i) {
		free(image[i]);
	}
	free(image);
}

/*
	Functia citeste o imagine alb-negru.
*/
void readPGMImage(FILE* file, int width, int height, Pixel** image) {
	int value = 0;

	for (int i = 0; i < width; ++i) {
		image[0][i].val = image[height - 1][i].val = 0;
	}

	int i = 0, j = width - 1;

	do {
		if (j == width - 1) {
			++i;
			j = 1;
			image[i][0].val = image[i][width - 1].val = 0;
		}
		value = fgetc(file);

		if (value == EOF) {
			break;
		}

		image[i][j].val = (unsigned char)value;
		++j;

	} while (1);
}

/*
	Functia citeste o imagine colorata.
*/
void readPNMImage(FILE* file, int width, int height, Pixel** image) {
	Color c, helpColor;
	int value = 0;
	helpColor.r = helpColor.g = helpColor.b = 0;

	for (int i = 0; i < width; ++i) {
		image[0][i].c = helpColor;
		image[height - 1][i].c = helpColor;
	}

	int i = 0, j = width - 1;

	do {
		if (j == width - 1) {
			++i;
			j = 1;
			image[i][0].c = image[i][width - 1].c = helpColor;
		}

		value = fgetc(file);
		if (value == EOF) {
			break;
		}

		c.r = (unsigned char)value;
		value = fgetc(file);
		c.g = (unsigned char)value;
		value = fgetc(file);
		c.b = (unsigned char)value;

		image[i][j].c = c;
		++j;

	} while (1);
}

/*
	Functia citeste o imagine din fisier si o salveaza in memorie.
	imageType = 5 => imagine al-negru
	imageType = 6 => imagine colorata
*/
Pixel** readImage(char* fileName, unsigned char* imageType, int* width, int* height,
		       unsigned char* maxVal) {

	Pixel** image = NULL;
	FILE* file;
	file = fopen(fileName, READ);

	// open file
	if (file == NULL) {
		printf("Could not open %s file.\n", fileName);
		MPI_Finalize();
		exit(-1);
	}

	// read image type
	char type[3];
	fscanf(file, "%s", type);
	if (strcmp(type, PGM) == 0) {
		*imageType = 5;
	} else if (strcmp(type, PNM) == 0){
		*imageType = 6;
	} else {
		printf("Incorrect image type.\n");
		MPI_Finalize();
		exit(-1);
	}
	fgetc(file);

	// read commented line
	while (fgetc(file) != NEWLINE_CODE) {}

	// read image dimensions
	int w, h;
	fscanf(file, "%d", &w);
	*width = w + 2;
	fscanf(file, "%d", &h);
	*height = h + 2;

	// alloc image memory
	image = (Pixel**)malloc(*height * sizeof(Pixel*));
	for (int i = 0; i < *height; ++i) {
		image[i] = (Pixel*)malloc(*width * sizeof(Pixel));
	}

	// read max value
	int max = 0;
	fscanf(file, "%d", &max);
	*maxVal = (unsigned char)max;

	fgetc(file);
	if (*imageType == 5) {
		readPGMImage(file, *width, *height, image);
	} else {
		readPNMImage(file, *width, *height, image);
	}

	// close file
	if (fclose(file)) {
		printf("An error occured when try to close the %s file.\n", fileName);
	}

	return image;
}

/*
	Salveaza o imagine alb-negru din memorie in fisier
*/
void writePGMImage(FILE* file, int width, int height, unsigned char maxVal, Pixel** image) {
	fprintf(file, "%s\n", PGM);

	fprintf(file, "%d %d\n", width - 2, height - 2);
	fprintf(file, "%u\n", maxVal);

	for (int i = 1; i < height - 1; ++i) {
		for (int j = 1; j < width - 1; ++j) {
			fputc(image[i][j].val, file);
		}
	}
}

/*
	Salveaza o imagine color din memorie in fisier
*/
void writePNMImage(FILE* file, int width, int height, unsigned char maxVal, Pixel** image) {
	fprintf(file, "%s\n", PNM);

	fprintf(file, "%d %d\n", width - 2, height - 2);
	fprintf(file, "%u\n", maxVal);

	for (int i = 1; i < height - 1; ++i) {
		for (int j = 1; j < width - 1; ++j) {
			fputc(image[i][j].c.r, file);
			fputc(image[i][j].c.g, file);
			fputc(image[i][j].c.b, file);
		}
	}
}

/*
	Salveaza imaginea din memorie in fisier
*/
void writeImage(char* fileName, unsigned char imageType, int width, int height,
		       unsigned char maxVal, Pixel** image) {
	FILE* file;
	file = fopen(fileName, WRITE);

	// open file
	if (file == NULL) {
		printf("Could not open %s file.\n", fileName);
		MPI_Finalize();
		freeMemory(height, image);
		exit(-1);
	}

	switch (imageType) {
		case 5:
			writePGMImage(file, width, height, maxVal, image);
			break;
		case 6:
			writePNMImage(file, width, height, maxVal, image);
			break;
		default: printf("Wrong image type.\n");
	}

	// close file
	if (fclose(file)) {
		printf("An error occured when try to close the %s file.\n", fileName);
	}
}

/*
	Creaza un tip de date corespunzator structurii Color.
	Tipul nou este utilizat in transferul de date intre task-uri.
*/
MPI_Datatype createColorType() {
	MPI_Datatype newType;
	MPI_Datatype types[3] = {MPI_UNSIGNED_CHAR, MPI_UNSIGNED_CHAR, MPI_UNSIGNED_CHAR};
	int blocklen[3] = {1, 1, 1};
	MPI_Aint disp[3];

	Color color;

	MPI_Get_address(&color.r, disp); 
	MPI_Get_address(&color.g, disp + 1);
	MPI_Get_address(&color.b, disp + 2);

	disp[1] -= disp[0];
	disp[2] -= disp[0];
	disp[0] = 0;

	MPI_Type_create_struct(3, blocklen, disp, types, &newType);
	MPI_Type_commit(&newType);

	return newType; 
}

/*
	Creaza un tip de date corespunzator structurii Pixel.
	Tipul nou este utilizat in transferul de date intre task-uri.
*/
void createPixelType(MPI_Datatype types[2]) {
	int blocklen[2] = {1, 1};
	MPI_Aint disp[2];
	MPI_Aint i, j;

	MPI_Datatype type[2];
	Pixel image[2];

	MPI_Get_address(image, &i);
	MPI_Get_address(image + 1, &j);
	disp[0] = 0;
	disp[1] = j - i;
	type[1] = MPI_UB;

	type[0] = MPI_UNSIGNED_CHAR;
	MPI_Type_create_struct(2, blocklen, disp, type, &types[0]);

	type[0] = createColorType();
	MPI_Type_create_struct(2, blocklen, disp, type, &types[1]);

	MPI_Type_commit(&types[0]);
	MPI_Type_commit(&types[1]);
}

/*
	Functia aplica filtrul corespunzator asupra imaginii primite ca parametru
*/
void applyFilter(float filterMatrix[3][3], unsigned char imageType, unsigned char maxVal,
			     Pixel** image, int start, int end, int width, Pixel** newImage) {
	float result;
	for (int i = start; i <= end; ++i) {
		for (int j = 1; j < width - 1; ++j) {
			if (imageType == 5) {
				result = filterMatrix[2][2] * image[i - 1][j - 1].val +
								  filterMatrix[2][1] * image[i - 1][j].val +
								  filterMatrix[2][0] * image[i - 1][j + 1].val +
								  filterMatrix[1][2] * image[i][j - 1].val +
								  filterMatrix[1][1] * image[i][j].val +
								  filterMatrix[1][0] * image[i][j + 1].val +
								  filterMatrix[0][2] * image[i + 1][j - 1].val +
								  filterMatrix[0][1] * image[i + 1][j].val +
								  filterMatrix[0][0] * image[i + 1][j + 1].val;
				result = (result < 0.0f) ? 0.0f : ((result > maxVal) ? maxVal : result);
				newImage[i - 1][j - 1].val = (unsigned char)result;
			} else if (imageType == 6) {
				result = filterMatrix[2][2] * image[i - 1][j - 1].c.r +
								  filterMatrix[2][1] * image[i - 1][j].c.r +
								  filterMatrix[2][0] * image[i - 1][j + 1].c.r +
								  filterMatrix[1][2] * image[i][j - 1].c.r +
								  filterMatrix[1][1] * image[i][j].c.r +
								  filterMatrix[1][0] * image[i][j + 1].c.r +
								  filterMatrix[0][2] * image[i + 1][j - 1].c.r +
								  filterMatrix[0][1] * image[i + 1][j].c.r +
								  filterMatrix[0][0] * image[i + 1][j + 1].c.r;
				result = (result < 0.0f) ? 0.0f : ((result > maxVal) ? maxVal : result);
				newImage[i - 1][j - 1].c.r = (unsigned char)result;

				result = filterMatrix[2][2] * image[i - 1][j - 1].c.g +
								  filterMatrix[2][1] * image[i - 1][j].c.g +
								  filterMatrix[2][0] * image[i - 1][j + 1].c.g +
								  filterMatrix[1][2] * image[i][j - 1].c.g +
								  filterMatrix[1][1] * image[i][j].c.g +
								  filterMatrix[1][0] * image[i][j + 1].c.g +
								  filterMatrix[0][2] * image[i + 1][j - 1].c.g +
								  filterMatrix[0][1] * image[i + 1][j].c.g +
								  filterMatrix[0][0] * image[i + 1][j + 1].c.g;
				result = (result < 0.0f) ? 0.0f : ((result > maxVal) ? maxVal : result);
				newImage[i - 1][j - 1].c.g = (unsigned char)result;

				result = filterMatrix[2][2] * image[i - 1][j - 1].c.b +
								  filterMatrix[2][1] * image[i - 1][j].c.b +
								  filterMatrix[2][0] * image[i - 1][j + 1].c.b +
								  filterMatrix[1][2] * image[i][j - 1].c.b +
								  filterMatrix[1][1] * image[i][j].c.b +
								  filterMatrix[1][0] * image[i][j + 1].c.b +
								  filterMatrix[0][2] * image[i + 1][j - 1].c.b +
								  filterMatrix[0][1] * image[i + 1][j].c.b +
								  filterMatrix[0][0] * image[i + 1][j + 1].c.b;
				result = (result < 0.0f) ? 0.0f : ((result > maxVal) ? maxVal : result);
				newImage[i - 1][j - 1].c.b = (unsigned char)result;
			}
		}
	}
}

/*
	Identifica filtrul care se doreste a se aplica imaginii date.
*/
void apply(char* filter, unsigned char imageType, unsigned char maxVal,
		   int start, int end, int width, Pixel** image, Pixel** newImage) {
	if (strcmp(filter, SMOOTH_FILTER) == 0) {
		applyFilter(SMOOTH, imageType, maxVal, image, start, end, width, newImage);
		return;
	}

	if (strcmp(filter, BLUR_FILTER) == 0) {
		applyFilter(BLUR, imageType, maxVal, image, start, end, width, newImage);
		return;
	}

	if (strcmp(filter, SHARPEN_FILTER) == 0) {
		applyFilter(SHARPEN, imageType, maxVal, image, start, end, width, newImage);
		return;
	}

	if (strcmp(filter, MEAN_FILTER) == 0) {
		applyFilter(MEAN, imageType, maxVal, image, start, end, width, newImage);
		return;
	}

	if (strcmp(filter, EMBOSS_FILTER) == 0) {
		applyFilter(EMBOSS, imageType, maxVal, image, start, end, width, newImage);
		return;
	}
}

/*
	Salveaza o copie a unei imagini alb-negru in memorie
*/
void copyPGMImage(Pixel** image, Pixel** imageCopy, int width, int height) {
	for (int i = 1; i < height - 1; ++i) {
		for (int j = 1; j < width - 1; ++j) {
			imageCopy[i][j].val = image[i - 1][j - 1].val;
		}
	}
}

/*
	Salveaza o copie a unei imagini color in memorie
*/
void copyPNMImage(Pixel** image, Pixel** imageCopy, int width, int height) {
	for (int i = 1; i < height - 1; ++i) {
		for (int j = 1; j < width - 1; ++j) {
			imageCopy[i][j].c.r = image[i - 1][j - 1].c.r;
			imageCopy[i][j].c.g = image[i - 1][j - 1].c.g;
			imageCopy[i][j].c.b = image[i - 1][j - 1].c.b;
		}
	}
}

/*
	Salveaza o copie a unei imagini in memorie
*/
void copyImage(Pixel** image, Pixel** imageCopy, unsigned char imageType,
			   int width, int height) {
	switch (imageType) {
		case 5:
			copyPGMImage(image, imageCopy, width, height);
			break;
		case 6:
			copyPNMImage(image, imageCopy, width, height);
			break;
		default: printf("Wrong image type.\n");
	}
}

int main(int argc, char* argv[]) {
	int rank, dim;
	int numberOfProcesses;
	int width, height;
	unsigned char imageType;
	unsigned char maxVal = 0;
	Pixel** image;

	MPI_Init(&argc, &argv);

	MPI_Datatype types[2];
	MPI_Datatype* type;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);

	createPixelType(types);

	if (rank == MASTER) {
		int start, end;
		int result;
		Pixel** imageHelper1;
		Pixel** imageHelper2;

		image = readImage(argv[1], &imageType, &width, &height, &maxVal);

		result = (int)ceil(height * 1.0 / numberOfProcesses * 1.0);
		if (result >= height) {
			result = height - 1;
		}

		imageHelper1 = (Pixel**)malloc((result + 1) * sizeof(Pixel*));
		for (int i = 0; i < result + 1; ++i) {
			imageHelper1[i] = (Pixel*)malloc(width * sizeof(Pixel));
		}

		for (int i = 0; i < result + 1; ++i) {
			for (int j = 0; j < width; ++j) {
				imageHelper1[i][j] = image[i][j];
			}
		}

		imageHelper2 = (Pixel**)malloc((result - 1) * sizeof(Pixel*));
		for (int i = 0; i < result - 1; ++i) {
			imageHelper2[i] = (Pixel*)malloc((width - 2) * sizeof(Pixel));
		}

		type = (imageType == 5) ? &types[0] : ((imageType == 6) ? &types[1] : NULL);
		if (type == NULL) {
			printf("Wrong image type.\n");
			MPI_Finalize();
			freeMemory(height, image);
			freeMemory(result + 1, imageHelper1);
			freeMemory(result - 1, imageHelper2);
			exit(-1);
		}

		for (int i = 1; i < numberOfProcesses; ++i) {
			start = result * i - 1;
			end = result * (i + 1);
			if (end >= height) {
				end = height - 1;
			}
			dim = end - start + 1;

			MPI_Send(&imageType, 1, MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD);
			MPI_Send(&width, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(&maxVal, 1, MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD);
			MPI_Send(&dim, 1, MPI_INT, i, 0, MPI_COMM_WORLD);

			for (int j = start; j <= end; ++j) {
				MPI_Send(image[j], width, *type, i, 0, MPI_COMM_WORLD);
			}
		}

		for (int k = 3; k < argc; ++k) {
			apply(argv[k], imageType, maxVal, 1, result - 1, width, imageHelper1, imageHelper2);

			if (numberOfProcesses != 1) {
				MPI_Send(imageHelper2[result - 2], width - 2, *type, 1, 0, MPI_COMM_WORLD);
				MPI_Recv(&imageHelper1[result][1], width - 2, *type, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}

			copyImage(imageHelper2, imageHelper1, imageType, width, result + 1);
		}

		if (argc > 3) {
			copyImage(imageHelper2, image, imageType, width, result + 1);
		}

		// get results from tasks
		for (int i = 1; i < numberOfProcesses; ++i) {
			start = result * i - 1;
			end = result * (i + 1);
			if (end >= height) {
				end = height - 1;
			}

			for (int j = start + 1; j < end; ++j) {
				MPI_Recv(&image[j][1], width - 2, *type, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}
		
		writeImage(argv[2], imageType, width, height, maxVal, image);

		freeMemory(height, image);
		freeMemory(result + 1, imageHelper1);
		freeMemory(result - 1, imageHelper2);
	} else {
		Pixel** newImage;

		MPI_Recv(&imageType, 1, MPI_UNSIGNED_CHAR, MASTER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&width, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&maxVal, 1, MPI_UNSIGNED_CHAR, MASTER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&dim, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		image = (Pixel**)malloc(dim * sizeof(Pixel*));
		for (int i = 0; i < dim; ++i) {
			image[i] = (Pixel*)malloc(width * sizeof(Pixel));
		}

		newImage = (Pixel**)malloc((dim - 2) * sizeof(Pixel*));
		for (int i = 0; i < dim - 2; ++i) {
			newImage[i] = (Pixel*)malloc((width - 2) * sizeof(Pixel));
		}

		type = (imageType == 5) ? &types[0] : ((imageType == 6) ? &types[1] : NULL);
		if (type == NULL) {
			printf("Wrong image type.\n");
			MPI_Finalize();
			freeMemory(dim, image);
			freeMemory(dim - 2, newImage);
			exit(-1);
		}

		// get data from master
		for (int i = 0; i < dim; ++i) {
			MPI_Recv(image[i], width, *type, MASTER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		// apply filters
		for (int k = 3; k < argc; ++k) {
			apply(argv[k], imageType, maxVal, 1, dim - 2, width, image, newImage);

			MPI_Recv(&image[0][1], width - 2, *type, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (rank < numberOfProcesses - 1) {
				MPI_Send(newImage[dim - 3], width - 2, *type, rank + 1, 0, MPI_COMM_WORLD);
			}

			MPI_Send(newImage[0], width - 2, *type, rank - 1, 0, MPI_COMM_WORLD);

			if (rank < numberOfProcesses - 1) {
				MPI_Recv(&image[dim - 1][1], width - 2, *type, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}

			copyImage(newImage, image, imageType, width, dim);
		}

		// send the results
		for (int i = 1; i < dim - 1; ++i) {
			MPI_Send(&image[i][1], width - 2, *type, MASTER, 0, MPI_COMM_WORLD);
		}

		freeMemory(dim, image);
		freeMemory(dim - 2, newImage);
	}
	
	MPI_Finalize();
	return 0;
}

