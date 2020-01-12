#ifndef UTILS_H_
#define UTILS_H_

#define READ "rb"
#define WRITE "wb"

#define PGM "P5"
#define PNM "P6"

#define SIZE 50
#define NEWLINE_CODE 10

#define SMOOTH_FILTER "smooth"
#define BLUR_FILTER "blur"
#define SHARPEN_FILTER "sharpen"
#define MEAN_FILTER "mean"
#define EMBOSS_FILTER "emboss"

#define MASTER 0

typedef struct {
	unsigned char r;
	unsigned char g;
	unsigned char b;
} Color;

typedef union {
	unsigned char val;
	Color c;
} Pixel;

float SMOOTH[3][3] = {
	{1 / 9.0f, 1 / 9.0f, 1 / 9.0f},
	{1 / 9.0f, 1 / 9.0f, 1 / 9.0f},
	{1 / 9.0f, 1 / 9.0f, 1 / 9.0f}
};

float BLUR[3][3] = {
	{1 / 16.0f, 1/ 8.0f, 1 / 16.0f},
	{1 / 8.0f, 1 / 4.0f, 1 / 8.0f},
	{1 / 16.0f, 1 / 8.0f, 1 / 16.0f}
};

float SHARPEN[3][3] = {
	{0, -2 / 3.0f, 0},
	{-2 / 3.0f, 11 / 3.0f, -2 / 3.0f},
	{0, -2 / 3.0f, 0}
};

float MEAN[3][3] = {
	{-1, -1, -1},
	{-1, 9, -1},
	{-1, -1, -1}
};

float EMBOSS[3][3] = {
	{0, 1, 0},
	{0, 0, 0},
	{0, -1, 0}
};

#endif // UTILS_H_