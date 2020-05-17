#include <mpi.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <png.h>

#define PNG_NO_SETJMP

const int MAX_ITER = 10000;

void write_png(const char *filename, const size_t width, const size_t height, const int *buffer);
void eval_correct(int *res_imag, int width, int height, double upper, double lower, double right, double left);

int main(int argc, char *argv[])
{
    MPI_Init(NULL, NULL);
#ifdef TIME
    double local_begin_time = MPI_Wtime();
    double local_end_time, min_begin_time, max_end_time;
#endif
    // int num_threads = strtol(argv[1], 0, 10);
    double left = strtod(argv[2], 0);
    double right = strtod(argv[3], 0);
    double lower = strtod(argv[4], 0);
    double upper = strtod(argv[5], 0);
    int width = strtol(argv[6], 0, 10);
    int height = strtol(argv[7], 0, 10);
    int local_height;
    const char *filename = argv[8];

    int *image = NULL;
    int *local_image = NULL;
    MPI_Status mpi_status;

    // local vars for mpi
    int my_rank, comm_size;
    int height_start, height_end;

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    int height_step = height / comm_size;
    if (my_rank == 0)
    {
#ifdef DEBUG
        printf("Debug Mode:\n");
#endif
#ifdef TIME
        printf("Record Time\n");
#endif
        height_start = 0;
        local_height = height - (comm_size - 1) * height_step;
        height_end = local_height;
        image = (int *)malloc(width * height * sizeof(int));
    }
    else
    {
        local_height = height_step;
        height_start = (height - (comm_size - 1) * height_step) + (my_rank - 1) * local_height;
        height_end = height_start + height_step;
    }

    local_image = (int *)malloc(width * height_step * sizeof(int));

    for (int j = height_start; j < height_end; ++j)
    {
        double y0 = j * ((upper - lower) / height) + lower;
        for (int i = 0; i < width; ++i)
        {
            double x0 = i * ((right - left) / width) + left;

            int repeats = 0;
            double x = 0;
            double y = 0;
            double length_squared = 0;
            for (; repeats < MAX_ITER && length_squared < 4; ++repeats)
            {
                double temp = x * x - y * y + x0;
                y = 2 * x * y + y0;
                x = temp;
                length_squared = x * x + y * y;
            }
            if (my_rank == 0)
            {
                image[j * width + i] = repeats;
            }
            else
                local_image[(j - height_start) * width + i] = repeats;
        }
    }

    if (my_rank == 0)
    {
        int active_proc_count = comm_size - 1;
        while (active_proc_count--)
        {
            MPI_Recv(local_image, height_step * width, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &mpi_status);
            memcpy(&image[(local_height + (mpi_status.MPI_SOURCE - 1) * height_step) * width], local_image, sizeof(int) * height_step * width);
        }
    }
    else
    {
        MPI_Send(local_image, height_step * width, MPI_INT, 0, 1, MPI_COMM_WORLD);
    }

    /* draw and cleanup */
    if (my_rank == 0)
    {
#ifdef DEBUG
        eval_correct(image, width, height, upper, lower, right, left);
#endif
        write_png(filename, width, height, image);
    }
    free(image);
    free(local_image);

#ifdef TIME
    local_end_time = MPI_Wtime();
    printf("Proc %d finished in %lf seconds\n", my_rank, local_end_time - local_begin_time);
    MPI_Reduce(&local_begin_time, &min_begin_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_end_time, &max_end_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
#endif

    MPI_Finalize();

#ifdef TIME
    if (my_rank == 0)
    {
        printf("mpi_static program finished in %lf seconds\n", max_end_time - min_begin_time);
    }
#endif
    return 0;
}

void write_png(const char *filename, const size_t width, const size_t height, const int *buffer)
{
    FILE *fp = fopen(filename, "wb");
    assert(fp);
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    assert(png_ptr);
    png_infop info_ptr = png_create_info_struct(png_ptr);
    assert(info_ptr);
    png_init_io(png_ptr, fp);
    png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_write_info(png_ptr, info_ptr);
    size_t row_size = 3 * width * sizeof(png_byte);
    png_bytep row = (png_bytep)malloc(row_size);
    for (int y = 0; y < height; ++y)
    {
        memset(row, 0, row_size);
        for (int x = 0; x < width; ++x)
        {
            int p = buffer[(height - 1 - y) * width + x];
            png_bytep color = row + x * 3;
            if (p != MAX_ITER)
            {
                if (p & 16)
                {
                    color[0] = 240;
                    color[1] = color[2] = p % 16 * 16;
                }
                else
                {
                    color[0] = p % 16 * 16;
                }
            }
        }
        png_write_row(png_ptr, row);
    }
    free(row);
    png_write_end(png_ptr, NULL);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
}

void eval_correct(int *res_imag, int width, int height, double upper, double lower, double right, double left)
{
    /* allocate memory for image */
    int *image = (int *)malloc(width * height * sizeof(int));
    assert(image);

    /* mandelbrot set */
    for (int j = 0; j < height; ++j)
    {
        double y0 = j * ((upper - lower) / height) + lower;
        for (int i = 0; i < width; ++i)
        {
            double x0 = i * ((right - left) / width) + left;

            int repeats = 0;
            double x = 0;
            double y = 0;
            double length_squared = 0;
            for (; repeats < MAX_ITER && length_squared < 4; ++repeats)
            {
                double temp = x * x - y * y + x0;
                y = 2 * x * y + y0;
                x = temp;
                length_squared = x * x + y * y;
            }
            image[j * width + i] = repeats;
        }
    }
    double error_count = 0;
    for (int i = 0; i < width * height; i++)
    {
        if (res_imag[i] != image[i])
        {
            error_count++;
        }
    }
    double error_rate = error_count / (width * height);
    printf("All errors: %lf, error rate compared with sequential program: %lf\n", error_count, error_rate);
    free(image);
}