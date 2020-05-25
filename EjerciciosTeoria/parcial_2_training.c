#define CACHE_LINE_SIZE 32

int *array; // pointer to start of array
int N;      // length of array, assumed to be multiple of the number of processors

void swap(int i, int j)
{
    int tmp = array[i];
    array[i] = array[j];
    array[j] = tmp;
}

// GEOMETRIC DATA DECOMPOSITION
void reverse()
{
    #pragma omp parallel
    {
        int P = omp_get_num_threads();      // number of threads executing this function
        int id = omp_get_thread_num();      // identifier of the thread executing this instance (0 .. P-1)
        int segmentLength = N / (2 * P);    // Divide entre dos porque es un swap entre el simetrico a partir de N/2.
        int segmentStart = id * segmentLength;
        for (int i = segmentStart; i < segmentStart + segmentLength; i++)
            swap(i, N - 1 - i);
    }
}

// BLOCK-CYCLIC GEOMETRIC DATA DECOMPOSITION V1
void compute()
{
    #pragma omp parallel
    {
        int P = omp_get_num_threads();
        int id= omp_get_thread_num();
        int block_size = CACHE_LINE_SIZE / sizeof(int);

        // loop jumping blocks cyclically
        for (int ii = id*block_size; ii < N; ii+=(P*block_size))
            // loop traversing each block
            for (int i = ii; ii < min(N, ii+block_size); i++)
                array[i] = foo(array[i], i);
    }
}

// BLOCK-CYCLIC GEOMETRIC DATA DECOMPOSITION V2
void compute()
{
    #pragma omp parallel
    {
        int P = omp_get_num_threads();
        int id= omp_get_thread_num();

        int segmentLength = N / (2 * P);
        int segmentStart1 = id * segmentLength;
        int segmentStart2 = N/2 + ((P-1-id) * segmentLength);

        for (int i = segmentStart1; i < segmentStart1+segmentLength; i++) {
            array[i] = foo(array[i], i);
            array[N - 1 - i] = foo(array[N - 1 - i], N - 1 - i);
        }
    }
}

void main()
{
    reverse();
    compute();
}