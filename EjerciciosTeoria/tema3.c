/*
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
*/

/*
10. Parallelise the following sequential code using task dependences in OpenMP, following a producer-
consumer execution model. Producer and consumer code should be in two different tasks.
*/

float sample[INPUT_SIZE+TAP1];
float coeff1[TAP1], coeff2[TAP2];
float data_out[INPUT_SIZE], final[INPUT_SIZE];

void main() {

    float sum;

    #pragma omp parallel for private(sum)
    for (int i=0; i<INPUT_SIZE; i++) {
        
        // Producer: Finite Impulse Response (FIR) filter
        sum = 0.0;
        #pragma omp taskloop firstprivate(i,sum) reduction(+: sum) // con taskgroup se hace wait para la siguiente suma
        for (int j = 0; j < TAP1; j++)
            sum += sample[i+j] * coeff1[j];
        
        #pragma omp taskwait
        data_out[i] = sum;
        
        // Consumer: apply correction function
        #pragma omp taskloop firstprivate(i) reduction(+: final[i]) depend(in:data_out[i])
        for (int j = 0; j < TAP2; j++)
            final[i] += correction(data_out[i], coeff2[j]);
    }
}

/*
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
*/

/*
12. We ask to parallelize the computation of the histogram of values appearing on a vector. The
histogram is another vector in which each position counts the number of elements in the input vector
that are in a certain value range. The following program shows a possible sequential implementation
for computing the histogram (vector frequency) of the input vector numbers:
*/

#define MAX_ELEM 1024*1024
#define HIST_SIZE 250

unsigned int numbers[MAX_ELEM];
unsigned int frequency[HIST_SIZE];

void ReadNumbers (int * input, int * size);

void FindBounds(int * input, int size, int * min, int * max) {
    for (int i=0; i<size; i++)
        if (input[i]>(*max)) (*max)=input[i];

    for (int i=0; i<size; i++)
        if (input[i]<(*min)) (*min)=input[i];
}

void FindFrequency(int * input, int size , int * histogram, int min, int max) {
    int tmp;
    for (int i=0; i<size; i++) {
        tmp = (input[i] - min) * (HIST_SIZE / (max - min - 1));
        histogram[tmp]++;
    }
}

void DrawHistogram(int * histogram, int minimum, int maximum);

void main() {
    
    int num_elem, max, min;
    
    ReadNumbers(numbers, &num_elem); // read input numbers
    
    max = min = numbers[0];

    FindBounds(numbers, num_elem, &min, &max); // returns the upper and lower

    // values for the histogram
    FindFrequency(numbers, num_elem, frequency, min, max); // compute histogram
    DrawHistogram(frequency, min, max); // print the histogram
}

/*
(a) Write a parallel version for function FindBounds ONLY using one OpenMP #pragma line, in
such a way that you minimize the possible parallel execution and synchronization overheads.
*/

void FindBounds(int * input, int size, int * min, int * max) {
    
    /* 
     * Modificando el codigo podemos hacer una version recursiva que vaya se quede siempre
     * con el elemento mas grande y mas pequenyo que vea en cada nodo (subiendo de la recursion)
    */

    for (int i=0; i<size; i++)
        if (input[i]>(*max)) (*max)=input[i];

    for (int i=0; i<size; i++)
        if (input[i]<(*min)) (*min)=input[i];
}

/*
(b) Write an alternative implementation for the following parallel version of function FindFrequency:
in which you improve the parallelism that is achieved, JUST using OpenMP pragmas.
*/
void FindFrequency(int * input, int size , int * histogram, int min, int max) {
    
    int tmp;

    #pragma omp parallel for private(tmp)
    for (int i=0; i<size; i++) {

        tmp = (input[i] - min) * (HIST_SIZE / (max - min - 1));
        
        #pragma omp atomic
        histogram[tmp]++;
    }
}

/*
(c) Write another alternative implementation for the same function based on the use of OpenMP
locks, in which you maximize the parallelism in the update of the histogram.
*/
void FindFrequency(int * input, int size , int * histogram, int min, int max) {
    
    omp_lock_t locks[histogram.size()];

    for (i = 0; i < histogram.size(); i++)
        omp_init_lock(&locks[i]);
    
    int tmp;

    #pragma omp parallel for private(tmp)
    for (int i=0; i<size; i++) {

        tmp = (input[i] - min) * (HIST_SIZE / (max - min - 1));
        
        omp_set_lock(&histogram[tmp])
        histogram[tmp]++;
        omp_unset_lock(&histogram[tmp])
    }
}

/*
(d) Finally, write a parallel version for function FindFrequency based on a divide and conquer
recursive strategy, making use of the synchronization and cut-off strategies that you consider
the most appropriate.
*/
#define CUT_OFF = 4;
#define MIN_SIZE = 8;

void rec_FindFrequency(int * input, int size , int * histogram, int min, int max, int k, int depth) {

    if (size <= MIN_SIZE)
        FindFrequency(input, size , histogram, min, max);

    else if (depth < CUT_OFF) {
        int n = size/2;
        #pragma omp task
        rec_FindFrequency(input, n, histogram, min, max, k, depth + 1);
        #pragma omp task
        rec_FindFrequency(input, size - n, histogram, min, max, k + n, depth + 1);

    } else {
        int n = size/2;
        rec_FindFrequency(input, n, histogram, min, max, k, depth);
        rec_FindFrequency(input, size - n, histogram, min, max, k + n, depth);
    }

}

void FindFrequency(int * input, int size , int * histogram, int min, int max) {
    
    int tmp;

    for (int i=0; i<size; i++) {

        tmp = (input[i] - min) * (HIST_SIZE / (max - min - 1));
        
        #pragma omp atomic
        histogram[tmp]++;
    }
}

// Llamada original
rec_FindFrequency(input,size,histogram,min,max,0,0)





/*
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
*/

/*
13. Assume a hash table implemented as a vector of chained lists, such as shown in the next figure and
type definition:
*/

#define SIZE_TABLE 1048576

typedef struct {
    int data;
    element * next;
} element;

typedef struct {
    omp_lock_t global_lock;
    element * entrada[SIZE_TABLE];
} HashTable;

HashTable table;

/*
Inside each list, the elements are stored ordered by their data field value. Also assume the following
code snippet to insert elements of the ToInsert vector of size num elem into the mentioned hash
table:
*/

#define MAX_ELEM 1024
int main() {

    int ToInsert[MAX_ELEM], num_elem, index;

    // ...
    omp_init_lock(&table.global_lock); // bloquea una entrada de la tabla de hash (toda la tira de elementos que contiene)

    #pragma omp parallel for private(index) schedule(static)
    for (i = 0; i < num_elem; i++) {
        index = hash_function(ToInsert[i], SIZE_TABLE);
        omp_set_lock (&table.global_lock);
        insert_elem (ToInsert[i], index);
        omp_unset_lock (&table.global_lock);
    }

    omp_destroy_lock(&table.global_lock);
    // ...
}

/*

a) to fasil , dibujar con calama respetando los bloqueos (preguntar cuando entra el set si estaba bloqueado)
b) pones un lock para cada entrada de la tabla
c) ?
d) pones un lock para cada elemento de la tabla (de todas las listas)

*/