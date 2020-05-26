#define N 1000

// BLOCK GEOMETRIC DATA DECOMPOSITION BY ROWS

    // Parallel version using omp forR
    sum = 0;
    #pragma omp parallel
    #pragma omp for private(j) schedule(static) reduction(+:sum)

    for(i=0; i<N; i++)
        for (j=0; j<N; j++)
            sum += f(Matrix_in[i][j])

    // Parallel version without using omp for and P multiple of N
    sum = 0;
    #pragma omp parallel private(i,j) reduction(+:sum)
    {
        int my_id = omp_get_thread_num();
        int BS = N / omp_get_num_threads();
        int i_start = my_id * BS;
        int i_end = i_start + BS;

        for (i = i_start; i < i_end; i++)
            for (j = 0; j < N; j++)
                sum += f(Matrix_in[i][j])
    }

    // Parallel version without using omp for and P NOT-multiple of N
    sum = 0;
    #pragma omp parallel private(i,j) reduction(+:sum)
    {
        int my_id   = omp_get_thread_num();
        int BS      = N / omp_get_num_threads();
        int i_start = my_id * BS;
        int i_end   = i_start + BS;

        if (my_id < something)  i_start = i_start + my_id;
        else                    i_start = i_start + something;

        if (my_id <something)   i_end = i_end+1;

        for (i = i_start; i < i_end; i++)
            for (j = 0; j < N; j++)
                sum += f(Matrix_in[i][j])
    }


// CYCLIC GEOMETRIC DATA DECOMPOSITION BY ROWS

    // Parallel version using omp for
    sum = 0;
    #pragma omp parallel
    #pragma omp for private(j) schedule(static, 1) reduction(+: sum)
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            sum += f(Matrix_in[i][j])

    // Parallel version without omp for   
    sum = 0;
    #pragma omp parallel private(i, j) reduction(+: sum)
    {
        int my_id   = omp_get_thread_num();
        int nt      = omp_get_num_threads();
        int i_start = my_id;
        int i_end   = N;
        
        for (i = i_start; i < i_end; i += nt)
            for (j = 0; j < N; j++)
                sum += f(Matrix_in[i][j])
    }

// BLOCK CYCLIC by Rows

    // Parallel version using omp for
    sum = 0;
    #pragma omp parallel
    #pragma omp for private(j) schedule(static,2) reduction(+:sum)
    for (i=0; i<N; i++)
        for (j=0; j<N; j++)
            sum+= f(Matrix_in[i][j])
    
    // Parallel version without omp for
    sum = 0;
    #pragma omp parallel private(ii,i,j) reduction(+:sum)
    {
        int B       = 1 /* Numero de filas de cada bloque*/
        int nt      = omp_get_num_threads();
        int my_id   = omp_get_thread_num();
        int i_start = my_id*B;
        int i_end   = N;

        for(ii=i_start; ii<i_end; ii+=nt*B)
            for (i=ii; i<ii+B; i++)
                for (j=0; j<N; j++)
                    sum+= f(Matrix_in[i][j])
    }


// BLOCK DECOMPOSITION

    // Parallel version without omp for
    sum = 0;
    #pragma omp parallel private(ii, i, j) reduction(+: sum) num_threads(8)
    {
        int my_id;
        int NC      = 4;                 // Numero de columnas por bloque
        int NF      = 2;                 // Numero de filas por bloque
        int block_i = my_id / NC;        // 0 or 1 indicates which row of blocks
        int block_j = my_id % NC;        // 0,1,2,3 indicates which column of blocks
        
        int BSi     = N / NF;            // Assume N is multiple of 2
        int BSj     = N / NC;            // Assume N is multiple of 4
        
        int i_start = block_i * Bsi;
        int i_end   = i_start + BSi;
        int j_start = block_j * BSj;
        int j_end   = j_start + BSj;
        
        for (i = i_start; i < i_end; i++)
            for (j = j_start; j < j_end; j++)
                sum += f(Matrix_in[i][j])
    }