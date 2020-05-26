
// Hash table BLOCK GEOMETRIC input data decomposition

#pragma omp parallel private(index, i)
{
	int myid = omp_get_thread_num();
	int numprocs = omp_get_num_threads();
	int lower = myid * (MAX_ELEM/numprocs);
	int upper = lower + (MAX_ELEM/numprocs);
	if (myid == numprocs-1) upper = MAX_ELEM;
	for (i = lower; i < upper; i++) {
		index = hash_function(ToInsert[i], SIZE_TABLE);
		omp_set_lock(&locks[index]);
		insert_elem (ToInsert[i], index);
		omp_unset_lock(&locks[index]);
	}
}

// Hash table BLOCK GEOMETRIC output data decomposition

#pragma omp parallel private(index, i)
{
	int myid = omp_get_thread_num();
	int numprocs = omp_get_num_threads();
	int lower = myid * (SIZE_TABLE/numprocs);
	int upper = lower + (SIZE_TABLE/numprocs);
	if (myid == numprocs-1) upper = SIZE_TABLE;
	for (i = 0; i < MAX_ELEM; i++) {
		index = hash_function(element[i], SIZE_TABLE);
		if ((index>=lower) && (index<upper)) insert_elem (element[i], index);
	}
}

// Hash table CYCLIC output
#pragma omp parallel private(index, i)
{
	// each thread “self-assigns” the range of elements to process
	int myid = omp_get_thread_num();
	int nt = omp_get_num_threads();
	// index%nt determines which is the thread that should
	// address this output
	// Assume nt=3
	// 0 1 2 3 4 5 6 7 8 9 10 - index
	// 0 1 2 0 1 2 0 1 2 0 1 - index%nt → my_id
	for (i = 0; i < MAX_ELEM; i++) {
		index = hash_function(ToInsert[i], SIZE_TABLE);
		if ((index%nt)==my_id) insert_elem (ToInsert[i], index);
	}
}

// Hash table BLOCK-CYCLIC output
#pragma omp parallel private(index, i)
{
	// each thread “self-assigns” the range of elements to process
	int myid = omp_get_thread_num();
	int nt = omp_get_num_threads();
	// Assume nt=3, block=2 -> blk
	// 0 1 2 3 4 5 6 7 8 9 10 11- index
	// 0 0 1 1 2 2 3 3 4 4 5 5 – (index/blk)
	// 0 0 1 1 2 2 0 0 1 1 2 2 – ((index/blk)%nt) → my_id
	for (i = 0; i < MAX_ELEM; i++) {
		index = hash_function(ToInsert[i], SIZE_TABLE);
		if ((index/blk)%nt==my_id) insert_elem (ToInsert[i], index);
	}
}


