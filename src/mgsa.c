/**
 * @file mgsa.c
 *
 * @author Sebastian Bauer
 *
 * @note if you want to compile as standalone, invoke
 *
 *  "gcc mgsa.c -DSTANDALONE `R CMD config --cppflags` `R CMD config --ldflags` -o mgsa"
 */
#include <stdio.h>

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* Enable debugging */
#define DEBUG

/* Define if file should be compiled as standalone (mainly for testing purposes) */
/* #define STANDALONE */

struct context
{
	/** @brief Number of sets */
	int number_of_sets;

	/** @brief Array of size "number_of_sets" indicating the state of the sets (active or not active) */
	int *sets_active;

	/** @brief Number of observable */
	int number_of_observable;

	/** @brief Array indicating whether a observable is active or not */
	int *observable;
};

/**
 * The work horse.
 *
 * @param sets pointer to the sets. Sets a made of observable entities.
 * @param sizes_of_sets specifies the length of each set
 * @param number_of_sets number of sets (length of sets)
 * @param n the number of observable entities.
 * @param o indices of entities which are "on" (0 based).
 * @param lo length of o
 */
static void do_mgsa_mcmc(int **sets, int *sizes_of_sets, int number_of_sets, int n, int *o, int lo)
{
	int i,j;
	struct context context;

#ifdef DEBUG

	for (i=0;i<number_of_sets;i++)
	{
		printf(" %d: ", sizes_of_sets[i]);
		for (j=0;j<sizes_of_sets[i];j++)
			printf("%d ", sets[i][j]);
		printf("\n");
	}

	printf(" n=%d\n",n);

	for (i=0;i<lo;i++)
	{
		printf(" o[%d]=%d\n",i,o[i]);
	}
#endif

	context.number_of_sets = number_of_sets;
	context.number_of_observable = n;
	if (!(context.observable = (int*)R_alloc(n,sizeof(context.observable[0]))))
		goto bailout;
	memset(context.observable,0,n * sizeof(context.observable[0]));
	for (i=0;i<lo;i++)
		context.observable[o[i]] = 1;
	if (!(context.sets_active = (int*)R_alloc(number_of_sets,sizeof(context.sets_active[0]))))
		goto bailout;
	memset(context.sets_active,0,number_of_sets * sizeof(context.sets_active[0]));


bailout:
	return;
}

/**
 * @param n the number of observable entities.
 * @param o (1 based)
 *
 * @note TODO: Check whether sets are real sets.
 */
SEXP mgsa_mcmc(SEXP sets, SEXP n, SEXP o, SEXP alpha, SEXP beta, SEXP p)
{
	int *xo,*no,lo;
	int **nsets, *lset, lsets;
	int i,j;

	if (LENGTH(n) != 1)
		error("Parameter 'n' needs to be atomic!");

	PROTECT(n = AS_INTEGER(n));
	PROTECT(o = AS_INTEGER(o));
	PROTECT(sets = AS_LIST(sets));

	/* Observations */
	xo = INTEGER_POINTER(o);
	lo = LENGTH(o);
	if (!(no =  (int*)R_alloc(lo,sizeof(no[0])))) /* R's garbage collection takes care of it */
		goto bailout;

	/* Turn 1 based observation indices into 0 based ones */
	for (i=0;i<lo;i++)
		no[i] = xo[i] - 1;

	/* Set associations. Turn 1 based observations indices into 0 based ones */
	lsets = LENGTH(sets);
	if (!(lset = (int*)R_alloc(lsets,sizeof(lset[0]))))
		goto bailout;
	if (!(nsets = (int**)R_alloc(lsets,sizeof(nsets[0]))))
		goto bailout;
	for (i=0;i<lsets;i++)
	{
		int *xset;

		SEXP set = VECTOR_ELT(sets, i);

		PROTECT(set = AS_INTEGER(set));
		lset[i] = LENGTH(set);
		xset = INTEGER_POINTER(set);
		if (!(nsets[i] = (int*)R_alloc(lset[i],sizeof(nsets[i][0]))))
		{
			UNPROTECT(1);
			goto bailout;
		}

		xset = INTEGER_POINTER(set);
		for (j=0;j<lset[i];j++)
			nsets[i][j] = xset[j];
		UNPROTECT(1);
	}

	do_mgsa_mcmc(nsets, lset, lsets, INTEGER_VALUE(n),no,lo);

bailout:
	UNPROTECT(3);
	return NULL_USER_OBJECT;
}

R_CallMethodDef callMethods[] = {
   {"mgsa_mcmc", (DL_FUNC)&mgsa_mcmc, 6},
   {NULL, NULL, 0}
};

/**
 * Called by R as the initialization routine.
 *
 * @param info
 */
void R_init_myLib(DllInfo *info)
{
   R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

/**
 * Called by R when the object is unloaded.
 *
 * @param info
 */
void R_unload_mylib(DllInfo *info)
{
  /* Release resources. */
}

#ifdef STANDALONE

int main(void)
{
	int t1[] = {0,1};
	int t2[] = {1,2};
	int o[] = {0,1};

	int *sets[] = {t1,t2};
	int sizes_of_sets[] = {sizeof(t1)/sizeof(t1[0]),sizeof(t2)/sizeof(t2[0])};

	do_mgsa_mcmc(sets, sizes_of_sets, sizeof(sets)/sizeof(sets[0]), 3, o, sizeof(o)/sizeof(o[0]));

}

#endif
