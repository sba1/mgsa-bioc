/**
 * @file mgsa.c
 *
 * @author Sebastian Bauer
 */
#include <stdio.h>

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


/**
 * The work-horse.
 *
 * @param sets pointer to the sets
 * @param lset specifies the length of each set
 * @param lsets number of sets (length of sets)
 * @param n the number of observable entities.
 * @param o indices of entities which are "on" (0 based).
 * @param lo length of o
 */
static void do_mgsa_mcmc(int **sets, int *lset, int lsets, int n, int *o, int lo)
{
	int i;

	printf("n=%d\n",n);

	for (i=0;i<lo;i++)
	{
		printf(" o[%d]=%d\n",i,o[i]);
	}
}

/**
 * @param n the number of observable entities.
 * @param o (1 based)
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

	/* Set associations  */
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
		if (!(nsets[i] = (int*)R_alloc(lset[i],sizeof(nsets[i][0]))))
		{
			UNPROTECT(1);
			goto bailout;
		}

		xset = INTEGER_POINTER(o);
		for (j=0;j<lset[i];j++)
			nsets[i][j] = lset[i];
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
