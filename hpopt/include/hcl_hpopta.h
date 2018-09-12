/*--------------------------------------------------------*/

/*
@file
@author Hamidreza Khaleghzadeh <hamidreza.khaleghzadeh@ucdconnect.ie>
@version 1.0
*/

/*--------------------------------------------------------*/

/**

* NOTE: Use the library hpopt for compilation.

**/

/*--------------------------------------------------------*/

#ifndef _HCL_HPOPTA_H_
#define _HCL_HPOPTA_H_

/*--------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

/*--------------------------------------------------------*/

/**
 * Implementation of the HPOPTA algorithm for heterogeneous
 * clusters
 * The goal of the algorithm is to determine optimal workload
 * distribution minimizing the execution time 
 * of the computations in its parallel execution.
 *
 * Refer to:
 * Khaleghzadeh, Reddy and Lastovetsky, 
 * A Novel Data-Partitioning Algorithm for Performance Optimization 
 * of Data-Parallel Applications on Heterogeneous HPC Platforms.
 *
 * @param verbosity Level of verbosity.
 * @param n The workload size.
 * @param p Number of processors.
 * @param npoints The number of points in each discrete time function in
 *								an array of size p
 * @param psizes A list of problem sizes in p discrete time functions. Values
 								 for the first function are followed by thoes of the
 									second one and so on.
 *               
 * @param etimes A list of execution times in p discrete time functions. Values
 								 for the first function are followed by thoes of the
 									second one and so on.
 *               
 * @param dOpt The optimal workload distribution output in an array of
 *             size npoints.
 * @param eOpt The optimal execution time.
 *
 * @return 0 on successful invocation, otherwise -1.
 */
extern int
hcl_hpopta(
    const int verbosity,
    const unsigned int n,
    const unsigned int p,
    const unsigned int* npoints,
    const unsigned int* psizes,
    const double* etimes,
    unsigned int* dOpt,
    double* eOpt
);

/*--------------------------------------------------------*/

#ifdef __cplusplus
}
#endif

/*--------------------------------------------------------*/

#endif /* _HCL_HPOPTA_H_ */

/*--------------------------------------------------------*/
