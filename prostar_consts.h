#ifndef PROSTAR_CONSTS_H
#define PROSTAR_CONSTS_H

/* Constants from Table 1 of Goñi et al. Genome Biology 2007 8:R263 doi:10.1186/gb-2007-8-12-r263  */
/* directly http://genomebiology.com/2007/8/12/R263/table/T1 */

/* The order in the vector is:  */
/*  (0) Twist                   */
/*  (1) Tilt                    */
/*  (2) Roll                    */
/*  (3) Shift                   */
/*  (4) Slide                   */
/*  (5) Rise                    */
                             
static const double prostar_mat[60] ={0.026, 0.038, 0.020, 1.69, 2.26, 7.65,
			              0.036, 0.038, 0.023, 1.32, 3.03, 8.93,
			              0.031, 0.037, 0.019, 1.46, 2.03, 7.08,
			              0.033, 0.036, 0.022, 1.03, 3.83, 9.07,
			              0.016, 0.025, 0.017, 1.07, 1.78, 6.38,
			              0.026, 0.042, 0.019, 1.43, 1.65, 8.04,
			              0.014, 0.026, 0.016, 1.08, 2.00, 6.23,
			              0.025, 0.038, 0.020, 1.32, 1.93, 8.56,
			              0.025, 0.036, 0.026, 1.20, 2.61, 9.53,
			              0.017, 0.018, 0.016, 0.72, 1.20, 6.23};

#endif
