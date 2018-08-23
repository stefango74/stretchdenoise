/*
 * ConstrainedLS.h
 *
 *  Created on: Apr 3, 2018
 *      Author: stef
 */

#ifndef SRC_CONSTRAINEDLS_H_
#define SRC_CONSTRAINEDLS_H_

class ConstrainedLS
{
public:
	ConstrainedLS(int rows, int cols, int rowsC, double *_H, double *_y, double *_C, double *_b, double *_x);
	ConstrainedLS(int rows, int cols, int rowsC, double *_H, double *_y, double *_C, double *_b, double *_x, double *_lb, double *_ub, double *_w);
};

#endif /* SRC_CONSTRAINEDLS_H_ */
