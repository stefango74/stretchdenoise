/*
 * ConstrainedLS.cpp
 *
 *  Created on: Mar 26, 2018
 *      Author: stef
 */

#include "ConstrainedLS.h"

#include <map>
#include <iostream>
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

/*
 * solve for x which minimizes ||y-Hx||_2^2, constraint to Cx=b
 */
void minimizeLSWithLinearConstraints(MatrixXd &H, VectorXd &y, MatrixXd &C, VectorXd &b, VectorXd &x)
{
	// transform using Lagrange multipliers, see page 4 in http://eeweb.poly.edu/iselesni/lecture_notes/least_squares/least_squares_SP.pdf
	// x = (H^T*H)^-1(H^T*y-C^T*(C*(H^T*H)^-1*C^T)^-1*(C*(H^T*H)^-1*H^T*y-b))
	// x = H2I*(H^T*y-C^T*(C*H2I*C^T)^-1*(C*H2I*H^T*y-b)), H2I = (H^T*H)^-1, H2IC = C*H2I, Hty = H^T*y, Ct = C^T
	MatrixXd Ht = H.transpose();
	MatrixXd H2I = (Ht*H).inverse();
	MatrixXd Hty = Ht*y;
	MatrixXd Ct = C.transpose();
	x = H2I*(Hty - Ct*((C*H2I*Ct).inverse())*(C*H2I*Hty - b));
}

/*
 * solve system with Lagrange multipliers and bounds iteratively, by making out-of-bound x_i constant
 */
ConstrainedLS::ConstrainedLS(int rows, int cols, int rowsC, double *_H, double *_y, double *_C, double *_b, double *_x, double *_lb, double *_ub, double *_w)
{
	int i, j, k;
	VectorXd x(rows);
	MatrixXd H(rows, cols);
	VectorXd y(rows);
	MatrixXd C(rowsC, cols);
	VectorXd b(rowsC);

	for (j = 0; j < rows; j++)
	{
		for (i = 0; i < cols; i++)
			H(j, i) = _H[j*cols + i];

		y(j) = _y[j];
	}

	for (j = 0; j < rowsC; j++)
	{
		for (i = 0; i < cols; i++)
			C(j, i) = _C[j*cols + i];

		b(j) = _b[j];
	}

	bool boundsOK;
	bool boundOK[cols];
	int colsNew = cols;
	map<int, int> colMap;	// map from current (shrunken) column indices to original

	for (j = 0; j < cols; j++)
	{
		boundOK[j] = true;
		colMap[j] = j;
	}

	do
	{
		boundsOK = true;
		minimizeLSWithLinearConstraints(H, y, C, b, x);

		for (j = 0; j < colsNew; j++)
			if ((x(j) < _lb[colMap[j]]) || (x(j) > _ub[colMap[j]]))
			{
				boundsOK = false;
				boundOK[colMap[j]] = false;

				// clamp and assign out-of-bound values
				if (x(j) < _lb[colMap[j]])
					x(j) = _lb[colMap[j]];

				if (x(j) > _ub[colMap[j]])
					x(j) = _ub[colMap[j]];

				_x[colMap[j]] = x(j);
			}

		if (!boundsOK)
		{
			// move out-of-bound x_i to constant
			colsNew = 0;

			for (j = 0; j < cols; j++)
				if (boundOK[j])
					colMap[colsNew++] = j;

			H = MatrixXd(rows, colsNew);
			C = MatrixXd(rowsC, colsNew);

			for (j = 0; j < rows; j++)
			{
				y(j) = _y[j];
				k = 0;

				for (i = 0; i < cols; i++)
				{
					if (boundOK[i])
						H(j, k++) = _H[j*cols + i];	// transfer matrix value
					else
						y(j) -= _x[i]*_H[j*cols + i];	// substract const value from right side
				}
			}

			for (j = 0; j < rowsC; j++)
			{
				b(j) = _b[j];
				k = 0;

				for (i = 0; i < cols; i++)
				{
					if (boundOK[i])
						C(j, k++) = _C[j*cols + i];	// transfer matrix value
					else
						b(j) -= _x[i]*_C[j*cols + i];	// substract const value from right side
				}
			}
		}

	} while (!boundsOK && (colsNew > 0));

	// assign remaining x values
	for (j = 0; j < colsNew; j++)
		_x[colMap[j]] = x(j);
}
