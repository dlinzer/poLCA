#include "R.h"

#define MAX_CLASSES 500  // Maximum number of latent classes

// function: hello
//    A little test
void hello(void) {
	Rprintf("Hello World\n");
}

// Function: ylik
//       A function to find the likelihood of an observation
//       given a vector of responses and a set of
//       choice probabilities for each latent class.
//
//       Args:
//            probs (double *):   A vector of choice probabilities
//                                over all items.
//            y (int *):          A vectorize matrix of choices
//            obs (int *):        Number of observations (rows of y)
//            items (int *):      Number of items (cols of y)
//            numChoices (int *): Vector given number of choices for
//                                each item.
//            classes (int):      Number of latent classes
//
//       Return:
//            llik (double *): likelihood vector for the observation
//
void ylik(double *probs, int *y, int *obs, int *items,
	  int *numChoices, int *classes, double *lik) {
	int i,j,k;
	const int citems = *items;
	const int cclasses = *classes;
	const int cobs = *obs;
	const double *firstprobs = probs;

	for (i=0;i<cobs;i++) {
		for (j=0;j<cclasses;j++) lik[j] = DBL_MAX;
		probs = (double *) firstprobs;
		for (k=0;k<citems;k++) {
			for (j=0;j<cclasses;j++) {
			   if (y[k]>0) lik[j] *= probs[y[k]-1];
			   probs += numChoices[k]; // move pointer to next item
			}
		}
		y += citems; // move pointer to next observation
		lik += cclasses; // move pointer to next observation
	}
}

// Function: postclass
//
//           A function to find the posterior distribution over the
//           latent classes for each observation.
//
//     Args:
//           prior (double *): Vector of prior class probs
//           probs (double *): An (N x J) matrix of choice probs
//           y (int *): An (N x J) matrix of choices
//           items (int *): Number of items
//           obs (int *): Number of observations
//           numChoices (int *): Vector given number of choices for
//                                each item.
//           classes (int):      Number of latent classes
//
//     Return:
//           post (double *): Vector of poster probs
//
void postclass(double *prior, double *probs, int *y,
               int *items, int *obs, int *numChoices,
               int *classes, double *posterior) {
	int i,j,totalChoices;
	double llik[MAX_CLASSES]; // Should probably calloc, limits num of classes0
	double denom;
	const int citems = *items;
	const int cobs = *obs;
	const int cclasses = *classes;
	int one = 1;

	totalChoices=0;
	for (i=0;i<citems;i++) totalChoices += numChoices[i];

	for (i=0;i<cobs;i++) {
		ylik(probs,y, (int *) &one,items,numChoices,classes,llik);
		denom = 0.0;
		for (j=0;j<cclasses;j++) denom+=prior[j]*llik[j];
		for (j=0;j<cclasses;j++) {
			posterior[j]=prior[j]*llik[j]/denom;
		}
		y+=citems; //Increment y pointer to next obs
		prior+=cclasses; // Increment proir pointer to next obs
		posterior+=cclasses;
	}
}

// Function: probhat
//     A function to return updates estimates the response probs
//     within each class given the data and posterior distributions
//     Just a bunch of conditional means.
//
//  Args:
//     y (int *): Matrix of response data
//     post (double *): Matrix of posterior estimates
//     items (int *): Number of items
//     obs (int *): Number of observations
//     numChoices (int *): Vector given number of choices for
//                         each item.
//     classes (int):      Number of latent classes
//
//   Returns:
//     probhat (double *): Estimates response probs in mat form
void probhat(int *y, double *post,
               int *items, int *obs, int *numChoices,
               int *classes, double *probhat) {
	double *denom;
	int i,j,k,cumChoices;
	const int citems = *items;
	const int cobs = *obs;
	const int cclasses = *classes;

	int totalChoices=0;
	for (i=0;i<citems;i++) totalChoices += numChoices[i];
	for (i=0;i<(totalChoices*cclasses);i++) probhat[i]=0.0;

 	denom = (double *) calloc((cclasses*citems),sizeof(double));
	for (i=0;i<(cclasses*citems);i++) denom[i]=0.0;

	for (i=0;i<cobs;i++) {
		for (j=0;j<cclasses;j++) {
			cumChoices=0;
			for (k=0;k<citems;k++) {
				if (y[k]>0) {
					probhat[j*numChoices[k]+cclasses*cumChoices+y[k]-1] += post[j];
					denom[j*citems+k] += post[j];
				}
				cumChoices += numChoices[k];
			}
		}
	    y+=citems;
	    post+=cclasses;
	}

	for (j=0;j<cclasses;j++) {
		cumChoices=0;
		for (k=0;k<citems;k++) {
			for (i=0;i<numChoices[k];i++) {
				probhat[j*numChoices[k]+cclasses*cumChoices+i] =
				             (double) probhat[j*numChoices[k]+cclasses*cumChoices+i]/
					         denom[j*citems+k];
			}
			cumChoices += numChoices[k];
		}
	}
	free(denom);
}



// Function: d2lldbeta2
//	Gives the first and second derivatives of the (cond) likelihood
//      w.r.t. the vector of betas.
//
//      Args:
//          rgivy (double *): Matrix of posteriors.
//			prior (double *): Matrix of priors.
//			x (double *): Matrix of indvars.
//          obs (int *): Number of observations.
//          classes (int *): Number of classes.
//          xcols (int *): Columns of x.
//      Return:
//          grad (double *): vector of first derivs.
//          hess (double *): Matrix of second derivs.
//
void d2lldbeta2(double *rgivy, double *prior, double *x, int *obs,
             int *classes, int *xcols, double *grad, double *hess) {
	int i,j,k,m,n,row,col,newrow,newcol;
	const int cobs = *obs;
	const int cclasses = *classes;
	const int cxcols = *xcols;
	const int crank =  cxcols*(cclasses-1);
	for (i=0;i<cobs;i++) {
		for (k=0;k<cxcols;k++) {
			//Gradient part
			for (j=1;j<cclasses;j++) {
				grad[(j-1)*cxcols+k] += (double) x[k]*(rgivy[j]-prior[j]);
			}
			//Hessian part
			for (m=0;m<=k;m++) {
				for (j=1;j<cclasses;j++) {
					col = (j-1)*cxcols + m;
					row = (j-1)*cxcols + k;
					//Diagonal block elements of the hessian
					hess[row*crank+col] += x[m]*x[k]*
					             ( -rgivy[j]*(1.0-rgivy[j]) +
					                prior[j]*(1.0-prior[j]) ) ;
					for (n=1;n<j;n++) {
						col = (n-1)*cxcols + m;
						//Subdiagonal elements of the hessian
						hess[row*crank+col] += x[m]*x[k]*
								  ( rgivy[j]*rgivy[n] - prior[j]*prior[n] );
					}
				}
			}
		}
		prior += cclasses;
		rgivy += cclasses;
		x += cxcols;
	}

	// Copy the upper elements of the symetric off-diag blocks.
	for (i=1;i<cclasses;i++) {
		for (j=i+1;j<cclasses;j++) {
			for (k=0;k<cxcols;k++) {
				for (n=k+1;n<cxcols;n++) {
					row = (j-1)*cxcols+n;
					col = (i-1)*cxcols+k;
					newrow = (j-1)*cxcols+k;
					newcol = (i-1)*cxcols+n;
					hess[newrow*crank+newcol]=hess[row*crank+col];
				}
			}
		}
	}

	// Copy lower to upper off diagional elements
	for (col=0;col<(cclasses-1)*cxcols;col++) {
		for (row=0;row<col;row++) {
			  hess[row*crank+col] = hess[col*crank+row];
		}
	}
}


