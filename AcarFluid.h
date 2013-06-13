/*
Copyright (c) 2013- Nervous System, inc

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "fftw3.h"
#include "ofVec2f.h"
#include <math.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <algorithm>

typedef double Real;
typedef ofVec2f vec2;

/*
	Fluid solver based on the spectral method in Level Set Driven Flows by Acar 2007

*/

class acarFluid {
public:
	acarFluid();
	void setup(int w, int h);
	void step();
	void allocate();
	void free();
	void reset();
	void smoothLevelSet(Real *in, Real *out);
	void solvePoisson(Real *in, Real *out);
	void init();
	void stroke(vector<vec2> &pts);
	//this really ought to be a single output grid
	void computeVelocity(Real *in, Real *outU, Real *outV);
	void advectLinear(Real * in, Real *out);
	const int getWidth() {return width;};
	const int getHeight() {return height;};
	Real getLevelSet(int i,int j) {return levelsetSmooth[i*height+j];};


	Real * levelset;
	
private:
	int width;
	int height;
	Real dt;
	Real dx, dxSq;
	//vorticity
	Real * vorticity;
	//velocity is the curl of vorticity
	Real * velocityU;
	Real * velocityV;
	//level set input, signed distance function, replaces vorticity I believe
	//smoothed by some function
	Real * levelsetSmooth;
	//fourier transform of the smoothed levelset
	Real * levelsetF;
	Real * levelsetTemp;
	//stream function, used to compute velocity
	Real * stream;
	//fourier transform of stream function
	Real * streamF;
	fftw_plan levelsetPlan;
	fftw_plan streamPlanI;
	//used in level set smoothing function
	Real sheetWidth;
	Real smoothingParam1;
	Real smoothingParam2;
};

