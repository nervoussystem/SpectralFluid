/*
Copyright (c) 2013- Nervous System, inc

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "AcarFluid.h"

inline int FloatToInt( float x )
{
	return _mm_cvtt_ss2si( _mm_load_ss( &x ) );
}

inline int FloatToInt( double x ) 
{
	return _mm_cvttsd_si32( _mm_load_sd( &x) );
}
acarFluid::acarFluid() {
	sheetWidth = 10.0;//0.1;//No idea what this should be
	smoothingParam1 = 3.14159265359/(2.0*sheetWidth*sheetWidth);
	smoothingParam2 = 3.14159265359/sheetWidth;
	dt = .01;
	dx = .1;
	dxSq = dx*dx;
}

void acarFluid::setup(int w, int h) {
	width = w;
	height = h;
	allocate();

	//setup the "plan"
	//consider using FFTW_DESTROY_INPUT for more efficiency, also consider FFTW_PATIENT 
	levelsetPlan = fftw_plan_r2r_2d(width,height,levelsetSmooth, levelsetF, FFTW_RODFT10, FFTW_RODFT10, 0);
	streamPlanI = fftw_plan_r2r_2d(width,height, streamF,stream, FFTW_RODFT01, FFTW_RODFT01, 0);
	reset();
}

void acarFluid::step() {
	//smooth level set?
	smoothLevelSet(levelset, levelsetSmooth);
	//solve for stream function using the level set
	fftw_execute( levelsetPlan);

	solvePoisson(levelsetF, streamF);
	fftw_execute(streamPlanI);
	//use stream function to compute velocity
	computeVelocity(stream, velocityU, velocityV);
	//use velocity to advect (level set?, smoothed level set?, vorticity?)
	advectLinear(levelset, levelsetTemp);
	for(int i=0;i<width*height;++i) levelset[i] = levelsetTemp[i];
	//repeat
}

/*
	smoothing function for level set input
	this could be changed so that it uses a faster function than sin
*/
void acarFluid::smoothLevelSet(Real *in, Real *out) {
	for(int i=0;i<width*height;++i) {
		if(abs(in[i]) < sheetWidth) {
			out[i] = smoothingParam1*sin(smoothingParam2*in[i]);
		} else {
			out[i] = 0;
		}
	}
}

/*
	solve the poisson equation in fourier space
	not sure how the "wave mode" variables should scale
*/
void acarFluid::solvePoisson(Real *in, Real *out) {
	int index = 0;
	Real wave1 = PI/width*dx;
	Real wave2 = PI/height*dx;
	wave1 *= wave1;
	wave2 *= wave2;
	for(int i=0;i<width;++i) {
		for(int j=0;j<height;++j) {
			index++;
			if(i==0 || j==0) {
				out[index] = 0;
			} else {
				out[index] = -in[index]/(wave1*i*i+wave2*j*j);
			}
		}
	}
}

/*
	compute the velocity field as the curl of the stream field
	u = dStream/dy
	v = dStream/dx
	should we introduce a length scale?
	we may need to divide by 4*width*height to normalize the effects of the FFT
*/
void acarFluid::computeVelocity(Real *in, Real * outU, Real * outV) {
	int starti = 1;
	int startj = 1;
	int endi = width-1;
	int endj = height-1;

	int index;
	Real normalizationConstant = 4.0*width*height;
	for(int i=starti;i<endi;++i) {
		for(int j=startj;j<endj;++j) {
			index = i*height+j;
			outU[index] = -0.5/dx*(in[index+1]-in[index-1])/normalizationConstant;
			outV[index] = 0.5/dx*(in[index+height]-in[index-height])/normalizationConstant;
		}
	}
}

void acarFluid::allocate() {
	vorticity = fftw_alloc_real(width*height);
	velocityU = fftw_alloc_real(width*height);
	velocityV = fftw_alloc_real(width*height);
	levelset = fftw_alloc_real(width*height);
	levelsetTemp = fftw_alloc_real(width*height);
	levelsetF = fftw_alloc_real(width*height);
	levelsetSmooth = fftw_alloc_real(width*height);
	stream = fftw_alloc_real(width*height);
	streamF = fftw_alloc_real(width*height);
}

void acarFluid::free() {
	fftw_free(vorticity);
	fftw_free(velocityU);
	fftw_free(velocityV);
	fftw_free(levelset);
	fftw_free(levelsetTemp);
	fftw_free(levelsetSmooth);
	fftw_free(levelsetF);
	fftw_free(stream);
	fftw_free(streamF);
}

void acarFluid::reset() {
	for(int i=0;i<width*height;++i) {
		levelset[i] = 0;
		velocityU[i] = 0;
		velocityV[i] = 0;
	}
}

/*
advection function from CinderFX
void Advect2D
( 
	RealT						aDissipation, 
	RealT						aDt, 
	const Grid2D<T>&			aSrc, 
	const Grid2D<Vec2<RealT> >&	aVel, 
	Grid2D<T>&					aDst,
	int							aBorder = 1
)
{
	// Range
	int iStart = aBorder;
	int iEnd   = aSrc.resX() - aBorder;
	int jStart = aBorder;
	int jEnd   = aSrc.resY() - aBorder;

	// Boundary
	const RealT xMin = (RealT)0.5;
	const RealT xMax = (RealT)aSrc.resX() - (RealT)1.5;
	const RealT yMin = (RealT)0.5;
	const RealT yMax = (RealT)aSrc.resY() - (RealT)1.5;

	// Process
	for( int j = jStart; j < jEnd; ++j ) {
		for( int i = iStart; i < iEnd; ++i ) {
			// Velocity
			const Vec2<RealT>& vel = aVel.at( i, j );

			// Previous
			RealT dx = aDt*vel.x;
			RealT dy = aDt*vel.y;
			RealT iPrev = i - dx;
			RealT jPrev = j - dy;
			iPrev = Clamp( iPrev, xMin, xMax );
			jPrev = Clamp( jPrev, yMin, yMax );

			// Advected value
			T advected = aSrc.bilinearSample( iPrev, jPrev );

			// Update
			aDst.at( i, j ) = aDissipation*advected;
		}
	}
}
*/

/*

template <typename RealT>
	DataT bilinearSample( RealT aX, RealT aY ) const {
		int x0 = FloatToInt( aX );
		int y0 = FloatToInt( aY );
		int x1 = x0 + 1;
		int y1 = y0 + 1;
		RealT a1 = aX - (RealT)x0;
		RealT b1 = aY - (RealT)y0;
		RealT a0 = (RealT)1 - a1;
		RealT b0 = (RealT)1 - b1;
		return b0*( a0*at( x0, y0 ) + a1*at( x1, y0 ) ) + 
			   b1*( a0*at( x0, y1 ) + a1*at( x1, y1 ) );
	}

	*/

void acarFluid::advectLinear(Real * in, Real *out) {
	int starti = 1;
	int startj = 1;
	int endi = width-1;
	int endj = height-1;

	Real maxX = width-1.5;
	Real maxY = height-1.5;
	int index;
	Real vx,vy;
	for(int i=starti;i<endi;++i) {
		for(int j=startj;j<endj;++j) {
			index = i*height+j;
			vx = velocityU[index];
			vy = velocityV[index];
			// Previous
			Real dx = dt*vx;
			Real dy = dt*vy;
			Real iPrev = i - dx;
			Real jPrev = j - dy;
			iPrev = std::min(std::max(iPrev, 0.5), maxX );
			jPrev = std::min(std::max( jPrev, 0.5), maxY );

			// Advected value
			int x0 = FloatToInt( iPrev );
			int y0 = FloatToInt( jPrev );
			int index2 = x0*height+y0;
			Real a1 = iPrev - (Real)x0;
			Real b1 = jPrev - (Real)y0;
			Real a0 = 1.0 - a1;
			Real b0 = 1.0 - b1;
			out[index] = b0*( a0*in[index2] + a1*in[index2+height] ) + 
					b1*( a0*in[index2+1] + a1*in[index2+height+1] );
		}
	}
}

void acarFluid::init() {
	int index = 0;
	for(Real i=1;i<width-1;++i) {
		for(Real j=1;j<height-1;++j) {
			index = i*height+j;
			levelset[index] = i-120.0+5.0*sin(20.0*PI*j/height);
		}
	}
}

void acarFluid::stroke(vector<vec2> &pts) {

}