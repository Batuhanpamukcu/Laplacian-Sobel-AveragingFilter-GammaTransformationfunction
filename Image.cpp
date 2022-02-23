/**********************************************************
 The original copy of the code can be found at http://web.eecs.utk.edu/~mkarakay/courses_files/testfiles.zip 
 and it is modified for ELM463/667
 
 * Image.cpp - the image library which implements
 *             the member functions defined in Image.h
 *
 * Author: Hairong Qi, ECE, University of Tennessee
 *
 * Created: 02/05/02
 *
 * Copyright (C) hqi@utk.edu
 **********************************************************/

#include "Image.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

/**
 * Default constructor.
 */ 
Image::Image() {
  image = NULL;
  nrows = 0;
  ncols = 0;
  maximum = 255;
  //createImage(0, 0);
}

/**
 * Constructor for grayscale images.
 * @param nrows Numbers of rows (height).
 * @param ncols Number of columns (width).
 * @return The created image.
 */
Image::Image(int nRows, int nCols) {
  if (nRows<=0 || nCols<=0) {
    cout << "Image: Index out of range.\n";
    exit(3);
  }
  image = NULL;
  createImage(nRows, nCols);
}

/**
 * Copy constructor. 
 * @param img Copy image.
 * @return The created image.
 */
Image::Image(const Image &img) {
  int rows, cols;

  image = NULL;
  nrows = img.getRow();
  ncols = img.getCol();
  createImage(nrows, ncols);             // allocate memory
  
  for (rows=0; rows < nrows; rows++)
    for (cols=0; cols < ncols; cols++)
		image[rows * ncols + cols] = img(rows, cols);
}

/**
 * Destructor.  Frees memory.
 */
Image::~Image() {
  if (image)
    delete [] image;       // free the image buffer
}






/**
 * Allocate memory for the image and initialize the content to be 0.
 */
void Image::createImage() {

  if (image != NULL)
    delete [] image;

  maximum = 255;

  image = (float *) new float [nrows * ncols];
  if (!image) {
    cout << "CREATEIMAGE: Out of memory.\n";
    exit(1);
  }

  initImage();
}


/**
 * Allocate memory for the image and initialize the content to be zero.
 * @param r Numbers of rows (height).
 * @param c Number of columns (width).
 */
void Image::createImage(int numberOfRows, int numberOfColumns) {
  
  if (image != NULL)
    delete [] image;

  nrows = numberOfRows;
  ncols = numberOfColumns;
  maximum = 255;

  image = (float *) new float [nrows * ncols];
  if (!image) {
    cout << "CREATEIMAGE: Out of memory.\n";
    exit(1);
  }

  initImage();
}

/**
 * Initialize the image.
 * @para init The value the image is initialized to. Default is 0.0.
 */
 
 
void Image::initImage(float initialValue) {
  int i;

  for (i = 0; i < nrows * ncols; i++){
    image[i] = initialValue;
    initialValue++;
    }
}

/*
 * Create new image
 */

// increase pixels
void Image::New(float initialValue) {
  int i;

  for (i = 0; i < nrows * ncols; i++){
    image[i] = initialValue;
    initialValue++;
    if (ncols-1 == i)
    initialValue = initialValue - ncols  + 2;
    else if (2*ncols-1 == i)
    initialValue = initialValue - ncols  + 2;
    else if (3*ncols-1 == i)
    initialValue = initialValue - ncols  + 2;
    else if (4*ncols-1 == i)
    initialValue = initialValue - ncols  + 2;
    else if (5*ncols-1 == i)
    initialValue = initialValue - ncols  + 2;
    else if (6*ncols-1 == i)
    initialValue = initialValue - ncols  + 2;
    else if (7*ncols-1 == i)
    initialValue = initialValue - ncols  + 2;
    else if (8*ncols-1 == i)
    initialValue = initialValue - ncols  + 2;
    else if (9*ncols-1 == i)
    initialValue = initialValue - ncols  + 2;
    }
}

/**
 * Returns the total number of rows in the image.
 * @return Total number of rows.
 * \ingroup getset
 */
int Image::getRow() const {
  return nrows;
}

/**
 * Returns the total number of columns in the image.
 * @return Total number of columns.
 * \ingroup getset
 */
int Image::getCol() const {
  return ncols;
}

/**
 * Returns the maximum pixel value of a gray-level image. 
 * @return The intensity of that pixel.
 * \ingroup getset
 */
float Image::getMaximum() const {
  int i, j;
  float maxi=-10000;

 
  for (i=0; i<nrows; i++)
    for (j=0; j<ncols; j++)
      if (maxi < image[i*ncols+j])
	maxi = image[i*ncols+j];
  
  return maxi;
}


/**
 * Returns the minimum pixel value of the image.
 * @return The minimum pixel value.
 * \ingroup getset
 */
float Image::getMinimum() const {
  int i, j;
  float mini=10000; 

  for (i=0; i<nrows; i++)
    for (j=0; j<ncols; j++)
      if (mini > image[i*ncols+j])
	mini = image[i*ncols+j];

  return mini;
}



/**
 * Returns the pixel value at rows, cols
 * @return The pixel value
 * \ingroup getset
 */
float Image::getPix(int rows, int cols) {
  return image[rows * ncols + cols];
}


/**
 * Returns the image. 
 * @return a gray-scale image
 * \ingroup getset
 */
Image Image::getImage() const {
  Image temp;
  int rows, cols;
  
  temp.createImage(nrows, ncols);   // temp is a gray-scale image
  for (rows = 0; rows < nrows; rows++)
    for (cols = 0; cols < ncols; cols++)
      temp(rows, cols) = image[rows * ncols + cols];
      
  return temp;
}

/**
 * Sets the total number of rows in an image.
 * @param r Total number of rows.
 * \ingroup getset
 */
void Image::setRow(int numberOfRows) {
  nrows = numberOfRows;
}

/**
 * Sets the total number of columns in an image.
 * @param c Total number of columns.
 * \ingroup getset
 */
void Image::setCol(int numberOfColumns) {
  ncols = numberOfColumns;
}


/**
 * Sets the pixel value at rows,cols.
 * @param row and col index.
 * \ingroup getset
 */
void Image::setPix(int rows, int cols, float value) {
  image[rows * ncols + cols] = value;
}


/**
 * Sets the image given a grayscale image. 
 * \ingroup getset
 */
void Image::setImage(Image &img) {
  int rows, cols;

  for (rows = 0; rows < nrows; rows++)
    for (cols = 0; cols < ncols; cols++)
      image[rows * ncols + cols] = img(rows, cols);
}


/**
 * Overloading () operator
 * \ingroup overload
 * @param i Row
 * @param j Column
 */
float & Image::operator()(int rows, int cols) const {
  return image[rows * ncols + cols];
}

/**
 * Overloading = operator.
 * \ingroup overload
 * @param img Image to copy.
 * @return Newly copied image.
 */
const Image Image::operator=(const Image& img) {
  int rows, cols;

  if (this == &img)
    return *this;

  nrows = img.getRow();
  ncols = img.getCol();
  createImage(nrows, ncols);             

  for (rows = 0; rows < nrows; rows++)
    for (cols = 0; cols < ncols; cols++)
	(*this)(rows, cols) = img(rows, cols);

  return *this;
}

/**
 * Overloading + operator.
 * \ingroup overload
 * @param img Image to add to specified image.
 * @return Addition of the two images.
 */
Image Image::operator+(const Image& img) const {
  int i, j, nr, nc;
  Image temp;

  nr = img.getRow();
  nc = img.getCol();

  if (nr != nrows || nc != ncols) {
    cout << "operator+: "
         << "Images are not of the same size or type, can't do addition\n";
    exit(3);
  }
  temp.createImage(nrows, ncols);
  
  for (i=0; i<nrows; i++)
    for (j=0; j<ncols; j++)
        temp(i,j) = image[i*ncols+j] + img(i,j);

  return temp;
}

/**
 * Overloading - operator.
 * \ingroup overload
 * @param img Image to subtract from specified image.
 * @return Subtraction of the two images.
 */
Image Image::operator-(const Image &img) const {
   int i, j, nr, nc;
  Image temp;

  nr = img.getRow();
  nc = img.getCol();

  if (nr != nrows || nc != ncols) {
    cout << "operator-: "
         << "Images are not of the same size or type, can't do subtraction\n";
    exit(3);
  }
  temp.createImage(nrows, ncols);             
  
  for (i=0; i<nrows; i++)
    for (j=0; j<ncols; j++)
        temp(i,j) = image[i*ncols+j] - img(i,j);

  return temp;
}

/**
 * Overloading * operator.  This function does pixel by pixel multiplication.
 * \ingroup overload
 * @param img Image to multiply with specified image.
 * @return Multiplication of the two images.
 */
Image Image::operator*(const Image &img) const {
  int i, j, nr, nc;
  Image temp;

  nr = img.getRow();
  nc = img.getCol();

  if (nr != nrows || nc != ncols) {
    cout << "operator*: "
         << "Images are not of the same size or type, can't do multiplication\n";
    exit(3);
  }
  temp.createImage(nrows, ncols);             
  
  for (i=0; i<nrows; i++)
    for (j=0; j<ncols; j++)
        temp(i,j) = image[i*ncols+j] * img(i,j);

  return temp;
}

/**
 * Overloading / operator.  This function does pixel by pixel division.
 * Specified image is the dividend.
 * \ingroup overload
 * @param img Image to be divided (divisor).
 * @return Quotient of the two images.
 */
Image Image::operator/(const Image &img) const {
  int i, j, nr, nc;
  Image temp;

  nr = img.getRow();
  nc = img.getCol();

  if (nr != nrows || nc != ncols) {
    cout << "operator/: "
         << "Images are not of the same size or type, can't do division\n";
    exit(3);
  }
  temp.createImage(nrows, ncols);             
  
  for (i=0; i<nrows; i++)
    for (j=0; j<ncols; j++)
        temp(i,j) = image[i*ncols+j] / ( img(i,j) + 0.001 );

  return temp;
}


/**
 * Overloading << operator.  Output the image to the specified destination.
 * \ingroup overload
 * @param out The specified output stream (or output destination).
 * @param img Image to be output.
 * @result Output image to the specified file destination.
 */
ostream & operator<<(ostream &out, Image &img) {
  int rows, cols;
  

    for (rows = 0; rows < img.getRow(); rows++) {
      for (cols = 0; cols < img.getCol(); cols++)
        out << setw(4) << img(rows, cols) << ' ';
      out << endl;
    }

  return out; 
}

/**
 * Overloading / operator.  The left operand is the image and the right
 * is the dividend (a double point number). Each pixel in the image is 
 * divided by the double point number.
 * \ingroup overload
 * @param img Image as the left operand.
 * @param val A double point number as the right operand.
 * @result Image divided by a double point number.
 */
Image operator/(Image &img, double val) {
  int i, j, nr, nc;
  Image temp;
  
  nr = img.getRow();
  nc = img.getCol();
  temp.createImage(nr, nc);
  
  for (i=0; i<nr; i++)
    for (j=0; j<nc; j++)
        temp(i,j) = img(i,j) / val;
  
  return temp;
}

/**
 * Overloading * operator.  The left operand is the image and the right
 * is a double point number. Each pixel in the image is multiplied by the
 * double point number.
 * \ingroup overload
 * @param img Image as the left operand.
 * @param s A double point number as the right operand.
 * @result Image multiplied by a double point scalar.
 */
Image operator*(Image &img, double s) {
  int i, j, nr, nc;
  Image temp;
  
  nr = img.getRow();
  nc = img.getCol();
  temp.createImage(nr, nc);
  
  for (i=0; i<nr; i++)
    for (j=0; j<nc; j++)
        temp(i,j) = img(i,j) * s;
  
  return temp;
}


/**
 * Overloading + operator.  The left operand is the image and the right
 * is a double point number. Each pixel in the image is added by the
 * double point number.
 * \ingroup overload
 * @param img Image as the left operand.
 * @param s A double point number as the right operand.
 * @result Image add a double point scalar.
 */
Image operator+(Image &img, double s) {
  int i, j, nr, nc;
  Image temp;
  
  nr = img.getRow();
  nc = img.getCol();
  temp.createImage(nr, nc);
  
  for (i=0; i<nr; i++)
    for (j=0; j<nc; j++)
        temp(i,j) = img(i,j) + s;
  
  return temp;
}  
  
/**
 * Overloading - operator.  The left operand is the image and the right
 * is a double point number. Each pixel in the image is subtracted by the
 * double point number.
 * \ingroup overload
 * @param img Image as the left operand.
 * @param s A double point number as the right operand.
 * @result Image subtract a double point scalar.
 */
Image operator-(Image &img, double s) {
  int i, j, nr, nc;
  Image temp;
  
  nr = img.getRow();
  nc = img.getCol();
  temp.createImage(nr, nc);
  
  for (i=0; i<nr; i++)
    for (j=0; j<nc; j++)
        temp(i,j) = img(i,j) - s;
  
  return temp;
} 

/**
 * Read image from a file                     
 * @param fname The name of the file 
 * @return An Image object
 */
  void Image::readImage(char *fname) {
  ifstream ifp;
  char dummy[80];
  unsigned char *img;
  int rows, cols;
  int nRows, nCols, nt, maxi;

  ifp.open(fname, ios::in | ios::binary);

  if (!ifp) {
    cout << "readImage: Can't read image: " << fname << endl;
    exit(1);
  }

  // identify image format
  ifp.getline(dummy, 80, '\n');

  if (dummy[0] == 'P' && dummy[1] == '5') 
     ;
  else {
    cout << "readImage: Can't identify image format." << endl;
    exit(1);
  }

  // skip the comments
  ifp.getline(dummy, 80, '\n');

  while (dummy[0] == '#') {
    ifp.getline(dummy, 80, '\n');
  }

  // read the row number and column number
  sscanf(dummy, "%d %d", &nCols, &nRows);

  // read the maximum pixel value
  ifp.getline(dummy, 80, '\n');
  sscanf(dummy, "%d", &maxi); 
  if (maxi > 255) {
    cout << "Don't know what to do: maximum value is over 255.\n";
    exit(1);
  }

  if (image != NULL)
  delete [] image;
  
  nrows = nRows;
  ncols = nCols;
  maximum = 255;
  
  // read the image data
  img = (unsigned char *) new unsigned char [nRows * nCols];
  if (!img) {
    cout << "READIMAGE: Out of memory.\n";
    exit(1);
  }
  image = (float *) new float [nRows * nCols];
  if (!image) {
    cout << "READIMAGE: Out of memory.\n";
    exit(1);
  }

    ifp.read((char *)img, (nRows * nCols * sizeof(unsigned char)));
    
    for (rows = 0; rows < nRows; rows++)
      for (cols = 0; cols < nCols; cols++)
          image[rows * nCols + cols] = (float) img[rows * nCols + cols];
      
  ifp.close();
  
  delete [] img;
}


/**
 * Write image buffer to a file.
 * @param fname The output file name.
 */
void Image::writeImage(char *fname, bool flag) {
  ofstream ofp;
  int i, j;
  int nRows, nCols, nt;
  unsigned char *img;

  ofp.open(fname, ios::out | ios::binary);

  if (!ofp) {
    cout << "writeImage: Can't write image: " << fname << endl;
    exit(1);
  }


  ofp << "P5" << endl;
  ofp << ncols << " " << nrows << endl;

 
  ofp << 255 << endl;
  
  

  // convert the image data type back to unsigned char
  img = (unsigned char *) new unsigned char [nrows * ncols];
  if (!img) {
    cout << "WRITEIMAGE: Out of memory.\n";
    exit(1);
  }

  float maxi = getMaximum();
  float mini = getMinimum();
  
  
    for (i = 0; i< nrows; i++)
      for (j = 0; j < ncols; j++) {
	  // rescale if the flag is set
	  if ((maxi != mini) && flag == true)
	    img[i * ncols + j] = (unsigned char)  ((image[i * ncols + j]-mini)/(float)(maxi-mini)*255.0); 
	  // any intensity that is larger than the maximum would be set as maximum
	  else if (image[i * ncols + j] > 255)
	    img[i * ncols + j] = 255;
	  else if (image[i * ncols + j] < 0)
	    img[i * ncols + j] = 0;
	  else
	    img[i * ncols + j] = (unsigned char)  image[i * ncols + j]; 
      }
      
    ofp.write((char *)img, (nrows * ncols * sizeof(unsigned char)));


  ofp.close();
  delete [] img;
}




// YOUR FUNCTIONS

/**
 * Returns the image. 
 * @return a gray-scale image
 * \ingroup getset
 */
Image Image::thresholdImage(float thresholdValue, float lowValue, float highValue) {
  Image temp;
  int rows, cols;
  
  temp.createImage(nrows, ncols);   // temp is a gray-scale image
  for (rows = 0; rows < nrows; rows++)
    for (cols = 0; cols < ncols; cols++)
      if (image[rows * ncols + cols] <= thresholdValue) 
	temp(rows, cols) = lowValue;
      else
	temp(rows, cols) = highValue;
      
      
  return temp;
}


  //END OF YOUR FUNCTIONS //
  
// MY FUNCTÝONS   
  
// NEGATIVE FUNCTÝON

/**
 * Convert the negative
 */
Image Image::NegativeFunction() {
	Image s;    // Create new variable of Image type
	int rows, cols;   // Rows and colums
	int L = 256;   // Maximum pixel value
	int r;   // Input pixel value
	
	s.createImage(nrows, ncols);   // Allocate memory for the image
	
	// Rotate pixels
	for (rows = 0; rows < nrows; rows++)
        for (cols = 0; cols < ncols; cols++){
        	// Take the negative for each pixel value
            r = Image::getPix(rows, cols);   // Get the input pixel number 
            s(rows, cols) = L - 1 - r;   // Each pixel convert to negative pixel value
		}
        
    return s;
}

// LOG TRANSFORMATION FUNCTÝON

/**
 * Take the log transformation
 */
Image Image::LogTransformationFunction() {
	Image s;     // Create new variable of Image type
	int rows, cols;     // Rows and colums
	int c = 1;     // Stable value
	int L = 256;     // Maximum pixel value 
	int r;    // Input pixel value
	
	s.createImage(nrows, ncols);    // Allocate memory for the image
	
	// Rotate pixels
	for (rows = 0; rows < nrows; rows++)
        for (cols = 0; cols < ncols; cols++){
            r = Image::getPix(rows, cols);	  // Get the input pixel number
            s(rows, cols) = (c*log(1+r))*(L/log(1+L));   // Take the log transformation
		}
        
    return s;
}

// GAMMA TRANSFORMATION FUNCTÝON

/**
 * Take the gamma transformation
 */
Image Image::GammaTransformationFunction(float gamma) {
	Image s;     // Create new variable of Image type
	int rows, cols;     // Rows and colums
	int c = 1;     // Stable value 
	int L = 256;     // Maximum pixel value
	int r;      // Input pixel value
	
	s.createImage(nrows, ncols);     // Allocate memory for the image
	
	// Rotate pixels
	for (rows = 0; rows < nrows; rows++)
        for (cols = 0; cols < ncols; cols++){
            r = Image::getPix(rows, cols);	 // Get the input pixel number
            s(rows, cols) = c*pow(r, gamma)*(L/pow(L, gamma));    // Take the gamma transformation
		}
        
    return s;
}

// HISTOGRAMEQUALIZATION FUNCTÝON

/**
 * Do HistogramEqualization
 */
 
Image Image::HistogramEqualization() {
    Image s;    // Create new variable of Image type for image
	int rows, cols;    // Rows and colums
	int r;    // Input pixel value
	int L = 255;    // Pixel value
	float max;     // Maximum pixel value
	double h[255] {0};   // Array that calculate how much pixels
	double P[255] {0};   // Array that calculate the pdf
	double T[255] {0};   // Array that calculate the cdf
	double S1[255] {0};    // Array that calculate the output
	int S2[255] {0};    // Array that calculate the output of round
	int j;   // Variable of loop
	
	s.createImage(nrows, ncols);   // Allocate memory for the image
	
	// Loop that calculate how much pixels
	for (rows = 0; rows < nrows; rows++)
        for (cols = 0; cols < ncols; cols++){
		r = Image::getPix(rows, cols);   // Get the input pixel value
		h[r] += 1;   // Calculate how many pixel's value are there
		}
	
	// Loop that calculate the pdf
	for (j = 0; j < L+1; j++)		
		P[j] = h[j] / (ncols*nrows);
	
	// Array that calculate the cdf
	T[0] = P[0]; // Assigning the initial value for cdf
	for (j = 1; j < L+1; j++)  
        T[j] = P[j] + T[j-1];    // Calculate the cdf
    
    max = Image::getMaximum();   // Find the maximum pixel value
    
    // Array that calculate the output and output of round
    for (j = 0; j < L+1; j++) {
        S1[j] =  T[j] * max;    // Calculate the output array
        S2[j] = round(S1[j]);     // Calculate the round output array
    }
    
    // Create the output image
    for (rows = 0; rows < nrows; rows++)
        for (cols = 0; cols < ncols; cols++){
		r = Image::getPix(rows, cols);   // Get the input pixel value
		s.setPix(rows, cols, S2[r]);   // Set the output pixel value
    }
    
    return s;
}

// FOURIERTRANSFORM FUNCTÝON


// NEW FIGURE
void Image::NewFigure(float initialValue) {
  int i;
  
  for (i = 0; i < nrows * ncols; i++){
    image[i] = initialValue;
    initialValue++;
  }
}

// LAPLACIAN
Image Image::Laplacian(Image& input){
    Image temp;   // Create new variable of Image type for image
    temp = input;
    int rows, cols;    // Rows and colums
	int k, l;    // For mask loop
    int MaskSize = 3;   // Mask Size
    
    // Laplacian Mask
    float mask[MaskSize][MaskSize] = {
        {-1, -1, -1},
        {-1,  8, -1},
        {-1, -1, -1}
    };
    
    // Implement mask to Image
    for(rows = 1; rows<nrows-1; rows++){
      for(cols = 1; cols<ncols-1; cols++){
      	temp(rows, cols) = 0;   // sum input and output
      	   for(k = 0; k<MaskSize; k++){
      	      for(l = 0; l<MaskSize; l++){
      	  	    temp(rows, cols) += input(rows - 1 + k, cols - 1 + l) * mask[k][l]; // Implement mask to Image
		  }
		}
      }
    }
  
    return temp;
}

// SHARPED
Image Image::Sharped(Image& input, Image& output){
	Image temp;   // Create new variable of Image type for image
	
	temp = input + output;   // sum input and output
	
	return temp;
}

// SOBEL
Image Image::Sobel(Image& input){
	Image temp, temp1, temp2;   // Create new variable of Image type for image
	temp = input;
	temp1 = input;
	temp2 = input;
	int rows, cols;   //  Rows and colums
	int k, l;   // For mask loop
	int MaskSize = 3;   // Mask Size
	
	// Mask for x
	float gx[MaskSize][MaskSize] = {
	    {-1, -2, -1},
	    {0,   0,  0},
	    {1,   2,  1}
	};
	
	// Implement mask to Image
	for(rows = 1; rows<nrows-1; rows++){
      for(cols = 1; cols<ncols-1; cols++){
      	temp(rows, cols) = 0;
      	   for(k = 0; k<MaskSize; k++){
      	      for(l = 0; l<MaskSize; l++){
      	  	    temp1(rows, cols) += input(rows - 1 + k, cols - 1 + l) * gx[k][l];   // Implement mask to Image
		  }
		}
      }
    }
	
	// Mask for y
	float gy[MaskSize][MaskSize] = {
	    {-1, 0, 1},
	    {-2, 0, 2},
	    {-1, 0, 1}
	};
	
	// Implement mask to Image
	for(rows = 1; rows<nrows-1; rows++){
      for(cols = 1; cols<ncols-1; cols++){
      	temp(rows, cols) = 0;
      	   for(k = 0; k<MaskSize; k++){
      	      for(l = 0; l<MaskSize; l++){
      	  	    temp2(rows, cols) += input(rows - 1 + k, cols - 1 + l) * gy[k][l];   // Implement mask to Image
		  }
		}
      }
    }
    
	// 	sum of amplitudes of gx and gy
	for(rows = 1; rows < nrows-1; rows++){
		for(cols = 1; cols < ncols-1; cols++){
			temp(rows, cols) = abs(temp1(rows, cols)) + abs(temp2(rows, cols));
		}
	}
    
	return temp;
}

// AVERAGINGFILTER
Image Image::AveragingFilter(Image& input){
	Image temp;  // Create new variable of Image type for image
	temp = input;
	
	int rows, cols;    //  Rows and colums
	int k, l;    // For mask loop
	int MaskSize = 5;  // Mask Size
	float Average = (1.0)/(MaskSize * MaskSize);  // Average
	
	// Avarege filter
	float i[MaskSize][MaskSize] = {
	    {1, 1, 1, 1, 1},
	    {1, 1, 1, 1, 1},
	    {1, 1, 1, 1, 1},
	    {1, 1, 1, 1, 1},
		{1, 1, 1, 1, 1}
	};
	
	// Implement mask to Image
	for(rows = 1; rows<nrows; rows++){
      for(cols = 1; cols<ncols; cols++){
      	temp(rows, cols) = 0;
      	   for(k = 0; k<MaskSize; k++){
      	      for(l = 0; l<MaskSize; l++){
      	  	    temp(rows, cols) += (input(rows - 1 + k, cols - 1 + l) * i[k][l]) * Average;  // Implement mask to Image
		  }
		}
      }
    }
	
	return temp;
}

//
    //END OF MY FUNCTIONS //
