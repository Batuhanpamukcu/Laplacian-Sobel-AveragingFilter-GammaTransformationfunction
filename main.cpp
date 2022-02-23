/**************************************
 The original copy of the code can be found at http://web.eecs.utk.edu/~mkarakay/courses_files/testfiles.zip 
 and it is modified for ELM463/667
 **************************************/

#include "Image.h"

using namespace std;

int main(int argc, char **argv)
{
	//SORU1
	//b
	Image inputImage, outputImage;  // Create new object of Image type for input image and output Image
  
    char inputFileName[] = "Fig0343(a).pgm";   // File Name for input image
    char outputFileName[] = "Fig0343(b).pgm";   // File Name for output image
    
    inputImage.readImage(inputFileName);   // Read the Image
    outputImage = inputImage.Laplacian(inputImage); // Run the Laplacian function
    outputImage.writeImage(outputFileName, true);  // Write the image
    
    //c
    Image outputImage2;  // Create new object of Image type for output Image
    
    char outputFileName2[] = "Fig0343(c).pgm";  // File Name for output image
    
    outputImage2 = outputImage.Sharped(inputImage, outputImage);  // Run the Sharped function
    outputImage2.writeImage(outputFileName2, true);  // Write the image
    
    //d
    Image inputImage3, outputImage3;   // Create new object of Image type for input image and output Image
    
    char outputFileName3[] = "Fig0343(d).pgm";   // File Name for output image
    
    inputImage3.readImage(inputFileName);   // Read the Image
    outputImage3 = inputImage.Sobel(inputImage);   // Run the Sobel function
    outputImage3.writeImage(outputFileName3, true);   // Write the image
    
    //e
    Image inputImage4, outputImage4;   // Create new object of Image type for input image and output Image
    
    char outputFileName4[] = "Fig0343(e).pgm";   // File Name for output image
    
    outputImage4 = outputImage3.AveragingFilter(outputImage3);   // Run the AveragingFilter function
    outputImage4.writeImage(outputFileName4, true);   // Write the image
    
    //f
    Image inputImage5, outputImage5;   // Create new object of Image type for input image and output Image
    
    char outputFileName5[] = "Fig0343(f).pgm";   // File Name for output image
    
    outputImage5 = outputImage2 * outputImage4;   // product outputImage2 and outputImage4
    outputImage5 = outputImage5 - outputImage5.getMinimum();   // Scaling
    outputImage5 = outputImage5 * (255/outputImage5.getMaximum());   // Scaling
    outputImage5.writeImage(outputFileName5, true);   // Write the image
    
    //g
    Image inputImage6, outputImage6;    // Create new object of Image type for input image and output Image
    
    char outputFileName6[] = "Fig0343(g).pgm";    // File Name for output image
    
    outputImage6 = outputImage6.Sharped(inputImage, outputImage5);  // Sum inputImage and outputImage5
	outputImage6 = outputImage6 - outputImage6.getMinimum();   // Scaling
    outputImage6 = outputImage6 * (255/outputImage6.getMaximum());    // Scaling
    outputImage6.writeImage(outputFileName6, true);   // Write the image
    
    //h
    Image inputImage7, outputImage7;   // Create new object of Image type for input image and output Image
    
    char outputFileName7[] = "Fig0343(h).pgm";    // File Name for output image
    
    outputImage7 = outputImage6.GammaTransformationFunction(0.5);   // Run the AveragingFilter function
    outputImage7.writeImage(outputFileName7, true);    // Write the image
}
