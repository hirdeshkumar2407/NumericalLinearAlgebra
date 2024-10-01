#include <Eigen/Dense>
#include <iostream>
#include <cstdlib>

// from https://github.com/nothings/stb/tree/master
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

using namespace Eigen;

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <image_path>" << std::endl;
    return 1;
  }

  const char* input_image_path = argv[1];

  // Load the image using stb_image
  int width, height, channels;
  // for greyscale images force to load only one channel
  unsigned char* image_data = stbi_load(input_image_path, &width, &height, &channels, 1);  
  if (!image_data) {
    std::cerr << "Error: Could not load image " << input_image_path
              << std::endl;
    return 1;
  }

  std::cout << "Image loaded: " << width << "x" << height << " with "
            << channels << " channels." << std::endl;

  // Prepare Eigen matrices 
  MatrixXd dark(height, width), light(height, width), rotate(width, height);


/*Purpose:
The raw image data is loaded as a 1D array of pixel values (image_data), but to apply transformations (like darkening or rotating), we need to fill the prepared Eigen matrices (dark, light, rotate) with these pixel values.
The loop reads each pixel from the image_data array and places it in the corresponding position in the matrices:
dark: stores the darkened version of the image.
light: stores the lightened version of the image (though it's not saved in this case).
rotate: stores a rotated version of the image.
Why needed:
The pixel values need to be transformed from raw data into a matrix form for easier manipulation. Each pixel's brightness is adjusted and mapped to a new value in the corresponding matrix.
Without filling the matrices, you wouldn't be able to perform complex operations like darkening or rotating the image.*/

  // Fill the matrices with image data
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      int index = (i * width + j) * channels;  // 1 channel (Greyscale) 3 channels (RGB)
      dark(i, j) = std::max(static_cast<double>(image_data[index]) - 50.,0.0) / 255.0;

      /*1. Darkening the Image (dark matrix)

dark(i, j) = std::max(static_cast<double>(image_data[index]) - 50., 0.0) / 255.0;
Logic:
Goal: The dark matrix stores the darkened version of the image by reducing the brightness of each pixel.
Explanation:
image_data[index]: This retrieves the pixel value from the raw image data (which ranges from 0 to 255, since it's grayscale).
static_cast<double>(image_data[index]) - 50.: This reduces the brightness of the pixel by subtracting 50 from its original value.
The static_cast<double> ensures the value is treated as a floating-point number.
The subtraction of 50 darkens the pixel.
std::max(..., 0.0): Ensures that the pixel value doesn't go below 0 (since pixel values can't be negative).
/ 255.0: Normalizes the pixel value to the range [0.0, 1.0], which is more convenient for mathematical operations in Eigen matrices.
Example:
If a pixel value is 100, the darkened value would be max(100 - 50, 0) = 50. After normalization: 50 / 255.0 ≈ 0.196.
If a pixel value is 30, the darkened value would be max(30 - 50, 0) = 0. After normalization: 0 / 255.0 = 0.
The dark matrix thus contains normalized pixel values where each pixel's brightness has been reduced, but no pixel goes below 0 (black).*/
      light(i, j) = std::min(static_cast<double>(image_data[index]) + 50.,255.) / 255.0;

      /*
Let's break down the logic behind the dark, light, and rotate transformations:

1. Darkening the Image (dark matrix)
cpp
Copy code
dark(i, j) = std::max(static_cast<double>(image_data[index]) - 50., 0.0) / 255.0;
Logic:
Goal: The dark matrix stores the darkened version of the image by reducing the brightness of each pixel.
Explanation:
image_data[index]: This retrieves the pixel value from the raw image data (which ranges from 0 to 255, since it's grayscale).
static_cast<double>(image_data[index]) - 50.: This reduces the brightness of the pixel by subtracting 50 from its original value.
The static_cast<double> ensures the value is treated as a floating-point number.
The subtraction of 50 darkens the pixel.
std::max(..., 0.0): Ensures that the pixel value doesn't go below 0 (since pixel values can't be negative).
/ 255.0: Normalizes the pixel value to the range [0.0, 1.0], which is more convenient for mathematical operations in Eigen matrices.
Example:
If a pixel value is 100, the darkened value would be max(100 - 50, 0) = 50. After normalization: 50 / 255.0 ≈ 0.196.
If a pixel value is 30, the darkened value would be max(30 - 50, 0) = 0. After normalization: 0 / 255.0 = 0.
The dark matrix thus contains normalized pixel values where each pixel's brightness has been reduced, but no pixel goes below 0 (black).

2. Lightening the Image (light matrix)
cpp
Copy code
light(i, j) = std::min(static_cast<double>(image_data[index]) + 50., 255.) / 255.0;
Logic:
Goal: The light matrix stores the lightened version of the image by increasing the brightness of each pixel.
Explanation:
static_cast<double>(image_data[index]) + 50.: This increases the brightness of the pixel by adding 50 to its original value.
std::min(..., 255.): Ensures that the pixel value doesn't exceed 255 (since 255 is the maximum brightness).
/ 255.0: Normalizes the pixel value to the range [0.0, 1.0], which is useful for further processing in Eigen matrices.
Example:
If a pixel value is 100, the lightened value would be min(100 + 50, 255) = 150. After normalization: 150 / 255.0 ≈ 0.588.
If a pixel value is 230, the lightened value would be min(230 + 50, 255) = 255. After normalization: 255 / 255.0 = 1.0.
The light matrix contains normalized pixel values where each pixel's brightness has been increased, but no pixel goes above 255 (white).*/
      rotate(width-j-1, i) = static_cast<double>(image_data[index]) / 255.0;

      /*. Rotating the Image (rotate matrix)
cpp
Copy code
rotate(width-j-1, i) = static_cast<double>(image_data[index]) / 255.0;
Logic:
Goal: The rotate matrix stores a rotated version of the image by rearranging the pixel positions. In this case, it rotates the image by 90 degrees counterclockwise.
Explanation:
image_data[index]: This retrieves the original pixel value.
static_cast<double>(image_data[index]) / 255.0: This normalizes the pixel value to the range [0.0, 1.0] as before.
rotate(width-j-1, i): This is where the rotation logic happens. Instead of placing the pixel at (i, j), the pixel is placed at a new position that rotates it by 90 degrees counterclockwise.
The new row index is width - j - 1, which effectively shifts the pixel from the right side of the image to the top in the rotated image.
The new column index is i, which keeps the row information but now represents a column in the rotated image.
Example:
Consider a pixel at position (i, j) in the original image. In the rotated image, this pixel will be placed at position (width - j - 1, i), rotating the image counterclockwise by 90 degrees.
For example, for a 3x3 image:

mathematica
Copy code
Original Image (positions):
(0,0) (0,1) (0,2)
(1,0) (1,1) (1,2)
(2,0) (2,1) (2,2)

After Rotation:
(2,0) (1,0) (0,0)
(2,1) (1,1) (0,1)
(2,2) (1,2) (0,2)
Thus, the rotate matrix contains the normalized pixel values of the original image rearranged to reflect a 90-degree counterclockwise rotation.*/
    }
  }

  
  // Free memory!!!
  stbi_image_free(image_data);

  Matrix<unsigned char, Dynamic, Dynamic, RowMajor> dark_image(height, width);

  /*Purpose:
The matrices (dark, rotate, etc.) hold pixel values in a normalized format (between 0.0 and 1.0). However, when saving the image back to disk, pixel values must be in the unsigned char format (0 to 255).
This mapping step uses Eigen's unaryExpr() function to convert each pixel value back from the normalized [0.0, 1.0] range to the standard [0, 255] range required by image files.
Why needed:
When saving the image to disk (with stbi_write_png), the pixel values must be in a specific format (0-255 integers). Since we manipulated them in a normalized format (as doubles), we need to convert them back to the correct range and data type.
Without this step, you can't save the image in the correct format, as image files require pixel data to be in the integer range of [0, 255].*/
  // Use Eigen's unaryExpr to map the grayscale values (0.0 to 1.0) to 0 to 255
  dark_image = dark.unaryExpr([](double val) -> unsigned char {
    return static_cast<unsigned char>(val * 255.0);
  });

  // Save the image using stb_image_write
  const std::string output_image_path1 = "dark_image.png";
  if (stbi_write_png(output_image_path1.c_str(), width, height, 1,
                     dark_image.data(), width) == 0) {
    std::cerr << "Error: Could not save grayscale image" << std::endl;

    return 1;
  }

  Matrix<unsigned char, Dynamic, Dynamic, RowMajor> rotate_image(width, height);
  // Use Eigen's unaryExpr to map the grayscale values (0.0 to 1.0) to 0 to 255
  rotate_image = rotate.unaryExpr([](double val) -> unsigned char {
    return static_cast<unsigned char>(val * 255.0);
  });

  // Save the image using stb_image_write
  const std::string output_image_path2 = "rotate_image.png";
  if (stbi_write_png(output_image_path2.c_str(), height, width, 1,
                     rotate_image.data(), height) == 0) {
    std::cerr << "Error: Could not save output image" << std::endl;
    
    return 1;
  }

  std::cout << "Images saved to " << output_image_path1 << " and " << output_image_path2 << std::endl;

  return 0;
}
