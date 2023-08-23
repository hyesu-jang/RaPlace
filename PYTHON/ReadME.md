# Basic steps of the Python program for "RaPlace"

First, after generating backward-warping Cartesian images from the polar image data, we obtain sinogram images using the Radon Transform on the Cartesian radar data. In the Python program, we implement it by calling the "radon" function in the skimage.transform library.

Then, we operate discrete FFT on the acquired sinograms, by using the function "numpy.fft.fft(variable_name, axis=0)", and utilize the FFT images as candidates.

Next, in order to determine the similarity score, we adopt the cross-correlation in frequency domain and mutable threshold, which are calculated from the auto-correlation value, for all previous keyframes when query images are published. In the Python program, we used the written function "fast_dft(Mq, Mi)" to implement the above functions. 

Finally, the candidate location is determined as the keyframe with the lowest distance, assessed from the similarity scores. And we can use the result to draw the P-R curve at one time.

# Instructions of some functions and variables

In this program, there are two ways to implement the Radon Transformation. If the function "radon(tmp, theta)" import from "skimage.transform" cannot be used, you can use the function "DiscreteRadonTransform(image, steps)" to implement the RT. But it should be noted that, the "theta" is an array, while the "steps" is a constant, which represents the maximum value of angle.

The variable "radar_data_dir" in the function "generateRadon" and "data_path" in the main function represent the the folder for storing radar images. The variable "gtpose" is got from loading "global_pose.csv".

# Package used

numpy, math, os, cv2, matplotlib.pyplot, radon (import from skimage.transform), ndimage (import from scipy)

**Make sure that you have installed these packages before running the code!!!**

# ATTENTION

During the operation of the program, there may be some "warning"s. If they do not affect the operation of the program, **please ignore the "warning"s**.

The possible warnings can be listed as follows:

<img width="713" alt="image" src="https://github.com/baorrr2020/RaPlace/assets/142761589/13227d17-de60-4d09-8dfa-db2d38f38770">

<img width="930" alt="image" src="https://github.com/baorrr2020/RaPlace/assets/142761589/883b414c-9e1f-40f3-89d3-f0bdf274e431">

<img width="547" alt="image" src="https://github.com/baorrr2020/RaPlace/assets/142761589/80644b9b-5b97-418a-a5b7-a1b33d607fef">

# Contact

If you have any questions about the python code, please contact me through sending e-mail to 
        
    mibh2020@163.com
