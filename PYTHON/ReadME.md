Basic steps of the Python program for "RaPlace":

First, after generating backward-warping Cartesian images from the polar image data, we obtain sinogram images using the Radon Transform on the Cartesian radar data. In the Python program, we implement it by calling the "radon" function in the skimage.transform library.

Then, we operate discrete FFT on the acquired sinograms, by using the function "numpy.fft.fft(variable_name, axis=0)", and utilize the FFT images as candidates.

Next, in order to determine the similarity score, we adopt the cross-correlation in frequency domain and mutable threshold, which are calculated from the auto-correlation value, for all previous keyframes when query images are published. In the Python program, we used the written function "fast_dft(Mq, Mi)" to implement the above functions. 

Finally, the candidate location is determined as the keyframe with the lowest distance, assessed from the similarity scores. And we can use the result to draw the P-R curve at one time.

In this program, there are two ways to implement the Radon Transformation. If the function "radon(tmp, theta)" import from "skimage.transform" cannot be used, you can use the function "DiscreteRadonTransform(image, steps)" to implement the RT. But it should be noted that, the "theta" is an array, while the "steps" is a constant, which represents the maximum value of angle.
