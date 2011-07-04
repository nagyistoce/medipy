import numpy

def addition(input1, input2, output) :
   """
   Add the two input images and store the result in the output image.

   :gui:
       input1 : Image
           First input

       input2 : Image
           Second input

       output : Image : output=True
           Output image
   """
   if input1.shape != input2.shape :
       raise Exception("Inputs must have the same shape")
   if input1.shape != output.shape :
       output.data = numpy.ndarray(shape=input1.shape, dtype=input1.dtype)
   output[:] = input1[:] + input2[:]

def subtraction(input1, input2, output) :
   """
   Store the result of input1 - input2 in the output image.

   :gui:
       input1 : Image
           First input

       input2 : Image
           Second input

       output : Image : output=True
           Output image
   """
   if input1.shape != input2.shape :
       raise Exception("Inputs must have the same shape")
   if input1.shape != output.shape :
       output.data = numpy.ndarray(shape=input1.shape)
   output.data = numpy.ndarray(shape=input1.shape, dtype=input1.dtype)
   output[:] = input1[:] - input2[:]

def multiplication(input, scalar, output) :
   """
   Store the result of input * scalar in the output image.

   :gui:
       input : Image
           Input

       scalar : Int
           Scalar

       output : Image : output=True
           Output image
   """
   if input.shape != output.shape :
       output.data = numpy.ndarray(shape=input.shape, dtype=input.dtype)
   output[:] = input[:] * scalar

def division(input, scalar, output) :
   """
   Store the result of input / scalar in the output image.

   :gui:
       input : Image
           Input

       scalar : Int
           Scalar

       output : Image : output=True
           Output image
   """
   if input.shape != output.shape :
       output.data = numpy.ndarray(shape=input.shape, dtype=input.dtype)
   output[:] = input[:] / scalar

def absolute_value(input, output):
    """ Compute the absolute value element_wise
        
        :gui:
            input : Image
                Input
            output : Image : output=True
                Output
    """
    if input.shape != output.shape :
       output.data = numpy.ndarray(shape=input.shape, dtype=input.dtype)
    output[:] = numpy.abs(input)
    output.copy_information(input)
    
def addition_scalar(input, scalar, output) :
  """
  Store the result of input + scalar in the ouput image

  :gui:
      input : Image
	  Input

      scalar : Int
	  Scalar

      output : Image : output=True
	  Output image
  """
  if input.shape != output.shape :
    output.data = numpy.ndarray(shape=input.shape, dtype=input.dtype)
  output[:] = input[:] + scalar


def substraction_scalar(input, scalar, output) :
  """
  Store the result of input - scalar in the ouput image

  :gui:
      input : Image
	  Input

      scalar : Int
	  Scalar

      output : Image : output=True
	  Output image
  """
  if input.shape != output.shape :
    output.data = numpy.ndarray(shape=input.shape, dtype=input.dtype)
  output[:] = input[:] - scalar

def multiplication_image(input1, input2, output) :
   """
   Store the result of input1 * input2 in the output image.

   :gui:
       input1 : Image
           First input

       input2 : Image
           Second input

       output : Image : output=True
           Output image
   """
   if input1.shape != input2.shape :
       raise Exception("Inputs must have the same shape")
   if input1.shape != output.shape :
       output.data = numpy.ndarray(shape=input1.shape)
   output.data = numpy.ndarray(shape=input1.shape, dtype=input1.dtype)
   output[:] = input1[:] * input2[:]

def division_image(input1, input2, output) :
   """
   Store the result of input1 / input2 in the output image.

   :gui:
       input1 : Image
           First input

       input2 : Image
           Second input

       output : Image : output=True
           Output image
   """
   if input1.shape != input2.shape :
    raise Exception("Inputs must have the same shape")
   if input1.shape != output.shape :
    output.data = numpy.ndarray(shape=input1.shape)
   output.data = numpy.ndarray(shape=input1.shape, dtype=input1.dtype)
   shape = input2.shape
   i=-1
   j=-1
   k=-1
   while i < shape[2]-1 :
      i+=1
      j=-1
      k=-1
      while j < shape[1]-1 :
	  j+=1
	  k=-1
	  while k < shape[0]-1 :
	      k+=1
	      if input2[i,j,k] != 0 :
		  output[i,j,k] = input1[i,j,k] / input2[i,j,k]
	      else :
		  output[i,j,k]=0


