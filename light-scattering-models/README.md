# Light Scattering Models

This repo hold different models for light scattering in a medium. We'll star off with a simple scattering model and move on to more complicated ones. A few examples of 
applications for these scattering models with be present here. The main code should be exportable for other applications.

### Compton
We will begin with compton scattering as it is a simple model to compute and doesn't require any complex calculations. The description of each function can be called in MATLAB by
writitng `help fcn_name`. Here is an example output from the Compton to help:
```
>> help compton_scatter
 the purpose of this function is to compute the compton scattering angle
 for a given wavelength of light. It is a very basic model for how light is
 scattered, which makes it the perfect location to begin. We will model the
 new wavelength as a function of input wavelength and scatter angle. If the
 user desires the scattered wavelength as a function of angle with a
 constant input wavelength, or vice versa, they can use this function. We
 can also input both parameters to create a heat map for this.
    Inputs: lambda:double, incident wavelength [m]
            theta:double, scatter angle [deg]
    Outputs:lambda_prime: double, sacttered wavelength of the light [m]

>> 
```
The output will give the input arguments required and the outputs along with their dimensions. Compton scattering is used to illustrate a simpler phenomena, Thompson is more
complex.