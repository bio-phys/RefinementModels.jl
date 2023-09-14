import Pkg

Pkg.activate("../../BioEn.jl")
import BioEn: Input, Utils 

Pkg.activate("../../RefinementModels.jl")

import RefinementModels: ContDoubleWell

path = "../data/ContDoubleWell" 

a = 3. 
x0 = 1.
half_width = 1.5

N=10000
Yo = [0.8] 
sigmas = [0.2]

if !isdir(path)
    println("Path \"$path\" does not exist. Exiting.")
    exit(-1)
end

M = length(Yo)

ContDoubleWell.distribution_max(a, x0)
fac = 5
x = ContDoubleWell.sample(fac*N, half_width=half_width, a=a, x0=x0)
len = length(x)
len < N && println("Not enough samples (N>$len)!. Exiting.")  && exit(-1)
sample = x[range(1,N)];

Y = Utils.normalize(Yo, sigmas)

yo = zeros(N, M)
yo[:,1] = sample
y = Utils.normalize(yo, sigmas);

w0 = ones(N)/N

Input.write_txt_normalized(path; Y_name="exp.txt", Y=Y, y_name="obs.txt", y=y, w0_name="w0.txt", w0=w0)

