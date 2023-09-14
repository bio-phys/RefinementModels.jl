import Pkg

Pkg.activate("../../BioEn.jl")
import BioEn: Input 

Pkg.activate("../../RefinementModels.jl")

import RefinementModels: Gauss

path = "../data/Gauss/"

N = convert(Int64, 1e1)
M = convert(Int64, 1e0)

mu_Y = 0.
sigma_Y = 1.

sigma_y = 2. 
offset_y = 1.

if !isdir(path)
    println("Path \"$path\" does not exist. Exiting.")
    exit(-1)
end

Y = Gauss.experimental(M; mu=mu_Y, sigma=sigma_Y);
y = Gauss.theoretical(N, Y; sigma=sigma_y, offset=offset_y);

w0 = ones(N)/N;

Input.write_txt_normalized(path; Y_name="exp.txt", Y=Y, y_name="obs.txt", y=y, w0_name="w0.txt", w0=w0)
