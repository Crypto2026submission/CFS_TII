We give three files: 
- attack_83.magma which attacks TII instance 83
- attack_89.magma which attacks TII instance 89
- attack_125.magma which attacks TII instance 125

attack_83.magma and attack_89.magma stop when a final GRS which is a power of the desired GRS is computed.
A standard [SS92] then finishes the job.

For attack_125.magma, the program ends when a valid support and Goppa polynomial have been recovered.

We also provide a file called Cmat.magma which generates a matrix code of quadratic relations for a desired alternant or 
Goppa code, and writes it into a file Cmat.txt. This is for the implementation in C. Then one needs to run 

python3 ./parse_cmat.py 

Once this is done, copy Cmat.magma into the C/ folder. The C program then performs the first system echelonization of Algorithm 2,
and write the obtained linear equations into stdout.
An example of Cmat.txt is already given in the C folder. It corresponds to a shortened Goppa code with r=6 and m=10.

To run the C implem, just type 

make
./main