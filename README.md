# Vascular Network Generator

This code generates a circular vascular network that inspired by the retinal network of a zebra fish. The idea of Fibonacci trees are used to generate the network. See the [design process](https://github.com/bobbielf2/vascular-net/blob/master/images/design_process.pdf).

Main author: Bowei (Bobbie) Wu

Also includes the following contributions and influces:

Bowei Wu - NURBS or B-splines polygon smoothing [Link](https://github.com/bobbielf2/Polygon-BSplines)  
Hai Zhu - polygon smoothing based on partition of unity (POU)   
Massimo Zanetti - smoothing cubic splines with periodic conditions [Link](https://www.mathworks.com/matlabcentral/fileexchange/59463-smoothing-cubic-splines-with-periodic-conditions)   
Nick Trefethen - random smooth function

### Example network

<img src=https://github.com/bobbielf2/vascular-net/raw/master/images/fibonacci_net.png width=200 height=200> <img src=https://github.com/bobbielf2/vascular-net/raw/master/images/retinal_network.png width=200 height=200>

* Left: network generated by this code.
* Right: retinal vascular network of a zebra fish, from [Alvarez et al., 2007]

## Reference

Alvarez et al. (2007). Genetic determinants of hyaloid and retinal vasculature in zebrafish. BMC developmental biology, 7(1):114.
