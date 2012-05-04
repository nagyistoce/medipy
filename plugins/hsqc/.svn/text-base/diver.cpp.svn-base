//#include <codecogs/maths/optimization/golden.h>
 
#include <iostream>
#include <iomanip>
#include <cmath>
 
// user-defined function
double f(double x) {
    return x * sin(x) - 2 * cos(x);
}
 
int main() 
{
    double x = Maths::Optimization::golden(f, 5, 6, 7);
    std::cout << "Function minimum is Y = " << std::setprecision(13) << f(x);
    std::cout << std::endl;
    std::cout << "for X = " << x << std::endl;
    return 0;
}