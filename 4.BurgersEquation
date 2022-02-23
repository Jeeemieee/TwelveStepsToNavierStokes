#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>


void displayVector (std::vector<double> &v);
void saveDataTextFile (std::vector<double> &coordinate, std::vector<double> &data, std::string &name);

int main() {
    //Input parameters
    int nx{101}; //Number of grid points
    double dx{2.0*M_PI/(nx-1.0)}; //Distance between grid points for a spatial domain of 2 units of length
    int nt{25}; //Number of time steps
    double sigma{0.2}; // time step size coefficient
    double nu{0.3}; //Kinematic viscosity
    double dt{sigma*pow(dx,2)/nu}; //Time step size
    double t{0};

    //Initial conditions
    std::vector<double> x (nx, 0); //Position vector
    for (int i{1}; i<nx; i++)
    {
        x[i] = i*dx;
    }
    std::vector<double> phi (x.size(),0);
    std::vector<double> phiPrime (x.size(),0);
    std::vector<double> u (x.size(),0);
    for (int i{0}; i<x.size(); i++)
    {
        phi[i] = exp(-pow((x[i]-4*t),2)/(4*nu*(t+1)))+exp(-pow((x[i]-4*t-2*M_PI),2)/(4*nu*(t+1)));
        phiPrime[i] = ((-2*(x[i]-4*t))/(4*nu*(t+1)))*exp(-pow((x[i]-4*t),2)/(4*nu*(t+1)))+((-2*(x[i]-4*t-2*M_PI))/(4*nu*(t+1)))*exp(-pow((x[i]-4*t-2*M_PI),2)/(4*nu*(t+1)));
        u[i] = -((2*nu)/(phi[i]))*phiPrime[i]+4;
    }

    //Temporary array
    std::vector<double> un (nx,1);//Velocity data previous time step
    for (int i{0}; i<nt; i++)
    {
        for (int j{0}; j<un.size(); j++)
        {
            un[j] = u[j];
        }
        for (int j{1}; j<nx-1; j++)
        {
            u[j] = un[j] - un[j]*(dt/dx)*(un[j]-un[j-1]) + nu*dt*pow(dx,-2)*(un[j+1]-2*un[j]+un[j-1]);
        }
    }


    std::string fileName{"4_BurgersEquationData.txt"};
    saveDataTextFile(x,u,fileName);
    return 0;
}
void displayVector (std::vector<double> &v)
{
    for(int i{0}; i < v.size();i++)
    {
        std::cout << v[i] << " ";
    }
}
void saveDataTextFile (std::vector<double> &coordinate, std::vector<double> &data, std::string &name)
{
    std::ofstream file_;
    file_.open(name);
    for (int i{0}; i<data.size(); i++)
    {
        file_ << coordinate[i] << "," << data[i];
        file_ << "\n";
    }
    file_.close();
}
