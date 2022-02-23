#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>


void displayVector (std::vector<double> &v);
void saveDataTextFile (std::vector<double> &coordinate, std::vector<double> &data, std::string &name);

int main() {
    //Input parameters
    int nx{41}; //Number of grid points
    double dx{2.0/(nx-1.0)}; //Distance between grid points for a spatial domain of 2 units of length
    int nt{25}; //Number of time steps
    double dt{0.025}; //Time step size
    int c{1}; //Wave speed

    //Initial conditions
    std::vector<double> u (nx,1); //Velocity vector
    std::fill_n(u.begin()+round(0.5/dx),round(1.0/dx + 1) - round(0.5/dx), 2.0);
    std::vector<double> x (nx, 0); //Position vector
    for (int i{1}; i<nx; i++)
    {
        x[i] = i*dx;
    }

    //Temporary array
    std::vector<double> un (nx,1);//Velocity data previous time step
    for (int i{0}; i<nt; i++)
    {
        for (int j{0}; j<un.size(); j++)
        {
            un[j] = u[j];
        }
        for (int j{1}; j<nx; j++)
        {//calculate velocities
            u[j] = un[j] - c*(dt/dx)*(un[j]-un[j-1]);
        }
    }
    std::string fileName{"1_LinearConvectionData.txt"};
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
