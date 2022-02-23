#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

std::vector<double> linspace(double start, double end, double steps);
void displayVector (std::vector<double> &v);
void displayVector (std::vector<int> &v);
std::vector<int> arrange(int start, int end, int stepSize);
void saveTimePositionData (std::vector<int> &time, std::vector<double> &x, std::vector<double> &y, std::string &fileName);
std::vector<std::vector<double>> matrixBoundaryValue (std::vector<std::vector<double>> &matrix, double boundaryValue);
void saveVelocityData (std::vector<std::vector<double>> &u, std::string &fileName);

int main() {
    //simulation variables
    int nx{51};
    int ny{51};
    int nt{200};
    double dx{2.0/(nx-1.0)};
    double dy{2.0/(ny-1.0)};
    double sigma{0.25};
    double nu{0.05};
    double dt{sigma*dx*dy*pow(nu,-1)};

    //storage variables
    int saveFactor{10}; //factor to save only a fraction of the total time steps
    std::string fileName{"7_DiffusionData.txt"};
    std::ofstream file_;//Produce the file with the defined name and close it, the functions will re-open the file, append with data and close it again.
    file_.open(fileName);
    file_.close();
    int n{0};

    //vector initialization
    std::vector<int> timeSteps = arrange(0,nt,saveFactor);
    std::vector<double> x = linspace(0,2,nx);
    std::vector<double> y = linspace(0,2,ny);
    saveTimePositionData(timeSteps,x,y,fileName);
    std::vector<std::vector<double>> u (ny,std::vector<double>(nx,1));
    std::vector<std::vector<double>> un (ny,std::vector<double>(nx,1));

    //data input
    for (int i{static_cast<int>(0.5/dy)};i<(round(1.0/dy + 1));i++)
    {//initial condition
        std::fill_n(u[i].begin()+round(0.5/dx),round(1.0/dx + 1) - round(0.5/dx), 2.0);
    }

    for (int t{0};t<nt;t++)
    {
        for (int i{0}; i<y.size();i++)
        {
            un[i] = u[i];
        }
        for (int j{1};j<y.size()-1;j++)
        {
            for (int i{1};i<x.size()-1;i++)
            {
                u[j][i] = un[j][i] + nu*dt*pow(dx,-2)*(un[j][i+1] - 2*un[j][i] + un[j][i-1]) + nu*dt*pow(dy,-2)*(un[j+1][i] - 2*un[j][i] + un[j-1][i]);
            }
        }
        u = matrixBoundaryValue(u,1.0);
        if ( t == timeSteps[n])
        {
            saveVelocityData(u,fileName);
            n++;
        }
    }

    return 0;
}
std::vector<double> linspace(double start, double end, double steps)
{
    double spacing{(end - start)/(steps-1)};
    std::vector<double> vector (steps, start);
    for (int i{1};i<steps;i++)
    {
        vector[i] += i*spacing;
    }

    return vector;
}
void displayVector (std::vector<double> &v)
{
    for(int i{0}; i < v.size();i++)
    {
        std::cout << v[i] << " ";
    }
}
void displayVector (std::vector<int> &v)
{
    for(int i{0}; i < v.size();i++)
    {
        std::cout << v[i] << " ";
    }
}
std::vector<int> arrange(int start, int end, int stepSize)
{
    std::vector<int> vector;
    for (int i{start}; i<end; i+=stepSize)
    {
        vector.push_back(start + i);
    }
    return vector;
}
void saveTimePositionData (std::vector<int> &time, std::vector<double> &x, std::vector<double> &y, std::string &fileName)
{
    std::ofstream file_;
    file_.open(fileName,std::ios::app);
    for (int i{0};i<time.size()-1;i++)
    {
        file_ << time[i] << ",";
    }
    file_ << time[time.size()-1] << "\n";
    for (int i{0};i<x.size()-1;i++)
    {
        file_ << x[i] << ",";
    }
    file_ << x[x.size()-1] << "\n";
    for (int i{0};i<y.size()-1;i++)
    {
        file_ << y[i] << ",";
    }
    file_ << y[y.size()-1] << "\n";
    file_.close();
}
std::vector<std::vector<double>> matrixBoundaryValue (std::vector<std::vector<double>> &matrix, double boundaryValue)
{//changes the boundary of the matrix to the specified value
    int M{0};
    M = matrix.size();
    int N{0};
    N = matrix[0].size();

    std::fill_n(matrix[0].begin(),N,boundaryValue);
    std::fill_n(matrix[M-1].begin(),N,boundaryValue);
    for (int i{1}; i<(M-1);i++)
    {
        matrix[i][0] = boundaryValue;
        matrix[i][N-1] = boundaryValue;
    }
    return matrix;
}
void saveVelocityData (std::vector<std::vector<double>> &u, std::string &fileName)
{
    std::ofstream file_;
    file_.open(fileName,std::ios::app);
    for (int i{0};i<u.size();i++)
    {
        for (int j{0};j<u[0].size()-1;j++)
        {
            file_ << u[i][j] << ",";
        }
        file_ << u[i][u[0].size()-1] << "\n";
    }
    file_.close();
}
