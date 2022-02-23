#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>

std::vector<double> linspace(double start, double end, double steps);
void displayVector (std::vector<double> &v);
void safeDataToTextFile (std::vector<double> &x, std::vector<double> &y, std::vector<std::vector<double>> &u, std::vector<std::vector<double>> &v, std::string &fileName);
std::vector<std::vector<double>> matrixBoundaryValue (std::vector<std::vector<double>> &matrix, double boundaryValue);

int main() {
    //variables
    int nx{81};
    int ny{81};
    int nt{100};
    double c{1};
    double dx{2.0/(nx-1.0)};
    double dy{2.0/(ny-1.0)};
    double sigma{0.2};
    double dt{sigma*dx};

    //vector initialization
    std::vector<double> x = linspace(0,2,nx);
    std::vector<double> y = linspace(0,2,ny);

    std::vector<std::vector<double>> u (ny,std::vector<double>(nx,1));
    std::vector<std::vector<double>> un (ny,std::vector<double>(nx,1));
    std::vector<std::vector<double>> v (ny,std::vector<double>(nx,1));
    std::vector<std::vector<double>> vn (ny,std::vector<double>(nx,1));


    //data input
    for (int i{static_cast<int>(0.5/dy)};i<(round(1.0/dy + 1));i++)
    {
        std::fill_n(u[i].begin()+round(0.5/dx),round(1.0/dx + 1) - round(0.5/dx), 2.0);
        std::fill_n(v[i].begin()+round(0.5/dx),round(1.0/dx + 1) - round(0.5/dx), 2.0);
    }

    for (int n{0}; n<(nt+1); n++)
    {
        for (int i{0}; i<y.size(); i++)
        {
            un[i] = u[i];
            vn[i] = v[i];
        }
        for (int j{1}; j<y.size(); j++)
        {
            for (int i{1}; i<x.size(); i++)
            {
                u[j][i] = un[j][i] - un[j][i]*(dt/dx)*(un[j][i]-un[j][i-1]) - vn[j][i]*(dt/dy)*(un[j][i]-un[j-1][i]);
                v[j][i] = vn[j][i] - un[j][i]*(dt/dx)*(vn[j][i]-vn[j][i-1]) - vn[j][i]*(dt/dy)*(vn[j][i]-vn[j-1][i]);
            }
        }
        u = matrixBoundaryValue(u,1.0);
        v = matrixBoundaryValue(v,1.0);
    }

    std::string fileName{"6_ConvectionData.txt"};
    safeDataToTextFile(x,y,u,v,fileName);

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

void safeDataToTextFile (std::vector<double> &t, std::vector<double> &x, std::vector<double> &y, std::vector<std::vector<double>> &u, std::vector<std::vector<double>> &v, std::string &fileName)
{
    std::ofstream file_;
    file_.open(fileName);
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
    for (int i{0};i<y.size();i++)
    {
        for (int j{0};j<x.size()-1;j++)
        {
            file_ << u[i][j] << ",";
        }
        file_ << u[i][x.size()-1] << "\n";
    }
    for (int i{0};i<y.size();i++)
    {
        for (int j{0};j<x.size()-1;j++)
        {
            file_ << v[i][j] << ",";
        }
        file_ << v[i][x.size()-1] << "\n";
    }
    file_.close();
}
