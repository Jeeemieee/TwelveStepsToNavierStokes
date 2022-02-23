#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

void displayVector (std::vector<double> &v);
void displayVector (std::vector<int> &v);
std::vector<double> linspace(double start, double end, double steps);
std::vector<std::vector<double>> applyBoundaryConditions(std::vector<std::vector<double>> &p, std::vector<double> &x, std::vector<double> &y, double &epsilon);
void displayMatrix (std::vector<std::vector<double>> &M);
std::vector<std::vector<double>> subtractMatrixABS(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &B);
double sumMatrix(std::vector<std::vector<double>> M);
void safeDataToTextFile (std::vector<double> &x, std::vector<double> &y, std::vector<std::vector<double>> &p, std::string &fileName);


int main() {
    //simulation variables
    int nx{50};
    int ny{50};
    double dx{2.0/(nx-1.0)};
    double dy{2.0/(ny-1.0)};

    //vector initialization
    std::vector<double> x = linspace(0,2,nx);
    std::vector<double> y = linspace(0,1,ny);
    std::vector<std::vector<double>> p (ny,std::vector<double>(nx,0.0));
    std::vector<std::vector<double>> pn (ny,std::vector<double>(nx,0.0));
    std::vector<std::vector<double>> b (ny,std::vector<double>(nx,0.0));

    //Boundary conditions
    double epsilon{pow(10,-5)};
    applyBoundaryConditions(p,x,y,epsilon);
    b[int(ny/4)][int(nx/4)] = 100;
    b[int(3*ny/4)][int(3*nx/4)] = -100;

    //Convergence criteria
    double convergenceValue{1.0};

    while (convergenceValue>=epsilon)
    {
        for (int i{0}; i<y.size();i++)
        {
            pn[i] = p[i];
        }
        for (int j{1}; j<(y.size()-1);j++)
        {
            for (int i{1}; i<(x.size()-1);i++)
            {
                p[j][i] = (pow(dy,2) * (p[j][i+1] + p[j][i-1]) + pow(dx,2) * (p[j+1][i] + p[j-1][i]) - b[j][i] * pow(dx,2) * pow(dy,2) )/(2.0*(pow(dy,2) + pow(dx,2)));

            }
        }
        applyBoundaryConditions(p,x,y,epsilon);

        convergenceValue = sumMatrix(subtractMatrixABS(p,pn))/ sumMatrix(pn);
    }

    std::string fileName{"10_PoissonEquationData.txt"};
    safeDataToTextFile(x,y,p,fileName);

    return 0;
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
void displayMatrix (std::vector<std::vector<double>> &M)
{
    for (int j{0}; j<M.size();j++)
    {
        for(int i{0}; i < M[0].size();i++)
        {
            std::cout << M[j][i] << " ";
        }
        std::cout << std::endl;
    }
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
std::vector<std::vector<double>> applyBoundaryConditions(std::vector<std::vector<double>> &p, std::vector<double> &x, std::vector<double> &y, double &epsilon)
{
    for (int i{0}; i<x.size();i++)
    {
        if (std::abs(x[i]-0.0)<=epsilon)
        {// p = 0 @ x = 0
            for(int j{0};j<y.size();j++)
            {
                p[j][i] = 0.0;
            }
        }
        else if (std::abs(x[i]-2.0)<=epsilon)
        {// p = 0 @ x = 2
            for(int j{0};j<y.size();j++)
            {
                p[j][i] = 0.0;
            }
        }
    }
    for (int j{0}; j<y.size();j++)
    {
        if (std::abs(y[j]-0.0)<=epsilon)
        {// p = 0 @ y = 0
            for(int i{0};i<x.size();i++)
            {
                p[j][i] = 0.0;
            }
        }
        else if (std::abs(y[j]-1.0)<=epsilon)
        {// p = 0 @ y = 1
            for(int i{0};i<x.size();i++)
            {
                p[j][i] = 0.0;
            }
        }
    }
    return p;
}
std::vector<std::vector<double>> subtractMatrixABS(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &B)
{
    std::vector<std::vector<double>> result (A.size(),std::vector<double>(A[0].size(),0.0));
    for (int j{0}; j<A.size();j++)
    {
        for (int i{0}; i<A[0].size();i++)
        {
            result[j][i] = std::abs(A[j][i]) - std::abs(B[j][i]);
        }
    }
    return result;
}
double sumMatrix(std::vector<std::vector<double>> M)
{
    double sum;
    for (int j{0}; j<M.size();j++)
    {
        for (int i{0}; i<M[0].size();i++)
        {
            sum += M[j][i];
        }
    }
    return sum;
}
void safeDataToTextFile (std::vector<double> &x, std::vector<double> &y, std::vector<std::vector<double>> &p, std::string &fileName)
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
            file_ << p[i][j] << ",";
        }
        file_ << p[i][x.size()-1] << "\n";
    }
    file_.close();
}
