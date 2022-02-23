#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>
#include <string>
#include <fstream>

std::vector<double> linspace(double start, double end, double steps);
std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> applyVelocityBoundaryConditions (std::vector<double> &x, std::vector<double> &y, std::vector<std::vector<double>> &u, std::vector<std::vector<double>> &v, double &convergenceValue );
std::vector<std::vector<double>> subtractMatrixABS(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &B);
double sumMatrix(std::vector<std::vector<double>> M);
std::vector<std::vector<double>> pressureIterator (std::vector<double> &x, std::vector<double> &y, std::vector<std::vector<double>> &u, std::vector<std::vector<double>> &v, std::vector<std::vector<double>> &p, double &pressureConvergenceValue, double &numberConvergenceValue, double &dx, double &dy, double &dt, double &rho);
void safeDataToTextFile (std::vector<double> &x, std::vector<double> &y, std::vector<std::vector<double>> &u, std::vector<std::vector<double>> &v, std::vector<std::vector<double>> &p, std::string &fileName);

int main() {
    //Input variables
    int nx{41};
    int ny{41};
    int nt{500};
    double dt{0.001};
    double pressureConvergenceValue {0.1};
    double numberConvergenceValue {pow(10,-5)};

    double rho{1};
    double nu{0.1};

    double dx{2.0/(nx-1.0)};
    double dy{2.0/(ny-1.0)};

    //Vector and Matrix initialization
    std::vector<double> x = linspace(0,2,nx);
    std::vector<double> y = linspace(0,2,ny);
    std::vector<std::vector<double>> u (ny,std::vector<double>(nx,0.0));
    std::vector<std::vector<double>> un = u;
    std::vector<std::vector<double>> v (ny,std::vector<double>(nx,0.0));
    std::vector<std::vector<double>> vn = v;
    std::vector<std::vector<double>> p (ny,std::vector<double>(nx,0.0));

    //Calculation
    for (int n{0}; n<nt; n++)
    {
        pressureIterator(x,y,u,v,p,pressureConvergenceValue,numberConvergenceValue,dx,dy,dt,rho);
        for (int j{0}; j<y.size();j++)
        {
            un[j] = u[j];
            vn[j] = v[j];
        }
        for (int j{1}; j<y.size()-1;j++)
        {
            for ( int i{1}; i<x.size()-1;i++)
            {
                u[j][i] = un[j][i] - un[j][i] * dt/dx * (un[j][i] - un[j][i-1]) - vn[j][i] * dt/dy * (un[j][i] - un[j-1][i]) - dt/(rho * 2 * dx) * (p[j][i+1] - p[j][i-1]) + nu * (dt * pow(dx,-2) * (un[j][i+1] - 2 * un[j][i] + un[j][i-1]) + dt * pow(dy,-2) * (un[j+1][i] - 2 * un[j][i] + un[j-1][i]));
                v[j][i] = vn[j][i] - un[j][i] * dt/dx * (vn[j][i] - vn[j][i-1]) - vn[j][i] * dt/dy * (vn[j][i] - vn[j-1][i]) - dt/(rho * 2 * dy) * (p[j+1][i] - p[j-1][i]) + nu * (dt * pow(dx,-2) * (vn[j][i+1] - 2 * vn[j][i] + vn[j][i-1]) + dt * pow(dy,-2) * (vn[j+1][i] - 2 * vn[j][i] + vn[j-1][i]));
            }
        }
        applyVelocityBoundaryConditions(x,y,u,v,numberConvergenceValue);
    }

    std::string fileName{"11_CavityFlow.txt"};
    safeDataToTextFile(x,y,u,v,p,fileName);

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
std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> applyVelocityBoundaryConditions (std::vector<double> &x, std::vector<double> &y, std::vector<std::vector<double>> &u, std::vector<std::vector<double>> &v, double &convergenceValue )
{
    for(int j{0}; j<y.size();j++)
    {
        for(int i{0}; i<(x.size());i++)
        {
            if ( std::abs(y[j] - 0.0) < convergenceValue )
            {
                u[j][i] = 0.0;
                v[j][i] = 0.0;
            }
            else if ( std::abs(x[i] - 0.0) < convergenceValue )
            {
                u[j][i] = 0.0;
                v[j][i] = 0.0;
            }
            else if ( std::abs(x[i] - 2.0) < convergenceValue )
            {
                u[j][i] = 0.0;
                v[j][i] = 0.0;
            }
            else if( std::abs(y[j] - 2.0) < convergenceValue )
            {
                u[j][i] = 1.0;
                v[j][i] = 0.0;
            }
        }
    }
    auto tuple = std::make_tuple(u,v);
    return tuple;
}
std::vector<std::vector<double>> applyPressureBoundaryConditions (std::vector<double> &x, std::vector<double> &y, std::vector<std::vector<double>> &p, double &convergenceValue )
{
    for(int j{0}; j<y.size();j++)
    {
        for(int i{0}; i<(x.size());i++)
        {
            if ( std::abs(x[i] - 2.0) < convergenceValue )
            {
                p[j][i] = p[j][i-1];
            }
            if ( std::abs(y[j] - 0.0) < convergenceValue )
            {
                p[j][i] = p[j+1][i];
            }
            if ( std::abs(x[i] - 0.0) < convergenceValue )
            {
                p[j][i] = p[j][i+1];
            }
            if( std::abs(y[j] - 2.0) < convergenceValue )
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
std::vector<std::vector<double>> pressureIterator (std::vector<double> &x, std::vector<double> &y, std::vector<std::vector<double>> &u, std::vector<std::vector<double>> &v, std::vector<std::vector<double>> &p, double &pressureConvergenceValue, double &numberConvergenceValue, double &dx, double &dy, double &dt, double &rho)
{
    double pressureConvergence{1.0};
    std::vector<std::vector<double>> pn = p;
    while (pressureConvergence>pressureConvergenceValue)
    {
        for (int j{0}; j<p.size(); j++)
        {
            pn[j] = p[j];
        }
        for (int j{1}; j<p.size()-1; j++)
        {
            for (int i{1}; i<p[0].size()-1; i++)
            {
                p[j][i] = ((pn[j][i+1] + pn[j][i-1]) * pow(dy,2) + (pn[j+1][i] + pn[j-1][i]) * pow(dx,2))/(2*(pow(dx,2) + pow(dy,2))) - (rho * pow(dx,2) * pow(dy,2))/(2*(pow(dx,2) + pow(dy,2))) * (pow(dt,-1) * ((u[j][i+1] - u[j][i-1])/(2*dx) + (v[j+1][i] - v[j-1][i])/(2*dy)) - pow(((u[j][i+1] - u[j][i-1])/(2*dx)),2) - 2 * (((u[j+1][i] - u[j-1][i])/(2*dy)) * ((v[j][i+1] - v[j][i-1])/(2*dx))) - pow(((v[j+1][i] - v[j-1][i])/(2*dy)),2));
            }
        }
        applyPressureBoundaryConditions(x,y,p,numberConvergenceValue);

        pressureConvergence = sumMatrix(subtractMatrixABS(p,pn))/ sumMatrix(pn);
    }
    return p;
}
void safeDataToTextFile (std::vector<double> &x, std::vector<double> &y, std::vector<std::vector<double>> &u, std::vector<std::vector<double>> &v, std::vector<std::vector<double>> &p, std::string &fileName)
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
