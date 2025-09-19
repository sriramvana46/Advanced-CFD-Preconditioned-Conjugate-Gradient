// Preconditioning using Jacobi
#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream> 
#define epsilon 1e-6
using namespace std;

void initialize(double T[17000])
{
    //Initializing interior cell centre Temperature values to be zero
    for(int i = 2; i < 130; i++)
    {
        int a = 130*i;
        for(int j = a-129; j < a-1; j++)
        {
            T[j] = 0.0;
        }
    }
}
void vector_b(double b[17000])
{
    for(int i = 3; i < 129; i++)
    {
        int a = 130*i;
        for(int j = a-128; j < a-2; j++)
        {
           b[j] = 0.0; //interior cells excluding boundary cells
        }
    }
    for(int i = 2; i < 128; i++)
    {
        b[130+i] = 0.0;//bottom boundary cells
    }
    for(int i = (130*128)+2; i < (130*129)-2; i++)
    {
        b[i] = 2.0; //Top boundary cells
    }
    for(int i = 2; i < 128; i++)
    {
        b[1+130*i] = 0.0; //Left boundary cells
    }
    for(int i = 2; i < 128; i++)
    {
        b[128+130*i] = 0.0; // Right boundary cells
    }
    b[131] = 0.0; // Bottom left corner cell
    b[258] = 0.0; //Bottom right corner cell
    b[1+130*128] = 2.0; //Top left corner cell
    b[128+130*128] = 2.0; //Top right corner cell

}

void compute_residual(double b[17000], double T[17000], double res[17000])
{

    for(int i = 3; i < 129; i++)
    {
        int a = 130*i;
        for(int j = a-128; j < a-2; j++)
        {
            res[j] = b[j] + T[j+1] + T[j-1] + T[j+130] + T[j-130] - 4*T[j]; // Residual calculation for interion cells excluding boundary cells
        }
    }
    for(int j = 132; j < 258; j++)
    {
        res[j] = b[j] + T[j-1] + T[j+1] + T[j+130] - 5*T[j]; //Residual calculation for bottom boundary cells
    }
    for(int i = (130*128)+2; i < (130*129)-2; i++)
    {
        res[i] = b[i] + T[i-1] + T[i+1] + T[i-130] - 5*T[i]; //Residual calculation for Top boundary cells
    }
    for(int i = 2; i < 128; i++)
    {
        res[1+130*i] = b[1+130*i] + T[2+130*i] + T[131+130*i] + T[-130+1+130*i] - 5*T[1+130*i]; //Residual calculation for Left boundary cells
    }
    for(int i = 2; i < 128; i++)
    {
        res[130*(i+1)-2] = b[130*(i+1)-2] + T[130*(i+1)-3] + T[130*(i+1)-132] + T[128+130*(i+1)] - 5*T[130*(i+1)-2]; //Residual calculation for Right boundary cells
    }
    res[131] = b[131] + T[132] + T[261] - 6*T[131]; // Bottom left corner cell
    res[258] = b[258] + T[257] + T[388] - 6*T[258]; //Bottom Right corner cell
    res[1+130*128] = b[1+130*128] + T[2+130*128] + T[1+130*127] - 6*T[1+130*128]; // Top left corner cell
    res[(130*129)-2] = b[(130*129)-2] + T[(130*129)-3] + T[(130*128)-2] - 6*T[(130*129)-2]; // Top right corner cell
}

void compute_do(double d[17000], double res[17000])
 {
    for(int i = 2; i < 130; i++)
    {
        int a = 130*i;
        for(int j = a-129; j < a-1; j++)
        {
            d[j] = res[j];
        }
    }
 }


void compute_Ad(double Ad[17000], double d[17000])
{
    for(int i = 3; i < 129; i++)
    {
        int a = 130*i;
        for(int j = a-128; j < a-2; j++)
        {
            Ad[j] = -d[j+1] - d[j-1] - d[j+130] - d[j-130] + 4*d[j]; //interion cells excluding boundary cells
        }
    }
    for(int j = 132; j < 258; j++)
    {
        Ad[j] = -d[j-1] - d[j+1] - d[j+130] + 5*d[j]; //bottom boundary cells
    }
    for(int i = (130*128)+2; i < (130*129)-2; i++)
    {
        Ad[i] = -d[i-1] - d[i+1] - d[i-130] + 5*d[i]; //Top boundary cells
    }
    for(int i = 2; i < 128; i++)
    {
        Ad[1+130*i] = -d[2+130*i] - d[131+130*i] - d[-130+1+130*i] + 5*d[1+130*i]; // Left boundary cells
    }
    for(int i = 2; i < 128; i++)
    {
        Ad[-2+130*(i+1)] = -d[-3+130*(i+1)] - d[-132+130*(i+1)] - d[128+130*(i+1)] + 5*d[-2+130*(i+1)]; //Right boundary cells
    }
    Ad[131] = -d[132] - d[261] + 6*d[131]; // Bottom left corner cell
    Ad[258] = -d[257] - d[388] + 6*d[258]; //Bottom Right corner cell
    Ad[1+130*128] = -d[2+130*128] - d[1+130*127] + 6*d[1+130*128]; // Top left corner cell
    Ad[(130*129)-2] = -d[(130*129)-3] - d[(130*128)-2] + 6*d[(130*129)-2]; // Top right corner cell
}
void compute_Pinverse_res(double Ad[17000], double d[17000])
{
    for(int i = 3; i < 129; i++)
    {
        int a = 130*i;
        for(int j = a-128; j < a-2; j++)
        {
            Ad[j] = 0.25*d[j]; //interion cells excluding boundary cells
        }
    }
    for(int j = 132; j < 258; j++)
    {
        Ad[j] = 0.2*d[j]; //bottom boundary cells
    }
    for(int i = (130*128)+2; i < (130*129)-2; i++)
    {
        Ad[i] = 0.2*d[i]; //Top boundary cells
    }
    for(int i = 2; i < 128; i++)
    {
        Ad[1+130*i] = 0.2*d[1+130*i]; // Left boundary cells
    }
    for(int i = 2; i < 128; i++)
    {
        Ad[-2+130*(i+1)] = 0.2*d[-2+130*(i+1)]; //Right boundary cells
    }
    Ad[131] = d[131]/6.0; // Bottom left corner cell
    Ad[258] = d[258]/6.0; //Bottom Right corner cell
    Ad[1+130*128] = d[1+130*128]/6.0; // Top left corner cell
    Ad[(130*129)-2] = d[(130*129)-2]/6.0; // Top right corner cell
}

double compute_vectormultiplication(double dt[17000], double Ad[17000])
{
    double sum = 0.0;
    for(int i = 2; i < 130; i++)
    {
        int a = 130*i;
        for(int j = a-129; j < a-1; j++)
        {
            sum+= dt[j]*Ad[j];
        }
    }
    return sum;
}

double compute_alpha(double Pinvres[17000],double res[17000], double d[17000], double Ad[17000])
{
    double alpha, dtAd, restPinvres;
    dtAd = compute_vectormultiplication(d,Ad);
    restPinvres = compute_vectormultiplication(res,Pinvres);
    alpha = (restPinvres)/(dtAd);
    return alpha;
}

void update_Temperature(double T[17000], double &alpha, double d[17000])
{
    for(int i = 2; i < 130; i++)
    {
        int a = 130*i;
        for(int j = a-129; j < a-1; j++)
        {
            T[j] = T[j] + alpha*d[j];
        }
    }
}

void update_residual(double res_upd[17000], double res[17000], double &alpha, double Ad[17000])
{
    for(int i = 2; i < 130; i++)
    {
        int a = 130*i;
        for(int j = a-129; j < a-1; j++)
        {
            res_upd[j] = res[j] - alpha*Ad[j];
        }
    }
}

double compute_beta(double res_upd[17000], double res[17000])
{
    double a,b,beta;
    double Pinvres[17000], Pinvres_upd[17000];
    compute_Pinverse_res(Pinvres, res);
    compute_Pinverse_res(Pinvres_upd, res_upd);
    a = compute_vectormultiplication(res,Pinvres);
    b = compute_vectormultiplication(res_upd,Pinvres_upd);
    beta = b/a;
    return beta;
}

void update_search_dir(double d[17000], double res_upd[17000], double &beta)
{
    for(int i = 2; i < 130; i++)
    {
        int a = 130*i;
        for(int j = a-129; j < a-1; j++)
        {
            d[j] = res_upd[j] + beta*d[j];
        }
    }
}


int main()
{
    double b[17000], Pinvres[17000];
    double T[17000], res[17000], res_upd[17000], res_o[17000], d[17000], Ad[17000], AT[17000];
    double  alpha, beta, Restpinvres;
    double restres = 1.0;
    int iter = 0;
    initialize(T);
    vector_b(b);
    compute_residual( b, T, res_o);
    compute_Pinverse_res(Pinvres, res_o);
    compute_do(d,Pinvres);
    compute_do(res,res_o);
    ofstream outfile("Preconditioning_Results.txt");
    do
    {
        compute_Ad(Ad,d);
        alpha = compute_alpha( Pinvres,res, d, Ad);
        update_Temperature(T,alpha,d);
        if(iter % 50 == 0)
        {
            compute_Ad(AT,T);
            compute_residual(b,T,res_upd);
        }
        else
        {
            update_residual(res_upd, res, alpha, Ad);
        }
        beta = compute_beta(res_upd,res);
        compute_do(res,res_upd);
        compute_Pinverse_res(Pinvres, res);
        update_search_dir(d, Pinvres,beta);
        iter++;
        restres = compute_vectormultiplication(res,res);
        outfile << "Iterations = " << iter << endl;
        outfile << "Res-norm = " << sqrt(restres) << endl;
        
        
        
    } while (sqrt(restres) > epsilon);
    outfile << "\nTemperature at the vertical centre line of the domain: " << endl;
    for(int i = 0; i<128; i++)
    {
        int a  =130*i;
        outfile << "T[" << 194 + a << "] = " << T[194+a] << endl;
    }
    outfile << "\nResidual 2-Norm = " << sqrt(restres) << endl;
    outfile << "Iterations = " << iter << endl;
    cout << "Results printed in Preconditioning_Results.txt" << endl;
return 0;
}