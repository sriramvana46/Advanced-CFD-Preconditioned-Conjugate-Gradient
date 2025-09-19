#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>
#define epsilon 1e-6

using namespace std;

void initialize_T(double T[17000])
{
    for(int i = 1; i <= 128*128; i++)
    {
        T[i] = 0.0;
    }
}

void vector_B(double b[17000])
{
    for(int i = 0; i <= 128*127; i++)
    {
        b[i] = 0.0;
    }
    for(int i = 128*127+1; i <= 128*128; i++)
    {
        b[i] = 2.0;
    }
}

double compute_dot_product(double A[17000], double B[17000])
{
    double sum = 0.0;
    for(int i = 1; i <= 128*128; i++)
    {
        sum = sum + A[i]*B[i];
    }
    return sum;
}

void compute_Ad(double Ad[17000], double Ap[17000], double Ae[17000], double Aw[17000], double An[17000], double As[17000], double d[17000])
{
    Ad[1] = Ap[1]*d[1] + Ae[1]*d[2] + An[1]*d[129];

    for(int i = 2; i<=128; i++)
    {
        Ad[i] = Aw[i]*d[i-1] + Ap[i]*d[i] + Ae[i]*d[i+1] + An[i]*d[i+128];
    }

    for(int i = 129; i <= 128*127; i++)
    {
        Ad[i] = As[i]*d[i-128] + Aw[i]*d[i-1] + Ap[i]*d[i] + Ae[i]*d[i+1] + An[i]*d[i+128];
    }

    for(int i = 128*127+1; i <= 128*128-1; i++)
    {
        Ad[i] = As[i]*d[i-128] + Aw[i]*d[i-1] + Ap[i]*d[i] + Ae[i]*d[i+1];
    }

    Ad[128*128] = As[128*128]*d[128*128-128] + Aw[128*128]*d[128*128-1] + Ap[128*128]*d[128*128];
}

void compute_residual(double res[17000], double b[17000], double AT[17000])
{
    for(int i = 1; i <= 128*128; i++)
    {
        res[i] = b[i] - AT[i];
    }
}
// Creating the Main diagonal of the A matrix
void vector_Ap(double Ap[17000])
{
    for(int i = 0; i < 128*128; i++)
    {
        Ap[i] = 4.0; // Initializing whole Array with 4.0
    }
    for(int i = 1; i<=126; i++)
    {
        Ap[i+1] = 5.0; // Corresponding to Bottom Boundary cells
        Ap[128*i+1] = 5.0; //Corresponding to the left boundary  cells
        Ap[128*(i+1)] = 5.0; // corresponding to the right boundary cells
        Ap[128*127 + i+1] = 5.0; // corresponding to the top boundary cells
    }
    Ap[1] = 6.0; // corresponding to bottom left corner cell
    Ap[128] = 6.0; // corresponding to bottom right corner cell
    Ap[128*127+1] = 6.0; // Corresponding to top left corner cell
    Ap[128*128] = 6.0; //corresponding to top right corner cell

}

void vector_Ae_Aw_An_As(double Ae[17000], double Aw[17000], double An[17000], double As[17000])
{
    for(int i = 0;i <= 17000; i++)
    {
        Ae[i] = -1.0;
        Aw[i] = -1.0;
        An[i] = -1.0;
        As[i] = -1.0;
    }
    for(int i = 1; i <= 126; i++)
    {
        Aw[1+128*i] = 0.0;
        Ae[128*(i+1)] = 0.0;
    }
}

void vector_Lw(double Lw[17000], double As[17000])
{
    for(int i = 0; i<=17000; i++)
    {
        Lw[i] = 0.0;
    }
    for(int i = 129; i <= 128*128; i++)
    {
        Lw[i] = As[i];
    }
}

void vector_Ls(double Ls[17000], double Aw[17000])
{
    for(int i = 0; i <= 17000; i++)
    {
        Ls[i] = 0.0;
    }
    for(int i = 2; i <= 128*128; i++)
    {
        Ls[i] = Aw[i];
    }
}

void compute_vector_Lp_Un_Ue(double Ap[17000], double Ae[17000], double An[17000], double Ue[17000], double Un[17000], double Lp[17000], double Ls[17000], double Lw[17000])
{
    Lp[1] = Ap[1];
    Un[1] = Ae[1]/Lp[1];
    Ue[1] = An[1]/Lp[1];

    for(int i = 2; i <=128; i++)
    {
        Lp[i] = Ap[i] - Ls[i]*Un[i-1];
        Un[i] = Ae[i]/Lp[i];
        Ue[i] = An[i]/Lp[i];
    }
    for(int i = 129; i<=128*127; i++)
    {
        Lp[i] = Ap[i] - Ls[i]*Un[i-1] - Lw[i]*Ue[i-128] ;
        Un[i] = Ae[i]/Lp[i];
        Ue[i] = An[i]/Lp[i];
    }
    for(int i = 128*127 + 1; i<=128*128-1; i++)
    {
        Lp[i] = Ap[i] - Ls[i]*Un[i-1] - Lw[i]*Ue[i-128];
        Un[i] = Ae[i]/Lp[i];
        Ue[i] = 0.0;
    }
    Lp[128*128] = Ap[128*128] - Ls[128*128]*Un[128*128-1] - Lw[128*128]*Ue[128*127];
    Un[128*128] = 0.0;
    Ue[128*128] = 0.0;
}

void construct_R(double R[17000], double res[17000], double Lp[17000], double Ls[17000], double Lw[17000])
{
    R[1] = res[1]/Lp[1];
    for(int i = 2; i <=128; i++)
    {
        R[i] = (res[i] - Ls[i]*R[i-1])/Lp[i];
    }
    for(int i = 129; i<=128*128; i++)
    {
        R[i] = (res[i] - Ls[i]*R[i-1] - Lw[i]*R[i-128])/Lp[i];
    }
}

void compute_searchdir(double del[17000], double R[17000], double Un[17000], double Ue[17000])
{
    del[128*128] = R[128*128];
    for(int i = 128*128-1; i >=128*127+1; i--)
    {
        del[i] = R[i] - Un[i]*del[i+1];
    }
    for(int i = 128*127; i>=1; i--)
    {
        del[i] = R[i] - Un[i]*del[i+1] - Ue[i]*del[i+128];
    }
}

double compute_alpha(double num, double denom)
{
   double alpha = num/denom;
   return alpha;
}

void update_Temp(double T[17000], double d[17000], double alpha)
{
    for(int i = 1; i <= 128*128; i++)
    {
        T[i] = T[i] + alpha*d[i];
    }
}

void update_residual(double res_upd[17000], double res[17000], double Ad[17000], double alpha)
{
    for(int i = 1; i <= 128*128; i++)
    {
        res_upd[i] = res[i] - alpha*Ad[i];
    }
}

double compute_beta(double num, double denom)
{
    double beta = num/denom;
    return beta;
}

void update_search_dir(double Pin_res[17000], double d[17000], double beta)
{
    for(int i = 1; i <= 128*128; i++)
    {
        d[i] = Pin_res[i] + beta*d[i];
    }
}

void update_vectors(double A[17000], double B[17000])
{
    for(int i = 1; i <= 128*128; i++)
    {
        A[i] = B[i];
    }
}

int main()
{
    double Ap[17000], Ae[17000], Aw[17000], An[17000], As[17000];
    double Lp[17000], Ls[17000], Lw[17000];
    double Up[17000], Un[17000], Ue[17000];
    double T[17000],b[17000], R[17000], res[17000], res_upd[17000], d[17000];
    double AT[17000], Ad[17000];
    double Pin_res[17000], Pin_res_upd[17000];
    double num, denom, alpha, beta;
    double restres = 1.0;
    int iter = 0;
    initialize_T(T);
    vector_B(b);
    vector_Ap(Ap);
    vector_Ae_Aw_An_As(Ae, Aw, An, As);
    vector_Lw(Lw, As);
    vector_Ls(Ls, Aw);
    compute_vector_Lp_Un_Ue(Ap, Ae, An, Ue, Un, Lp, Ls, Lw);
    compute_Ad(AT, Ap, Ae, Aw, An, As, T);
    compute_residual(res, b, AT);
    construct_R(R, res, Lp, Ls, Lw);
    compute_searchdir(d, R, Un, Ue);
    ofstream outfile("ILU_Results.txt");

    do
    {
        compute_Ad(Ad, Ap, Ae, Aw, An, As, d);
        construct_R(R, res, Lp, Ls, Lw);
        compute_searchdir(Pin_res, R, Un, Ue);
        num = compute_dot_product(res, Pin_res);
        denom = compute_dot_product(d, Ad);
        alpha = compute_alpha(num, denom);
        update_Temp(T, d, alpha);
        if(iter % 50 == 0)
        {
            compute_Ad(AT, Ap, Ae, Aw, An, As, T);
            compute_residual(res_upd,b,AT);
        }
        else
        {
            update_residual(res_upd, res, Ad, alpha);
        }
        construct_R(R, res_upd, Lp, Ls, Lw);
        compute_searchdir(Pin_res_upd, R, Un, Ue);
        num = compute_dot_product(res_upd, Pin_res_upd);
        denom = compute_dot_product(res, Pin_res);
        beta = compute_beta(num, denom);
        update_search_dir(Pin_res_upd, d, beta);
        update_vectors(res, res_upd);
        iter++;
        restres = compute_dot_product(res, res);
        outfile << "Iterations = " << iter << endl;
        outfile << "Res-norm = " << sqrt(restres) << endl;

    } while (sqrt(restres) > epsilon);

    outfile << "\nTemperature at the vertical centre line of the domain: " << endl;
    for(int i = 0; i < 128; i++)
    {
        outfile << "T[" << 65+128*i << "]\t=\t " << T[65+128*i] << endl;
    }
    outfile << "\nResidual 2-Norm = " << sqrt(restres) << endl;
    outfile << "Iterations = " << iter << endl;
    cout << "Results printed in ILU_Results.txt" << endl;
    return 0;
}
