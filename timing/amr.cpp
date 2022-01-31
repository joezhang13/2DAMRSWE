#include <mpi.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <string>
#include <cstring>
#include <vector>
#include <numeric>
#include "flux.h"
using namespace std;

const double g = 9.8;
const double rho = 1000;
const double L = 10;
const double CFL = 0.45;
const double GHcr = 0.1;
const double normal[4][2] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};

double initial(double x, double y)
{
    double h;
    h = 6 + 3*exp(-8*(x-5)*(x-5) - 5*(y-5)*(y-5));

    return h;
}

int main(int argc, char** argv) {
    // MPI initialization
    int P;
    int ID;
    MPI_Init(&argc , &argv);
    MPI_Comm_size(MPI_COMM_WORLD , &P);
    MPI_Comm_rank(MPI_COMM_WORLD , &ID);
    int N = 400;    //number of initial grids in each dimension
    double Tend = 0.5;    //running time
    int nout = 10;    //number of output
    //int nt = 1000;    //number of time steps
    int M = ceil(double(N)/double(P));
    int m;    //number of strips of initial grids in each processor
    if (M*(P - 1) > N - 1)
        M--;
    if (ID < P - 1)
        m = M;
    else
        m = N - M*ID; 

    // Mesh information
    int Ncell = N*N;    //number of coarse cells
    double cL = L/N;
    double area = cL*cL;
    int Ni = 2*N*(N-1);    //number of coarse interior edges
    int Nb = 4*N;    //number of coarse boundary edges
    int Ncellf = 4*Ncell;    //number of refined cells
    double cLf = cL/2;
    double areaf = cLf*cLf;
    int NiF = 2*(2*N)*(2*N-1);    //number of refined interior edges
    int NbF = 4*(2*N);    //number of refined boundary edges

    vector<vector<int> > C2C(Ncell, vector<int>(0));     //coarse cells to coarse cells
    vector<vector<int> > C2Cn(Ncell, vector<int>(0));    //normal vectors corresponding to C2C
    vector<vector<int> > C2F(Ncell, vector<int>(0));     //coarse cells to fine cells
    vector<vector<int> > C2Fn(Ncell, vector<int>(0));    //normal vectors corresponding to C2F
    vector<vector<int> > F2F(Ncellf, vector<int>(0));    //fine cells to fine cells
    vector<vector<int> > F2Fn(Ncellf, vector<int>(0));   //normal vectors corresponding to F2F
    vector<vector<int> > C2Bn(Ncell, vector<int>(0));    //normal vectors corresponding to coarse boundary cells
    vector<vector<int> > F2Bn(Ncellf, vector<int>(0));   //normal vectors corresponding to fine boundary cells
    vector<vector<int> > CF(Ncell, vector<int>(4));      //mapping from coarse to fine
    vector<int> FC(Ncellf);                              //mapping from fine to coarse

    for (int i = 0; i < Ncell; i++)
    {
        for (int j = 0; j < 4; j++)
            CF[i][j] = 4*i + j;
    }
    for (int i = 0; i < Ncell; i++)
    {
        for (int j = 0; j < 4; j++)
            FC[4*i+j] = i;
    }
    
    for (int i = 0; i < N; i++)        //vertical edges
    {
        for (int j = 0; j < N - 1; j++)
        {
            int iCell = i*N + j;
            C2C[iCell].push_back(iCell+1);
            C2Cn[iCell].push_back(0);
            C2F[iCell].push_back(CF[iCell+1][0]);
            C2Fn[iCell].push_back(0);
            C2F[iCell].push_back(CF[iCell+1][2]);
            C2Fn[iCell].push_back(0);
            C2F[iCell+1].push_back(CF[iCell][1]);
            C2Fn[iCell+1].push_back(1);
            C2F[iCell+1].push_back(CF[iCell][3]);
            C2Fn[iCell+1].push_back(1);
            F2F[CF[iCell][1]].push_back(CF[iCell+1][0]);
            F2Fn[CF[iCell][1]].push_back(0);
            F2F[CF[iCell][3]].push_back(CF[iCell+1][2]);
            F2Fn[CF[iCell][3]].push_back(0);
        }
    }
    for (int i = 0; i < N - 1; i++)        //horizontal edges
    {
        for (int j = 0; j < N; j++)
        {
            int iCell = i*N + j;
            C2C[iCell].push_back(iCell+N);
            C2Cn[iCell].push_back(2);
            C2F[iCell].push_back(CF[iCell+N][0]);
            C2Fn[iCell].push_back(2);
            C2F[iCell].push_back(CF[iCell+N][1]);
            C2Fn[iCell].push_back(2);
            C2F[iCell+N].push_back(CF[iCell][2]);
            C2Fn[iCell+N].push_back(3);
            C2F[iCell+N].push_back(CF[iCell][3]);
            C2Fn[iCell+N].push_back(3);
            F2F[CF[iCell][2]].push_back(CF[iCell+N][0]);
            F2Fn[CF[iCell][2]].push_back(2);
            F2F[CF[iCell][3]].push_back(CF[iCell+N][1]);
            F2Fn[CF[iCell][3]].push_back(2);
        }
    }
    for (int i = 0; i < Ncell; i++)       //edges inside the refined cells
    {
        F2F[CF[i][0]].push_back(CF[i][1]);
        F2Fn[CF[i][0]].push_back(0);
        F2F[CF[i][2]].push_back(CF[i][3]);
        F2Fn[CF[i][2]].push_back(0);
        F2F[CF[i][0]].push_back(CF[i][2]);
        F2Fn[CF[i][0]].push_back(2);
        F2F[CF[i][1]].push_back(CF[i][3]);
        F2Fn[CF[i][1]].push_back(2);
    }
    for (int i = 0; i < N; i++)            //boundary edges
    {
        int ibdown = i;
        int ibup = N*(N-1) + i;
        int ibleft = i*N;
        int ibright = N - 1 + i*N;
        C2Bn[ibdown].push_back(3);
        C2Bn[ibup].push_back(2);
        C2Bn[ibleft].push_back(1);
        C2Bn[ibright].push_back(0);
        F2Bn[CF[ibdown][0]].push_back(3);
        F2Bn[CF[ibdown][1]].push_back(3);
        F2Bn[CF[ibup][2]].push_back(2);
        F2Bn[CF[ibup][3]].push_back(2);
        F2Bn[CF[ibleft][0]].push_back(1);
        F2Bn[CF[ibleft][2]].push_back(1);
        F2Bn[CF[ibright][1]].push_back(0);
        F2Bn[CF[ibright][3]].push_back(0);
    }
    
    // Initial condition
    int Nc = m*N;
    int Nf = 4*Nc;
    int istart = ID*M*N;
    vector<int> I(Nc);
    vector<vector<double> > Uc(3, vector<double>(Nc));
    vector<vector<double> > Uf(3, vector<double>(Nf));
    vector<vector<double> > Rc(3, vector<double>(Nc));
    vector<vector<double> > Rf(3, vector<double>(Nf));
    double x;
    double y;
    double h;
    int iic;
    int iif;
    double smax = 0;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < N; j++)
        {
            iic = i*N + j;
            x = (j+0.5)*cL;
            y = (i+0.5+M*ID)*cL;
            h = initial(x,y);
            Uc[0][iic] = h;
            Uc[1][iic] = 0;
            Uc[2][iic] = 0;
            smax = max(smax, sqrt(g*h));
            for (int k = 0; k < 4; k++)
            {
                iif = CF[iic][k];
                Uf[0][iif] = h;
                Uf[1][iif] = 0;
                Uf[2][iif] = 0;
            }
        }
    }
    double dtl = CFL*cLf/smax;    //use CFL condition to determine the time step size
    double dt;
    MPI_Allreduce(&dtl, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    // FV solver
    int nt = ceil(int(Tend/dt));
    dt = double(Tend/nt);
    vector<double> hc(Nc);
    vector<vector<double> > gradh(2, vector<double>(Nc));    //refinement criteria (|grad h|)
    int igstart = istart;

    //time loop
    double t0;
    double tend;
    double t;
    if (ID == 0)
    {
        t0 = MPI_Wtime();
    }
    double tup0;
    double tupend;
    double tup = 0;
    double ti0;
    double tiend;
    double ti = 0;
    double tlb0;
    double tlbend;
    double tlb = 0;
    double ts = 0;
    int iout = 0;
    for (int n = 0; n < nt; n++)
    {
        if (ID == 0)
        {
            ti0 = MPI_Wtime();
        }
        int igap = istart - igstart;    //length of the ghost cells in the front (before repartition)
        vector<vector<int> > recvFrom(P, vector<int> (0));
        vector<vector<int> > recvStart(P, vector<int> (0));
        vector<vector<int> > recvLength(P, vector<int> (0));
        vector<int> sendTo(0);
        vector<int> sendStart(0);
        vector<int> sendLength(0);
        //obtain values of ghost cells from other processors
        vector<int> istartAll(P);
        if (ID != P - 1)
        {
            MPI_Send(&istart, 1, MPI_INT, P - 1, 0, MPI_COMM_WORLD);
        }
        if (ID == P - 1)
        {
            for (int i = 0; i < P - 1; i++)
            {
                MPI_Recv(&istartAll[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            istartAll[P-1] = istart;
        }
        MPI_Bcast(&istartAll[0], P, MPI_INT, P - 1, MPI_COMM_WORLD);
        int gfront;
        int gback;
        int gfrontLength = min(N, istart);
        int gbackLength = min(N, N*N-istart-Nc);
        for (int i = 1; i < P; i++)
        {
            gfront = max(istartAll[i] - N, 0);
            for (int j = 0; j < i; j++)
            {
                if (gfront < istartAll[j+1])
                {
                    recvFrom[i].push_back(j);
                    recvStart[i].push_back(max(gfront-istartAll[j], 0));
                    recvLength[i].push_back(min(istartAll[j+1]-gfront, istartAll[j+1]-istartAll[j]));
                }
            }
        }
        for (int i = 0; i < P-1; i++)
        {
            gback = min(istartAll[i+1]+N-1, N*N-1);
            for (int j = i + 1; j < P; j++)
            {
                if (gback > istartAll[j] - 1)
                {
                    recvFrom[i].push_back(j);
                    recvStart[i].push_back(0);
                    recvLength[i].push_back(min(gback-istartAll[j]+1, istartAll[j]-istartAll[j-1]));
                }
            }
        }
        for (int i = 0; i < P; i++)
        {
            if (i != ID)
            {
                for (int j = 0; j < recvFrom[i].size(); j++)
                {
                    if (recvFrom[i][j] == ID)
                    {
                        sendTo.push_back(i);
                        sendStart.push_back(recvStart[i][j]);
                        sendLength.push_back(recvLength[i][j]);
                    }
                }
            }
        }
        hc.resize(gfrontLength+Nc+gbackLength);
        int nSendhc = sendTo.size();
        int nRecvhc = recvFrom[ID].size();
        vector<MPI_Request> reqSendhc(nSendhc);
        vector<MPI_Request> reqRecvhc(nRecvhc);
        int frontidx = 0;
        int backidx = gfrontLength+Nc;
        int nreq = 0;
        for (int i = 0; i < sendTo.size(); i++)
        {
            MPI_Isend(&Uc[0][sendStart[i]+igap], sendLength[i], MPI_DOUBLE, sendTo[i], 0, 
                            MPI_COMM_WORLD, &reqSendhc[nreq]);
            nreq++;
        }
        nreq = 0;
        for (int i = 0; i < recvFrom[ID].size(); i++)
        {
            if (recvFrom[ID][i] < ID)
            {
                MPI_Irecv(&hc[frontidx], recvLength[ID][i], MPI_DOUBLE, recvFrom[ID][i], 0, 
                                MPI_COMM_WORLD, &reqRecvhc[nreq]);
                frontidx += recvLength[ID][i];
            }
            else
            {
                MPI_Irecv(&hc[backidx], recvLength[ID][i], MPI_DOUBLE, recvFrom[ID][i], 0, 
                                MPI_COMM_WORLD, &reqRecvhc[nreq]);
                backidx += recvLength[ID][i];
            }
            nreq++;
        }
        copy(&Uc[0][igap], &Uc[0][igap+Nc], &hc[gfrontLength]);
        MPI_Waitall(nSendhc, reqSendhc.data(), MPI_STATUS_IGNORE);
        MPI_Waitall(nRecvhc, reqRecvhc.data(), MPI_STATUS_IGNORE);
        //calculate the refinement indicator
        int itemp;
        int ig;
        for (int i = 0; i < 2; i++)
            gradh[i].resize(Nc);
        for (int i = 0; i < Nc; i++)
        {
            itemp = i + istart;    //global index
            ig = i + gfrontLength;    //local index
            if (itemp%N==0)    //left boundary
                gradh[0][i] = (hc[ig+1] - hc[ig])/cL;
            else if (itemp%N==N-1)    //right boundary
                gradh[0][i] = (hc[ig] - hc[ig-1])/cL;
            else
                gradh[0][i] = (hc[ig+1] - hc[ig-1])/(2*cL);
            if (itemp<N)    //lower boundary
                gradh[1][i] = (hc[ig+N] - hc[ig])/cL;
            else if (itemp>=N*(N-1))    //upper boundary
                gradh[1][i] = (hc[ig] - hc[ig-N])/cL;
            else
                gradh[1][i] = (hc[ig+N] - hc[ig-N])/(2*cL);
        }

        //update the indicator function
        double gradhmag;
        I.resize(Nc);
        for (int i = 0; i < Nc; i++)
        {
            gradhmag = sqrt(gradh[0][i]*gradh[0][i] + gradh[1][i]*gradh[1][i]);
            if (gradhmag > GHcr)
                I[i] = 1;
            else
                I[i] = 0;
        }

        if (ID == 0)
        {
            tiend = MPI_Wtime();
            ti += tiend - ti0;
        }

        //load balancing and repartition
        if (ID == 0)
        {
            tlb0 = MPI_Wtime();
        }

        vector<int> istartOld(P);
        vector<int> iendOld(P);
        vector<int> istartNew(P);
        vector<int> iendNew(P);
        vector<MPI_Request> reqSendbarOld(P-1);
        vector<MPI_Request> reqRecvbarOld(P-1);
        if (ID != P - 1)
        {
            MPI_Isend(&istart, 1, MPI_INT, P - 1, 0, MPI_COMM_WORLD, &reqSendbarOld[ID]);
        }
        if (ID == P - 1)
        {
            for (int i = 0; i < P - 1; i++)
            {
                MPI_Irecv(&istartOld[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &reqRecvbarOld[i]);
            }
            istartOld[P-1] = istart;
            MPI_Waitall(P-1, reqRecvbarOld.data(), MPI_STATUS_IGNORE);
        }
        MPI_Bcast(&istartOld[0], P, MPI_INT, P - 1, MPI_COMM_WORLD);
        int sumLoc = 0;
        for (int i = 0; i < Nc; i++)
        {
            if (I[i])
                sumLoc += 4;
            else
                sumLoc += 1;
        }
        int sumScan = 0;
        int sumAll = 0;
        MPI_Scan(&sumLoc, &sumScan, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        sumAll = sumScan;
        MPI_Bcast(&sumAll, 1, MPI_INT, P-1, MPI_COMM_WORLD);
        int sumFront = sumScan - sumLoc;
        int avg = sumAll / P;
        vector<MPI_Request> reqSendbarNew(P);
        vector<MPI_Request> reqRecvbarNew(P);
        int barNum = int(ceil(float(sumFront)/float(avg)));
        if (ID == P - 1)
        {   
            for (int i = 0; i < P; i++)
            {
                MPI_Irecv(&istartNew[i], 1, MPI_INT, MPI_ANY_SOURCE, i, MPI_COMM_WORLD, &reqRecvbarNew[i]);
            }
        }
        for (int i = 0; i < Nc; i++)  
        {
            if (I[i])
                sumFront += 4;
            else
                sumFront += 1;
            if (sumFront > avg * barNum)
            {  
                istartNew[barNum] = istart + i + 1; 
                MPI_Isend(&istartNew[barNum], 1, MPI_INT, P - 1, barNum, MPI_COMM_WORLD, &reqSendbarNew[barNum]);
                barNum++; 
            }
            if (barNum > P - 1)
                break;
        }
        if (ID == P - 1)
        {
            MPI_Waitall(P , reqRecvbarNew.data(), MPI_STATUS_IGNORE);
            istartNew[0] = 0;
        }
        MPI_Bcast(&istartNew[0], P, MPI_INT, P - 1, MPI_COMM_WORLD);
        for (int i = 0; i < P - 1; i++)
        {
            iendOld[i] = istartOld[i+1] - 1;
            iendNew[i] = istartNew[i+1] - 1;
        }
        iendOld[P-1] = N*N - 1;
        iendNew[P-1] = N*N - 1;
        Nc = iendNew[ID] - istartNew[ID] + 1;
        Nf = Nc*4;
        istart = istartNew[ID];
        for (int i = 0; i < P - 1; i++)
        {
            iendNew[i] = min(iendNew[i]+N, N*N - 1);
            istartNew[i+1] = max(istartNew[i+1]-N, 0);
        }
        igstart = istartNew[ID];
        for (int i = 0; i < P; i++)
        {
            recvFrom[i].clear();
            recvStart[i].clear();
            recvLength[i].clear();
        }
        sendTo.clear();
        sendStart.clear();
        sendLength.clear();
        int kstart;
        int kend;
        int kk;
        for (int i = 0; i < P; i++)
        {
            kstart = istartNew[i];
            kend = iendNew[i];
            kk = kstart;
            for (int j = 0; j < P; j++)
            {
                if (iendOld[j]>=kstart && iendOld[j]<kend)
                {
                    recvFrom[i].push_back(j);
                    recvStart[i].push_back(kk - istartOld[j]);
                    recvLength[i].push_back(iendOld[j] - kk + 1);
                    kk = iendOld[j] + 1;
                }
                if (iendOld[j]>=kend)
                {
                    recvFrom[i].push_back(j);
                    recvStart[i].push_back(kk - istartOld[j]);
                    recvLength[i].push_back(kend - kk + 1);
                    break;
                }
            }
        }
        for (int i = 0; i < P; i++)
        {
            for (int j = 0; j < recvFrom[i].size(); j++)
            {
                if (recvFrom[i][j] == ID)
                {
                    sendTo.push_back(i);
                    sendStart.push_back(recvStart[i][j]);
                    sendLength.push_back(recvLength[i][j]);
                }
            }
        }
        //data transfer between processors (indicator function I)
        int ulength = iendNew[ID] - istartNew[ID] + 1;
        vector<int> Itemp(I);
        I.clear();
        I.resize(ulength);
        int nSendI = sendTo.size();
        int nRecvI = recvFrom[ID].size();
        vector<MPI_Request> reqSendI(nSendI);
        vector<MPI_Request> reqRecvI(nRecvI);
        nreq = 0;
        for (int j = 0; j < sendTo.size(); j++)
        {
            MPI_Isend(&Itemp[sendStart[j]], sendLength[j], MPI_INT, sendTo[j], 0,
                MPI_COMM_WORLD, &reqSendI[nreq]);
            nreq++;
        }
        int ilength = 0;
        nreq = 0;
        for (int j = 0; j < recvFrom[ID].size(); j++)
        {
            MPI_Irecv(&I[ilength], recvLength[ID][j], MPI_INT, recvFrom[ID][j], 0,
                    MPI_COMM_WORLD, &reqRecvI[nreq]);
            nreq++;
            ilength += recvLength[ID][j];
        }
        MPI_Waitall(nSendI, reqSendI.data(), MPI_STATUS_IGNORE);
        MPI_Waitall(nRecvI, reqRecvI.data(), MPI_STATUS_IGNORE);
        //data transfer between processors (Uc and Uf combined)
        int nSend = sendTo.size()*3;
        int nRecv = recvFrom[ID].size()*3;
        vector<MPI_Request> reqSend(nSend);
        vector<MPI_Request> reqRecv(nRecv);
        vector<vector<vector<double> > > Usend(sendTo.size(), vector<vector<double> >(3, vector<double>(0)));
        vector<vector<vector<double> > > Urecv(recvFrom[ID].size(), vector<vector<double> >(3, vector<double>(0)));
        int idxc;
        nreq = 0;
        for (int i = 0; i < sendTo.size(); i++)
        {
            for (int j = 0; j < sendLength[i]; j++)
            {
                idxc = sendStart[i] + j;
                if (Itemp[idxc])
                {
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 4; l++)
                            Usend[i][k].push_back(Uf[k][(idxc+igap)*4+l]);
                    }
                }
                else
                {
                    for (int k = 0; k < 3; k++)
                        Usend[i][k].push_back(Uc[k][idxc+igap]);
                }
            }
            for (int j = 0; j < 3; j++)
            {
                MPI_Isend(&Usend[i][j][0], Usend[i][j].size(), MPI_DOUBLE, sendTo[i], j+1,
                    MPI_COMM_WORLD, &reqSend[nreq]);
                nreq++;
            }
        }
        nreq = 0;
        int iidx = 0;
        for (int i = 0; i < recvFrom[ID].size(); i++)
        {
            int urlength = 0;
            for (int j = 0; j < recvLength[ID][i]; j++)
            {
                if (I[iidx])
                    urlength += 4;
                else
                    urlength++;
                iidx++;
            }
            for (int j = 0; j < 3; j++)
            {
                Urecv[i][j].resize(urlength);
                MPI_Irecv(&Urecv[i][j][0], urlength, MPI_DOUBLE, recvFrom[ID][i], j+1, 
                            MPI_COMM_WORLD, &reqRecv[nreq]);
                nreq++;
            }
        }
        MPI_Waitall(nSend, reqSend.data(), MPI_STATUS_IGNORE);
        MPI_Waitall(nRecv, reqRecv.data(), MPI_STATUS_IGNORE);
        //reallocate the transferred combined data into Uc and Uf
        for (int i = 0; i < 3; i++)
        {
            Uc[i].clear();
            Uf[i].clear();
            Uc[i].resize(ulength);
            Uf[i].resize(ulength*4);
        }
        idxc = 0;
        for (int i = 0; i < recvFrom[ID].size(); i++)
        {
            iidx = 0;
            for (int j = 0; j < recvLength[ID][i]; j++)
            {
                if (I[idxc])
                {
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 4; l++)
                        {
                            Uf[k][idxc*4+l] = Urecv[i][k][iidx+l];
                        }
                    }
                    iidx += 4;
                }
                else
                {
                    for (int k = 0; k < 3; k++)
                    {
                        Uc[k][idxc] = Urecv[i][k][iidx];
                    }
                    iidx++;
                }
                idxc++;
            }
        }

        if (ID == 0)
        {
            tlbend = MPI_Wtime();
            tlb += tlbend - tlb0;
        }

        if (ID == 0)
        {
            tup0 = MPI_Wtime();
        }

        //initialize the residual vectors
        for (int i = 0; i < 3; i++)
        {
            Rc[i].assign(Nc, 0);
        }
        for (int i = 0; i < 3; i++)
        {
            Rf[i].assign(Nf, 0);
        }

        //loop over interior edges
        double nvec[2];
        double F[3];
        double stemp;
        double uL[3];
        double uR[3];
        int elem1;
        int elem2;
        int elemL;
        int elemR;
        int iin;
        bool inRangeL;
        bool inRangeR;
        for (int i = 0; i < ulength; i++)
        {
            elem1 = i + igstart;
            inRangeL = (elem1>=istart) && (elem1<istart+Nc);
            //F2F
            if (I[i])
            {
                for (int j = 0; j < 4; j++)
                {
                    elemL = CF[i+igstart][j];
                    for (int k = 0; k < F2F[elemL].size(); k++)
                    {
                        elemR = F2F[elemL][k];
                        elem2 = FC[elemR];
                        inRangeR = (elem2>=istart) && (elem2<istart+Nc);
                        if (inRangeL || inRangeR)
                        {
                            if (I[FC[elemR]-igstart])
                            {
                                iin = F2Fn[elemL][k];
                                for (int l = 0; l < 3; l++)
                                {
                                    uL[l] = Uf[l][elemL-igstart*4];
                                    uR[l] = Uf[l][elemR-igstart*4];
                                }
                                memcpy(nvec, normal[iin], sizeof(nvec));
                                flux(uL, uR, nvec, F, &stemp);
                                if (inRangeL)
                                {
                                    for (int l = 0; l < 3; l++)
                                        Rf[l][elemL-istart*4] += F[l]*cLf;
                                }
                                if (inRangeR)
                                {
                                    for (int l = 0; l < 3; l++)
                                        Rf[l][elemR-istart*4] -= F[l]*cLf;
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                elemL = i + igstart;
                //C2C
                for (int j = 0; j < C2C[elemL].size(); j++)
                {
                    elemR = C2C[elemL][j];
                    elem2 = elemR;
                    inRangeR = (elem2>=istart) && (elem2<istart+Nc);
                    if (inRangeL || inRangeR)
                    {
                        if (I[elemR-igstart] == 0)
                        {
                            iin = C2Cn[elemL][j];
                            for (int l = 0; l < 3; l++)
                            {
                                uL[l] = Uc[l][elemL-igstart];
                                uR[l] = Uc[l][elemR-igstart];
                            }
                            memcpy(nvec, normal[iin], sizeof(nvec));
                            flux(uL, uR, nvec, F, &stemp);
                            if (inRangeL)
                            {
                                for (int l = 0; l < 3; l++)
                                    Rc[l][elemL-istart] += F[l]*cL;
                            }
                            if (inRangeR)
                            {
                                for (int l = 0; l < 3; l++)
                                    Rc[l][elemR-istart] -= F[l]*cL;
                            }
                        }
                    }
                }
                //C2F
                for (int j = 0; j < C2F[elemL].size(); j++)
                {
                    elemR = C2F[elemL][j];
                    elem2 = FC[elemR];
                    inRangeR = (elem2>=istart) && (elem2<istart+Nc);
                    if (inRangeL || inRangeR)
                    {
                        if (I[FC[elemR]-igstart])
                        {
                            iin = C2Fn[elemL][j];
                            for (int l = 0; l < 3; l++)
                            {
                                uL[l] = Uc[l][elemL-igstart];
                                uR[l] = Uf[l][elemR-igstart*4];
                            }
                            memcpy(nvec, normal[iin], sizeof(nvec));
                            flux(uL, uR, nvec, F, &stemp);
                            if (inRangeL)
                            {
                                for (int l = 0; l < 3; l++)
                                    Rc[l][elemL-istart] += F[l]*cLf;
                            }
                            if (inRangeR)
                            {
                                for (int l = 0; l < 3; l++)
                                    Rf[l][elemR-istart*4] -= F[l]*cLf;
                            }
                        }
                    }
                }
            }
        }
        //loop over boundary edges
        int elemB;
        double hp;
        for (int i = 0; i < Nc; i++)
        {
            if (I[i+istart-igstart])
            {
                //F2B
                for (int j = 0; j < 4; j++)
                {
                    elemB = CF[i+istart][j];
                    for (int k = 0; k < F2Bn[elemB].size(); k++)
                    {
                        iin = F2Bn[elemB][k];
                        memcpy(nvec, normal[iin], sizeof(nvec));
                        hp = Uf[0][elemB-igstart*4];
                        F[0] = 0;
                        F[1] = 0.5*g*hp*hp*nvec[0];
                        F[2] = 0.5*g*hp*hp*nvec[1];
                        for (int l = 0; l < 3; l++)
                            Rf[l][elemB-istart*4] += F[l]*cLf;
                    }
                }
            }
            else
            {
                //C2B
                elemB = i + istart;
                for (int j = 0; j < C2Bn[elemB].size(); j++)
                {
                    iin = C2Bn[elemB][j];
                    memcpy(nvec, normal[iin], sizeof(nvec));
                    hp = Uc[0][elemB-igstart];
                    F[0] = 0;
                    F[1] = 0.5*g*hp*hp*nvec[0];
                    F[2] = 0.5*g*hp*hp*nvec[1];
                    for (int l = 0; l < 3; l++)
                        Rc[l][elemB-istart] += F[l]*cL;
                }
            }
        }

        //update cell states
        int itempc;
        int itempf;
        for (int i = 0; i < Nc; i++)
        {
            itempc = i + istart - igstart; 
            if (I[itempc])
            {
                for (int j = 0; j < 3; j++)
                {
                    Uc[j][itempc] = 0;
                }
                for (int j = 0; j < 4; j++)
                {
                    itempf = 4*itempc + j;
                    for (int k = 0; k < 3; k++)
                    {
                        Uf[k][itempf] -= dt/areaf*Rf[k][4*i+j];
                        Uc[k][itempc] += Uf[k][itempf];
                    }
                }
                for (int j = 0; j < 3; j++)
                {
                    Uc[j][itempc] /= 4;
                }
            }
            else
            {
                for (int j = 0; j < 3; j++)
                {
                    Uc[j][itempc] -= dt/area*Rc[j][i];
                }
                for (int j = 0; j < 4; j++)
                {
                    itempf = 4*itempc + j;
                    for (int k = 0; k < 3; k++)
                    {
                        Uf[k][itempf] = Uc[k][itempc];
                    }
                }
            }
        }

        if (ID == 0)
        {
            tupend = MPI_Wtime();
            tup += tupend - tup0;
        }
        
    }

    if (ID == 0)
    {
        tend = MPI_Wtime();
        t = tend - t0;
        cout << "Elapsed time: " << t << "s" << endl;
        cout << "Indicator calculating time: " << ti << "s" << endl;
        cout << "Load balancing time: " << tlb << "s" << endl;
        cout << "Updating time: " << tup << "s" << endl;
        cout << "nt = " << nt << endl;
    }

    /*
    int itempc;
    int itempf;
    ofstream outfile;
    for (int iid = 0; iid < P; iid++)
    {
        if (ID == iid)
        {
            outfile.open("output.dat", ios::app);
            //if (ID == 0)
            //{
            //    outfile << "time step: " << nt - 1 << endl;
            //}
            for (int i = 0; i < Nc; i++)
            {
                itempc = i + istart - igstart; 
                if (I[itempc])
                {
                    for (int j = 0; j < 4; j++)
                    {
                        itempf = itempc*4 + j;
                        outfile << I[itempc] << " " << (istart+i)*4+j << " " << 
                                Uf[0][itempf] << " " << Uf[1][itempf] << " " << Uf[2][itempf] <<endl;
                    }
                }
                else
                {
                    outfile << I[itempc] << " " << istart+i << " " << 
                            Uc[0][itempc] << " " << Uc[1][itempc] << " " << Uc[2][itempc] <<endl;
                }
            }
            outfile.close();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    */

    // Finalize MPI
    MPI_Finalize();

    return 0;
}