#include "CouplingGraph.h"

//析构
CouplingGraph::~CouplingGraph(){
    for(int i=0; i<positions; i++)
        delete [] dist[i];
    delete [] dist;
}

//初始化dist和path
void CouplingGraph::initDist(){
    dist = new int*[positions];
    path = new int*[positions];
    for(int i=0; i<positions; i++)
    {
        dist[i] = new int[positions];
        path[i] = new int[positions];
    }

    for(int i=0; i<positions; i++)
        for(int j=0; j<positions; j++)
        {
            if(i==j)
                dist[i][j]=0;
            else
                dist[i][j]=MAX_DIST;
            path[i][j] = j;
        }


    //Initialize distance matrix to a adjacency matrix
    //把邻接表转为邻接矩阵
    for(int i=0;i<positions;i++){
        for(unsigned int j=0;j<adjList[i].size();j++){
            int k=adjList[i][j];
            dist[i][k]=1;
        }
    }

    computeDist();
}


//建图
void CouplingGraph::buildQX20(){
    fstream fin;
    fin.open("./device/QX20",ios::in);
    int n,m;
    fin>>n>>m;
    positions=n;
    adjList.resize(positions);
    deg.resize(positions);
    while(m--){
        int u,v;
        fin>>u>>v;
        adjList[u].push_back(v);
        adjList[v].push_back(u);
        deg[u]+=1;
        deg[v]+=1;
    }
    initDist();
    fin.close();
    cout<<"[INFO] build tokyo Over"<<endl;
}

void CouplingGraph::buildQX5()
{
    fstream fin;
    fin.open("./device/QX5",ios::in);
    int n,m;
    fin>>n>>m;
    positions=n;
    adjList.resize(positions);
    deg.resize(positions);
    while(m--){
        int u,v;
        fin>>u>>v;
        adjList[u].push_back(v);
        // adjList[v].push_back(u);
        deg[u]+=1;
        deg[v]+=1;
    }
    initDist();
    fin.close();
    cout<<"[INFO] build QX5 Over"<<endl;
}

void CouplingGraph::build_IBMQ_GUADALUPE()
{
    fstream fin;
    fin.open("./device/IBMQ_GUADALUPE",ios::in);
    int n,m;
    fin>>n>>m;
    positions=n;
    adjList.resize(positions);
    deg.resize(positions);
    while(m--){
        int u,v;
        fin>>u>>v;
        adjList[u].push_back(v);
        adjList[v].push_back(u);
        deg[u]+=1;
        deg[v]+=1;
    }
    initDist();
    fin.close();
    
    cout<<"[INFO] build IBMQ GUADALUPE Over"<<endl;
}

void CouplingGraph::build_IBM_Johannesburg(){
    fstream fin;
    fin.open("./device/IBM_Johannesburg",ios::in);
    int n,m;
    fin>>n>>m;
    positions=n;
    adjList.resize(positions);
    deg.resize(positions);
    while(m--){
        int u,v;
        fin>>u>>v;
        adjList[u].push_back(v);
        adjList[v].push_back(u);
        deg[u]+=1;
        deg[v]+=1;
    }
    initDist();
    fin.close();
    
    cout<<"[INFO] build IBM_Johannesburg Over"<<endl;
}
void CouplingGraph::build_Google_Sycamore(){
    fstream fin;
    fin.open("./device/Google_Sycamore",ios::in);
    int n,m;
    fin>>n>>m;
    positions=n;
    adjList.resize(positions);
    deg.resize(positions);
    while(m--){
        int u,v;
        fin>>u>>v;
        adjList[u].push_back(v);
        adjList[v].push_back(u);
        deg[u]+=1;
        deg[v]+=1;
    }
    initDist();
    fin.close();
    
    cout<<"[INFO] build Google Sycamore Over"<<endl;
}
void CouplingGraph::build_IBM_Rochester()
{
    fstream fin;
    fin.open("./device/IBM_Rochester",ios::in);
    int n,m;
    fin>>n>>m;
    positions=n;
    adjList.resize(positions);
    deg.resize(positions);
    while(m--){
        int u,v;
        fin>>u>>v;
        adjList[u].push_back(v);
        adjList[v].push_back(u);
        deg[u]+=1;
        deg[v]+=1;
    }
    initDist();
    fin.close();
    cout<<"[INFO] build IBM Rochester Over"<<endl;
}
void CouplingGraph::build_IBM_Penguin_v1()
{
    fstream fin;
    fin.open("./device/Penguin_v1",ios::in);
    int n,m;
    fin>>n>>m;
    positions=n;
    adjList.resize(positions);
    deg.resize(positions);
    while(m--){
        int u,v;
        fin>>u>>v;
        adjList[u].push_back(v);
        adjList[v].push_back(u);
        deg[u]+=1;
        deg[v]+=1;
    }
    initDist();
    fin.close();
}
void CouplingGraph::build_IBM_Penguin_v2()
{
    fstream fin;
    fin.open("./device/Penguin_v2",ios::in);
    int n,m;
    fin>>n>>m;
    positions=n;
    adjList.resize(positions);
    deg.resize(positions);
    while(m--){
        int u,v;
        fin>>u>>v;
        adjList[u].push_back(v);
        adjList[v].push_back(u);
        deg[u]+=1;
        deg[v]+=1;
    }
    initDist();
    fin.close();
}
void CouplingGraph::build_nXn_grid(int n=17)
{
    positions=n*n;
    adjList.resize(positions);
    // deg.resize(positions);

    vector<pair<int,int> >dxy={{1,0},{0,1},{1,1},{1,-1}};//完全网格

    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            int u=i*n+j;
            for(auto d:dxy){
                int x=i+d.first;
                int y=j+d.second;
                if(x>=0&&x<n&&y>=0&&y<n){
                    int v=x*n+y;
                    adjList[u].push_back(v);
                    adjList[v].push_back(u);
                    // deg[u]+=1;
                    // deg[v]+=1;
                    // cout<<u<<" "<<v<<endl;
                }
            }
        }
    }
    initDist();
    cout<<"[INFO] build "<<n<<"X"<<n<<"_grid Over"<<endl;
}
void CouplingGraph::buildQ16()
{
    positions=16;
    adjList.resize(positions);

    adjList[0].push_back(1);
    adjList[0].push_back(4);
    adjList[0].push_back(5);

    adjList[1].push_back(0);
    adjList[1].push_back(2);
    adjList[1].push_back(4);
    adjList[1].push_back(5);

    adjList[2].push_back(1);
    adjList[2].push_back(3);
    adjList[2].push_back(6);
    adjList[2].push_back(7);

    adjList[3].push_back(2);
    adjList[3].push_back(6);
    adjList[3].push_back(7);

    adjList[4].push_back(0);
    adjList[4].push_back(1);
    adjList[4].push_back(5);
    adjList[4].push_back(8);

    adjList[5].push_back(0);
    adjList[5].push_back(1);
    adjList[5].push_back(4);
    adjList[5].push_back(6);
    adjList[5].push_back(9);
    adjList[5].push_back(10);

    adjList[6].push_back(2);
    adjList[6].push_back(3);
    adjList[6].push_back(5);
    adjList[6].push_back(7);
    adjList[6].push_back(9);
    adjList[6].push_back(10);

    adjList[7].push_back(2);
    adjList[7].push_back(3);
    adjList[7].push_back(6);
    adjList[7].push_back(11);

    adjList[8].push_back(4);
    adjList[8].push_back(9);
    adjList[8].push_back(12);
    adjList[8].push_back(13);

    adjList[9].push_back(5);
    adjList[9].push_back(6);
    adjList[9].push_back(8);
    adjList[9].push_back(10);
    adjList[9].push_back(12);
    adjList[9].push_back(13);

    adjList[10].push_back(5);
    adjList[10].push_back(6);
    adjList[10].push_back(9);
    adjList[10].push_back(11);
    adjList[10].push_back(14);
    adjList[10].push_back(15);

    adjList[11].push_back(7);
    adjList[11].push_back(10);
    adjList[11].push_back(14);
    adjList[11].push_back(15);

    adjList[12].push_back(8);
    adjList[12].push_back(9);
    adjList[12].push_back(13);

    adjList[13].push_back(8);
    adjList[13].push_back(9);
    adjList[13].push_back(12);
    adjList[13].push_back(14);

    adjList[14].push_back(10);
    adjList[14].push_back(11);
    adjList[14].push_back(13);
    adjList[14].push_back(15);

    adjList[15].push_back(10);
    adjList[15].push_back(11);
    adjList[15].push_back(14);
    //Initialize distance matrix to a adjacency matrix
    initDist();
}

void CouplingGraph::buildQ10(){
    positions=10;
    adjList.resize(positions);

    adjList[0].push_back(1);
    adjList[0].push_back(3);

    adjList[1].push_back(0);
    adjList[1].push_back(2);
    adjList[1].push_back(4);
    adjList[1].push_back(5);

    adjList[2].push_back(1);
    adjList[2].push_back(4);
    adjList[2].push_back(5);

    adjList[3].push_back(0);
    adjList[3].push_back(4);
    adjList[3].push_back(6);
    adjList[3].push_back(7);

    adjList[4].push_back(1);
    adjList[4].push_back(2);
    adjList[4].push_back(3);
    adjList[4].push_back(5);
    adjList[4].push_back(6);
    adjList[4].push_back(7);

    adjList[5].push_back(1);
    adjList[5].push_back(2);
    adjList[5].push_back(4);
    adjList[5].push_back(8);
    adjList[5].push_back(9);

    adjList[6].push_back(3);
    adjList[6].push_back(4);
    adjList[6].push_back(7);

    adjList[7].push_back(3);
    adjList[7].push_back(4);
    adjList[7].push_back(6);
    adjList[7].push_back(8);

    adjList[8].push_back(5);
    adjList[8].push_back(7);
    adjList[8].push_back(9);

    adjList[9].push_back(5);
    adjList[9].push_back(8);

    //Initialize distance matrix to a adjacency matrix
    initDist();
}

void CouplingGraph::buildQ9(){
    positions=9;
    adjList.resize(positions);

    adjList[0].push_back(1);
    adjList[0].push_back(3);

    adjList[1].push_back(0);
    adjList[1].push_back(2);
    adjList[1].push_back(4);
    adjList[1].push_back(5);

    adjList[2].push_back(1);
    adjList[2].push_back(4);
    adjList[2].push_back(5);

    adjList[3].push_back(0);
    adjList[3].push_back(4);
    adjList[3].push_back(6);
    adjList[3].push_back(7);

    adjList[4].push_back(1);
    adjList[4].push_back(2);
    adjList[4].push_back(3);
    adjList[4].push_back(5);
    adjList[4].push_back(6);
    adjList[4].push_back(7);

    adjList[5].push_back(1);
    adjList[5].push_back(2);
    adjList[5].push_back(4);
    adjList[5].push_back(8);

    adjList[6].push_back(3);
    adjList[6].push_back(4);
    adjList[6].push_back(7);

    adjList[7].push_back(3);
    adjList[7].push_back(4);
    adjList[7].push_back(6);
    adjList[7].push_back(8);

    adjList[8].push_back(5);
    adjList[8].push_back(7);

    //Initialize distance matrix to a adjacency matrix
    initDist();
}


void CouplingGraph::buildSquare(){
    positions=6;
    adjList.resize(positions);
    adjList[0].push_back(1);
    adjList[0].push_back(2);
    adjList[1].push_back(0);
    adjList[1].push_back(3);
    adjList[2].push_back(0);
    adjList[2].push_back(3);
    adjList[2].push_back(4);
    adjList[3].push_back(1);
    adjList[3].push_back(2);
    adjList[3].push_back(5);
    adjList[4].push_back(2);
    adjList[4].push_back(5);
    adjList[5].push_back(3);
    adjList[5].push_back(4);

    //Initialize distance matrix to a adjacency matrix
    initDist();
}

//Floyd-Warshall Algorithm
 void CouplingGraph::computeDist(){
     for(int k = 0 ; k < positions ; k ++)
    {
        for(int i = 0 ; i < positions ; i ++)
        {
            for(int j = 0 ; j < positions ; j ++)
            {
                if(dist[i][j] > dist[i][k] + dist[k][j])
                {
                    dist[i][j] = dist[i][k] + dist[k][j];
                    path[i][j] = path[i][k];
                }
            }
        }
    }
 }

//找路径
 vector<int> CouplingGraph::getPath(int i, int j){
    vector<int> pathBetween;
    pathBetween.push_back(i);
    int k = path[i][j];
    while(k!=j){
        pathBetween.push_back(k);
        k = path[k][j];
    }
    pathBetween.push_back(j);
    return pathBetween;
 }

//打印
 void CouplingGraph::showInfo(){
    for(unsigned int i=0; i<adjList.size(); i++)
    {
        cout<<i<<":";
        for(unsigned int j=0;j<adjList[i].size();j++)
            cout<<adjList[i][j]<<",";
        cout<<std::endl;
    }
 }
