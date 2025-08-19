/*
 * Amplitrop: A program for global tropical sampling
 * 
 * Author: Michael Borinsky
 * Affiliation: Perimeter Institute for Theoretical Physics & University of Waterloo
 * Year: 2025
 *
 * License: MIT License (see LICENSE file for details)
 * If you use this code in your research, please cite the paper
 * 'Tropicalized quantum field theory and global tropical sampling' 
 * by Michael Borinsky
 */

class graph
{
public:
    array< array<int, MAXK>, MAXV > vertices = {};
    array< pair<int, int>, MAXE > edges = {};
    array< int, MAXH > legs = {};

    int V=0;
    int E=0;
    int n=0;
    int k=0;
};

ostream& operator<<(ostream& os, const graph& gg) 
{
    array< int, MAXH > vertMap;

    for( int v = 0; v < gg.V; v++ )
    {
        for( int i = 0; i < gg.k; i++ )
        {
            int h = gg.vertices[v][i];
            vertMap[h] = v;
        }
    }

    for( int e = 0; e < gg.E; e++ )
    {
        int h1, h2;
        tie(h1, h2) = gg.edges[e];
        int v1 = vertMap[h1];
        int v2 = vertMap[h2];

        os << "( " << v1 << ", " << v2 << " ), ";
    }
    os << "[ ";
    for( int l = 0; l < gg.n; l++ )
    {
        int h = gg.legs[l];
        os << vertMap[h] << " ";
    }
    os << "]";
    os << endl;
    return os;
}


void tree_dfs( 
    int u, 
    int p, 
    array< pair< int, int >, MAXV >& parent, 
    array< int, MAXV >& depth, 
    const array< array< pair< int, int >, MAXK>, MAXV >& adj,
    int k
)
{
    for( int i = 0; i < k; i++ )
    {
        int v, e;
        tie(v, e) = adj[u][i];

        if( v == -1 || v == p )
            continue;

        parent[v] = make_pair( u, e );
        depth[v] = depth[u] + 1;

        tree_dfs( v, u, parent, depth, adj, k );
    }
}

void get_tree_path_vec( 
    Eigen::Ref<Eigen::VectorXd> edge_vec, 
    int u, 
    int v,
    const array< pair< int, int >, MAXV >& parent, 
    const array< int, MAXV >& depth
)
{
    edge_vec.setZero();

    while( u != v )
    {
        if( depth[u] > depth[v] )
        {
            int e;
            tie(u, e) = parent[u];

            edge_vec[e] = 1;
        }
        else
        {
            int e;
            tie(v, e) = parent[v];

            edge_vec[e] = -1;
        }
    }
}

void fundamental_basis( 
    Eigen::MatrixXd& mat, 
    const graph& gg, 
    const array< pair< int, int >, MAXE >& type 
)
{
    array< int, MAXH > hEMap;
    array< int, MAXH > partner;
    partner.fill(-1);

    for( int e = 0; e < gg.E; e++ )
    {
        int h1, h2;
        tie(h1, h2) = gg.edges[e];

        hEMap[h1] = e;
        hEMap[h2] = e;

        if( type[e] != make_pair(-1, -1) )
            continue;

        partner[h1] = h2;
        partner[h2] = h1;
    }
    
    array< int, MAXH > hVMap;
    for( int v = 0; v < gg.V; v++ )
    {
        for( int i = 0; i < gg.k; i++ )
        {
            int h = gg.vertices[v][i];

            hVMap[h] = v;
        }
    }

    array< array< pair< int, int >, MAXK >, MAXV > adj;
    for( int v = 0; v < gg.V; v++ )
    {
        for( int i = 0; i < gg.k; i++ )
        {
            int h = gg.vertices[v][i];

            if( partner[h] == -1 )
                adj[v][i] = make_pair(-1, -1);
            else
                adj[v][i] = make_pair(hVMap[partner[h]], hEMap[h]);
        }
    }

    array< pair< int, int >, MAXV > parent;
    array< int, MAXV > depth = {};
    tree_dfs( 0, -1, parent, depth, adj, gg.k );

    int nL = 0;
    for( int e = 0; e < gg.E; e++ )
    {
        if( type[e] == make_pair(-1, -1) )
            continue;

        int h1, h2;
        tie(h1, h2) = gg.edges[e];

        int u = hVMap[h1];
        int v = hVMap[h2];

        get_tree_path_vec( mat.col(nL), u, v, parent, depth );
        mat.col(nL)[e] = 1;
        nL++;
    }
}

bool next_combination( array<int, MAXE>& comb, int nC, int nE )
{
    int i = nC - 1;
    while (i >= 0 && comb[i] == nE - nC + i) {
        --i;
    }
    if (i < 0) return false;

    ++comb[i];
    for (int j = i + 1; j < nC; ++j) {
        comb[j] = comb[j - 1] + 1;
    }
    return true;
}


bool dfs_has_bridge(
    int u, 
    int parent,
    const array< array< int, MAXK>, MAXV >& adj, 
    int k, 
    array<bool, MAXV>& visited, 
    array<int, MAXV>& disc, 
    array<int, MAXV>& low, 
    int& timer
    )
{
    visited[u] = true;
    disc[u] = low[u] = ++timer;

    for( int i = 0; i < k; i++ )
    {
        int v = adj[u][i];

        if( !visited[v] )
        {
            if( dfs_has_bridge( v, u, adj, k, visited, disc, low, timer ) )
                return true;

            low[u] = min( low[u], low[v] );

            if( low[v] > disc[u] )
                return true;
        }
        else if( v != parent )
            low[u] = min(low[u], disc[v]);
    }

    return false;
}

bool is_k4_D4_primitive( const graph& gg )
{
    assert( gg.k == 4 );

    int nC = gg.E - 3;
    array<int, MAXE> comb;
    for( int e = 0; e < nC; e++ )
        comb[e] = e;
    
    array< int, MAXH > hVMap;

    for( int v = 0; v < gg.V; v++ )
    {
        for( int i = 0; i < gg.k; i++ )
        {
            int h = gg.vertices[v][i];
            hVMap[h] = v;
        }
    }

    array< array<int, MAXV>, MAXV > adjMat = {};
    
    for( int e = 0; e < gg.E; e++ )
    {
        int h1, h2;
        tie(h1, h2) = gg.edges[e];
        int v1 = hVMap[h1];
        int v2 = hVMap[h2];

        if( v1 == v2 )
            return false;
        
        if( v2 < v1 ) swap( v1, v2 );

        adjMat[v1][v2]++;

        if( adjMat[v1][v2] >= 2 )
            return false;
    }

    do
    {
        array< int, MAXH > partner;
        partner.fill(-1);

        for( int e = 0; e < nC; e++ )
        {
            int h1, h2;
            tie(h1, h2) = gg.edges[comb[e]];
            partner[h1] = h2;
            partner[h2] = h1;
        }

        bool too_many_legs = false;

        array< array<int, MAXK>, MAXV > adj;
        for( int v = 0; v < gg.V; v++ )
        {
            int num_legs = 0;
            for( int i = 0; i < gg.k; i++ )
            {
                int h = gg.vertices[v][i];

                if( partner[h] == -1 )
                {
                    adj[v][i] = v;
                    num_legs++;
                }
                else
                    adj[v][i] = hVMap[partner[h]];
            }
            
            if( num_legs > 2 )
            {
                too_many_legs = true;
                break;
            }
        }

        if( too_many_legs )
            continue;

        array<bool, MAXV> visited = {};
        array<int, MAXV> disc;
        array<int, MAXV> low;
        int timer = 0;

        if( dfs_has_bridge( 0, -1, adj, gg.k, visited, disc, low, timer ) )
            return false;
    }
    while( next_combination( comb, nC, gg.E ) );

    return true;
}

