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
 * by Michael Borinsky https://arxiv.org/abs/2508.14263
 */

double omega( int L, int n, double D, int k )
{
    assert( ((L-1)*k + n) % (k-2) == 0 );

    int E = ((L-1)*k + n)/(k-2);

    return E - L*D/2;
}

void fill_Z_B_tables( 
    long long maxL, 
    long long maxn, 
    double D,
    int k,
    array< array< double, MAXN+2*MAXL+1 >, MAXL+1 > &Z,
    array< array< double, MAXN+2*MAXL+1 >, MAXL+1 > &B
 )
{
    Z[0][k] = 1.0/tgamma(k+1);

    for( int L = 0; L <= maxL; L++ )
    {
        for( int n = 2; n <= 2*(maxL-L)+maxn; n++ )
        {
            B[L][n] = (n-1)*n*Z[L][n];

            for( int LL = 0; LL <= L; LL++ )
                for( int nn = 0; nn <= n-2; nn++ )
                    B[L][n] += (nn+1)*(nn+2)*Z[LL][nn+2]*B[L-LL][n-nn];

            if( L < maxL )
            {
                double w = omega(L+1, n-2, D, k);
                if( w > 1e-6 )
                    Z[L+1][n-2] = B[L][n]/(2*w);
            }
        }
    }

/*
    for( int L = 0; L <= maxL; L++ )
    {
        cout << L << ": ";
        for( int n = 2; n <= 2*(maxL-L)+maxn; n++ )
        {
            cout << "( " << n << " : " << B[L][n] << " | " << Z[L][n] << "), ";
        }
        cout << endl;
    }
*/
}


template< class Generator >
bool sample_LL_nn( 
    int L, 
    int n, 
    int& LL,
    int& nn,
    const array< array< double, MAXN+2*MAXL+1 >, MAXL+1 > &Z,
    const array< array< double, MAXN+2*MAXL+1 >, MAXL+1 > &B,
    Generator& rng 
)
{
    double R = true_random::uniform(rng) * B[L][n];

    R -= (n-1)*n*Z[L][n];

    if( R <= 0 )
        return true;

    for( LL = 0; LL <= L; LL++ )
    {
        for( nn = 0; nn <= n-2; nn++ )
        {
            R -= (nn+1)*(nn+2)*Z[LL][nn+2]*B[L-LL][n-nn];

            if( R <= 0 )
                return false;
        }
    }

    return false;
}

template< class Generator >
void sample_beaded_graph( 
    int L, 
    int n, 
    int k,
    graph& gg, 
    array<int, MAXE>& parent,
    array<pair<int, int>, MAXE>& type,
    const array< array< double, MAXN+2*MAXL+1 >, MAXL+1 > &Z,
    const array< array< double, MAXN+2*MAXL+1 >, MAXL+1 > &B,
    Generator& rng 
)
{
    int LL, nn;
    if( sample_LL_nn( L, n, LL, nn, Z, B, rng ) ) // returns true for case *
    {
        sample_1PI_graph( L, n, k, gg, parent, type, Z, B, rng );
        return;
    }

    sample_1PI_graph( LL, nn+2, k, gg, parent, type, Z, B, rng );

    graph gra;
    array<int, MAXE> gra_parent;
    array<pair<int, int>, MAXE> gra_type;
    sample_beaded_graph( L-LL, n-nn, k, gra, gra_parent, gra_type, Z, B, rng );

    int N = 2*gg.E + gg.n;
    int NE = gg.E;
    for( int v = 0; v < gra.V; v++ )
    {
        array<int, MAXK> vert;
        for( int i = 0; i < k; i++ )
            vert[i] = gra.vertices[v][i] + N;

        gg.vertices[gg.V++] = vert;
    }

    for( int e = 0; e < gra.E; e++ )
    {
        int h1, h2;
        tie(h1, h2) = gra.edges[e];

        int ee = gg.E++;
        gg.edges[ee] = make_pair(h1+N,h2+N);

        if( gra_parent[e] != -1 )
            parent[ee] = gra_parent[e] + NE;
        else
            parent[ee] = -1;

        type[ee] = gra_type[e];
    }

    int i1, i2;
    i1 = gg.legs[gg.n-1];
    i2 = gra.legs[gra.n-2];

    int ee = gg.E++;
    gg.edges[ee] = make_pair(i1, i2+N);
    parent[ee] = -1;

    i1 = gg.legs[gg.n-2];
    i2 = gra.legs[gra.n-1];
    
    gg.n -= 2;
    gra.n -= 2;
    
    for( int l = 0; l < gra.n; l++ )
    {
        int h = gra.legs[l];
        gg.legs[gg.n++] = h + N;
    }

    shuffle( begin(gg.legs), begin(gg.legs) + gg.n, rng );

    gg.legs[gg.n++] = i1;
    gg.legs[gg.n++] = i2 + N;
}

template< class Generator >
void sample_1PI_graph( 
    int L, 
    int n, 
    int k,
    graph& gg, 
    array<int, MAXE>& parent,
    array<pair<int, int>, MAXE>& type,
    const array< array< double, MAXN+2*MAXL+1 >, MAXL+1 > &Z,
    const array< array< double, MAXN+2*MAXL+1 >, MAXL+1 > &B,
    Generator& rng 
)
{
    if( k == n && L == 0)
    {
        gg.E = 0;
        gg.V = 1;
        gg.n = k;
        gg.k = k;
 
        array<int, MAXK> vertex;
        for( int i = 0; i < k; i++ )
        {
            vertex[i] = i;
            gg.legs[i] = i;
        }

        gg.vertices = { vertex };
        type.fill( make_pair(-1, -1) );
        
        return;
    }

    sample_beaded_graph( L-1, n+2, k, gg, parent, type, Z, B, rng );

    int ee = gg.E++;
    for( int e = 0; e < ee; e++ )
    {
        if( parent[e] == -1 )
            parent[e] = ee;
    }

    int h1 = gg.legs[gg.n-2];
    int h2 = gg.legs[gg.n-1];
    gg.n -= 2;

    parent[ee] = -1;
    type[ee] = make_pair( L, n );
    gg.edges[ee] = make_pair( h1, h2 );
}

template< class Generator >
double sample_edge_lengths_tree( 
    array< double, MAXE >& X,
    int e, 
    double shift,
    double D,
    int k,
    const array<array<int, MAXE>, MAXE>& adj,
    const array<int, MAXE>& deg,
    const array<pair<int, int>, MAXE>& type,
    Generator& rng
)
{
    if( deg[e] == 0 )
    {
        X[e] = shift + log(1. - true_random::uniform(rng));
        return 0.0;
    }

    int L, n;
    tie(L, n) = type[e];
    assert( L != -1 && n != -1 );

    double w = omega( L, n, D, k );
    double val = shift + log(1. - true_random::uniform(rng)) / w;
    X[e] = val;
    double log_psi_tr = val;

    for( int i = 0; i < deg[e]; i++ )
        log_psi_tr += sample_edge_lengths_tree( X, adj[e][i], val, D, k, adj, deg, type, rng );

    return log_psi_tr;
}


template< class Generator >
double sample_edge_lengths(
    array< double, MAXE >& X,
    double D,
    int k,
    const graph& gg,
    const array<int, MAXE>& parent,
    const array<pair<int, int>, MAXE>& type,
    Generator& rng
)
{
    array<array<int, MAXE>, MAXE> adj;
    array<int, MAXE> deg = {};
    int root = -1;
    
    for( int e = 0; e < gg.E; e++ )
    {
        int u = parent[e];
        if( u == -1 )
        {   
            assert( root == -1 );
            root = e;
            continue;
        }

        adj[u][deg[u]++] = e;
    }
    assert( root != -1 );

    X[root] = 0;

    double log_psi_tr = 0;
    for( int i = 0; i < deg[root]; i++ )
        log_psi_tr += sample_edge_lengths_tree( X, adj[root][i], 0, D, k, adj, deg, type, rng );

    return log_psi_tr;
}

