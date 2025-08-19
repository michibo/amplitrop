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


#include <cxxopts.hpp>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include <xsum/xsum/xsum.hpp>

#include <chrono>

#include <Eigen/Dense>
#include <array>
 
#include <cmath>
#include <iostream>

#include "random.hpp"

using namespace std;


constexpr int MAXK = 4;   // maximal k in phi^k
constexpr int MINK = 3;   // minimal k in phi^k
constexpr int MAXL = 100; // maximal loop order
constexpr int MAXN = 4;   // maximal multiplicity (at maximal loop order)


constexpr int MAXE = (MAXK*(MAXL-1)+MAXN) / (MINK-2);
constexpr int MAXV = (2*(MAXL-1)+MAXN) / (MINK-2);
constexpr int MAXH = 2*MAXE+MAXN;

#include "graph.hpp"
#include "sampler.hpp"


int main(int argc, char** argv) 
{
    int L, n, k, D;
    long long N, seed;
    bool check_primitive;

    try
    {
        cxxopts::Options options("amplitrop", "Proof of concept implementation for global tropical sampling computations of perturbation theoretic expansions in phi^k theory");

        options.add_options()
          ("N", "Number of Samples", cxxopts::value<long long>())
          ("n", "Number of legs", cxxopts::value<int>())
          ("k", "k in phi^k-theory", cxxopts::value<int>())
          ("D", "Dimension D", cxxopts::value<int>())
          ("L", "Loop order", cxxopts::value<int>())
          ("S", "Seed", cxxopts::value<long long>()->default_value("31415926"))
          ("prim", "Primitive", cxxopts::value<bool>()->default_value("false"))
        ;

        auto result = options.parse(argc, argv);

        if( result.count("help") )
        {
            std::cout << options.help() << std::endl;
            return 0;
        }

        L = result["L"].as<int>();
        n = result["n"].as<int>();
        k = result["k"].as<int>();
        D = result["D"].as<int>();

        N = result["N"].as<long long>();
        seed = result["S"].as<long long>();

        check_primitive = result["prim"].as<bool>();

        if( L > MAXL || n - 2*(MAXL-L) > MAXN || k < MINK || k > MAXK )
        {
            cout << "This program was compiled with max loop order " << MAXL << ", max leg number " << MAXN << " at that loop order and k <= " << MAXK << " and k >= " << MINK << ". Increase the respective values in the code and recompile." << endl;
            return 1;
        }
    }
    catch (const cxxopts::exceptions::parsing& e) 
    {
        std::cerr << "Error parsing options: " << e.what() << std::endl;
        std::cerr << "Run with --help for usage information." << std::endl;
        return 1;
    }
    catch (const cxxopts::exceptions::specification& e)
    {
        std::cerr << "Error in option specification: " << e.what() << std::endl;
        return 1;
    }
    catch (const std::exception& e)
    {
        std::cerr << "An unexpected error occurred: " << e.what() << std::endl;
        return 1;
    }

    array< array< double, MAXN+2*MAXL+1 >, MAXL+1 > Z = {};
    array< array< double, MAXN+2*MAXL+1 >, MAXL+1 > B = {};

    fill_Z_B_tables( L, n, D, k, Z, B );

    true_random::xoshiro256 rand( seed );

    int max_threads = omp_get_max_threads();
    vector<true_random::xoshiro256> rngs;

    for( int i = 0; i < max_threads; i++ )
    {
        rngs.push_back(rand);
        rand.jump();
    }

    int M = 5;
    vector< vector< xsum::xsum_large > > moment_sum;
    for( int i = 0; i < M; i++ )
        moment_sum.push_back( vector< xsum::xsum_large >( max_threads ) );

    int nV = (2*(L-1)+n) / (k-2);
    int nE = (k*(L-1)+n) / (k-2);

    Eigen::MatrixXd laplace(L, L);
    Eigen::MatrixXd matCycle(nE, L);
        
    Eigen::VectorXd Xvec(nE);
    Eigen::LDLT<Eigen::MatrixXd> ldlt;

    size_t Nprim = 0;

    auto start = chrono::high_resolution_clock::now();

    #pragma omp parallel for reduction(+:Nprim) schedule(dynamic) firstprivate(laplace,matCycle,Xvec,ldlt)
    for( long long i = 0; i < N; i++ )
    {
        int tid = omp_get_thread_num();
        graph gg;
        array<int, MAXE> parent;
        array<pair<int, int>, MAXE> type;

        sample_1PI_graph( L, n, k, gg, parent, type, Z, B, rngs[tid] );

        if( check_primitive )
        {
            // complete graph gg -> ggc
            graph ggc = gg;
            int nH = 2*ggc.E + ggc.n;
            array<int, MAXK> vert;

            for( int i = 0; i < ggc.n; i++ )
            {
                int h = nH + i;
                ggc.edges[ggc.E++] = make_pair(h, ggc.legs[i]);
                vert[i] = h;
            }

            ggc.vertices[ggc.V++] = vert;
            ggc.legs = {};

            if( !is_k4_D4_primitive( ggc ) )
                continue;

            Nprim++;
        }
/*
        cout << gg << endl;
        for( int e = 0; e < gg.E; e++ )
            cout << parent[e] << ", ";
        cout << endl;
*/
        fundamental_basis( matCycle, gg, type );

        array<double, MAXE> X;
        double log_psi_tr = sample_edge_lengths( X, D, k, gg, parent, type, rngs[tid] );

        double norm = exp( - log_psi_tr / L );

        double maxX = 0.;
        double phi = 0.;
        for( int e = 0; e < gg.E; e++ )
        {
//            cout << X[e] << ", ";
            double val = exp( X[e] );

            if( val > maxX )
                maxX = val;

            phi += val;

            Xvec[e] = val * norm;
        }
//        cout << endl;

        log_psi_tr = 0;

        laplace = matCycle.transpose() * Xvec.asDiagonal() * matCycle;

        ldlt.compute( laplace );

        if (ldlt.info() != Eigen::Success) {
            cerr << "LDLT decomposition failed" << endl;
        }

        double log_psi = ldlt.vectorD().array().log().sum();

        double w = omega( L, n, D, k );
        double amp = tgamma(w+1) * exp( (log_psi_tr - log_psi) * D/2.0 ) * pow( maxX / phi, w );

        if( !isfinite(amp) )
        {
            cout << "NOT FINITE " << amp << endl;

            cout << Xvec << endl;

            for (int i = 0; i < laplace.rows(); ++i) {
                for (int j = 0; j < laplace.cols(); ++j)
                    cout << laplace(i,j) << " ";
                cout << "\n";
            }
        }

        double amp_power = 1.0;
        for( int i = 0; i < M; i++ )
        {
            moment_sum[i][tid].add( amp_power );
            amp_power *= amp;
        }
    }

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;

    vector< xsum::xsum_large > moment( M );
    for( int i = 0; i < M; i++ )
        for( int j = 0; j < max_threads; j++ )
            moment[i].add( moment_sum[i][j] );

    double m1 = moment[1].round() / N;
    double m2 = moment[2].round() / N;  
    double m3 = moment[3].round() / N;
    double m4 = moment[4].round() / N;

    double avg = m1;
    double var = m2 - m1 * m1;
    double err = sqrt( var / N );

    double mu3 = m3 - 3 * m1 * m2 + 2 * m1 * m1 * m1;
    double mu4 = m4 - 4 * m1 * m3 + 6 * m1 * m1 * m2 - 3 * m1 * m1 * m1 * m1;

    double skew = mu3 / pow( var, 1.5 );

    double kur = mu4 / (var * var);

    json o;

    o["N"] = N;
    o["L"] = L;
    o["n"] = n;

    o["M1"] = m1;
    o["M2"] = m2;
    o["M3"] = m3;
    o["M4"] = m4;

    double leg_factor = tgamma(n+1);

    if( !check_primitive )
    {
        o["Tropical Amplitude"] = leg_factor*Z[L][n];
        o["Amplitude"] = (leg_factor*Z[L][n]*avg);
        o["Amplitude Err"] = (leg_factor*Z[L][n]*err);
        o["Relative Amplitude Err"] = (err/avg);
    }
    else
    {
        double prim_prop = (Nprim/(double)N);
        double prim_prop_var = Nprim/(double)N - prim_prop*prim_prop;
        double prim_prop_err = sqrt( prim_prop_var / N );

        o["Period"] = (leg_factor*B[L-1][n+2]/2*avg);
        o["Period Err"] = (leg_factor*B[L-1][n+2]/2*err);
        o["Rel Period Err"] = (err/avg);

        o["Pos Hepp"] = leg_factor*B[L-1][n+2]/2;
        o["Prim Hepp"] = leg_factor*(B[L-1][n+2]/2*prim_prop);
        o["Prim Hepp Err"] = leg_factor*(B[L-1][n+2]/2*prim_prop_err);
        o["Rel Hepp Err"] = (prim_prop_err/prim_prop);
    }

    o["Skewness"] = skew;
    o["Kurtosis"] = kur;
    o["Seed"] = seed;

    o["Seconds"] = elapsed.count();
    o["Threads"] = max_threads;
    
    cout << o.dump() << endl;
        
    return 0;
}

