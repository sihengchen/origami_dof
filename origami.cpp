#include "origami.hpp"
#include <iostream>
#include <math.h>
#include <vector>
#include <cassert>
#include <random>
#include <algorithm>
#include <fstream>
#include <omp.h>
#include <cstdio>
#include "/Users/Siheng/Downloads/SuiteSparse/include/SuiteSparseQR.hpp"

/*
Codes for the paper:
Chen, S., & Mahadevan, L. (2019). Rigidity percolation and geometric information
in floppy origami. Proceedings of the National Academy of Sciences, 116(17), 8119-8124.

For research and educational purpose only.

Please contact siheng_chen@g.harvard.edu for any questions about the codes.

origami::gen_DoF is the main function calculating the rank of the rigidity matrix (thus the DoF)
*/

origami::origami(long input_n_row_quad, long input_n_col_quad, long input_n_cst_coplanar, long input_start_i, long input_start_j){
  n_row_quad = input_n_row_quad;
  n_col_quad = input_n_col_quad;
  n_cst_coplanar = input_n_cst_coplanar;
  start_i = input_start_i;
  start_j = input_start_j;

  n_quad = n_row_quad * n_col_quad;
  n_row_node = n_row_quad+1;
  n_col_node = n_col_quad+1;
  n_node = n_row_node * n_col_node;
  n_cst_edges = n_quad * 2 + n_row_quad + n_col_quad + n_quad;
  row_num = 0;
  shuffle_done = 0;
  n_nodes_in_largest_cluster_square = 0;

  gen_random_cst_list();


}
void origami::add_A_entry(long r, long c, double x)
{
    ((long*)A->i)[A->nnz] = r;
    ((long*)A->j)[A->nnz] = c;
    ((double*)A->x)[A->nnz] = x;
    (A->nnz)++;
}

/*// Get the coordinates of the four corner nodes
int * genCoord(int i, int j, int * coord, double gammma, double theta, double L)
{
  int start_x_index = i/2, start_y_index = j/2;
  int x_odd_flag = i % 2, y_odd_flag = j % 2;
  double start_x_coord = 1;
  return coord;
}*/

// When constructing the rigidity matrix, only the nij is useful
//       ----3---->              /     /`.           both diagonal towards upper direction (towards large y)
//      /        /              /        ``
//     4(up)    2 (up)         5           6
//    /        /              /             ``
//   -----1--->              /               ``
double * origami::genVectors(int quad_i, int quad_j, double * vec_quad)
{
  // (i, j) is the quad at ith column (i along x) and jth row (j along y).
  // The first quad is 1, so it start with (1, 1)
  int x_odd_flag = (quad_i+start_i) % 2, y_odd_flag = (quad_j+start_j) % 2;
  vec_quad[12] = L;  vec_quad[13] = S;
  vec_quad[15] = -L; vec_quad[16] = S;
  if (y_odd_flag == 1) {
    vec_quad[3] = V; vec_quad[4] = S; vec_quad[5] = 0;
    vec_quad[9] = V; vec_quad[10] = S; vec_quad[11] = 0;
    vec_quad[12] += V;
    vec_quad[15] += V;
  }
  else {
    vec_quad[3] = -V; vec_quad[4] = S; vec_quad[5] = 0;
    vec_quad[9] = -V; vec_quad[10] = S; vec_quad[11] = 0;
    vec_quad[12] -= V;
    vec_quad[15] -= V;
  }
  if (x_odd_flag == 1) {
    vec_quad[0] = L; vec_quad[1] = 0; vec_quad[2] = H;
    vec_quad[6] = L; vec_quad[7] = 0; vec_quad[8] = H;
                                      vec_quad[14] = H;
                                      vec_quad[17] = -H;
  }
  else {
    vec_quad[0] = L; vec_quad[1] = 0; vec_quad[2] = -H;
    vec_quad[6] = L; vec_quad[7] = 0; vec_quad[8] = -H;
                                      vec_quad[14] = -H;
                                      vec_quad[17] = H;
  }
  return vec_quad;
}

// Normalize a vector. Return the norm.
double origami::normalize_vec(std::vector<double> vec_to_norm, int vec_size)
{
  double square_sum = 0;
  for (int i=0;i<vec_size;i++){
    square_sum += vec_to_norm[i]*vec_to_norm[i];
  }
  return sqrt(square_sum);
}

//For a given quad, add all the constraints
// i along the x, j along y.
// i from 0 to n_col_quad (number of columns); j from 0 to n_row_quad
// Calculation: j*n_col_node+i is the current
void origami::gen_edge_Constraint_Triplets(double * vec_quad)
{

  std::vector<double> vec_to_norm = {0,0,0,0,0,0};
  double vec_norm;
  for (long i=0;i<n_col_quad;i++){
    for (long j=0;j<n_row_quad;j++){
      long f = j * n_col_node + i; // First node at the bottom left corner;
      //std::cout << "f" << f << std::endl;
      //i, j start from (0,0), so we use (i+1, j+1) for genVectors
      vec_quad = genVectors(i+1, j+1, vec_quad);
      vec_to_norm = {vec_quad[0], vec_quad[1], vec_quad[2], -vec_quad[0], -vec_quad[1], -vec_quad[2]};
      vec_norm = normalize_vec(vec_to_norm, 6);
      //std::cout << "Norm" << vec_norm << std::endl;
      add_A_entry(row_num, 3*f,   vec_quad[0]/vec_norm);
      add_A_entry(row_num, 3*f+1, vec_quad[1]/vec_norm);
      add_A_entry(row_num, 3*f+2, vec_quad[2]/vec_norm);
      add_A_entry(row_num, 3*(f+1),   -vec_quad[0]/vec_norm);
      add_A_entry(row_num, 3*(f+1)+1, -vec_quad[1]/vec_norm);
      add_A_entry(row_num, 3*(f+1)+2, -vec_quad[2]/vec_norm);
      row_num++;

      vec_to_norm = {-vec_quad[9], -vec_quad[10], -vec_quad[11], vec_quad[9], vec_quad[10], vec_quad[11]};
      vec_norm = normalize_vec(vec_to_norm, 6);
      add_A_entry(row_num, 3*(f+n_col_node),   -vec_quad[9]/vec_norm);
      add_A_entry(row_num, 3*(f+n_col_node)+1, -vec_quad[10]/vec_norm);
      add_A_entry(row_num, 3*(f+n_col_node)+2, -vec_quad[11]/vec_norm);
      add_A_entry(row_num, 3*(f),   vec_quad[9]/vec_norm);
      add_A_entry(row_num, 3*(f)+1, vec_quad[10]/vec_norm);
      add_A_entry(row_num, 3*(f)+2, vec_quad[11]/vec_norm);
      row_num++;

      if (i == n_col_quad-1) {
        vec_to_norm = {vec_quad[3], vec_quad[4], vec_quad[5], -vec_quad[3], -vec_quad[4], -vec_quad[5]};
        vec_norm = normalize_vec(vec_to_norm, 6);
        add_A_entry(row_num, 3*(f+1),   vec_quad[3]/vec_norm);
        add_A_entry(row_num, 3*(f+1)+1, vec_quad[4]/vec_norm);
        add_A_entry(row_num, 3*(f+1)+2, vec_quad[5]/vec_norm);
        add_A_entry(row_num, 3*(f+n_col_node+1),   -vec_quad[3]/vec_norm);
        add_A_entry(row_num, 3*(f+n_col_node+1)+1, -vec_quad[4]/vec_norm);
        add_A_entry(row_num, 3*(f+n_col_node+1)+2, -vec_quad[5]/vec_norm);
        row_num++;
      }
      if (j == n_row_quad-1){
        vec_to_norm = {-vec_quad[6], -vec_quad[7], -vec_quad[8], vec_quad[6], vec_quad[7], vec_quad[8]};
        vec_norm = normalize_vec(vec_to_norm, 6);
        add_A_entry(row_num, 3*(f+n_col_node+1),   -vec_quad[6]/vec_norm);
        add_A_entry(row_num, 3*(f+n_col_node+1)+1, -vec_quad[7]/vec_norm);
        add_A_entry(row_num, 3*(f+n_col_node+1)+2, -vec_quad[8]/vec_norm);
        add_A_entry(row_num, 3*(f+n_col_node),   vec_quad[6]/vec_norm);
        add_A_entry(row_num, 3*(f+n_col_node)+1, vec_quad[7]/vec_norm);
        add_A_entry(row_num, 3*(f+n_col_node)+2, vec_quad[8]/vec_norm);
        row_num++;
      }

      vec_to_norm = {vec_quad[12],vec_quad[13],vec_quad[14]};
      double norm5 = normalize_vec(vec_to_norm, 3);
      vec_to_norm = {vec_quad[15],vec_quad[16],vec_quad[17]};
      double norm6 = normalize_vec(vec_to_norm, 3);

      double r = ((double) rand() / (RAND_MAX));
      //std::cout << "RAND" << r << std::endl;
      if ((norm5 < norm6 && CURR_FOLD_OPTION == FOLD_CLOSE_PAIRS) || (r < 0.5 && CURR_FOLD_OPTION == FOLD_RANDOM)) {
        vec_to_norm = {vec_quad[12], vec_quad[13], vec_quad[14], -vec_quad[12], -vec_quad[13], -vec_quad[14]};
        vec_norm = normalize_vec(vec_to_norm, 6);
        add_A_entry(row_num, 3*(f),   vec_quad[12]/vec_norm);
        add_A_entry(row_num, 3*(f)+1, vec_quad[13]/vec_norm);
        add_A_entry(row_num, 3*(f)+2, vec_quad[14]/vec_norm);
        add_A_entry(row_num, 3*(f+n_col_node+1),   -vec_quad[12]/vec_norm);
        add_A_entry(row_num, 3*(f+n_col_node+1)+1, -vec_quad[13]/vec_norm);
        add_A_entry(row_num, 3*(f+n_col_node+1)+2, -vec_quad[14]/vec_norm);
        row_num++;
      }
      else {
        vec_to_norm = {vec_quad[15], vec_quad[16], vec_quad[17], -vec_quad[15], -vec_quad[16], -vec_quad[17]};
        vec_norm = normalize_vec(vec_to_norm, 6);
        add_A_entry(row_num, 3*(f+1),   vec_quad[15]/vec_norm);
        add_A_entry(row_num, 3*(f+1)+1, vec_quad[16]/vec_norm);
        add_A_entry(row_num, 3*(f+1)+2, vec_quad[17]/vec_norm);
        add_A_entry(row_num, 3*(f+n_col_node),   -vec_quad[15]/vec_norm);
        add_A_entry(row_num, 3*(f+n_col_node)+1, -vec_quad[16]/vec_norm);
        add_A_entry(row_num, 3*(f+n_col_node)+2, -vec_quad[17]/vec_norm);
        row_num++;
      }

    }
  }

  //return row_num;
}


/*
points:
4  3
1  2

a3-a1, b3-b1, c3-c1
a2-a1, b2-b1, c2-c1
a4-a1, b4-b1, c4-c1

(a3-a1)((b2-b1)(c4-c1)-(b4-b1)(c2-c1)) +
(b3-b1)((a4-a1)(c2-c1)-(a2-a1)(c4-c1)) +
(c3-c1)((a2-a1)(b4-b1)-(a4-a1)(b2-b1))

a2-a1, b2-b1, c2-c1: vec_quad[0,1,2];
a3-a1, b3-b1, c3-c1: vec_quad[12,13,14];
a4-a1, b4-b1, c4-c1: vec_quad[9,10,11]
*/
void origami::gen_Coplanar_Constraint_Triplets(std::vector<int> constraint_list, double * vec_quad) {
  std::vector<double> vec_to_norm = {0,0,0,0,0,0,0,0,0,0,0,0};
  double vec_norm;
  for (long t=0; t<n_cst_coplanar; t++){
    long quad_num = constraint_list[t]; // # of the quad that has coplanar constraint
    long quad_j = floor(double(quad_num) / double(n_col_quad)); // ith column (i is along x)
    long quad_i = quad_num % n_col_quad; // jth row (j is along y)
    vec_quad = genVectors(quad_i+1, quad_j+1, vec_quad);
    long f = quad_j * n_col_node + quad_i;
    double a2a1 = vec_quad[0], b2b1=vec_quad[1], c2c1=vec_quad[2];
    double a3a1 = vec_quad[12], b3b1=vec_quad[13], c3c1=vec_quad[14];
    double a4a1 = vec_quad[9], b4b1=vec_quad[10], c4c1=vec_quad[11];
    /*for (int tt = 0;tt<18;tt++)
    {
      std::cout << tt << " " << vec_quad[tt] << std::endl;
    }*/
    //Derivative cst / a1,b1,c1
    vec_to_norm[0] =  b4b1*c2c1 + b2b1*c3c1 - b4b1*c3c1 - b2b1*c4c1 + b3b1*c4c1 - b3b1*c2c1;
    vec_to_norm[1] = -a4a1*c2c1 + a4a1*c3c1 - a2a1*c3c1 + a3a1*c2c1 - a3a1*c4c1 + a2a1*c4c1;
    vec_to_norm[2] =  a4a1*b2b1 + a2a1*b3b1 - a4a1*b3b1 - a2a1*b4b1 + a3a1*b4b1 - a3a1*b2b1;
    vec_to_norm[3] =  b4b1*c3c1 - b3b1*c4c1;
    vec_to_norm[4] = -a4a1*c3c1 + a3a1*c4c1;
    vec_to_norm[5] =  a4a1*b3b1 - a3a1*b4b1;
    vec_to_norm[6] = -b4b1*c2c1 + b2b1*c4c1;
    vec_to_norm[7] =  a4a1*c2c1 - a2a1*c4c1;
    vec_to_norm[8] = -a4a1*b2b1 + a2a1*b4b1;
    vec_to_norm[9] =  b3b1*c2c1 - b2b1*c3c1;
    vec_to_norm[10] =-a3a1*c2c1 + a2a1*c3c1;
    vec_to_norm[11] = a3a1*b2b1 - a2a1*b3b1;

    vec_norm = normalize_vec(vec_to_norm, 12);

    add_A_entry(row_num, 3*(f),   vec_to_norm[0]/vec_norm);
    add_A_entry(row_num, 3*(f)+1, vec_to_norm[1]/vec_norm);
    add_A_entry(row_num, 3*(f)+2, vec_to_norm[2]/vec_norm);
    add_A_entry(row_num, 3*(f+1),   vec_to_norm[3]/vec_norm);
    add_A_entry(row_num, 3*(f+1)+1, vec_to_norm[4]/vec_norm);
    add_A_entry(row_num, 3*(f+1)+2, vec_to_norm[5]/vec_norm);
    add_A_entry(row_num, 3*(f+n_col_node+1),   vec_to_norm[6]/vec_norm);
    add_A_entry(row_num, 3*(f+n_col_node+1)+1, vec_to_norm[7]/vec_norm);
    add_A_entry(row_num, 3*(f+n_col_node+1)+2, vec_to_norm[8]/vec_norm);
    add_A_entry(row_num, 3*(f+n_col_node),   vec_to_norm[9]/vec_norm);
    add_A_entry(row_num, 3*(f+n_col_node)+1, vec_to_norm[10]/vec_norm);
    add_A_entry(row_num, 3*(f+n_col_node)+2, vec_to_norm[11]/vec_norm);
    row_num++;
    //a1, a2, a3
  }
  //return row_num;

}

// Generate DoF of a given origami
SuiteSparse_long origami::gen_DoF(std::vector <int> constraint_list)
{
  cholmod_common Common, * com;
  com = &Common;
  cholmod_l_start(com);


  A = cholmod_l_allocate_triplet(n_cst_edges+n_cst_coplanar, 3*n_node, n_cst_edges * 6 +n_cst_coplanar * 12, 0, CHOLMOD_REAL, com);

  double * vec_quad = new double[18];
  //Generate all edge constraints
  gen_edge_Constraint_Triplets(vec_quad);
  //Generate all coplanar constraints
  gen_Coplanar_Constraint_Triplets(constraint_list, vec_quad);
  delete [] vec_quad;

  // Calculate Rank
  SuiteSparse_long econ = 0;//min(rgd_Matrix->nrow, rgd_Matrix->ncol);
  //cholmod_l_print_triplet(A, "A", com);
  cholmod_sparse * rgd_Matrix = cholmod_l_triplet_to_sparse(A, n_cst_edges * 6 +n_cst_coplanar * 12, com);
  SuiteSparse_long rank = SuiteSparseQR <double> (SPQR_ORDERING_GIVEN, SPQR_DEFAULT_TOL, econ, rgd_Matrix, NULL, NULL, com);

  // Cleaning
  cholmod_l_free_triplet(&A, com);
  cholmod_l_free_sparse(&rgd_Matrix, com);
  cholmod_l_finish(com);
  SuiteSparse_long dof = 3 * n_node - rank - 6;
  return dof;
}


// Randomly generate n_cst_coplanar coplanar constraints.
// Result is updated in constraint_all and constraint_all_extended
// constraint_all_extended will be updated by find_Rigid_Quads() later
void origami::gen_random_cst_list() {
  std::vector<int> constraint_list(n_quad);
  std::iota (std::begin(constraint_list), std::end(constraint_list), 0);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::shuffle(constraint_list.begin(), constraint_list.end(), gen);
  constraint_list.resize(n_cst_coplanar);

  constraint_all = constraint_list;
  //Make sure both are sorted initially
  std::sort (constraint_all.begin(),constraint_all.end());
  constraint_all_extended = constraint_all;
  shuffle_done = 1;
}
