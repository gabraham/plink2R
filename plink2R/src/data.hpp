#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigen>
#include <Eigen/SVD>
#include <Eigen/Sparse>

#define IMPUTE_NONE 0
#define IMPUTE_AVG 1
#define IMPUTE_RANDOM 2

//#include <boost/filesystem.hpp>
//#include <boost/iostreams/device/mapped_file.hpp>

//#include "util.hpp"

#define PACK_DENSITY 4
#define PLINK_NA 3

#define PHENO_BINARY_12 0
#define PHENO_CONTINUOUS 1
#define PLINK_PHENO_MISSING -9

// The BED file magic numbers
#define PLINK_OFFSET 3

#define COVAR_ACTION_TRAIN_TEST 0
#define COVAR_ACTION_TRAIN_ONLY 1

#define COVAR_ACTION_TRAIN_TEST_STR "traintest"
#define COVAR_ACTION_TRAIN_ONLY_STR "train"

/* 3 is 11 in binary, we need a 2 bit mask for each of the 4 positions */
#define MASK0 3	  /* 3 << 2 * 0 */
#define MASK1 12  /* 3 << 2 * 1 */
#define MASK2 48  /* 3 << 2 * 2 */
#define MASK3 192 /* 3 << 2 * 3 */

#define BUFSIZE 100
#define DATA_MODE_TRAIN 1
#define DATA_MODE_TEST 2

using namespace Eigen;

typedef Array<bool, Dynamic, Dynamic> ArrayXXb;
typedef Array<bool, Dynamic, 1> ArrayXb;

class Data {
   public:
      
      MatrixXd X, X2, Y;
      MatrixXd Xtrain, Xtest; // only used if data is small enough
      MatrixXd X2train, X2test;
      MatrixXd Ytrain, Ytest;
      unsigned int N, p, K;
      unsigned long long len, filesize;
      unsigned int np, nsnps, ncovar;
      ArrayXb mask;
      char *geno_filename;
      //boost::iostreams::mapped_file_source geno_fin;
      int geno_fin_fd;
      unsigned char* data;
      //std::vector<unsigned int> covar_ignore_pred_idx;
      std::vector<unsigned int> covar_actions;
      unsigned int mode;
      VectorXd ones, zeros;
      VectorXd geno;
      VectorXd *geno_ptr;
      const char *bedfile, *famfile, *bimfile;
      bool verbose;
      
      Data(const char* bedfile, const char* famfile, bool verbose);
      ~Data();
      void read_bed(int impute);
      MatrixXd read_plink_pheno(const char *filename, unsigned int firstcol);
      void load_snp_double(unsigned int j, double *geno);
      VectorXd load_snp(unsigned int j);

      std::string tolower(const std::string& v);
      VectorXd get_coordinate(unsigned int j);
      VectorXd get_snp(unsigned int j);
};


