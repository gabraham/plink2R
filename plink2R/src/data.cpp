
#include "data.hpp"

#include <fcntl.h>
#include <stdexcept>

Data::Data(const char* bedfile, const char* famfile, bool verbose)
{
   srand48(time(NULL));
   N = 0;
   p = 0;
   K = 0;
   nsnps = 0;
   ncovar = 0;
   this->bedfile = bedfile;
   this->famfile = famfile;
   this->verbose = verbose;
   if(verbose)
      std::cout << "bedfile: " << bedfile 
	 << " famfile: " << famfile << std::endl;
   this->Y = read_plink_pheno(famfile, 6);
   this->mask = VectorXd::Ones(N).array() == 1.0;
}

Data::~Data()
{
}

/* 
 *                   plink BED           sparsnp
 * minor homozyous:  00 => numeric 0     10 => numeric 2
 * heterozygous:     10 => numeric 2     01 => numeric 1
 * major homozygous: 11 => numeric 3     00 => numeric 0
 * missing:          01 => numeric 1     11 => numeric 3
 *
 *
 * http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml says,
 * The bytes in plink are read backwards HGFEDCBA, not GHEFCDAB, but we read
 * them forwards as a character (a proper byte)
 *
 * By default, plink usage dosage of the *major* allele, since allele A1 is
 * usually the minor allele and the code "1" refers to the second allele A2,
 * so that "11" is A2/A2 or major/major.
 *
 * We always use minor allele dosage, to be consistent with the output from
 * plink --recodeA which used minor allele dosage by default.
 *
 * out: array of genotypes
 * in: array of packed genotypes (bytes)
 * n: number of bytes in input
 * 
 */
void decode_plink(unsigned char *out,
      const unsigned char *in, const unsigned int n)
{
   unsigned int i, k;
   unsigned char tmp, geno;
   unsigned int a1, a2;

   for(i = 0 ; i < n ; ++i)
   {
      tmp = in[i];
      k = PACK_DENSITY * i;
      
      /* geno is interpreted as a char, however a1 and a2 are bits for allele 1 and
       * allele 2. The final genotype is the sum of the alleles, except for 01
       * which denotes missing.
       */
      geno = (tmp & MASK0);
      a1 = !(geno & 1);
      a2 = !(geno >> 1);
      out[k] = (geno == 1) ? 3 : a1 + a2;
      k++;

      geno = (tmp & MASK1) >> 2; 
      a1 = !(geno & 1);
      a2 = !(geno >> 1);
      out[k] = (geno == 1) ? 3 : a1 + a2;
      k++;

      geno = (tmp & MASK2) >> 4; 
      a1 = !(geno & 1);
      a2 = !(geno >> 1);
      out[k] = (geno == 1) ? 3 : a1 + a2;
      k++;

      geno = (tmp & MASK3) >> 6; 
      a1 = !(geno & 1);
      a2 = !(geno >> 1);
      out[k] = (geno == 1) ? 3 : a1 + a2;
   }
}

// Expects PLINK BED in SNP-major format
void Data::read_bed(int impute)
{
   if(verbose)
      std::cout << ">>> Reading BED file '" << this->bedfile << "'" << std::endl;
   std::ifstream in(this->bedfile, std::ios::in | std::ios::binary);

   if(!in)
   {
      std::cerr << "[read_bed] Error reading file " 
	 << this->bedfile << std::endl;
      throw std::runtime_error("io error");
   }

   in.seekg(0, std::ifstream::end);

   // file size in bytes, ignoring first 3 bytes (2byte magic number + 1byte mode)
   len = (unsigned int)in.tellg() - 3;

   // size of packed data, in bytes, per SNP
   np = (unsigned int)ceil((double)N / PACK_DENSITY);
   nsnps = len / np;
   in.seekg(3, std::ifstream::beg);

   unsigned char* tmp = new unsigned char[np];

   // Allocate more than the sample size since data must take up whole bytes
   unsigned char* tmp2 = new unsigned char[np * PACK_DENSITY];
   X = MatrixXd(N, nsnps);

   if(verbose)
      std::cout << ">>> Detected BED file: " << this->bedfile <<
      " with " << len << " bytes, " << N << " samples, " << nsnps 
      << " SNPs." << std::endl;
   VectorXd tmp3(N);

   double* avg = new double[nsnps]; 

   // The Rcpp code in read_plink will take care of converting PLINK_NA to
   // NA_REAL
   if(impute == IMPUTE_NONE)
   {
      if(verbose)
	 std::cout << "impute: IMPUTE_NONE" << std::endl;
   
      for(unsigned int j = 0 ; j < nsnps ; j++)
      {
         // read raw genotypes
         in.read((char*)tmp, sizeof(char) * np);

         // decode the genotypes
         decode_plink(tmp2, tmp, np);

         for(unsigned int i = 0 ; i < N ; i++)
            tmp3(i) = (double)tmp2[i];
         X.col(j) = tmp3;
      }
   }
   // impute by average
   else if(impute == IMPUTE_AVG)
   {
      if(verbose)
	 std::cout << "impute: IMPUTE_AVG" << std::endl;

      for(unsigned int j = 0 ; j < nsnps ; j++)
      {
         // read raw genotypes
         in.read((char*)tmp, sizeof(char) * np);

         // decode the genotypes
         decode_plink(tmp2, tmp, np);

         // Compute average per SNP, excluding missing values
         avg[j] = 0;
         unsigned int ngood = 0;
         for(unsigned int i = 0 ; i < N ; i++)
         {
            double s = (double)tmp2[i];
            if(s != PLINK_NA)
            {
               avg[j] += s;
               ngood++;
            }
         }
         avg[j] /= ngood;

         // Impute using average per SNP
         for(unsigned int i = 0 ; i < N ; i++)
         {
            double s = (double)tmp2[i];
            if(s != PLINK_NA)
               tmp3(i) = s;
            else
               tmp3(i) = avg[j];
         }

         X.col(j) = tmp3;
      }
   }
   // impute by sampling genotypes
   else if(impute == IMPUTE_RANDOM)
   {
      if(verbose)
	 std::cout << "impute: IMPUTE_RANDOM" << std::endl;

      for(unsigned int j = 0 ; j < nsnps ; j++)
      {
         // read raw genotypes
         in.read((char*)tmp, sizeof(char) * np);

         // decode the genotypes
         decode_plink(tmp2, tmp, np);

	 double proportions[3];
	 unsigned int ngood = 0, sum0 = 0, sum1 = 0, sum2 = 0;
	 char x;

	 // first pass to get genotype proportions
         for(unsigned int i = 0 ; i < N ; i++)
	 {
	    x = tmp2[i];
	    if(x != PLINK_NA)
	    {
	       sum0 += (x == 0);
	       sum1 += (x == 1);
	       sum2 += (x == 2);
	       ngood++;
	       tmp3(i) = (double)x;
	    }
	 }

	 if(ngood < N)
	 {
	    proportions[0] = (double)sum0 / ngood;
	    proportions[1] = (double)sum1 / ngood;
	    proportions[2] = (double)sum2 / ngood;
	    
	    double cumsum[3] = {
	       proportions[0],
	       proportions[0] + proportions[1],
	       1
	    };

	    // second pass to sample the missing genotypes
            for(unsigned int i = 0 ; i < N ; i++)
	    {
	       x = tmp2[i];
	       if(x == PLINK_NA)
	       {
	          double r = drand48();
	          if(r < cumsum[0])
	             tmp3(i) = 0.0;
	          else if(r < cumsum[1])
	             tmp3(i) = 1.0;
	          else
	             tmp3(i) = 2.0;
	       }
	    }
	 }

         X.col(j) = tmp3;
      }
   }

   p = X.cols();

   delete[] tmp;
   delete[] tmp2;
   delete[] avg;

   in.close();
}

// Reads PLINK phenotype files:
// FID IID pheno1 pheno2 ...
// Need to be able to read continuous phenotypes
// 
// firstcol is one-based, 3 for pheno file, 6 for FAM file (ignoring gender),
// 5 for FAM file (with gender)
MatrixXd Data::read_plink_pheno(const char *filename, unsigned int firstcol)
{
   std::ifstream in(filename, std::ios::in);

   if(!in)
   {
      std::cerr << "[read_plink_pheno] Error reading file " 
	 << filename << std::endl;
      throw std::string("io error");
   }
   std::vector<std::string> lines;
   
   while(in)
   {
      std::string line;
      std::getline(in, line);
      if(!in.eof())
	 lines.push_back(line);
   }

   if(verbose)
      std::cout << ">>> Detected pheno file " <<
	 filename << ", " << lines.size() << " samples";

   in.close();

   unsigned int numtok = 0, numfields;

   MatrixXd Z(0, 0);

   for(unsigned int i = 0 ; i < lines.size() ; i++)
   {
      std::stringstream ss(lines[i]);
      std::string s;
      std::vector<std::string> tokens;

      while(ss >> s)
	 tokens.push_back(s);

      numtok = tokens.size();
      numfields = numtok - firstcol + 1;
      if(i == 0)
	 Z.resize(lines.size(), numfields);

      VectorXd y(numfields);
      for(unsigned int j = 0 ; j < numfields ; j++)
	 y(j) = std::atof(tokens[j + firstcol - 1].c_str());
      Z.row(i) = y;
   }

   if(verbose)
      std::cout << ", " << Z.cols() << " columns (ex. FAM+INDIV IDs)" << std::endl;
   
   N = Z.rows();

   return Z;
}

std::string Data::tolower(const std::string& v)
{
   std::string r = v;
   std::transform(r.begin(), r.end(), r.begin(), ::tolower);
   return r;
}

// Assumes data is SNP-major ordered
VectorXd Data::get_snp(unsigned int j)
{
   geno = load_snp(j);
   return geno;
}

VectorXd Data::load_snp(unsigned int j)
{
   double geno_dat[N];
   load_snp_double(j, geno_dat);
   VectorXd geno = Map<VectorXd>(geno_dat, N, 1);
   return geno;
}

// Read a SNP from disk, for the correct set (training or test),
// impute missing and standardise
void Data::load_snp_double(unsigned int j, double *geno)
{
   unsigned int i;

   // Not found in cache, load from disk
   unsigned char* tmp2 = new unsigned char[np * PACK_DENSITY];
   decode_plink(tmp2, data + PLINK_OFFSET + j * np, np);

   // Compute average over non-missing genotypes
   unsigned int k = 0, ngood = 0;
   double sum = 0;
   for(i = 0 ; i < N ; i++)
   {
      if(mask(i))
      {
	 double v = (double)tmp2[i];
         geno[k] = v;
	 if(v != PLINK_NA)
	 {
	    ngood++;
	    sum += v;
	 }
	 k++;
      }
   }

   delete[] tmp2;

   // Compute standard dev over non-missing genotypes
   double mean = sum / ngood;
   double sum2 = 0, v;
   for(i = 0 ; i < N ; i++)
   {
      v = geno[i];
      if(v != PLINK_NA)
      {
	 v -= mean;
	 sum2 += v * v;
      }
   }
   double sd = sqrt(sum2 / (ngood - 1));

   if(ngood == N)
   {
      for(i = 0 ; i < N ; i++)
	 geno[i] = (geno[i] - mean) / sd;
   }
   else
   {
      // Impute missing genotypes and standardise to zero-mean unit-variance
      for(i = 0 ; i < N ; i++)
      {
         v = geno[i];
         if(v == PLINK_NA)
            geno[i] = 0;
         else
            geno[i] = (v - mean) / sd;
      }
   }
}

