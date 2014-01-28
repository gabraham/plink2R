
#include "data.hpp"

#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <stdexcept>

Data::Data(const char* bedfile, const char* famfile, bool verbose)
{
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

   //mmap_bed(bedfile);
}

Data::~Data()
{
   //geno_fin.close();
   munmap(data, filesize);
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
void Data::read_bed(bool impute)
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

   if(impute)
   {
      for(unsigned int i = 0 ; i < nsnps ; i++)
      {
         // read raw genotypes
         in.read((char*)tmp, sizeof(char) * np);

         // decode the genotypes
         decode_plink(tmp2, tmp, np);

         // Compute average per SNP, excluding missing values
         avg[i] = 0;
         unsigned int ngood = 0;
         for(unsigned int j = 0 ; j < N ; j++)
         {
            double s = (double)tmp2[j];
            if(s != PLINK_NA)
            {
               avg[i] += s;
               ngood++;
            }
         }
         avg[i] /= ngood;

         // Impute using average per SNP
         for(unsigned int j = 0 ; j < N ; j++)
         {
            double s = (double)tmp2[j];
            if(s != PLINK_NA)
               tmp3(j) = s;
            else
               tmp3(j) = avg[i];
         }

         X.col(i) = tmp3;
      }
   }
   else
   {
      for(unsigned int i = 0 ; i < nsnps ; i++)
      {
         // read raw genotypes
         in.read((char*)tmp, sizeof(char) * np);

         // decode the genotypes
         decode_plink(tmp2, tmp, np);

         for(unsigned int j = 0 ; j < N ; j++)
	    tmp3(j) = (double)tmp2[j];
         X.col(i) = tmp3;
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

void Data::mmap_bed(const char *filename)
{
  // boost::iostreams::mapped_file_params params;
   //params.path = std::string(filename, strlen(filename));
   //params.offset = 0;
   //params.length = -1;
   //geno_fin = boost::iostreams::mapped_file_source(params);

   geno_fin_fd = open(filename, O_RDONLY);
   if(geno_fin_fd == -1)
      throw std::runtime_error(
	 "can't open file for reading");

   struct stat64 stat_buf;
   int rc = stat64(filename, &stat_buf);
   filesize = stat_buf.st_size;

   if(verbose)
      std::cout << "filesize: " << filesize << std::endl;
   len = filesize - PLINK_OFFSET;
   data = (unsigned char*)mmap(NULL, filesize, PROT_READ, MAP_SHARED, geno_fin_fd, 0);

   if(N == 0)
      throw std::runtime_error(
	 "haven't read a FAM/PHENO file so don't know what sample size is");

   if(verbose)
      std::cout << filename << " len: " << len << " bytes" << std::endl;
   np = (unsigned int)ceil((double)N / PACK_DENSITY);
   nsnps = (unsigned int)(len / np);

   if(verbose)
   {
      std::cout << filename << " np: " << np << std::endl;
      std::cout << filename << " nsnps: " << nsnps << std::endl;
   }
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
   double mean = sum / k;
   double sum2 = 0, v;
   for(i = N - 1 ; i != -1 ; --i)
   {
      v = geno[i];
      if(v != PLINK_NA)
      {
	 v -= mean;
	 sum2 += v * v;
      }
   }
   double sd = sqrt(sum2 / (ngood - 1));
   double mean_sd = mean / sd;

   if(ngood == N)
   {
      for(i = N - 1 ; i != -1 ; --i)
	 geno[i] = (geno[i] - mean) / sd;
   }
   else
   {
      // Impute missing genotypes and standardise to zero-mean unit-variance
      for(i = N - 1 ; i != -1 ; --i)
      {
         v = geno[i];
         if(v == PLINK_NA)
            geno[i] = mean_sd;
         else
            geno[i] = (v - mean) / sd;
      }
   }
}

