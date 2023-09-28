#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define map(x)  ( (x) > 0 ? 1 : -1 )
#define demap(x) ( (x) > 0 ? 1 : 0 )

const short 
        g = 0b10011, // {11,15} polynomial
        taps = 0b1100,
        n = 15,
        k = 11; 

int syndromes[n];
int codewords[1 << k];

void awgn();
void hamming_encoder();
void hamming_decoder_hard();
void hamming_decoder_soft();
void generate_syndromes();
void generate_codewords();

int main(int argc, char *argv[])
{
   if(argc != 3)
   {
      printf("\nUsage: BPSK Eb/No(dB) data-block-size \n\n");
      exit(1);
   }

   float EbNodB, EcNodB, EcNo, sigma, sigsqrd, R=(float)k/n;  // change to R = (float)k/n if add coding
   int i, blk_errs=0, flag;
   long blk_cnt=0, bit_errs=0;
   int N, blk_errs_max = 150;   // 150 block errors should provide sufficient reliability
                                // in the measurement of error probabilities

   srand(time(NULL));  // choose random seed for rand()

   EbNodB = atof(argv[1]);
   N = atoi(argv[2]);      // for convenience, N should be even

   int Nc = ceil((float) N/k) * n;
   // printf("%d\n", Nc);
   float x[N], r[Nc];
   uint8_t y[N], y2[N];
   uint8_t m[N], c[Nc];

   EcNodB = EbNodB + 10*log10(R);  // general formula *if* there was coding
   EcNo = pow(10.0, EcNodB/10.0);
   sigsqrd = 0.5/EcNo; // var = No/2, Ec=1
   sigma = sqrt(sigsqrd);

   generate_syndromes();
   generate_codewords();

   while(blk_errs < blk_errs_max)
   {
      memset(m, 0, N);
      memset(c, 0, Nc);
      memset(r, 0, Nc);
      memset(y, 0, N);
      blk_cnt++;
      // get a block of data
      for(i=0; i<N; i++) m[i] = (rand()>>5) & 0x1;
   

      hamming_encoder(m, c, N);
      // send over additive white Gaussian noise (AWGN) channel (which also
      // converts (0, 1) to (-1, +1))
      awgn(sigma, c, r, Nc);

      // hamming_decoder_hard(r, y, Nc);
      hamming_decoder_soft(r, y, Nc);

      /*
         Debug printing
      */

      // printf("\nm: ");
      // for (int i = 0; i < N; i++) printf("%d", m[i]);
      // printf("\nc: ");

      // for (int i = 0; i < Nc; i++) printf("%d", c[i]);
      // printf("\nr: ");

      // for (int i = 0; i < Nc; i++) printf("%d", demap(r[i]));
      // printf("\ny: ");

      // for (int i = 0; i < N; i++) printf("%d", y[i]);
      // printf("\n");
      // printf("\n");

      // getchar();

      // count bit errors and block errors
      flag = 0;
      for(i=0; i<N; i++)
         if(demap(y[i]) != m[i]) { bit_errs++; flag = 1; }
      if(flag) blk_errs++; 
   }

   printf("bit-error rate = %10.2e   block-error rate = %10.2e \n",
             (float)bit_errs/(N*blk_cnt), (float)blk_errs/blk_cnt);

   // NOTE: When the bit-error rate Pb is Pb << 1 (usual case),
   // then block-error rate Pblk ~= N*Pb. The reason is that 
   // Pblk = 1 - (1 - Pb)^N ~= 1 - (1 - N*Pb) = N*Pb.

   return 0;
}

// marsaglia polar method
void awgn (sigma, x, y, len)
float sigma, y[];
uint8_t x[];
int len;
{
	int k;
	double u, v, s, a, b;
	
	for(k = 0; k < len; k+=2)
	{
            do
            {
		        u = (rand()/((double)RAND_MAX)) * 2.0 - 1.0;
		        v = (rand()/((double)RAND_MAX)) * 2.0 - 1.0;
                s = u*u + v*v;
            } while(s >= 1.0 || s == 0.0);
            s = sqrt(-2.0*log(s)/s);
	    y[k] = map(x[k]) + sigma*u*s;
	    y[k+1] = map(x[k+1]) + sigma*v*s;
	}
}

void hamming_encoder(m, c, message_len) 
uint8_t m[], c[];
int message_len;
{
   
    for (int idx_c = 0, idx_m = 0; idx_m < ceil((float) message_len/k)*k; idx_c += n, idx_m += k) {
        int enc_reg = 0;
        for (int i = 0; i < n; i++) {
            // for (int i = 3; i >= 0; i--) printf("%d", (enc_reg >> i) & 0b1);
            // printf("\n");
            uint8_t parity = enc_reg & 0b1;
            uint8_t mb = idx_m + i >= message_len ? 0 : m[idx_m + i]; // zero pad
            uint8_t f = parity ^ mb;
            if (i >= k) f = 0;
            enc_reg >>= 1;
            enc_reg ^= (f ? taps : 0);
            c[idx_c + i] = i < k ? mb : parity;
            // for(int i = 0; i < k; i++) printf("%d", m[idx_m + i]);
            // printf("\n");
            // for(int i = 0; i < n; i++) printf("%d", c[idx_c + i]);
            // printf("\n");

            // for (int i = 3; i >= 0; i--) printf("%d", (enc_reg >> i) & 0b1);
            // printf("\n\n");
            
            // getchar();
        }
    }
}

void hamming_decoder_hard(r, y, message_len) 
float r[];
uint8_t y[];
int message_len;
{
   for (int idx = 0, idx_y = 0; idx < message_len; idx += n, idx_y += k) {
        int dec_reg = 0;
        for (int i = 0; i < n; i++) {
            uint8_t parity = dec_reg & 0b1;
            uint8_t mb = demap(r[idx + i]);
            uint8_t f = parity ^ mb;
            dec_reg >>= 1;
            dec_reg ^= (f ? taps : 0);
            if (i < k) y[idx_y + i] = mb;
        }
        // decoder register holds remainder
        // nonzero remainder = error occurred
        if (dec_reg) for (int s = 0; s<n; s++) {
            if (syndromes[s] == dec_reg) {
               y[idx_y + s] ^= 1;
               break;
            }
        }
    }
}

void hamming_decoder_soft(r, y, message_len) 
float r[];
uint8_t y[];
int message_len;
{
   for (int idx_r = 0, idx_y = 0; idx_r < message_len; idx_r += n, idx_y += k) {
      float d = n * n;
      int codeword_id = 0;
      for (int i = 0; i < (1 << k); i++) {
         float temp = 0;
         for (int j = 0; j < n; j++) {
            uint8_t c = (codewords[i] >> j) & 0b1;
            temp += pow(map(c) - r[idx_r + j], 2);
         } if (temp <= d) { d = temp; codeword_id = i; }
      } 
      for (int i = 0; i < k; i++) {
         y[idx_y + i] = (codewords[codeword_id] >> i) & 0b1;
      }
   }
}

// modified encoder to create error syndromes
void generate_syndromes() 
{
   for (int i = 0; i < n; i++) {
      uint8_t m[n] = {}, c[2*n-k] = {};
      m[i] = 1;
      long enc_reg = 0;
      for (int i = 0; i < 2*n-k; i++) {
            uint8_t parity = enc_reg & 0b1;
            uint8_t mb = i > n ? 0 : m[i];
            uint8_t f = parity ^ mb;
            if (i >= n) f = 0;
            enc_reg >>= 1;
            enc_reg ^= (f ? taps : 0);
            c[i] = parity;
      }
      for (int j = n; j < 2*n-k; j++)
         syndromes[i] += c[j] << (j - n);
   }
}

// create all possible codewords givn {n, k} hamming code
void generate_codewords() 
{
   int c_i;
   uint8_t m[k], c[n];
   for (int i = 0; i < (1 << k); i++) {
      memset(m, 0, k);
      for (int j = 0; j < k; j++) m[j] = (i >> j) & 0b1;
      hamming_encoder(m, c, k);
      c_i = 0;
      for (int j = 0; j < n; j++) c_i |= c[j] << j; 
      codewords[i] = c_i;
   }
}
