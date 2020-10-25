#include <stdio.h>  
#include <stdlib.h>
#include <unistd.h>  
#include <math.h>


// global vars
u_int32_t sequence[8192];
u_int32_t pattern[8192];
int ts;
int tp;
int ns;
int np;
long long int charComp_bs = 0;
long long int charComp_kmp = 0;
long long int charComp_rk = 0;
unsigned long runtime_bs = 0;
unsigned long runtime_kmp = 0;
unsigned long runtime_rk = 0;

int rabinkarpSearch()
{
	struct timeval tvs1;
  	struct timeval tve1;
  	unsigned long reportedTime = 0;
  	unsigned long long int charComp = 0;
  	unsigned long long int hashComp = 0;

	int m = tp*np;
	int n = ts*ns;
	int q = 8191; // a prime larger than m
	int c = fmod(pow(2,m-1),q);
	int fp = 0;
	int ft = 0;
	int s;

	// preprocessing
	for( int i = 0; i <= m-1; i++)
	{
		fp = (2*fp + ((pattern[i + (i >> 4)] >> ((i & 0x0F) << 1)) & 0x03)) % q;
		ft = (2*ft + ((sequence[i + (i >> 4)] >> ((i & 0x0F) << 1)) & 0x03)) % q;
	}

	//matching
	int flag;
	for( s = 0; s <= n-m; s++)
	{
		hashComp++;
		if( fp == ft)
		{
			flag = 1;
			for( int j = 0; j < m && flag == 1; j++)
			{	
				charComp++;	
				if( ((pattern[j + (j >> 4)] >> ((j & 0x0F) << 1)) & 0x03) != ((sequence[(s+j) + ((s+j) >> 4)] >> (((s+j) & 0x0F) << 1)) & 0x03))
				{
					flag = 0;
					break;
				}
			}
			if( flag == 1)
			{
				gettimeofday(&tve1, NULL);
				reportedTime = tve1.tv_usec - tvs1.tv_usec;
				printf("\n-----RABIN-KARP SEARCH-----\nFound pattern at position %d.\nPerformed %lld character comparisons and %lld hash value comparisons.\n=>Total number of comparisons: %lld\nRuntime was %ld ms.\n", s+1, charComp, hashComp, charComp+hashComp,reportedTime);
				runtime_rk = reportedTime;
				charComp_rk = charComp+hashComp;
				return s;
			}
		}
		ft = ((ft - ((sequence[s + (s >> 4)] >> ((s & 0x0F) << 1)) & 0x03)*c)*2 + ((sequence[(s+m) + ((s+m) >> 4)] >> (((s+m) & 0x0F) << 1)) & 0x03))%q;
	}

	gettimeofday(&tve1, NULL);
	reportedTime = tve1.tv_usec - tvs1.tv_usec;
	printf("\n-----RABIN-KARP SEARCH-----\nPattern NOT FOUND.\nPerformed %lld character comparisons and %lld hash value comparisons.\n=>Total number of comparisons: %lld\nRuntime was %ld ms.\n", s+1, charComp, hashComp, charComp+hashComp,reportedTime);
	runtime_rk = reportedTime;
	charComp_rk = charComp+hashComp;
	return -1;			

}

int* failureFnc()
{
	int m = np*tp;
	u_int32_t* failFnc = calloc(sizeof(u_int32_t), m);
	failFnc[0] = 0;
	int i = 1;
	int j = 0;

	while( i < m)
	{
		if( ((pattern[i + (i >> 4)] >> ((i & 0x0F) << 1)) & 0x03) == ((pattern[j + (j >> 4)] >> ((j & 0x0F) << 1)) & 0x03))
		{
			// matched j+1 characters
			failFnc[i] = j+1;
			i++;
			j++;
		}
		else if( j > 0)
		{
			// use failure fnc to shift P
			j = failFnc[j-1];
		}
		else
		{
			// no match
			failFnc[i] = 0;
			i++;
		}
	}

	return failFnc;
}

int kmpSearch()
{
	struct timeval tvs2;
  	struct timeval tve2;
  	unsigned long reportedTime = 0;
  	unsigned long long int charComp = 0;
  	int m = np*tp;
	int i = 0;
	int j = 0;

	gettimeofday(&tvs2, NULL); 
	int* failFnc = failureFnc();

	while( i < ns*ts)
	{
		charComp++;
		
		if( ((sequence[i + (i >> 4)] >> ((i & 0x0F) << 1)) & 0x03) == ((pattern[j + (j >> 4)] >> ((j & 0x0F) << 1)) & 0x03))
		{	
			if( j == m-1)
			{
				// match
				gettimeofday(&tve2, NULL);
				reportedTime = tve2.tv_usec - tvs2.tv_usec;
				printf("\n-----KMP SEARCH-----\nFound pattern at position %d.\nPerformed %lld character comparisons.\nRuntime was %ld ms.\n", i-j+1, charComp, reportedTime);
				runtime_kmp = reportedTime;
				charComp_kmp = charComp;
				return i-j;
			}
			else
			{
				i++;
				j++;
			}
		}
		else
		{
			if( j > 0)
			{
				j = failFnc[j-1];
			}
			else
			{
				i++;
				j = 0;
			}
		}
	}
	// no match
	gettimeofday(&tve2, NULL);
	reportedTime = tve2.tv_usec - tvs2.tv_usec;
	printf("\n-----KMP SEARCH-----\nPattern NOT FOUND.\nPerformed %lld character comparisons.\nRuntime was %ld ms.\n", charComp, reportedTime);
	runtime_kmp = reportedTime;
	charComp_kmp = charComp;
	return -1;
}

int bruteForceSearch()
{
	struct timeval tvs3;
  	struct timeval tve3;
  	unsigned long reportedTime = 0;
	unsigned long long int charComp = 0;
	int index = 0;
	int flag = 0;

	gettimeofday(&tvs3, NULL); 
	for( int i = 0; i < ns*ts && flag == 0; i++)
	{
		flag = 1;
		for( int j = 0; j < np*tp; j++)
		{	
			charComp++;
			
			// NON-MATCH
			if( ((pattern[j + (j >> 4)] >> ((j & 0x0F) << 1)) & 0x03) != ((sequence[(i+j) + ((i+j) >> 4)] >> (((i+j) & 0x0F) << 1)) & 0x03))
			{
				flag = 0;
				break;
			}
		}
		// PATTERN FOUND
		if( flag == 1)
		{	
			gettimeofday(&tve3, NULL);
			reportedTime = tve3.tv_usec - tvs3.tv_usec;
			printf("\n----BRUTE FORCE SEARCH-----\nFound pattern at position %d.\nPerformed %lld character comparisons.\nRuntime was %ld ms.\n", i+1, charComp, reportedTime);
			runtime_bs = reportedTime;
			charComp_bs = charComp;
			return i+1;
		}
		
	}
	gettimeofday(&tve3, NULL);
	reportedTime = tve3.tv_usec - tvs3.tv_usec;
	printf("\n----BRUTE FORCE SEARCH-----\nPattern not found.\nPerformed %lld character comparisons.\nRuntime was %ld ms.\n", charComp, reportedTime);
	runtime_bs = reportedTime;
	charComp_bs = charComp;
	return -1;
}

int main(int argc, char *argv[]) 
{
	int option;
	FILE* sequence_fp;
	FILE* pattern_fp;
	char* sequence_fn;
	char* pattern_fn;

	while(( option = getopt( argc, argv, "i:p:")) != -1)
	{
		switch(option)
		{
			case 'i':
				//printf("Given Option: %c\n", option);
				//printf("filename: %s\n", optarg);
				sequence_fn = optarg;
				sequence_fp = fopen(optarg,"r");
            	break;
			case 'p':
				//printf("Given Option: %c\n", option);
				//printf("filename: %s\n", optarg);
				pattern_fn = optarg;
				pattern_fp = fopen(optarg,"r");
            	break;
		}
	}
	
	// read SEQUENCE FILE
	char c, c1, c2;

	ts = 0;
	tp = 0;
	ns = 0;
	np = 0;

	int stuff;
	fscanf(sequence_fp, "%d\n", &stuff);

	// get t and n
	int i = 0;
	while(c != EOF)
    {   
        while((c = getc(sequence_fp)) != '\n' && c != EOF)
        {
			if (ts == 1) ns++;
        }
        ts++;
    }
    ts--;
   // printf("ts: %d\n", ts); 
   // printf("ns: %d\n", ns);

    fclose(sequence_fp);
    sequence_fp = fopen(sequence_fn, "r");
    fscanf(sequence_fp, "%s\n", &c);

    // encode sequences
    /*
	* A:=00 
	* T:=01 
	* G:=10 
	* C:=11
	*/
    i = 0;
    int j = 0;

    while((c = getc(sequence_fp)) != '\n' && c != EOF);
	for(i = 0; i < ns*ts; i++)
	{
		int bit_encoding = 0;

        switch(getc(sequence_fp))
        {
            case 'A': bit_encoding = 0; break;
            case 'T': bit_encoding = 1; break;
            case 'G': bit_encoding = 2; break;
            case 'C': bit_encoding = 3; break;
            case '\n': bit_encoding = -1; j++; i--;break;
        }
        if(bit_encoding != -1)
        {
            sequence[i + i / 16] += bit_encoding << ((i % 16) * 2);
        }

        if(j >= ts) break;
	}

    // read PATTERN FILE
    
    // get t and n
	i = 0;
	while(c1 != EOF)
    {   
        while((c1 = getc(pattern_fp)) != '\n' && c1 != EOF)
        {
			if (tp == 1) np++;
        }
        tp++;
    }
    tp--;
   // printf("tp: %d\n", tp); 
   // printf("np: %d\n", np);

    fclose(pattern_fp);
    pattern_fp = fopen(pattern_fn, "r");
    fscanf(pattern_fp, "%s\n", &c1);

    // encode sequences
    i = 0;
    j = 0;

    while((c1 = getc(pattern_fp)) != '\n' && c1 != EOF);
	for(i = 0; i < np*tp; i++)
	{
		int bit_encoding = 0;

        switch(getc(pattern_fp))
        {
            case 'A': bit_encoding = 0; break;
            case 'T': bit_encoding = 1; break;
            case 'G': bit_encoding = 2; break;
            case 'C': bit_encoding = 3; break;
            case '\n': bit_encoding = -1; j++; i--;break;
        }
        if(bit_encoding != -1)
        {
            pattern[i + i / 16] += bit_encoding << ((i % 16) * 2);
        }

        if(j >= tp) break;
	}

	// BRUTE FORCE SEARCH
	bruteForceSearch();

	// SEARCH W/ KNUTH-MORRIS-PRATT 
	kmpSearch();

	// SEARCH W/ RABIN-KARP
	rabinkarpSearch();

	if( charComp_bs <= charComp_kmp && charComp_bs <= charComp_rk) printf("\nBest algorithm was BRUTE FORCE.");
	else if( charComp_kmp <= charComp_bs && charComp_kmp <= charComp_rk) printf("\nBest algorithm was KNUTH-MORRIS-PRATT.");
	else if( charComp_rk <= charComp_bs && charComp_rk <= charComp_kmp) printf("\nBest algorithm was RABIN-KARP.");
	printf(" (in terms of character comparisons)\n");

	if( runtime_bs <= runtime_kmp && runtime_bs <= runtime_rk) printf("\nBest algorithm was BRUTE FORCE.");
	else if( runtime_kmp <= runtime_bs && runtime_kmp <= runtime_rk) printf("\nBest algorithm was KNUTH-MORRIS-PRATT.");
	else if( runtime_rk <= runtime_bs && runtime_rk <= runtime_kmp) printf("\nBest algorithm was RABIN-KARP.");
	printf(" (in terms of runtime)\n");

	return 0;
}