/*
 *   Sieve of Eratosthenes
 *
 *   Programmed by Michael J. Quinn
 *
 *   Last modification: 7 September 2001
 */

#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MIN(a,b)  ((a)<(b)?(a):(b))

int main (int argc, char *argv[])
{
   unsigned long int    count;        /* Local prime count */
   double elapsed_time; /* Parallel execution time */
   unsigned long int    first;        /* Index of first multiple */
   int   local_first;
   unsigned long int    global_count = 0; /* Global prime count */
   unsigned long long int    high_value;   /* Highest value on this proc */
   unsigned long int    i;
   int    id;           /* Process ID number */
   unsigned long int    index;        /* Index of current prime */
   unsigned long long int    low_value;    /* Lowest value on this proc */
   char  *marked;       /* Portion of 2,...,'n' */
   char  *local_prime_marked;
   unsigned long long int    n;            /* Sieving from 2, ..., 'n' */
   int    p;            /* Number of processes */
   unsigned long int    proc0_size;   /* Size of proc 0's subarray */
   unsigned long int    prime;
   unsigned long int  local_prime;        /* Current prime */
   unsigned long int    size;         /* Elements in 'marked' */
   unsigned long int  local_prime_size;
   unsigned long long int    low_value_first;
   unsigned long int temp;

   unsigned long int first2,prime2;
     MPI_Init (&argc, &argv);

   /* Start the timer */

   MPI_Comm_rank (MPI_COMM_WORLD, &id);
   MPI_Comm_size (MPI_COMM_WORLD, &p);
   MPI_Barrier(MPI_COMM_WORLD);
   elapsed_time = -MPI_Wtime();

   if (argc != 2) {
      if (!id) printf ("Command line: %s <m>\n", argv[0]);
      MPI_Finalize();
      exit (1);
   }

   n = atoll(argv[1]);

   /* Figure out this process's share of the array, as
      well as the integers represented by the first and
      last array elements */



   /* Add you code here  */
/* Figure out this process's share of the array, as
 *              well as the integers represented by the first and
 *                              last array elements */



        /* Add you code here  */
         low_value = 2 + id * (n - 1) / p;
         low_value_first = 2 + 0 * (n - 1) / p;
         high_value = 1 + (id + 1) * (n - 1) / p;
         size = high_value - low_value + 1;

         /* Bail out if all the primes used for sieving are
          *               not all held by process 0 */

         proc0_size = (n - 1) / p;

         if ((2 + proc0_size) < (int) sqrt((double) n)) {
                  if (!id) printf("Too many processes\n");
                  MPI_Finalize();
                  exit(1);
         }

         /* Allocate this process's share of the array. */

         marked = (char *) malloc(size);
        char *marked_prime = (char *)malloc((sqrt(n))*sizeof(char));

         if (marked == NULL) {
                  printf("Cannot allocate enough memory\n");
                  MPI_Finalize();
                  exit(1);
         }

         for (i = 0; i < size; i++) marked[i] = 0;
         if (!id) index = 0;
         for(i=0; i< (sqrt(n)); i++)
                marked_prime[i] = 0;
                if(1)
                index = 0;
                prime = 3;
                prime2 = 5;
                do
                {
                    if (prime * prime > low_value)
                        first = prime * prime - low_value;
                    else
                    {
                        if (!(low_value % prime))
                            first = 0;
                        else
                            first = prime - (low_value % prime);
                    }
                    if (prime2 * prime2 > low_value)
                        first2 = prime2 * prime2 - low_value;
                    else
                    {
                        if (!(low_value % prime2))
                            first2 = 0;
                        else
                            first2 = prime2 - (low_value % prime2);
                    }

                    int j = 0;
                    register long long int B = 100000;
                    int condition = 0;
                    if (condition)
                    {

                        for (i = 0; i < size; i += B)
                        {
                            for (j = first; j < MIN(size, i + B); j += prime)
                            {
                                temp = j + low_value;
                                if ((temp) % 2 == 0)
                                {
                                }
                                else
                                {
                                    marked[(j) / 2] = 1;
                                }
                            }
                            for (j = first2; j < MIN(size, i + B); j += prime2)
                            {
                                temp = j + low_value;
                                if ((temp) % 2 == 0)
                                {
                                }
                                else
                                {
                                    marked[(j) / 2] = 1;
                                }
                            }
                        }
                    }
                    if (!condition)
                    {
                        for (i = first; i < size; i += prime)
                        {
                            temp = i + low_value;
                            if ((temp) % 2 == 0)
                            {
                            }
                            else
                            {
                                marked[(i) / 2] = 1;
                            }
                        }
                    }

                    long long int first_0 = prime * prime - low_value_first;
                    for (i = first_0; i < (sqrt(n)); i += prime)
                    {
                        temp = i + low_value_first;
                        if ((temp) % 2 == 0)
                        {
                        }
                        else
                        {
                            marked_prime[(i) / 2] = 1;
                        }
                    }
                    first_0 = prime2 * prime2 - low_value_first;
                    for (i = first_0; i < (sqrt(n)); i += prime2)
                    {
                        temp = i + low_value_first;
                        if ((temp) % 2 == 0)
                        {
                        }
                        else
                        {
                            marked_prime[(i) / 2] = 1;
                        }
                    }

                    if (1)
                    {
                        while (marked_prime[++index])
                            ;
                        prime = index * 2 + 3;
                        if (prime2 * prime2 > n)
                        {
                            prime2 = prime;
                        }
                        prime2 = index * 2 + 3;
                    }
                } while (prime * prime <= n);

         count = 0;
         for (i = 0; i < size; i++){
                temp = i+low_value;
                if((temp%2) == 0){
                }else{
                     if(!marked[(i)/2])
                                count++;
                }
         }

         MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM,
                                                 0, MPI_COMM_WORLD);


   /* Stop the timer */

   elapsed_time += MPI_Wtime();















   /* Print the results */

   if (!id) {
      printf("The total number of prime: %ld, total time: %10.6f, total node %d\n", global_count + 1, elapsed_time, p);

   }
   MPI_Finalize ();
   return 0;
}