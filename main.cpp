#include <iostream>
#include <pthread.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sstream>



using namespace std;

void exit(string message, int exitcode);
double *vector_double(long m);
long *vector_long(long m);
long **matrix_long(long m, long n);
void *graph_thread ( void *arg );



////Global variable
   pthread_mutex_t mutex; //mutex variable for file writing
   pthread_mutex_t mutex_average; //mutex variable for averaging
   int *threadind;	//Thread index
   long N=0; //Number of nodes
   double *C_average; //Averaged clustering coefficient
   double *neighbours_average; //Averaged node degree
   double *triang_average; //Averaged number of triangles


//int main(int nprocs,int ngraphs) {
int main(int argc, char *argv[]) {

	
   int ngraphs=0; //Number of graphs to be generated (read from command line argument)
   int nprocs=0; //Number of threads to use (read from command line argument)
   int i,j;	//Generic integers
   pthread_t *pth; //thread pointers for parallel execution
   pthread_attr_t pth_attr;
   clock_t start, finish; // time variables
   ofstream myfile; //Output file handle

   
if ( argc != 4 ) // 3 arguments are the default (ngraphs, nodes, nprocs)
    {
        exit("Not appropriate amount of arguments", 1);
    }
    else 
    {
        ngraphs = atoi(argv[1]); //setting number of graphs to be generated
        N = atol(argv[2]); //setting number of nodes per graph
        nprocs = atoi(argv[3]);  //setting number of threads to use
	}

	//Allocating memory for variables, and setting the array to zero
	C_average=vector_double(N);
	neighbours_average=vector_double(N);
	triang_average=vector_double(N);

    //Pointer allocation for threads
    pth = ( pthread_t * ) malloc ( nprocs * sizeof ( pthread_t ) );


    
    //Start the clock!
    start = time(NULL);
    printf("Number of graphs: %i \n",ngraphs);
    printf("Number of points: %li \n",N);
	printf("Number of threads: %i \n",nprocs);
    
    //set the attr to be joinable thread
    pthread_attr_init(&pth_attr);
    pthread_attr_setdetachstate(&pth_attr, PTHREAD_CREATE_JOINABLE);

    //Initializing the mutex variables
    pthread_mutex_init(&mutex, NULL);
    pthread_mutex_init(&mutex_average, NULL);
    
//    //Index array
    threadind = ( int * ) malloc ( ngraphs * sizeof ( int ) );
    
    j=0;
    while(j<(ngraphs-ngraphs%nprocs)){ //Will use multithreads while sensible
    //Starting threads
         for ( i = 0; i < nprocs; i++ ){
    	    threadind[j]=j;
   	        pthread_create(&pth[i],&pth_attr,graph_thread,&threadind[j]);
			j+=1;
         }
	
	//Waiting for threads to end
 	    for ( i = 0; i < nprocs; i++ ){
		    (void) pthread_join(pth[i], NULL);
	    }
	}
	while(j<ngraphs){ //Will repeat until N graphs are ready
    //Starting thread
    	    threadind[j]=j;
   	        pthread_create(&pth[0],&pth_attr,graph_thread,&threadind[j]);
			j+=1;	
	//Waiting for threads to end
		    (void) pthread_join(pth[0], NULL);
	}
	
	//destroying mutexes
    pthread_mutex_destroy(&mutex);
    pthread_mutex_destroy(&mutex_average);

	//Writing data to hard drive
    myfile.open("Averaged_data.txt");
    for ( i = 0; i < N; i++ ){
		   myfile <<i<<"\t"<<neighbours_average[i]/(double) ngraphs <<"\t"<<triang_average[i]/(double) ngraphs <<"\t"<<C_average[i]/(double) ngraphs <<"\n";
		}
    myfile.flush();
    myfile.close();
    	
	//Freeing memory
	free(C_average);
	free(neighbours_average);
	free(triang_average);

	//Exiting
	finish = time(NULL);
    printf("Total elapsed time: %i minutes %i seconds\n",(int)((long)(finish - start)/60),(int)((long)(finish - start) % 60));
    printf("Press any key to continue!");
    getchar();
    return 0;
    
}

//---------------------------------------------------------------------    
// Graph creating routione

void *graph_thread( void *arg){
  int threadID; //ID of the thread, handed down from main
  long i,j,k;	//Generic long integers for stepping
  double x;	//Generic doubles
  long new1,new2,new3; //Indices of new connections
  long Sum;	//Sum of connections
  long **links; //Matrix containing the links' endpoints
  long *neighbours; //Number of neighbours for each point
  long *triang; // Number of triangles
  double *Cluster; //Clustering coefficients
  ofstream myfile; //Output file handle
  
//  //Assigning threadID from input
  threadID=*((int*)arg);
  
      //Initializing random number generator
    srand(time(NULL)+threadID);
;
    //Throwing away the first 10000 numbers, as it is customary
    for (i = 0; i < 10000; i++) {
        x = rand();
    }

  
  //Allocating memory for variables, and setting the array to zero
  links=matrix_long(3*(N-4)+4,2);
  neighbours=vector_long(N);
  triang=vector_long(N);
  Cluster=vector_double(N);

  
  //Initializing links
  links[0][0]=0;links[0][1]=1;
  links[1][0]=0;links[1][1]=2;  
  links[2][0]=0;links[2][1]=3; 
  links[3][0]=1;links[3][1]=2; 
  
  //Initializing neighbours vector
  neighbours[0]=3;
  neighbours[1]=2;
  neighbours[2]=2;
  neighbours[3]=1;
  
  //Initializing S
  Sum=4;
  
  //Intializing number of triangles
  triang[0]=1;triang[1]=1;triang[2]=1;
  
//Creating graph
   for ( i = 4; i < N; i++ ){ 
   //for each new point	
	   //resetting connection points to forbidden value (self)
	   new1=i;new2=i;new3=i;
	   
	   //Assigning connections
	   while(new3==i){ //repeat until the last connection is also set	   
	   //getting new random number, renormalized to S
	   x=(double)rand()/((double)RAND_MAX + 1)*(double)Sum*2;
	   //Cycle to determine where should the new connection be
	   j=-1;
   	   while(x>0){
  		   j++;
	   	   x-=(double)neighbours[j];
	   	  }
	   	//See if j is a good index    
	   		if(new1==i){ //first index is not yet chosen
   			       new1=j;
			   }
			   else{
			   		if(new2==i && new1!=j){ //second index is not yet chosen, and index1 is not j
							   new2=j;
							   }
   				    else{
						 if(new1!=j && new2!=j){//index1 and index2 are not j
						 			new3=j;
									 }
						 	   }		
			   }
      }	
      
     //Sorting new conections as new1<new2<new3
     if(new1>new2){k=new1;new1=new2;new2=k;}
     if(new2>new3){k=new2;new2=new3;new3=k;}
     if(new1>new2){k=new1;new1=new2;new2=k;}
	 	 			   	   
     //New connections have been selected. Time to update everything.
	//Updating links array (indices are in increasing order, so that links[n][0]<links[n][1] )
	links[Sum][0]=new1;links[Sum][1]=i;
	links[Sum+1][0]=new2;links[Sum+1][1]=i;
	links[Sum+2][0]=new3;links[Sum+2][1]=i;
	//Updating neighbours vector
	neighbours[i]=3;
	neighbours[new1]+=1;neighbours[new2]+=1;neighbours[new3]+=1;
	
	//Updating number of triangles
	for ( k = 0; k < Sum; k++ ){ //for each link
	    //Connects any two of the three points
		if((links[k][0]==new1)&&(links[k][1]==new2)){
			triang[new1]+=1;
			triang[new2]+=1;
			triang[i]+=1;
		}
		if((links[k][0]==new1)&&(links[k][1]==new3)){
			triang[new1]+=1;
			triang[new3]+=1;
			triang[i]+=1;
		}
		if((links[k][0]==new2)&&(links[k][1]==new3)){
			triang[new2]+=1;
			triang[new3]+=1;
			triang[i]+=1;
		}
	}
	//Updating sum of nodes
	Sum+=3;
			
	}
     
    //Determining clustering coefficient
    for ( i = 0; i < 3; i++ ){
		  Cluster[i]=triang[i]*2/(double)(neighbours[i]*(neighbours[i]-1));
		}
		// it is possible for the 3. point to have only 1 neighbor, thus psecial care is to be taken here
	if(neighbours[3]!=1){Cluster[3]=triang[3]*2/(double)(neighbours[3]*(neighbours[3]-1));}
	else{Cluster[3]=0;}
    for ( i = 4; i < N; i++ ){
		  Cluster[i]=triang[i]*2/(double)(neighbours[i]*(neighbours[i]-1));
		}
	
		

    pthread_mutex_lock(&mutex);
    myfile.open("C_nozero_data.txt",ios::out | ios::app);
    for ( i = 0; i < N; i++ ){
		if(Cluster[i]!=0){
		   myfile <<i<<"\t"<<neighbours[i]<<"\t"<< triang[i]<<"\t"<<Cluster[i]<<"\n";}
		}
    myfile.flush();
    myfile.close();
    pthread_mutex_unlock(&mutex);

	//Adding results to global averaged data
	    pthread_mutex_lock(&mutex_average);
 		for ( i = 0; i < N; i++ ){
			C_average[i]+=Cluster[i];
			neighbours_average[i]+=(double)neighbours[i];
			triang_average[i]+=(double)triang[i];
		}
	    pthread_mutex_unlock(&mutex_average);

	//Freeing variables
       free(triang);
       free(neighbours);
       free(Cluster);
    for (i = 0; i < 3*(N-4)+4; i++){
              free(links[i]);}
	free(links);

	printf("Graph %i ready\n",threadID);
    pthread_exit((void*) threadID);		
}
    

    




//--------------------------------------------------------------------------------
//General routines, useful for everything

//Writing exit messages
void exit(string message, int exitcode = 0) {
    cout << "Application exited with error: " << endl << message << endl;
    exit(exitcode);
}

//Allocating vector of long integers, size: m
long *vector_long(long m) {
    long *res = (long *) calloc(m, sizeof (long));
    return (res);
}

//Allocating matrix of long integers, size: m x n
long **matrix_long(long m, long n) {
    long **res = (long **) calloc(m, sizeof (long *));
    long i;
    for (i = 0; i < m; i++)
        res[i] = (long *) calloc(n, sizeof (long));
    return (res);
}

//Allocating vector of doubles, size: m
double *vector_double(long m) {

    double *res = (double *) calloc(m, sizeof (double));
    return (res);
}



