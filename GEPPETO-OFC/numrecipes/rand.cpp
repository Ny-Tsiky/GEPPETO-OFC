//rand.cpp - uniform random numbers

//include this file to use rand2()
//rand2() returns a random unsigned in range [0 .. 2^32-1]
//pseudorandom sequence is defined by x[n]= x[n-24] + x[n-55]  (mod 2^32)
//see Knuth vol 2 page 27

/*
#include <fstream.h>

static unsigned randdata[55+1];        	//the last 55 numbers in the sequence, plus randsrc offset
static unsigned *randend= randdata+54; 	//pointer to end of array
static unsigned *randsrc, *randdst;    	//pointers to x[n-24], x[n-55]
static char randfile[]= "rand.dat";    	//file randdata[] is saved in
static bool write_rand= true;		//set false to reuse same random sequence next run

static struct rand_init {
    rand_init();
    ~rand_init();
} Rand_Init;

rand_init::rand_init()
//read state of random number generator from randfile
//to create initial randfile copy any old file (say this one) to "rand.dat"
{
    ifstream in(randfile,ios::binary);
    if (!in) cerr << "can't read " << randfile << '\n', exit(0);
    in.read((char*)randdata,sizeof(randdata));
    int j= randdata[55]%55;                 //mod 55 in case read bogus data
    randsrc= randdata+j;
    randdst= randdata+(j+31)%55;
}

rand_init::~rand_init()
//save state of random number generator to randfile
{
    if (!write_rand) return;		    //will use same sequence next run
    randdata[55]= randsrc-randdata;         //save randsrc offset in randdata[55]
    ofstream out(randfile,ios::binary);
    out.write((char*)randdata,sizeof(randdata)); //save state of random generator
}

unsigned rand2()
//return a random unsigned in range [0 .. 2^32-1]
{
    randdst= (randdst>randdata)? randdst-1: randend;
    randsrc= (randsrc>randdata)? randsrc-1: randend;
    return (*randdst+= *randsrc);
}
*/




static unsigned randdata[55+1];        	//the last 55 numbers in the sequence, plus randsrc offset
static unsigned *randend= randdata+54; 	//pointer to end of array
static unsigned *randsrc, *randdst;    	//pointers to x[n-24], x[n-55]
static char randfile[]= "rand.dat";    	//file randdata[] is saved in
static bool write_rand= true;		//set false to reuse same random sequence next run
FILE *randfi;

static struct rand_init {
  rand_init();
  ~rand_init();
} Rand_Init;

rand_init::rand_init()
//read state of random number generator from randfile
//to create initial randfile copy any old file (say this one) to "rand.dat"
{
  randfi = fopen(randfile,"rb");
  if (!randfi) fprintf(stderr,"can't read %s \n",randfile), exit(0);
  fread((char*)randdata,sizeof(randdata),1,randfi);
  fclose(randfi);
  int j= randdata[55]%55;                 //mod 55 in case read bogus data
  randsrc= randdata+j;
  randdst= randdata+(j+31)%55;
}

rand_init::~rand_init()
//save state of random number generator to randfile
{
  if (!write_rand) return;		    //will use same sequence next run
  randdata[55]= randsrc-randdata;         //save randsrc offset in randdata[55]
  randfi = fopen(randfile,"wb");
  fwrite((char*)randdata,sizeof(randdata),1,randfi);
  fclose(randfi);
}

unsigned rand2()
//return a random unsigned in range [0 .. 2^32-1]
{
  randdst= (randdst>randdata)? randdst-1: randend;
  randsrc= (randsrc>randdata)? randsrc-1: randend;
  return (*randdst+= *randsrc);
}
