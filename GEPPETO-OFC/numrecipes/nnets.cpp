#include <stdlib.h>
#include <math.h>
#include "security.h"
#include "nnets.h"

float 
Sigmoid(float t)
{
	return(tanh(t));
}

float Squash(float t)
{
  return((fabs(t)<CEILING)?t:CEILING);
}

float
DerSigmoid(float t)
{
  double d=tanh(t); 
  return(1-sqr(d));
}

float DerSquash(float t)
{
  return((fabs(t)<CEILING)?1:0);
}

double 
Rand(double min, double max)   /* Pseudo-pseudorandom number generator */
{
	return( ((rand()&127)/128.0)*(max-min) +min);
}



Layer
CreateLayer(int layer_size)
{
  float* f;

  f = (float*) calloc(layer_size, sizeof(float));
  if (f!=NULL) return(f); 
  else {
    fprintf(stderr,"Error allocating memory for layer\n");
    exit(-1);
  }
}


void
FreeLayer(Layer layer)
{
  free(layer);
}


Weights
CreateWeights(int dest_size, int source_size)
{
  float** w;
  int i;

  w = (float**) calloc(source_size, sizeof(float*));
  if (!w) {
    fprintf(stderr,"Error allocating memory for weights\n");
    exit(-2);
  }
  for (i=0; i<dest_size; i++) {
    w[i] = (float*)  calloc(source_size, sizeof(float));
    if (!w[i]){
      fprintf(stderr,"Error allocating memory for weights\n");
      exit(-2);
    }
  }
  return(w);
}

Net*
CreateNet(int in_size, int hid_size, int out_size)
{
  Net* net;

  net = (Net*) malloc(sizeof(struct Network));
  net->in=CreateLayer(in_size);
  net->hid=CreateLayer(hid_size);
  net->out=CreateLayer(out_size);
  net->w=CreateWeights(hid_size,in_size);
  net->v=CreateWeights(out_size,hid_size);
  net->size.in=in_size+1;   /* 1 for bias */
  net->size.hid=hid_size;
  net->size.out=out_size;

  return(net);
}

void
InitNet(Net* net)
{
  int i,j;

  for (i=0; i<net->size.hid; i++) 
    for (j=0; j<net->size.in; j++) net->w[i][j]=Rand(-1,1)/net->size.in/10;

  for (i=0; i<net->size.out; i++) 
    for (j=0; j<net->size.hid; j++) net->v[i][j]=Rand(-1,1)/net->size.in/10;
}

float
BackProp(Layer errors, Net* net)
{
  int i,j,k;
  float err=0;
  Layer backerrors=CreateLayer(net->size.hid);
  int 
    Nin=net->size.in,
    Nout=net->size.out,
    Nhid=net->size.hid;
  

  for (j=0; j<Nout; j++) {
    //    errors[j]*=DerSigmoid(net->out[j]);
    err+=sqr(errors[j]);
  }

  for (k=0; k<Nhid; k++) {
    backerrors[k]=0;
    for (j=0; j<Nout; j++) 
      backerrors[k]+=net->v[j][k]*errors[j];
    backerrors[k]*=DerSigmoid(net->hid[k]);
  }
  
  for (j=0; j<Nhid; j++) 
    for (k=0; k<Nin; k++) 
      net->w[j][k]+=ETA*backerrors[j]*net->in[k];

  for (j=0; j<Nout; j++)
    for (k=0; k<Nhid; k++) 
      net->v[j][k]+=ETA*errors[j]*net->hid[k];

  FreeLayer(backerrors);
  return(err);
}


void 
Prop(Layer inputs, Net* net)
{
  int i,j;
  int 
    Nin=net->size.in,
    Nout=net->size.out,
    Nhid=net->size.hid;

  for (j=0; j<Nin; j++) net->in[j]=inputs[j];
  net->in[Nin]=1;

  for (i=0; i<Nhid; i++) {
    net->hid[i]=0;
    for (j=0; j<Nin; j++)
      net->hid[i]+=net->w[i][j]*net->in[j];
    net->hid[i]=Sigmoid(net->hid[i]);
  }
  for (i=0; i<Nout; i++) {
    net->out[i]=0;
    for (j=0; j<Nin; j++)
      net->out[i]+=net->v[i][j]*net->hid[j];
    // net->out[i]=Sigmoid(net->out[i]);
  }
}


void 
SaveWeights(Net* net)
{
	FILE *fweights;
	
	fweights=OpenFile("/tmp/bw","w");
	fwrite(net->w, sizeof(float),net->size.in*net->size.hid,fweights);
	fwrite(net->v, sizeof(float),net->size.out*net->size.hid,fweights);
	fclose(fweights);
}


void 
LoadWeights(Net* net)
{
	FILE *fweights;
	
	fweights=OpenFile("/tmp/bw","r");
	fread(net->w, sizeof(float),net->size.in*net->size.hid,fweights);
	fread(net->v, sizeof(float),net->size.out*net->size.hid,fweights);
	fclose(fweights);
}



