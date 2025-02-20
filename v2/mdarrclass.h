//////////////////////////////////////////////////////
// header classes for multidimensional arrays ////////
//////////////////////////////////////////////////////
// R. Hillary CSIRO 2023 /////////////////////////////
//////////////////////////////////////////////////////

class arr2d
{

public:
 
  int dd1,dd2;
  double **t; 
  double **arr2(int,int);
  double& operator()(int,int);
  void del();

};

double** arr2d::arr2(int d1,int d2)
{
  dd1 = d1;
  dd2 = d2;
  t = new double*[d1];
  for(int ii=0;ii<d1;ii++) t[ii] = new double[d2];

  return t;

}

double& arr2d::operator()(int i1,int i2) 
{
  return(t[i1][i2]);
}

void arr2d::del()
{
   
  for(int ii=0;ii<dd1;ii++) delete[] t[ii];

  delete[] t;

}

class arr3d
{

public:
 
  int dd1,dd2,dd3; 
  double ***t; 
  double ***arr3(int,int,int);
  double& operator()(int,int,int);
  void del(); 

};

double*** arr3d::arr3(int d1,int d2,int d3)
{
   
  dd1 = d1;
  dd2 = d2; 
  dd3 = d3;
  t = new double**[d1];
  for(int ii=0;ii<d1;ii++) {

    t[ii] = new double*[d2];
    for(int jj=0;jj<d2;jj++) t[ii][jj] = new double[d3];

  }

  return t;

}

double& arr3d::operator()(int i1,int i2,int i3) 
{
  return(t[i1][i2][i3]);
}

void arr3d::del()
{

  for(int ii=0;ii<dd1;ii++) 
    for(int jj=0;jj<dd2;jj++) delete[] t[ii][jj];

  for(int ii=0;ii<dd1;ii++) delete[] t[ii];
 
  delete[] t;

}

class arr4d
{

public:
 
  int dd1,dd2,dd3,dd4; 
  double ****t; 
  double ****arr4(int,int,int,int);
  double& operator()(int,int,int,int);
  void del(); 

};

double**** arr4d::arr4(int d1,int d2,int d3,int d4)
{

  dd1 = d1;
  dd2 = d2; 
  dd3 = d3;
  dd4 = d4;
  t = new double***[d1];
  for(int ii=0;ii<d1;ii++) {

    t[ii] = new double**[d2];
    for(int jj=0;jj<d2;jj++) {

      t[ii][jj] = new double*[d3];
      for(int kk=0;kk<d3;kk++) t[ii][jj][kk] = new double[d4];

    }
  }

  return t;
}

double& arr4d::operator()(int i1,int i2,int i3,int i4) 
{
  return(t[i1][i2][i3][i4]);
}

void arr4d::del()
{
   
  for(int ii=0;ii<dd1;ii++) 
    for(int jj=0;jj<dd2;jj++)  
      for(int kk=0;kk<dd3;kk++) delete[] t[ii][jj][kk];  

  for(int ii=0;ii<dd1;ii++) 
    for(int jj=0;jj<dd2;jj++) delete[] t[ii][jj];

  for(int ii=0;ii<dd1;ii++) delete[] t[ii]; 

  delete[] t;

}

class arr5d
{

public:
 
  int dd1,dd2,dd3,dd4,dd5; 
  double *****t; 
  double *****arr5(int,int,int,int,int);
  double& operator()(int,int,int,int,int);
  void del(); 

};

double***** arr5d::arr5(int d1,int d2,int d3,int d4,int d5)
{

  dd1 = d1;
  dd2 = d2; 
  dd3 = d3;
  dd4 = d4;
  dd5 = d5;
  t = new double****[d1];
  for(int ii=0;ii<d1;ii++) {

    t[ii] = new double***[d2];
    for(int jj=0;jj<d2;jj++) {

      t[ii][jj] = new double**[d3];
      for(int kk=0;kk<d3;kk++) {

        t[ii][jj][kk] = new double*[d4];
        for(int ll=0;ll<d4;ll++) t[ii][jj][kk][ll] = new double[d5];
      
      }
    }
  }

  return t;
}

double& arr5d::operator()(int i1,int i2,int i3,int i4,int i5) 
{
  return(t[i1][i2][i3][i4][i5]);
}

void arr5d::del()
{
   
  for(int ii=0;ii<dd1;ii++) 
    for(int jj=0;jj<dd2;jj++)  
      for(int kk=0;kk<dd3;kk++)
        for(int ll=0;ll<dd4;ll++) delete[] t[ii][jj][kk][ll];  

  for(int ii=0;ii<dd1;ii++) 
    for(int jj=0;jj<dd2;jj++)  
      for(int kk=0;kk<dd3;kk++) delete[] t[ii][jj][kk]; 

  for(int ii=0;ii<dd1;ii++) 
    for(int jj=0;jj<dd2;jj++) delete[] t[ii][jj];

  for(int ii=0;ii<dd1;ii++) delete[] t[ii]; 

  delete[] t;

}
