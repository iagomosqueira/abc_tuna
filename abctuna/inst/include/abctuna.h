/**
 * @author      : Iago MOSQUEIRA (WMR) (iago.mosqueira@wur.nl)
 * @file        : abctuna
 * @created     : Tuesday Feb 28, 2023 07:30:35 CET
 */

#ifndef ABCTUNA_H

#define ABCTUNA_H

// header for multidimensional arrays

// function decs

double **arr2(int,int);
double ***arr3(int,int,int);
double ****arr4(int,int,int,int);

// function defs

double **arr2(int d1,int d2) 
{
 
  double **t = new double*[d1];
  for(int ii=0;ii<d1;ii++) t[ii] = new double[d2];

  return t;

}

double ***arr3(int d1,int d2,int d3)
{
  double ***t = new double**[d1];
  for(int ii=0;ii<d1;ii++) {

    t[ii] = new double*[d2];
    for(int jj=0;jj<d2;jj++) t[ii][jj] = new double[d3];

  }

  return t;
}

double ****arr4(int d1,int d2,int d3,int d4)
{

  double ****t = new double***[d1];
  for(int ii=0;ii<d1;ii++) {

    t[ii] = new double**[d2];
    for(int jj=0;jj<d2;jj++) {

      t[ii][jj] = new double*[d3];
      for(int kk=0;kk<d3;kk++) t[ii][jj][kk] = new double[d4];

    }
  }

  return t;
}

class arr2d
{

public:
 
  double **t; 
  double **arr2(int,int);
  double& operator()(int,int);

};

double** arr2d::arr2(int d1,int d2)
{
  t = new double*[d1];
  for(int ii=0;ii<d1;ii++) t[ii] = new double[d2];

  return t;

}

double& arr2d::operator()(int i1,int i2) 
{
  return(t[i1][i2]);
}

class arr3d
{

public:
 
  double ***t; 
  double ***arr3(int,int,int);
  double& operator()(int,int,int);

};

double*** arr3d::arr3(int d1,int d2,int d3)
{
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

class arr4d
{

public:
 
  double ****t; 
  double ****arr4(int,int,int,int);
  double& operator()(int,int,int,int);

};

double**** arr4d::arr4(int d1,int d2,int d3,int d4)
{

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


#endif /* end of include guard ABCTUNA_H */

