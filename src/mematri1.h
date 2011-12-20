// 
// constructors
//
template <class DT> 
mematrix<DT>::mematrix(int nr, int nc)
{
	if (nr<=0)
	{
		//fprintf(stderr,"mematrix(): nr <= 0\n");
		error("mematrix(): nr <= 0");
	}
	if (nc<=0)
	{
		//fprintf(stderr,"mematrix(): nc <= 0\n");
		error("mematrix(): nc <= 0");
	}
	nrow = nr;
	ncol = nc;
	nelements = nr*nc;
	data = new (nothrow) DT [ncol*nrow];
	if (!data) 
	{
		//fprintf(stderr,"mematrix(nr,nc): cannot allocate memory (%d,%d)\n",nrow,ncol);
		error("mematrix(nr,nc): cannot allocate memory");
	}
//	fprintf(stderr,"mematrix(nr,nc): can allocate memory (%d,%d)\n",nrow,ncol);
}
template <class DT> 
mematrix<DT>::mematrix(const mematrix <DT> & M)
{
	ncol = M.ncol;
	nrow = M.nrow;
	nelements = M.nelements;
	data = new (nothrow) DT [M.ncol*M.nrow];
	if (!data)
	{
		//fprintf(stderr,"mematrix const(mematrix): cannot allocate memory (%d,%d)\n",M.nrow,M.ncol);
		error("mematrix const(mematrix): cannot allocate memory");
	}
//	fprintf(stderr,"mematrix const(mematrix): can allocate memory (%d,%d)\n",M.nrow,M.ncol);
	for(int i=0 ; i<M.ncol*M.nrow;i++) data[i] = M.data[i];
}	
// 
// operators
//
template <class DT> 
mematrix<DT> &mematrix<DT>::operator=(const mematrix<DT> &M)
{
	if(this != &M)
	{
		delete [] data;
		data = new (nothrow) DT [M.ncol*M.nrow];
		if (!data) 
		{
			//fprintf(stderr,"mematrix=: cannot allocate memory (%d,%d)\n",M.nrow,M.ncol);
			delete [] data;
			error("mematrix=: cannot allocate memory");
		}
		ncol = M.ncol;
		nrow = M.nrow;
		nelements = M.nelements;
		for(int i=0 ; i<M.ncol*M.nrow;i++) data[i] = M.data[i];
//		fprintf(stderr,"mematrix=: can allocate memory (%d,%d)\n",M.nrow,M.ncol);
	}
	return *this;
}
template <class DT> 
DT &mematrix<DT>::operator[](int i)
{
	if (i<0 || i>=(ncol*nrow)) 
	{
		//fprintf(stderr,"mematrix[]: %d out of bounds (0,%d)\n",i,nrow*ncol-1);
		error("mematrix[]: out of bounds");
	}
	return data[i];
}
template <class DT> 
mematrix<DT> mematrix<DT>::operator+(DT toadd)
{
	mematrix<DT> temp(nrow,ncol);
	for (int i=0;i<nelements;i++) temp.data[i] = data[i] + toadd;
	return temp;
}
template <class DT> 
mematrix<DT> mematrix<DT>::operator+(mematrix<DT> &M)
{
	if (ncol != M.ncol || nrow != M.nrow)
	{
		//fprintf(stderr,"mematrix+: matrices not equal in size (%d,%d) and (%d,%d)",nrow,ncol,M.nrow,M.ncol);
		error("mematrix+: matrices not equal in size");
	}
	mematrix<DT> temp(nrow,ncol);
	for (int i=0;i<nelements;i++) temp.data[i] = data[i] + M.data[i];
	return temp;
}
template <class DT> 
mematrix<DT> mematrix<DT>::operator-(DT toadd)
{
	mematrix<DT> temp(nrow,ncol);
	for (int i=0;i<nelements;i++) temp.data[i] = data[i] - toadd;
	return temp;
}
template <class DT> 
mematrix<DT> mematrix<DT>::operator-(mematrix<DT> &M)
{
	if (ncol != M.ncol || nrow != M.nrow)
	{
		//fprintf(stderr,"mematrix-: matrices not equal in size (%d,%d) and (%d,%d)",nrow,ncol,M.nrow,M.ncol);
		error("mematrix-: matrices not equal in size");
	}
	mematrix<DT> temp(nrow,ncol);
	for (int i=0;i<nelements;i++) temp.data[i] = data[i] - M.data[i];
	return temp;
}
template <class DT> 
mematrix<DT> mematrix<DT>::operator*(DT toadd)
{
	mematrix<DT> temp(nrow,ncol);
	for (int i=0;i<nelements;i++) temp.data[i] = data[i] * toadd;
	return temp;
}
template <class DT> 
mematrix<DT> mematrix<DT>::operator*(mematrix<DT> &M)
{
	if (ncol != M.nrow)
	{
		//fprintf(stderr,"mematrix*: ncol != nrow (%d,%d) and (%d,%d)",nrow,ncol,M.nrow,M.ncol);
		error("mematrix*: ncol != nrow");
	}
	mematrix<DT> temp(nrow,M.ncol);
	for (int j=0;j<temp.nrow;j++)
	{
		for (int i=0;i<temp.ncol;i++)
		{
			DT sum = 0;
			for (int j1=0;j1<ncol;j1++)
				sum += data[j*ncol+j1]*M.data[j1*M.ncol+i];
			temp[j*temp.ncol+i] = sum;
		}
	}
	return temp;
}

// 
// operations
//
template <class DT> 
void mematrix<DT>::reinit(int nr, int nc)
{
	if (nelements>0) delete [] data;
	if (nr<=0)
	{
		//fprintf(stderr,"mematrix(): nr <= 0\n");
		error("mematrix(): nr <= 0");
	}
	if (nc<=0)
	{
		//fprintf(stderr,"mematrix(): nc <= 0\n");
		error("mematrix(): nc <= 0");
	}
	nrow = nr;
	ncol = nc;
	nelements = nr*nc;
	data = new (nothrow) DT [ncol*nrow];
	if (!data) 
	{
		//fprintf(stderr,"mematrix(nr,nc): cannot allocate memory (%d,%d)\n",nrow,ncol);
		error("mematrix(nr,nc): cannot allocate memory");
	}
}
template <class DT> 
DT mematrix<DT>::get(int nr, int nc)
{
	if (nc<0 || nc>ncol) 
	{
		//fprintf(stderr,"mematrix::get: column out of range: %d not in (0,%d)\n",nc,ncol);
		error("mematrix::get: column out of range");
	}
	if (nr <0 || nr>nrow) 
	{
		//printf("mematrix::get: row out of range: %d not in (0,%d)\n",nr,nrow);
		error("mematrix::get: row out of range");
	}
	DT temp = data[nr*ncol+nc];
	return temp;
}
template <class DT> 
void mematrix<DT>::put(DT value, int nr, int nc)
{
	if (nc<0 || nc>ncol) 
	{
		//fprintf(stderr,"mematrix::put: column out of range: %d not in (0,%d)\n",nc,ncol);
		error("mematrix::put: column out of range");
	}
	if (nr <0 || nr>nrow) 
	{
		//printf("mematrix::put: row out of range: %d not in (0,%d)\n",nr,nrow);
		error("mematrix::put: row out of range");
	}
	data[nr*ncol+nc] = value;
}
template <class DT> 
DT mematrix<DT>::column_mean(int nc)
{
	if (nc>= ncol || nc <0)
	{
		//fprintf(stderr,"colmM bad column\n");
		error("colmM bad column");
	}
	DT out = 0.0;
	for (int i =0;i<nrow;i++) out+= DT(data[i*ncol+nc]);
	out /= DT(nrow);
	return out;
}
template <class DT> 
void mematrix<DT>::print(void)
{
	/**
	cout << "nrow=" << nrow << "; ncol=" << ncol << "; nelements=" << nelements << "\n";
	for (int i=0;i<nrow;i++) {
		cout << "nr=" << i << ":\t";
		for (int j=0;j<ncol;j++)
			cout << data[i*ncol+j] << "\t";
		cout << "\n";
	}
	**/
	Rprintf("mematrix::print called... but not defined :(\n");
}
template <class DT> 
void mematrix<DT>::delete_column(int delcol)
{
	if (delcol > ncol || delcol < 0) 
	{
		//fprintf(stderr,"mematrix::delete_column: column out of range\n");
		error("mematrix::delete_column: column out of range");
	}
	mematrix<DT> temp = *this;
	if (nelements>0) delete [] data;
	ncol--;
	nelements = ncol*nrow;
	data = new (nothrow) DT [ncol*nrow];
	if (!data) 
	{
		//fprintf(stderr,"mematrix::delete_column: cannot allocate memory (%d,%d)\n",nrow,ncol);
		delete [] data;
		error("mematrix::delete_column: cannot allocate memory");
	}
	int newcol=0;
	for (int nr=0;nr<temp.nrow;nr++) {
		newcol = 0;
		for (int nc=0;nc<temp.ncol;nc++)
			if (nc != delcol) data[nr*ncol+(newcol++)] = temp[nr*temp.ncol+nc];
	}

}
template <class DT> 
void mematrix<DT>::delete_row(int delrow)
{
	if (delrow > nrow || delrow < 0) 
	{
		//fprintf(stderr,"mematrix::delete_row: row out of range\n");
		error("mematrix::delete_row: row out of range");
	}
	mematrix<DT> temp = *this;
	if (nelements>0) delete [] data;
	nrow--;
	nelements = ncol*nrow;
	data = new (nothrow) DT [ncol*nrow];
	if (!data) 
	{
		//fprintf(stderr,"mematrix::delete_row: cannot allocate memory (%d,%d)\n",nrow,ncol);
		delete [] data;
		error("mematrix::delete_row: cannot allocate memory ");
	}
	int newrow=0;
	for (int nc=0;nc<temp.ncol;nc++) {
		newrow = 0;
		for (int nr=0;nr<temp.nrow;nr++)
			if (nr != delrow) data[nr*ncol+(newrow++)] = temp[nr*temp.ncol+nc];
	}

}

// 
// other functions
//
template <class DT> 
mematrix<DT> transpose(mematrix <DT> &M)
{
	mematrix<DT> temp(M.ncol,M.nrow);
	for (int i=0;i<temp.nrow;i++)
		for (int j=0;j<temp.ncol;j++)
			temp.data[i*temp.ncol+j] = M.data[j*M.ncol+i];
	return temp;
}

template <class DT> 
mematrix<DT> reorder(mematrix <DT> &M, mematrix <int> order)
{
	if (M.nrow != order.nrow)
	{
		//fprintf(stderr,"reorder: M & order have differet # of rows\n");
		error("reorder: M & order have differet # of rows");
	}
	mematrix<DT> temp(M.nrow,M.ncol);
	for (int i=0;i<temp.nrow;i++)
		for (int j=0;j<temp.ncol;j++)
			temp.data[order[i]*temp.ncol+j] = M.data[i*M.ncol+j];
	return temp;
}

template <class DT> 
mematrix<DT> productMatrDiag(mematrix <DT> &M, mematrix <DT> &D)
{
	if (M.ncol != D.nrow)
	{
		//fprintf(stderr,"productMatrDiag: wrong dimenstions");
		error("productMatrDiag: wrong dimenstions");
	}
	mematrix<DT> temp(M.nrow,M.ncol);
	for (int i=0;i<temp.nrow;i++)
		for (int j=0;j<temp.ncol;j++)
			temp.data[i*temp.ncol+j] = M.data[i*M.ncol+j]*D.data[j];
//			temp.put(M.get(i,j)*D.get(j,0),i,j);
	return temp;
}
 
template <class DT> 
mematrix<double> todouble(mematrix <DT> &M)
{
	mematrix<double> temp(M.nrow,M.ncol);
	for (int i=0;i<temp.nelements;i++)
			temp.data[i] = double(M.data[i]);
	return temp;
}

// written by Mike Dinolfo 12/98
// modified Yurii Aulchenko 2008-04-22
template <class DT> 
mematrix<DT> invert(mematrix <DT> &M)
  {
	if (M.ncol != M.nrow) 
	{
		//fprintf(stderr,"invert: only square matrices possible\n");
		error("invert: only square matrices possible");
	}
	if (M.ncol == 1) 
	{
		mematrix<DT> temp(1,1);
		temp[0] = 1./M[0];
	}
	for (int i=0;i<M.ncol;i++) 
		if (M.data[i*M.ncol+i]==0) 
		{
			//fprintf(stderr,"invert: zero elements in diagonal\n");
			error("invert: zero elements in diagonal");
		}
	int actualsize = M.ncol;
	int maxsize = M.ncol;
	mematrix<DT> temp = M;
    for (int i=1; i < actualsize; i++) temp.data[i] /= temp.data[0]; // normalize row 0
    for (int i=1; i < actualsize; i++)  { 
      for (int j=i; j < actualsize; j++)  { // do a column of L
        DT sum = 0.0;
        for (int k = 0; k < i; k++)  
            sum += temp.data[j*maxsize+k] * temp.data[k*maxsize+i];
        temp.data[j*maxsize+i] -= sum;
        }
      if (i == actualsize-1) continue;
      for (int j=i+1; j < actualsize; j++)  {  // do a row of U
        DT sum = 0.0;
        for (int k = 0; k < i; k++)
            sum += temp.data[i*maxsize+k]*temp.data[k*maxsize+j];
        temp.data[i*maxsize+j] = 
           (temp.data[i*maxsize+j]-sum) / temp.data[i*maxsize+i];
        }
      }
    for ( int i = 0; i < actualsize; i++ )  // invert L
      for ( int j = i; j < actualsize; j++ )  {
        DT x = 1.0;
        if ( i != j ) {
          x = 0.0;
          for ( int k = i; k < j; k++ ) 
              x -= temp.data[j*maxsize+k]*temp.data[k*maxsize+i];
          }
        temp.data[j*maxsize+i] = x / temp.data[j*maxsize+j];
        }
    for ( int i = 0; i < actualsize; i++ )   // invert U
      for ( int j = i; j < actualsize; j++ )  {
        if ( i == j ) continue;
        DT sum = 0.0;
        for ( int k = i; k < j; k++ )
            sum += temp.data[k*maxsize+j]*( (i==k) ? 1.0 : temp.data[i*maxsize+k] );
        temp.data[i*maxsize+j] = -sum;
        }
    for ( int i = 0; i < actualsize; i++ )   // final inversion
      for ( int j = 0; j < actualsize; j++ )  {
        DT sum = 0.0;
        for ( int k = ((i>j)?i:j); k < actualsize; k++ )  
            sum += ((j==k)?1.0:temp.data[j*maxsize+k])*temp.data[k*maxsize+i];
        temp.data[j*maxsize+i] = sum;
        }
    return temp;
    }

