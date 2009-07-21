void Rprint(mematrix<double> matrix)
{
	Rprintf("nrow=%d; ncol=%d; nelements=%d\n",matrix.nrow,matrix.ncol,matrix.nelements);
	for (int i=0;i<matrix.nrow;i++) {
		Rprintf("nr=%d:\t",i);
		for (int j=0;j<matrix.ncol;j++)
			Rprintf("%f\t",matrix.data[i*matrix.ncol+j]);
		Rprintf("\n");
	}
}

