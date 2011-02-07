

#ifdef __cplusplus
extern "C" {
#endif



	bool getDataReal(double *inData, unsigned int inDataRowLength, double *outData, unsigned int datasize,
			int step,
			unsigned long int index, unsigned int margin);
	bool getDataNew(AbstractMatrix *inData, double *outData, unsigned int datasize,
			int step,
			unsigned long int index, unsigned int margin);
	void getDataOld(char const *inData, unsigned int inDataRowLength, double *outData, unsigned int datasize,
				int step,
				unsigned long int index, unsigned int margin);





#ifdef __cplusplus
}
#endif
