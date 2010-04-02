#ifndef SMV_BARTLETT_TEST_H
#define SMV_BARTLETT_TEST_H


#include <list>

extern "C" {

		
class my_small_vector
	{
	
	public:
		my_small_vector(double * vector_, unsigned long number_)
			{
			vector=vector_;
			number=number_;
			}
//		~my_small_vector(void)
//			{
//			delete[] vector;
//			}

	double * vector; 
	unsigned long number; //amount of cells in vector
};



double bartlett_test(std::list<my_small_vector> * samples);

double get_mean(my_small_vector vec);
double var(my_small_vector vec);


}

#endif
