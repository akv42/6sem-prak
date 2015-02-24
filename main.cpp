#define PUAS

#include "fft.h"
#include "puasson.h"

#include <iostream>

using namespace std;

int main(int argc, char **argv)
{
#ifdef FFT
	fft_test();
#endif
#ifdef PUAS
	puas();
#endif
	cout << "Nothing to do!" << endl;
	return 0;
}