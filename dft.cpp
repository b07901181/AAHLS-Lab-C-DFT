#include<math.h>
#include "dft.h"

void dft(DTYPE real_sample[1024], DTYPE imag_sample[1024],DTYPE real_op[1024],DTYPE imag_op[1024])
{
	//Write your code here
#pragma HLS INTERFACE mode=axis port=real_sample
#pragma HLS INTERFACE mode=axis port=imag_sample
#pragma HLS INTERFACE mode=axis port=real_op
#pragma HLS INTERFACE mode=axis port=imag_op

	DTYPE R1u[SIZE/2], R2u[SIZE/2], I1u[SIZE/2], I2u[SIZE/2];
	DTYPE R1d[SIZE/2], R2d[SIZE/2], I1d[SIZE/2], I2d[SIZE/2];
#pragma HLS ARRAY_PARTITION variable=R1u,R2u,I1u,I2u,R1d,R2d,I1d,I2d type=cyclic factor=8

	transfer_data_in(real_sample, imag_sample, R1u, R1d, I1u, I1d);
	fft_stage(R1u, R1d, I1u, I1d, 1, R2u, R2d, I2u, I2d);
	fft_stage(R2u, R2d, I2u, I2d, 2, R1u, R1d, I1u, I1d);
	fft_stage(R1u, R1d, I1u, I1d, 3, R2u, R2d, I2u, I2d);
	fft_stage(R2u, R2d, I2u, I2d, 4, R1u, R1d, I1u, I1d);
	fft_stage(R1u, R1d, I1u, I1d, 5, R2u, R2d, I2u, I2d);
	fft_stage(R2u, R2d, I2u, I2d, 6, R1u, R1d, I1u, I1d);
	fft_stage(R1u, R1d, I1u, I1d, 7, R2u, R2d, I2u, I2d);
	fft_stage(R2u, R2d, I2u, I2d, 8, R1u, R1d, I1u, I1d);
	fft_stage(R1u, R1d, I1u, I1d, 9, R2u, R2d, I2u, I2d);
	fft_stage(R2u, R2d, I2u, I2d, 10, R1u, R1d, I1u, I1d);
	/*
	bit_reversal(R1u, R1d, I1u, I1d, R2u, R2d, I2u, I2d);
	transfer_data_out(R2u, R2d, I2u, I2d, real_op, imag_op);
	*/
	bit_reversal_out(R1u, R1d, I1u, I1d, real_op, imag_op);
}

void transfer_data_in(DTYPE Rin[SIZE], DTYPE Iin[SIZE],
		DTYPE Ruout[SIZE/2], DTYPE Rdout[SIZE/2], DTYPE Iuout[SIZE/2], DTYPE Idout[SIZE/2]){
	for(int i=0; i<SIZE/2; i++){
#pragma HLS PIPELINE II=1
		Ruout[i] = Rin[i];
		Iuout[i] = Iin[i];
	}
	for(int i=0; i<SIZE/2; i++){
#pragma HLS PIPELINE II=1
		Rdout[i] = Rin[i+SIZE/2];
		Idout[i] = Iin[i+SIZE/2];
	}
}
/*
void transfer_data_out(DTYPE Ruin[SIZE/2], DTYPE Rdin[SIZE/2], DTYPE Iuin[SIZE/2], DTYPE Idin[SIZE/2],
		DTYPE Rout[SIZE], DTYPE Iout[SIZE]){
	for(int i=0; i<SIZE/2; i++){
#pragma HLS PIPELINE II=1
		Rout[i] = Ruin[i];
		Iout[i] = Iuin[i];
	}
	for(int i=0; i<SIZE/2; i++){
#pragma HLS PIPELINE II=1
		Rout[i+SIZE/2] = Rdin[i];
		Iout[i+SIZE/2] = Idin[i];
	}
}
*/
void bit_reversal_out(DTYPE Ruin[SIZE/2], DTYPE Rdin[SIZE/2], DTYPE Iuin[SIZE/2], DTYPE Idin[SIZE/2],
		DTYPE Rout[SIZE], DTYPE Iout[SIZE]){
	int SIZE2 = SIZE >> 1;
	for(int i=0; i<SIZE2; i++){
#pragma HLS PIPELINE II=1
		int bitr = BITR[i];
		if(bitr < SIZE2){
			Rout[i] = Ruin[bitr];
			Iout[i] = Iuin[bitr];
		}
		else{
			Rout[i] = Rdin[bitr - SIZE2];
			Iout[i] = Idin[bitr - SIZE2];
		}
	}
	for(int i=0; i<SIZE/2; i++){
#pragma HLS PIPELINE II=1
		int bitr = BITR[i + SIZE2];
		if(bitr < SIZE2){
			Rout[i+SIZE/2] = Ruin[bitr];
			Iout[i+SIZE/2] = Iuin[bitr];
		}
		else{
			Rout[i+SIZE/2] = Rdin[bitr - SIZE2];
			Iout[i+SIZE/2] = Idin[bitr - SIZE2];
		}
	}
}

void fft_stage(DTYPE Ruin[SIZE/2], DTYPE Rdin[SIZE/2], DTYPE Iuin[SIZE/2], DTYPE Idin[SIZE/2],
		int stage, DTYPE Ruout[SIZE/2], DTYPE Rdout[SIZE/2], DTYPE Iuout[SIZE/2], DTYPE Idout[SIZE/2]){
	int SIZE4 = SIZE >> 2;
	int angle = 0;
	/*
	for(int i=0; i<SIZE4; i++){
#pragma HLS PIPELINE II=1
#pragma HLS UNROLL factor=8
		angle = (i>>(stage-1))<<(stage-1);
		Ruout[i<<1] = Ruin[i] + Rdin[i];
		Iuout[i<<1] = Iuin[i] + Idin[i];
		Ruout[(i<<1)+1] = (Ruin[i] - Rdin[i])*W_real[angle]
							+ (Idin[i] - Iuin[i])*W_imag[angle];
		Iuout[(i<<1)+1] = (Iuin[i] - Idin[i])*W_real[angle]
							+ (Ruin[i] - Rdin[i])*W_imag[angle];
	}
	for(int i=0; i<SIZE4; i++){
#pragma HLS PIPELINE II=1
#pragma HLS UNROLL factor=8
		angle = ((i + SIZE4)>>(stage-1))<<(stage-1);
		Rdout[i<<1] = Ruin[i + SIZE4] + Rdin[i + SIZE4];
		Idout[i<<1] = Iuin[i + SIZE4] + Idin[i + SIZE4];
		Rdout[(i<<1)+1] = (Ruin[i + SIZE4] - Rdin[i + SIZE4])*W_real[angle]
							+ (Idin[i + SIZE4] - Iuin[i + SIZE4])*W_imag[angle];
		Idout[(i<<1)+1] = (Iuin[i + SIZE4] - Idin[i + SIZE4])*W_real[angle]
							+ (Ruin[i + SIZE4] - Rdin[i + SIZE4])*W_imag[angle];
	}*/
	// the following code perform II=1 even with no array_partition and unroll
	for(int i=0; i<SIZE/2; i++){
#pragma HLS PIPELINE II=1
#pragma HLS UNROLL factor=16
		if(i%2==0){
			Ruout[i] = Ruin[i>>1] + Rdin[i>>1];
			Iuout[i] = Iuin[i>>1] + Idin[i>>1];
		}
		else{
			angle = (i>>stage)<<(stage-1);
			Ruout[i] = (Ruin[i>>1] - Rdin[i>>1])*W_real[angle]
						+ (Idin[i>>1] - Iuin[i>>1])*W_imag[angle];
			Iuout[i] = (Iuin[i>>1] - Idin[i>>1])*W_real[angle]
						+ (Ruin[i>>1] - Rdin[i>>1])*W_imag[angle];
		}
	}
	for(int i=0; i<SIZE/2; i++){
#pragma HLS PIPELINE II=1
#pragma HLS UNROLL factor=16
		if(i%2==0){
			Rdout[i] = Ruin[(i>>1) +SIZE4] + Rdin[(i>>1) +SIZE4];
			Idout[i] = Iuin[(i>>1) +SIZE4] + Idin[(i>>1) +SIZE4];
		}
		else{
			angle = ((i+SIZE/2)>>stage)<<(stage-1);
			Rdout[i] = (Ruin[(i>>1) +SIZE4] - Rdin[(i>>1) +SIZE4])*W_real[angle]
						+ (Idin[(i>>1) +SIZE4] - Iuin[(i>>1) +SIZE4])*W_imag[angle];
			Idout[i] = (Iuin[(i>>1) +SIZE4] - Idin[(i>>1) +SIZE4])*W_real[angle]
						+ (Ruin[(i>>1) +SIZE4] - Rdin[(i>>1) +SIZE4])*W_imag[angle];
		}
	}
}
/*
void bit_reversal(DTYPE Ruin[SIZE/2], DTYPE Rdin[SIZE/2], DTYPE Iuin[SIZE/2], DTYPE Idin[SIZE/2],
		DTYPE Ruout[SIZE/2], DTYPE Rdout[SIZE/2], DTYPE Iuout[SIZE/2], DTYPE Idout[SIZE/2]){
	for(int i=0; i<SIZE/2; i++){
#pragma HLS PIPELINE II=1
		int bitr = BITR[i];
		if(bitr < SIZE/2){
			Ruout[bitr] = Ruin[i];
			Iuout[bitr] = Iuin[i];
		}
		else{
			Rdout[bitr-SIZE/2] = Ruin[i];
			Idout[bitr-SIZE/2] = Iuin[i];
		}
	}
	for(int i=0; i<SIZE/2; i++){
#pragma HLS PIPELINE II=1
			int bitr = BITR[i+SIZE/2];
			if(bitr < SIZE/2){
				Ruout[bitr] = Rdin[i];
				Iuout[bitr] = Idin[i];
			}
			else{
				Rdout[bitr-SIZE/2] = Rdin[i];
				Idout[bitr-SIZE/2] = Idin[i];
			}
		}
}
*/
