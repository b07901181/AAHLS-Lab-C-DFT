
typedef float DTYPE;
#define SIZE 1024 		/* SIZE OF DFT */
void dft(DTYPE real_sample[1024], DTYPE imag_sample[1024],DTYPE real_op[1024], DTYPE imag_op[1024]);

void transfer_data_in(DTYPE Rin[SIZE], DTYPE Iin[SIZE],
		DTYPE Ruout[SIZE/2], DTYPE Rdout[SIZE/2], DTYPE Iuout[SIZE/2], DTYPE Idout[SIZE/2]);
/*
void transfer_data_out(DTYPE Ruin[SIZE/2], DTYPE Rdin[SIZE/2], DTYPE Iuin[SIZE/2], DTYPE Idin[SIZE/2],
		DTYPE Rout[SIZE], DTYPE Iout[SIZE]);
*/
void fft_stage(DTYPE Ruin[SIZE/2], DTYPE Rdin[SIZE/2], DTYPE Iuin[SIZE/2], DTYPE Idin[SIZE/2],
		int stage, DTYPE Ruout[SIZE/2], DTYPE Rdout[SIZE/2], DTYPE Iuout[SIZE/2], DTYPE Idout[SIZE/2]);
/*
void bit_reversal(DTYPE Ruin[SIZE/2], DTYPE Rdin[SIZE/2], DTYPE Iuin[SIZE/2], DTYPE Idin[SIZE/2],
		DTYPE Ruout[SIZE/2], DTYPE Rdout[SIZE/2], DTYPE Iuout[SIZE/2], DTYPE Idout[SIZE/2]);
*/
void bit_reversal_out(DTYPE Ruin[SIZE/2], DTYPE Rdin[SIZE/2], DTYPE Iuin[SIZE/2], DTYPE Idin[SIZE/2],
		DTYPE Rout[SIZE], DTYPE Iout[SIZE]);


const DTYPE W_real[]={1.0, 0.99998116, 0.9999247, 0.9998306, 0.9996988, 0.9995294, 0.99932235, 0.99907774, 0.99879545, 0.99847555, 0.9981181, 0.99772304, 0.99729043, 0.9968203, 0.9963126, 0.9957674, 0.9951847, 0.9945646, 0.993907, 0.9932119, 0.99247956, 0.99170977, 0.99090266, 0.9900582, 0.9891765, 0.9882576, 0.9873014, 0.9863081, 0.98527765, 0.9842101, 0.9831055, 0.9819639, 0.98078525, 0.9795698, 0.9783174, 0.97702813, 0.9757021, 0.97433937, 0.97293997, 0.9715039, 0.97003126, 0.9685221, 0.96697646, 0.96539444, 0.96377605, 0.9621214, 0.9604305, 0.95870346, 0.95694035, 0.9551412, 0.953306, 0.951435, 0.94952816, 0.9475856, 0.9456073, 0.94359344, 0.94154406, 0.9394592, 0.937339, 0.9351835, 0.9329928, 0.93076694, 0.9285061, 0.9262102, 0.9238795, 0.92151403, 0.9191139, 0.9166791, 0.9142098, 0.91170603, 0.909168, 0.9065957, 0.9039893, 0.9013488, 0.8986745, 0.89596623, 0.8932243, 0.89044875, 0.88763964, 0.8847971, 0.8819213, 0.8790122, 0.8760701, 0.873095, 0.87008697, 0.86704624, 0.86397284, 0.86086696, 0.8577286, 0.854558, 0.8513552, 0.84812033, 0.8448536, 0.841555, 0.8382247, 0.8348629, 0.8314696, 0.82804507, 0.8245893, 0.8211025, 0.8175848, 0.8140363, 0.81045717, 0.8068476, 0.8032075, 0.79953724, 0.7958369, 0.79210657, 0.7883464, 0.78455657, 0.7807372, 0.7768885, 0.77301043, 0.76910335, 0.76516724, 0.7612024, 0.7572088, 0.7531868, 0.7491364, 0.74505776, 0.7409511, 0.7368166, 0.7326543, 0.72846437, 0.7242471, 0.72000253, 0.71573085, 0.7114322, 0.70710677, 0.70275474, 0.69837624, 0.69397146, 0.68954057, 0.6850837, 0.680601, 0.6760927, 0.671559, 0.66699994, 0.6624158, 0.6578067, 0.65317285, 0.6485144, 0.64383155, 0.63912445, 0.6343933, 0.62963825, 0.6248595, 0.6200572, 0.6152316, 0.6103828, 0.60551107, 0.60061646, 0.5956993, 0.5907597, 0.58579785, 0.58081394, 0.57580817, 0.57078075, 0.5657318, 0.56066155, 0.55557024, 0.55045795, 0.545325, 0.54017144, 0.53499764, 0.52980363, 0.52458966, 0.519356, 0.51410276, 0.50883013, 0.50353837, 0.49822766, 0.4928982, 0.48755017, 0.48218378, 0.47679922, 0.47139674, 0.4659765, 0.46053872, 0.45508358, 0.44961134, 0.44412214, 0.43861625, 0.43309382, 0.42755508, 0.42200026, 0.41642955, 0.41084316, 0.4052413, 0.3996242, 0.39399204, 0.38834503, 0.38268343, 0.37700742, 0.3713172, 0.36561298, 0.35989505, 0.35416353, 0.34841868, 0.34266073, 0.33688986, 0.3311063, 0.3253103, 0.31950203, 0.31368175, 0.30784965, 0.30200595, 0.2961509, 0.29028466, 0.28440753, 0.2785197, 0.27262136, 0.26671275, 0.2607941, 0.25486565, 0.24892761, 0.24298018, 0.2370236, 0.2310581, 0.22508392, 0.21910124, 0.21311031, 0.20711137, 0.20110464, 0.19509032, 0.18906866, 0.18303989, 0.17700422, 0.17096189, 0.16491312, 0.15885815, 0.15279719, 0.14673047, 0.14065824, 0.1345807, 0.1284981, 0.12241068, 0.11631863, 0.110222206, 0.10412163, 0.09801714, 0.091908954, 0.08579731, 0.07968244, 0.07356457, 0.06744392, 0.061320737, 0.055195246, 0.049067676, 0.04293826, 0.036807224, 0.030674804, 0.024541229, 0.01840673, 0.012271538, 0.0061358847, 6.123234e-17, -0.0061358847, -0.012271538, -0.01840673, -0.024541229, -0.030674804, -0.036807224, -0.04293826, -0.049067676, -0.055195246, -0.061320737, -0.06744392, -0.07356457, -0.07968244, -0.08579731, -0.091908954, -0.09801714, -0.10412163, -0.110222206, -0.11631863, -0.12241068, -0.1284981, -0.1345807, -0.14065824, -0.14673047, -0.15279719, -0.15885815, -0.16491312, -0.17096189, -0.17700422, -0.18303989, -0.18906866, -0.19509032, -0.20110464, -0.20711137, -0.21311031, -0.21910124, -0.22508392, -0.2310581, -0.2370236, -0.24298018, -0.24892761, -0.25486565, -0.2607941, -0.26671275, -0.27262136, -0.2785197, -0.28440753, -0.29028466, -0.2961509, -0.30200595, -0.30784965, -0.31368175, -0.31950203, -0.3253103, -0.3311063, -0.33688986, -0.34266073, -0.34841868, -0.35416353, -0.35989505, -0.36561298, -0.3713172, -0.37700742, -0.38268343, -0.38834503, -0.39399204, -0.3996242, -0.4052413, -0.41084316, -0.41642955, -0.42200026, -0.42755508, -0.43309382, -0.43861625, -0.44412214, -0.44961134, -0.45508358, -0.46053872, -0.4659765, -0.47139674, -0.47679922, -0.48218378, -0.48755017, -0.4928982, -0.49822766, -0.50353837, -0.50883013, -0.51410276, -0.519356, -0.52458966, -0.52980363, -0.53499764, -0.54017144, -0.545325, -0.55045795, -0.55557024, -0.56066155, -0.5657318, -0.57078075, -0.57580817, -0.58081394, -0.58579785, -0.5907597, -0.5956993, -0.60061646, -0.60551107, -0.6103828, -0.6152316, -0.6200572, -0.6248595, -0.62963825, -0.6343933, -0.63912445, -0.64383155, -0.6485144, -0.65317285, -0.6578067, -0.6624158, -0.66699994, -0.671559, -0.6760927, -0.680601, -0.6850837, -0.68954057, -0.69397146, -0.69837624, -0.70275474, -0.70710677, -0.7114322, -0.71573085, -0.72000253, -0.7242471, -0.72846437, -0.7326543, -0.7368166, -0.7409511, -0.74505776, -0.7491364, -0.7531868, -0.7572088, -0.7612024, -0.76516724, -0.76910335, -0.77301043, -0.7768885, -0.7807372, -0.78455657, -0.7883464, -0.79210657, -0.7958369, -0.79953724, -0.8032075, -0.8068476, -0.81045717, -0.8140363, -0.8175848, -0.8211025, -0.8245893, -0.82804507, -0.8314696, -0.8348629, -0.8382247, -0.841555, -0.8448536, -0.84812033, -0.8513552, -0.854558, -0.8577286, -0.86086696, -0.86397284, -0.86704624, -0.87008697, -0.873095, -0.8760701, -0.8790122, -0.8819213, -0.8847971, -0.88763964, -0.89044875, -0.8932243, -0.89596623, -0.8986745, -0.9013488, -0.9039893, -0.9065957, -0.909168, -0.91170603, -0.9142098, -0.9166791, -0.9191139, -0.92151403, -0.9238795, -0.9262102, -0.9285061, -0.93076694, -0.9329928, -0.9351835, -0.937339, -0.9394592, -0.94154406, -0.94359344, -0.9456073, -0.9475856, -0.94952816, -0.951435, -0.953306, -0.9551412, -0.95694035, -0.95870346, -0.9604305, -0.9621214, -0.96377605, -0.96539444, -0.96697646, -0.9685221, -0.97003126, -0.9715039, -0.97293997, -0.97433937, -0.9757021, -0.97702813, -0.9783174, -0.9795698, -0.98078525, -0.9819639, -0.9831055, -0.9842101, -0.98527765, -0.9863081, -0.9873014, -0.9882576, -0.9891765, -0.9900582, -0.99090266, -0.99170977, -0.99247956, -0.9932119, -0.993907, -0.9945646, -0.9951847, -0.9957674, -0.9963126, -0.9968203, -0.99729043, -0.99772304, -0.9981181, -0.99847555, -0.99879545, -0.99907774, -0.99932235, -0.9995294, -0.9996988, -0.9998306, -0.9999247, -0.99998116};

const DTYPE W_imag[]={-0.0, -0.0061358847, -0.012271538, -0.01840673, -0.024541229, -0.030674804, -0.036807224, -0.04293826, -0.049067676, -0.055195246, -0.061320737, -0.06744392, -0.07356457, -0.07968244, -0.08579731, -0.091908954, -0.09801714, -0.10412163, -0.110222206, -0.11631863, -0.12241068, -0.1284981, -0.1345807, -0.14065824, -0.14673047, -0.15279719, -0.15885815, -0.16491312, -0.17096189, -0.17700422, -0.18303989, -0.18906866, -0.19509032, -0.20110464, -0.20711137, -0.21311031, -0.21910124, -0.22508392, -0.2310581, -0.2370236, -0.24298018, -0.24892761, -0.25486565, -0.2607941, -0.26671275, -0.27262136, -0.2785197, -0.28440753, -0.29028466, -0.2961509, -0.30200595, -0.30784965, -0.31368175, -0.31950203, -0.3253103, -0.3311063, -0.33688986, -0.34266073, -0.34841868, -0.35416353, -0.35989505, -0.36561298, -0.3713172, -0.37700742, -0.38268343, -0.38834503, -0.39399204, -0.3996242, -0.4052413, -0.41084316, -0.41642955, -0.42200026, -0.42755508, -0.43309382, -0.43861625, -0.44412214, -0.44961134, -0.45508358, -0.46053872, -0.4659765, -0.47139674, -0.47679922, -0.48218378, -0.48755017, -0.4928982, -0.49822766, -0.50353837, -0.50883013, -0.51410276, -0.519356, -0.52458966, -0.52980363, -0.53499764, -0.54017144, -0.545325, -0.55045795, -0.55557024, -0.56066155, -0.5657318, -0.57078075, -0.57580817, -0.58081394, -0.58579785, -0.5907597, -0.5956993, -0.60061646, -0.60551107, -0.6103828, -0.6152316, -0.6200572, -0.6248595, -0.62963825, -0.6343933, -0.63912445, -0.64383155, -0.6485144, -0.65317285, -0.6578067, -0.6624158, -0.66699994, -0.671559, -0.6760927, -0.680601, -0.6850837, -0.68954057, -0.69397146, -0.69837624, -0.70275474, -0.70710677, -0.7114322, -0.71573085, -0.72000253, -0.7242471, -0.72846437, -0.7326543, -0.7368166, -0.7409511, -0.74505776, -0.7491364, -0.7531868, -0.7572088, -0.7612024, -0.76516724, -0.76910335, -0.77301043, -0.7768885, -0.7807372, -0.78455657, -0.7883464, -0.79210657, -0.7958369, -0.79953724, -0.8032075, -0.8068476, -0.81045717, -0.8140363, -0.8175848, -0.8211025, -0.8245893, -0.82804507, -0.8314696, -0.8348629, -0.8382247, -0.841555, -0.8448536, -0.84812033, -0.8513552, -0.854558, -0.8577286, -0.86086696, -0.86397284, -0.86704624, -0.87008697, -0.873095, -0.8760701, -0.8790122, -0.8819213, -0.8847971, -0.88763964, -0.89044875, -0.8932243, -0.89596623, -0.8986745, -0.9013488, -0.9039893, -0.9065957, -0.909168, -0.91170603, -0.9142098, -0.9166791, -0.9191139, -0.92151403, -0.9238795, -0.9262102, -0.9285061, -0.93076694, -0.9329928, -0.9351835, -0.937339, -0.9394592, -0.94154406, -0.94359344, -0.9456073, -0.9475856, -0.94952816, -0.951435, -0.953306, -0.9551412, -0.95694035, -0.95870346, -0.9604305, -0.9621214, -0.96377605, -0.96539444, -0.96697646, -0.9685221, -0.97003126, -0.9715039, -0.97293997, -0.97433937, -0.9757021, -0.97702813, -0.9783174, -0.9795698, -0.98078525, -0.9819639, -0.9831055, -0.9842101, -0.98527765, -0.9863081, -0.9873014, -0.9882576, -0.9891765, -0.9900582, -0.99090266, -0.99170977, -0.99247956, -0.9932119, -0.993907, -0.9945646, -0.9951847, -0.9957674, -0.9963126, -0.9968203, -0.99729043, -0.99772304, -0.9981181, -0.99847555, -0.99879545, -0.99907774, -0.99932235, -0.9995294, -0.9996988, -0.9998306, -0.9999247, -0.99998116, -1.0, -0.99998116, -0.9999247, -0.9998306, -0.9996988, -0.9995294, -0.99932235, -0.99907774, -0.99879545, -0.99847555, -0.9981181, -0.99772304, -0.99729043, -0.9968203, -0.9963126, -0.9957674, -0.9951847, -0.9945646, -0.993907, -0.9932119, -0.99247956, -0.99170977, -0.99090266, -0.9900582, -0.9891765, -0.9882576, -0.9873014, -0.9863081, -0.98527765, -0.9842101, -0.9831055, -0.9819639, -0.98078525, -0.9795698, -0.9783174, -0.97702813, -0.9757021, -0.97433937, -0.97293997, -0.9715039, -0.97003126, -0.9685221, -0.96697646, -0.96539444, -0.96377605, -0.9621214, -0.9604305, -0.95870346, -0.95694035, -0.9551412, -0.953306, -0.951435, -0.94952816, -0.9475856, -0.9456073, -0.94359344, -0.94154406, -0.9394592, -0.937339, -0.9351835, -0.9329928, -0.93076694, -0.9285061, -0.9262102, -0.9238795, -0.92151403, -0.9191139, -0.9166791, -0.9142098, -0.91170603, -0.909168, -0.9065957, -0.9039893, -0.9013488, -0.8986745, -0.89596623, -0.8932243, -0.89044875, -0.88763964, -0.8847971, -0.8819213, -0.8790122, -0.8760701, -0.873095, -0.87008697, -0.86704624, -0.86397284, -0.86086696, -0.8577286, -0.854558, -0.8513552, -0.84812033, -0.8448536, -0.841555, -0.8382247, -0.8348629, -0.8314696, -0.82804507, -0.8245893, -0.8211025, -0.8175848, -0.8140363, -0.81045717, -0.8068476, -0.8032075, -0.79953724, -0.7958369, -0.79210657, -0.7883464, -0.78455657, -0.7807372, -0.7768885, -0.77301043, -0.76910335, -0.76516724, -0.7612024, -0.7572088, -0.7531868, -0.7491364, -0.74505776, -0.7409511, -0.7368166, -0.7326543, -0.72846437, -0.7242471, -0.72000253, -0.71573085, -0.7114322, -0.70710677, -0.70275474, -0.69837624, -0.69397146, -0.68954057, -0.6850837, -0.680601, -0.6760927, -0.671559, -0.66699994, -0.6624158, -0.6578067, -0.65317285, -0.6485144, -0.64383155, -0.63912445, -0.6343933, -0.62963825, -0.6248595, -0.6200572, -0.6152316, -0.6103828, -0.60551107, -0.60061646, -0.5956993, -0.5907597, -0.58579785, -0.58081394, -0.57580817, -0.57078075, -0.5657318, -0.56066155, -0.55557024, -0.55045795, -0.545325, -0.54017144, -0.53499764, -0.52980363, -0.52458966, -0.519356, -0.51410276, -0.50883013, -0.50353837, -0.49822766, -0.4928982, -0.48755017, -0.48218378, -0.47679922, -0.47139674, -0.4659765, -0.46053872, -0.45508358, -0.44961134, -0.44412214, -0.43861625, -0.43309382, -0.42755508, -0.42200026, -0.41642955, -0.41084316, -0.4052413, -0.3996242, -0.39399204, -0.38834503, -0.38268343, -0.37700742, -0.3713172, -0.36561298, -0.35989505, -0.35416353, -0.34841868, -0.34266073, -0.33688986, -0.3311063, -0.3253103, -0.31950203, -0.31368175, -0.30784965, -0.30200595, -0.2961509, -0.29028466, -0.28440753, -0.2785197, -0.27262136, -0.26671275, -0.2607941, -0.25486565, -0.24892761, -0.24298018, -0.2370236, -0.2310581, -0.22508392, -0.21910124, -0.21311031, -0.20711137, -0.20110464, -0.19509032, -0.18906866, -0.18303989, -0.17700422, -0.17096189, -0.16491312, -0.15885815, -0.15279719, -0.14673047, -0.14065824, -0.1345807, -0.1284981, -0.12241068, -0.11631863, -0.110222206, -0.10412163, -0.09801714, -0.091908954, -0.08579731, -0.07968244, -0.07356457, -0.06744392, -0.061320737, -0.055195246, -0.049067676, -0.04293826, -0.036807224, -0.030674804, -0.024541229, -0.01840673, -0.012271538, -0.0061358847};

const int BITR[]={0, 512, 256, 768, 128, 640, 384, 896, 64, 576, 320, 832, 192, 704, 448, 960, 32, 544, 288, 800, 160, 672, 416, 928, 96, 608, 352, 864, 224, 736,
480, 992, 16, 528, 272, 784, 144, 656, 400, 912, 80, 592, 336, 848, 208, 720, 464, 976, 48, 560, 304, 816, 176, 688, 432, 944, 112, 624, 368, 880, 240, 752, 496,
1008, 8, 520, 264, 776, 136, 648, 392, 904, 72, 584, 328, 840, 200, 712, 456, 968, 40, 552, 296, 808, 168, 680, 424, 936, 104, 616, 360, 872, 232, 744, 488, 1000,
24, 536, 280, 792, 152, 664, 408, 920, 88, 600, 344, 856, 216, 728, 472, 984, 56, 568, 312, 824, 184, 696, 440, 952, 120, 632, 376, 888, 248, 760, 504, 1016, 4,
516, 260, 772, 132, 644, 388, 900, 68, 580, 324, 836, 196, 708, 452, 964, 36, 548, 292, 804, 164, 676, 420, 932, 100, 612, 356, 868, 228, 740, 484, 996, 20, 532,
276, 788, 148, 660, 404, 916, 84, 596, 340, 852, 212, 724, 468, 980, 52, 564, 308, 820, 180, 692, 436, 948, 116, 628, 372, 884, 244, 756, 500, 1012, 12, 524, 268,
780, 140, 652, 396, 908, 76, 588, 332, 844, 204, 716, 460, 972, 44, 556, 300, 812, 172, 684, 428, 940, 108, 620, 364, 876, 236, 748, 492, 1004, 28, 540, 284, 796,
156, 668, 412, 924, 92, 604, 348, 860, 220, 732, 476, 988, 60, 572, 316, 828, 188, 700, 444, 956, 124, 636, 380, 892, 252, 764, 508, 1020, 2, 514, 258, 770, 130,
642, 386, 898, 66, 578, 322, 834, 194, 706, 450, 962, 34, 546, 290, 802, 162, 674, 418, 930, 98, 610, 354, 866, 226, 738, 482, 994, 18, 530, 274, 786, 146, 658,
402, 914, 82, 594, 338, 850, 210, 722, 466, 978, 50, 562, 306, 818, 178, 690, 434, 946, 114, 626, 370, 882, 242, 754, 498, 1010, 10, 522, 266, 778, 138, 650, 394,
906, 74, 586, 330, 842, 202, 714, 458, 970, 42, 554, 298, 810, 170, 682, 426, 938, 106, 618, 362, 874, 234, 746, 490, 1002, 26, 538, 282, 794, 154, 666, 410, 922,
90, 602, 346, 858, 218, 730, 474, 986, 58, 570, 314, 826, 186, 698, 442, 954, 122, 634, 378, 890, 250, 762, 506, 1018, 6, 518, 262, 774, 134, 646, 390, 902, 70,
582, 326, 838, 198, 710, 454, 966, 38, 550, 294, 806, 166, 678, 422, 934, 102, 614, 358, 870, 230, 742, 486, 998, 22, 534, 278, 790, 150, 662, 406, 918, 86, 598,
342, 854, 214, 726, 470, 982, 54, 566, 310, 822, 182, 694, 438, 950, 118, 630, 374, 886, 246, 758, 502, 1014, 14, 526, 270, 782, 142, 654, 398, 910, 78, 590, 334,
846, 206, 718, 462, 974, 46, 558, 302, 814, 174, 686, 430, 942, 110, 622, 366, 878, 238, 750, 494, 1006, 30, 542, 286, 798, 158, 670, 414, 926, 94, 606, 350, 862,
222, 734, 478, 990, 62, 574, 318, 830, 190, 702, 446, 958, 126, 638, 382, 894, 254, 766, 510, 1022, 1, 513, 257, 769, 129, 641, 385, 897, 65, 577, 321, 833, 193,
705, 449, 961, 33, 545, 289, 801, 161, 673, 417, 929, 97, 609, 353, 865, 225, 737, 481, 993, 17, 529, 273, 785, 145, 657, 401, 913, 81, 593, 337, 849, 209, 721,
465, 977, 49, 561, 305, 817, 177, 689, 433, 945, 113, 625, 369, 881, 241, 753, 497, 1009, 9, 521, 265, 777, 137, 649, 393, 905, 73, 585, 329, 841, 201, 713, 457,
969, 41, 553, 297, 809, 169, 681, 425, 937, 105, 617, 361, 873, 233, 745, 489, 1001, 25, 537, 281, 793, 153, 665, 409, 921, 89, 601, 345, 857, 217, 729, 473, 985,
57, 569, 313, 825, 185, 697, 441, 953, 121, 633, 377, 889, 249, 761, 505, 1017, 5, 517, 261, 773, 133, 645, 389, 901, 69, 581, 325, 837, 197, 709, 453, 965, 37, 549,
293, 805, 165, 677, 421, 933, 101, 613, 357, 869, 229, 741, 485, 997, 21, 533, 277, 789, 149, 661, 405, 917, 85, 597, 341, 853, 213, 725, 469, 981, 53, 565, 309, 821,
181, 693, 437, 949, 117, 629, 373, 885, 245, 757, 501, 1013, 13, 525, 269, 781, 141, 653, 397, 909, 77, 589, 333, 845, 205, 717, 461, 973, 45, 557, 301, 813, 173,
685, 429, 941, 109, 621, 365, 877, 237, 749, 493, 1005, 29, 541, 285, 797, 157, 669, 413, 925, 93, 605, 349, 861, 221, 733, 477, 989, 61, 573, 317, 829, 189, 701,
445, 957, 125, 637, 381, 893, 253, 765, 509, 1021, 3, 515, 259, 771, 131, 643, 387, 899, 67, 579, 323, 835, 195, 707, 451, 963, 35, 547, 291, 803, 163, 675, 419, 931,
99, 611, 355, 867, 227, 739, 483, 995, 19, 531, 275, 787, 147, 659, 403, 915, 83, 595, 339, 851, 211, 723, 467, 979, 51, 563, 307, 819, 179, 691, 435, 947, 115, 627,
371, 883, 243, 755, 499, 1011, 11, 523, 267, 779, 139, 651, 395, 907, 75, 587, 331, 843, 203, 715, 459, 971, 43, 555, 299, 811, 171, 683, 427, 939, 107, 619, 363, 875,
235, 747, 491, 1003, 27, 539, 283, 795, 155, 667, 411, 923, 91, 603, 347, 859, 219, 731, 475, 987, 59, 571, 315, 827, 187, 699, 443, 955, 123, 635, 379, 891, 251, 763,
507, 1019, 7, 519, 263, 775, 135, 647, 391, 903, 71, 583, 327, 839, 199, 711, 455, 967, 39, 551, 295, 807, 167, 679, 423, 935, 103, 615, 359, 871, 231, 743, 487, 999,
23, 535, 279, 791, 151, 663, 407, 919, 87, 599, 343, 855, 215, 727, 471, 983, 55, 567, 311, 823, 183, 695, 439, 951, 119, 631, 375, 887, 247, 759, 503, 1015, 15, 527,
271, 783, 143, 655, 399, 911, 79, 591, 335, 847, 207, 719, 463, 975, 47, 559, 303, 815, 175, 687, 431, 943, 111, 623, 367, 879, 239, 751, 495, 1007, 31, 543, 287, 799,
159, 671, 415, 927, 95, 607, 351, 863, 223, 735, 479, 991, 63, 575, 319, 831, 191, 703, 447, 959, 127, 639, 383, 895, 255, 767, 511, 1023};