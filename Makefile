fft00: Fft00.cpp
	g++ -lm Fft00.cpp -o fft00 -D FFT00_TEST

fftanalyze00: Fftanalyze00.cpp
	g++ -lm Fftanalyze00.cpp Fft00.cpp -o fftanalyze00 -D FFT_TEST

hilbert00: Hilbert00.cpp
	g++ -lm Hilbert00.cpp Fft00.cpp -o hilbert00 -D HILBERT_TEST

dwvd00: Dwvd00.cpp Hilbert00.cpp Fft00.cpp
	g++ -lm Dwvd00.cpp Hilbert00.cpp Fft00.cpp -o dwvd00 -D DWVD_TEST

dwvd00sp: Dwvd00.cpp Hilbert00.cpp Fft00.cpp
	g++ -lm Dwvd00.cpp Hilbert00.cpp Fft00.cpp -o dwvd00sp -D DWVDSP_TEST

gabor00: Gabor00.cpp
	g++ -lm Gabor00.cpp -o gabor00 -D GABOR_TEST

orthomra00: OrthoMra00.cpp
	g++ -lm OrthoMra00.cpp OrthoWavelet00.cpp Matrix.cpp Interpolation00.cpp Fft00.cpp -o orthomra00 -D ORTHO_MRA_TEST

orthoimra00: OrthoMra00.cpp
	g++ -lm OrthoMra00.cpp OrthoWavelet00.cpp Matrix.cpp Interpolation00.cpp Fft00.cpp -o orthoimra00 -D ORTHO_MRA_I_TEST

graphorthomra00: GraphOrthoMra00.cpp
	g++ -lm GraphOrthoMra00.cpp OrthoMra00.cpp OrthoWavelet00.cpp Matrix.cpp Interpolation00.cpp Fft00.cpp -o graphorthomra00 -D GRAPH_ORTHO_TEST

graphorthoimra00: GraphOrthoMra00.cpp
	g++ -lm GraphOrthoMra00.cpp OrthoMra00.cpp OrthoWavelet00.cpp Matrix.cpp Interpolation00.cpp Fft00.cpp -o graphorthoimra00 -D GRAPH_ORTHO_I_TEST

noisecancel00: GraphOrthoMra00.cpp
	g++ -lm GraphOrthoMra00.cpp OrthoMra00.cpp OrthoWavelet00.cpp Matrix.cpp Interpolation00.cpp Fft00.cpp -o noisecancel00 -D GRAPH_ORTHO_NOISE_TEST

sharpening00: GraphOrthoMra00.cpp
	g++ -lm GraphOrthoMra00.cpp OrthoMra00.cpp OrthoWavelet00.cpp Matrix.cpp Interpolation00.cpp Fft00.cpp -o sharpening00 -D GRAPH_ORTHO_SHARP_TEST

contour00: GraphOrthoMra00.cpp
	g++ -lm GraphOrthoMra00.cpp OrthoMra00.cpp OrthoWavelet00.cpp Matrix.cpp Interpolation00.cpp Fft00.cpp -o contour00 -D GRAPH_ORTHO_CONTOUR_TEST