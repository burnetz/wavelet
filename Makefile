fft00: Fft00.cpp
	g++ -lm Fft00.cpp -o fft00 -D FFT00_TEST

fftanalyze00: Fftanalyze00.cpp
	g++ -lm Fftanalyze00.cpp Fft00.cpp -o fftanalyze00 -D FFT_TEST

hilbert00: Hilbert00.cpp
	g++ -lm Hilbert00.cpp Fft00.cpp -o hilbert00 -D HILBERT_TEST

dwvd00: Dwvd00.cpp Hilbert00.cpp Fft00.cpp
	g++ -lm Dwvd00.cpp Hilbert00.cpp Fft00.cpp -o dwvd00 -D DWVD_TEST