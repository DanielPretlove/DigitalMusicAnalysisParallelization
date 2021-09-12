using System;
using System.Diagnostics;
using System.Numerics;
using System.Threading;
using System.Threading.Tasks;
using System.Windows.Media.TextFormatting;

namespace DigitalMusicAnalysis
{
    public class timefreq
    {
        public float[][] timeFreqData;
        public int wSamp;
        public Complex[] twiddles;
        public int NumOfProcessors = Environment.ProcessorCount;
        public int N;
        public Complex[] X; 
        public float[][] Y;
        public float fftMax = 0;
        public timefreq(float[] x, int windowSamp)
        {
            double pi = 3.14159265;
            Complex i = Complex.ImaginaryOne;
            this.wSamp = windowSamp;
            twiddles = new Complex[wSamp];
            Parallel.For(0, wSamp, new ParallelOptions { MaxDegreeOfParallelism = NumOfProcessors }, ii =>
            {
                double a = 2 * pi * ii / (double)wSamp;
                twiddles[ii] = Complex.Pow(Complex.Exp(-i), (float)a);
            });

            timeFreqData = new float[wSamp / 2][];

            int nearest = (int)Math.Ceiling((double)x.Length / (double)wSamp);
            nearest = nearest * wSamp;

            Complex[] compX = new Complex[nearest];
            for (int kk = 0; kk < nearest; kk++)
            {
                if (kk < x.Length)
                {
                    compX[kk] = x[kk];
                }
                else
                {
                    compX[kk] = Complex.Zero;
                }
            }


            int cols = 2 * nearest / wSamp;

            for (int jj = 0; jj < wSamp / 2; jj++)
            {
                timeFreqData[jj] = new float[cols];
            }

            timeFreqData = stft(compX, wSamp);

        }

        float[][] stft(Complex[] x, int wSamp)
        {

            N = x.Length;
            this.X = x;

            Y = new float[wSamp / 2][];

            Parallel.For(0, wSamp / 2, new ParallelOptions { MaxDegreeOfParallelism = NumOfProcessors }, ll =>
            {
                Y[ll] = new float[2 * (int)Math.Floor((double)N / (double)wSamp)];
            });

            Thread[] num_of_threads = new Thread[NumOfProcessors];

            for (int thread = 0; thread < NumOfProcessors; thread++)
            {
                num_of_threads[thread] = new Thread(stftThreads);
                num_of_threads[thread].Start(thread);
            }

            for (int thread = 0; thread < NumOfProcessors; thread++)
            {
                num_of_threads[thread].Join();
            }


            for (int ii = 0; ii < 2 * Math.Floor((double)N / (double)wSamp) - 1; ii++)
            {
                for (int kk = 0; kk < wSamp / 2; kk++)
                {
                    Y[kk][ii] /= fftMax;
                }
            }

            return Y;
        }

        Complex[] fft(Complex[] a)
        {
            Complex[] A = bitReverseCopy(a);
            double n = A.Length;
            for(double s  = 1; s < Math.Log(n); s++)
            {
                double m = Math.Pow(2, s);
                double temp = (-2 * Math.PI) / m;
                double Omega_m = Math.Exp(-2 * Math.PI * Math.Sqrt(-1) / m);

                for(double k = 0; k < n-1; k+=m)
                {
                    double Omega = 1;
                    for(int j = 0; j < m / 2 - 1; j++)
                    {
                        Complex t = Omega * A[(int)k + j + (int)m / 2];
                        Complex u = A[(int)k + j];
                        A[(int)k + j] = u + t;
                        A[(int)k + j + (int)m / 2] = u - t;
                        Omega = Omega * Omega_m;
                    }
                }
            }
            return A;
        }

        public static string Reverse(string s)
        {
            char[] charArray = s.ToCharArray();
            Array.Reverse(charArray);
            return new string(charArray);
        }

        Complex[] bitReverseCopy(Complex[] a)
        {
            Complex[] A = a;
            int n = a.Length;
            for(int k = 0; k < n; k++)
            {
                int x = 6;
                int kk = 3;
                string binaryString = Convert.ToString(x, 2).PadLeft(kk, '0');
                int reversed = Convert.ToInt32(Reverse(binaryString), 2);
            }

            return A;
        }

        public void stftThreads(object data)
        {
            int threadId = (int)data;
            int size = (2 * (int)Math.Floor(N / (double)wSamp) - 1) / NumOfProcessors;
            int StartThread = threadId * size;
            int EndThread = Math.Min (StartThread + size, (2 * (int)Math.Floor(N / (double)wSamp) - 1));
            Complex[] temp = new Complex[wSamp];
            Complex[] tempFFT = new Complex[wSamp];

            for (int ii = StartThread; ii < EndThread; ii++)
            {

                for (int jj = 0; jj < wSamp; jj++)
                {
                    temp[jj] = X[ii * (wSamp / 2) + jj];
                }

                tempFFT = fft(temp);

                for (int kk = 0; kk < wSamp / 2; kk++)
                {
                    Y[kk][ii] = (float)Complex.Abs(tempFFT[kk]);

                    if (Y[kk][ii] > fftMax)
                    {
                        fftMax = Y[kk][ii];
                    }
                }
            }
        }

    }
}
