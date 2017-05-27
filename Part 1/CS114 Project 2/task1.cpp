#include <cmath>  
#include <random>
#include <iostream>
using namespace std;

const int N = 100000;
const double pi = 3.1415926;

double mag(double Ax, double Ay, double Az);
double dotprod(double Ax, double Ay, double Az, double Bx, double By, double Bz);
double eval_Ihat(double x);
double mean(double x[]);
double std_dev(double Ihat[], double I);

int main()
{
	double *x = new double[N];
	double *Ihat = new double[N];
	double I, std, s;
	random_device rd;
	mt19937 gen(rd());
	//Task 1-1
	/*uniform_real_distribution<double> dist(-2, 2);
	for (int j = 0; j < 10; j++)
	{
		for (int i = 0; i < N; i++)
		{
			x[i] = dist(gen);
			Ihat[i] = eval_Ihat(x[i]) / .25;
		}
		I = mean(Ihat);
		std = std_dev(Ihat, I);
		s = (2 * std) / sqrt(N);
		cout << j << "      ";
		cout << "I = " << I << "         ";
		cout << "s = " << s << "         ";
		cout << "[ " << I - s << ", " << I + s << " ]\n";
	}*/
	//Task 1-2
	/*uniform_real_distribution<double> dist(0, 1);
	double p, lambda = 0;
	for (int l = 0; l < 3; l++)
	{
		switch (l) {
		case 0: lambda = .1;
			break;
		case 1: lambda = 1;
			break;
		case 2: lambda = 10;
			break;
		}
		cout << "Lambda = " << lambda << endl;
		for (int j = 0; j < 10; j++)
		{
			for (int i = 0; i < N; i++)
			{
				x[i] = 1 - (log(dist(gen)) / lambda);
				Ihat[i] = eval_Ihat(x[i]);
				p = lambda * exp(-1 * lambda * (x[i] - 1));
				Ihat[i] = Ihat[i] / p;
			}
			I = mean(Ihat);
			std = std_dev(Ihat, I);
			s = (2 * std) / sqrt(N);
			cout << j << "      ";
			cout << "I = " << I << "         ";
			cout << "s = " << s << "         ";
			cout << "[ " << I - s << ", " << I + s << " ]\n";
		}
		cout << "-------------------------------------------------\n";
	}*/
	//Task 2-1
	/*uniform_real_distribution<double> dist(0, 1); 
	double p = (1 / pow(pi, 2));
	for (int j = 0; j < 10; j++)
	{
		for (int i = 0; i < N; i++)
		{
			Ihat[i] = sin((pi / 2) * dist(gen));
			//Ihat[i] += sin((pi * 2) * dist(gen));
			Ihat[i] = Ihat[i] / (p);
		}
		I = mean(Ihat);
		std = std_dev(Ihat, I);
		s = (2 * std) / sqrt(N);
		cout << j << "      ";
		cout << "I = " << I << "         ";
		cout << "s = " << s << "         ";
		cout << "[ " << I - s << ", " << I + s << " ]\n";
	}*/
	//Task 2-2
	/*uniform_real_distribution<double> dist(0, 1);
	double p = 1 / (2 * pi);
	double theta, phi;
	for (int j = 0; j < 10; j++)
	{
		for (int i = 0; i < N; i++)
		{
			theta = acos(dist(gen));
			phi = 2 * pi * dist(gen);
			Ihat[i] = pow((cos(theta) * sin(theta) * cos(phi)), 2);
			Ihat[i] /= p;
		}
		I = mean(Ihat);
		std = std_dev(Ihat, I);
		s = (2 * std) / sqrt(N);
		cout << j << "      ";
		cout << "I = " << I << "         ";
		cout << "s = " << s << "         ";
		cout << "[ " << I - s << ", " << I + s << " ]\n";
	}*/
	//Task 3
	double Wx, Wy, Wz, Rx, Ry, Rz = 0,
		   lx, ly, lz, dis, ints, theta, phi;
	const double L = 100, Cx = 1, Cy = 1, Cz = 5, r = 1;
	double a, b, c, D; // D = b^2 - 4ac
	double p = 1 / (2 * pi);
	uniform_real_distribution<double> Rdist(-.5, .5);
	uniform_real_distribution<double> Wdist(0, 1);
	for (int j = 0; j < 10; j++)
	{
		for (int i = 0; i < N; i++)
		{
			theta = acos(Wdist(gen));
			phi = 2 * pi * Wdist(gen);
			Wx = sin(theta) * cos(phi); 
			Wy = sin(theta) * sin(phi); 
			Wz = cos(theta);
			Rx = Rdist(gen); 
			Ry = Rdist(gen);
			lx = Rx - Cx; 
			ly = Ry - Cy; 
			lz = Rz - Cz;
			a = dotprod(Wx, Wy, Wz, Wx, Wy, Wz);
			b = 2 * dotprod(Wx, Wy, Wz, lx, ly, lz);
			c = dotprod(lx, ly, lz, lx, ly, lz) - pow(r, 2);
			dis = b * b - 4 * a * c;
			if (dis < 0) ints = 0;
			else ints = 1;
			Ihat[i] = L * ints * cos(theta * pi /180);
			Ihat[i] /= p;
		}
		I = mean(Ihat);
		std = std_dev(Ihat, I);
		s = (2 * std) / sqrt(N);
		cout << j << "      ";
		cout << "I = " << I << "         ";
		cout << "s = " << s << "         ";
		cout << "[ " << I - s << ", " << I + s << " ]\n";
	}
	cin.get();
	return 0;
}
double mag(double Ax, double Ay, double Az)
{
	double sum = sqrt(pow(Ax, 2) + pow(Ay, 2) + pow(Az, 2));
	return sum;
}
double dotprod(double Ax, double Ay, double Az, double Bx, double By, double Bz)
{
	double sum = Ax * Bx + Ay * By + Az * Bz;
	return sum;
}
double eval_Ihat(double x)
{
	double I = exp(-pow(x, 2) / 2);
	return I;
}

double mean(double x[])
{
	double m = 0;
	for (int i = 0; i < N; i++)
		m += x[i];
	m /= N;
	return m;
}

double std_dev(double Ihat[], double I)
{
	double std = 0;
	for (int i = 0; i < N; i++)
		std += pow((Ihat[i] - I), 2);
	std /= (N - 1);
	std = sqrt(std);
	return std;
}