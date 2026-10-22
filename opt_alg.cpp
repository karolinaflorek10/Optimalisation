#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	// Zmienne wej�ciowe:
	// ff - wska�nik do funkcji celu
	// N - liczba zmiennych funkcji celu
	// lb, ub - dolne i g�rne ograniczenie
	// epslion - zak��dana dok�adno�� rozwi�zania
	// Nmax - maksymalna liczba wywo�a� funkcji celu
	// ud1, ud2 - user data
	try
	{
		solution Xopt;
		Xopt = rand_mat(N);									// losujemy macierz Nx1 stosuj�c rozk�ad jednostajny na przedziale [0,1]
		for (int i = 0; i < N; ++i)							// przeskalowywujemy rozwi�zanie do przedzia�u [lb, ub]
			Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
		Xopt.fit_fun(ff, ud1, ud2);							// obliczamy warto�� funkcji celu
		Xopt.ud = (Xopt.y);									// warto�� funckji celu zapami�tujemy w macierzy ud
		if (Xopt.y < epsilon)								// sprawdzamy 1. kryterium stopu
		{
			Xopt.flag = 1;
			return Xopt;
		}
		solution X;
		while (true)
		{
			X = rand_mat(N);									// losujemy macierz Nx1 stosuj�c rozk�ad jednostajny na przedziale [0,1]
			for (int i = 0; i < N; ++i)							// przeskalowywujemy rozwi�zanie do przedzia�u [lb, ub]
				X.x(i) = (ub(i) - lb(i)) * X.x(i) + lb(i);
			X.fit_fun(ff, ud1, ud2);							// obliczmy warto�� funkcji celu
			if (X.y < Xopt.y)
			{
				Xopt = X;
				if (Xopt.y < epsilon)							// sprawdzmy 1. kryterium stopu
				{
					Xopt.flag = 1;								// flaga = 1 ozancza znalezienie rozwi�zanie z zadan� dok�adno�ci�
					break;
				}
			}
			if (solution::f_calls > Nmax)						// sprawdzmy 2. kryterium stopu
			{
				Xopt.flag = 0;									// flaga = 0 ozancza przekroczenie maksymalnej liczby wywo�a� funkcji celu
				break;
			}
			Xopt.ud.add_row(Xopt.y);							// warto�� funckji dodajemy jako kolejny wiersz macierzy ud
		}
		Xopt.ud.add_row(Xopt.y);
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2] { 0, 0 };
		solution X0, X1;
		X0.x = x0;
		X0.fit_fun(ff, ud1, ud2);
		X1.x = x0 + d;
		X1.fit_fun(ff, ud1, ud2);
		
		if (m2d(X1.y) == m2d(X0.y))
		{
			p[0] = m2d(X0.x);
			p[1] = m2d(X1.x);
			return p;
		}
		if (m2d(X1.y) > m2d(X0.y))
		{
			d = -d;
			X1.x = m2d(X0.x) + d;
			X1.fit_fun(ff, ud1, ud2);
			if (m2d(X1.y) >= m2d(X0.y))
			{
				p[0] = m2d(X1.x);
				p[1] = m2d(X0.x) - d;
				return p;
			}
		}
		
		int i = 1;
		solution X_prev = X0, X_curr = X1, X_next;
		while (true)
		{
			if (solution::f_calls > Nmax)
				throw string("Przekroczono Nmax");
				
			double x_val = x0 + pow(alpha, i) * d;
			X_next.x = x_val;
			X_next.fit_fun(ff, ud1, ud2);
			
			if (m2d(X_curr.y) <= m2d(X_next.y))
			{
				if (d > 0)
				{
					p[0] = m2d(X_prev.x);
					p[1] = m2d(X_next.x);
				}
				else
				{
					p[0] = m2d(X_next.x);
					p[1] = m2d(X_prev.x);
				}
				break;
			}
			X_prev = X_curr;
			X_curr = X_next;
			i++;
		}
		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		int k = 1;
		double sqrt5 = sqrt(5.0);
		
		while (true)
		{
			double Fk = (pow((1 + sqrt5) / 2, k) - pow((1 - sqrt5) / 2, k)) / sqrt5;
			if (Fk > (b - a) / epsilon)
				break;
			k++;
		}
		
		int n = k;
		double a_i = a, b_i = b;
		
		double* F = new double[n + 1];
		F[0] = 1; F[1] = 1;
		for (int i = 2; i <= n; ++i)
			F[i] = F[i - 1] + F[i - 2];
			
		double c = b_i - (F[n - 1] / F[n]) * (b_i - a_i);
		double d = a_i + b_i - c;
		
		solution Xc, Xd;
		Xc.x = c; Xc.fit_fun(ff, ud1, ud2);
		Xd.x = d; Xd.fit_fun(ff, ud1, ud2);
		
		for (int i = 0; i < n - 3; ++i)
		{
			if (m2d(Xc.y) < m2d(Xd.y))
			{
				b_i = d;
				d = c;
				Xd = Xc;
				c = b_i - (F[n - i - 2] / F[n - i - 1]) * (b_i - a_i);
				Xc.x = c; Xc.fit_fun(ff, ud1, ud2);
			}
			else
			{
				a_i = c;
				c = d;
				Xc = Xd;
				d = a_i + b_i - c;
				Xd.x = d; Xd.fit_fun(ff, ud1, ud2);
			}
		}
		delete[] F;
		
		Xopt = Xc;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		double c = (a + b) / 2.0;
		Xopt.x = a;
		Xopt.fit_fun(ff, ud1, ud2);
		double fa = m2d(Xopt.y);
		
		Xopt.x = b;
		Xopt.fit_fun(ff, ud1, ud2);
		double fb = m2d(Xopt.y);
		
		Xopt.x = c;
		Xopt.fit_fun(ff, ud1, ud2);
		double fc = m2d(Xopt.y);
		
		double d = c, d_prev = c, fd = fc;
		int i = 0;
		
		Xopt.ud = b - a;
		
		while (true) 
		{
			double l = fa * (b * b - c * c) + fb * (c * c - a * a) + fc * (a * a - b * b);
			double m = fa * (b - c) + fb * (c - a) + fc * (a - b);
			
			if (m <= 0.0) 
				throw string("m <= 0");
				
			d_prev = d;
			d = 0.5 * l / m;
			
			if (d <= a || d >= b || d == c) 
				throw string("d poza przedzialem");
				
			Xopt.x = d;
			Xopt.fit_fun(ff, ud1, ud2);
			fd = m2d(Xopt.y);
			
			if (d < c) 
			{
				if (fd < fc) 
				{
					b = c;
					fb = fc;
					c = d;
					fc = fd;
				}
				else 
				{
					a = d;
					fa = fd;
				}
			}
			else 
			{
				if (fd < fc) 
				{
					a = c;
					fa = fc;
					c = d;
					fc = fd;
				}
				else 
				{
					b = d;
					fb = fd;
				}
			}
			
			i++;
			Xopt.ud.add_row(b - a);
			
			if (solution::f_calls > Nmax) 
				throw string("Przekroczono Nmax");
				
			if ((b - a) < epsilon || (i > 1 && fabs(d - d_prev) < gamma)) 
				break;
		}
		
		Xopt.x = d;
		Xopt.y = fd;

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
