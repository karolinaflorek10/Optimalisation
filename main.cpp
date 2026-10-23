/*********************************************
Kod stanowi uzupe�nienie materia��w do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udost�pniony na licencji CC BY-SA 3.0
Autor: dr in�. �ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G�rniczo-Hutnicza
Data ostatniej modyfikacji: 15.09.2026
*********************************************/

#include"opt_alg.h"
#include <iomanip>

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();
string to_csv_str(double val)
{
	if (isnan(val))
		return "Brak zbieznosci";
	ostringstream oss;
	oss << setprecision(20) << val; 
	string s = oss.str();
	replace(s.begin(), s.end(), '.', ',');
	return s;
}

int main()
{
	try
	{
		lab1();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	return 0;
}

void lab0()
{
	double epsilon;												// dok�adno��
	int Nmax;													// maksymalna liczba wywo�a� funkcji celu
	matrix lb, ub;												// dolne oraz g�rne ograniczenie
	solution opt;												// rozwi�zanie optymalne znalezione przez algorytm
	if (0) // Funkcja testowa
	{
		epsilon = 1e-2;
		Nmax = 10000;
		lb = matrix(2, 1, -5);
		ub = matrix(2, 1, 5);
		matrix a(2, 1);											// dok�adne rozwi�zanie optymalne
		a(0) = -1;
		a(1) = 2;
		opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);			// wywo�anie procedury optymalizacji
		cout << opt << endl << endl;							// wypisanie wyniku
		solution::clear_calls();								// wyzerowanie licznik�w
		int* n = get_size(opt.ud);
		cout << n[0] << '\t' << n[1] << endl << endl;
	}
	if (1) // Wahad�o
	{
		Nmax = 1000;
		epsilon = 1e-2;
		lb = 0, ub = 5;
		double teta_opt = 1;									// maksymalne wychylenie wahad�a
		opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);		// wywo�anie procedury optymalizacji
		cout << opt << endl << endl;							// wypisanie wyniku
		solution::clear_calls();								// wyzerowanie licznik�w

		// Zapis zmian warto�ci funckji celu do pliku
		ofstream Sout("optymalicja_lab0.csv");					// definiujemy strumie� do pliku .csv
		Sout << opt.ud;											// zapisujemy wyniki w pliku
		Sout.close();											// zamykamy strumie�

		// Zapis symulacji do pliku csv
		matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz�tkowe
			MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });	// MT zawiera moment si�y dzia�aj�cy na wahad�o oraz czas dzia�ania
		matrix* Y = solve_ode(df0, 0, 0.1, 20, Y0, NAN, MT);	// rozwi�zujemy r�wnanie r�niczkowe
		Sout.open("symulacja_lab0.csv");						// otwieramy plik .csv
		Sout << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
		Sout.close();											// zamykamy strumie�
		Y[0].~matrix();											// usuwamy z pami�ci rozwi�zanie RR
		Y[1].~matrix();
	}
}

void lab1()
{

	const double eps_fib = 1e-4;
	const double eps_lag = 1e-4;
	const double gamma_lag = 1e-6;
	const double krok_d = 5.0;
	const int limit_nmax = 10000;
	
	double wsp_alfy[3] = { 1.25, 2.1, 2.7 };
	solution optymalne;
	
	ofstream plik_tabela1("tabela1.csv");
	plik_tabela1 << "Współczynnik ekspansji;Lp.;x(0);a;b;Liczba wywołań funkcji celu;x*;y*;Liczba wywołań funkcji celu;Minimum lokalne/globalne;x*;y*;Liczba wywołań funkcji celu;Minimum lokalne/globalne" << endl;
	
	ofstream plik_tabela2("tabela2.csv");
	plik_tabela2 << "Współczynnik ekspansji;Średnia f_calls eksp.;Średnia f_calls fib;Minima globalne fib;Minima lokalne fib;Średnia f_calls lag;Minima globalne lag;Minima lokalne lag;Brak zbieżności lag" << endl;
	
	for (int idx_a = 0; idx_a < 3; ++idx_a)
	{
		int numer_porzadkowy = 1;
		
		double suma_f_eksp = 0.0;
		double suma_f_fib = 0.0;
		int fib_globalne = 0;
		int fib_lokalne = 0;
		
		double suma_f_lag = 0.0;
		int lag_globalne = 0;
		int lag_lokalne = 0;
		int lag_brak_zbieznosci = 0;
		int lag_udane_proby = 0;
		
		for (int proba = 0; proba < 100; ++proba)
		{
			double x0_start = 200.0 * m2d(rand_mat()) - 100.0;
			double* przedzial_eksp = expansion(ff1T, x0_start, krok_d, wsp_alfy[idx_a], limit_nmax);
			double wywolania_eksp = solution::f_calls;
			suma_f_eksp += wywolania_eksp;
			solution::clear_calls();
			
			optymalne = fib(ff1T, przedzial_eksp[0], przedzial_eksp[1], eps_fib);
			double fib_x_val = m2d(optymalne.x);
			double fib_y_val = m2d(optymalne.y);
			double wywolania_fib = solution::f_calls;
			suma_f_fib += wywolania_fib;
			
			string fib_typ_min = (fib_y_val < -0.5) ? "Minimum globalne" : "Minimum lokalne";
			if (fib_y_val < -0.5)
				fib_globalne++;
			else
				fib_lokalne++;
			solution::clear_calls();
			
			string lag_x_str = "Brak zbieznosci";
			string lag_y_str = "Brak zbieznosci";
			double wywolania_lag = 0;
			string lag_typ_min = "Brak zbieżności";
			
			try
			{
				optymalne = lag(ff1T, przedzial_eksp[0], przedzial_eksp[1], eps_lag, gamma_lag, limit_nmax);
				double lx = m2d(optymalne.x);
				double ly = m2d(optymalne.y);
				wywolania_lag = solution::f_calls;
				
				if (!isnan(lx) && !isnan(ly))
				{
					lag_x_str = to_csv_str(lx);
					lag_y_str = to_csv_str(ly);
					lag_udane_proby++;
					suma_f_lag += wywolania_lag;
					
					if (ly < -0.5)
					{
						lag_typ_min = "Minimum globalne";
						lag_globalne++;
					}
					else
					{
						lag_typ_min = "Minimum lokalne";
						lag_lokalne++;
					}
				}
				else
				{
					lag_brak_zbieznosci++;
				}
			}
			catch (string&) 
			{
				wywolania_lag = solution::f_calls;
				lag_brak_zbieznosci++;
			}
			catch (const char*&) 
			{
				wywolania_lag = solution::f_calls;
				lag_brak_zbieznosci++;
			}
			catch (...) 
			{
				wywolania_lag = solution::f_calls;
				lag_brak_zbieznosci++;
			}
			solution::clear_calls();

			plik_tabela1 << to_csv_str(wsp_alfy[idx_a]) << ";" 
			             << numer_porzadkowy++ << ";" 
			             << to_csv_str(x0_start) << ";" 
			             << to_csv_str(przedzial_eksp[0]) << ";" 
			             << to_csv_str(przedzial_eksp[1]) << ";" 
			             << wywolania_eksp << ";" 
			             << to_csv_str(fib_x_val) << ";" 
			             << to_csv_str(fib_y_val) << ";" 
			             << wywolania_fib << ";" 
			             << fib_typ_min << ";" 
			             << lag_x_str << ";" 
			             << lag_y_str << ";" 
			             << wywolania_lag << ";" 
			             << lag_typ_min << endl;
			
			delete[] przedzial_eksp;
		}
		
		double sr_lag_calls = (lag_udane_proby > 0) ? (suma_f_lag / lag_udane_proby) : 0.0;
		plik_tabela2 << to_csv_str(wsp_alfy[idx_a]) << ";"
		             << to_csv_str(suma_f_eksp / 100.0) << ";"
		             << to_csv_str(suma_f_fib / 100.0) << ";"
		             << fib_globalne << ";"
		             << fib_lokalne << ";"
		             << to_csv_str(sr_lag_calls) << ";"
		             << lag_globalne << ";"
		             << lag_lokalne << ";"
		             << lag_brak_zbieznosci << endl;
	}
	plik_tabela1.close();
	plik_tabela2.close();

	double real_x0 = 0.005;
	double* real_przedzial = expansion(ff1R, real_x0, 0.001, 2.0, limit_nmax);
	double real_ce = solution::f_calls;
	solution::clear_calls();

	optymalne = fib(ff1R, real_przedzial[0], real_przedzial[1], 1e-5);
	matrix DA_fib_mat = optymalne.x;
	double real_fy = m2d(optymalne.y);
	double real_cf = solution::f_calls;
	solution::clear_calls();

	matrix DA_lag_mat = DA_fib_mat;
	double real_ly = real_fy;
	double real_cl = 0;
	try
	{
		optymalne = lag(ff1R, real_przedzial[0], real_przedzial[1], 1e-5, 1e-6, limit_nmax);
		DA_lag_mat = optymalne.x;
		real_ly = m2d(optymalne.y);
		real_cl = solution::f_calls;
	}
	catch (...)
	{
		real_cl = solution::f_calls;
	}
	solution::clear_calls();

	ofstream plik_tabela3("tabela3.csv");
	plik_tabela3 << "metoda (1 = Fibonacci, 2 = Lagrange);x0 [cm2];a [cm2];b [cm2];f_calls eksp.;D_A [cm2];y;f_calls" << endl;
	plik_tabela3 << "1;" << to_csv_str(real_x0 * 1e4) << ";" << to_csv_str(real_przedzial[0] * 1e4) << ";" << to_csv_str(real_przedzial[1] * 1e4) << ";" << real_ce << ";" << to_csv_str(m2d(DA_fib_mat) * 1e4) << ";" << to_csv_str(real_fy) << ";" << real_cf << endl;
	plik_tabela3 << "2;" << to_csv_str(real_x0 * 1e4) << ";" << to_csv_str(real_przedzial[0] * 1e4) << ";" << to_csv_str(real_przedzial[1] * 1e4) << ";" << real_ce << ";" << to_csv_str(m2d(DA_lag_mat) * 1e4) << ";" << to_csv_str(real_ly) << ";" << real_cl << endl;
	plik_tabela3.close();
	delete[] real_przedzial;

	matrix stan_poczatkowy(3, new double[3] { 5.0, 1.0, 20.0 });
	matrix* wynik_ode_fib = solve_ode(df1, 0.0, 1.0, 2000, stan_poczatkowy, NAN, matrix(m2d(DA_fib_mat)));
	matrix* wynik_ode_lag = solve_ode(df1, 0.0, 1.0, 2000, stan_poczatkowy, NAN, matrix(m2d(DA_lag_mat)));
	
	ofstream plik_symulacja("symulacja.csv");
	plik_symulacja << "t;V_A fib;V_A lag;V_B fib;V_B lag;T_B fib;T_B lag" << endl;
	
	int kroki = get_len(wynik_ode_fib[0]);
	for (int k = 0; k < kroki; ++k)
	{
		plik_symulacja << to_csv_str(wynik_ode_fib[0](k)) << ";"
		               << to_csv_str(wynik_ode_fib[1](k, 0)) << ";"
		               << to_csv_str(wynik_ode_lag[1](k, 0)) << ";"
		               << to_csv_str(wynik_ode_fib[1](k, 1)) << ";"
		               << to_csv_str(wynik_ode_lag[1](k, 1)) << ";"
		               << to_csv_str(wynik_ode_fib[1](k, 2)) << ";"
		               << to_csv_str(wynik_ode_lag[1](k, 2)) << endl;
	}
	plik_symulacja.close();
	
	delete[] wynik_ode_fib;
	delete[] wynik_ode_lag;

	ofstream plik_wykres("wykres.csv");
	plik_wykres << "x;y" << endl;
	for (double x_val = -100.0; x_val <= 100.0; x_val += 0.5)
	{
		solution punkt;
		punkt.x = x_val;
		punkt.fit_fun(ff1T);
		plik_wykres << to_csv_str(x_val) << ";" << to_csv_str(m2d(punkt.y)) << endl;
	}
	plik_wykres.close();
}


void lab2()
{
// Deklaracja parametrów optymalizacji:
    // s       - początkowa długość kroku
    // alphaHJ - współczynnik zmniejszania kroku w metodzie Hooke'a-Jeevesa
    // alphaR  - współczynnik zwiększania kroku w metodzie Rosenbrocka
    // beta    - współczynnik zmniejszania i odwracania kroku w metodzie Rosenbrocka
    // epsilon - dokładność zakończenia optymalizacji
    double s, alphaHJ, alphaR, beta, epsilon;
    int Nmax; // Maksymalna liczba wywołań funkcji celu

    // Obiekt przechowujący wynik optymalizacji:
    // opt.x - znaleziony punkt, opt.y - wartość funkcji celu
    solution opt;
    // x0 - wektor punktu początkowego
    // s0 - wektor początkowych kroków dla metody Rosenbrocka
    matrix x0, s0;

    if (0) // Funkcja testowa
    {
        s = 0.5; // Początkowa długość kroku
        alphaHJ = 0.5; // Hooke-Jeeves: po nieudanym poszukiwaniu krok zmniejsza się o połowę
        alphaR = 2; // Rosenbrock: po udanym ruchu krok zwiększa się dwukrotnie
        beta = 0.5; // Rosenbrock: po nieudanym ruchu krok zmniejsza się o połowę
                    // i zmienia znak, czyli kierunek
        epsilon = 1e-3; // Dokładność optymalizacji: 0.001
        Nmax = 1000; // Limit wywołań funkcji celu
		s0 = matrix(2, 1, s); //Utworzenie wektora 2x1 zawierajacego poczatkowe kroki:
		// s0 = [0.5; 0.5]
		
		//Losowanie punktu początkowego x0 = [x1; x2]
		//
}
}
void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}
