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
		lab3();
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
    const double eps = 1e-3;
    const int Nmax = 1000;
    const double alphaHJ = 0.5;
    const double alphaR = 2.0;
    const double beta = 0.5;

    // =========================================================
    // 5.a Testowa funkcja celu
    // =========================================================
    double steps[3] = { 0.05, 0.1, 0.2 };

    // Tabela 1: nagłówki 1:1 z szablonu Excela
    ofstream plik_tab1("tabela1.csv");
    plik_tab1 << "Długość kroku;Lp.;x1(0);x2(0);"
              << "x1*;x2*;y*;Liczba wywołań funkcji celu;Minimum globalne [TAK/NIE];"
              << "x1*;x2*;y*;Liczba wywołań funkcji celu;Minimum globalne [TAK/NIE]\n";

    // Tabela 2: wartości średnie dla minimów globalnych
    ofstream plik_tab2("tabela2.csv");
    plik_tab2 << "Długość kroku;"
              << "x1*;x2*;y*;Liczba wywołań funkcji celu;Liczba minimów globalnych;"
              << "x1*;x2*;y*;Liczba wywołań funkcji celu;Liczba minimów globalnych\n";

    for (int is = 0; is < 3; ++is)
    {
        double s = steps[is];
        double s_tab[2] = { s, s };
        matrix s0(2, s_tab);

        // Sumy do wyliczenia średnich z prób zakończonych minimum globalnym
        double sum_x1_hj = 0, sum_x2_hj = 0, sum_y_hj = 0, sum_calls_hj = 0;
        int glob_hj = 0;

        double sum_x1_r = 0, sum_x2_r = 0, sum_y_r = 0, sum_calls_r = 0;
        int glob_r = 0;

        for (int i = 0; i < 100; ++i)
        {
            matrix r = rand_mat(2, 1);
            double x0_arr[2] = { 2.0 * r(0, 0) - 1.0, 2.0 * r(1, 0) - 1.0 };
            matrix x0(2, x0_arr);

            // Hooke-Jeeves
            solution::clear_calls();
            solution opt_hj = HJ(ff2T, x0, s, alphaHJ, eps, Nmax);
            double calls_hj = solution::f_calls;
            double y_hj = m2d(opt_hj.y);
            bool is_glob_hj = (y_hj < 0.05);

            if (is_glob_hj)
            {
                glob_hj++;
                sum_x1_hj += opt_hj.x(0, 0);
                sum_x2_hj += opt_hj.x(1, 0);
                sum_y_hj += y_hj;
                sum_calls_hj += calls_hj;
            }

            // Rosenbrock
            solution::clear_calls();
            solution opt_r = Rosen(ff2T, x0, s0, alphaR, beta, eps, Nmax);
            double calls_r = solution::f_calls;
            double y_r = m2d(opt_r.y);
            bool is_glob_r = (y_r < 0.05);

            if (is_glob_r)
            {
                glob_r++;
                sum_x1_r += opt_r.x(0, 0);
                sum_x2_r += opt_r.x(1, 0);
                sum_y_r += y_r;
                sum_calls_r += calls_r;
            }
            solution::clear_calls();

            plik_tab1 << to_csv_str(s) << ";"
                      << (i + 1) << ";"
                      << to_csv_str(x0(0, 0)) << ";"
                      << to_csv_str(x0(1, 0)) << ";"
                      << to_csv_str(opt_hj.x(0, 0)) << ";"
                      << to_csv_str(opt_hj.x(1, 0)) << ";"
                      << to_csv_str(y_hj) << ";"
                      << calls_hj << ";"
                      << (is_glob_hj ? "TAK" : "NIE") << ";"
                      << to_csv_str(opt_r.x(0, 0)) << ";"
                      << to_csv_str(opt_r.x(1, 0)) << ";"
                      << to_csv_str(y_r) << ";"
                      << calls_r << ";"
                      << (is_glob_r ? "TAK" : "NIE") << "\n";
        }

        // Zapis wiersza podsumowania do Tabeli 2
        plik_tab2 << to_csv_str(s) << ";"
                  << (glob_hj > 0 ? to_csv_str(sum_x1_hj / glob_hj) : "-") << ";"
                  << (glob_hj > 0 ? to_csv_str(sum_x2_hj / glob_hj) : "-") << ";"
                  << (glob_hj > 0 ? to_csv_str(sum_y_hj / glob_hj) : "-") << ";"
                  << (glob_hj > 0 ? to_csv_str(sum_calls_hj / glob_hj) : "-") << ";"
                  << glob_hj << ";"
                  << (glob_r > 0 ? to_csv_str(sum_x1_r / glob_r) : "-") << ";"
                  << (glob_r > 0 ? to_csv_str(sum_x2_r / glob_r) : "-") << ";"
                  << (glob_r > 0 ? to_csv_str(sum_y_r / glob_r) : "-") << ";"
                  << (glob_r > 0 ? to_csv_str(sum_calls_r / glob_r) : "-") << ";"
                  << glob_r << "\n";
    }

    plik_tab1.close();
    plik_tab2.close();

    // =========================================================
    // Arkusz: Wykres 
    // =========================================================
    {
        // Losujemy punkt startowy w przedziale [-1, 1] z pełną precyzją
        matrix r_plot = rand_mat(2, 1);
        double x0_plot[2] = { 2.0 * r_plot(0, 0) - 1.0, 2.0 * r_plot(1, 0) - 1.0 };
        matrix x0_p(2, x0_plot);

        double s_plot = 0.5; // Krok 0.5 jak w konspekcie i na wzorcu
        double s0_p_arr[2] = { s_plot, s_plot };
        matrix s0_p(2, s0_p_arr);

        // 1. Hooke-Jeeves: rejestracja punktu bazowego XB w KAŻDEJ iteracji
        vector<pair<double, double>> path_hj;
        solution XB(x0_p), XB_old, X;
        XB.fit_fun(ff2T);
        path_hj.push_back({ XB.x(0, 0), XB.x(1, 0) }); // iteracja 0

        double s_curr = s_plot;
        solution::clear_calls();

        while (true)
        {
            X = HJ_trial(ff2T, XB, s_curr);
            if (X.y < XB.y)
            {
                while (true)
                {
                    XB_old = XB;
                    XB = X;
                    path_hj.push_back({ XB.x(0, 0), XB.x(1, 0) });

                    double next_x[2] = { 2.0 * XB.x(0, 0) - XB_old.x(0, 0),
                                         2.0 * XB.x(1, 0) - XB_old.x(1, 0) };
                    X.x = matrix(2, next_x);
                    X.fit_fun(ff2T);
                    X = HJ_trial(ff2T, X, s_curr);

                    if (X.y >= XB.y || solution::f_calls > Nmax)
                        break;
                }
            }
            else
            {
                s_curr *= alphaHJ;
                // Zapisujemy punkt bazowy również przy zmniejszaniu kroku (jak na wzorcu)
                path_hj.push_back({ XB.x(0, 0), XB.x(1, 0) });
            }

            if (s_curr < eps || solution::f_calls > Nmax)
                break;
        }

        // 2. Rosenbrock: rejestracja punktu bazowego XB po każdym obiegu kierunków
        vector<pair<double, double>> path_r;
        XB = solution(x0_p);
        XB.fit_fun(ff2T);
        path_r.push_back({ XB.x(0, 0), XB.x(1, 0) }); // iteracja 0

        int n_dim = 2;
        matrix l(n_dim, 1), p(n_dim, 1), s_vec(s0_p), D = ident_mat(n_dim);
        solution::clear_calls();

        while (true)
        {
            for (int j = 0; j < n_dim; ++j)
            {
                X.x = XB.x + s_vec(j, 0) * D[j];
                X.fit_fun(ff2T);
                if (X.y < XB.y)
                {
                    XB = X;
                    l(j, 0) += s_vec(j, 0);
                    s_vec(j, 0) *= alphaR;
                }
                else
                {
                    p(j, 0) += 1.0;
                    s_vec(j, 0) *= -beta;
                }
            }
            path_r.push_back({ XB.x(0, 0), XB.x(1, 0) });

            bool change = true;
            for (int j = 0; j < n_dim; ++j)
                if (l(j, 0) == 0.0 || p(j, 0) == 0.0)
                {
                    change = false;
                    break;
                }

            if (change)
            {
                matrix Q(n_dim, n_dim), v(n_dim, 1);
                for (int i = 0; i < n_dim; ++i)
                    for (int j = 0; j <= i; ++j)
                        Q(i, j) = l(i, 0);
                Q = D * Q;
                v = Q[0] / norm(Q[0]);
                D.set_col(v, 0);

                for (int i = 1; i < n_dim; ++i)
                {
                    matrix temp(n_dim, 1);
                    for (int j = 0; j < i; ++j)
                        temp = temp + trans(Q[i]) * D[j] * D[j];
                    v = (Q[i] - temp) / norm(Q[i] - temp);
                    D.set_col(v, i);
                }
                s_vec = s0_p;
                l = matrix(n_dim, 1);
                p = matrix(n_dim, 1);
            }

            double max_s = abs(s_vec(0, 0));
            if (max_s < abs(s_vec(1, 0)))
                max_s = abs(s_vec(1, 0));

            if (max_s < eps || solution::f_calls > Nmax)
                break;
        }

        // Zapis do pliku arkusz_wykres.csv z 6 miejscami po przecinku
        ofstream plik_w("arkusz_wykres.csv");
        plik_w << fixed << setprecision(6);
        plik_w << "Nr iteracji;x1*(HJ);x2*(HJ);x1*(Rosen);x2*(Rosen)\n";
        size_t max_iter = max(path_hj.size(), path_r.size());
        for (size_t it = 0; it < max_iter; ++it)
        {
            plik_w << it << ";";
            if (it < path_hj.size())
            {
                string s1 = to_string(path_hj[it].first);
                string s2 = to_string(path_hj[it].second);
                // wymuszenie polskiego przecinka
                replace(s1.begin(), s1.end(), '.', ',');
                replace(s2.begin(), s2.end(), '.', ',');
                // format do 6 cyfr po przecinku
                plik_w << s1.substr(0, s1.find(',') + 7) << ";" << s2.substr(0, s2.find(',') + 7) << ";";
            }
            else
                plik_w << ";;";

            if (it < path_r.size())
            {
                string s1 = to_string(path_r[it].first);
                string s2 = to_string(path_r[it].second);
                replace(s1.begin(), s1.end(), '.', ',');
                replace(s2.begin(), s2.end(), '.', ',');
                plik_w << s1.substr(0, s1.find(',') + 7) << ";" << s2.substr(0, s2.find(',') + 7) << "\n";
            }
            else
                plik_w << ";\n";
        }
        plik_w.close();
        solution::clear_calls();
    }

    // =========================================================
    // 5.b Problem rzeczywisty (Tabela 3 i Symulacja)
    // =========================================================
    double s_real = 2.0;
    double s0_real_arr[2] = { s_real, s_real };
    matrix s0_real(2, s0_real_arr);

    matrix r_real = rand_mat(2, 1);
    double x0_real_arr[2] = { 20.0 * r_real(0, 0), 20.0 * r_real(1, 0) };
    matrix x0_real(2, x0_real_arr);

    // HJ
    solution::clear_calls();
    solution real_hj = HJ(ff2R, x0_real, s_real, alphaHJ, eps, Nmax);
    double calls_real_hj = solution::f_calls;
    solution::clear_calls();

    // Rosenbrock
    solution real_r = Rosen(ff2R, x0_real, s0_real, alphaR, beta, eps, Nmax);
    double calls_real_r = solution::f_calls;
    solution::clear_calls();

    // Tabela 3
    ofstream plik_tab3("tabela3.csv");
    plik_tab3 << "Długość kroku;k1*;k2*;Q*;Liczba wywołań funkcji celu;k1*;k2*;Q*;Liczba wywołań funkcji celu\n";
    plik_tab3 << to_csv_str(s_real) << ";"
              << to_csv_str(real_hj.x(0, 0)) << ";"
              << to_csv_str(real_hj.x(1, 0)) << ";"
              << to_csv_str(m2d(real_hj.y)) << ";"
              << calls_real_hj << ";"
              << to_csv_str(real_r.x(0, 0)) << ";"
              << to_csv_str(real_r.x(1, 0)) << ";"
              << to_csv_str(m2d(real_r.y)) << ";"
              << calls_real_r << "\n";
    plik_tab3.close();

    // Arkusz Symulacja: symulacja dla obu metod
    double y0_tab[2] = { 0.0, 0.0 };
    matrix Y0(2, y0_tab);
    double yref_tab[2] = { 3.14, 0.0 };
    matrix Y_ref(2, yref_tab);

    matrix* Y_sim_hj = solve_ode(df2, 0.0, 0.1, 100.0, Y0, Y_ref, real_hj.x);
    matrix* Y_sim_r  = solve_ode(df2, 0.0, 0.1, 100.0, Y0, Y_ref, real_r.x);

    ofstream plik_sym("symulacja.csv");
    plik_sym << "t;alpha_HJ;omega_HJ;alpha_Rosen;omega_Rosen\n";
    int* sim_sz = get_size(Y_sim_hj[0]);
    int n_steps = sim_sz[0];
    for (int t = 0; t < n_steps; ++t)
    {
        plik_sym << to_csv_str(Y_sim_hj[0](t, 0)) << ";"
                 << to_csv_str(Y_sim_hj[1](t, 0)) << ";"
                 << to_csv_str(Y_sim_hj[1](t, 1)) << ";"
                 << to_csv_str(Y_sim_r[1](t, 0)) << ";"
                 << to_csv_str(Y_sim_r[1](t, 1)) << "\n";
    }
    plik_sym.close();
    delete[] sim_sz;
    delete[] Y_sim_hj;
    delete[] Y_sim_r;

    cout << "Wygenerowano pliki\n";
}

void lab3()
{

    matrix x0; // Punkt początkowy optymalizacji

    // c_ex - początkowy współczynnik zewnętrznej funkcji kary
    // c_in - początkowy współczynnik wewnętrznej funkcji kary
    // dc_ex - współczynnik zwiększania kary zewnętrznej
    // dc_in - współczynnik zmniejszania kary wewnętrznej
    // epsilon - dokładność optymalizacji

    double c_ex, c_in, dc_ex, dc_in, epsilon;

    int Nmax;     // Maksymalna liczba wywołań funkcji celu
    solution opt; // Obiekt przechowujący rozwiązanie

    // ==========================================
    // CZĘŚĆ 1 - FUNKCJA TESTOWA
    // ==========================================

    if (0) // Funkcja testowa
    {
        // Parametr ograniczenia:
        //
        // x1^2 + x2^2 <= a^2
        //
        // W konspekcie należy przyjąć kolejno:
        // a = 4, a = 4.4934 oraz a = 5.
        //
        // Tutaj rozpatrujemy przypadek a = 4.

        matrix a(4);

        // Początkowa wartość współczynnika kary zewnętrznej.
        c_ex = 1;

        // Początkowa wartość współczynnika kary wewnętrznej.
        c_in = 10;

        // W metodzie zewnętrznej zwiększamy karę
        // po każdej optymalizacji.
        //
        // Przykład: 1 -> 2 -> 4 -> 8 -> ...

        dc_ex = 2;

        // W metodzie wewnętrznej zmniejszamy
        // współczynnik funkcji bariery.
        //
        // Przykład: 10 -> 5 -> 2.5 -> 1.25 -> ...

        dc_in = 0.5;

        epsilon = 1e-4;
        Nmax = 10000;

        // Losujemy punkt początkowy.
        //
        // rand_mat(2,1) zwraca dwie liczby losowe
        // z przedziału [0,1].
        //
        // Po przekształceniu otrzymujemy współrzędne
        // z przedziału [1,6].

        do
        {
            x0 = 5 * rand_mat(2, 1) + 1;
        }
        // Punkt musi znajdować się
        // wewnątrz koła o promieniu a.
        //
        // norm(x0) = sqrt(x1^2 + x2^2)

        while (norm(x0) > a);

        cout << x0 << endl << endl;

        // ZEWNĘTRZNA FUNKCJA KARY
        //
        // ff3Ta - funkcja celu z karą zewnętrzną.
        // a     - parametr ograniczenia.

        opt = pen(ff3Ta, x0, c_ex, dc_ex, epsilon, Nmax, a);

        cout << opt << endl;

        // Wyświetlamy odległość znalezionego punktu
        // od początku układu współrzędnych.
        //
        // Ta wartość jest potrzebna w Tabeli 1 (kolumna r*).

        cout << norm(opt.x) << endl << endl;

        // Zerujemy licznik wywołań funkcji celu
        // przed uruchomieniem kolejnej metody.

        solution::clear_calls();

        // WEWNĘTRZNA FUNKCJA KARY
        //
        // Używamy tego samego punktu początkowego x0.

        opt = pen(ff3Tb, x0, c_in, dc_in, epsilon, Nmax, a);

        cout << opt << endl;
        cout << norm(opt.x) << endl << endl;

        solution::clear_calls();
    }

    // ==========================================
    // CZĘŚĆ 2 - PROBLEM RZECZYWISTY
    // ==========================================

    if (0) // Rzut piłka
    {
        // Dla problemu rzeczywistego
        // stosujemy wyłącznie zewnętrzną funkcję kary.

        c_ex = 1;
        dc_ex = 2;

        epsilon = 1e-4;
        Nmax = 10000;

        // Losowanie dwóch parametrów:
        //
        // x0(0) - początkowa prędkość pozioma v0x
        // x0(1) - prędkość kątowa omega
        //
        // Oba parametry należą do przedziału [-10,10].
        x0 = matrix(2, 1);

        x0(0) = 20 * m2d(rand_mat()) - 10;
        x0(1) = 20 * m2d(rand_mat()) - 10;

        cout << x0 << endl << endl;

        // Optymalizacja parametrów lotu piłki.
        //
        // ff3R - funkcja celu z zewnętrzną funkcją kary.

        opt = pen(ff3R, x0, c_ex, dc_ex, epsilon, Nmax);

        // W funkcji ff3R minimalizujemy ujemny zasięg.
        //
        // Zmieniamy znak wartości funkcji celu,
        // aby wyświetlić dodatnią wartość.

        opt.y = -opt.y;

        cout << opt << endl;

        solution::clear_calls();

        // SYMULACJA DLA ZNALEZIONYCH PARAMETRÓW
        //
        // Wektor warunków początkowych:
        //
        // Y0(0) = x(0)     = 0 m
        // Y0(1) = vx(0)    = opt.x(0)
        // Y0(2) = y(0)     = 100 m
        // Y0(3) = vy(0)    = 0 m/s

        matrix Y0 = matrix(4, new double[4] { 0, opt.x(0), 100, 0 });

        // Rozwiązujemy równania ruchu.
        //
        // t0 = 0 s
        // dt = 0.01 s
        // tend = 7 s
        // opt.x(1) - optymalna prędkość kątowa

        matrix* Y = solve_ode(df3, 0, 0.01, 7, Y0, NAN, opt.x(1));

        // Szukamy:
        //
        // i50 - indeksu punktu, w którym y jest najbliższe 50 m
        // i0  - indeksu punktu, w którym y jest najbliższe 0 m

        int n = get_len(Y[0]), i50 = 0, i0 = 0;

        for (int i = 0; i < n; ++i)
        {
            // Szukamy punktu na wysokości około 50 m.

            if (abs(Y[1](i, 2) - 50) < abs(Y[1](i50, 2) - 50))
                i50 = i;

            // Szukamy punktu uderzenia w podłoże,
            // czyli położenia y ≈ 0.

            if (abs(Y[1](i, 2)) < abs(Y[1](i0, 2)))
                i0 = i;
        }

        // Współrzędna pozioma przy y ≈ 50 m.
        // Powinna należeć do przedziału [3,7] m.

        cout << endl << "x50 = " << Y[1](i50, 0) << endl;

        // Zasięg poziomy piłki w momencie uderzenia w podłoże.

        cout << endl << "xend = " << Y[1](i0, 0) << endl;

        // Zapisujemy symulację do pliku CSV.
        //
        // Kolejne kolumny:
        // t, x, vx, y, vy

        ofstream Sout("symulacja_lab3.csv");

        Sout << hcat(Y[0], Y[1]);

        Sout.close();

        Y[0].~matrix();
        Y[1].~matrix();
    }

    // ==========================================
    // CZĘŚĆ 3 - TEST MODELU PIŁKI
    // ==========================================

    if (1) // Test modelu
    {
        // Parametry testowe:
        //
        // v0x = 5 m/s
        // omega = 10 rad/s
        //
        // Nie wykonujemy tutaj optymalizacji.
        // Sprawdzamy wyłącznie model matematyczny.

        matrix Y0 = matrix(4, new double[4] { 0, 5, 100, 0 });

        matrix* Y = solve_ode(df3, 0, 0.01, 7, Y0, NAN, 10);

        // Zapisujemy trajektorię do pliku CSV.

        ofstream Sout("symulacja_lab3.csv");

        Sout << hcat(Y[0], Y[1]);

        Sout.close();

        Y[0].~matrix();
        Y[1].~matrix();
    }
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
